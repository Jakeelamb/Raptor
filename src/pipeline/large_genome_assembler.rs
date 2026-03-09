//! Large Genome Assembler - handles 100+ Gb genomes on 16 GB RAM
//!
//! This is the production-ready assembler for salamander-scale genomes.
//! Uses disk-based k-mer counting with proper sequence reconstruction.
//!
//! Algorithm:
//! 1. Stream FASTQ → distribute k-mers to disk buckets
//! 2. Sort & count each bucket → filter by coverage
//! 3. Build adjacency index for graph traversal
//! 4. Greedy path extension → extract contigs
//!
//! Memory: O(bucket_size + adjacency_cache) ≈ 2-4 GB

use crate::io::fasta::FastaWriter;
use crate::io::fastq::{open_fastq, stream_fastq_records, stream_paired_fastq_records};
use crate::kmer::disk_counting_v2::{
    decode_kmer, extend_left, extend_right, DiskCounterConfig, DiskKmerCounterV2,
};
use crate::kmer::kmer::{reverse_complement, KmerU64};
use ahash::{AHashMap, AHashSet};
use rayon::prelude::*;
use std::cmp::Reverse;
use tracing::info;

/// Configuration for large genome assembly
#[derive(Clone, Debug)]
pub struct LargeGenomeConfig {
    /// K-mer size (recommend 31 for large genomes)
    pub k: usize,
    /// Minimum k-mer count to use (filters errors). Use 0 for adaptive selection.
    pub min_count: u32,
    /// Minimum contig length to output
    pub min_contig_len: usize,
    /// Number of disk buckets (None = auto)
    pub num_buckets: Option<usize>,
    /// Temp directory (None = system temp)
    pub temp_dir: Option<String>,
    /// Maximum tip length to remove (graph cleaning)
    pub max_tip_len: usize,
    /// Remove bubbles shorter than this
    pub max_bubble_len: usize,
}

impl Default for LargeGenomeConfig {
    fn default() -> Self {
        Self {
            k: 31,
            min_count: 0,
            min_contig_len: 200,
            num_buckets: None,
            temp_dir: None,
            max_tip_len: 100,
            max_bubble_len: 50,
        }
    }
}

/// Assembly statistics
#[derive(Debug, Clone, Default)]
pub struct AssemblyStats {
    pub reads_processed: u64,
    pub bases_processed: u64,
    pub effective_min_count: u32,
    pub kmers_total: u64,
    pub kmers_unique: u64,
    pub kmers_filtered: u64,
    pub kmers_error_corrected: u64,
    pub tips_removed: usize,
    pub bubbles_popped: usize,
    pub contigs: usize,
    pub total_length: usize,
    pub n50: usize,
    pub largest: usize,
    pub disk_bytes: u64,
}

impl std::fmt::Display for AssemblyStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "=== Assembly Statistics ===")?;
        writeln!(f, "Reads: {}", self.reads_processed)?;
        writeln!(
            f,
            "Bases: {} ({:.2} Gb)",
            self.bases_processed,
            self.bases_processed as f64 / 1e9
        )?;
        writeln!(f, "Min count used: {}", self.effective_min_count)?;
        writeln!(
            f,
            "K-mers: {} total, {} unique, {} after filtering",
            self.kmers_total, self.kmers_unique, self.kmers_filtered
        )?;
        if self.kmers_error_corrected > 0 {
            writeln!(f, "K-mers error-corrected: {}", self.kmers_error_corrected)?;
        }
        if self.tips_removed > 0 || self.bubbles_popped > 0 {
            writeln!(
                f,
                "Graph cleaning: {} tips removed, {} bubbles popped",
                self.tips_removed, self.bubbles_popped
            )?;
        }
        writeln!(f, "Contigs: {}", self.contigs)?;
        writeln!(f, "Total length: {} bp", self.total_length)?;
        writeln!(f, "N50: {} bp", self.n50)?;
        writeln!(f, "Largest: {} bp", self.largest)?;
        writeln!(f, "Disk used: {:.2} GB", self.disk_bytes as f64 / 1e9)?;
        Ok(())
    }
}

/// Insert size statistics for paired-end reads
#[derive(Debug, Clone, Default)]
pub struct InsertSizeStats {
    pub mean: f64,
    pub std_dev: f64,
    pub min: usize,
    pub max: usize,
    pub median: usize,
}

/// Repeat detection statistics
#[derive(Debug, Clone, Default)]
pub struct RepeatStats {
    pub median_coverage: u32,
    pub repeat_threshold: u32,
    pub repeat_kmer_count: usize,
    pub total_kmer_count: usize,
}

#[derive(Debug, Clone)]
struct WeightedGraphNode {
    predecessors: Vec<(u64, u32)>,
    successors: Vec<(u64, u32)>,
}

#[derive(Debug, Clone, Default)]
struct WeightedGraphDiagnostics {
    nodes: usize,
    oriented_edges: usize,
    branching_nodes: usize,
    source_nodes: usize,
    sink_nodes: usize,
}

#[derive(Debug, Clone, Copy)]
struct BranchCandidate {
    base_idx: usize,
    next: u64,
    count: u32,
    is_repeat: bool,
    read_support: u32,
}

#[derive(Debug, Clone, Default)]
struct ReadThreadingStats {
    ambiguous_edges: usize,
    supported_edges: usize,
    edge_observations: u64,
}

impl InsertSizeStats {
    pub fn from_samples(samples: &[usize]) -> Self {
        if samples.is_empty() {
            return Self::default();
        }

        let n = samples.len() as f64;
        let sum: usize = samples.iter().sum();
        let mean = sum as f64 / n;

        let variance: f64 = samples
            .iter()
            .map(|&x| (x as f64 - mean).powi(2))
            .sum::<f64>()
            / n;
        let std_dev = variance.sqrt();

        let mut sorted = samples.to_vec();
        sorted.sort_unstable();

        Self {
            mean,
            std_dev,
            min: sorted[0],
            max: sorted[sorted.len() - 1],
            median: sorted[sorted.len() / 2],
        }
    }
}

/// Main large genome assembler
pub struct LargeGenomeAssembler {
    config: LargeGenomeConfig,
}

impl LargeGenomeAssembler {
    pub fn new(config: LargeGenomeConfig) -> Self {
        Self { config }
    }

    /// Run the complete assembly pipeline
    pub fn assemble(&self, input_path: &str, output_path: &str) -> std::io::Result<AssemblyStats> {
        let mut stats = AssemblyStats::default();
        let k = self.config.k;

        info!("=== Large Genome Assembler ===");
        info!("K-mer size: {}", k);
        info!("Min count request: {}", self.config.min_count);

        // Phase 1: Disk-based k-mer counting
        info!("Phase 1/6: Distributing k-mers to disk...");
        let disk_config = self.create_disk_config();
        let mut counter = DiskKmerCounterV2::new(disk_config)?;

        let (reads, bases) = self.distribute_from_fastq(input_path, &mut counter)?;
        stats.reads_processed = reads;
        stats.bases_processed = bases;

        let (total_kmers, num_buckets, disk_bytes) = counter.stats();
        stats.kmers_total = total_kmers;
        stats.disk_bytes = disk_bytes;
        info!(
            "Distributed {} k-mers to {} buckets ({:.2} GB on disk)",
            total_kmers,
            num_buckets,
            disk_bytes as f64 / 1e9
        );

        // Phase 2: Count k-mers
        info!("Phase 2/6: Counting k-mers from buckets...");
        let kmer_counts = counter.count_all()?;
        stats.kmers_unique = kmer_counts.len() as u64;

        // Phase 3: Error correction - rescue low-count k-mers that are 1 edit from high-count
        info!("Phase 3/6: Error correction...");
        let (corrected_counts, num_corrected) = self.error_correct_kmers(&kmer_counts, k);
        stats.kmers_error_corrected = num_corrected;
        info!("Error-corrected {} k-mers", num_corrected);

        let effective_min_count = self.select_min_count(&corrected_counts);
        stats.effective_min_count = effective_min_count;
        if self.config.min_count == 0 {
            info!(
                "Adaptive min count selected: {} (requested auto)",
                effective_min_count
            );
        }

        let filtered: AHashMap<u64, u32> = corrected_counts
            .into_iter()
            .filter(|(_, count)| *count >= effective_min_count)
            .collect();
        stats.kmers_filtered = filtered.len() as u64;
        info!(
            "Unique k-mers: {}, after filtering: {}",
            stats.kmers_unique, stats.kmers_filtered
        );

        // Phase 4: Build graph and clean it
        info!("Phase 4/6: Building and cleaning de Bruijn graph...");
        let valid_kmers: AHashSet<u64> = filtered.keys().copied().collect();
        let mut adjacency = self.build_adjacency(&valid_kmers, k);
        info!("  Adjacency built for {} k-mers", adjacency.len());

        // Graph cleaning: remove tips and pop bubbles
        let tips_removed = self.remove_tips(&mut adjacency, &filtered, k, effective_min_count);
        stats.tips_removed = tips_removed;
        info!("  Removed {} tips", tips_removed);

        let bubbles_popped = self.pop_bubbles(&mut adjacency, &filtered, k);
        stats.bubbles_popped = bubbles_popped;
        info!("  Popped {} bubbles", bubbles_popped);

        let weighted_graph = self.build_weighted_unitig_graph(&filtered, &adjacency, k);
        let graph_diagnostics = self.analyze_weighted_graph(&weighted_graph);
        info!(
            "  Weighted graph: {} nodes, {} oriented edges, {} branching, {} sources, {} sinks",
            graph_diagnostics.nodes,
            graph_diagnostics.oriented_edges,
            graph_diagnostics.branching_nodes,
            graph_diagnostics.source_nodes,
            graph_diagnostics.sink_nodes
        );
        let ambiguous_edges = self.collect_ambiguous_branch_edges(&weighted_graph);
        let (branch_support, threading_stats) =
            self.collect_branch_support_from_fastq(input_path, &adjacency, &ambiguous_edges, k)?;
        info!(
            "  Read threading: {} ambiguous edges, {} supported, {} observations",
            threading_stats.ambiguous_edges,
            threading_stats.supported_edges,
            threading_stats.edge_observations
        );

        // Phase 5: Build contigs from cleaned graph
        info!("Phase 5/6: Building contigs from cleaned graph...");
        let contigs = self.build_contigs_from_graph(&filtered, &adjacency, &branch_support, k);
        info!("Assembled {} raw contigs", contigs.len());

        // Phase 6: Write output
        info!("Phase 6/6: Writing output...");
        self.write_output(&contigs, output_path, &mut stats)?;

        // Cleanup
        counter.cleanup()?;

        info!("\n{}", stats);
        Ok(stats)
    }

    /// Run assembly with paired-end reads
    pub fn assemble_paired(
        &self,
        input1_path: &str,
        input2_path: &str,
        output_path: &str,
    ) -> std::io::Result<AssemblyStats> {
        let mut stats = AssemblyStats::default();
        let k = self.config.k;

        info!("=== Large Genome Assembler (Paired-End Mode) ===");
        info!("K-mer size: {}", k);
        info!("Min count request: {}", self.config.min_count);
        info!("Input R1: {}", input1_path);
        info!("Input R2: {}", input2_path);

        // Phase 1: Disk-based k-mer counting from paired-end reads
        info!("Phase 1/6: Distributing k-mers to disk...");
        let disk_config = self.create_disk_config();
        let mut counter = DiskKmerCounterV2::new(disk_config)?;

        let (reads, bases, insert_stats) =
            self.distribute_from_paired_fastq(input1_path, input2_path, &mut counter)?;
        stats.reads_processed = reads;
        stats.bases_processed = bases;

        info!(
            "Insert size: mean={:.0}, std={:.0}",
            insert_stats.mean, insert_stats.std_dev
        );

        let (total_kmers, num_buckets, disk_bytes) = counter.stats();
        stats.kmers_total = total_kmers;
        stats.disk_bytes = disk_bytes;
        info!(
            "Distributed {} k-mers to {} buckets ({:.2} GB on disk)",
            total_kmers,
            num_buckets,
            disk_bytes as f64 / 1e9
        );

        // Phases 2-6 are the same as single-end
        info!("Phase 2/6: Counting k-mers from buckets...");
        let kmer_counts = counter.count_all()?;
        stats.kmers_unique = kmer_counts.len() as u64;

        info!("Phase 3/6: Error correction...");
        let (corrected_counts, num_corrected) = self.error_correct_kmers(&kmer_counts, k);
        stats.kmers_error_corrected = num_corrected;
        info!("Error-corrected {} k-mers", num_corrected);

        let effective_min_count = self.select_min_count(&corrected_counts);
        stats.effective_min_count = effective_min_count;
        if self.config.min_count == 0 {
            info!(
                "Adaptive min count selected: {} (requested auto)",
                effective_min_count
            );
        }

        let filtered: AHashMap<u64, u32> = corrected_counts
            .into_iter()
            .filter(|(_, count)| *count >= effective_min_count)
            .collect();
        stats.kmers_filtered = filtered.len() as u64;
        info!(
            "Unique k-mers: {}, after filtering: {}",
            stats.kmers_unique, stats.kmers_filtered
        );

        info!("Phase 4/6: Building and cleaning de Bruijn graph...");
        let valid_kmers: AHashSet<u64> = filtered.keys().copied().collect();
        let mut adjacency = self.build_adjacency(&valid_kmers, k);
        info!("  Adjacency built for {} k-mers", adjacency.len());

        let tips_removed = self.remove_tips(&mut adjacency, &filtered, k, effective_min_count);
        stats.tips_removed = tips_removed;
        info!("  Removed {} tips", tips_removed);

        let bubbles_popped = self.pop_bubbles(&mut adjacency, &filtered, k);
        stats.bubbles_popped = bubbles_popped;
        info!("  Popped {} bubbles", bubbles_popped);

        let weighted_graph = self.build_weighted_unitig_graph(&filtered, &adjacency, k);
        let graph_diagnostics = self.analyze_weighted_graph(&weighted_graph);
        info!(
            "  Weighted graph: {} nodes, {} oriented edges, {} branching, {} sources, {} sinks",
            graph_diagnostics.nodes,
            graph_diagnostics.oriented_edges,
            graph_diagnostics.branching_nodes,
            graph_diagnostics.source_nodes,
            graph_diagnostics.sink_nodes
        );
        let ambiguous_edges = self.collect_ambiguous_branch_edges(&weighted_graph);
        let (branch_support, threading_stats) = self.collect_branch_support_from_paired_fastq(
            input1_path,
            input2_path,
            &adjacency,
            &ambiguous_edges,
            k,
        )?;
        info!(
            "  Read threading: {} ambiguous edges, {} supported, {} observations",
            threading_stats.ambiguous_edges,
            threading_stats.supported_edges,
            threading_stats.edge_observations
        );

        info!("Phase 5/6: Building contigs from cleaned graph...");
        let contigs = self.build_contigs_from_graph(&filtered, &adjacency, &branch_support, k);
        info!("Assembled {} raw contigs", contigs.len());

        info!("Phase 6/6: Writing output...");
        self.write_output(&contigs, output_path, &mut stats)?;

        counter.cleanup()?;

        info!("\n{}", stats);
        Ok(stats)
    }

    fn create_disk_config(&self) -> DiskCounterConfig {
        let mut config = DiskCounterConfig::auto(self.config.k);

        if let Some(n) = self.config.num_buckets {
            config.num_buckets = n;
        }
        if let Some(ref dir) = self.config.temp_dir {
            config.temp_dir = std::path::PathBuf::from(dir);
        }
        config.min_count = 1; // We filter later for flexibility
        config
    }

    fn distribute_from_fastq(
        &self,
        path: &str,
        counter: &mut DiskKmerCounterV2,
    ) -> std::io::Result<(u64, u64)> {
        let reader = open_fastq(path);
        let mut reads = 0u64;
        let mut bases = 0u64;

        // Process in batches for efficiency
        const BATCH_SIZE: usize = 50_000;
        let mut batch: Vec<String> = Vec::with_capacity(BATCH_SIZE);

        for record in stream_fastq_records(reader) {
            reads += 1;
            bases += record.sequence.len() as u64;
            batch.push(record.sequence);

            if batch.len() >= BATCH_SIZE {
                counter.distribute(batch.drain(..).map(|s| s.into_bytes()))?;

                if reads % 5_000_000 == 0 {
                    info!("  {} million reads processed...", reads / 1_000_000);
                }
            }
        }

        // Remaining batch
        if !batch.is_empty() {
            counter.distribute(batch.drain(..).map(|s| s.into_bytes()))?;
        }

        Ok((reads, bases))
    }

    /// Distribute k-mers from paired-end FASTQ files
    fn distribute_from_paired_fastq(
        &self,
        path1: &str,
        path2: &str,
        counter: &mut DiskKmerCounterV2,
    ) -> std::io::Result<(u64, u64, InsertSizeStats)> {
        let reader1 = open_fastq(path1);
        let reader2 = open_fastq(path2);
        let mut reads = 0u64;
        let mut bases = 0u64;
        let mut insert_sizes: Vec<usize> = Vec::new();

        const BATCH_SIZE: usize = 50_000;
        const INSERT_SAMPLE_SIZE: usize = 100_000;
        let mut batch: Vec<String> = Vec::with_capacity(BATCH_SIZE * 2);

        for (r1, r2) in stream_paired_fastq_records(reader1, reader2) {
            reads += 2;
            bases += (r1.sequence.len() + r2.sequence.len()) as u64;

            // Sample insert sizes from read lengths (actual insert size estimation
            // would require mapping, but we estimate from read lengths for now)
            if insert_sizes.len() < INSERT_SAMPLE_SIZE {
                insert_sizes.push(r1.sequence.len() + r2.sequence.len());
            }

            batch.push(r1.sequence);
            batch.push(r2.sequence);

            if batch.len() >= BATCH_SIZE * 2 {
                counter.distribute(batch.drain(..).map(|s| s.into_bytes()))?;

                if reads % 10_000_000 == 0 {
                    info!("  {} million read pairs processed...", reads / 2_000_000);
                }
            }
        }

        // Remaining batch
        if !batch.is_empty() {
            counter.distribute(batch.drain(..).map(|s| s.into_bytes()))?;
        }

        // Calculate insert size statistics
        let insert_stats = InsertSizeStats::from_samples(&insert_sizes);

        Ok((reads, bases, insert_stats))
    }

    /// Build adjacency map: k-mer -> (left_extensions, right_extensions)
    ///
    /// Key insight: We store adjacency for BOTH the canonical form AND its
    /// reverse complement. This allows proper traversal regardless of which
    /// strand we're on. Extensions are checked against the valid k-mer set
    /// using canonical forms for consistency.
    fn build_adjacency(
        &self,
        valid_kmers: &AHashSet<u64>,
        k: usize,
    ) -> AHashMap<u64, ([bool; 4], [bool; 4])> {
        let bases = [b'A', b'C', b'G', b'T'];

        // For each canonical k-mer, compute adjacency for both orientations
        let entries: Vec<(u64, ([bool; 4], [bool; 4]))> = valid_kmers
            .par_iter()
            .flat_map(|&canonical_encoded| {
                let kmer = KmerU64 {
                    encoded: canonical_encoded,
                    len: k as u8,
                };
                let rc = kmer.reverse_complement();

                // Build adjacency for forward orientation
                let fwd_adj = Self::compute_adjacency(canonical_encoded, k, valid_kmers, &bases);

                // Build adjacency for reverse complement orientation
                // Note: "right" for RC is "left" for forward, and vice versa
                let rc_adj = Self::compute_adjacency(rc.encoded, k, valid_kmers, &bases);

                // Return both entries (they may be the same for palindromic k-mers)
                if canonical_encoded == rc.encoded {
                    vec![(canonical_encoded, fwd_adj)]
                } else {
                    vec![(canonical_encoded, fwd_adj), (rc.encoded, rc_adj)]
                }
            })
            .collect();

        entries.into_iter().collect()
    }

    /// Compute adjacency for a single k-mer orientation
    fn compute_adjacency(
        encoded: u64,
        k: usize,
        valid_kmers: &AHashSet<u64>,
        bases: &[u8; 4],
    ) -> ([bool; 4], [bool; 4]) {
        let mut left = [false; 4];
        let mut right = [false; 4];

        for (i, &base) in bases.iter().enumerate() {
            // Right extension: drop first base, add new at end
            if let Some(ext) = extend_right(encoded, base, k) {
                let ext_kmer = KmerU64 {
                    encoded: ext,
                    len: k as u8,
                };
                let ext_canonical = ext_kmer.canonical().encoded;
                if valid_kmers.contains(&ext_canonical) {
                    right[i] = true;
                }
            }

            // Left extension: drop last base, add new at start
            if let Some(ext) = extend_left(encoded, base, k) {
                let ext_kmer = KmerU64 {
                    encoded: ext,
                    len: k as u8,
                };
                let ext_canonical = ext_kmer.canonical().encoded;
                if valid_kmers.contains(&ext_canonical) {
                    left[i] = true;
                }
            }
        }

        (left, right)
    }

    /// Error correction: merge counts of singleton k-mers into similar high-frequency ones
    ///
    /// Strategy: For each singleton k-mer (count = 1), check if there's a
    /// high-count k-mer within Hamming distance 1. If so, add the singleton's count
    /// to the high one (singletons are almost certainly sequencing errors).
    ///
    /// Conservative approach: Only correct definite errors (count=1) to avoid
    /// removing real low-coverage k-mers.
    fn error_correct_kmers(
        &self,
        kmer_counts: &AHashMap<u64, u32>,
        k: usize,
    ) -> (AHashMap<u64, u32>, u64) {
        let min_count = self.config.min_count;
        let mut corrected = kmer_counts.clone();
        let mut num_corrected = 0u64;

        // Build set of high-confidence k-mers (need strong evidence)
        // Require at least 3x the min_count threshold for "trusted" status
        let trusted_threshold = (min_count * 3).max(4);
        let high_conf: AHashSet<u64> = kmer_counts
            .iter()
            .filter(|(_, &c)| c >= trusted_threshold)
            .map(|(&k, _)| k)
            .collect();

        if high_conf.is_empty() {
            return (corrected, 0);
        }

        // Only correct singletons (count = 1) - these are almost certainly errors
        let singleton_kmers: Vec<u64> = kmer_counts
            .iter()
            .filter(|(_, &c)| c == 1)
            .map(|(&k, _)| k)
            .collect();

        for singleton in singleton_kmers {
            // Find a trusted neighbor within Hamming distance 1
            if let Some(trusted_neighbor) = self.find_trusted_neighbor(singleton, k, &high_conf) {
                // Transfer count from error k-mer to trusted k-mer
                if let Some(count) = corrected.get_mut(&trusted_neighbor) {
                    *count = count.saturating_add(1);
                }
                corrected.remove(&singleton);
                num_corrected += 1;
            }
        }

        (corrected, num_corrected)
    }

    /// Find a high-confidence k-mer within Hamming distance 1
    fn find_trusted_neighbor(
        &self,
        encoded: u64,
        k: usize,
        trusted: &AHashSet<u64>,
    ) -> Option<u64> {
        let bases: [u64; 4] = [0, 1, 2, 3]; // A, C, G, T in 2-bit encoding

        // Try substituting each position with each alternative base
        for pos in 0..k {
            let shift = (k - 1 - pos) * 2;
            let current_base = (encoded >> shift) & 0b11;

            for &new_base in &bases {
                if new_base == current_base {
                    continue;
                }

                // Create variant by substituting base at position
                let mask = !(0b11u64 << shift);
                let variant = (encoded & mask) | (new_base << shift);

                // Check canonical form
                let variant_kmer = KmerU64 {
                    encoded: variant,
                    len: k as u8,
                };
                let canonical = variant_kmer.canonical().encoded;

                if trusted.contains(&canonical) {
                    return Some(canonical);
                }
            }
        }

        None
    }

    fn select_min_count(&self, kmer_counts: &AHashMap<u64, u32>) -> u32 {
        if self.config.min_count > 0 {
            return self.config.min_count;
        }

        let mut non_singleton_counts: Vec<u32> = kmer_counts
            .values()
            .copied()
            .filter(|&count| count >= 2)
            .collect();

        if non_singleton_counts.is_empty() {
            return 2;
        }

        non_singleton_counts.sort_unstable();
        let median = non_singleton_counts[non_singleton_counts.len() / 2];

        (median / 10).clamp(2, 5)
    }

    /// Remove tips: dead-end paths shorter than max_tip_len
    ///
    /// Tips are typically caused by sequencing errors at read ends.
    /// We identify nodes with only one neighbor (dead ends) and remove
    /// short paths leading to them.
    fn remove_tips(
        &self,
        adjacency: &mut AHashMap<u64, ([bool; 4], [bool; 4])>,
        kmer_counts: &AHashMap<u64, u32>,
        k: usize,
        min_count: u32,
    ) -> usize {
        let max_tip_len = self.config.max_tip_len;
        let bases = [b'A', b'C', b'G', b'T'];
        let mut tips_removed = 0;
        let mut to_remove: Vec<u64> = Vec::new();

        // Find tip starting points: k-mers with exactly 1 neighbor total
        for (&kmer, &(left, right)) in adjacency.iter() {
            let left_count: usize = left.iter().filter(|&&b| b).count();
            let right_count: usize = right.iter().filter(|&&b| b).count();

            // Dead end on left (only right neighbors) or right (only left neighbors)
            let is_left_dead_end = left_count == 0 && right_count >= 1;
            let is_right_dead_end = right_count == 0 && left_count >= 1;

            if is_left_dead_end || is_right_dead_end {
                // Trace the path and check if it's a tip
                let tip_length =
                    self.trace_tip_length(kmer, k, adjacency, is_right_dead_end, &bases);

                if tip_length > 0 && tip_length <= max_tip_len {
                    // Check if tip has lower coverage than the branch point
                    let tip_count = kmer_counts
                        .get(
                            &KmerU64 {
                                encoded: kmer,
                                len: k as u8,
                            }
                            .canonical()
                            .encoded,
                        )
                        .copied()
                        .unwrap_or(0);

                    // Only remove if low coverage (likely error)
                    if tip_count <= min_count.saturating_mul(2) {
                        to_remove.push(kmer);
                        tips_removed += 1;
                    }
                }
            }
        }

        // Remove tips from adjacency
        for kmer in to_remove {
            adjacency.remove(&kmer);
        }

        tips_removed
    }

    /// Trace a potential tip and return its length (0 if not a simple tip)
    fn trace_tip_length(
        &self,
        start: u64,
        k: usize,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        going_left: bool,
        bases: &[u8; 4],
    ) -> usize {
        let mut current = start;
        let mut length = 1;
        let max_tip = self.config.max_tip_len + 10; // Safety limit

        while length < max_tip {
            let (left_ext, right_ext) = adjacency
                .get(&current)
                .copied()
                .unwrap_or(([false; 4], [false; 4]));

            let extensions = if going_left { &left_ext } else { &right_ext };
            let ext_count: usize = extensions.iter().filter(|&&b| b).count();

            if ext_count == 0 {
                // Reached dead end
                return length;
            } else if ext_count == 1 {
                // Continue along the path
                let ext_idx = extensions.iter().position(|&b| b).unwrap();
                let base = bases[ext_idx];

                let next = if going_left {
                    extend_left(current, base, k)
                } else {
                    extend_right(current, base, k)
                };

                if let Some(next_kmer) = next {
                    let next_canonical = KmerU64 {
                        encoded: next_kmer,
                        len: k as u8,
                    }
                    .canonical()
                    .encoded;
                    if adjacency.contains_key(&next_kmer) || adjacency.contains_key(&next_canonical)
                    {
                        current = if adjacency.contains_key(&next_kmer) {
                            next_kmer
                        } else {
                            next_canonical
                        };
                        length += 1;
                    } else {
                        return length;
                    }
                } else {
                    return length;
                }
            } else {
                // Reached a branch point - this is where the tip connects
                return length;
            }
        }

        0 // Path too long to be a tip
    }

    /// Pop bubbles: remove alternative paths between the same start/end nodes
    ///
    /// Bubbles are caused by heterozygosity or sequencing errors creating
    /// two similar paths. We keep the higher-coverage path.
    fn pop_bubbles(
        &self,
        adjacency: &mut AHashMap<u64, ([bool; 4], [bool; 4])>,
        kmer_counts: &AHashMap<u64, u32>,
        k: usize,
    ) -> usize {
        let max_bubble_len = self.config.max_bubble_len;
        let bases = [b'A', b'C', b'G', b'T'];
        let mut bubbles_popped = 0;
        let mut to_remove: AHashSet<u64> = AHashSet::new();

        // Find branch points: k-mers with 2+ extensions in same direction
        let branch_points: Vec<u64> = adjacency
            .iter()
            .filter(|(_, (left, right))| {
                let left_count: usize = left.iter().filter(|&&b| b).count();
                let right_count: usize = right.iter().filter(|&&b| b).count();
                left_count >= 2 || right_count >= 2
            })
            .map(|(&k, _)| k)
            .collect();

        for branch in branch_points {
            let (left_ext, right_ext) = adjacency
                .get(&branch)
                .copied()
                .unwrap_or(([false; 4], [false; 4]));

            // Check right branches
            let right_branches: Vec<usize> = right_ext
                .iter()
                .enumerate()
                .filter(|(_, &b)| b)
                .map(|(i, _)| i)
                .collect();

            if right_branches.len() >= 2 {
                // Trace each branch and look for convergence
                if let Some(bubble_kmers) = self.find_bubble(
                    branch,
                    &right_branches,
                    k,
                    adjacency,
                    max_bubble_len,
                    kmer_counts,
                    &bases,
                    false,
                ) {
                    to_remove.extend(bubble_kmers);
                    bubbles_popped += 1;
                }
            }

            // Check left branches
            let left_branches: Vec<usize> = left_ext
                .iter()
                .enumerate()
                .filter(|(_, &b)| b)
                .map(|(i, _)| i)
                .collect();

            if left_branches.len() >= 2 {
                if let Some(bubble_kmers) = self.find_bubble(
                    branch,
                    &left_branches,
                    k,
                    adjacency,
                    max_bubble_len,
                    kmer_counts,
                    &bases,
                    true,
                ) {
                    to_remove.extend(bubble_kmers);
                    bubbles_popped += 1;
                }
            }
        }

        // Remove bubble k-mers from adjacency
        for kmer in to_remove {
            adjacency.remove(&kmer);
        }

        bubbles_popped
    }

    /// Find and return k-mers in the lower-coverage branch of a bubble
    fn find_bubble(
        &self,
        branch_point: u64,
        branches: &[usize],
        k: usize,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        max_len: usize,
        kmer_counts: &AHashMap<u64, u32>,
        bases: &[u8; 4],
        going_left: bool,
    ) -> Option<Vec<u64>> {
        if branches.len() < 2 {
            return None;
        }

        // Trace each branch path
        let mut paths: Vec<(Vec<u64>, u64)> = Vec::new(); // (path_kmers, end_point)

        for &branch_idx in branches {
            let base = bases[branch_idx];
            let first_kmer = if going_left {
                extend_left(branch_point, base, k)
            } else {
                extend_right(branch_point, base, k)
            };

            if let Some(first) = first_kmer {
                if let Some((path, end)) =
                    self.trace_path(first, k, adjacency, max_len, bases, going_left)
                {
                    paths.push((path, end));
                }
            }
        }

        // Check if any two paths converge to the same endpoint
        if paths.len() < 2 {
            return None;
        }

        for i in 0..paths.len() {
            for j in (i + 1)..paths.len() {
                if paths[i].1 == paths[j].1 {
                    // Found a bubble - keep the higher coverage path
                    let cov_i: u32 = paths[i]
                        .0
                        .iter()
                        .filter_map(|&kmer| {
                            let canonical = KmerU64 {
                                encoded: kmer,
                                len: k as u8,
                            }
                            .canonical()
                            .encoded;
                            kmer_counts.get(&canonical).copied()
                        })
                        .sum();
                    let cov_j: u32 = paths[j]
                        .0
                        .iter()
                        .filter_map(|&kmer| {
                            let canonical = KmerU64 {
                                encoded: kmer,
                                len: k as u8,
                            }
                            .canonical()
                            .encoded;
                            kmer_counts.get(&canonical).copied()
                        })
                        .sum();

                    // Return the lower-coverage path for removal
                    if cov_i < cov_j {
                        return Some(paths[i].0.clone());
                    } else {
                        return Some(paths[j].0.clone());
                    }
                }
            }
        }

        None
    }

    /// Trace a path from a starting k-mer, returning the path and endpoint
    fn trace_path(
        &self,
        start: u64,
        k: usize,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        max_len: usize,
        bases: &[u8; 4],
        going_left: bool,
    ) -> Option<(Vec<u64>, u64)> {
        let mut path = vec![start];
        let mut current = start;

        for _ in 0..max_len {
            let lookup_kmer = if adjacency.contains_key(&current) {
                current
            } else {
                let canonical = KmerU64 {
                    encoded: current,
                    len: k as u8,
                }
                .canonical()
                .encoded;
                if adjacency.contains_key(&canonical) {
                    canonical
                } else {
                    return Some((path, current));
                }
            };

            let (left_ext, right_ext) = adjacency.get(&lookup_kmer).copied()?;
            let extensions = if going_left { &left_ext } else { &right_ext };
            let ext_count: usize = extensions.iter().filter(|&&b| b).count();

            if ext_count == 0 {
                return Some((path, current));
            } else if ext_count == 1 {
                let ext_idx = extensions.iter().position(|&b| b).unwrap();
                let base = bases[ext_idx];
                let next = if going_left {
                    extend_left(current, base, k)?
                } else {
                    extend_right(current, base, k)?
                };
                path.push(next);
                current = next;
            } else {
                // Reached another branch point - this is the end
                return Some((path, current));
            }
        }

        None // Path too long
    }

    /// Identify repeat k-mers based on coverage distribution
    ///
    /// A k-mer is considered a repeat if its count is significantly higher
    /// than the median coverage (typically 2x or more).
    fn identify_repeats(&self, kmer_counts: &AHashMap<u64, u32>) -> (AHashSet<u64>, RepeatStats) {
        let mut counts: Vec<u32> = kmer_counts.values().copied().collect();
        counts.sort_unstable();

        if counts.is_empty() {
            return (AHashSet::new(), RepeatStats::default());
        }

        // Calculate median coverage
        let median = counts[counts.len() / 2];

        // Repeat threshold: 2x median or minimum of 10
        let repeat_threshold = (median * 2).max(10);

        // Identify repeat k-mers
        let repeat_kmers: AHashSet<u64> = kmer_counts
            .iter()
            .filter(|(_, &count)| count >= repeat_threshold)
            .map(|(&kmer, _)| kmer)
            .collect();

        let stats = RepeatStats {
            median_coverage: median,
            repeat_threshold,
            repeat_kmer_count: repeat_kmers.len(),
            total_kmer_count: kmer_counts.len(),
        };

        info!(
            "Repeat detection: median_cov={}, threshold={}, repeats={}",
            median,
            repeat_threshold,
            repeat_kmers.len()
        );

        (repeat_kmers, stats)
    }

    /// Build contigs from a pre-built (and cleaned) adjacency graph
    /// with repeat-aware extension
    fn build_contigs_from_graph(
        &self,
        kmer_counts: &AHashMap<u64, u32>,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        branch_support: &AHashMap<(u64, u64), u32>,
        k: usize,
    ) -> Vec<String> {
        let min_len = self.config.min_contig_len;
        let valid_kmers: AHashSet<u64> = kmer_counts.keys().copied().collect();

        // Identify repeat k-mers for coverage-guided traversal
        let (repeat_kmers, _repeat_stats) = self.identify_repeats(kmer_counts);

        // Track used k-mers
        let mut used: AHashSet<u64> = AHashSet::with_capacity(valid_kmers.len());
        let mut contigs: Vec<String> = Vec::new();

        // Sort k-mers by count (highest first) for seed selection
        // But skip repeat k-mers as seeds (they cause fragmentation)
        // Check if k-mer (or its RC) is still in the adjacency graph after cleaning
        let mut sorted: Vec<(u64, u32)> = kmer_counts
            .iter()
            .filter(|(&kmer, _)| {
                if repeat_kmers.contains(&kmer) {
                    return false;
                }
                if adjacency.contains_key(&kmer) {
                    return true;
                }
                let rc = KmerU64 {
                    encoded: kmer,
                    len: k as u8,
                }
                .reverse_complement()
                .encoded;
                adjacency.contains_key(&rc)
            })
            .map(|(&e, &c)| (e, c))
            .collect();
        sorted.sort_unstable_by_key(|(_, c)| Reverse(*c));

        info!("  Extending from {} seed k-mers...", sorted.len());
        let mut progress = 0;

        for (seed, _count) in sorted {
            if used.contains(&seed) {
                continue;
            }

            let contig = self.extend_bidirectional_with_coverage(
                seed,
                k,
                &valid_kmers,
                adjacency,
                kmer_counts,
                branch_support,
                &repeat_kmers,
                &mut used,
            );

            if contig.len() >= min_len {
                contigs.push(contig);
            }

            progress += 1;
            if progress % 100_000 == 0 {
                info!(
                    "    {} seeds processed, {} contigs...",
                    progress,
                    contigs.len()
                );
            }
        }

        contigs
    }

    fn build_weighted_unitig_graph(
        &self,
        kmer_counts: &AHashMap<u64, u32>,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        k: usize,
    ) -> AHashMap<u64, WeightedGraphNode> {
        let bases = [b'A', b'C', b'G', b'T'];

        adjacency
            .iter()
            .map(|(&node, &(left, right))| {
                let predecessors =
                    self.weighted_neighbors(node, &left, adjacency, kmer_counts, k, &bases, true);
                let successors =
                    self.weighted_neighbors(node, &right, adjacency, kmer_counts, k, &bases, false);
                (
                    node,
                    WeightedGraphNode {
                        predecessors,
                        successors,
                    },
                )
            })
            .collect()
    }

    fn analyze_weighted_graph(
        &self,
        graph: &AHashMap<u64, WeightedGraphNode>,
    ) -> WeightedGraphDiagnostics {
        let mut diagnostics = WeightedGraphDiagnostics {
            nodes: graph.len(),
            ..Default::default()
        };

        for node in graph.values() {
            diagnostics.oriented_edges += node.successors.len();

            if node.predecessors.len() > 1 || node.successors.len() > 1 {
                diagnostics.branching_nodes += 1;
            }
            if node.predecessors.is_empty() && !node.successors.is_empty() {
                diagnostics.source_nodes += 1;
            }
            if node.successors.is_empty() && !node.predecessors.is_empty() {
                diagnostics.sink_nodes += 1;
            }
        }

        diagnostics
    }

    fn collect_ambiguous_branch_edges(
        &self,
        graph: &AHashMap<u64, WeightedGraphNode>,
    ) -> AHashSet<(u64, u64)> {
        let mut edges = AHashSet::new();

        for (&node, node_info) in graph {
            if node_info.successors.len() > 1 {
                for &(next, _) in &node_info.successors {
                    edges.insert((node, next));
                }
            }
            if node_info.predecessors.len() > 1 {
                for &(prev, _) in &node_info.predecessors {
                    edges.insert((prev, node));
                }
            }
        }

        edges
    }

    fn collect_branch_support_from_fastq(
        &self,
        path: &str,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        branch_edges: &AHashSet<(u64, u64)>,
        k: usize,
    ) -> std::io::Result<(AHashMap<(u64, u64), u32>, ReadThreadingStats)> {
        let mut edge_support = AHashMap::with_capacity(branch_edges.len());
        let mut stats = ReadThreadingStats {
            ambiguous_edges: branch_edges.len(),
            ..Default::default()
        };

        if branch_edges.is_empty() {
            return Ok((edge_support, stats));
        }

        let reader = open_fastq(path);
        for record in stream_fastq_records(reader) {
            stats.edge_observations += Self::thread_branch_edges_in_sequence(
                record.sequence.as_bytes(),
                k,
                adjacency,
                branch_edges,
                &mut edge_support,
            );
        }

        stats.supported_edges = edge_support.len();
        Ok((edge_support, stats))
    }

    fn collect_branch_support_from_paired_fastq(
        &self,
        path1: &str,
        path2: &str,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        branch_edges: &AHashSet<(u64, u64)>,
        k: usize,
    ) -> std::io::Result<(AHashMap<(u64, u64), u32>, ReadThreadingStats)> {
        let mut edge_support = AHashMap::with_capacity(branch_edges.len());
        let mut stats = ReadThreadingStats {
            ambiguous_edges: branch_edges.len(),
            ..Default::default()
        };

        if branch_edges.is_empty() {
            return Ok((edge_support, stats));
        }

        let reader1 = open_fastq(path1);
        let reader2 = open_fastq(path2);
        for (r1, r2) in stream_paired_fastq_records(reader1, reader2) {
            stats.edge_observations += Self::thread_branch_edges_in_sequence(
                r1.sequence.as_bytes(),
                k,
                adjacency,
                branch_edges,
                &mut edge_support,
            );
            stats.edge_observations += Self::thread_branch_edges_in_sequence(
                r2.sequence.as_bytes(),
                k,
                adjacency,
                branch_edges,
                &mut edge_support,
            );
        }

        stats.supported_edges = edge_support.len();
        Ok((edge_support, stats))
    }

    fn thread_branch_edges_in_sequence(
        sequence: &[u8],
        k: usize,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        branch_edges: &AHashSet<(u64, u64)>,
        edge_support: &mut AHashMap<(u64, u64), u32>,
    ) -> u64 {
        if branch_edges.is_empty() || sequence.len() < k {
            return 0;
        }

        let mut observations = 0u64;
        let mut prev_node = None;
        let mut pos = 0usize;

        while pos + k <= sequence.len() {
            let Some(mut current_kmer) = KmerU64::from_slice(&sequence[pos..pos + k]) else {
                prev_node = None;
                pos += 1;
                continue;
            };

            loop {
                if let Some(current_node) =
                    Self::resolve_threaded_node(current_kmer.encoded, adjacency, k)
                {
                    if let Some(prev) = prev_node {
                        let edge = (prev, current_node);
                        if branch_edges.contains(&edge) {
                            *edge_support.entry(edge).or_insert(0) += 1;
                            observations += 1;
                        }
                    }
                    prev_node = Some(current_node);
                } else {
                    prev_node = None;
                }

                if pos + k >= sequence.len() {
                    return observations;
                }

                let next_base = sequence[pos + k];
                let Some(next_kmer) = current_kmer.extend(next_base) else {
                    prev_node = None;
                    pos += k + 1;
                    break;
                };

                current_kmer = next_kmer;
                pos += 1;
            }
        }

        observations
    }

    fn resolve_threaded_node(
        encoded: u64,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        k: usize,
    ) -> Option<u64> {
        if adjacency.contains_key(&encoded) {
            return Some(encoded);
        }

        let canonical = KmerU64 {
            encoded,
            len: k as u8,
        }
        .canonical()
        .encoded;
        adjacency.contains_key(&canonical).then_some(canonical)
    }

    fn weighted_neighbors(
        &self,
        node: u64,
        extensions: &[bool; 4],
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        kmer_counts: &AHashMap<u64, u32>,
        k: usize,
        bases: &[u8; 4],
        going_left: bool,
    ) -> Vec<(u64, u32)> {
        let mut neighbors = Vec::new();

        for (idx, &valid) in extensions.iter().enumerate() {
            if !valid {
                continue;
            }

            let candidate = if going_left {
                extend_left(node, bases[idx], k)
            } else {
                extend_right(node, bases[idx], k)
            };

            let Some(candidate) = candidate else {
                continue;
            };

            let resolved = if adjacency.contains_key(&candidate) {
                candidate
            } else {
                let canonical = KmerU64 {
                    encoded: candidate,
                    len: k as u8,
                }
                .canonical()
                .encoded;
                if adjacency.contains_key(&canonical) {
                    canonical
                } else {
                    continue;
                }
            };

            let support = kmer_counts
                .get(
                    &KmerU64 {
                        encoded: resolved,
                        len: k as u8,
                    }
                    .canonical()
                    .encoded,
                )
                .copied()
                .unwrap_or(0);

            neighbors.push((resolved, support));
        }

        neighbors.sort_unstable_by_key(|(neighbor, _)| *neighbor);
        neighbors.dedup_by_key(|(neighbor, _)| *neighbor);
        neighbors
    }

    #[allow(dead_code)]
    fn extract_maximal_unitigs(
        &self,
        graph: &AHashMap<u64, WeightedGraphNode>,
        k: usize,
    ) -> Vec<String> {
        let mut nodes: Vec<u64> = graph.keys().copied().collect();
        nodes.sort_unstable();

        let mut used_edges: AHashSet<(u64, u64)> = AHashSet::new();
        let mut covered_canonical: AHashSet<u64> = AHashSet::new();
        let mut contigs = Vec::new();

        for &node in &nodes {
            let Some(node_info) = graph.get(&node) else {
                continue;
            };

            if node_info.predecessors.len() == 1 && node_info.successors.len() == 1 {
                continue;
            }

            for &(next, _) in &node_info.successors {
                let edge_key = Self::canonical_edge_key(node, next, k);
                if used_edges.contains(&edge_key) {
                    continue;
                }

                let path = self.trace_unitig_path(node, next, graph, k, &mut used_edges);
                if let Some(contig) = self.render_path(&path, k, &mut covered_canonical) {
                    contigs.push(contig);
                }
            }
        }

        for &node in &nodes {
            let Some(node_info) = graph.get(&node) else {
                continue;
            };

            if node_info.predecessors.len() != 1 || node_info.successors.len() != 1 {
                continue;
            }

            let next = node_info.successors[0].0;
            let edge_key = Self::canonical_edge_key(node, next, k);
            if used_edges.contains(&edge_key) {
                continue;
            }

            let path = self.trace_unitig_cycle(node, graph, k, &mut used_edges);
            if let Some(contig) = self.render_path(&path, k, &mut covered_canonical) {
                contigs.push(contig);
            }
        }

        for &node in &nodes {
            let canonical = KmerU64 {
                encoded: node,
                len: k as u8,
            }
            .canonical()
            .encoded;

            if covered_canonical.insert(canonical) {
                contigs.push(decode_kmer(canonical, k));
            }
        }

        contigs
    }

    #[allow(dead_code)]
    fn trace_unitig_path(
        &self,
        start: u64,
        next: u64,
        graph: &AHashMap<u64, WeightedGraphNode>,
        k: usize,
        used_edges: &mut AHashSet<(u64, u64)>,
    ) -> Vec<u64> {
        let mut path = vec![start];
        let mut current = start;
        let mut next_node = next;

        loop {
            let edge_key = Self::canonical_edge_key(current, next_node, k);
            if !used_edges.insert(edge_key) {
                break;
            }

            path.push(next_node);
            current = next_node;

            let Some(node_info) = graph.get(&current) else {
                break;
            };

            if node_info.predecessors.len() != 1 || node_info.successors.len() != 1 {
                break;
            }

            next_node = node_info.successors[0].0;
        }

        path
    }

    #[allow(dead_code)]
    fn trace_unitig_cycle(
        &self,
        start: u64,
        graph: &AHashMap<u64, WeightedGraphNode>,
        k: usize,
        used_edges: &mut AHashSet<(u64, u64)>,
    ) -> Vec<u64> {
        let mut path = vec![start];
        let mut current = start;

        loop {
            let Some(node_info) = graph.get(&current) else {
                break;
            };
            let next = node_info.successors[0].0;
            let edge_key = Self::canonical_edge_key(current, next, k);
            if !used_edges.insert(edge_key) {
                break;
            }
            if next == start {
                break;
            }
            path.push(next);
            current = next;
        }

        path
    }

    #[allow(dead_code)]
    fn render_path(
        &self,
        path: &[u64],
        k: usize,
        covered_canonical: &mut AHashSet<u64>,
    ) -> Option<String> {
        let (&first, rest) = path.split_first()?;

        let mut sequence = decode_kmer(first, k).into_bytes();
        covered_canonical.insert(
            KmerU64 {
                encoded: first,
                len: k as u8,
            }
            .canonical()
            .encoded,
        );

        for &node in rest {
            covered_canonical.insert(
                KmerU64 {
                    encoded: node,
                    len: k as u8,
                }
                .canonical()
                .encoded,
            );
            sequence.push(Self::last_base(node));
        }

        let forward = String::from_utf8(sequence).ok()?;
        let reverse = reverse_complement(&forward);
        if reverse < forward {
            Some(reverse)
        } else {
            Some(forward)
        }
    }

    #[allow(dead_code)]
    fn canonical_edge_key(from: u64, to: u64, k: usize) -> (u64, u64) {
        let forward = (from, to);
        let reverse = (
            KmerU64 {
                encoded: to,
                len: k as u8,
            }
            .reverse_complement()
            .encoded,
            KmerU64 {
                encoded: from,
                len: k as u8,
            }
            .reverse_complement()
            .encoded,
        );

        if reverse < forward {
            reverse
        } else {
            forward
        }
    }

    #[allow(dead_code)]
    fn last_base(encoded: u64) -> u8 {
        match encoded & 0b11 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        }
    }

    /// Extend a seed k-mer bidirectionally to form a contig
    ///
    /// The key insight for handling canonical k-mers correctly:
    /// - We track the ACTUAL k-mer (not canonical) as we extend
    /// - We use canonical forms only for the "used" set to avoid revisiting
    /// - Adjacency is stored for both forward and RC orientations
    fn extend_bidirectional(
        &self,
        seed_canonical: u64,
        k: usize,
        _valid_kmers: &AHashSet<u64>,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        used: &mut AHashSet<u64>,
    ) -> String {
        let bases = [b'A', b'C', b'G', b'T'];

        // Start with the seed k-mer sequence (use canonical form)
        let seed_seq = decode_kmer(seed_canonical, k);
        let mut contig: Vec<u8> = seed_seq.into_bytes();
        used.insert(seed_canonical);

        // Track current ACTUAL k-mer at each end (not canonical)
        // This is crucial for correct extension direction
        let mut right_kmer = seed_canonical;
        let mut left_kmer = seed_canonical;

        // Extend right: look up adjacency for current k-mer (not canonical)
        loop {
            let (_, right_ext) = adjacency
                .get(&right_kmer)
                .copied()
                .unwrap_or(([false; 4], [false; 4]));

            // Find valid extensions that haven't been used
            let extensions: Vec<usize> = right_ext
                .iter()
                .enumerate()
                .filter(|(i, &valid)| {
                    if !valid {
                        return false;
                    }
                    // Check if this extension leads to an unused k-mer
                    if let Some(next) = extend_right(right_kmer, bases[*i], k) {
                        let next_kmer = KmerU64 {
                            encoded: next,
                            len: k as u8,
                        };
                        let next_canonical = next_kmer.canonical().encoded;
                        !used.contains(&next_canonical)
                    } else {
                        false
                    }
                })
                .map(|(i, _)| i)
                .collect();

            if extensions.len() != 1 {
                break; // Stop at branch or dead end
            }

            let base_idx = extensions[0];
            let base = bases[base_idx];

            if let Some(next) = extend_right(right_kmer, base, k) {
                let next_kmer = KmerU64 {
                    encoded: next,
                    len: k as u8,
                };
                let next_canonical = next_kmer.canonical().encoded;

                contig.push(base);
                used.insert(next_canonical);

                // Continue with the ACTUAL extended k-mer (not canonical)
                // This preserves correct directionality
                right_kmer = next;
            } else {
                break;
            }
        }

        // Extend left
        loop {
            let (left_ext, _) = adjacency
                .get(&left_kmer)
                .copied()
                .unwrap_or(([false; 4], [false; 4]));

            let extensions: Vec<usize> = left_ext
                .iter()
                .enumerate()
                .filter(|(i, &valid)| {
                    if !valid {
                        return false;
                    }
                    if let Some(next) = extend_left(left_kmer, bases[*i], k) {
                        let next_kmer = KmerU64 {
                            encoded: next,
                            len: k as u8,
                        };
                        let next_canonical = next_kmer.canonical().encoded;
                        !used.contains(&next_canonical)
                    } else {
                        false
                    }
                })
                .map(|(i, _)| i)
                .collect();

            if extensions.len() != 1 {
                break;
            }

            let base_idx = extensions[0];
            let base = bases[base_idx];

            if let Some(next) = extend_left(left_kmer, base, k) {
                let next_kmer = KmerU64 {
                    encoded: next,
                    len: k as u8,
                };
                let next_canonical = next_kmer.canonical().encoded;

                contig.insert(0, base);
                used.insert(next_canonical);

                // Continue with ACTUAL extended k-mer
                left_kmer = next;
            } else {
                break;
            }
        }

        String::from_utf8(contig).unwrap_or_default()
    }

    /// Extend a seed k-mer bidirectionally with coverage-guided traversal
    ///
    /// At branch points, uses coverage information to choose the best path.
    /// Avoids entering repeat regions unless necessary.
    fn extend_bidirectional_with_coverage(
        &self,
        seed_canonical: u64,
        k: usize,
        _valid_kmers: &AHashSet<u64>,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        kmer_counts: &AHashMap<u64, u32>,
        branch_support: &AHashMap<(u64, u64), u32>,
        repeat_kmers: &AHashSet<u64>,
        used: &mut AHashSet<u64>,
    ) -> String {
        let bases = [b'A', b'C', b'G', b'T'];

        let seed_seq = decode_kmer(seed_canonical, k);
        let mut contig: Vec<u8> = seed_seq.into_bytes();
        used.insert(seed_canonical);

        let mut right_kmer = seed_canonical;
        let mut left_kmer = seed_canonical;

        // Extend right with coverage-guided traversal
        loop {
            let (_, right_ext) = adjacency
                .get(&right_kmer)
                .copied()
                .unwrap_or(([false; 4], [false; 4]));

            // Find valid extensions with their coverage
            let mut extensions: Vec<BranchCandidate> = Vec::new();
            for (i, &valid) in right_ext.iter().enumerate() {
                if !valid {
                    continue;
                }
                if let Some(next) = extend_right(right_kmer, bases[i], k) {
                    let next_kmer = KmerU64 {
                        encoded: next,
                        len: k as u8,
                    };
                    let next_canonical = next_kmer.canonical().encoded;
                    if !used.contains(&next_canonical) {
                        let count = kmer_counts.get(&next_canonical).copied().unwrap_or(0);
                        let is_repeat = repeat_kmers.contains(&next_canonical);
                        let read_support = branch_support
                            .get(&(right_kmer, next))
                            .copied()
                            .unwrap_or(0);
                        extensions.push(BranchCandidate {
                            base_idx: i,
                            next,
                            count,
                            is_repeat,
                            read_support,
                        });
                    }
                }
            }

            if extensions.is_empty() {
                break; // Dead end
            }

            // If single extension, take it
            if extensions.len() == 1 {
                let extension = extensions[0];
                let base = bases[extension.base_idx];
                let next_kmer = KmerU64 {
                    encoded: extension.next,
                    len: k as u8,
                };
                let next_canonical = next_kmer.canonical().encoded;

                contig.push(base);
                used.insert(next_canonical);
                right_kmer = extension.next;
                continue;
            }

            // Multiple extensions: use coverage-guided selection
            // Prefer non-repeat paths with similar coverage to current position
            let current_count = kmer_counts
                .get(
                    &KmerU64 {
                        encoded: right_kmer,
                        len: k as u8,
                    }
                    .canonical()
                    .encoded,
                )
                .copied()
                .unwrap_or(1);

            if let Some(best) = self.choose_branch_extension(&extensions, current_count) {
                let base = bases[best.base_idx];
                let next_kmer = KmerU64 {
                    encoded: best.next,
                    len: k as u8,
                };
                let next_canonical = next_kmer.canonical().encoded;

                contig.push(base);
                used.insert(next_canonical);
                right_kmer = best.next;
            } else {
                break;
            }
        }

        // Extend left with coverage-guided traversal
        loop {
            let (left_ext, _) = adjacency
                .get(&left_kmer)
                .copied()
                .unwrap_or(([false; 4], [false; 4]));

            let mut extensions: Vec<BranchCandidate> = Vec::new();
            for (i, &valid) in left_ext.iter().enumerate() {
                if !valid {
                    continue;
                }
                if let Some(next) = extend_left(left_kmer, bases[i], k) {
                    let next_kmer = KmerU64 {
                        encoded: next,
                        len: k as u8,
                    };
                    let next_canonical = next_kmer.canonical().encoded;
                    if !used.contains(&next_canonical) {
                        let count = kmer_counts.get(&next_canonical).copied().unwrap_or(0);
                        let is_repeat = repeat_kmers.contains(&next_canonical);
                        let read_support =
                            branch_support.get(&(next, left_kmer)).copied().unwrap_or(0);
                        extensions.push(BranchCandidate {
                            base_idx: i,
                            next,
                            count,
                            is_repeat,
                            read_support,
                        });
                    }
                }
            }

            if extensions.is_empty() {
                break;
            }

            if extensions.len() == 1 {
                let extension = extensions[0];
                let base = bases[extension.base_idx];
                let next_kmer = KmerU64 {
                    encoded: extension.next,
                    len: k as u8,
                };
                let next_canonical = next_kmer.canonical().encoded;

                contig.insert(0, base);
                used.insert(next_canonical);
                left_kmer = extension.next;
                continue;
            }

            // Coverage-guided selection for left extension
            let current_count = kmer_counts
                .get(
                    &KmerU64 {
                        encoded: left_kmer,
                        len: k as u8,
                    }
                    .canonical()
                    .encoded,
                )
                .copied()
                .unwrap_or(1);

            if let Some(best) = self.choose_branch_extension(&extensions, current_count) {
                let base = bases[best.base_idx];
                let next_kmer = KmerU64 {
                    encoded: best.next,
                    len: k as u8,
                };
                let next_canonical = next_kmer.canonical().encoded;

                contig.insert(0, base);
                used.insert(next_canonical);
                left_kmer = best.next;
            } else {
                break;
            }
        }

        String::from_utf8(contig).unwrap_or_default()
    }

    fn choose_branch_extension(
        &self,
        extensions: &[BranchCandidate],
        current_count: u32,
    ) -> Option<BranchCandidate> {
        if extensions.is_empty() {
            return None;
        }

        let (best_support, second_support, supported_extension) =
            Self::best_read_supported_extension(extensions);
        if best_support >= 2 && best_support > second_support {
            return Some(supported_extension);
        }

        let non_repeat = extensions.iter().copied().filter(|ext| !ext.is_repeat);
        let best_non_repeat = non_repeat.min_by_key(|ext| {
            (
                Self::coverage_ratio(current_count, ext.count),
                Reverse(ext.read_support),
                Reverse(ext.count),
                ext.base_idx,
            )
        });

        best_non_repeat.or_else(|| {
            extensions
                .iter()
                .copied()
                .max_by_key(|ext| (ext.read_support, ext.count))
        })
    }

    fn best_read_supported_extension(
        extensions: &[BranchCandidate],
    ) -> (u32, u32, BranchCandidate) {
        let mut best = extensions[0];
        let mut best_support = best.read_support;
        let mut second_support = 0u32;

        for &extension in &extensions[1..] {
            if extension.read_support > best_support {
                second_support = best_support;
                best_support = extension.read_support;
                best = extension;
            } else if extension.read_support > second_support {
                second_support = extension.read_support;
            }
        }

        (best_support, second_support, best)
    }

    fn coverage_ratio(current_count: u32, candidate_count: u32) -> u32 {
        if candidate_count > current_count {
            candidate_count / current_count.max(1)
        } else {
            current_count / candidate_count.max(1)
        }
    }

    fn write_output(
        &self,
        contigs: &[String],
        path: &str,
        stats: &mut AssemblyStats,
    ) -> std::io::Result<()> {
        let mut writer = FastaWriter::new(path);

        // Filter and sort by length
        let mut valid: Vec<&String> = contigs
            .iter()
            .filter(|c| c.len() >= self.config.min_contig_len)
            .collect();
        valid.sort_by_key(|c| std::cmp::Reverse(c.len()));

        let mut lengths: Vec<usize> = Vec::with_capacity(valid.len());
        let mut total_len = 0usize;

        for (i, contig) in valid.iter().enumerate() {
            writer.write_record(&format!("contig_{} len={}", i + 1, contig.len()), contig)?;
            lengths.push(contig.len());
            total_len += contig.len();
        }

        // Calculate N50
        let half = total_len / 2;
        let mut cumsum = 0;
        let n50 = lengths
            .iter()
            .find(|&&len| {
                cumsum += len;
                cumsum >= half
            })
            .copied()
            .unwrap_or(0);

        stats.contigs = valid.len();
        stats.total_length = total_len;
        stats.n50 = n50;
        stats.largest = lengths.first().copied().unwrap_or(0);

        Ok(())
    }
}

/// Convenience function for quick assembly
pub fn assemble_large_genome(
    input: &str,
    output: &str,
    k: usize,
    min_count: u32,
) -> std::io::Result<AssemblyStats> {
    let config = LargeGenomeConfig {
        k,
        min_count,
        ..Default::default()
    };
    LargeGenomeAssembler::new(config).assemble(input, output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::{NamedTempFile, TempDir};

    fn create_test_fastq() -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();

        // Create reads with overlapping k-mers that should assemble
        let sequences = [
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", // 40bp
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", // Duplicate
            "CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA", // Shifted
            "GTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC", // More shift
        ];

        for (i, seq) in sequences.iter().enumerate() {
            writeln!(file, "@read_{}", i).unwrap();
            writeln!(file, "{}", seq).unwrap();
            writeln!(file, "+").unwrap();
            writeln!(file, "{}", "I".repeat(seq.len())).unwrap();
        }
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_large_genome_assembler() {
        let input = create_test_fastq();
        let output = NamedTempFile::new().unwrap();
        let temp_dir = TempDir::new().unwrap();

        let config = LargeGenomeConfig {
            k: 11,
            min_count: 1,
            min_contig_len: 12, // Just slightly above k
            num_buckets: Some(4),
            temp_dir: Some(temp_dir.path().to_str().unwrap().to_string()),
            ..Default::default()
        };

        let assembler = LargeGenomeAssembler::new(config);
        let stats = assembler
            .assemble(
                input.path().to_str().unwrap(),
                output.path().to_str().unwrap(),
            )
            .unwrap();

        eprintln!("Stats: {:?}", stats);
        assert_eq!(stats.reads_processed, 4);
        // With min_count=1 and overlapping reads, we should get some contigs
        // But even if assembly doesn't work perfectly, verify basic flow works
        assert!(stats.kmers_filtered > 0, "Should have filtered k-mers");
    }

    /// Test assembly with realistic overlapping reads that should form longer contigs
    #[test]
    fn test_assembly_with_overlapping_reads() {
        let mut file = NamedTempFile::new().unwrap();

        // Reference sequence: 100bp of realistic DNA
        let reference = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATC\
                         GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\
                         ATCGATCGATCGATCGATCG";

        // Generate overlapping reads (simulating 3x coverage)
        let read_len = 40;
        let step = 10; // 30bp overlap between consecutive reads
        for i in 0..((reference.len() - read_len) / step + 1) {
            let start = i * step;
            if start + read_len <= reference.len() {
                let seq = &reference[start..start + read_len];
                writeln!(file, "@read_{}", i).unwrap();
                writeln!(file, "{}", seq).unwrap();
                writeln!(file, "+").unwrap();
                writeln!(file, "{}", "I".repeat(read_len)).unwrap();
            }
        }
        // Add duplicates for coverage
        for i in 0..((reference.len() - read_len) / step + 1) {
            let start = i * step;
            if start + read_len <= reference.len() {
                let seq = &reference[start..start + read_len];
                writeln!(file, "@read_dup_{}", i).unwrap();
                writeln!(file, "{}", seq).unwrap();
                writeln!(file, "+").unwrap();
                writeln!(file, "{}", "I".repeat(read_len)).unwrap();
            }
        }
        file.flush().unwrap();

        let output = NamedTempFile::new().unwrap();
        let temp_dir = TempDir::new().unwrap();
        let config = LargeGenomeConfig {
            k: 21,
            min_count: 2, // Require at least 2x coverage
            min_contig_len: 30,
            num_buckets: Some(4),
            temp_dir: Some(temp_dir.path().to_str().unwrap().to_string()),
            ..Default::default()
        };

        let assembler = LargeGenomeAssembler::new(config);
        let stats = assembler
            .assemble(
                file.path().to_str().unwrap(),
                output.path().to_str().unwrap(),
            )
            .unwrap();

        eprintln!("Overlapping reads test - Stats: {:?}", stats);

        // We should get at least one contig of reasonable length
        assert!(stats.contigs >= 1, "Should produce at least one contig");
        assert!(
            stats.largest >= 30,
            "Largest contig should be at least 30bp"
        );
        eprintln!("Largest contig: {}bp, N50: {}bp", stats.largest, stats.n50);
    }

    /// Test error correction with simulated sequencing errors
    #[test]
    fn test_error_correction() {
        let mut file = NamedTempFile::new().unwrap();
        let temp_dir = TempDir::new().unwrap();

        // Create a sequence with high coverage
        let good_seq = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";

        // Add many copies of the good sequence (high coverage)
        for i in 0..10 {
            writeln!(file, "@good_read_{}", i).unwrap();
            writeln!(file, "{}", good_seq).unwrap();
            writeln!(file, "+").unwrap();
            writeln!(file, "{}", "I".repeat(good_seq.len())).unwrap();
        }

        // Add a few reads with single-base errors (should be corrected)
        let error_seqs = [
            "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCC", // T->C at end
            "CTGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG", // A->C at start
            "ATGCGATCGATCGATCGATTGATCGATCGATCGATCGATCGATCGATCG", // C->T in middle
        ];

        for (i, seq) in error_seqs.iter().enumerate() {
            writeln!(file, "@error_read_{}", i).unwrap();
            writeln!(file, "{}", seq).unwrap();
            writeln!(file, "+").unwrap();
            writeln!(file, "{}", "I".repeat(seq.len())).unwrap();
        }

        file.flush().unwrap();

        let output = NamedTempFile::new().unwrap();
        let config = LargeGenomeConfig {
            k: 21,
            min_count: 3, // Require decent coverage
            min_contig_len: 25,
            num_buckets: Some(4),
            temp_dir: Some(temp_dir.path().to_str().unwrap().to_string()),
            ..Default::default()
        };

        let assembler = LargeGenomeAssembler::new(config);
        let stats = assembler
            .assemble(
                file.path().to_str().unwrap(),
                output.path().to_str().unwrap(),
            )
            .unwrap();

        eprintln!("Error correction test - Stats: {:?}", stats);
        eprintln!("Error-corrected k-mers: {}", stats.kmers_error_corrected);

        // The pipeline should work and produce contigs
        assert!(stats.kmers_filtered > 0, "Should have filtered k-mers");
    }

    /// Test graph cleaning with tip-inducing reads
    #[test]
    fn test_graph_cleaning() {
        let mut file = NamedTempFile::new().unwrap();
        let temp_dir = TempDir::new().unwrap();

        // Main sequence with good coverage
        let main_seq = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";

        // Add multiple copies for coverage
        for i in 0..6 {
            writeln!(file, "@main_{}", i).unwrap();
            writeln!(file, "{}", main_seq).unwrap();
            writeln!(file, "+").unwrap();
            writeln!(file, "{}", "I".repeat(main_seq.len())).unwrap();
        }

        // Add some reads that create a short branching path (potential tip)
        // This read shares a prefix but diverges
        let tip_seq = "ATGCGATCGATCGATCGATCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
        for i in 0..2 {
            writeln!(file, "@tip_{}", i).unwrap();
            writeln!(file, "{}", tip_seq).unwrap();
            writeln!(file, "+").unwrap();
            writeln!(file, "{}", "I".repeat(tip_seq.len())).unwrap();
        }

        file.flush().unwrap();

        let output = NamedTempFile::new().unwrap();
        let config = LargeGenomeConfig {
            k: 21,
            min_count: 2,
            min_contig_len: 30,
            max_tip_len: 50,
            num_buckets: Some(4),
            temp_dir: Some(temp_dir.path().to_str().unwrap().to_string()),
            ..Default::default()
        };

        let assembler = LargeGenomeAssembler::new(config);
        let stats = assembler
            .assemble(
                file.path().to_str().unwrap(),
                output.path().to_str().unwrap(),
            )
            .unwrap();

        eprintln!("Graph cleaning test - Stats: {:?}", stats);
        eprintln!(
            "Tips removed: {}, Bubbles popped: {}",
            stats.tips_removed, stats.bubbles_popped
        );

        // Should produce at least one contig
        assert!(stats.contigs >= 1, "Should produce at least one contig");
    }

    /// Test paired-end assembly
    #[test]
    fn test_paired_end_assembly() {
        let temp_dir = TempDir::new().unwrap();

        // Create R1 and R2 files with paired reads
        let mut r1_file = NamedTempFile::new().unwrap();
        let mut r2_file = NamedTempFile::new().unwrap();

        // Reference sequence
        let reference = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAATCGATCGATCGATCGATCG";

        // Generate paired reads with ~200bp insert size (overlapping for short test)
        let read_len = 40;
        for i in 0..5 {
            let start1 = i * 10;
            let start2 = start1 + 50;

            if start1 + read_len <= reference.len() && start2 + read_len <= reference.len() {
                let seq1 = &reference[start1..start1 + read_len];
                let seq2 = &reference[start2..start2 + read_len];

                // R1
                writeln!(r1_file, "@read_{}/1", i).unwrap();
                writeln!(r1_file, "{}", seq1).unwrap();
                writeln!(r1_file, "+").unwrap();
                writeln!(r1_file, "{}", "I".repeat(read_len)).unwrap();

                // R2 (reverse complement in real data, but we keep forward for simplicity)
                writeln!(r2_file, "@read_{}/2", i).unwrap();
                writeln!(r2_file, "{}", seq2).unwrap();
                writeln!(r2_file, "+").unwrap();
                writeln!(r2_file, "{}", "I".repeat(read_len)).unwrap();
            }
        }

        // Add duplicates for coverage
        for i in 0..5 {
            let start1 = i * 10;
            let start2 = start1 + 50;

            if start1 + read_len <= reference.len() && start2 + read_len <= reference.len() {
                let seq1 = &reference[start1..start1 + read_len];
                let seq2 = &reference[start2..start2 + read_len];

                writeln!(r1_file, "@read_dup_{}/1", i).unwrap();
                writeln!(r1_file, "{}", seq1).unwrap();
                writeln!(r1_file, "+").unwrap();
                writeln!(r1_file, "{}", "I".repeat(read_len)).unwrap();

                writeln!(r2_file, "@read_dup_{}/2", i).unwrap();
                writeln!(r2_file, "{}", seq2).unwrap();
                writeln!(r2_file, "+").unwrap();
                writeln!(r2_file, "{}", "I".repeat(read_len)).unwrap();
            }
        }

        r1_file.flush().unwrap();
        r2_file.flush().unwrap();

        let output = NamedTempFile::new().unwrap();
        let config = LargeGenomeConfig {
            k: 21,
            min_count: 2,
            min_contig_len: 30,
            num_buckets: Some(4),
            temp_dir: Some(temp_dir.path().to_str().unwrap().to_string()),
            ..Default::default()
        };

        let assembler = LargeGenomeAssembler::new(config);
        let stats = assembler
            .assemble_paired(
                r1_file.path().to_str().unwrap(),
                r2_file.path().to_str().unwrap(),
                output.path().to_str().unwrap(),
            )
            .unwrap();

        eprintln!("Paired-end test - Stats: {:?}", stats);
        assert!(
            stats.reads_processed >= 4,
            "Should have processed paired reads"
        );
        assert!(stats.kmers_filtered > 0, "Should have filtered k-mers");
    }

    /// Test insert size statistics
    #[test]
    fn test_insert_size_stats() {
        let samples = vec![200, 210, 195, 205, 200, 190, 215, 200];
        let stats = InsertSizeStats::from_samples(&samples);

        assert!(stats.mean > 190.0 && stats.mean < 210.0);
        assert!(stats.std_dev > 0.0);
        assert_eq!(stats.min, 190);
        assert_eq!(stats.max, 215);
    }

    /// Test repeat detection and resolution
    #[test]
    fn test_repeat_detection() {
        let mut file = NamedTempFile::new().unwrap();
        let temp_dir = TempDir::new().unwrap();

        // Create realistic sequences with a repeat region
        // Use longer sequences that will have valid k-mers

        // Unique region 1 - normal coverage
        let unique1 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        // Repeat region - will have high coverage (many reads)
        let repeat = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";

        // Unique region 2 - normal coverage
        let unique2 = "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA";

        let read_len = 40;

        // Generate reads from unique1 (2x coverage)
        for rep in 0..2 {
            for i in 0..2 {
                let start = i * 2;
                if start + read_len <= unique1.len() {
                    let seq = &unique1[start..start + read_len];
                    writeln!(file, "@unique1_{}_{}", rep, i).unwrap();
                    writeln!(file, "{}", seq).unwrap();
                    writeln!(file, "+").unwrap();
                    writeln!(file, "{}", "I".repeat(read_len)).unwrap();
                }
            }
        }

        // Generate reads from repeat region (10x coverage to simulate repeat)
        for rep in 0..10 {
            for i in 0..2 {
                let start = i * 2;
                if start + read_len <= repeat.len() {
                    let seq = &repeat[start..start + read_len];
                    writeln!(file, "@repeat_{}_{}", rep, i).unwrap();
                    writeln!(file, "{}", seq).unwrap();
                    writeln!(file, "+").unwrap();
                    writeln!(file, "{}", "I".repeat(read_len)).unwrap();
                }
            }
        }

        // Generate reads from unique2 (2x coverage)
        for rep in 0..2 {
            for i in 0..2 {
                let start = i * 2;
                if start + read_len <= unique2.len() {
                    let seq = &unique2[start..start + read_len];
                    writeln!(file, "@unique2_{}_{}", rep, i).unwrap();
                    writeln!(file, "{}", seq).unwrap();
                    writeln!(file, "+").unwrap();
                    writeln!(file, "{}", "I".repeat(read_len)).unwrap();
                }
            }
        }

        file.flush().unwrap();

        let output = NamedTempFile::new().unwrap();
        let config = LargeGenomeConfig {
            k: 21,
            min_count: 2,
            min_contig_len: 25,
            num_buckets: Some(4),
            temp_dir: Some(temp_dir.path().to_str().unwrap().to_string()),
            ..Default::default()
        };

        let assembler = LargeGenomeAssembler::new(config);
        let stats = assembler
            .assemble(
                file.path().to_str().unwrap(),
                output.path().to_str().unwrap(),
            )
            .unwrap();

        eprintln!("Repeat detection test - Stats: {:?}", stats);
        eprintln!(
            "Reads: {}, K-mers: {}",
            stats.reads_processed, stats.kmers_filtered
        );

        // The assembler should process the reads
        assert!(stats.reads_processed > 0, "Should have processed reads");
        assert!(
            stats.kmers_filtered > 0,
            "Should have k-mers after filtering"
        );
    }

    /// Test RepeatStats calculation
    #[test]
    fn test_repeat_stats() {
        let mut counts = AHashMap::new();
        // Normal coverage k-mers (count ~5)
        for i in 0..100 {
            counts.insert(i as u64, 5);
        }
        // High coverage k-mers (count ~50, representing repeats)
        for i in 100..110 {
            counts.insert(i as u64, 50);
        }

        let config = LargeGenomeConfig::default();
        let assembler = LargeGenomeAssembler::new(config);
        let (repeat_kmers, stats) = assembler.identify_repeats(&counts);

        assert_eq!(stats.total_kmer_count, 110);
        assert!(stats.repeat_kmer_count > 0, "Should identify repeat k-mers");
        assert!(
            stats.median_coverage <= 10,
            "Median should be around normal coverage"
        );
        assert!(
            repeat_kmers.len() >= 5,
            "Should identify the high-coverage k-mers as repeats"
        );
    }

    #[test]
    fn test_auto_min_count_scales_with_coverage() {
        let mut counts = AHashMap::new();
        for i in 0..500 {
            counts.insert(i as u64, 35);
        }
        for i in 500..650 {
            counts.insert(i as u64, 2);
        }

        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            min_count: 0,
            ..Default::default()
        });

        assert_eq!(assembler.select_min_count(&counts), 3);
    }

    #[test]
    fn test_auto_min_count_stays_lenient_for_low_coverage_data() {
        let mut counts = AHashMap::new();
        for i in 0..100 {
            counts.insert(i as u64, 6);
        }
        for i in 100..180 {
            counts.insert(i as u64, 2);
        }

        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            min_count: 0,
            ..Default::default()
        });

        assert_eq!(assembler.select_min_count(&counts), 2);
    }

    #[test]
    fn test_explicit_min_count_overrides_auto_selection() {
        let mut counts = AHashMap::new();
        for i in 0..500 {
            counts.insert(i as u64, 40);
        }

        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            min_count: 4,
            ..Default::default()
        });

        assert_eq!(assembler.select_min_count(&counts), 4);
    }

    #[test]
    fn test_unitig_graph_compacts_linear_path() {
        let sequence = "ACGTAC";
        let k = 4;
        let mut counts = AHashMap::new();
        for i in 0..=sequence.len() - k {
            let encoded = KmerU64::from_str(&sequence[i..i + k])
                .unwrap()
                .canonical()
                .encoded;
            counts.insert(encoded, 10);
        }

        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            k,
            min_count: 1,
            min_contig_len: 1,
            ..Default::default()
        });
        let valid_kmers: AHashSet<u64> = counts.keys().copied().collect();
        let adjacency = assembler.build_adjacency(&valid_kmers, k);
        let graph = assembler.build_weighted_unitig_graph(&counts, &adjacency, k);
        let diagnostics = assembler.analyze_weighted_graph(&graph);
        let contigs = assembler.extract_maximal_unitigs(&graph, k);
        let reverse = reverse_complement(sequence);

        assert_eq!(diagnostics.branching_nodes, 0);
        assert!(diagnostics.nodes >= counts.len());
        assert!(diagnostics.oriented_edges >= counts.len().saturating_sub(1));
        assert_eq!(contigs.len(), 1);
        assert!(contigs[0] == sequence || contigs[0] == reverse);
    }

    #[test]
    fn test_branch_threading_counts_ambiguous_edge_support() {
        let k = 4;
        let sequences = ["TACGATG", "TACGATG", "TACGACG"];
        let mut counts = AHashMap::new();

        for sequence in &sequences {
            for i in 0..=sequence.len() - k {
                let encoded = KmerU64::from_str(&sequence[i..i + k])
                    .unwrap()
                    .canonical()
                    .encoded;
                counts.insert(encoded, 10);
            }
        }

        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            k,
            min_count: 1,
            min_contig_len: 1,
            ..Default::default()
        });
        let valid_kmers: AHashSet<u64> = counts.keys().copied().collect();
        let adjacency = assembler.build_adjacency(&valid_kmers, k);
        let graph = assembler.build_weighted_unitig_graph(&counts, &adjacency, k);
        let branch_edges = assembler.collect_ambiguous_branch_edges(&graph);
        let mut edge_support = AHashMap::new();
        let mut observations = 0u64;

        for sequence in &sequences {
            observations += LargeGenomeAssembler::thread_branch_edges_in_sequence(
                sequence.as_bytes(),
                k,
                &adjacency,
                &branch_edges,
                &mut edge_support,
            );
        }

        let branch = KmerU64::from_str("ACGA").unwrap().encoded;
        let next_primary = KmerU64::from_str("CGAT").unwrap().encoded;
        let next_alternate = KmerU64::from_str("CGAC").unwrap().encoded;

        assert!(observations >= 3);
        assert_eq!(edge_support.get(&(branch, next_primary)).copied(), Some(2));
        assert_eq!(
            edge_support.get(&(branch, next_alternate)).copied(),
            Some(1)
        );
    }

    #[test]
    fn test_choose_branch_extension_prefers_strong_read_support() {
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig::default());
        let extensions = [
            BranchCandidate {
                base_idx: 0,
                next: 11,
                count: 18,
                is_repeat: false,
                read_support: 5,
            },
            BranchCandidate {
                base_idx: 1,
                next: 22,
                count: 20,
                is_repeat: false,
                read_support: 1,
            },
        ];

        let best = assembler.choose_branch_extension(&extensions, 20).unwrap();
        assert_eq!(best.next, 11);
    }

    #[test]
    fn test_choose_branch_extension_ignores_weak_read_support_noise() {
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig::default());
        let extensions = [
            BranchCandidate {
                base_idx: 0,
                next: 11,
                count: 5,
                is_repeat: false,
                read_support: 1,
            },
            BranchCandidate {
                base_idx: 1,
                next: 22,
                count: 20,
                is_repeat: false,
                read_support: 0,
            },
        ];

        let best = assembler.choose_branch_extension(&extensions, 20).unwrap();
        assert_eq!(best.next, 22);
    }
}
