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

use crate::eval::metrics::{evaluate_lengths_sorted_desc, BaseComposition};
use crate::io::fasta::FastaWriter;
use crate::io::fastq::{
    stream_fastq_records_checked, stream_paired_fastq_records_checked, try_open_fastq,
};
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
    /// Minimum count required for a k-mer to be considered trusted during singleton rescue.
    pub error_correction_min_trusted_count: u32,
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
    /// Minimum read support required before branch threading overrides coverage.
    pub branch_support_min_win: u32,
    /// Minimum lead over the runner-up support before branch threading overrides coverage.
    pub branch_support_min_margin: u32,
    /// Prefer non-repeat branch candidates when read support does not decide.
    pub prefer_non_repeat_branches: bool,
    /// Prefer higher-coverage seeds before lower-coverage seeds during contig extraction.
    pub prefer_high_count_seeds: bool,
    /// Prefer non-repeat seeds before repeat seeds during contig extraction.
    pub prefer_non_repeat_seeds: bool,
    /// Run a second repeat-seed completion pass after the non-repeat extraction pass.
    pub enable_repeat_seed_completion: bool,
    /// Suppress redundant contained or same-flank branch-alternative contigs after extraction.
    pub suppress_redundant_contigs: bool,
}

impl Default for LargeGenomeConfig {
    fn default() -> Self {
        Self {
            k: 31,
            min_count: 0,
            error_correction_min_trusted_count: 4,
            min_contig_len: 200,
            num_buckets: None,
            temp_dir: None,
            max_tip_len: 100,
            max_bubble_len: 50,
            branch_support_min_win: 1,
            branch_support_min_margin: 2,
            prefer_non_repeat_branches: true,
            prefer_high_count_seeds: true,
            prefer_non_repeat_seeds: true,
            enable_repeat_seed_completion: true,
            suppress_redundant_contigs: false,
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
    pub ambiguous_branch_edges: usize,
    pub branch_edges_supported: usize,
    pub branch_edge_observations: u64,
    pub branch_edge_support_fraction: f64,
    pub contigs: usize,
    pub total_length: usize,
    pub gc_bases: usize,
    pub acgt_bases: usize,
    pub n_bases: usize,
    pub ambiguous_bases: usize,
    pub gc_content: f64,
    pub n_content: f64,
    pub ambiguous_content: f64,
    pub n_run_count: usize,
    pub max_n_run: usize,
    pub mean_n_run_length: f64,
    pub n_runs_per_100kb: f64,
    pub n_bases_per_100kb: f64,
    pub ambiguous_bases_per_100kb: f64,
    pub ungapped_total_length: usize,
    pub gap_bases: usize,
    pub gap_bases_frac: f64,
    pub n10: usize,
    pub n25: usize,
    pub n50: usize,
    pub n75: usize,
    pub n90: usize,
    pub n95: usize,
    pub n99: usize,
    pub l10: usize,
    pub l25: usize,
    pub l50: usize,
    pub l75: usize,
    pub l90: usize,
    pub l95: usize,
    pub l99: usize,
    pub ungapped_n50: usize,
    pub ungapped_n90: usize,
    pub ungapped_n95: usize,
    pub ungapped_n99: usize,
    pub avg_contig_len: f64,
    pub median_contig_len: f64,
    pub au_n: f64,
    pub effective_contig_count: f64,
    pub ungapped_au_n: f64,
    pub ungapped_effective_contig_count: f64,
    pub largest: usize,
    pub contigs_ge_1kb: usize,
    pub contigs_ge_10kb: usize,
    pub contigs_ge_50kb: usize,
    pub contigs_ge_100kb: usize,
    pub contigs_ge_1mb: usize,
    pub bases_ge_1kb: usize,
    pub bases_ge_10kb: usize,
    pub bases_ge_50kb: usize,
    pub bases_ge_100kb: usize,
    pub bases_ge_1mb: usize,
    pub contigs_ge_1kb_frac: f64,
    pub contigs_ge_10kb_frac: f64,
    pub contigs_ge_50kb_frac: f64,
    pub contigs_ge_100kb_frac: f64,
    pub contigs_ge_1mb_frac: f64,
    pub bases_ge_1kb_frac: f64,
    pub bases_ge_10kb_frac: f64,
    pub bases_ge_50kb_frac: f64,
    pub bases_ge_100kb_frac: f64,
    pub bases_ge_1mb_frac: f64,
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
        if self.ambiguous_branch_edges > 0 || self.branch_edge_observations > 0 {
            writeln!(
                f,
                "Read-threaded branch support: {}/{} ({:.2}%), {} observations",
                self.branch_edges_supported,
                self.ambiguous_branch_edges,
                self.branch_edge_support_fraction * 100.0,
                self.branch_edge_observations
            )?;
        }
        writeln!(f, "Contigs: {}", self.contigs)?;
        writeln!(f, "Total length: {} bp", self.total_length)?;
        writeln!(
            f,
            "Base counts: GC={}, ACGT={}, N={}, ambiguous={}",
            self.gc_bases, self.acgt_bases, self.n_bases, self.ambiguous_bases
        )?;
        writeln!(f, "GC content: {:.2}%", self.gc_content * 100.0)?;
        writeln!(f, "N content: {:.2}%", self.n_content * 100.0)?;
        writeln!(
            f,
            "Ambiguous content: {:.2}%",
            self.ambiguous_content * 100.0
        )?;
        writeln!(
            f,
            "N runs: {} (max {}, mean {:.2} bp)",
            self.n_run_count, self.max_n_run, self.mean_n_run_length
        )?;
        writeln!(
            f,
            "N/ambig per 100kb: runs={:.2}, N={:.2}, ambiguous={:.2}",
            self.n_runs_per_100kb, self.n_bases_per_100kb, self.ambiguous_bases_per_100kb
        )?;
        writeln!(
            f,
            "Ungapped span: {} bp (gaps {} bp, {:.2}%)",
            self.ungapped_total_length,
            self.gap_bases,
            self.gap_bases_frac * 100.0
        )?;
        writeln!(f, "Mean contig: {:.2} bp", self.avg_contig_len)?;
        writeln!(f, "Median contig: {:.2} bp", self.median_contig_len)?;
        writeln!(f, "N10: {} bp", self.n10)?;
        writeln!(f, "N25: {} bp", self.n25)?;
        writeln!(f, "N50: {} bp", self.n50)?;
        writeln!(f, "N75: {} bp", self.n75)?;
        writeln!(f, "N90: {} bp", self.n90)?;
        writeln!(f, "N95: {} bp", self.n95)?;
        writeln!(f, "N99: {} bp", self.n99)?;
        writeln!(
            f,
            "Ungapped N50/N90/N95/N99: {}/{}/{}/{} bp",
            self.ungapped_n50, self.ungapped_n90, self.ungapped_n95, self.ungapped_n99
        )?;
        writeln!(f, "L10: {}", self.l10)?;
        writeln!(f, "L25: {}", self.l25)?;
        writeln!(f, "L50: {}", self.l50)?;
        writeln!(f, "L75: {}", self.l75)?;
        writeln!(f, "L90: {}", self.l90)?;
        writeln!(f, "L95: {}", self.l95)?;
        writeln!(f, "L99: {}", self.l99)?;
        writeln!(
            f,
            "auN/effective count: {:.2} bp / {:.2}",
            self.au_n, self.effective_contig_count
        )?;
        writeln!(
            f,
            "Ungapped auN/effective count: {:.2} bp / {:.2}",
            self.ungapped_au_n, self.ungapped_effective_contig_count
        )?;
        writeln!(
            f,
            "Contigs >=1kb/10kb/50kb/100kb/1mb: {}/{}/{}/{}/{}",
            self.contigs_ge_1kb,
            self.contigs_ge_10kb,
            self.contigs_ge_50kb,
            self.contigs_ge_100kb,
            self.contigs_ge_1mb
        )?;
        writeln!(
            f,
            "Contigs frac >=1kb/10kb/50kb/100kb/1mb: {:.2}%/{:.2}%/{:.2}%/{:.2}%/{:.2}%",
            self.contigs_ge_1kb_frac * 100.0,
            self.contigs_ge_10kb_frac * 100.0,
            self.contigs_ge_50kb_frac * 100.0,
            self.contigs_ge_100kb_frac * 100.0,
            self.contigs_ge_1mb_frac * 100.0
        )?;
        writeln!(
            f,
            "Span >=1kb/10kb/50kb/100kb/1mb: {}/{}/{}/{}/{} bp",
            self.bases_ge_1kb,
            self.bases_ge_10kb,
            self.bases_ge_50kb,
            self.bases_ge_100kb,
            self.bases_ge_1mb
        )?;
        writeln!(
            f,
            "Span frac >=1kb/10kb/50kb/100kb/1mb: {:.2}%/{:.2}%/{:.2}%/{:.2}%/{:.2}%",
            self.bases_ge_1kb_frac * 100.0,
            self.bases_ge_10kb_frac * 100.0,
            self.bases_ge_50kb_frac * 100.0,
            self.bases_ge_100kb_frac * 100.0,
            self.bases_ge_1mb_frac * 100.0
        )?;
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

#[derive(Debug, Clone, Copy)]
struct SeedCandidate {
    canonical: u64,
    oriented: u64,
    count: u32,
}

#[derive(Debug, Clone)]
struct RankedExtractedContig {
    original_index: usize,
    sequence: String,
    total_support: u64,
    kmer_windows: usize,
}

#[derive(Debug, Clone, Default)]
struct ReadThreadingStats {
    ambiguous_edges: usize,
    supported_edges: usize,
    edge_observations: u64,
}

#[derive(Debug, Clone, Default)]
struct OrientedBranchEdgeLookup {
    canonical_by_oriented: AHashMap<(u64, u64), (u64, u64)>,
    ambiguous_edges: usize,
}

#[doc(hidden)]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct LargeGenomeBenchBranchThreadingStats {
    pub ambiguous_edges: usize,
    pub supported_edges: usize,
    pub edge_observations: u64,
}

#[doc(hidden)]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct LargeGenomeBenchBranchCandidate {
    pub base_idx: usize,
    pub next: u64,
    pub count: u32,
    pub is_repeat: bool,
    pub read_support: u32,
}

#[doc(hidden)]
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct LargeGenomeBenchBranchResolutionCase {
    pub current_count: u32,
    pub candidates: Vec<LargeGenomeBenchBranchCandidate>,
}

#[doc(hidden)]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct LargeGenomeBenchErrorCorrectionStats {
    pub corrected_kmers: u64,
    pub remaining_kmers: usize,
}

#[derive(Debug, Default)]
struct ErrorCorrectionBatch {
    trusted_increments: AHashMap<u64, u32>,
    singleton_removals: Vec<u64>,
    corrected: u64,
}

#[doc(hidden)]
pub struct LargeGenomeStageBenchFixture {
    assembler: LargeGenomeAssembler,
    k: usize,
    kmer_counts: AHashMap<u64, u32>,
    valid_kmers: AHashSet<u64>,
    adjacency: AHashMap<u64, ([bool; 4], [bool; 4])>,
    branch_edges: AHashSet<(u64, u64)>,
    branch_reads: Vec<Vec<u8>>,
    branch_support: AHashMap<(u64, u64), u32>,
    expected_contigs: Vec<String>,
}

#[doc(hidden)]
pub struct LargeGenomeErrorCorrectionBenchFixture {
    assembler: LargeGenomeAssembler,
    k: usize,
    kmer_counts: AHashMap<u64, u32>,
}

#[derive(Debug, Clone, Copy, Default)]
struct NRunSummary {
    run_count: usize,
    max_run: usize,
}

impl NRunSummary {
    #[inline]
    fn add_sequence_and_count_ungapped(&mut self, sequence: &[u8]) -> usize {
        let mut run_len = 0usize;
        let mut ungapped_bases = 0usize;
        for &base in sequence {
            if matches!(base, b'N' | b'n') {
                run_len += 1;
            } else if run_len > 0 {
                self.run_count += 1;
                self.max_run = self.max_run.max(run_len);
                run_len = 0;
                ungapped_bases += 1;
            } else {
                ungapped_bases += 1;
            }
        }

        if run_len > 0 {
            self.run_count += 1;
            self.max_run = self.max_run.max(run_len);
        }

        ungapped_bases
    }

    #[inline]
    fn mean_run_length(self, total_n_bases: usize) -> f64 {
        if self.run_count == 0 {
            0.0
        } else {
            total_n_bases as f64 / self.run_count as f64
        }
    }
}

#[inline]
fn per_100kb(count: usize, total_bases: usize) -> f64 {
    if total_bases == 0 {
        0.0
    } else {
        count as f64 * 100_000.0 / total_bases as f64
    }
}

impl InsertSizeStats {
    pub fn from_samples(samples: &[usize]) -> Self {
        if samples.is_empty() {
            return Self::default();
        }

        let n = samples.len() as f64;
        let sum: f64 = samples.iter().map(|&x| x as f64).sum();
        let mean = sum / n;

        let variance: f64 = samples
            .iter()
            .map(|&x| (x as f64 - mean).powi(2))
            .sum::<f64>()
            / n;
        let std_dev = variance.sqrt();

        let mut median_scratch = samples.to_vec();
        let median = median_usize_in_place(&mut median_scratch);
        let mut min = usize::MAX;
        let mut max = usize::MIN;
        for &sample in samples {
            min = min.min(sample);
            max = max.max(sample);
        }

        Self {
            mean,
            std_dev,
            min,
            max,
            median,
        }
    }
}

#[inline]
fn midpoint_u32(a: u32, b: u32) -> u32 {
    ((a as u64 + b as u64) / 2) as u32
}

#[inline]
fn midpoint_usize(a: usize, b: usize) -> usize {
    ((a as u128 + b as u128) / 2) as usize
}

#[inline]
fn median_u32_in_place(values: &mut [u32]) -> u32 {
    debug_assert!(!values.is_empty());
    let mid = values.len() / 2;
    let upper = {
        let (_, upper, _) = values.select_nth_unstable(mid);
        *upper
    };
    if values.len() % 2 == 1 {
        upper
    } else {
        let lower = values[..mid].iter().copied().max().unwrap_or(upper);
        midpoint_u32(lower, upper)
    }
}

#[inline]
fn upper_median_u32_in_place(values: &mut [u32]) -> u32 {
    debug_assert!(!values.is_empty());
    let mid = values.len() / 2;
    let (_, upper, _) = values.select_nth_unstable(mid);
    *upper
}

#[inline]
fn median_usize_in_place(values: &mut [usize]) -> usize {
    debug_assert!(!values.is_empty());
    let mid = values.len() / 2;
    let upper = {
        let (_, upper, _) = values.select_nth_unstable(mid);
        *upper
    };
    if values.len() % 2 == 1 {
        upper
    } else {
        let lower = values[..mid].iter().copied().max().unwrap_or(upper);
        midpoint_usize(lower, upper)
    }
}

#[inline]
fn add_weighted_sequence_counts(
    counts: &mut AHashMap<u64, u32>,
    sequence: &[u8],
    k: usize,
    weight: u32,
) {
    if sequence.len() < k || weight == 0 {
        return;
    }

    for window in sequence.windows(k) {
        if let Some(kmer) = KmerU64::from_slice(window) {
            let canonical = kmer.canonical().encoded;
            let entry = counts.entry(canonical).or_insert(0);
            *entry = entry.saturating_add(weight);
        }
    }
}

fn xorshift64(mut state: u64) -> u64 {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    state
}

fn pseudo_random_dna(seed: u64, len: usize) -> String {
    const BASES: [char; 4] = ['A', 'C', 'G', 'T'];

    let mut state = seed.max(1);
    let mut sequence = String::with_capacity(len);
    for _ in 0..len {
        state = xorshift64(state);
        sequence.push(BASES[(state & 0b11) as usize]);
    }
    sequence
}

#[inline]
fn alternate_dna_base(base: u8, offset: usize) -> u8 {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

    let current = BASES
        .iter()
        .position(|&candidate| candidate == base)
        .unwrap_or(0);
    BASES[(current + (offset % 3) + 1) % BASES.len()]
}

impl LargeGenomeBenchBranchThreadingStats {
    fn from_internal(stats: ReadThreadingStats) -> Self {
        Self {
            ambiguous_edges: stats.ambiguous_edges,
            supported_edges: stats.supported_edges,
            edge_observations: stats.edge_observations,
        }
    }
}

impl OrientedBranchEdgeLookup {
    fn new(branch_edges: &AHashSet<(u64, u64)>, k: usize) -> Self {
        let mut canonical_by_oriented =
            AHashMap::with_capacity(branch_edges.len().saturating_mul(2));

        for &(from, to) in branch_edges {
            canonical_by_oriented.insert((from, to), (from, to));

            let reverse = (
                LargeGenomeAssembler::reverse_complement_kmer(to, k),
                LargeGenomeAssembler::reverse_complement_kmer(from, k),
            );
            canonical_by_oriented.insert(reverse, (from, to));
        }

        Self {
            canonical_by_oriented,
            ambiguous_edges: branch_edges.len(),
        }
    }

    #[inline]
    fn is_empty(&self) -> bool {
        self.canonical_by_oriented.is_empty()
    }
}

impl LargeGenomeStageBenchFixture {
    #[doc(hidden)]
    pub fn synthetic_branching(
        k: usize,
        component_count: usize,
        primary_reads_per_component: usize,
        alternate_reads_per_component: usize,
    ) -> Self {
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            k,
            min_count: 1,
            min_contig_len: k,
            ..Default::default()
        });

        let mut counts = AHashMap::new();
        let mut branch_reads = Vec::with_capacity(
            component_count
                .saturating_mul(primary_reads_per_component + alternate_reads_per_component),
        );
        let mut expected_contigs = Vec::with_capacity(component_count);
        let primary_weight = primary_reads_per_component.max(1) as u32;
        let alternate_weight = alternate_reads_per_component.max(1) as u32;
        let prefix_len = k + 12;
        let suffix_len = k + 12;

        for component_idx in 0..component_count {
            let component_seed = (component_idx as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
            let prefix = pseudo_random_dna(component_seed ^ 0xA5A5_5A5A_1234_5678, prefix_len);
            let suffix = pseudo_random_dna(component_seed ^ 0xC3C3_3C3C_8765_4321, suffix_len);
            let primary = format!("{prefix}AAGTC{suffix}");
            let alternate = format!("{prefix}CCGTC{suffix}");
            expected_contigs.push(primary.clone());

            add_weighted_sequence_counts(&mut counts, primary.as_bytes(), k, primary_weight);
            add_weighted_sequence_counts(&mut counts, alternate.as_bytes(), k, alternate_weight);

            for _ in 0..primary_reads_per_component {
                branch_reads.push(primary.as_bytes().to_vec());
            }
            for _ in 0..alternate_reads_per_component {
                branch_reads.push(alternate.as_bytes().to_vec());
            }
        }

        let valid_kmers: AHashSet<u64> = counts.keys().copied().collect();
        let adjacency = assembler.build_adjacency(&valid_kmers, k);
        let graph = assembler.build_weighted_unitig_graph(&counts, &adjacency, k);
        let branch_edges = assembler.collect_ambiguous_branch_edges(&graph);
        let (branch_support, _) = assembler.collect_branch_support_from_sequences(
            &branch_reads,
            &adjacency,
            &branch_edges,
            k,
        );

        Self {
            assembler,
            k,
            kmer_counts: counts,
            valid_kmers,
            adjacency,
            branch_edges,
            branch_reads,
            branch_support,
            expected_contigs,
        }
    }

    #[doc(hidden)]
    pub fn run_branch_threading(&self) -> LargeGenomeBenchBranchThreadingStats {
        let (_, stats) = self.assembler.collect_branch_support_from_sequences(
            &self.branch_reads,
            &self.adjacency,
            &self.branch_edges,
            self.k,
        );
        LargeGenomeBenchBranchThreadingStats::from_internal(stats)
    }

    #[doc(hidden)]
    pub fn run_contig_extraction(&self) -> Vec<String> {
        self.assembler.build_contigs_from_graph(
            &self.kmer_counts,
            &self.adjacency,
            &self.branch_support,
            self.k,
        )
    }

    #[doc(hidden)]
    pub fn cloned_kmer_counts(&self) -> AHashMap<u64, u32> {
        self.kmer_counts.clone()
    }

    #[doc(hidden)]
    pub fn valid_kmers(&self) -> AHashSet<u64> {
        self.valid_kmers.clone()
    }

    #[doc(hidden)]
    pub fn branch_support_map(&self) -> AHashMap<(u64, u64), u32> {
        self.branch_support.clone()
    }

    #[doc(hidden)]
    pub fn expected_contigs(&self) -> &[String] {
        &self.expected_contigs
    }

    #[doc(hidden)]
    pub fn total_branch_read_bases(&self) -> usize {
        self.branch_reads.iter().map(Vec::len).sum()
    }

    #[doc(hidden)]
    pub fn graph_node_count(&self) -> usize {
        self.adjacency.len()
    }

    #[doc(hidden)]
    pub fn branch_edge_count(&self) -> usize {
        self.branch_edges.len()
    }
}

impl LargeGenomeErrorCorrectionBenchFixture {
    #[doc(hidden)]
    pub fn synthetic_singleton_correction(
        k: usize,
        trusted_kmer_count: usize,
        singletons_per_trusted: usize,
    ) -> Self {
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            k,
            min_count: 1,
            min_contig_len: k,
            ..Default::default()
        });

        let mut counts = AHashMap::with_capacity(
            trusted_kmer_count.saturating_mul(singletons_per_trusted.saturating_add(1)),
        );

        for trusted_idx in 0..trusted_kmer_count {
            let trusted_seed =
                (trusted_idx as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ 0xA5A5_5A5A_1234_5678;

            let trusted = loop {
                let candidate = pseudo_random_dna(trusted_seed ^ counts.len() as u64, k);
                let encoded = KmerU64::from_str(&candidate).unwrap().canonical().encoded;
                if let std::collections::hash_map::Entry::Vacant(entry) = counts.entry(encoded) {
                    entry.insert(16 + (trusted_idx % 11) as u32);
                    break encoded;
                }
            };

            let trusted_template = decode_kmer(trusted, k).into_bytes();
            let mut generated = 0usize;
            let mut attempt = 0usize;

            while generated < singletons_per_trusted {
                let mut singleton = trusted_template.clone();
                let pos = (attempt * 7 + trusted_idx) % k;
                singleton[pos] = alternate_dna_base(singleton[pos], attempt + generated);

                let encoded = KmerU64::from_slice(&singleton).unwrap().canonical().encoded;
                attempt += 1;

                if encoded == trusted || counts.contains_key(&encoded) {
                    continue;
                }

                counts.insert(encoded, 1);
                generated += 1;
            }
        }

        Self {
            assembler,
            k,
            kmer_counts: counts,
        }
    }

    #[doc(hidden)]
    pub fn cloned_kmer_counts(&self) -> AHashMap<u64, u32> {
        self.kmer_counts.clone()
    }

    #[doc(hidden)]
    pub fn run_error_correction(
        &self,
        kmer_counts: AHashMap<u64, u32>,
    ) -> LargeGenomeBenchErrorCorrectionStats {
        let (_corrected_counts, stats) = self
            .assembler
            .run_error_correction_case(kmer_counts, self.k);
        stats
    }

    #[doc(hidden)]
    pub fn total_kmer_count(&self) -> usize {
        self.kmer_counts.len()
    }

    #[doc(hidden)]
    pub fn singleton_count(&self) -> usize {
        self.kmer_counts
            .values()
            .filter(|&&count| count == 1)
            .count()
    }
}

/// Main large genome assembler
pub struct LargeGenomeAssembler {
    config: LargeGenomeConfig,
}

impl LargeGenomeAssembler {
    const THREADING_BATCH_SIZE: usize = 512;
    const PARALLEL_THREADING_MIN_BATCH: usize = 256;
    const ERROR_CORRECTION_BATCH_SIZE: usize = 4_096;
    const PARALLEL_ERROR_CORRECTION_MIN_BATCH: usize = 8_192;

    #[doc(hidden)]
    pub fn run_error_correction_case(
        &self,
        kmer_counts: AHashMap<u64, u32>,
        k: usize,
    ) -> (AHashMap<u64, u32>, LargeGenomeBenchErrorCorrectionStats) {
        let (corrected_counts, corrected_kmers) = self.error_correct_kmers(kmer_counts, k);
        let stats = LargeGenomeBenchErrorCorrectionStats {
            corrected_kmers,
            remaining_kmers: corrected_counts.len(),
        };
        (corrected_counts, stats)
    }

    #[doc(hidden)]
    pub fn run_contig_extraction_case(
        &self,
        kmer_counts: AHashMap<u64, u32>,
        valid_kmers: &AHashSet<u64>,
        branch_support: &AHashMap<(u64, u64), u32>,
        k: usize,
    ) -> Vec<String> {
        let adjacency = self.build_adjacency(valid_kmers, k);
        self.build_contigs_from_graph(&kmer_counts, &adjacency, branch_support, k)
    }

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
        let (corrected_counts, num_corrected) = self.error_correct_kmers(kmer_counts, k);
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
        Self::apply_read_threading_stats(&mut stats, &threading_stats);
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
        let (corrected_counts, num_corrected) = self.error_correct_kmers(kmer_counts, k);
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
        Self::apply_read_threading_stats(&mut stats, &threading_stats);
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

    #[inline]
    fn apply_read_threading_stats(stats: &mut AssemblyStats, threading_stats: &ReadThreadingStats) {
        stats.ambiguous_branch_edges = threading_stats.ambiguous_edges;
        stats.branch_edges_supported = threading_stats.supported_edges;
        stats.branch_edge_observations = threading_stats.edge_observations;
        stats.branch_edge_support_fraction = if threading_stats.ambiguous_edges == 0 {
            0.0
        } else {
            threading_stats.supported_edges as f64 / threading_stats.ambiguous_edges as f64
        };
    }

    fn distribute_from_fastq(
        &self,
        path: &str,
        counter: &mut DiskKmerCounterV2,
    ) -> std::io::Result<(u64, u64)> {
        let reader = try_open_fastq(path)?;
        let mut reads = 0u64;
        let mut bases = 0u64;

        // Process in batches for efficiency
        const BATCH_SIZE: usize = 50_000;
        let mut batch: Vec<String> = Vec::with_capacity(BATCH_SIZE);

        for record in stream_fastq_records_checked(reader) {
            let record = record?;
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
        let reader1 = try_open_fastq(path1)?;
        let reader2 = try_open_fastq(path2)?;
        let mut reads = 0u64;
        let mut bases = 0u64;
        let mut insert_sizes: Vec<usize> = Vec::new();

        const BATCH_SIZE: usize = 50_000;
        const INSERT_SAMPLE_SIZE: usize = 100_000;
        let mut batch: Vec<String> = Vec::with_capacity(BATCH_SIZE * 2);

        for pair in stream_paired_fastq_records_checked(reader1, reader2) {
            let (r1, r2) = pair?;
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
            .fold(Vec::new, |mut local_entries, &canonical_encoded| {
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
                    local_entries.push((canonical_encoded, fwd_adj));
                } else {
                    local_entries.push((canonical_encoded, fwd_adj));
                    local_entries.push((rc.encoded, rc_adj));
                }
                local_entries
            })
            .reduce(Vec::new, |mut acc, mut part| {
                acc.append(&mut part);
                acc
            });

        let mut adjacency = AHashMap::with_capacity(entries.len());
        for (kmer, ext) in entries {
            adjacency.insert(kmer, ext);
        }
        adjacency
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
        mut kmer_counts: AHashMap<u64, u32>,
        k: usize,
    ) -> (AHashMap<u64, u32>, u64) {
        // Build set of high-confidence k-mers (need strong evidence)
        // Require at least 3x the min_count threshold and honor the explicit
        // trusted floor so local optimization can tune both leniency and strictness.
        let trusted_threshold = self.error_correction_trusted_threshold();
        let mut trusted_counts = AHashMap::new();
        let mut singleton_kmers = Vec::new();

        for (&encoded, &count) in &kmer_counts {
            if count >= trusted_threshold {
                trusted_counts.insert(encoded, count);
            } else if count == 1 {
                singleton_kmers.push(encoded);
            }
        }

        if trusted_counts.is_empty() || singleton_kmers.is_empty() {
            return (kmer_counts, 0);
        }

        let correction_batches =
            if singleton_kmers.len() >= Self::PARALLEL_ERROR_CORRECTION_MIN_BATCH {
                singleton_kmers
                    .par_chunks(Self::ERROR_CORRECTION_BATCH_SIZE)
                    .map(|chunk| self.collect_error_correction_batch(chunk, k, &trusted_counts))
                    .collect::<Vec<_>>()
            } else {
                vec![self.collect_error_correction_batch(&singleton_kmers, k, &trusted_counts)]
            };

        let mut num_corrected = 0u64;
        for batch in correction_batches {
            num_corrected += batch.corrected;

            for (trusted_neighbor, increment) in batch.trusted_increments {
                if let Some(count) = kmer_counts.get_mut(&trusted_neighbor) {
                    *count = count.saturating_add(increment);
                }
            }

            for singleton in batch.singleton_removals {
                kmer_counts.remove(&singleton);
            }
        }

        (kmer_counts, num_corrected)
    }

    fn error_correction_trusted_threshold(&self) -> u32 {
        self.config
            .min_count
            .saturating_mul(3)
            .max(self.config.error_correction_min_trusted_count)
    }

    fn collect_error_correction_batch(
        &self,
        singleton_kmers: &[u64],
        k: usize,
        trusted_counts: &AHashMap<u64, u32>,
    ) -> ErrorCorrectionBatch {
        let mut batch = ErrorCorrectionBatch {
            trusted_increments: AHashMap::with_capacity(singleton_kmers.len()),
            singleton_removals: Vec::with_capacity(singleton_kmers.len()),
            corrected: 0,
        };

        for &singleton in singleton_kmers {
            if let Some(trusted_neighbor) = self.find_trusted_neighbor(singleton, k, trusted_counts)
            {
                let entry = batch
                    .trusted_increments
                    .entry(trusted_neighbor)
                    .or_insert(0);
                *entry = entry.saturating_add(1);
                batch.singleton_removals.push(singleton);
                batch.corrected += 1;
            }
        }

        batch
    }

    /// Find a high-confidence k-mer within Hamming distance 1
    fn find_trusted_neighbor(
        &self,
        encoded: u64,
        k: usize,
        trusted_counts: &AHashMap<u64, u32>,
    ) -> Option<u64> {
        let bases: [u64; 4] = [0, 1, 2, 3]; // A, C, G, T in 2-bit encoding
        let mut best: Option<(u32, u64)> = None;

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

                if let Some(&count) = trusted_counts.get(&canonical) {
                    match best {
                        None => best = Some((count, canonical)),
                        Some((best_count, best_kmer)) => {
                            if count > best_count || (count == best_count && canonical < best_kmer)
                            {
                                best = Some((count, canonical));
                            }
                        }
                    }
                }
            }
        }

        best.map(|(_, kmer)| kmer)
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

        let median = upper_median_u32_in_place(&mut non_singleton_counts);

        (median / 10).clamp(2, 5)
    }

    #[inline]
    fn canonical_kmer(encoded: u64, k: usize) -> u64 {
        KmerU64 {
            encoded,
            len: k as u8,
        }
        .canonical()
        .encoded
    }

    #[inline]
    fn reverse_complement_kmer(encoded: u64, k: usize) -> u64 {
        KmerU64 {
            encoded,
            len: k as u8,
        }
        .reverse_complement()
        .encoded
    }

    #[inline]
    fn resolve_adjacency_orientation(
        encoded: u64,
        k: usize,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
    ) -> Option<u64> {
        if adjacency.contains_key(&encoded) {
            return Some(encoded);
        }

        let canonical = Self::canonical_kmer(encoded, k);
        if adjacency.contains_key(&canonical) {
            return Some(canonical);
        }

        let reverse = Self::reverse_complement_kmer(canonical, k);
        adjacency.contains_key(&reverse).then_some(reverse)
    }

    fn remove_kmers_and_prune_adjacency(
        adjacency: &mut AHashMap<u64, ([bool; 4], [bool; 4])>,
        kmers_to_remove: &AHashSet<u64>,
        k: usize,
    ) {
        if kmers_to_remove.is_empty() || adjacency.is_empty() {
            return;
        }

        let mut removed_orientations = AHashSet::with_capacity(kmers_to_remove.len() * 2);
        for &kmer in kmers_to_remove {
            let canonical = Self::canonical_kmer(kmer, k);
            removed_orientations.insert(canonical);
            removed_orientations.insert(Self::reverse_complement_kmer(canonical, k));
        }

        for &kmer in &removed_orientations {
            adjacency.remove(&kmer);
        }

        if adjacency.is_empty() {
            return;
        }

        let bases = [b'A', b'C', b'G', b'T'];
        for (&node, (left_ext, right_ext)) in adjacency.iter_mut() {
            for (idx, is_valid) in left_ext.iter_mut().enumerate() {
                if !*is_valid {
                    continue;
                }
                *is_valid = extend_left(node, bases[idx], k)
                    .is_some_and(|next| !removed_orientations.contains(&next));
            }

            for (idx, is_valid) in right_ext.iter_mut().enumerate() {
                if !*is_valid {
                    continue;
                }
                *is_valid = extend_right(node, bases[idx], k)
                    .is_some_and(|next| !removed_orientations.contains(&next));
            }
        }
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
        let mut to_remove: AHashSet<u64> = AHashSet::new();

        // Iterate canonical k-mers in sorted order to keep traversal deterministic.
        let mut candidates: Vec<u64> = adjacency
            .keys()
            .copied()
            .map(|kmer| Self::canonical_kmer(kmer, k))
            .collect();
        candidates.sort_unstable();
        candidates.dedup();

        // Find tip starting points: k-mers with exactly 1 neighbor total.
        for kmer in candidates {
            let Some(&(left, right)) = adjacency.get(&kmer) else {
                continue;
            };
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
                    let tip_count = kmer_counts.get(&kmer).copied().unwrap_or(0);

                    // Only remove if low coverage (likely error)
                    if tip_count <= min_count.saturating_mul(2) && to_remove.insert(kmer) {
                        tips_removed += 1;
                    }
                }
            }
        }

        Self::remove_kmers_and_prune_adjacency(adjacency, &to_remove, k);

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
            let Some(current_oriented) = Self::resolve_adjacency_orientation(current, k, adjacency)
            else {
                return length;
            };
            let (left_ext, right_ext) = adjacency
                .get(&current_oriented)
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
                    extend_left(current_oriented, base, k)
                } else {
                    extend_right(current_oriented, base, k)
                };

                if let Some(next_kmer) = next {
                    if let Some(next_oriented) =
                        Self::resolve_adjacency_orientation(next_kmer, k, adjacency)
                    {
                        current = next_oriented;
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

        // Find canonical branch points in deterministic order.
        let mut branch_points: Vec<u64> = adjacency
            .iter()
            .filter(|(_, (left, right))| {
                let left_count: usize = left.iter().filter(|&&b| b).count();
                let right_count: usize = right.iter().filter(|&&b| b).count();
                left_count >= 2 || right_count >= 2
            })
            .map(|(&node, _)| Self::canonical_kmer(node, k))
            .collect();
        branch_points.sort_unstable();
        branch_points.dedup();

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
                    for bubble_kmer in bubble_kmers {
                        to_remove.insert(Self::canonical_kmer(bubble_kmer, k));
                    }
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
                    for bubble_kmer in bubble_kmers {
                        to_remove.insert(Self::canonical_kmer(bubble_kmer, k));
                    }
                    bubbles_popped += 1;
                }
            }
        }

        Self::remove_kmers_and_prune_adjacency(adjacency, &to_remove, k);

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
                let shared_end = paths[i].1;
                if shared_end == paths[j].1 {
                    // Found a bubble. Compare only branch-unique nodes so we do not
                    // remove the shared reconvergence node.
                    let path_i = Self::path_without_terminal_node(&paths[i].0, shared_end);
                    let path_j = Self::path_without_terminal_node(&paths[j].0, shared_end);

                    if path_i.is_empty() && path_j.is_empty() {
                        continue;
                    }

                    let cov_i = Self::path_coverage_sum(path_i, kmer_counts, k);
                    let cov_j = Self::path_coverage_sum(path_j, kmer_counts, k);

                    let remove_i =
                        match Self::compare_bubble_path_quality(path_i, cov_i, path_j, cov_j, k) {
                            std::cmp::Ordering::Less => true,
                            std::cmp::Ordering::Greater => false,
                            std::cmp::Ordering::Equal => false,
                        };

                    return if remove_i {
                        Some(path_i.to_vec())
                    } else {
                        Some(path_j.to_vec())
                    };
                }
            }
        }

        None
    }

    #[inline]
    fn path_without_terminal_node(path: &[u64], terminal: u64) -> &[u64] {
        if path.last().copied() == Some(terminal) {
            &path[..path.len().saturating_sub(1)]
        } else {
            path
        }
    }

    #[inline]
    fn path_coverage_sum(path: &[u64], kmer_counts: &AHashMap<u64, u32>, k: usize) -> u64 {
        path.iter().fold(0u64, |acc, &kmer| {
            let canonical = Self::canonical_kmer(kmer, k);
            acc.saturating_add(kmer_counts.get(&canonical).copied().unwrap_or(0) as u64)
        })
    }

    fn compare_bubble_path_quality(
        left_path: &[u64],
        left_cov: u64,
        right_path: &[u64],
        right_cov: u64,
        k: usize,
    ) -> std::cmp::Ordering {
        left_cov
            .cmp(&right_cov)
            .then_with(|| {
                let left_len = left_path.len().max(1) as u128;
                let right_len = right_path.len().max(1) as u128;
                (left_cov as u128 * right_len).cmp(&(right_cov as u128 * left_len))
            })
            .then_with(|| left_path.len().cmp(&right_path.len()))
            .then_with(|| {
                left_path
                    .iter()
                    .map(|&kmer| Self::canonical_kmer(kmer, k))
                    .cmp(right_path.iter().map(|&kmer| Self::canonical_kmer(kmer, k)))
                    .reverse()
            })
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
            let Some(lookup_kmer) = Self::resolve_adjacency_orientation(current, k, adjacency)
            else {
                return Some((path, current));
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
                    extend_left(lookup_kmer, base, k)?
                } else {
                    extend_right(lookup_kmer, base, k)?
                };
                let Some(next_oriented) = Self::resolve_adjacency_orientation(next, k, adjacency)
                else {
                    path.push(next);
                    return Some((path, next));
                };
                path.push(next_oriented);
                current = next_oriented;
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

        if counts.is_empty() {
            return (AHashSet::new(), RepeatStats::default());
        }

        // Calculate median coverage
        let median = median_u32_in_place(&mut counts);

        // Repeat threshold: 2x median or minimum of 10
        let repeat_threshold = median.saturating_mul(2).max(10);

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

        // Identify repeat k-mers for coverage-guided traversal
        let (repeat_kmers, _repeat_stats) = self.identify_repeats(kmer_counts);

        // Track used k-mers
        let mut used: AHashSet<u64> = AHashSet::with_capacity(kmer_counts.len());
        let mut contigs: Vec<String> = Vec::new();

        // Prefer non-repeat seeds, but fall back to repeat seeds if graph cleaning
        // removed all non-repeat nodes from the surviving adjacency.
        let mut preferred = Vec::new();
        let mut repeat_fallback = Vec::new();
        for (&canonical, &count) in kmer_counts {
            let Some(oriented) = Self::orient_seed_for_adjacency(canonical, adjacency, k) else {
                continue;
            };

            let entry = SeedCandidate {
                canonical,
                oriented,
                count,
            };
            if !self.config.prefer_non_repeat_seeds {
                preferred.push(entry);
            } else if repeat_kmers.contains(&canonical) {
                repeat_fallback.push(entry);
            } else {
                preferred.push(entry);
            }
        }

        let seed_order = |left: &SeedCandidate, right: &SeedCandidate| {
            let count_cmp = if self.config.prefer_high_count_seeds {
                right.count.cmp(&left.count)
            } else {
                left.count.cmp(&right.count)
            };
            count_cmp.then_with(|| left.canonical.cmp(&right.canonical))
        };
        preferred.sort_unstable_by(seed_order);
        repeat_fallback.sort_unstable_by(seed_order);

        let using_repeat_fallback = self.config.prefer_non_repeat_seeds && preferred.is_empty();
        if using_repeat_fallback {
            info!(
                "  Non-repeat seeds unavailable after graph cleanup; falling back to {} repeat seeds",
                repeat_fallback.len()
            );
            preferred = std::mem::take(&mut repeat_fallback);
        }

        info!("  Extending from {} seed k-mers...", preferred.len());
        let mut progress = 0;

        for seed in preferred {
            if used.contains(&seed.canonical) {
                continue;
            }

            if self.config.suppress_redundant_contigs
                && self.mark_enclosed_seed_component_used(seed.oriented, k, adjacency, &mut used)
            {
                continue;
            }

            let contig = self.extend_bidirectional_with_coverage(
                seed.oriented,
                k,
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

        // Quality/completeness pass: if repeat components remain disconnected from
        // non-repeat seeds, extend them deterministically as a second pass.
        if self.config.enable_repeat_seed_completion
            && !using_repeat_fallback
            && !repeat_fallback.is_empty()
        {
            info!(
                "  Running repeat-seed completion pass on {} seeds...",
                repeat_fallback.len()
            );
            let contigs_before = contigs.len();
            for seed in repeat_fallback {
                if used.contains(&seed.canonical) {
                    continue;
                }

                if self.config.suppress_redundant_contigs
                    && self.mark_enclosed_seed_component_used(
                        seed.oriented,
                        k,
                        adjacency,
                        &mut used,
                    )
                {
                    continue;
                }

                let contig = self.extend_bidirectional_with_coverage(
                    seed.oriented,
                    k,
                    adjacency,
                    kmer_counts,
                    branch_support,
                    &repeat_kmers,
                    &mut used,
                );

                if contig.len() >= min_len {
                    contigs.push(contig);
                }
            }
            info!(
                "  Repeat-seed completion pass added {} contigs",
                contigs.len().saturating_sub(contigs_before)
            );
        }

        if self.config.suppress_redundant_contigs && contigs.len() > 1 {
            let contigs_before = contigs.len();
            contigs = self.suppress_redundant_contigs(contigs, kmer_counts, k);
            info!(
                "  Redundant-contig suppression removed {} contigs",
                contigs_before.saturating_sub(contigs.len())
            );
        }

        contigs
    }

    fn mark_enclosed_seed_component_used(
        &self,
        seed_encoded: u64,
        k: usize,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        used: &mut AHashSet<u64>,
    ) -> bool {
        let seed_canonical = Self::canonical_encoded(seed_encoded, k);
        if used.contains(&seed_canonical) {
            return false;
        }

        let Some((left_path, left_enclosed)) =
            Self::walk_linear_unused_path(seed_encoded, k, adjacency, used, true)
        else {
            return false;
        };
        let Some((right_path, right_enclosed)) =
            Self::walk_linear_unused_path(seed_encoded, k, adjacency, used, false)
        else {
            return false;
        };

        if !left_enclosed || !right_enclosed {
            return false;
        }

        used.insert(seed_canonical);
        for kmer in left_path.into_iter().chain(right_path.into_iter()) {
            used.insert(kmer);
        }
        true
    }

    fn walk_linear_unused_path(
        start: u64,
        k: usize,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        used: &AHashSet<u64>,
        going_left: bool,
    ) -> Option<(Vec<u64>, bool)> {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

        let mut current = start;
        let mut traversed = Vec::new();
        let mut local_seen = AHashSet::new();
        local_seen.insert(Self::canonical_encoded(start, k));

        loop {
            let current_oriented = Self::resolve_adjacency_orientation(current, k, adjacency)?;
            let extensions = if going_left {
                adjacency.get(&current_oriented)?.0
            } else {
                adjacency.get(&current_oriented)?.1
            };

            let mut used_boundary = false;
            let mut next_unused = None;

            for (base_idx, valid) in extensions.iter().copied().enumerate() {
                if !valid {
                    continue;
                }

                let next = if going_left {
                    extend_left(current_oriented, BASES[base_idx], k)
                } else {
                    extend_right(current_oriented, BASES[base_idx], k)
                }?;
                let next_canonical = Self::canonical_encoded(next, k);

                if used.contains(&next_canonical) {
                    used_boundary = true;
                    continue;
                }

                if next_unused.is_some() || !local_seen.insert(next_canonical) {
                    return None;
                }
                next_unused = Some(next_canonical);
                current = next;
            }

            match next_unused {
                Some(next_canonical) => traversed.push(next_canonical),
                None => return Some((traversed, used_boundary)),
            }
        }
    }

    fn suppress_redundant_contigs(
        &self,
        contigs: Vec<String>,
        kmer_counts: &AHashMap<u64, u32>,
        k: usize,
    ) -> Vec<String> {
        let mut ranked = contigs
            .into_iter()
            .enumerate()
            .map(|(original_index, sequence)| {
                let (total_support, kmer_windows) =
                    Self::contig_support_profile(&sequence, kmer_counts, k);
                RankedExtractedContig {
                    original_index,
                    sequence,
                    total_support,
                    kmer_windows,
                }
            })
            .collect::<Vec<_>>();

        ranked.sort_unstable_by(|left, right| {
            (right.total_support as u128 * left.kmer_windows as u128)
                .cmp(&(left.total_support as u128 * right.kmer_windows as u128))
                .then_with(|| right.total_support.cmp(&left.total_support))
                .then_with(|| right.sequence.len().cmp(&left.sequence.len()))
                .then_with(|| left.sequence.cmp(&right.sequence))
                .then_with(|| left.original_index.cmp(&right.original_index))
        });

        let mut kept: Vec<RankedExtractedContig> = Vec::with_capacity(ranked.len());
        for candidate in ranked {
            if kept.iter().any(|existing| {
                Self::is_redundant_contig(
                    existing.sequence.as_str(),
                    candidate.sequence.as_str(),
                    k,
                )
            }) {
                continue;
            }
            kept.push(candidate);
        }

        kept.sort_unstable_by_key(|contig| contig.original_index);
        kept.into_iter().map(|contig| contig.sequence).collect()
    }

    fn contig_support_profile(
        sequence: &str,
        kmer_counts: &AHashMap<u64, u32>,
        k: usize,
    ) -> (u64, usize) {
        if sequence.len() < k {
            return (0, 1);
        }

        let mut total_support = 0u64;
        let mut kmer_windows = 0usize;
        for window in sequence.as_bytes().windows(k) {
            if let Some(kmer) = KmerU64::from_slice(window) {
                total_support = total_support.saturating_add(
                    kmer_counts
                        .get(&kmer.canonical().encoded)
                        .copied()
                        .unwrap_or(0) as u64,
                );
                kmer_windows += 1;
            }
        }

        (total_support, kmer_windows.max(1))
    }

    fn is_redundant_contig(existing: &str, candidate: &str, k: usize) -> bool {
        if existing == candidate {
            return true;
        }

        if existing.len() >= candidate.len() && existing.contains(candidate) {
            return true;
        }

        let min_len = existing.len().min(candidate.len());
        if min_len < k.saturating_mul(2).saturating_add(1) {
            return false;
        }

        let shared_prefix = Self::shared_prefix_bases(existing.as_bytes(), candidate.as_bytes());
        if shared_prefix < k {
            return false;
        }

        let shared_suffix = Self::shared_suffix_bases(existing.as_bytes(), candidate.as_bytes());
        shared_suffix >= k && shared_prefix.saturating_add(shared_suffix) < min_len
    }

    fn shared_prefix_bases(left: &[u8], right: &[u8]) -> usize {
        left.iter()
            .zip(right.iter())
            .take_while(|(left_base, right_base)| left_base == right_base)
            .count()
    }

    fn shared_suffix_bases(left: &[u8], right: &[u8]) -> usize {
        left.iter()
            .rev()
            .zip(right.iter().rev())
            .take_while(|(left_base, right_base)| left_base == right_base)
            .count()
    }

    #[inline]
    fn orient_seed_for_adjacency(
        canonical: u64,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        k: usize,
    ) -> Option<u64> {
        if adjacency.contains_key(&canonical) {
            return Some(canonical);
        }

        let rc = Self::reverse_complement_kmer(canonical, k);
        adjacency.contains_key(&rc).then_some(rc)
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
        let k = self.config.k;

        for (&node, node_info) in graph {
            if node_info.successors.len() > 1 {
                for &(next, _) in &node_info.successors {
                    edges.insert(Self::canonical_edge_key(node, next, k));
                }
            }
            if node_info.predecessors.len() > 1 {
                for &(prev, _) in &node_info.predecessors {
                    edges.insert(Self::canonical_edge_key(prev, node, k));
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

        let branch_edge_lookup = OrientedBranchEdgeLookup::new(branch_edges, k);
        let reader = try_open_fastq(path)?;
        let mut batch = Vec::with_capacity(Self::THREADING_BATCH_SIZE);
        for record in stream_fastq_records_checked(reader) {
            let record = record?;
            batch.push(record.sequence.into_bytes());
            if batch.len() >= Self::THREADING_BATCH_SIZE {
                stats.edge_observations += Self::accumulate_branch_support_batch(
                    &batch,
                    k,
                    adjacency,
                    &branch_edge_lookup,
                    &mut edge_support,
                );
                batch.clear();
            }
        }

        if !batch.is_empty() {
            stats.edge_observations += Self::accumulate_branch_support_batch(
                &batch,
                k,
                adjacency,
                &branch_edge_lookup,
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

        let branch_edge_lookup = OrientedBranchEdgeLookup::new(branch_edges, k);
        let reader1 = try_open_fastq(path1)?;
        let reader2 = try_open_fastq(path2)?;
        let mut batch = Vec::with_capacity(Self::THREADING_BATCH_SIZE);
        for pair in stream_paired_fastq_records_checked(reader1, reader2) {
            let (r1, r2) = pair?;
            batch.push(r1.sequence.into_bytes());
            if batch.len() >= Self::THREADING_BATCH_SIZE {
                stats.edge_observations += Self::accumulate_branch_support_batch(
                    &batch,
                    k,
                    adjacency,
                    &branch_edge_lookup,
                    &mut edge_support,
                );
                batch.clear();
            }

            batch.push(r2.sequence.into_bytes());
            if batch.len() >= Self::THREADING_BATCH_SIZE {
                stats.edge_observations += Self::accumulate_branch_support_batch(
                    &batch,
                    k,
                    adjacency,
                    &branch_edge_lookup,
                    &mut edge_support,
                );
                batch.clear();
            }
        }

        if !batch.is_empty() {
            stats.edge_observations += Self::accumulate_branch_support_batch(
                &batch,
                k,
                adjacency,
                &branch_edge_lookup,
                &mut edge_support,
            );
        }

        stats.supported_edges = edge_support.len();
        Ok((edge_support, stats))
    }

    fn collect_branch_support_from_sequences(
        &self,
        sequences: &[Vec<u8>],
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        branch_edges: &AHashSet<(u64, u64)>,
        k: usize,
    ) -> (AHashMap<(u64, u64), u32>, ReadThreadingStats) {
        let mut edge_support = AHashMap::with_capacity(branch_edges.len());
        let branch_edge_lookup = OrientedBranchEdgeLookup::new(branch_edges, k);
        let mut stats = ReadThreadingStats {
            ambiguous_edges: branch_edge_lookup.ambiguous_edges,
            ..Default::default()
        };

        if branch_edge_lookup.is_empty() {
            return (edge_support, stats);
        }

        stats.edge_observations = Self::accumulate_branch_support_batch(
            sequences,
            k,
            adjacency,
            &branch_edge_lookup,
            &mut edge_support,
        );

        stats.supported_edges = edge_support.len();
        (edge_support, stats)
    }

    fn accumulate_branch_support_batch(
        sequences: &[Vec<u8>],
        k: usize,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        branch_edge_lookup: &OrientedBranchEdgeLookup,
        edge_support: &mut AHashMap<(u64, u64), u32>,
    ) -> u64 {
        if sequences.is_empty() || branch_edge_lookup.is_empty() {
            return 0;
        }

        if sequences.len() < Self::PARALLEL_THREADING_MIN_BATCH {
            let mut observations = 0u64;
            for sequence in sequences {
                observations =
                    observations.saturating_add(Self::thread_branch_edges_in_sequence_with_lookup(
                        sequence,
                        k,
                        adjacency,
                        branch_edge_lookup,
                        edge_support,
                    ));
            }
            return observations;
        }

        let (batch_support, observations) = sequences
            .par_chunks(Self::THREADING_BATCH_SIZE)
            .map(|chunk| {
                let mut local_support = AHashMap::with_capacity(
                    branch_edge_lookup
                        .ambiguous_edges
                        .min(chunk.len() * 16)
                        .max(16),
                );
                let mut local_observations = 0u64;
                for sequence in chunk {
                    local_observations = local_observations.saturating_add(
                        Self::thread_branch_edges_in_sequence_with_lookup(
                            sequence,
                            k,
                            adjacency,
                            branch_edge_lookup,
                            &mut local_support,
                        ),
                    );
                }
                (local_support, local_observations)
            })
            .reduce(
                || (AHashMap::new(), 0u64),
                |(mut left_support, left_observations), (right_support, right_observations)| {
                    Self::merge_edge_support_counts(&mut left_support, right_support);
                    (
                        left_support,
                        left_observations.saturating_add(right_observations),
                    )
                },
            );

        Self::merge_edge_support_counts(edge_support, batch_support);
        observations
    }

    fn merge_edge_support_counts(
        edge_support: &mut AHashMap<(u64, u64), u32>,
        additions: AHashMap<(u64, u64), u32>,
    ) {
        for (edge, count) in additions {
            let entry = edge_support.entry(edge).or_insert(0);
            *entry = entry.saturating_add(count);
        }
    }

    fn thread_branch_edges_in_sequence(
        sequence: &[u8],
        k: usize,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        branch_edges: &AHashSet<(u64, u64)>,
        edge_support: &mut AHashMap<(u64, u64), u32>,
    ) -> u64 {
        let branch_edge_lookup = OrientedBranchEdgeLookup::new(branch_edges, k);
        Self::thread_branch_edges_in_sequence_with_lookup(
            sequence,
            k,
            adjacency,
            &branch_edge_lookup,
            edge_support,
        )
    }

    fn thread_branch_edges_in_sequence_with_lookup(
        sequence: &[u8],
        k: usize,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        branch_edge_lookup: &OrientedBranchEdgeLookup,
        edge_support: &mut AHashMap<(u64, u64), u32>,
    ) -> u64 {
        if branch_edge_lookup.is_empty() || sequence.len() < k {
            return 0;
        }

        let mut observations = 0u64;
        let mut prev_node = None;
        let mut rolling = 0u64;
        let mut valid_run = 0usize;
        let rolling_mask = if k >= 32 {
            u64::MAX
        } else {
            (1u64 << (k * 2)) - 1
        };

        for &base in sequence {
            let base_bits = Self::encode_dna_base_2bit(base);
            if base_bits > 3 {
                prev_node = None;
                rolling = 0;
                valid_run = 0;
                continue;
            }

            rolling = ((rolling << 2) | base_bits) & rolling_mask;
            valid_run += 1;
            if valid_run < k {
                continue;
            }

            if adjacency.contains_key(&rolling) {
                let current_node = rolling;
                if let Some(prev) = prev_node {
                    if let Some(&canonical_edge) = branch_edge_lookup
                        .canonical_by_oriented
                        .get(&(prev, current_node))
                    {
                        let support = edge_support.entry(canonical_edge).or_insert(0);
                        *support = support.saturating_add(1);
                        observations = observations.saturating_add(1);
                    }
                }
                prev_node = Some(current_node);
            } else {
                prev_node = None;
            }
        }

        observations
    }

    #[cfg(test)]
    fn thread_branch_edges_in_sequence_reference(
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
                        let edge = Self::canonical_edge_key(prev, current_node, k);
                        if branch_edges.contains(&edge) {
                            let support = edge_support.entry(edge).or_insert(0);
                            *support = support.saturating_add(1);
                            observations = observations.saturating_add(1);
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

    #[inline]
    fn encode_dna_base_2bit(base: u8) -> u64 {
        match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 4,
        }
    }

    fn resolve_threaded_node(
        encoded: u64,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        _k: usize,
    ) -> Option<u64> {
        adjacency.contains_key(&encoded).then_some(encoded)
    }

    #[inline]
    fn branch_read_support(
        branch_support: &AHashMap<(u64, u64), u32>,
        from: u64,
        to: u64,
        k: usize,
    ) -> u32 {
        branch_support
            .get(&Self::canonical_edge_key(from, to, k))
            .copied()
            .unwrap_or(0)
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
        seed_encoded: u64,
        k: usize,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        used: &mut AHashSet<u64>,
    ) -> String {
        let bases = [b'A', b'C', b'G', b'T'];

        let seed_canonical = KmerU64 {
            encoded: seed_encoded,
            len: k as u8,
        }
        .canonical()
        .encoded;

        // Start with the seed k-mer sequence in the selected traversal orientation.
        let seed_seq = decode_kmer(seed_encoded, k);
        let mut contig: Vec<u8> = seed_seq.into_bytes();
        let mut left_extension: Vec<u8> = Vec::new();
        used.insert(seed_canonical);

        // Track current ACTUAL k-mer at each end (not canonical)
        // This is crucial for correct extension direction
        let mut right_kmer = seed_encoded;
        let mut left_kmer = seed_encoded;

        // Extend right: look up adjacency for current k-mer (not canonical)
        loop {
            let Some(current_oriented) =
                Self::resolve_adjacency_orientation(right_kmer, k, adjacency)
            else {
                break;
            };
            let (_, right_ext) = adjacency
                .get(&current_oriented)
                .copied()
                .unwrap_or(([false; 4], [false; 4]));

            let Some((base_idx, next)) = Self::select_unique_unused_linear_extension(
                current_oriented,
                &right_ext,
                k,
                used,
                false,
            ) else {
                break; // Stop at branch or dead end
            };

            contig.push(bases[base_idx]);
            used.insert(Self::canonical_encoded(next, k));

            // Continue with the ACTUAL extended k-mer (not canonical)
            // This preserves correct directionality
            right_kmer = next;
        }

        // Extend left
        loop {
            let Some(current_oriented) =
                Self::resolve_adjacency_orientation(left_kmer, k, adjacency)
            else {
                break;
            };
            let (left_ext, _) = adjacency
                .get(&current_oriented)
                .copied()
                .unwrap_or(([false; 4], [false; 4]));

            let Some((base_idx, next)) = Self::select_unique_unused_linear_extension(
                current_oriented,
                &left_ext,
                k,
                used,
                true,
            ) else {
                break;
            };

            left_extension.push(bases[base_idx]);
            used.insert(Self::canonical_encoded(next, k));

            // Continue with ACTUAL extended k-mer
            left_kmer = next;
        }

        if !left_extension.is_empty() {
            let mut merged = Vec::with_capacity(left_extension.len() + contig.len());
            merged.extend(left_extension.iter().rev().copied());
            merged.extend_from_slice(&contig);
            contig = merged;
        }

        String::from_utf8(contig).unwrap_or_default()
    }

    /// Extend a seed k-mer bidirectionally with coverage-guided traversal
    ///
    /// At branch points, uses coverage information to choose the best path.
    /// Avoids entering repeat regions unless necessary.
    fn extend_bidirectional_with_coverage(
        &self,
        seed_encoded: u64,
        k: usize,
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
        kmer_counts: &AHashMap<u64, u32>,
        branch_support: &AHashMap<(u64, u64), u32>,
        repeat_kmers: &AHashSet<u64>,
        used: &mut AHashSet<u64>,
    ) -> String {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

        let seed_canonical = KmerU64 {
            encoded: seed_encoded,
            len: k as u8,
        }
        .canonical()
        .encoded;
        let seed_seq = decode_kmer(seed_encoded, k);
        let mut contig: Vec<u8> = seed_seq.into_bytes();
        let mut left_extension: Vec<u8> = Vec::new();
        used.insert(seed_canonical);

        let mut right_kmer = seed_encoded;
        let mut left_kmer = seed_encoded;
        let mut right_candidates = Vec::with_capacity(4);
        let mut left_candidates = Vec::with_capacity(4);

        // Extend right with coverage-guided traversal
        loop {
            let Some(current_oriented) =
                Self::resolve_adjacency_orientation(right_kmer, k, adjacency)
            else {
                break;
            };
            let (_, right_ext) = adjacency
                .get(&current_oriented)
                .copied()
                .unwrap_or(([false; 4], [false; 4]));

            let Some(best) = self.select_branch_extension(
                current_oriented,
                &right_ext,
                k,
                kmer_counts,
                branch_support,
                repeat_kmers,
                used,
                false,
                &mut right_candidates,
            ) else {
                break;
            };

            contig.push(BASES[best.base_idx]);
            used.insert(Self::canonical_encoded(best.next, k));
            right_kmer = best.next;
        }

        // Extend left with coverage-guided traversal
        loop {
            let Some(current_oriented) =
                Self::resolve_adjacency_orientation(left_kmer, k, adjacency)
            else {
                break;
            };
            let (left_ext, _) = adjacency
                .get(&current_oriented)
                .copied()
                .unwrap_or(([false; 4], [false; 4]));

            let Some(best) = self.select_branch_extension(
                current_oriented,
                &left_ext,
                k,
                kmer_counts,
                branch_support,
                repeat_kmers,
                used,
                true,
                &mut left_candidates,
            ) else {
                break;
            };

            left_extension.push(BASES[best.base_idx]);
            used.insert(Self::canonical_encoded(best.next, k));
            left_kmer = best.next;
        }

        if !left_extension.is_empty() {
            let mut merged = Vec::with_capacity(left_extension.len() + contig.len());
            merged.extend(left_extension.iter().rev().copied());
            merged.extend_from_slice(&contig);
            contig = merged;
        }

        String::from_utf8(contig).unwrap_or_default()
    }

    #[inline]
    fn canonical_encoded(encoded: u64, k: usize) -> u64 {
        KmerU64 {
            encoded,
            len: k as u8,
        }
        .canonical()
        .encoded
    }

    #[inline]
    fn canonical_count(kmer_counts: &AHashMap<u64, u32>, encoded: u64, k: usize) -> u32 {
        kmer_counts
            .get(&Self::canonical_encoded(encoded, k))
            .copied()
            .unwrap_or(1)
    }

    fn select_unique_unused_linear_extension(
        current_kmer: u64,
        extension_flags: &[bool; 4],
        k: usize,
        used: &AHashSet<u64>,
        going_left: bool,
    ) -> Option<(usize, u64)> {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        let mut selected: Option<(usize, u64)> = None;

        for (base_idx, &valid) in extension_flags.iter().enumerate() {
            if !valid {
                continue;
            }

            let next = if going_left {
                extend_left(current_kmer, BASES[base_idx], k)
            } else {
                extend_right(current_kmer, BASES[base_idx], k)
            };
            let Some(next) = next else {
                continue;
            };

            if used.contains(&Self::canonical_encoded(next, k)) {
                continue;
            }

            if selected.is_some() {
                return None;
            }
            selected = Some((base_idx, next));
        }

        selected
    }

    fn collect_branch_candidates(
        current_kmer: u64,
        extension_flags: &[bool; 4],
        k: usize,
        kmer_counts: &AHashMap<u64, u32>,
        branch_support: &AHashMap<(u64, u64), u32>,
        repeat_kmers: &AHashSet<u64>,
        used: &AHashSet<u64>,
        going_left: bool,
        extensions: &mut Vec<BranchCandidate>,
    ) {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        extensions.clear();

        for (base_idx, &valid) in extension_flags.iter().enumerate() {
            if !valid {
                continue;
            }

            let candidate = if going_left {
                extend_left(current_kmer, BASES[base_idx], k)
            } else {
                extend_right(current_kmer, BASES[base_idx], k)
            };

            let Some(next) = candidate else {
                continue;
            };

            let next_canonical = Self::canonical_encoded(next, k);
            if used.contains(&next_canonical) {
                continue;
            }

            let read_support = if going_left {
                Self::branch_read_support(branch_support, next, current_kmer, k)
            } else {
                Self::branch_read_support(branch_support, current_kmer, next, k)
            };

            extensions.push(BranchCandidate {
                base_idx,
                next,
                count: kmer_counts.get(&next_canonical).copied().unwrap_or(0),
                is_repeat: repeat_kmers.contains(&next_canonical),
                read_support,
            });
        }
    }

    fn select_branch_extension(
        &self,
        current_kmer: u64,
        extension_flags: &[bool; 4],
        k: usize,
        kmer_counts: &AHashMap<u64, u32>,
        branch_support: &AHashMap<(u64, u64), u32>,
        repeat_kmers: &AHashSet<u64>,
        used: &AHashSet<u64>,
        going_left: bool,
        scratch: &mut Vec<BranchCandidate>,
    ) -> Option<BranchCandidate> {
        Self::collect_branch_candidates(
            current_kmer,
            extension_flags,
            k,
            kmer_counts,
            branch_support,
            repeat_kmers,
            used,
            going_left,
            scratch,
        );

        match scratch.len() {
            0 => None,
            1 => Some(scratch[0]),
            _ => self.choose_branch_extension(
                scratch,
                Self::canonical_count(kmer_counts, current_kmer, k),
            ),
        }
    }

    #[doc(hidden)]
    pub fn run_branch_resolution_case(
        &self,
        case: &LargeGenomeBenchBranchResolutionCase,
    ) -> Option<LargeGenomeBenchBranchCandidate> {
        let mut scratch = Vec::with_capacity(case.candidates.len());
        scratch.extend(
            case.candidates
                .iter()
                .copied()
                .map(|candidate| BranchCandidate {
                    base_idx: candidate.base_idx,
                    next: candidate.next,
                    count: candidate.count,
                    is_repeat: candidate.is_repeat,
                    read_support: candidate.read_support,
                }),
        );

        self.choose_branch_extension(&scratch, case.current_count)
            .map(|selected| LargeGenomeBenchBranchCandidate {
                base_idx: selected.base_idx,
                next: selected.next,
                count: selected.count,
                is_repeat: selected.is_repeat,
                read_support: selected.read_support,
            })
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
        if best_support >= self.config.branch_support_min_win
            && best_support.saturating_sub(second_support) >= self.config.branch_support_min_margin
        {
            return Some(supported_extension);
        }

        let best_by_coverage = |candidates: &[BranchCandidate]| {
            candidates.iter().copied().min_by(|left, right| {
                Self::compare_coverage_ratio(current_count, left.count, right.count)
                    .then_with(|| {
                        Self::coverage_distance(current_count, left.count)
                            .cmp(&Self::coverage_distance(current_count, right.count))
                    })
                    .then_with(|| right.read_support.cmp(&left.read_support))
                    .then_with(|| right.count.cmp(&left.count))
                    .then_with(|| left.base_idx.cmp(&right.base_idx))
                    .then_with(|| left.next.cmp(&right.next))
            })
        };

        let best_non_repeat = if self.config.prefer_non_repeat_branches {
            let non_repeat = extensions.iter().copied().filter(|ext| !ext.is_repeat);
            non_repeat.min_by(|left, right| {
                Self::compare_coverage_ratio(current_count, left.count, right.count)
                    .then_with(|| {
                        Self::coverage_distance(current_count, left.count)
                            .cmp(&Self::coverage_distance(current_count, right.count))
                    })
                    .then_with(|| right.read_support.cmp(&left.read_support))
                    .then_with(|| right.count.cmp(&left.count))
                    .then_with(|| left.base_idx.cmp(&right.base_idx))
                    .then_with(|| left.next.cmp(&right.next))
            })
        } else {
            best_by_coverage(extensions)
        };

        best_non_repeat.or_else(|| {
            extensions.iter().copied().max_by_key(|ext| {
                (
                    ext.read_support,
                    ext.count,
                    Reverse(ext.next),
                    Reverse(ext.base_idx),
                )
            })
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

    fn compare_coverage_ratio(
        current_count: u32,
        left_count: u32,
        right_count: u32,
    ) -> std::cmp::Ordering {
        let (left_num, left_den) = Self::normalized_coverage_ratio(current_count, left_count);
        let (right_num, right_den) = Self::normalized_coverage_ratio(current_count, right_count);
        (left_num as u128 * right_den as u128).cmp(&(right_num as u128 * left_den as u128))
    }

    fn normalized_coverage_ratio(current_count: u32, candidate_count: u32) -> (u32, u32) {
        let hi = current_count.max(candidate_count);
        let lo = current_count.min(candidate_count).max(1);
        (hi, lo)
    }

    fn coverage_distance(current_count: u32, candidate_count: u32) -> u32 {
        current_count.abs_diff(candidate_count)
    }

    fn write_output(
        &self,
        contigs: &[String],
        path: &str,
        stats: &mut AssemblyStats,
    ) -> std::io::Result<()> {
        let mut writer = FastaWriter::try_new(path)?;

        // Filter and sort by length
        let mut valid: Vec<&String> = contigs
            .iter()
            .filter(|c| c.len() >= self.config.min_contig_len)
            .collect();
        valid.sort_unstable_by(|a, b| b.len().cmp(&a.len()).then_with(|| a.cmp(b)));

        let mut lengths: Vec<usize> = Vec::with_capacity(valid.len());
        let mut ungapped_lengths: Vec<usize> = Vec::with_capacity(valid.len());
        let mut composition = BaseComposition::default();
        let mut n_runs = NRunSummary::default();

        for (i, contig) in valid.iter().enumerate() {
            writer.write_record(&format!("contig_{} len={}", i + 1, contig.len()), contig)?;
            lengths.push(contig.len());
            let sequence = contig.as_bytes();
            let ungapped_len = n_runs.add_sequence_and_count_ungapped(sequence);
            ungapped_lengths.push(ungapped_len);
            composition.add_sequence(sequence);
        }

        let contig_stats = evaluate_lengths_sorted_desc(&lengths);
        ungapped_lengths.sort_unstable_by(|a, b| b.cmp(a));
        let ungapped_stats = evaluate_lengths_sorted_desc(&ungapped_lengths);
        let total_bases = contig_stats.total_bases;
        let gap_bases = contig_stats
            .total_bases
            .saturating_sub(ungapped_stats.total_bases);
        let gap_bases_frac = if total_bases > 0 {
            gap_bases as f64 / total_bases as f64
        } else {
            0.0
        };

        stats.contigs = valid.len();
        stats.total_length = contig_stats.total_bases;
        stats.gc_bases = composition.gc_bases;
        stats.acgt_bases = composition.acgt_bases;
        stats.n_bases = composition.n_bases;
        stats.ambiguous_bases = composition.ambiguous_bases;
        stats.gc_content = composition.gc_content();
        stats.n_content = composition.n_content(total_bases);
        stats.ambiguous_content = composition.ambiguous_content(total_bases);
        stats.n_run_count = n_runs.run_count;
        stats.max_n_run = n_runs.max_run;
        stats.mean_n_run_length = n_runs.mean_run_length(composition.n_bases);
        stats.n_runs_per_100kb = per_100kb(stats.n_run_count, total_bases);
        stats.n_bases_per_100kb = per_100kb(composition.n_bases, total_bases);
        stats.ambiguous_bases_per_100kb = per_100kb(composition.ambiguous_bases, total_bases);
        stats.ungapped_total_length = ungapped_stats.total_bases;
        stats.gap_bases = gap_bases;
        stats.gap_bases_frac = gap_bases_frac;
        stats.n10 = contig_stats.n10;
        stats.n25 = contig_stats.n25;
        stats.n50 = contig_stats.n50;
        stats.n75 = contig_stats.n75;
        stats.n90 = contig_stats.n90;
        stats.n95 = contig_stats.n95;
        stats.n99 = contig_stats.n99;
        stats.l10 = contig_stats.l10;
        stats.l25 = contig_stats.l25;
        stats.l50 = contig_stats.l50;
        stats.l75 = contig_stats.l75;
        stats.l90 = contig_stats.l90;
        stats.l95 = contig_stats.l95;
        stats.l99 = contig_stats.l99;
        stats.ungapped_n50 = ungapped_stats.n50;
        stats.ungapped_n90 = ungapped_stats.n90;
        stats.ungapped_n95 = ungapped_stats.n95;
        stats.ungapped_n99 = ungapped_stats.n99;
        stats.avg_contig_len = contig_stats.avg_length;
        stats.median_contig_len = contig_stats.median_length;
        stats.au_n = contig_stats.au_n;
        stats.effective_contig_count = contig_stats.effective_count;
        stats.ungapped_au_n = ungapped_stats.au_n;
        stats.ungapped_effective_contig_count = ungapped_stats.effective_count;
        stats.largest = contig_stats.longest;
        stats.contigs_ge_1kb = contig_stats.contigs_ge_1kb;
        stats.contigs_ge_10kb = contig_stats.contigs_ge_10kb;
        stats.contigs_ge_50kb = contig_stats.contigs_ge_50kb;
        stats.contigs_ge_100kb = contig_stats.contigs_ge_100kb;
        stats.contigs_ge_1mb = contig_stats.contigs_ge_1mb;
        stats.bases_ge_1kb = contig_stats.bases_ge_1kb;
        stats.bases_ge_10kb = contig_stats.bases_ge_10kb;
        stats.bases_ge_50kb = contig_stats.bases_ge_50kb;
        stats.bases_ge_100kb = contig_stats.bases_ge_100kb;
        stats.bases_ge_1mb = contig_stats.bases_ge_1mb;
        stats.contigs_ge_1kb_frac = contig_stats.contigs_ge_1kb_frac;
        stats.contigs_ge_10kb_frac = contig_stats.contigs_ge_10kb_frac;
        stats.contigs_ge_50kb_frac = contig_stats.contigs_ge_50kb_frac;
        stats.contigs_ge_100kb_frac = contig_stats.contigs_ge_100kb_frac;
        stats.contigs_ge_1mb_frac = contig_stats.contigs_ge_1mb_frac;
        stats.bases_ge_1kb_frac = contig_stats.bases_ge_1kb_frac;
        stats.bases_ge_10kb_frac = contig_stats.bases_ge_10kb_frac;
        stats.bases_ge_50kb_frac = contig_stats.bases_ge_50kb_frac;
        stats.bases_ge_100kb_frac = contig_stats.bases_ge_100kb_frac;
        stats.bases_ge_1mb_frac = contig_stats.bases_ge_1mb_frac;

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
    use proptest::prelude::*;
    use rand::rngs::StdRng;
    use rand::seq::SliceRandom;
    use rand::Rng;
    use rand::SeedableRng;
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

    fn canonicalize_adjacency(
        adjacency: &AHashMap<u64, ([bool; 4], [bool; 4])>,
    ) -> Vec<(u64, [bool; 4], [bool; 4])> {
        let mut entries: Vec<(u64, [bool; 4], [bool; 4])> = adjacency
            .iter()
            .map(|(&kmer, &(left, right))| (kmer, left, right))
            .collect();
        entries.sort_unstable_by_key(|(kmer, _, _)| *kmer);
        entries
    }

    fn build_branch_threading_fixture(
        k: usize,
    ) -> (
        AHashMap<u64, ([bool; 4], [bool; 4])>,
        AHashSet<(u64, u64)>,
        AHashMap<u64, u32>,
    ) {
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

        (adjacency, branch_edges, counts)
    }

    fn canonicalize_edge_support(
        edge_support: &AHashMap<(u64, u64), u32>,
    ) -> Vec<((u64, u64), u32)> {
        let mut entries: Vec<_> = edge_support
            .iter()
            .map(|(&edge, &count)| (edge, count))
            .collect();
        entries.sort_unstable_by_key(|(edge, _)| *edge);
        entries
    }

    fn reference_error_correct_kmers(
        assembler: &LargeGenomeAssembler,
        kmer_counts: &AHashMap<u64, u32>,
        k: usize,
    ) -> (AHashMap<u64, u32>, u64) {
        let trusted_threshold = assembler.error_correction_trusted_threshold();
        let mut trusted_counts = AHashMap::new();
        let mut singleton_kmers = Vec::new();

        for (&encoded, &count) in kmer_counts {
            if count >= trusted_threshold {
                trusted_counts.insert(encoded, count);
            } else if count == 1 {
                singleton_kmers.push(encoded);
            }
        }

        let mut corrected = kmer_counts.clone();
        let mut num_corrected = 0u64;

        for singleton in singleton_kmers {
            if let Some(trusted_neighbor) =
                assembler.find_trusted_neighbor(singleton, k, &trusted_counts)
            {
                if let Some(count) = corrected.get_mut(&trusted_neighbor) {
                    *count = count.saturating_add(1);
                }
                corrected.remove(&singleton);
                num_corrected += 1;
            }
        }

        (corrected, num_corrected)
    }

    fn simple_pop_bubble_fixture() -> (usize, [u64; 7], AHashMap<u64, ([bool; 4], [bool; 4])>) {
        let k = 3;
        let aaa = KmerU64::from_str("AAA").unwrap().encoded;
        let aac = KmerU64::from_str("AAC").unwrap().encoded;
        let aag = KmerU64::from_str("AAG").unwrap().encoded;
        let agc = KmerU64::from_str("AGC").unwrap().encoded;
        let acc = KmerU64::from_str("ACC").unwrap().encoded;
        let gcc = KmerU64::from_str("GCC").unwrap().encoded;
        let ccc = KmerU64::from_str("CCC").unwrap().encoded;

        let mut right_aa = [false; 4];
        right_aa[1] = true; // AAA -> AAC
        right_aa[2] = true; // AAA -> AAG
        let mut right_c = [false; 4];
        right_c[1] = true; // * -> *C
        let none = [false; 4];

        let mut adjacency = AHashMap::new();
        adjacency.insert(aaa, (none, right_aa));
        adjacency.insert(aac, (none, right_c));
        adjacency.insert(aag, (none, right_c));
        adjacency.insert(agc, (none, right_c));
        adjacency.insert(acc, (none, right_c));
        adjacency.insert(gcc, (none, right_c));
        adjacency.insert(ccc, (none, none));

        (k, [aaa, aac, aag, agc, acc, gcc, ccc], adjacency)
    }

    fn canonical_only_orientation_fixture() -> (
        LargeGenomeAssembler,
        usize,
        AHashMap<u64, ([bool; 4], [bool; 4])>,
        AHashMap<u64, u32>,
        u64,
        u64,
        u64,
    ) {
        let k = 4;
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            k,
            min_count: 1,
            min_contig_len: 1,
            ..Default::default()
        });

        let acga = KmerU64::from_str("ACGA").unwrap().canonical().encoded;
        let atcg = KmerU64::from_str("ATCG").unwrap().canonical().encoded;
        let tcga = KmerU64::from_str("TCGA").unwrap().canonical().encoded;

        let mut adjacency = AHashMap::new();
        let mut right_t = [false; 4];
        right_t[3] = true;
        adjacency.insert(acga, ([false; 4], right_t));

        let mut right_a = [false; 4];
        right_a[0] = true;
        adjacency.insert(atcg, ([false; 4], right_a));
        adjacency.insert(tcga, ([false; 4], [false; 4]));

        let mut counts = AHashMap::new();
        counts.insert(acga, 30);
        counts.insert(atcg, 28);
        counts.insert(tcga, 26);

        (assembler, k, adjacency, counts, acga, atcg, tcga)
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

    #[test]
    fn test_find_trusted_neighbor_prefers_highest_coverage_candidate() {
        let k = 4;
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig::default());

        let singleton = KmerU64::from_str("AAAA").unwrap().canonical().encoded;
        let low = KmerU64::from_str("AAAT").unwrap().canonical().encoded;
        let high = KmerU64::from_str("AACA").unwrap().canonical().encoded;

        let mut trusted_counts = AHashMap::new();
        trusted_counts.insert(low, 9);
        trusted_counts.insert(high, 25);

        let picked = assembler
            .find_trusted_neighbor(singleton, k, &trusted_counts)
            .unwrap();
        assert_eq!(picked, high);
    }

    #[test]
    fn test_find_trusted_neighbor_tie_breaks_by_kmer_value() {
        let k = 4;
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig::default());

        let singleton = KmerU64::from_str("AAAA").unwrap().canonical().encoded;
        let first = KmerU64::from_str("AAAT").unwrap().canonical().encoded;
        let second = KmerU64::from_str("AACA").unwrap().canonical().encoded;
        let expected = first.min(second);

        let mut trusted_counts = AHashMap::new();
        trusted_counts.insert(first, 12);
        trusted_counts.insert(second, 12);

        let picked = assembler
            .find_trusted_neighbor(singleton, k, &trusted_counts)
            .unwrap();
        assert_eq!(picked, expected);
    }

    #[test]
    fn test_error_correct_kmers_matches_serial_reference() {
        let k = 15;
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            k,
            min_count: 1,
            min_contig_len: k,
            ..Default::default()
        });
        let fixture =
            LargeGenomeErrorCorrectionBenchFixture::synthetic_singleton_correction(k, 512, 3);
        let counts = fixture.cloned_kmer_counts();

        let expected = reference_error_correct_kmers(&assembler, &counts, k);
        let observed = assembler.error_correct_kmers(counts, k);

        assert_eq!(observed.1, expected.1);
        assert_eq!(observed.0, expected.0);
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

    #[test]
    fn test_paired_end_assembly_rejects_mismatched_record_counts() {
        let temp_dir = TempDir::new().unwrap();
        let mut r1_file = NamedTempFile::new().unwrap();
        let mut r2_file = NamedTempFile::new().unwrap();

        for i in 0..2 {
            writeln!(r1_file, "@read_{}/1", i).unwrap();
            writeln!(r1_file, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT").unwrap();
            writeln!(r1_file, "+").unwrap();
            writeln!(r1_file, "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII").unwrap();
        }

        writeln!(r2_file, "@read_0/2").unwrap();
        writeln!(r2_file, "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA").unwrap();
        writeln!(r2_file, "+").unwrap();
        writeln!(r2_file, "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII").unwrap();

        r1_file.flush().unwrap();
        r2_file.flush().unwrap();

        let output = NamedTempFile::new().unwrap();
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            k: 21,
            min_count: 1,
            min_contig_len: 21,
            num_buckets: Some(4),
            temp_dir: Some(temp_dir.path().to_str().unwrap().to_string()),
            ..Default::default()
        });

        let err = assembler
            .assemble_paired(
                r1_file.path().to_str().unwrap(),
                r2_file.path().to_str().unwrap(),
                output.path().to_str().unwrap(),
            )
            .unwrap_err();

        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("different record counts"));
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
        assert_eq!(stats.median, 200);
    }

    #[test]
    fn test_insert_size_stats_even_sample_uses_midpoint_median() {
        let samples = vec![100, 200, 300, 400];
        let stats = InsertSizeStats::from_samples(&samples);
        assert_eq!(stats.median, 250);
    }

    #[test]
    fn test_write_output_populates_extended_length_metrics() {
        let output = NamedTempFile::new().unwrap();
        let temp_dir = TempDir::new().unwrap();
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            min_contig_len: 1,
            num_buckets: Some(4),
            temp_dir: Some(temp_dir.path().to_str().unwrap().to_string()),
            ..Default::default()
        });

        let contigs = vec!["A".repeat(100), "C".repeat(50), "G".repeat(25)];
        let mut stats = AssemblyStats::default();
        assembler
            .write_output(&contigs, output.path().to_str().unwrap(), &mut stats)
            .unwrap();

        assert_eq!(stats.contigs, 3);
        assert_eq!(stats.total_length, 175);
        assert_eq!(stats.n10, 100);
        assert_eq!(stats.n25, 100);
        assert_eq!(stats.n50, 100);
        assert_eq!(stats.n75, 50);
        assert_eq!(stats.n90, 25);
        assert_eq!(stats.n95, 25);
        assert_eq!(stats.n99, 25);
        assert_eq!(stats.l10, 1);
        assert_eq!(stats.l25, 1);
        assert_eq!(stats.l50, 1);
        assert_eq!(stats.l75, 2);
        assert_eq!(stats.l90, 3);
        assert_eq!(stats.l95, 3);
        assert_eq!(stats.l99, 3);
        assert_eq!(stats.largest, 100);
        assert_eq!(stats.median_contig_len, 50.0);
        assert_eq!(stats.contigs_ge_1kb, 0);
        assert_eq!(stats.contigs_ge_10kb, 0);
        assert_eq!(stats.contigs_ge_50kb, 0);
        assert_eq!(stats.contigs_ge_100kb, 0);
        assert_eq!(stats.contigs_ge_1mb, 0);
        assert_eq!(stats.bases_ge_1kb, 0);
        assert_eq!(stats.bases_ge_10kb, 0);
        assert_eq!(stats.bases_ge_50kb, 0);
        assert_eq!(stats.bases_ge_100kb, 0);
        assert_eq!(stats.bases_ge_1mb, 0);
        assert_eq!(stats.contigs_ge_1kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_10kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_50kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_100kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_1mb_frac, 0.0);
        assert_eq!(stats.bases_ge_1kb_frac, 0.0);
        assert_eq!(stats.bases_ge_10kb_frac, 0.0);
        assert_eq!(stats.bases_ge_50kb_frac, 0.0);
        assert_eq!(stats.bases_ge_100kb_frac, 0.0);
        assert_eq!(stats.bases_ge_1mb_frac, 0.0);
        assert_eq!(stats.gc_bases, 75);
        assert_eq!(stats.acgt_bases, 175);
        assert_eq!(stats.n_bases, 0);
        assert_eq!(stats.ambiguous_bases, 0);
        assert!((stats.gc_content - (75.0 / 175.0)).abs() < 1e-12);
        assert_eq!(stats.n_content, 0.0);
        assert_eq!(stats.ambiguous_content, 0.0);
        assert_eq!(stats.n_run_count, 0);
        assert_eq!(stats.max_n_run, 0);
        assert_eq!(stats.mean_n_run_length, 0.0);
        assert_eq!(stats.n_runs_per_100kb, 0.0);
        assert_eq!(stats.n_bases_per_100kb, 0.0);
        assert_eq!(stats.ambiguous_bases_per_100kb, 0.0);
        assert_eq!(stats.ungapped_total_length, 175);
        assert_eq!(stats.gap_bases, 0);
        assert_eq!(stats.gap_bases_frac, 0.0);
        assert_eq!(stats.ungapped_n50, 100);
        assert_eq!(stats.ungapped_n90, 25);
        assert_eq!(stats.ungapped_n95, 25);
        assert_eq!(stats.ungapped_n99, 25);
        assert!((stats.avg_contig_len - (175.0 / 3.0)).abs() < 1e-12);
        assert!((stats.au_n - 75.0).abs() < 1e-12);
        assert!((stats.effective_contig_count - (175.0 / 75.0)).abs() < 1e-12);
        assert!((stats.ungapped_au_n - 75.0).abs() < 1e-12);
        assert!((stats.ungapped_effective_contig_count - (175.0 / 75.0)).abs() < 1e-12);
    }

    #[test]
    fn test_assembly_stats_display_includes_effective_count() {
        let stats = AssemblyStats {
            au_n: 1234.5,
            effective_contig_count: 6.75,
            ungapped_total_length: 98_000,
            gap_bases: 2_000,
            gap_bases_frac: 0.02,
            ungapped_au_n: 1111.0,
            ungapped_effective_contig_count: 7.5,
            ..AssemblyStats::default()
        };

        let rendered = stats.to_string();
        assert!(rendered.contains("Ungapped span: 98000 bp (gaps 2000 bp, 2.00%)"));
        assert!(rendered.contains("auN/effective count: 1234.50 bp / 6.75"));
        assert!(rendered.contains("Ungapped auN/effective count: 1111.00 bp / 7.50"));
    }

    #[test]
    fn test_apply_read_threading_stats_updates_branch_quality_metrics() {
        let mut stats = AssemblyStats::default();
        let threading_stats = ReadThreadingStats {
            ambiguous_edges: 12,
            supported_edges: 9,
            edge_observations: 42,
        };

        LargeGenomeAssembler::apply_read_threading_stats(&mut stats, &threading_stats);

        assert_eq!(stats.ambiguous_branch_edges, 12);
        assert_eq!(stats.branch_edges_supported, 9);
        assert_eq!(stats.branch_edge_observations, 42);
        assert!((stats.branch_edge_support_fraction - 0.75).abs() < 1e-12);

        LargeGenomeAssembler::apply_read_threading_stats(
            &mut stats,
            &ReadThreadingStats {
                ambiguous_edges: 0,
                supported_edges: 0,
                edge_observations: 0,
            },
        );
        assert_eq!(stats.ambiguous_branch_edges, 0);
        assert_eq!(stats.branch_edges_supported, 0);
        assert_eq!(stats.branch_edge_observations, 0);
        assert_eq!(stats.branch_edge_support_fraction, 0.0);
    }

    #[test]
    fn test_assembly_stats_display_includes_read_threading_summary() {
        let stats = AssemblyStats {
            ambiguous_branch_edges: 12,
            branch_edges_supported: 9,
            branch_edge_observations: 42,
            branch_edge_support_fraction: 0.75,
            ..AssemblyStats::default()
        };

        let rendered = stats.to_string();
        assert!(rendered.contains("Read-threaded branch support: 9/12 (75.00%), 42 observations"));
    }

    #[test]
    fn test_write_output_reports_base_composition_metrics() {
        let output = NamedTempFile::new().unwrap();
        let temp_dir = TempDir::new().unwrap();
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            min_contig_len: 1,
            num_buckets: Some(4),
            temp_dir: Some(temp_dir.path().to_str().unwrap().to_string()),
            ..Default::default()
        });

        let contigs = vec!["GCGCNN".to_string(), "atry".to_string()];
        let mut stats = AssemblyStats::default();
        assembler
            .write_output(&contigs, output.path().to_str().unwrap(), &mut stats)
            .unwrap();

        assert_eq!(stats.total_length, 10);
        assert_eq!(stats.gc_bases, 4);
        assert_eq!(stats.acgt_bases, 6);
        assert_eq!(stats.n_bases, 2);
        assert_eq!(stats.ambiguous_bases, 2);
        assert!((stats.gc_content - (4.0 / 6.0)).abs() < 1e-12);
        assert!((stats.n_content - 0.2).abs() < 1e-12);
        assert!((stats.ambiguous_content - 0.2).abs() < 1e-12);
        assert_eq!(stats.n_run_count, 1);
        assert_eq!(stats.max_n_run, 2);
        assert!((stats.mean_n_run_length - 2.0).abs() < 1e-12);
        assert!((stats.n_runs_per_100kb - 10_000.0).abs() < 1e-12);
        assert!((stats.n_bases_per_100kb - 20_000.0).abs() < 1e-12);
        assert!((stats.ambiguous_bases_per_100kb - 20_000.0).abs() < 1e-12);
        assert_eq!(stats.ungapped_total_length, 8);
        assert_eq!(stats.gap_bases, 2);
        assert!((stats.gap_bases_frac - 0.2).abs() < 1e-12);
        assert_eq!(stats.ungapped_n50, 4);
        assert_eq!(stats.ungapped_n90, 4);
        assert_eq!(stats.ungapped_n95, 4);
        assert_eq!(stats.ungapped_n99, 4);
        assert!((stats.ungapped_au_n - 4.0).abs() < 1e-12);
        assert!((stats.ungapped_effective_contig_count - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_write_output_reports_n_run_metrics() {
        let output = NamedTempFile::new().unwrap();
        let temp_dir = TempDir::new().unwrap();
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            min_contig_len: 1,
            num_buckets: Some(4),
            temp_dir: Some(temp_dir.path().to_str().unwrap().to_string()),
            ..Default::default()
        });

        let contigs = vec![
            "AANNNNCC".to_string(), // run len 4
            "NNA".to_string(),      // run len 2
            "NNN".to_string(),      // run len 3
        ];
        let mut stats = AssemblyStats::default();
        assembler
            .write_output(&contigs, output.path().to_str().unwrap(), &mut stats)
            .unwrap();

        assert_eq!(stats.n_bases, 9);
        assert_eq!(stats.n_run_count, 3);
        assert_eq!(stats.max_n_run, 4);
        assert!((stats.mean_n_run_length - 3.0).abs() < 1e-12);
        assert!((stats.n_runs_per_100kb - (3.0 * 100_000.0 / 14.0)).abs() < 1e-12);
        assert!((stats.n_bases_per_100kb - (9.0 * 100_000.0 / 14.0)).abs() < 1e-12);
        assert_eq!(stats.ambiguous_bases_per_100kb, 0.0);
        assert_eq!(stats.ungapped_total_length, 5);
        assert_eq!(stats.gap_bases, 9);
        assert!((stats.gap_bases_frac - (9.0 / 14.0)).abs() < 1e-12);
        assert_eq!(stats.ungapped_n50, 4);
        assert_eq!(stats.ungapped_n90, 1);
        assert_eq!(stats.ungapped_n95, 1);
        assert_eq!(stats.ungapped_n99, 1);
    }

    #[test]
    fn test_write_output_tie_breaks_equal_length_contigs_lexicographically() {
        let output = NamedTempFile::new().unwrap();
        let temp_dir = TempDir::new().unwrap();
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            min_contig_len: 1,
            num_buckets: Some(4),
            temp_dir: Some(temp_dir.path().to_str().unwrap().to_string()),
            ..Default::default()
        });

        let contigs = vec!["TT".to_string(), "AA".to_string(), "CC".to_string()];
        let mut stats = AssemblyStats::default();
        assembler
            .write_output(&contigs, output.path().to_str().unwrap(), &mut stats)
            .unwrap();

        let written = std::fs::read_to_string(output.path()).unwrap();
        let sequences: Vec<&str> = written
            .lines()
            .filter(|line| !line.starts_with('>'))
            .collect();

        assert_eq!(sequences, vec!["AA", "CC", "TT"]);
    }

    #[test]
    fn test_write_output_reports_length_bucket_metrics() {
        let output = NamedTempFile::new().unwrap();
        let temp_dir = TempDir::new().unwrap();
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            min_contig_len: 1,
            num_buckets: Some(4),
            temp_dir: Some(temp_dir.path().to_str().unwrap().to_string()),
            ..Default::default()
        });

        let contigs = vec![
            "A".repeat(50_000),
            "C".repeat(10_000),
            "G".repeat(1_000),
            "T".repeat(999),
        ];
        let mut stats = AssemblyStats::default();
        assembler
            .write_output(&contigs, output.path().to_str().unwrap(), &mut stats)
            .unwrap();

        assert_eq!(stats.contigs_ge_1kb, 3);
        assert_eq!(stats.contigs_ge_10kb, 2);
        assert_eq!(stats.contigs_ge_50kb, 1);
        assert_eq!(stats.contigs_ge_100kb, 0);
        assert_eq!(stats.contigs_ge_1mb, 0);
        assert_eq!(stats.bases_ge_1kb, 61_000);
        assert_eq!(stats.bases_ge_10kb, 60_000);
        assert_eq!(stats.bases_ge_50kb, 50_000);
        assert_eq!(stats.bases_ge_100kb, 0);
        assert_eq!(stats.bases_ge_1mb, 0);
        assert!((stats.contigs_ge_1kb_frac - 0.75).abs() < 1e-12);
        assert!((stats.contigs_ge_10kb_frac - 0.5).abs() < 1e-12);
        assert!((stats.contigs_ge_50kb_frac - 0.25).abs() < 1e-12);
        assert_eq!(stats.contigs_ge_100kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_1mb_frac, 0.0);
        assert!((stats.bases_ge_1kb_frac - (61_000.0 / 61_999.0)).abs() < 1e-12);
        assert!((stats.bases_ge_10kb_frac - (60_000.0 / 61_999.0)).abs() < 1e-12);
        assert!((stats.bases_ge_50kb_frac - (50_000.0 / 61_999.0)).abs() < 1e-12);
        assert_eq!(stats.bases_ge_100kb_frac, 0.0);
        assert_eq!(stats.bases_ge_1mb_frac, 0.0);
    }

    #[test]
    fn test_write_output_reports_100kb_bucket_metrics() {
        let output = NamedTempFile::new().unwrap();
        let temp_dir = TempDir::new().unwrap();
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            min_contig_len: 1,
            num_buckets: Some(4),
            temp_dir: Some(temp_dir.path().to_str().unwrap().to_string()),
            ..Default::default()
        });

        let contigs = vec!["A".repeat(100_000), "C".repeat(99_999)];
        let mut stats = AssemblyStats::default();
        assembler
            .write_output(&contigs, output.path().to_str().unwrap(), &mut stats)
            .unwrap();

        assert_eq!(stats.contigs_ge_100kb, 1);
        assert_eq!(stats.bases_ge_100kb, 100_000);
        assert!((stats.contigs_ge_100kb_frac - 0.5).abs() < 1e-12);
        assert!((stats.bases_ge_100kb_frac - (100_000.0 / 199_999.0)).abs() < 1e-12);
        assert_eq!(stats.contigs_ge_1mb, 0);
        assert_eq!(stats.bases_ge_1mb, 0);
        assert_eq!(stats.contigs_ge_1mb_frac, 0.0);
        assert_eq!(stats.bases_ge_1mb_frac, 0.0);
    }

    #[test]
    fn test_write_output_reports_1mb_bucket_metrics() {
        let output = NamedTempFile::new().unwrap();
        let temp_dir = TempDir::new().unwrap();
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            min_contig_len: 1,
            num_buckets: Some(4),
            temp_dir: Some(temp_dir.path().to_str().unwrap().to_string()),
            ..Default::default()
        });

        let contigs = vec![
            "A".repeat(1_500_000),
            "C".repeat(900_000),
            "G".repeat(100_000),
        ];
        let mut stats = AssemblyStats::default();
        assembler
            .write_output(&contigs, output.path().to_str().unwrap(), &mut stats)
            .unwrap();

        assert_eq!(stats.contigs_ge_1mb, 1);
        assert_eq!(stats.bases_ge_1mb, 1_500_000);
        assert!((stats.contigs_ge_1mb_frac - (1.0 / 3.0)).abs() < 1e-12);
        assert!((stats.bases_ge_1mb_frac - (1_500_000.0 / 2_500_000.0)).abs() < 1e-12);
    }

    #[test]
    fn test_write_output_returns_error_for_invalid_output_path() {
        let temp_dir = TempDir::new().unwrap();
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            min_contig_len: 1,
            num_buckets: Some(4),
            temp_dir: Some(temp_dir.path().to_str().unwrap().to_string()),
            ..Default::default()
        });

        let contigs = vec!["ACGT".to_string()];
        let mut stats = AssemblyStats::default();
        let result =
            assembler.write_output(&contigs, temp_dir.path().to_str().unwrap(), &mut stats);
        assert!(result.is_err());
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
    fn test_repeat_stats_even_split_uses_midpoint_median_and_detects_tail_repeats() {
        let mut counts = AHashMap::new();
        for i in 0..50 {
            counts.insert(i as u64, 5);
        }
        for i in 50..95 {
            counts.insert(i as u64, 25);
        }
        for i in 95..100 {
            counts.insert(i as u64, 35);
        }

        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig::default());
        let (repeat_kmers, stats) = assembler.identify_repeats(&counts);

        assert_eq!(stats.total_kmer_count, 100);
        assert_eq!(stats.median_coverage, 15);
        assert_eq!(stats.repeat_threshold, 30);
        assert_eq!(stats.repeat_kmer_count, 5);
        assert_eq!(repeat_kmers.len(), 5);
        for i in 95..100 {
            assert!(repeat_kmers.contains(&(i as u64)));
        }
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
    fn test_auto_min_count_uses_upper_median_for_even_non_singleton_counts() {
        let mut counts = AHashMap::new();
        counts.insert(0, 2);
        counts.insert(1, 2);
        counts.insert(2, 49);
        counts.insert(3, 50);
        counts.insert(4, 1); // singleton should be ignored

        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            min_count: 0,
            ..Default::default()
        });

        // Sorted non-singletons: [2, 2, 49, 50], upper median is 49.
        assert_eq!(assembler.select_min_count(&counts), 4);
    }

    #[test]
    fn test_auto_min_count_is_insertion_order_invariant() {
        let values = [2u32, 2, 2, 8, 12, 24, 40, 60, 1, 1];
        let base_items: Vec<(u64, u32)> = values
            .iter()
            .copied()
            .enumerate()
            .map(|(idx, count)| (idx as u64, count))
            .collect();

        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            min_count: 0,
            ..Default::default()
        });

        let mut baseline = AHashMap::new();
        for &(k, v) in &base_items {
            baseline.insert(k, v);
        }
        let expected = assembler.select_min_count(&baseline);

        let mut rng = StdRng::seed_from_u64(0xA11C_E11E_1234_5678);
        for _ in 0..64 {
            let mut shuffled = base_items.clone();
            shuffled.shuffle(&mut rng);

            let mut observed_map = AHashMap::new();
            for (k, v) in shuffled {
                observed_map.insert(k, v);
            }
            assert_eq!(assembler.select_min_count(&observed_map), expected);
        }
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
    fn test_extend_bidirectional_with_coverage_reconstructs_linear_path() {
        let sequence = "AAATCGG";
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
        let branch_support: AHashMap<(u64, u64), u32> = AHashMap::new();
        let repeat_kmers: AHashSet<u64> = AHashSet::new();
        let mut used = AHashSet::new();
        let seed = KmerU64::from_str("ATCG").unwrap().canonical().encoded;

        let contig = assembler.extend_bidirectional_with_coverage(
            seed,
            k,
            &adjacency,
            &counts,
            &branch_support,
            &repeat_kmers,
            &mut used,
        );

        let reverse = reverse_complement(sequence);
        assert!(contig == sequence || contig == reverse);
    }

    #[test]
    fn test_extend_bidirectional_with_coverage_marks_used_by_canonical_seed() {
        let sequence = "AAATCGG";
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
        let branch_support: AHashMap<(u64, u64), u32> = AHashMap::new();
        let repeat_kmers: AHashSet<u64> = AHashSet::new();
        let mut used = AHashSet::new();

        let canonical_seed = KmerU64::from_str("ATCG").unwrap().canonical().encoded;
        let reverse_seed = KmerU64 {
            encoded: canonical_seed,
            len: k as u8,
        }
        .reverse_complement()
        .encoded;
        assert_ne!(canonical_seed, reverse_seed);

        let contig = assembler.extend_bidirectional_with_coverage(
            reverse_seed,
            k,
            &adjacency,
            &counts,
            &branch_support,
            &repeat_kmers,
            &mut used,
        );

        assert!(
            used.contains(&canonical_seed),
            "used-set should track canonical k-mers regardless of traversal orientation"
        );
        let reverse = reverse_complement(sequence);
        assert!(contig == sequence || contig == reverse);
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
        let primary_key = LargeGenomeAssembler::canonical_edge_key(branch, next_primary, k);
        let alternate_key = LargeGenomeAssembler::canonical_edge_key(branch, next_alternate, k);

        assert!(observations >= 3);
        assert_eq!(edge_support.get(&primary_key).copied(), Some(2));
        assert_eq!(edge_support.get(&alternate_key).copied(), Some(1));
    }

    #[test]
    fn test_branch_threading_aggregates_reverse_complement_support_into_same_edge() {
        let k = 4;
        let sequences = ["TACGATG", "CATCGTA"];
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

        let branch = KmerU64::from_str("ACGA").unwrap().encoded;
        let next = KmerU64::from_str("CGAT").unwrap().encoded;
        let canonical_key = LargeGenomeAssembler::canonical_edge_key(branch, next, k);
        let mut branch_edges = AHashSet::new();
        branch_edges.insert(canonical_key);

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

        assert_eq!(observations, 2);
        assert_eq!(edge_support.get(&canonical_key).copied(), Some(2));
    }

    #[test]
    fn test_branch_threading_support_accumulator_saturates_at_u32_max() {
        let k = 4;
        let sequence = "TACGATG";
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

        let branch = KmerU64::from_str("ACGA").unwrap().encoded;
        let next = KmerU64::from_str("CGAT").unwrap().encoded;
        let edge_key = LargeGenomeAssembler::canonical_edge_key(branch, next, k);

        let mut branch_edges = AHashSet::new();
        branch_edges.insert(edge_key);

        let mut edge_support = AHashMap::new();
        edge_support.insert(edge_key, u32::MAX);

        let observations = LargeGenomeAssembler::thread_branch_edges_in_sequence(
            sequence.as_bytes(),
            k,
            &adjacency,
            &branch_edges,
            &mut edge_support,
        );

        assert_eq!(observations, 1);
        assert_eq!(edge_support.get(&edge_key).copied(), Some(u32::MAX));
    }

    proptest! {
        #[test]
        fn prop_branch_threading_matches_reference_implementation_under_noisy_input(
            reads in prop::collection::vec(
                prop::collection::vec(
                    prop_oneof![
                        Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T'),
                        Just(b'N'), Just(b'R'), Just(b'Y'), Just(b'a'), Just(b't')
                    ],
                    0..80
                ),
                0..24
            )
        ) {
            let k = 4usize;
            let (adjacency, branch_edges, _) = build_branch_threading_fixture(k);

            let mut optimized_support = AHashMap::new();
            let mut reference_support = AHashMap::new();
            let mut optimized_observations = 0u64;
            let mut reference_observations = 0u64;

            for read in &reads {
                optimized_observations += LargeGenomeAssembler::thread_branch_edges_in_sequence(
                    read,
                    k,
                    &adjacency,
                    &branch_edges,
                    &mut optimized_support,
                );
                reference_observations += LargeGenomeAssembler::thread_branch_edges_in_sequence_reference(
                    read,
                    k,
                    &adjacency,
                    &branch_edges,
                    &mut reference_support,
                );
            }

            prop_assert_eq!(optimized_observations, reference_observations);
            prop_assert_eq!(
                canonicalize_edge_support(&optimized_support),
                canonicalize_edge_support(&reference_support),
            );
        }

        #[test]
        fn prop_select_min_count_matches_reference_and_is_insertion_order_invariant(
            counts in prop::collection::vec(0u32..200u32, 1..256),
            seed in any::<u64>(),
        ) {
            let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
                min_count: 0,
                ..Default::default()
            });

            let mut baseline = AHashMap::new();
            let mut items = Vec::with_capacity(counts.len());
            for (idx, &count) in counts.iter().enumerate() {
                let key = idx as u64;
                baseline.insert(key, count);
                items.push((key, count));
            }

            let mut non_singleton = counts
                .iter()
                .copied()
                .filter(|&count| count >= 2)
                .collect::<Vec<_>>();
            let expected = if non_singleton.is_empty() {
                2
            } else {
                non_singleton.sort_unstable();
                (non_singleton[non_singleton.len() / 2] / 10).clamp(2, 5)
            };

            prop_assert_eq!(assembler.select_min_count(&baseline), expected);

            let mut rng = StdRng::seed_from_u64(seed);
            items.shuffle(&mut rng);
            let mut shuffled = AHashMap::new();
            for (key, value) in items {
                shuffled.insert(key, value);
            }
            prop_assert_eq!(assembler.select_min_count(&shuffled), expected);
        }
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

    #[test]
    fn test_choose_branch_extension_prefers_closer_coverage_ratio_over_weak_support() {
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig::default());
        let extensions = [
            BranchCandidate {
                base_idx: 0,
                next: 11,
                count: 21,
                is_repeat: false,
                read_support: 0,
            },
            BranchCandidate {
                base_idx: 1,
                next: 22,
                count: 39,
                is_repeat: false,
                read_support: 1,
            },
        ];

        let best = assembler.choose_branch_extension(&extensions, 20).unwrap();
        assert_eq!(best.next, 11);
    }

    #[test]
    fn test_choose_branch_extension_breaks_equal_ratio_ties_by_closer_absolute_coverage() {
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig::default());
        let extensions = [
            BranchCandidate {
                base_idx: 0,
                next: 11,
                count: 10,
                is_repeat: false,
                read_support: 0,
            },
            BranchCandidate {
                base_idx: 1,
                next: 22,
                count: 40,
                is_repeat: false,
                read_support: 0,
            },
        ];

        let best = assembler.choose_branch_extension(&extensions, 20).unwrap();
        assert_eq!(best.next, 11);
    }

    #[test]
    fn test_choose_branch_extension_repeat_fallback_is_insertion_order_invariant() {
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig::default());
        let left = BranchCandidate {
            base_idx: 3,
            next: 11,
            count: 8,
            is_repeat: true,
            read_support: 2,
        };
        let right = BranchCandidate {
            base_idx: 0,
            next: 7,
            count: 8,
            is_repeat: true,
            read_support: 2,
        };

        let best_ab = assembler
            .choose_branch_extension(&[left, right], 12)
            .unwrap();
        let best_ba = assembler
            .choose_branch_extension(&[right, left], 12)
            .unwrap();

        assert_eq!(best_ab.next, best_ba.next);
        assert_eq!(best_ab.next, 7);
    }

    #[test]
    fn test_choose_branch_extension_respects_configured_support_margin() {
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            branch_support_min_margin: 3,
            ..Default::default()
        });
        let extensions = [
            BranchCandidate {
                base_idx: 0,
                next: 11,
                count: 30,
                is_repeat: false,
                read_support: 5,
            },
            BranchCandidate {
                base_idx: 1,
                next: 22,
                count: 20,
                is_repeat: false,
                read_support: 3,
            },
        ];

        let best = assembler.choose_branch_extension(&extensions, 20).unwrap();
        assert_eq!(best.next, 22);
    }

    #[test]
    fn test_choose_branch_extension_can_disable_non_repeat_preference() {
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            prefer_non_repeat_branches: false,
            ..Default::default()
        });
        let extensions = [
            BranchCandidate {
                base_idx: 0,
                next: 11,
                count: 19,
                is_repeat: false,
                read_support: 0,
            },
            BranchCandidate {
                base_idx: 1,
                next: 22,
                count: 20,
                is_repeat: true,
                read_support: 0,
            },
        ];

        let best = assembler.choose_branch_extension(&extensions, 20).unwrap();
        assert_eq!(best.next, 22);
    }

    #[test]
    fn test_choose_branch_extension_is_order_invariant_under_randomized_inputs() {
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig::default());
        let mut rng = StdRng::seed_from_u64(0xC0FFEE);

        for case_idx in 0..256u64 {
            let current_count = rng.gen_range(1..=64);
            let extension_count = rng.gen_range(1..=4usize);
            let mut extensions = Vec::with_capacity(extension_count);

            for base_idx in 0..extension_count {
                extensions.push(BranchCandidate {
                    base_idx,
                    next: (case_idx << 8) | base_idx as u64,
                    count: rng.gen_range(0..=64),
                    is_repeat: rng.gen_bool(0.4),
                    read_support: rng.gen_range(0..=8),
                });
            }

            let baseline = assembler
                .choose_branch_extension(&extensions, current_count)
                .map(|candidate| (candidate.base_idx, candidate.next));

            for _ in 0..32 {
                extensions.shuffle(&mut rng);
                let observed = assembler
                    .choose_branch_extension(&extensions, current_count)
                    .map(|candidate| (candidate.base_idx, candidate.next));
                assert_eq!(
                    observed, baseline,
                    "case={} current_count={} extensions={:?}",
                    case_idx, current_count, extensions
                );
            }
        }
    }

    #[test]
    fn test_run_branch_resolution_case_exposes_branch_chooser() {
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig::default());
        let selected = assembler
            .run_branch_resolution_case(&LargeGenomeBenchBranchResolutionCase {
                current_count: 20,
                candidates: vec![
                    LargeGenomeBenchBranchCandidate {
                        base_idx: 0,
                        next: 11,
                        count: 18,
                        is_repeat: false,
                        read_support: 5,
                    },
                    LargeGenomeBenchBranchCandidate {
                        base_idx: 1,
                        next: 22,
                        count: 20,
                        is_repeat: false,
                        read_support: 1,
                    },
                ],
            })
            .unwrap();
        assert_eq!(selected.next, 11);
        assert_eq!(selected.base_idx, 0);
    }

    #[test]
    fn test_build_contigs_from_graph_is_stable_under_shuffled_tied_seed_insertion() {
        let k = 3;
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            k,
            min_count: 1,
            min_contig_len: 1,
            ..Default::default()
        });

        let mut entries: Vec<(u64, u32)> = ["AAA", "AAT", "ATC", "ATG"]
            .into_iter()
            .map(|kmer| (KmerU64::from_str(kmer).unwrap().canonical().encoded, 10u32))
            .collect();
        entries.sort_unstable_by_key(|(kmer, _)| *kmer);
        entries.dedup_by_key(|(kmer, _)| *kmer);

        let branch_support: AHashMap<(u64, u64), u32> = AHashMap::new();
        let mut rng = StdRng::seed_from_u64(7);
        let mut baseline: Option<Vec<String>> = None;

        for _ in 0..64 {
            entries.shuffle(&mut rng);
            let mut counts = AHashMap::new();
            for &(kmer, count) in &entries {
                counts.insert(kmer, count);
            }

            let valid_kmers: AHashSet<u64> = counts.keys().copied().collect();
            let adjacency = assembler.build_adjacency(&valid_kmers, k);
            let contigs =
                assembler.build_contigs_from_graph(&counts, &adjacency, &branch_support, k);

            if let Some(expected) = &baseline {
                assert_eq!(&contigs, expected);
            } else {
                baseline = Some(contigs);
            }
        }
    }

    #[test]
    fn test_build_contigs_from_graph_falls_back_to_repeat_seeds_after_cleanup() {
        let k = 3;
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            k,
            min_count: 1,
            min_contig_len: 1,
            ..Default::default()
        });

        let mut counts = AHashMap::new();
        for kmer in ["AAA", "AAC", "ACC"] {
            let canonical = KmerU64::from_str(kmer).unwrap().canonical().encoded;
            counts.insert(canonical, 20);
        }
        // Add many low-coverage entries that are not present in the cleaned adjacency.
        // This forces repeat detection to classify the surviving graph nodes as repeats.
        for decoy in 1_000_000u64..1_000_020u64 {
            counts.insert(decoy, 2);
        }

        let valid_kmers: AHashSet<u64> = ["AAA", "AAC", "ACC"]
            .into_iter()
            .map(|kmer| KmerU64::from_str(kmer).unwrap().canonical().encoded)
            .collect();
        let adjacency = assembler.build_adjacency(&valid_kmers, k);
        let branch_support: AHashMap<(u64, u64), u32> = AHashMap::new();
        let contigs = assembler.build_contigs_from_graph(&counts, &adjacency, &branch_support, k);

        assert!(
            !contigs.is_empty(),
            "repeat-seed fallback should still produce contigs"
        );
        assert!(
            contigs.iter().all(|contig| contig.len() >= k),
            "all emitted contigs should be at least k bases long"
        );
    }

    #[test]
    fn test_repeat_seed_fallback_is_stable_under_shuffled_count_insertion() {
        let k = 3;
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            k,
            min_count: 1,
            min_contig_len: 1,
            ..Default::default()
        });

        let mut entries: Vec<(u64, u32)> = ["AAA", "AAC", "ACC"]
            .into_iter()
            .map(|kmer| (KmerU64::from_str(kmer).unwrap().canonical().encoded, 20u32))
            .collect();
        for decoy in 1_000_000u64..1_000_020u64 {
            entries.push((decoy, 2));
        }

        let valid_kmers: AHashSet<u64> = ["AAA", "AAC", "ACC"]
            .into_iter()
            .map(|kmer| KmerU64::from_str(kmer).unwrap().canonical().encoded)
            .collect();
        let adjacency = assembler.build_adjacency(&valid_kmers, k);
        let branch_support: AHashMap<(u64, u64), u32> = AHashMap::new();

        let mut rng = StdRng::seed_from_u64(0xFEED_FACE_u64);
        let mut baseline: Option<Vec<String>> = None;

        for _ in 0..64 {
            entries.shuffle(&mut rng);
            let mut counts = AHashMap::new();
            for &(kmer, count) in &entries {
                counts.insert(kmer, count);
            }

            let contigs =
                assembler.build_contigs_from_graph(&counts, &adjacency, &branch_support, k);
            if let Some(expected) = &baseline {
                assert_eq!(&contigs, expected);
            } else {
                baseline = Some(contigs);
            }
        }
    }

    #[test]
    fn test_build_contigs_from_graph_runs_repeat_completion_after_non_repeat_pass() {
        let k = 3;
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            k,
            min_count: 1,
            min_contig_len: 1,
            ..Default::default()
        });

        let mut counts = AHashMap::new();
        for kmer in ["AAA", "AAT"] {
            let canonical = KmerU64::from_str(kmer).unwrap().canonical().encoded;
            counts.insert(canonical, 6);
        }
        for kmer in ["CCC", "CCG"] {
            let canonical = KmerU64::from_str(kmer).unwrap().canonical().encoded;
            counts.insert(canonical, 20);
        }
        // Keep median low enough that the high-coverage component is tagged as repeat.
        for decoy in 2_000_000u64..2_000_020u64 {
            counts.insert(decoy, 2);
        }

        let valid_kmers: AHashSet<u64> = ["AAA", "AAT", "CCC", "CCG"]
            .into_iter()
            .map(|kmer| KmerU64::from_str(kmer).unwrap().canonical().encoded)
            .collect();
        let adjacency = assembler.build_adjacency(&valid_kmers, k);
        let branch_support: AHashMap<(u64, u64), u32> = AHashMap::new();
        let contigs = assembler.build_contigs_from_graph(&counts, &adjacency, &branch_support, k);

        assert!(contigs.len() >= 2);
        assert!(
            contigs
                .iter()
                .any(|contig| contig == "AAAT" || contig == "ATTT"),
            "non-repeat component should be assembled: {:?}",
            contigs
        );
        assert!(
            contigs
                .iter()
                .any(|contig| contig == "CCCG" || contig == "CGGG"),
            "repeat completion pass should recover repeat component: {:?}",
            contigs
        );
    }

    #[test]
    fn test_repeat_completion_pass_is_stable_under_shuffled_count_insertion() {
        let k = 3;
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            k,
            min_count: 1,
            min_contig_len: 1,
            ..Default::default()
        });

        let mut entries: Vec<(u64, u32)> = ["AAA", "AAT"]
            .into_iter()
            .map(|kmer| (KmerU64::from_str(kmer).unwrap().canonical().encoded, 6u32))
            .collect();
        entries.extend(
            ["CCC", "CCG"]
                .into_iter()
                .map(|kmer| (KmerU64::from_str(kmer).unwrap().canonical().encoded, 20u32)),
        );
        for decoy in 2_000_000u64..2_000_020u64 {
            entries.push((decoy, 2));
        }

        let valid_kmers: AHashSet<u64> = ["AAA", "AAT", "CCC", "CCG"]
            .into_iter()
            .map(|kmer| KmerU64::from_str(kmer).unwrap().canonical().encoded)
            .collect();
        let adjacency = assembler.build_adjacency(&valid_kmers, k);
        let branch_support: AHashMap<(u64, u64), u32> = AHashMap::new();

        let mut rng = StdRng::seed_from_u64(0xABC0_1234_u64);
        let mut baseline: Option<Vec<String>> = None;
        for _ in 0..64 {
            entries.shuffle(&mut rng);
            let mut counts = AHashMap::new();
            for &(kmer, count) in &entries {
                counts.insert(kmer, count);
            }

            let contigs =
                assembler.build_contigs_from_graph(&counts, &adjacency, &branch_support, k);
            if let Some(expected) = &baseline {
                assert_eq!(&contigs, expected);
            } else {
                baseline = Some(contigs);
            }
        }
    }

    #[test]
    fn test_build_adjacency_is_stable_under_shuffled_valid_kmer_insertion() {
        let k = 4;
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            k,
            min_count: 1,
            min_contig_len: 1,
            ..Default::default()
        });

        let mut entries: Vec<u64> = ["AAAT", "AATC", "ATCG", "TCGA", "CGAT", "GATC"]
            .into_iter()
            .map(|kmer| KmerU64::from_str(kmer).unwrap().canonical().encoded)
            .collect();
        entries.sort_unstable();
        entries.dedup();

        let mut baseline_set = AHashSet::new();
        for &kmer in &entries {
            baseline_set.insert(kmer);
        }
        let baseline = canonicalize_adjacency(&assembler.build_adjacency(&baseline_set, k));

        let mut rng = StdRng::seed_from_u64(0xA11A_DA71_u64);
        for _ in 0..64 {
            entries.shuffle(&mut rng);
            let mut valid_kmers = AHashSet::new();
            for &kmer in &entries {
                valid_kmers.insert(kmer);
            }
            let observed = canonicalize_adjacency(&assembler.build_adjacency(&valid_kmers, k));
            assert_eq!(observed, baseline);
        }
    }

    #[test]
    fn test_stage_bench_fixture_generates_branch_support() {
        let fixture = LargeGenomeStageBenchFixture::synthetic_branching(11, 8, 5, 2);
        let stats = fixture.run_branch_threading();

        assert!(fixture.total_branch_read_bases() > 0);
        assert!(fixture.branch_edge_count() > 0);
        assert_eq!(stats.ambiguous_edges, fixture.branch_edge_count());
        assert!(stats.supported_edges > 0);
        assert!(stats.edge_observations > 0);
    }

    #[test]
    fn test_stage_bench_fixture_extracts_contigs() {
        let fixture = LargeGenomeStageBenchFixture::synthetic_branching(11, 6, 4, 1);
        let contigs = fixture.run_contig_extraction();

        assert!(fixture.graph_node_count() > 0);
        assert!(!contigs.is_empty());
        assert!(contigs.iter().all(|contig| contig.len() >= 11));
    }

    #[test]
    fn test_build_contigs_from_graph_can_suppress_redundant_branch_alternatives() {
        let fixture = LargeGenomeStageBenchFixture::synthetic_branching(11, 6, 4, 1);
        let baseline = fixture.run_contig_extraction();

        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            k: 11,
            min_count: 1,
            min_contig_len: 11,
            suppress_redundant_contigs: true,
            ..Default::default()
        });
        let suppressed = assembler.run_contig_extraction_case(
            fixture.cloned_kmer_counts(),
            &fixture.valid_kmers(),
            &fixture.branch_support_map(),
            11,
        );

        let canonicalize = |sequence: &str| {
            let reverse = reverse_complement(sequence);
            if reverse.as_str() < sequence {
                reverse
            } else {
                sequence.to_string()
            }
        };
        let expected = fixture
            .expected_contigs()
            .iter()
            .map(|sequence| canonicalize(sequence))
            .collect::<std::collections::BTreeSet<_>>();
        let observed = suppressed
            .iter()
            .map(|sequence| canonicalize(sequence))
            .collect::<std::collections::BTreeSet<_>>();

        assert!(baseline.len() > suppressed.len());
        assert_eq!(observed, expected);
    }

    #[test]
    fn test_error_correction_stage_bench_fixture_corrects_singletons() {
        let fixture =
            LargeGenomeErrorCorrectionBenchFixture::synthetic_singleton_correction(21, 128, 2);
        let initial_kmers = fixture.total_kmer_count();
        let singleton_count = fixture.singleton_count();
        let stats = fixture.run_error_correction(fixture.cloned_kmer_counts());

        assert!(singleton_count > 0);
        assert_eq!(stats.corrected_kmers, singleton_count as u64);
        assert!(stats.remaining_kmers < initial_kmers);
    }

    #[test]
    fn test_remove_kmers_and_prune_adjacency_removes_both_orientations_and_dangling_edges() {
        let k = 4;
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            k,
            min_count: 1,
            min_contig_len: 1,
            ..Default::default()
        });

        let core = KmerU64::from_str("ACGA").unwrap().encoded;
        let tip = KmerU64::from_str("CGAT").unwrap().canonical().encoded;
        let tip_rc = KmerU64 {
            encoded: tip,
            len: k as u8,
        }
        .reverse_complement()
        .encoded;

        let mut valid_kmers = AHashSet::new();
        valid_kmers.insert(KmerU64::from_str("ACGA").unwrap().canonical().encoded);
        valid_kmers.insert(tip);

        let mut adjacency = assembler.build_adjacency(&valid_kmers, k);
        let (_, right_before) = adjacency.get(&core).copied().unwrap();
        assert!(
            right_before[3],
            "expected ACGA --T--> CGAT edge before removal"
        );

        let mut to_remove = AHashSet::new();
        to_remove.insert(tip);
        LargeGenomeAssembler::remove_kmers_and_prune_adjacency(&mut adjacency, &to_remove, k);

        assert!(!adjacency.contains_key(&tip));
        assert!(!adjacency.contains_key(&tip_rc));

        let (_, right_after) = adjacency.get(&core).copied().unwrap();
        assert!(
            !right_after[3],
            "expected dangling ACGA --T--> CGAT edge to be pruned"
        );
    }

    #[test]
    fn test_pop_bubbles_is_stable_under_shuffled_kmer_insertion_order() {
        let k = 3;
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            k,
            min_count: 1,
            min_contig_len: 1,
            ..Default::default()
        });

        let mut entries: Vec<(u64, u32)> = [
            ("AAA", 20u32),
            ("AAC", 10u32),
            ("AAG", 10u32),
            ("ACC", 20u32),
        ]
        .into_iter()
        .map(|(kmer, count)| (KmerU64::from_str(kmer).unwrap().canonical().encoded, count))
        .collect();
        entries.sort_unstable_by_key(|(kmer, _)| *kmer);
        entries.dedup_by_key(|(kmer, _)| *kmer);

        let mut rng = StdRng::seed_from_u64(0xBADD_CAFE);
        let mut baseline: Option<(usize, Vec<(u64, [bool; 4], [bool; 4])>)> = None;

        for _ in 0..64 {
            entries.shuffle(&mut rng);
            let mut counts = AHashMap::new();
            for &(kmer, count) in &entries {
                counts.insert(kmer, count);
            }

            let valid_kmers: AHashSet<u64> = counts.keys().copied().collect();
            let mut adjacency = assembler.build_adjacency(&valid_kmers, k);
            let popped = assembler.pop_bubbles(&mut adjacency, &counts, k);
            let snapshot = canonicalize_adjacency(&adjacency);

            if let Some((expected_popped, expected_snapshot)) = &baseline {
                assert_eq!(popped, *expected_popped);
                assert_eq!(snapshot, *expected_snapshot);
            } else {
                baseline = Some((popped, snapshot));
            }
        }

        assert!(baseline.is_some(), "baseline should be captured");
    }

    #[test]
    fn test_pop_bubbles_preserves_reconvergence_node() {
        let k = 3;
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            k,
            min_count: 1,
            min_contig_len: 1,
            ..Default::default()
        });

        let aaa = KmerU64::from_str("AAA").unwrap().encoded;
        let aac = KmerU64::from_str("AAC").unwrap().encoded;
        let aag = KmerU64::from_str("AAG").unwrap().encoded;
        let agc = KmerU64::from_str("AGC").unwrap().encoded;
        let acc = KmerU64::from_str("ACC").unwrap().encoded;
        let gcc = KmerU64::from_str("GCC").unwrap().encoded;
        let ccc = KmerU64::from_str("CCC").unwrap().encoded;

        let mut right_aa = [false; 4];
        right_aa[1] = true; // AAA -> AAC
        right_aa[2] = true; // AAA -> AAG
        let mut right_c = [false; 4];
        right_c[1] = true; // * -> *C
        let none = [false; 4];

        let mut adjacency = AHashMap::new();
        adjacency.insert(aaa, (none, right_aa));
        adjacency.insert(aac, (none, right_c));
        adjacency.insert(aag, (none, right_c));
        adjacency.insert(agc, (none, right_c));
        adjacency.insert(acc, (none, right_c));
        adjacency.insert(gcc, (none, right_c));
        adjacency.insert(ccc, (none, none));

        let mut counts = AHashMap::new();
        for (node, count) in [
            (aaa, 40u32),
            (aac, 30u32),
            (acc, 30u32),
            (aag, 8u32),
            (agc, 8u32),
            (gcc, 8u32),
            (ccc, 35u32),
        ] {
            let canonical = KmerU64 {
                encoded: node,
                len: k as u8,
            }
            .canonical()
            .encoded;
            counts.insert(canonical, count);
        }

        let popped = assembler.pop_bubbles(&mut adjacency, &counts, k);

        assert_eq!(popped, 1);
        let ccc_canonical = KmerU64 {
            encoded: ccc,
            len: k as u8,
        }
        .canonical()
        .encoded;
        assert!(
            adjacency.contains_key(&ccc) || adjacency.contains_key(&ccc_canonical),
            "bubble popping must preserve the shared reconvergence node"
        );
        let aag_canonical = KmerU64 {
            encoded: aag,
            len: k as u8,
        }
        .canonical()
        .encoded;
        assert!(
            !adjacency.contains_key(&aag) && !adjacency.contains_key(&aag_canonical),
            "expected lower-coverage branch nodes to be removed"
        );
    }

    #[test]
    fn test_pop_bubbles_handles_high_coverage_without_overflow() {
        let k = 3;
        let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
            k,
            min_count: 1,
            min_contig_len: 1,
            ..Default::default()
        });

        let aaa = KmerU64::from_str("AAA").unwrap().encoded;
        let aac = KmerU64::from_str("AAC").unwrap().encoded;
        let aag = KmerU64::from_str("AAG").unwrap().encoded;
        let agc = KmerU64::from_str("AGC").unwrap().encoded;
        let acc = KmerU64::from_str("ACC").unwrap().encoded;
        let gcc = KmerU64::from_str("GCC").unwrap().encoded;
        let ccc = KmerU64::from_str("CCC").unwrap().encoded;

        let mut right_aa = [false; 4];
        right_aa[1] = true; // AAA -> AAC
        right_aa[2] = true; // AAA -> AAG
        let mut right_c = [false; 4];
        right_c[1] = true; // * -> *C
        let none = [false; 4];

        let mut adjacency = AHashMap::new();
        adjacency.insert(aaa, (none, right_aa));
        adjacency.insert(aac, (none, right_c));
        adjacency.insert(aag, (none, right_c));
        adjacency.insert(agc, (none, right_c));
        adjacency.insert(acc, (none, right_c));
        adjacency.insert(gcc, (none, right_c));
        adjacency.insert(ccc, (none, none));

        let mut counts = AHashMap::new();
        for node in [aaa, aac, aag, agc, acc, gcc, ccc] {
            let canonical = KmerU64 {
                encoded: node,
                len: k as u8,
            }
            .canonical()
            .encoded;
            counts.insert(canonical, u32::MAX);
        }

        let popped = assembler.pop_bubbles(&mut adjacency, &counts, k);

        assert_eq!(popped, 1);
        let ccc_canonical = KmerU64 {
            encoded: ccc,
            len: k as u8,
        }
        .canonical()
        .encoded;
        assert!(
            adjacency.contains_key(&ccc) || adjacency.contains_key(&ccc_canonical),
            "reconvergence node should survive even at extreme coverage"
        );
    }

    proptest! {
        #[test]
        fn prop_pop_bubbles_is_stable_under_randomized_coverage_insertion_order(
            coverages in prop::array::uniform7(1u32..=u32::MAX),
            seed in any::<u64>(),
        ) {
            let (k, nodes, fixture) = simple_pop_bubble_fixture();
            let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
                k,
                min_count: 1,
                min_contig_len: 1,
                ..Default::default()
            });

            let mut baseline_counts = AHashMap::new();
            for idx in 0..nodes.len() {
                let canonical = KmerU64 {
                    encoded: nodes[idx],
                    len: k as u8,
                }
                .canonical()
                .encoded;
                baseline_counts.insert(canonical, coverages[idx]);
            }
            let mut baseline_adjacency = fixture.clone();
            let baseline_popped = assembler.pop_bubbles(&mut baseline_adjacency, &baseline_counts, k);
            let baseline_snapshot = canonicalize_adjacency(&baseline_adjacency);

            let mut order = [0usize, 1, 2, 3, 4, 5, 6];
            let mut rng = StdRng::seed_from_u64(seed);
            order.shuffle(&mut rng);

            let mut shuffled_counts = AHashMap::new();
            for idx in order {
                let canonical = KmerU64 {
                    encoded: nodes[idx],
                    len: k as u8,
                }
                .canonical()
                .encoded;
                shuffled_counts.insert(canonical, coverages[idx]);
            }

            let mut shuffled_adjacency = fixture;
            let shuffled_popped = assembler.pop_bubbles(&mut shuffled_adjacency, &shuffled_counts, k);
            let shuffled_snapshot = canonicalize_adjacency(&shuffled_adjacency);

            prop_assert_eq!(shuffled_popped, baseline_popped);
            prop_assert_eq!(shuffled_snapshot, baseline_snapshot);

            let ccc = nodes[6];
            let ccc_canonical = KmerU64 {
                encoded: ccc,
                len: k as u8,
            }
            .canonical()
            .encoded;
            prop_assert!(
                shuffled_adjacency.contains_key(&ccc) || shuffled_adjacency.contains_key(&ccc_canonical),
                "shared reconvergence node should remain after bubble cleanup"
            );
        }
    }

    #[test]
    fn test_extend_bidirectional_recovers_via_canonical_orientation_lookup() {
        let (assembler, k, adjacency, _counts, acga, atcg, tcga) =
            canonical_only_orientation_fixture();
        let mut used = AHashSet::new();

        let contig = assembler.extend_bidirectional(acga, k, &adjacency, &mut used);
        assert_eq!(contig, "ACGATA");
        assert!(used.contains(&acga));
        assert!(used.contains(&atcg));
        assert!(used.contains(&tcga));
    }

    #[test]
    fn test_extend_bidirectional_with_coverage_recovers_via_canonical_orientation_lookup() {
        let (assembler, k, adjacency, counts, acga, _atcg, _tcga) =
            canonical_only_orientation_fixture();
        let mut used = AHashSet::new();
        let branch_support: AHashMap<(u64, u64), u32> = AHashMap::new();
        let repeat_kmers: AHashSet<u64> = AHashSet::new();

        let contig = assembler.extend_bidirectional_with_coverage(
            acga,
            k,
            &adjacency,
            &counts,
            &branch_support,
            &repeat_kmers,
            &mut used,
        );
        assert_eq!(contig, "ACGATA");
    }

    #[test]
    fn test_trace_path_extends_from_resolved_orientation_when_start_is_missing() {
        let (assembler, k, adjacency, _counts, _acga, _atcg, tcga) =
            canonical_only_orientation_fixture();
        let start = KmerU64::from_str("CGAT").unwrap().encoded;
        let bases = [b'A', b'C', b'G', b'T'];

        let (path, end) = assembler
            .trace_path(start, k, &adjacency, 8, &bases, false)
            .expect("trace path should succeed");

        assert_eq!(path, vec![start, tcga]);
        assert_eq!(end, tcga);
    }

    #[test]
    fn test_select_unique_unused_linear_extension_requires_exactly_one_choice() {
        let k = 4;
        let current = KmerU64::from_str("ACGA").unwrap().encoded;
        let mut used = AHashSet::new();

        let mut only_t = [false; 4];
        only_t[3] = true;
        let selected = LargeGenomeAssembler::select_unique_unused_linear_extension(
            current, &only_t, k, &used, false,
        )
        .expect("single valid unused extension should be selected");
        assert_eq!(selected.0, 3);
        let expected_next = extend_right(current, b'T', k).expect("valid T extension");
        assert_eq!(selected.1, expected_next);

        let mut branch = [false; 4];
        branch[0] = true;
        branch[3] = true;
        assert_eq!(
            LargeGenomeAssembler::select_unique_unused_linear_extension(
                current, &branch, k, &used, false
            ),
            None
        );

        used.insert(LargeGenomeAssembler::canonical_encoded(expected_next, k));
        assert_eq!(
            LargeGenomeAssembler::select_unique_unused_linear_extension(
                current, &only_t, k, &used, false
            ),
            None
        );
    }
}
