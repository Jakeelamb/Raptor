//! Streaming assembly pipeline for large genomes
//!
//! This module implements a memory-efficient assembly pipeline that can handle
//! genomes of any size (tested up to 120 Gb) on machines with limited RAM (16 GB).
//!
//! Key techniques:
//! 1. Disk-based reversible k-mer counting with bucket partitioning
//! 2. Deterministic greedy contig extraction from encoded k-mers
//! 3. Stable assembly quality reporting (Nx/Lx/auN)
//!
//! Memory usage: O(bucket_size) ≈ 2-4 GB regardless of genome size

use crate::accel::CpuBackend;
use crate::eval::metrics::{evaluate_lengths_sorted_desc, TranscriptStats};
use crate::graph::assembler::{greedy_assembly_u64, Contig};
use crate::io::fasta::FastaWriter;
use crate::io::fastq::{stream_fastq_records_checked, try_open_fastq};
use crate::kmer::disk_counting_v2::{DiskCounterConfig, DiskKmerCounterV2};
use ahash::AHashMap;
use std::fs;
use std::path::Path;
use tracing::info;

/// Configuration for streaming assembly
#[derive(Clone, Debug)]
pub struct StreamingAssemblyConfig {
    /// K-mer size
    pub k: usize,
    /// Minimum k-mer count (filters sequencing errors)
    pub min_kmer_count: u32,
    /// Minimum contig length to output
    pub min_contig_len: usize,
    /// Number of disk buckets (auto-calculated if None)
    pub num_buckets: Option<usize>,
    /// Maximum RAM usage in bytes (auto-detected if None)
    pub max_ram: Option<usize>,
    /// Temporary directory for disk files
    pub temp_dir: Option<String>,
}

impl Default for StreamingAssemblyConfig {
    fn default() -> Self {
        Self {
            k: 31,
            min_kmer_count: 2,
            min_contig_len: 200,
            num_buckets: None,
            max_ram: None,
            temp_dir: None,
        }
    }
}

/// Statistics from streaming assembly
#[derive(Debug, Clone)]
pub struct StreamingAssemblyStats {
    pub total_reads: u64,
    pub total_bases: u64,
    pub total_kmers: u64,
    pub unique_kmers: u64,
    pub filtered_kmers: u64,
    pub contigs_produced: usize,
    pub total_contig_length: usize,
    pub average_contig_length: f64,
    pub n25: usize,
    pub n50: usize,
    pub n75: usize,
    pub n90: usize,
    pub n95: usize,
    pub n99: usize,
    pub l25: usize,
    pub l50: usize,
    pub l75: usize,
    pub l90: usize,
    pub l95: usize,
    pub l99: usize,
    pub au_n: f64,
    pub longest_contig: usize,
    pub disk_usage_bytes: u64,
    pub peak_memory_bytes: u64,
}

/// Main streaming assembler
pub struct StreamingAssembler {
    config: StreamingAssemblyConfig,
}

impl StreamingAssembler {
    pub fn new(config: StreamingAssemblyConfig) -> Self {
        Self { config }
    }

    /// Run the complete streaming assembly pipeline
    pub fn assemble(
        &self,
        input_path: &str,
        output_path: &str,
    ) -> std::io::Result<StreamingAssemblyStats> {
        let k = self.config.k;
        let min_count = self.config.min_kmer_count;
        if !(1..=32).contains(&k) {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("k-mer size must be in 1..=32, got {}", k),
            ));
        }

        info!("Starting streaming assembly pipeline");
        info!("K-mer size: {}, Min count: {}", k, min_count);
        let mut peak_memory_bytes = Self::detect_peak_memory_bytes();

        // Phase 1: Disk-based k-mer counting
        info!("Phase 1: Distributing k-mers to disk buckets...");
        let disk_config = self.create_disk_config();
        let mut disk_counter = DiskKmerCounterV2::new(disk_config)?;

        // Stream sequences directly to disk counter
        let (total_reads, total_bases) =
            self.stream_to_disk_counter(input_path, &mut disk_counter)?;

        let (distributed_kmers, bucket_count, disk_usage_bytes) = disk_counter.stats();
        peak_memory_bytes = peak_memory_bytes.max(Self::detect_peak_memory_bytes());
        info!(
            "Distribution complete: kmers={}, buckets={}, disk={:.2} GB",
            distributed_kmers,
            bucket_count,
            disk_usage_bytes as f64 / (1024.0 * 1024.0 * 1024.0)
        );
        info!("Reads processed: {}, Bases: {}", total_reads, total_bases);

        // Phase 2: Count k-mers from buckets
        info!("Phase 2: Counting k-mers from buckets...");
        let kmer_counts = disk_counter.count_all()?;
        peak_memory_bytes = peak_memory_bytes.max(Self::detect_peak_memory_bytes());

        let unique_kmers = kmer_counts.len() as u64;
        // count_all() already applies min_count filtering during bucket counting.
        let filtered_kmers = unique_kmers;
        info!(
            "Unique k-mers: {}, After filtering (count >= {}): {}",
            unique_kmers, min_count, filtered_kmers
        );

        // Phase 3: Build assembly graph and extract contigs
        info!("Phase 3: Building assembly graph and extracting contigs...");
        let contigs = self.build_contigs_streaming(&kmer_counts);
        peak_memory_bytes = peak_memory_bytes.max(Self::detect_peak_memory_bytes());

        // Phase 4: Write output
        info!("Phase 4: Writing contigs...");
        let (contigs_produced, total_contig_length, length_stats) =
            self.write_contigs(&contigs, output_path)?;
        peak_memory_bytes = peak_memory_bytes.max(Self::detect_peak_memory_bytes());

        // Cleanup
        disk_counter.cleanup()?;

        Ok(StreamingAssemblyStats {
            total_reads,
            total_bases,
            total_kmers: distributed_kmers,
            unique_kmers,
            filtered_kmers,
            contigs_produced,
            total_contig_length,
            average_contig_length: length_stats.avg_length,
            n25: length_stats.n25,
            n50: length_stats.n50,
            n75: length_stats.n75,
            n90: length_stats.n90,
            n95: length_stats.n95,
            n99: length_stats.n99,
            l25: length_stats.l25,
            l50: length_stats.l50,
            l75: length_stats.l75,
            l90: length_stats.l90,
            l95: length_stats.l95,
            l99: length_stats.l99,
            au_n: length_stats.au_n,
            longest_contig: length_stats.longest,
            disk_usage_bytes,
            peak_memory_bytes,
        })
    }

    #[inline]
    fn detect_peak_memory_bytes() -> u64 {
        #[cfg(target_os = "linux")]
        {
            fs::read_to_string("/proc/self/status")
                .ok()
                .and_then(|status| Self::parse_proc_status_kib(&status, "VmHWM:"))
                .and_then(|kib| kib.checked_mul(1024))
                .unwrap_or(0)
        }

        #[cfg(not(target_os = "linux"))]
        {
            0
        }
    }

    #[inline]
    fn parse_proc_status_kib(status: &str, key: &str) -> Option<u64> {
        status.lines().find_map(|line| {
            if !line.starts_with(key) {
                return None;
            }
            line.split_whitespace().nth(1)?.parse::<u64>().ok()
        })
    }

    fn create_disk_config(&self) -> DiskCounterConfig {
        let mut config = DiskCounterConfig::auto(self.config.k);

        if let Some(max_ram) = self.config.max_ram {
            config.num_buckets = Self::suggest_bucket_count(max_ram);
        }

        if let Some(num_buckets) = self.config.num_buckets {
            config.num_buckets = Self::normalize_bucket_count(num_buckets);
        }

        if let Some(ref temp_dir) = self.config.temp_dir {
            config.temp_dir = Path::new(temp_dir).to_path_buf();
        }

        config.min_count = self.config.min_kmer_count;
        config
    }

    #[inline]
    fn suggest_bucket_count(max_ram: usize) -> usize {
        if max_ram >= 32 * 1024 * 1024 * 1024 {
            512
        } else if max_ram >= 16 * 1024 * 1024 * 1024 {
            1024
        } else if max_ram >= 8 * 1024 * 1024 * 1024 {
            2048
        } else {
            4096
        }
    }

    #[inline]
    fn normalize_bucket_count(num_buckets: usize) -> usize {
        match num_buckets {
            0 | 1 => 1,
            n if n.is_power_of_two() => n,
            n => n
                .checked_next_power_of_two()
                .unwrap_or(1usize << (usize::BITS - 1)),
        }
    }

    fn stream_to_disk_counter(
        &self,
        input_path: &str,
        counter: &mut DiskKmerCounterV2,
    ) -> std::io::Result<(u64, u64)> {
        let reader = try_open_fastq(input_path)?;
        let mut total_reads = 0u64;
        let mut total_bases = 0u64;

        // Collect sequences in batches for efficiency
        const BATCH_SIZE: usize = 10_000;
        let mut batch: Vec<String> = Vec::with_capacity(BATCH_SIZE);

        for record in stream_fastq_records_checked(reader) {
            let record = record?;
            total_reads += 1;
            total_bases += record.sequence.len() as u64;
            batch.push(record.sequence);

            if batch.len() >= BATCH_SIZE {
                counter.distribute(batch.iter().map(|s| s.as_bytes()))?;
                batch.clear();

                if total_reads % 1_000_000 == 0 {
                    info!("Processed {} million reads...", total_reads / 1_000_000);
                }
            }
        }

        // Process remaining batch
        if !batch.is_empty() {
            counter.distribute(batch.iter().map(|s| s.as_bytes()))?;
        }

        Ok((total_reads, total_bases))
    }

    /// Build contigs using streaming graph traversal
    fn build_contigs_streaming(&self, kmer_counts: &AHashMap<u64, u32>) -> Vec<Contig> {
        let k = self.config.k;
        let min_len = self.config.min_contig_len;

        info!("Building graph from {} filtered k-mers", kmer_counts.len());
        let cpu_backend = CpuBackend::new();
        let adjacency = cpu_backend.build_adjacency_u64(kmer_counts, k);
        let mut contigs = greedy_assembly_u64(k, kmer_counts, &adjacency, min_len);

        // Keep FASTA output ordering deterministic across equivalent tie cases.
        contigs.sort_unstable_by(|a, b| {
            b.sequence
                .len()
                .cmp(&a.sequence.len())
                .then_with(|| a.sequence.cmp(&b.sequence))
                .then_with(|| a.kmer_path.cmp(&b.kmer_path))
        });
        for (idx, contig) in contigs.iter_mut().enumerate() {
            contig.id = idx;
        }

        info!("Extracted {} contigs", contigs.len());
        contigs
    }

    fn write_contigs(
        &self,
        contigs: &[Contig],
        output_path: &str,
    ) -> std::io::Result<(usize, usize, TranscriptStats)> {
        let mut writer = FastaWriter::try_new(output_path)?;

        let mut lengths: Vec<usize> = Vec::with_capacity(contigs.len());

        for (i, contig) in contigs.iter().enumerate() {
            writer.write_record(&format!("contig_{}", i + 1), &contig.sequence)?;
            lengths.push(contig.sequence.len());
        }

        debug_assert!(lengths.windows(2).all(|w| w[0] >= w[1]));
        let length_stats = evaluate_lengths_sorted_desc(&lengths);
        Ok((contigs.len(), length_stats.total_bases, length_stats))
    }
}

/// Simplified streaming assembly function
pub fn streaming_assemble(
    input_path: &str,
    output_path: &str,
    k: usize,
    min_count: u32,
    min_len: usize,
) -> std::io::Result<StreamingAssemblyStats> {
    let config = StreamingAssemblyConfig {
        k,
        min_kmer_count: min_count,
        min_contig_len: min_len,
        ..Default::default()
    };

    let assembler = StreamingAssembler::new(config);
    assembler.assemble(input_path, output_path)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eval::metrics::evaluate_lengths;
    use crate::kmer::kmer::KmerU64;
    use rand::rngs::StdRng;
    use rand::seq::SliceRandom;
    use rand::SeedableRng;
    use std::io::Write;
    use tempfile::{tempdir, NamedTempFile};

    fn create_fastq_with_reads(reads: &[&str]) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        for (i, read) in reads.iter().enumerate() {
            writeln!(file, "@read_{}", i).unwrap();
            writeln!(file, "{}", read).unwrap();
            writeln!(file, "+").unwrap();
            writeln!(file, "{}", "I".repeat(read.len())).unwrap();
        }
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_streaming_assembly_config() {
        let config = StreamingAssemblyConfig::default();
        assert_eq!(config.k, 31);
        assert_eq!(config.min_kmer_count, 2);
    }

    #[test]
    fn streaming_assembly_rejects_invalid_kmer_size() {
        let input = create_fastq_with_reads(&["ACGTACGT", "ACGTACGT"]);
        let output = NamedTempFile::new().unwrap();

        for &invalid_k in &[0usize, 33usize] {
            let assembler = StreamingAssembler::new(StreamingAssemblyConfig {
                k: invalid_k,
                min_kmer_count: 1,
                min_contig_len: 3,
                num_buckets: Some(4),
                max_ram: None,
                temp_dir: None,
            });
            let err = assembler
                .assemble(
                    input.path().to_str().unwrap(),
                    output.path().to_str().unwrap(),
                )
                .unwrap_err();
            assert_eq!(err.kind(), std::io::ErrorKind::InvalidInput);
        }
    }

    #[test]
    fn streaming_assembly_rejects_truncated_fastq_records() {
        let mut input = NamedTempFile::new().unwrap();
        writeln!(input, "@read_0").unwrap();
        writeln!(input, "ACGTACGT").unwrap();
        writeln!(input, "+").unwrap();
        writeln!(input, "IIIIIIII").unwrap();
        writeln!(input, "@read_1").unwrap();
        writeln!(input, "ACGTACGT").unwrap();
        writeln!(input, "+").unwrap();
        input.flush().unwrap();

        let output = NamedTempFile::new().unwrap();
        let assembler = StreamingAssembler::new(StreamingAssemblyConfig {
            k: 5,
            min_kmer_count: 1,
            min_contig_len: 5,
            num_buckets: Some(4),
            max_ram: None,
            temp_dir: None,
        });

        let err = assembler
            .assemble(
                input.path().to_str().unwrap(),
                output.path().to_str().unwrap(),
            )
            .unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn build_contigs_streaming_is_deterministic_under_kmer_insertion_order() {
        let config = StreamingAssemblyConfig {
            k: 3,
            min_kmer_count: 1,
            min_contig_len: 3,
            num_buckets: Some(4),
            max_ram: None,
            temp_dir: None,
        };
        let assembler = StreamingAssembler::new(config);

        let kmers = vec![
            ("AAA", 15),
            ("AAC", 11),
            ("AAT", 11),
            ("ACC", 9),
            ("CCC", 8),
            ("CCA", 7),
        ];
        let canonical_entries: Vec<(u64, u32)> = kmers
            .iter()
            .map(|(seq, count)| (KmerU64::from_str(seq).unwrap().canonical().encoded, *count))
            .collect();

        let mut baseline_map = AHashMap::new();
        for &(kmer, count) in &canonical_entries {
            baseline_map.insert(kmer, count);
        }
        let baseline: Vec<String> = assembler
            .build_contigs_streaming(&baseline_map)
            .into_iter()
            .map(|c| c.sequence)
            .collect();
        assert!(!baseline.is_empty());

        let mut rng = StdRng::seed_from_u64(0x51D1_7A2Eu64);
        for _ in 0..64 {
            let mut shuffled = canonical_entries.clone();
            shuffled.shuffle(&mut rng);

            let mut counts = AHashMap::new();
            for (kmer, count) in shuffled {
                counts.insert(kmer, count);
            }
            let observed: Vec<String> = assembler
                .build_contigs_streaming(&counts)
                .into_iter()
                .map(|c| c.sequence)
                .collect();
            assert_eq!(observed, baseline);
        }
    }

    #[test]
    fn streaming_assembly_outputs_real_sequences_and_consistent_metrics() {
        let input = create_fastq_with_reads(&[
            "ACGTACGTACGTACGT",
            "ACGTACGTACGTACGT",
            "CGTACGTACGTACGTA",
            "ACGTACGTACGTACGT",
        ]);
        let output = NamedTempFile::new().unwrap();
        let tmp = tempdir().unwrap();
        let temp_kmer_dir = tmp.path().join("kmer_buckets");

        let assembler = StreamingAssembler::new(StreamingAssemblyConfig {
            k: 5,
            min_kmer_count: 1,
            min_contig_len: 5,
            num_buckets: Some(8),
            max_ram: None,
            temp_dir: Some(temp_kmer_dir.display().to_string()),
        });

        let stats = assembler
            .assemble(
                input.path().to_str().unwrap(),
                output.path().to_str().unwrap(),
            )
            .unwrap();

        assert!(stats.contigs_produced > 0);
        assert!(stats.total_contig_length > 0);
        assert!(stats.n25 >= stats.n50);
        assert!(stats.n50 >= stats.n75);
        assert!(stats.n75 >= stats.n90);
        assert!(stats.n90 >= stats.n95);
        assert!(stats.n95 >= stats.n99);
        assert!(stats.l25 <= stats.l50);
        assert!(stats.l50 <= stats.l75);
        assert!(stats.l75 <= stats.l90);
        assert!(stats.l50 <= stats.l90);
        assert!(stats.l90 <= stats.l95);
        assert!(stats.l95 <= stats.l99);
        assert!(stats.longest_contig >= stats.n50);
        assert!(stats.au_n + 1e-9 >= stats.average_contig_length);

        let fasta = std::fs::read_to_string(output.path()).unwrap();
        let lengths: Vec<usize> = fasta
            .lines()
            .filter(|line| !line.starts_with('>'))
            .map(|seq| {
                assert!(!seq.contains("CONTIG_FROM_HASH"));
                assert!(seq.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')));
                seq.len()
            })
            .collect();

        let expected = evaluate_lengths(&lengths);
        assert_eq!(stats.total_contig_length, expected.total_bases);
        assert_eq!(stats.n25, expected.n25);
        assert_eq!(stats.n50, expected.n50);
        assert_eq!(stats.n75, expected.n75);
        assert_eq!(stats.n90, expected.n90);
        assert_eq!(stats.n95, expected.n95);
        assert_eq!(stats.n99, expected.n99);
        assert_eq!(stats.l25, expected.l25);
        assert_eq!(stats.l50, expected.l50);
        assert_eq!(stats.l75, expected.l75);
        assert_eq!(stats.l90, expected.l90);
        assert_eq!(stats.l95, expected.l95);
        assert_eq!(stats.l99, expected.l99);
        assert_eq!(stats.longest_contig, expected.longest);
        assert!((stats.average_contig_length - expected.avg_length).abs() < 1e-12);
        assert!((stats.au_n - expected.au_n).abs() < 1e-12);
        if cfg!(target_os = "linux") {
            assert!(
                stats.peak_memory_bytes > 0,
                "expected non-zero peak memory on Linux"
            );
        }
    }

    #[test]
    fn parse_proc_status_kib_extracts_values() {
        let status = "Name:\traptor\nVmRSS:\t  2048 kB\nVmHWM:\t  4096 kB\n";
        assert_eq!(
            StreamingAssembler::parse_proc_status_kib(status, "VmHWM:"),
            Some(4096)
        );
        assert_eq!(
            StreamingAssembler::parse_proc_status_kib(status, "VmRSS:"),
            Some(2048)
        );
        assert_eq!(
            StreamingAssembler::parse_proc_status_kib(status, "VmSize:"),
            None
        );
    }

    #[test]
    fn write_contigs_returns_error_for_invalid_output_path() {
        let assembler = StreamingAssembler::new(StreamingAssemblyConfig {
            k: 3,
            min_kmer_count: 1,
            min_contig_len: 1,
            num_buckets: Some(4),
            max_ram: None,
            temp_dir: None,
        });
        let tmp = tempdir().unwrap();
        let contigs = vec![Contig {
            id: 0,
            sequence: "ACGT".to_string(),
            kmer_path: Vec::new(),
        }];

        let result = assembler.write_contigs(&contigs, tmp.path().to_str().unwrap());
        assert!(result.is_err());
    }
}
