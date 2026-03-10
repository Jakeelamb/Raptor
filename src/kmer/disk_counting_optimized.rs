//! Optimized Disk-based k-mer counting with LZ4 compression
//!
//! Improvements over disk_counting_v2:
//! - LZ4 compression for bucket files (2-4x smaller)
//! - Memory-mapped I/O for faster reading
//! - Parallel sequence distribution
//! - Streaming compression/decompression

use crate::kmer::kmer::KmerU64;
use ahash::AHashMap;
use lz4_flex::{compress_prepend_size, decompress_size_prepended};
use memmap2::Mmap;
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

/// Configuration for optimized disk-based k-mer counting
#[derive(Clone, Debug)]
pub struct OptimizedDiskConfig {
    pub k: usize,
    pub num_buckets: usize,
    pub min_count: u32,
    pub temp_dir: PathBuf,
    pub compression_enabled: bool,
    pub parallel_distribution: bool,
    pub chunk_size: usize,
}

impl Default for OptimizedDiskConfig {
    fn default() -> Self {
        Self {
            k: 31,
            num_buckets: 1024,
            min_count: 2,
            temp_dir: std::env::temp_dir().join("raptor_kmer_opt"),
            compression_enabled: true,
            parallel_distribution: true,
            chunk_size: 100_000, // Sequences per parallel chunk
        }
    }
}

impl OptimizedDiskConfig {
    pub fn auto(k: usize) -> Self {
        let ram = Self::detect_ram();
        let num_buckets = if ram >= 32 * 1024 * 1024 * 1024 {
            512
        } else if ram >= 16 * 1024 * 1024 * 1024 {
            1024
        } else if ram >= 8 * 1024 * 1024 * 1024 {
            2048
        } else {
            4096
        };

        Self {
            k,
            num_buckets,
            ..Default::default()
        }
    }

    fn detect_ram() -> usize {
        if let Ok(contents) = fs::read_to_string("/proc/meminfo") {
            for line in contents.lines() {
                if line.starts_with("MemTotal:") {
                    if let Some(kb) = line.split_whitespace().nth(1) {
                        if let Ok(kb) = kb.parse::<usize>() {
                            return kb * 1024;
                        }
                    }
                }
            }
        }
        16 * 1024 * 1024 * 1024
    }
}

/// Statistics from disk counting
#[derive(Debug, Clone, Default)]
pub struct DiskCountingStats {
    pub total_kmers: u64,
    pub unique_kmers: u64,
    pub disk_bytes_raw: u64,
    pub disk_bytes_compressed: u64,
    pub compression_ratio: f64,
}

/// Optimized disk-based k-mer counter
pub struct OptimizedDiskCounter {
    config: OptimizedDiskConfig,
    bucket_paths: Vec<PathBuf>,
    bucket_counts: Vec<u64>,
    total_kmers: u64,
    compressed_sizes: Vec<u64>,
}

#[inline]
fn encode_base_2bit(base: u8) -> Option<u64> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

impl OptimizedDiskCounter {
    pub fn new(config: OptimizedDiskConfig) -> std::io::Result<Self> {
        if config.k == 0 || config.k > 32 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!(
                    "invalid k-mer size {}: supported range is 1..=32 for u64 encoding",
                    config.k
                ),
            ));
        }
        if config.num_buckets == 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "num_buckets must be greater than zero",
            ));
        }
        if config.chunk_size == 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "chunk_size must be greater than zero",
            ));
        }

        fs::create_dir_all(&config.temp_dir)?;
        Ok(Self {
            config,
            bucket_paths: Vec::new(),
            bucket_counts: Vec::new(),
            total_kmers: 0,
            compressed_sizes: Vec::new(),
        })
    }

    /// Distribute k-mers to buckets with optional parallelism
    pub fn distribute<I, S>(&mut self, sequences: I) -> std::io::Result<()>
    where
        I: Iterator<Item = S>,
        S: AsRef<[u8]>,
    {
        let k = self.config.k;
        let num_buckets = self.config.num_buckets;
        let mask = (num_buckets - 1) as u64;

        if self.config.parallel_distribution {
            self.distribute_parallel(sequences, k, num_buckets, mask)
        } else {
            self.distribute_sequential(sequences, k, num_buckets, mask)
        }
    }

    #[inline]
    fn bucket_sequence(
        seq: &[u8],
        k: usize,
        num_buckets: usize,
        bucket_mask: u64,
        use_mask_bucket: bool,
        rolling_mask: u64,
        bucket_data: &mut [Vec<u64>],
    ) -> u64 {
        if seq.len() < k {
            return 0;
        }

        let mut rolling = 0u64;
        let mut valid_run = 0usize;
        let mut total = 0u64;

        for &base in seq {
            if let Some(base_bits) = encode_base_2bit(base) {
                rolling = ((rolling << 2) | base_bits) & rolling_mask;
                valid_run += 1;

                if valid_run >= k {
                    let canonical = KmerU64 {
                        encoded: rolling,
                        len: k as u8,
                    }
                    .canonical()
                    .encoded;
                    let bucket_id = if use_mask_bucket {
                        (canonical & bucket_mask) as usize
                    } else {
                        (canonical % num_buckets as u64) as usize
                    };
                    bucket_data[bucket_id].push(canonical);
                    total += 1;
                }
            } else {
                rolling = 0;
                valid_run = 0;
            }
        }

        total
    }

    fn distribute_sequential<I, S>(
        &mut self,
        sequences: I,
        k: usize,
        num_buckets: usize,
        bucket_mask: u64,
    ) -> std::io::Result<()>
    where
        I: Iterator<Item = S>,
        S: AsRef<[u8]>,
    {
        // Accumulate k-mers in memory per bucket
        let mut bucket_data: Vec<Vec<u64>> = (0..num_buckets).map(|_| Vec::new()).collect();
        let use_mask_bucket = num_buckets.is_power_of_two();
        let rolling_mask = if k >= 32 {
            u64::MAX
        } else {
            (1u64 << (k * 2)) - 1
        };
        let mut total = 0u64;

        for seq in sequences {
            let seq_total = Self::bucket_sequence(
                seq.as_ref(),
                k,
                num_buckets,
                bucket_mask,
                use_mask_bucket,
                rolling_mask,
                &mut bucket_data,
            );
            total = total.checked_add(seq_total).ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "total k-mer count overflowed u64 during distribution",
                )
            })?;
        }

        // Write buckets with optional compression
        self.write_buckets(bucket_data)?;
        self.total_kmers = total;

        Ok(())
    }

    fn count_parallel_chunk(
        sequences: &[Vec<u8>],
        k: usize,
        num_buckets: usize,
        bucket_mask: u64,
        use_mask_bucket: bool,
        rolling_mask: u64,
    ) -> std::io::Result<(Vec<Vec<u64>>, u64)> {
        if sequences.is_empty() {
            return Ok(((0..num_buckets).map(|_| Vec::new()).collect(), 0));
        }

        let worker_chunks = rayon::current_num_threads().max(1);
        let sub_chunk_size = (sequences.len() / worker_chunks).max(1);

        let partial_results: Vec<(Vec<Vec<u64>>, u64)> = sequences
            .par_chunks(sub_chunk_size)
            .map(|chunk| -> std::io::Result<(Vec<Vec<u64>>, u64)> {
                let mut bucket_data: Vec<Vec<u64>> = (0..num_buckets).map(|_| Vec::new()).collect();
                let mut total = 0u64;
                for seq in chunk {
                    let seq_total = Self::bucket_sequence(
                        seq,
                        k,
                        num_buckets,
                        bucket_mask,
                        use_mask_bucket,
                        rolling_mask,
                        &mut bucket_data,
                    );
                    total = total.checked_add(seq_total).ok_or_else(|| {
                        std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            "chunk k-mer count overflowed u64 during parallel distribution",
                        )
                    })?;
                }
                Ok((bucket_data, total))
            })
            .collect::<std::io::Result<Vec<_>>>()?;

        let mut merged: Vec<Vec<u64>> = (0..num_buckets).map(|_| Vec::new()).collect();
        let mut total = 0u64;
        for (partial_buckets, partial_total) in partial_results {
            total = total.checked_add(partial_total).ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "chunk merge k-mer count overflowed u64 during parallel distribution",
                )
            })?;
            for (bucket_idx, bucket) in partial_buckets.into_iter().enumerate() {
                merged[bucket_idx].extend(bucket);
            }
        }

        Ok((merged, total))
    }

    fn distribute_parallel<I, S>(
        &mut self,
        sequences: I,
        k: usize,
        num_buckets: usize,
        bucket_mask: u64,
    ) -> std::io::Result<()>
    where
        I: Iterator<Item = S>,
        S: AsRef<[u8]>,
    {
        let chunk_size = self.config.chunk_size;
        let use_mask_bucket = num_buckets.is_power_of_two();
        let rolling_mask = if k >= 32 {
            u64::MAX
        } else {
            (1u64 << (k * 2)) - 1
        };

        let mut merged: Vec<Vec<u64>> = (0..num_buckets).map(|_| Vec::new()).collect();
        let mut total = 0u64;
        let mut buffered: Vec<Vec<u8>> = Vec::with_capacity(chunk_size.max(1));

        for sequence in sequences {
            buffered.push(sequence.as_ref().to_vec());
            if buffered.len() >= chunk_size {
                let (chunk_buckets, chunk_total) = Self::count_parallel_chunk(
                    &buffered,
                    k,
                    num_buckets,
                    bucket_mask,
                    use_mask_bucket,
                    rolling_mask,
                )?;
                total = total.checked_add(chunk_total).ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        "total k-mer count overflowed u64 during parallel distribution",
                    )
                })?;
                for (i, bucket) in chunk_buckets.into_iter().enumerate() {
                    merged[i].extend(bucket);
                }
                buffered.clear();
            }
        }
        if !buffered.is_empty() {
            let (chunk_buckets, chunk_total) = Self::count_parallel_chunk(
                &buffered,
                k,
                num_buckets,
                bucket_mask,
                use_mask_bucket,
                rolling_mask,
            )?;
            total = total.checked_add(chunk_total).ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "total k-mer count overflowed u64 during parallel distribution",
                )
            })?;
            for (i, bucket) in chunk_buckets.into_iter().enumerate() {
                merged[i].extend(bucket);
            }
        }

        // Write buckets
        self.write_buckets(merged)?;
        self.total_kmers = total;

        Ok(())
    }

    fn write_buckets(&mut self, bucket_data: Vec<Vec<u64>>) -> std::io::Result<()> {
        let compress = self.config.compression_enabled;
        let temp_dir = &self.config.temp_dir;

        let results: Vec<(PathBuf, u64, u64)> = bucket_data
            .into_par_iter()
            .enumerate()
            .map(|(i, data)| -> std::io::Result<(PathBuf, u64, u64)> {
                let count = data.len() as u64;
                let ext = if compress { "lz4" } else { "bin" };
                let path = temp_dir.join(format!("bucket_{:05}.{}", i, ext));

                if count == 0 {
                    // Create empty file
                    File::create(&path)?;
                    return Ok((path, 0, 0));
                }

                // Convert to bytes
                let raw_bytes: Vec<u8> = data.iter().flat_map(|&v| v.to_le_bytes()).collect();

                let compressed_size = if compress {
                    let compressed = compress_prepend_size(&raw_bytes);
                    let mut file = BufWriter::new(File::create(&path)?);
                    file.write_all(&compressed)?;
                    file.flush()?;
                    compressed.len() as u64
                } else {
                    let mut file = BufWriter::new(File::create(&path)?);
                    file.write_all(&raw_bytes)?;
                    file.flush()?;
                    raw_bytes.len() as u64
                };

                Ok((path, count, compressed_size))
            })
            .collect::<std::io::Result<Vec<_>>>()?;

        self.bucket_paths = results.iter().map(|(p, _, _)| p.clone()).collect();
        self.bucket_counts = results.iter().map(|(_, c, _)| *c).collect();
        self.compressed_sizes = results.iter().map(|(_, _, s)| *s).collect();

        Ok(())
    }

    /// Count k-mers using memory-mapped I/O
    pub fn count_all(&self) -> std::io::Result<AHashMap<u64, u32>> {
        let min_count = self.config.min_count;
        let compress = self.config.compression_enabled;

        // Process buckets in parallel
        let bucket_results: Vec<AHashMap<u64, u32>> = self
            .bucket_paths
            .par_iter()
            .enumerate()
            .filter(|(i, _)| self.bucket_counts[*i] > 0)
            .map(|(_, path)| {
                if compress {
                    Self::count_bucket_compressed(path, min_count)
                } else {
                    Self::count_bucket_mmap(path, min_count)
                }
            })
            .collect::<std::io::Result<Vec<_>>>()?;

        // Merge results
        let total_unique: usize = bucket_results.iter().map(|m| m.len()).sum();
        let mut merged = AHashMap::with_capacity(total_unique);
        for bucket in bucket_results {
            merged.extend(bucket);
        }

        Ok(merged)
    }

    fn count_bucket_mmap(path: &Path, min_count: u32) -> std::io::Result<AHashMap<u64, u32>> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };

        if mmap.is_empty() {
            return Ok(AHashMap::new());
        }
        if mmap.len() % 8 != 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "bucket file '{}' has invalid byte length {} (not divisible by 8)",
                    path.display(),
                    mmap.len()
                ),
            ));
        }

        let num_kmers = mmap.len() / 8;
        let mut kmers: Vec<u64> = Vec::with_capacity(num_kmers);

        for chunk in mmap.chunks_exact(8) {
            let bytes: [u8; 8] = chunk.try_into().unwrap();
            kmers.push(u64::from_le_bytes(bytes));
        }

        Self::count_sorted_kmers(kmers, min_count)
    }

    fn count_bucket_compressed(path: &Path, min_count: u32) -> std::io::Result<AHashMap<u64, u32>> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };

        if mmap.is_empty() {
            return Ok(AHashMap::new());
        }

        // Decompress
        let decompressed = decompress_size_prepended(&mmap)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        if decompressed.len() % 8 != 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "compressed bucket '{}' decoded to invalid byte length {} (not divisible by 8)",
                    path.display(),
                    decompressed.len()
                ),
            ));
        }

        let num_kmers = decompressed.len() / 8;
        let mut kmers: Vec<u64> = Vec::with_capacity(num_kmers);

        for chunk in decompressed.chunks_exact(8) {
            let bytes: [u8; 8] = chunk.try_into().unwrap();
            kmers.push(u64::from_le_bytes(bytes));
        }

        Self::count_sorted_kmers(kmers, min_count)
    }

    fn count_sorted_kmers(
        mut kmers: Vec<u64>,
        min_count: u32,
    ) -> std::io::Result<AHashMap<u64, u32>> {
        if kmers.is_empty() {
            return Ok(AHashMap::new());
        }

        // Sort for counting
        kmers.sort_unstable();

        // Count consecutive
        let mut counts = AHashMap::with_capacity(kmers.len() / 5);
        let mut current = kmers[0];
        let mut count = 1u32;

        for &kmer in &kmers[1..] {
            if kmer == current {
                count = count.checked_add(1).ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "k-mer count overflow for encoded value {} (exceeds u32)",
                            current
                        ),
                    )
                })?;
            } else {
                if count >= min_count {
                    counts.insert(current, count);
                }
                current = kmer;
                count = 1;
            }
        }

        if count >= min_count {
            counts.insert(current, count);
        }

        Ok(counts)
    }

    /// Get detailed statistics
    pub fn stats(&self) -> DiskCountingStats {
        let raw_bytes = self.total_kmers * 8;
        let compressed_bytes: u64 = self.compressed_sizes.iter().sum();

        DiskCountingStats {
            total_kmers: self.total_kmers,
            unique_kmers: 0, // Filled after counting
            disk_bytes_raw: raw_bytes,
            disk_bytes_compressed: compressed_bytes,
            compression_ratio: if compressed_bytes > 0 {
                raw_bytes as f64 / compressed_bytes as f64
            } else {
                1.0
            },
        }
    }

    /// Get basic stats tuple for compatibility
    pub fn stats_tuple(&self) -> (u64, usize, u64) {
        let disk_bytes: u64 = self.compressed_sizes.iter().sum();
        (self.total_kmers, self.bucket_paths.len(), disk_bytes)
    }

    /// Cleanup temporary files
    pub fn cleanup(&self) -> std::io::Result<()> {
        for path in &self.bucket_paths {
            let _ = fs::remove_file(path);
        }
        let _ = fs::remove_dir(&self.config.temp_dir);
        Ok(())
    }
}

impl Drop for OptimizedDiskCounter {
    fn drop(&mut self) {
        let _ = self.cleanup();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn run_counter(config: OptimizedDiskConfig, seqs: &[&[u8]]) -> AHashMap<u64, u32> {
        let mut counter = OptimizedDiskCounter::new(config).expect("counter should construct");
        counter
            .distribute(seqs.iter().copied())
            .expect("distribution should succeed");
        counter.count_all().expect("counting should succeed")
    }

    #[test]
    fn test_optimized_counter_compressed() {
        let config = OptimizedDiskConfig {
            k: 11,
            num_buckets: 4,
            min_count: 1,
            temp_dir: std::env::temp_dir().join("raptor_test_opt_compress"),
            compression_enabled: true,
            parallel_distribution: false,
            ..Default::default()
        };

        let mut counter = OptimizedDiskCounter::new(config).unwrap();

        let seqs = vec!["ACGTACGTACGTACGT", "ACGTACGTACGTACGT", "TGCATGCATGCATGCA"];

        counter
            .distribute(seqs.iter().map(|s| s.as_bytes()))
            .unwrap();

        let stats = counter.stats();
        assert!(stats.total_kmers > 0);
        assert!(stats.compression_ratio >= 1.0, "Compression should help");

        let counts = counter.count_all().unwrap();
        assert!(!counts.is_empty());
    }

    #[test]
    fn test_optimized_counter_uncompressed() {
        let config = OptimizedDiskConfig {
            k: 11,
            num_buckets: 4,
            min_count: 1,
            temp_dir: std::env::temp_dir().join("raptor_test_opt_raw"),
            compression_enabled: false,
            parallel_distribution: false,
            ..Default::default()
        };

        let mut counter = OptimizedDiskCounter::new(config).unwrap();

        let seqs = vec!["ACGTACGTACGTACGT", "ACGTACGTACGTACGT", "TGCATGCATGCATGCA"];

        counter
            .distribute(seqs.iter().map(|s| s.as_bytes()))
            .unwrap();

        let counts = counter.count_all().unwrap();
        assert!(!counts.is_empty());
    }

    #[test]
    fn test_parallel_distribution() {
        let config = OptimizedDiskConfig {
            k: 11,
            num_buckets: 4,
            min_count: 1,
            temp_dir: std::env::temp_dir().join("raptor_test_opt_parallel"),
            compression_enabled: true,
            parallel_distribution: true,
            chunk_size: 2,
            ..Default::default()
        };

        let mut counter = OptimizedDiskCounter::new(config).unwrap();

        let seqs: Vec<&str> = (0..10).map(|_| "ACGTACGTACGTACGT").collect();

        counter
            .distribute(seqs.iter().map(|s| s.as_bytes()))
            .unwrap();

        let counts = counter.count_all().unwrap();
        assert!(!counts.is_empty());
    }

    #[test]
    fn parallel_distribution_matches_sequential_with_non_power_of_two_buckets() {
        let temp = TempDir::new().expect("create temp dir");
        let seqs: Vec<&[u8]> = vec![
            b"ACGTACGTACGTACGT",
            b"NNNNACGTACGTNNNN",
            b"TGCATGCATGCATGCA",
            b"ACGTACGTACGTACGT",
            b"TTTTCCCCAAAAGGGG",
            b"GATANNNNCATG",
            b"CCCCCCCCCCCCCCCC",
        ];

        let sequential = run_counter(
            OptimizedDiskConfig {
                k: 7,
                num_buckets: 3,
                min_count: 1,
                temp_dir: temp.path().join("seq_non_pow2"),
                compression_enabled: false,
                parallel_distribution: false,
                chunk_size: 2,
            },
            &seqs,
        );

        let parallel = run_counter(
            OptimizedDiskConfig {
                k: 7,
                num_buckets: 3,
                min_count: 1,
                temp_dir: temp.path().join("par_non_pow2"),
                compression_enabled: false,
                parallel_distribution: true,
                chunk_size: 2,
            },
            &seqs,
        );

        assert_eq!(parallel, sequential);
    }

    #[test]
    fn compressed_and_uncompressed_modes_produce_identical_counts() {
        let temp = TempDir::new().expect("create temp dir");
        let seqs: Vec<&[u8]> = vec![
            b"ACGTACGTACGTACGT",
            b"ACGTACGTACGTACGT",
            b"TGCATGCATGCATGCA",
            b"GATTACAGATTACA",
            b"TTTTTTTTTTTTTTTT",
        ];

        let compressed = run_counter(
            OptimizedDiskConfig {
                k: 9,
                num_buckets: 8,
                min_count: 1,
                temp_dir: temp.path().join("compressed"),
                compression_enabled: true,
                parallel_distribution: true,
                chunk_size: 2,
            },
            &seqs,
        );

        let uncompressed = run_counter(
            OptimizedDiskConfig {
                k: 9,
                num_buckets: 8,
                min_count: 1,
                temp_dir: temp.path().join("uncompressed"),
                compression_enabled: false,
                parallel_distribution: true,
                chunk_size: 2,
            },
            &seqs,
        );

        assert_eq!(compressed, uncompressed);
    }

    #[test]
    fn parallel_distribution_is_input_order_invariant() {
        let temp = TempDir::new().expect("create temp dir");
        let seqs_a: Vec<&[u8]> = vec![
            b"ACGTACGTACGTACGT",
            b"NNNNACGTACGTNNNN",
            b"TGCATGCATGCATGCA",
            b"TTTTCCCCAAAAGGGG",
            b"GATANNNNCATG",
            b"CCCCCCCCCCCCCCCC",
        ];
        let seqs_b: Vec<&[u8]> = seqs_a.iter().copied().rev().collect();

        let counts_a = run_counter(
            OptimizedDiskConfig {
                k: 7,
                num_buckets: 8,
                min_count: 1,
                temp_dir: temp.path().join("order_a"),
                compression_enabled: true,
                parallel_distribution: true,
                chunk_size: 3,
            },
            &seqs_a,
        );

        let counts_b = run_counter(
            OptimizedDiskConfig {
                k: 7,
                num_buckets: 8,
                min_count: 1,
                temp_dir: temp.path().join("order_b"),
                compression_enabled: true,
                parallel_distribution: true,
                chunk_size: 3,
            },
            &seqs_b,
        );

        assert_eq!(counts_a, counts_b);
    }

    #[test]
    fn optimized_counter_rejects_invalid_config() {
        let temp_base = std::env::temp_dir();

        let invalid_k = OptimizedDiskConfig {
            k: 0,
            temp_dir: temp_base.join("raptor_test_opt_invalid_k"),
            ..Default::default()
        };
        assert!(OptimizedDiskCounter::new(invalid_k).is_err());

        let invalid_bucket_count = OptimizedDiskConfig {
            num_buckets: 0,
            temp_dir: temp_base.join("raptor_test_opt_invalid_bucket_count"),
            ..Default::default()
        };
        assert!(OptimizedDiskCounter::new(invalid_bucket_count).is_err());

        let invalid_chunk_size = OptimizedDiskConfig {
            chunk_size: 0,
            temp_dir: temp_base.join("raptor_test_opt_invalid_chunk_size"),
            ..Default::default()
        };
        assert!(OptimizedDiskCounter::new(invalid_chunk_size).is_err());
    }

    #[test]
    fn optimized_counter_rejects_truncated_bucket_file() {
        let config = OptimizedDiskConfig {
            k: 11,
            num_buckets: 1,
            min_count: 1,
            temp_dir: std::env::temp_dir().join("raptor_test_opt_truncated_bucket"),
            compression_enabled: false,
            parallel_distribution: false,
            ..Default::default()
        };
        let mut counter = OptimizedDiskCounter::new(config).expect("counter should construct");
        let path = counter.config.temp_dir.join("bucket_00000.bin");
        fs::create_dir_all(&counter.config.temp_dir).expect("temp dir should exist");
        fs::write(&path, [1u8, 2, 3]).expect("write truncated bucket payload");
        counter.bucket_paths = vec![path];
        counter.bucket_counts = vec![1];

        let err = counter
            .count_all()
            .expect_err("truncated bucket must return an error");
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    }
}
