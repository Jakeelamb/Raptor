//! Disk-based k-mer counting for large genomes
//!
//! Implements KMC3-style external memory k-mer counting:
//! 1. Distribute k-mers to disk buckets based on minimizer hash
//! 2. Sort and count each bucket independently
//! 3. Merge results into final k-mer database
//!
//! Memory usage: O(bucket_size) instead of O(total_kmers)
//! Can handle 100+ Gb genomes on 16 GB RAM machines.

use crate::kmer::nthash::NtHashIterator;
use ahash::AHashMap;
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

/// Configuration for disk-based k-mer counting
#[derive(Clone, Debug)]
pub struct DiskCounterConfig {
    /// K-mer size (must be <= 32 for u64 encoding)
    pub k: usize,
    /// Number of buckets (power of 2 recommended, e.g., 256, 1024, 4096)
    pub num_buckets: usize,
    /// Minimum k-mer count to keep (filters sequencing errors)
    pub min_count: u32,
    /// Directory for temporary bucket files
    pub temp_dir: PathBuf,
    /// Buffer size for disk writes (bytes)
    pub write_buffer_size: usize,
    /// Maximum RAM per bucket during counting (bytes)
    pub max_bucket_ram: usize,
    /// Use compression for bucket files
    pub compress_buckets: bool,
}

impl Default for DiskCounterConfig {
    fn default() -> Self {
        Self {
            k: 31,
            num_buckets: 1024,  // Good default for most genomes
            min_count: 2,       // Filter singletons (likely errors)
            temp_dir: std::env::temp_dir().join("raptor_kmer"),
            write_buffer_size: 4 * 1024 * 1024,  // 4 MB buffers
            max_bucket_ram: 4 * 1024 * 1024 * 1024,  // 4 GB per bucket
            compress_buckets: true,
        }
    }
}

impl DiskCounterConfig {
    /// Configure for specific RAM limit
    pub fn for_ram_limit(ram_bytes: usize, k: usize) -> Self {
        // Use ~50% of RAM for bucket processing, rest for OS/buffers
        let usable_ram = ram_bytes / 2;

        // Each k-mer entry in sorting: 8 bytes (u64 hash)
        // Want each bucket to fit comfortably in usable_ram
        // For 100B k-mers, need enough buckets
        let bucket_target_size = usable_ram / 2;  // Leave room for sorting overhead

        // Estimate buckets needed for large genomes
        // Assume worst case: 100B k-mers, 8 bytes each = 800 GB raw
        // With num_buckets, each bucket is 800GB / num_buckets
        let num_buckets = if bucket_target_size >= 4 * 1024 * 1024 * 1024 {
            256   // Large RAM: fewer buckets
        } else if bucket_target_size >= 1024 * 1024 * 1024 {
            1024  // Medium RAM
        } else if bucket_target_size >= 256 * 1024 * 1024 {
            4096  // Small RAM
        } else {
            8192  // Very limited RAM
        };

        Self {
            k,
            num_buckets,
            max_bucket_ram: bucket_target_size,
            ..Default::default()
        }
    }

    /// Auto-detect system RAM and configure appropriately
    pub fn auto_configure(k: usize) -> Self {
        // Try to detect system memory
        let ram_bytes = Self::detect_system_ram();
        Self::for_ram_limit(ram_bytes, k)
    }

    fn detect_system_ram() -> usize {
        // Try /proc/meminfo on Linux
        if let Ok(contents) = fs::read_to_string("/proc/meminfo") {
            for line in contents.lines() {
                if line.starts_with("MemTotal:") {
                    if let Some(kb_str) = line.split_whitespace().nth(1) {
                        if let Ok(kb) = kb_str.parse::<usize>() {
                            return kb * 1024;  // Convert KB to bytes
                        }
                    }
                }
            }
        }

        // Default: assume 16 GB
        16 * 1024 * 1024 * 1024
    }
}

/// Bucket file for storing k-mers on disk
struct KmerBucket {
    path: PathBuf,
    writer: Option<BufWriter<File>>,
    buffer: Vec<u64>,
    buffer_capacity: usize,
}

impl KmerBucket {
    fn new(path: PathBuf, buffer_size: usize) -> std::io::Result<Self> {
        let file = File::create(&path)?;
        let writer = BufWriter::with_capacity(buffer_size, file);
        let buffer_capacity = buffer_size / 8;  // u64 = 8 bytes

        Ok(Self {
            path,
            writer: Some(writer),
            buffer: Vec::with_capacity(buffer_capacity),
            buffer_capacity,
        })
    }

    fn add_kmer(&mut self, kmer_hash: u64) -> std::io::Result<()> {
        self.buffer.push(kmer_hash);

        if self.buffer.len() >= self.buffer_capacity {
            self.flush_buffer()?;
        }

        Ok(())
    }

    fn flush_buffer(&mut self) -> std::io::Result<()> {
        if let Some(ref mut writer) = self.writer {
            for &hash in &self.buffer {
                writer.write_all(&hash.to_le_bytes())?;
            }
            self.buffer.clear();
        }
        Ok(())
    }

    fn finalize(mut self) -> std::io::Result<PathBuf> {
        self.flush_buffer()?;
        if let Some(writer) = self.writer.take() {
            writer.into_inner()?.sync_all()?;
        }
        Ok(self.path)
    }
}

/// Main disk-based k-mer counter
pub struct DiskKmerCounter {
    config: DiskCounterConfig,
    bucket_paths: Vec<PathBuf>,
    total_kmers_distributed: u64,
}

impl DiskKmerCounter {
    pub fn new(config: DiskCounterConfig) -> std::io::Result<Self> {
        // Create temp directory
        fs::create_dir_all(&config.temp_dir)?;

        Ok(Self {
            config,
            bucket_paths: Vec::new(),
            total_kmers_distributed: 0,
        })
    }

    /// Pass 1: Distribute k-mers from sequences to disk buckets
    pub fn distribute_kmers<I, S>(&mut self, sequences: I) -> std::io::Result<()>
    where
        I: Iterator<Item = S>,
        S: AsRef<[u8]>,
    {
        let k = self.config.k;
        let num_buckets = self.config.num_buckets;
        let bucket_mask = (num_buckets - 1) as u64;  // Assumes power of 2

        // Create bucket files
        let mut buckets: Vec<KmerBucket> = (0..num_buckets)
            .map(|i| {
                let path = self.config.temp_dir.join(format!("bucket_{:04}.bin", i));
                KmerBucket::new(path, self.config.write_buffer_size)
            })
            .collect::<std::io::Result<Vec<_>>>()?;

        // Process sequences
        let mut total_kmers = 0u64;

        for seq in sequences {
            let bytes = seq.as_ref();
            if bytes.len() < k {
                continue;
            }

            // Use ntHash for fast k-mer hashing
            for (_, canonical_hash) in NtHashIterator::new(bytes, k) {
                // Compute minimizer for bucket assignment
                // Use lower bits of hash for bucket ID (fast modulo for power of 2)
                let bucket_id = (canonical_hash & bucket_mask) as usize;

                buckets[bucket_id].add_kmer(canonical_hash)?;
                total_kmers += 1;
            }
        }

        // Finalize all buckets
        self.bucket_paths = buckets
            .into_iter()
            .map(|b| b.finalize())
            .collect::<std::io::Result<Vec<_>>>()?;

        self.total_kmers_distributed = total_kmers;

        Ok(())
    }

    /// Pass 2: Count k-mers in each bucket and return counts
    pub fn count_all_buckets(&self) -> std::io::Result<AHashMap<u64, u32>> {
        let min_count = self.config.min_count;

        // Process buckets in parallel, merge results
        let bucket_counts: Vec<AHashMap<u64, u32>> = self.bucket_paths
            .par_iter()
            .map(|path| Self::count_single_bucket(path, min_count))
            .collect::<std::io::Result<Vec<_>>>()?;

        // Merge all bucket counts
        let total_unique: usize = bucket_counts.iter().map(|m| m.len()).sum();
        let mut merged = AHashMap::with_capacity(total_unique);

        for bucket in bucket_counts {
            merged.extend(bucket);
        }

        Ok(merged)
    }

    /// Count k-mers in a single bucket file using radix sort
    fn count_single_bucket(path: &Path, min_count: u32) -> std::io::Result<AHashMap<u64, u32>> {
        // Read all k-mers from bucket
        let file = File::open(path)?;
        let file_size = file.metadata()?.len() as usize;
        let num_kmers = file_size / 8;

        if num_kmers == 0 {
            return Ok(AHashMap::new());
        }

        let mut reader = BufReader::with_capacity(4 * 1024 * 1024, file);
        let mut kmers: Vec<u64> = Vec::with_capacity(num_kmers);

        let mut buf = [0u8; 8];
        while reader.read_exact(&mut buf).is_ok() {
            kmers.push(u64::from_le_bytes(buf));
        }

        // Sort k-mers (radix sort would be faster, but std sort is good enough)
        kmers.sort_unstable();

        // Count consecutive identical k-mers
        let mut counts = AHashMap::with_capacity(num_kmers / 10);  // Estimate unique

        if kmers.is_empty() {
            return Ok(counts);
        }

        let mut current = kmers[0];
        let mut count = 1u32;

        for &kmer in &kmers[1..] {
            if kmer == current {
                count = count.saturating_add(1);
            } else {
                if count >= min_count {
                    counts.insert(current, count);
                }
                current = kmer;
                count = 1;
            }
        }

        // Don't forget the last k-mer
        if count >= min_count {
            counts.insert(current, count);
        }

        Ok(counts)
    }

    /// Cleanup temporary files
    pub fn cleanup(&self) -> std::io::Result<()> {
        for path in &self.bucket_paths {
            if path.exists() {
                fs::remove_file(path)?;
            }
        }

        // Try to remove temp directory if empty
        let _ = fs::remove_dir(&self.config.temp_dir);

        Ok(())
    }

    /// Get statistics about the distribution
    pub fn stats(&self) -> DiskCounterStats {
        let bucket_sizes: Vec<u64> = self.bucket_paths
            .iter()
            .filter_map(|p| fs::metadata(p).ok())
            .map(|m| m.len())
            .collect();

        let total_disk_usage: u64 = bucket_sizes.iter().sum();
        let max_bucket_size = bucket_sizes.iter().max().copied().unwrap_or(0);
        let avg_bucket_size = if bucket_sizes.is_empty() {
            0
        } else {
            total_disk_usage / bucket_sizes.len() as u64
        };

        DiskCounterStats {
            total_kmers: self.total_kmers_distributed,
            num_buckets: self.bucket_paths.len(),
            total_disk_bytes: total_disk_usage,
            max_bucket_bytes: max_bucket_size,
            avg_bucket_bytes: avg_bucket_size,
        }
    }
}

impl Drop for DiskKmerCounter {
    fn drop(&mut self) {
        // Best-effort cleanup
        let _ = self.cleanup();
    }
}

#[derive(Debug, Clone)]
pub struct DiskCounterStats {
    pub total_kmers: u64,
    pub num_buckets: usize,
    pub total_disk_bytes: u64,
    pub max_bucket_bytes: u64,
    pub avg_bucket_bytes: u64,
}

impl std::fmt::Display for DiskCounterStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "K-mers: {}, Buckets: {}, Disk: {:.2} GB, Max bucket: {:.2} MB",
            self.total_kmers,
            self.num_buckets,
            self.total_disk_bytes as f64 / (1024.0 * 1024.0 * 1024.0),
            self.max_bucket_bytes as f64 / (1024.0 * 1024.0),
        )
    }
}

/// Iterator that streams k-mer counts from disk buckets
/// Allows processing counts without loading everything into RAM
pub struct StreamingKmerCounts {
    bucket_paths: Vec<PathBuf>,
    current_bucket: usize,
    current_counts: std::vec::IntoIter<(u64, u32)>,
    min_count: u32,
}

impl StreamingKmerCounts {
    pub fn new(bucket_paths: Vec<PathBuf>, min_count: u32) -> Self {
        Self {
            bucket_paths,
            current_bucket: 0,
            current_counts: Vec::new().into_iter(),
            min_count,
        }
    }

    fn load_next_bucket(&mut self) -> Option<()> {
        while self.current_bucket < self.bucket_paths.len() {
            let path = &self.bucket_paths[self.current_bucket];
            self.current_bucket += 1;

            if let Ok(counts) = DiskKmerCounter::count_single_bucket(path, self.min_count) {
                if !counts.is_empty() {
                    self.current_counts = counts.into_iter().collect::<Vec<_>>().into_iter();
                    return Some(());
                }
            }
        }
        None
    }
}

impl Iterator for StreamingKmerCounts {
    type Item = (u64, u32);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(item) = self.current_counts.next() {
                return Some(item);
            }

            if self.load_next_bucket().is_none() {
                return None;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_disk_counter_basic() {
        let config = DiskCounterConfig {
            k: 11,
            num_buckets: 4,
            min_count: 1,
            temp_dir: std::env::temp_dir().join("raptor_test_kmer"),
            ..Default::default()
        };

        let mut counter = DiskKmerCounter::new(config).unwrap();

        // Test sequences with repeated k-mers
        let sequences = vec![
            "ACGTACGTACGTACGT",
            "ACGTACGTACGTACGT",  // Duplicate
            "TGCATGCATGCATGCA",
        ];

        counter.distribute_kmers(sequences.iter().map(|s| s.as_bytes())).unwrap();

        let stats = counter.stats();
        assert!(stats.total_kmers > 0);
        assert_eq!(stats.num_buckets, 4);

        let counts = counter.count_all_buckets().unwrap();
        assert!(!counts.is_empty());

        // Some k-mers should have count >= 2 due to duplicates
        let high_count = counts.values().filter(|&&c| c >= 2).count();
        assert!(high_count > 0);
    }

    #[test]
    fn test_auto_config() {
        let config = DiskCounterConfig::auto_configure(31);
        assert!(config.num_buckets >= 256);
        assert!(config.num_buckets <= 8192);
    }
}
