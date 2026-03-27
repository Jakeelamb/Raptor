//! Disk-based k-mer counting V2 - with sequence reconstruction
//!
//! This version stores encoded k-mers (not hashes) so we can reconstruct
//! the actual sequences for assembly. Uses KmerU64 encoding which is
//! reversible and canonical.
//!
//! Key insight: For k≤32, a 2-bit encoded k-mer fits in u64 (same size as hash!)
//! but is fully reversible back to the original sequence.

use crate::kmer::kmer::KmerU64;
use ahash::AHashMap;
use memmap2::Mmap;
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

#[inline]
fn kmer_mask(k: usize) -> u64 {
    if k >= 32 {
        u64::MAX
    } else {
        (1u64 << (k * 2)) - 1
    }
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

#[inline]
fn count_sorted_kmers(sorted_kmers: &[u64], min_count: u32) -> AHashMap<u64, u32> {
    let mut counts = AHashMap::with_capacity(sorted_kmers.len() / 4);
    let mut idx = 0usize;

    while idx < sorted_kmers.len() {
        let current = sorted_kmers[idx];
        idx += 1;
        let mut count = 1u32;

        while idx < sorted_kmers.len() && sorted_kmers[idx] == current {
            count = count.saturating_add(1);
            idx += 1;
        }

        if count >= min_count {
            counts.insert(current, count);
        }
    }

    counts
}

/// Configuration for disk-based k-mer counting
#[derive(Clone, Debug)]
pub struct DiskCounterConfig {
    pub k: usize,
    pub num_buckets: usize,
    pub min_count: u32,
    pub temp_dir: PathBuf,
    pub write_buffer_size: usize,
}

impl Default for DiskCounterConfig {
    fn default() -> Self {
        Self {
            k: 31,
            num_buckets: 1024,
            min_count: 2,
            temp_dir: std::env::temp_dir().join("raptor_kmer_v2"),
            write_buffer_size: 4 * 1024 * 1024,
        }
    }
}

impl DiskCounterConfig {
    /// Auto-configure based on available RAM
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
        16 * 1024 * 1024 * 1024 // Default 16 GB
    }
}

/// A bucket file storing encoded k-mers
struct KmerBucket {
    path: PathBuf,
    writer: BufWriter<File>,
    count: u64,
}

impl KmerBucket {
    fn new(path: PathBuf, buffer_size: usize) -> std::io::Result<Self> {
        let file = File::create(&path)?;
        Ok(Self {
            path,
            writer: BufWriter::with_capacity(buffer_size, file),
            count: 0,
        })
    }

    fn open_append(path: PathBuf, buffer_size: usize) -> std::io::Result<Self> {
        let existing_count = fs::metadata(&path).map(|m| m.len() / 8).unwrap_or(0);
        let file = fs::OpenOptions::new().append(true).open(&path)?;
        Ok(Self {
            path,
            writer: BufWriter::with_capacity(buffer_size, file),
            count: existing_count,
        })
    }

    fn add(&mut self, encoded: u64) -> std::io::Result<()> {
        self.writer.write_all(&encoded.to_le_bytes())?;
        self.count += 1;
        Ok(())
    }

    fn finalize(self) -> std::io::Result<(PathBuf, u64)> {
        let mut writer = self.writer;
        writer.flush()?;
        drop(writer);
        Ok((self.path, self.count))
    }
}

/// Disk-based k-mer counter with sequence reconstruction support
pub struct DiskKmerCounterV2 {
    config: DiskCounterConfig,
    bucket_paths: Vec<PathBuf>,
    bucket_counts: Vec<u64>,
    open_buckets: Option<Vec<KmerBucket>>,
    total_kmers: u64,
}

impl DiskKmerCounterV2 {
    pub fn new(config: DiskCounterConfig) -> std::io::Result<Self> {
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
        if config.write_buffer_size == 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "write_buffer_size must be greater than zero",
            ));
        }

        fs::create_dir_all(&config.temp_dir)?;
        let num_buckets = config.num_buckets;
        let bucket_paths = (0..num_buckets)
            .map(|i| config.temp_dir.join(format!("bucket_{:05}.bin", i)))
            .collect();
        Ok(Self {
            config,
            bucket_paths,
            bucket_counts: vec![0; num_buckets],
            open_buckets: None,
            total_kmers: 0,
        })
    }

    fn ensure_buckets_open(&mut self) -> std::io::Result<()> {
        if self.open_buckets.is_some() {
            return Ok(());
        }

        let append_existing = self.bucket_counts.iter().any(|&count| count > 0);
        let buckets = self
            .bucket_paths
            .iter()
            .cloned()
            .map(|path| {
                if append_existing {
                    KmerBucket::open_append(path, self.config.write_buffer_size)
                } else {
                    KmerBucket::new(path, self.config.write_buffer_size)
                }
            })
            .collect::<std::io::Result<Vec<_>>>()?;
        self.open_buckets = Some(buckets);
        Ok(())
    }

    fn finalize_open_buckets(&mut self) -> std::io::Result<()> {
        let Some(buckets) = self.open_buckets.take() else {
            return Ok(());
        };

        for (bucket_idx, bucket) in buckets.into_iter().enumerate() {
            let (_, count) = bucket.finalize()?;
            self.bucket_counts[bucket_idx] = count;
        }

        Ok(())
    }

    /// Pass 1: Distribute encoded k-mers to disk buckets.
    /// Can be called multiple times — appends to existing bucket files.
    pub fn distribute<I, S>(&mut self, sequences: I) -> std::io::Result<()>
    where
        I: Iterator<Item = S>,
        S: AsRef<[u8]>,
    {
        let k = self.config.k;
        let num_buckets = self.config.num_buckets;
        let use_mask_bucket = num_buckets.is_power_of_two();
        let bucket_mask = (num_buckets.saturating_sub(1)) as u64;
        let rolling_mask = kmer_mask(k);

        self.ensure_buckets_open()?;
        let buckets = self
            .open_buckets
            .as_mut()
            .expect("buckets must be initialized before distribution");

        let mut batch_total = 0u64;

        for seq in sequences {
            let bytes = seq.as_ref();
            if bytes.len() < k {
                continue;
            }

            // Rolling extraction avoids O(n*k) repeated re-encoding per window.
            let mut rolling = 0u64;
            let mut valid_run = 0usize;

            for &base in bytes {
                if let Some(base_bits) = encode_base_2bit(base) {
                    rolling = ((rolling << 2) | base_bits) & rolling_mask;
                    valid_run += 1;

                    if valid_run >= k {
                        let canonical_encoded = canonical(rolling, k);
                        let bucket_id = if use_mask_bucket {
                            (canonical_encoded & bucket_mask) as usize
                        } else {
                            (canonical_encoded % num_buckets as u64) as usize
                        };
                        buckets[bucket_id].add(canonical_encoded)?;
                        batch_total += 1;
                    }
                } else {
                    rolling = 0;
                    valid_run = 0;
                }
            }
        }

        for (bucket_idx, bucket) in buckets.iter().enumerate() {
            self.bucket_counts[bucket_idx] = bucket.count;
        }
        self.total_kmers += batch_total;

        Ok(())
    }

    /// Pass 2: Count k-mers and return map of encoded -> count
    pub fn count_all(&mut self) -> std::io::Result<AHashMap<u64, u32>> {
        self.finalize_open_buckets()?;
        let min_count = self.config.min_count;

        // Process buckets in parallel
        let bucket_results: Vec<AHashMap<u64, u32>> = self
            .bucket_paths
            .par_iter()
            .enumerate()
            .filter(|(i, _)| self.bucket_counts[*i] > 0)
            .map(|(_, path)| Self::count_bucket(path, min_count))
            .collect::<std::io::Result<Vec<_>>>()?;

        // Merge results
        let total_unique: usize = bucket_results.iter().map(|m| m.len()).sum();
        let mut merged = AHashMap::with_capacity(total_unique);
        for bucket in bucket_results {
            merged.extend(bucket);
        }

        Ok(merged)
    }

    fn count_bucket(path: &Path, min_count: u32) -> std::io::Result<AHashMap<u64, u32>> {
        let file = File::open(path)?;
        let file_size = file.metadata()?.len() as usize;
        if file_size % 8 != 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "bucket file '{}' has invalid byte length {} (not divisible by 8)",
                    path.display(),
                    file_size
                ),
            ));
        }
        let num_kmers = file_size / 8;

        if num_kmers == 0 {
            return Ok(AHashMap::new());
        }

        let mmap = unsafe { Mmap::map(&file)? };
        let mut kmers = Vec::with_capacity(num_kmers);

        for chunk in mmap.chunks_exact(8) {
            let bytes: [u8; 8] = chunk.try_into().expect("chunk size must be 8");
            kmers.push(u64::from_le_bytes(bytes));
        }

        kmers.sort_unstable();
        Ok(count_sorted_kmers(&kmers, min_count))
    }

    /// Get statistics
    pub fn stats(&self) -> (u64, usize, u64) {
        let disk_bytes = self.bucket_counts.iter().fold(0u64, |acc, count| {
            acc.saturating_add(count.saturating_mul(8))
        });
        (self.total_kmers, self.bucket_paths.len(), disk_bytes)
    }

    /// Cleanup temporary files
    pub fn cleanup(&mut self) -> std::io::Result<()> {
        self.open_buckets.take();
        for path in &self.bucket_paths {
            let _ = fs::remove_file(path);
        }
        let _ = fs::remove_dir(&self.config.temp_dir);
        Ok(())
    }
}

impl Drop for DiskKmerCounterV2 {
    fn drop(&mut self) {
        let _ = self.cleanup();
    }
}

/// Decode a u64-encoded k-mer back to a string
pub fn decode_kmer(encoded: u64, k: usize) -> String {
    KmerU64 {
        encoded,
        len: k as u8,
    }
    .decode()
}

/// Get the suffix of an encoded k-mer (last k-1 bases)
pub fn get_suffix(encoded: u64, k: usize) -> u64 {
    let mask = (1u64 << ((k - 1) * 2)) - 1;
    encoded & mask
}

/// Get the prefix of an encoded k-mer (first k-1 bases)
pub fn get_prefix(encoded: u64, _k: usize) -> u64 {
    encoded >> 2
}

/// Extend encoded k-mer with a new base at the end
pub fn extend_right(encoded: u64, base: u8, k: usize) -> Option<u64> {
    let base_bits = match base {
        b'A' | b'a' => 0u64,
        b'C' | b'c' => 1u64,
        b'G' | b'g' => 2u64,
        b'T' | b't' => 3u64,
        _ => return None,
    };

    let mask = kmer_mask(k);
    Some(((encoded << 2) | base_bits) & mask)
}

/// Extend encoded k-mer with a new base at the beginning
pub fn extend_left(encoded: u64, base: u8, k: usize) -> Option<u64> {
    let base_bits = match base {
        b'A' | b'a' => 0u64,
        b'C' | b'c' => 1u64,
        b'G' | b'g' => 2u64,
        b'T' | b't' => 3u64,
        _ => return None,
    };

    Some((encoded >> 2) | (base_bits << ((k - 1) * 2)))
}

/// Get canonical form of encoded k-mer
pub fn canonical(encoded: u64, k: usize) -> u64 {
    let kmer = KmerU64 {
        encoded,
        len: k as u8,
    };
    kmer.canonical().encoded
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_roundtrip() {
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACG"; // 31-mer
        let kmer = KmerU64::from_str(seq).unwrap();
        let decoded = kmer.decode();
        assert_eq!(seq, decoded);
    }

    #[test]
    fn test_extend_right() {
        let kmer = KmerU64::from_str("ACGT").unwrap();
        let extended = extend_right(kmer.encoded, b'A', 4).unwrap();
        let result = KmerU64 {
            encoded: extended,
            len: 4,
        };
        assert_eq!(result.decode(), "CGTA");
    }

    #[test]
    fn test_extend_right_handles_k32_without_shift_overflow() {
        let input = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let kmer = KmerU64::from_str(input).unwrap();
        let extended = extend_right(kmer.encoded, b'T', 32).unwrap();
        let result = KmerU64 {
            encoded: extended,
            len: 32,
        };
        assert_eq!(result.decode(), "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT");
    }

    #[test]
    fn count_sorted_kmers_applies_min_count_filter() {
        let counts = count_sorted_kmers(&[1, 1, 1, 7, 7, 9], 2);
        assert_eq!(counts.get(&1), Some(&3));
        assert_eq!(counts.get(&7), Some(&2));
        assert!(!counts.contains_key(&9));
    }

    #[test]
    fn test_disk_counter() {
        let config = DiskCounterConfig {
            k: 11,
            num_buckets: 4,
            min_count: 1,
            temp_dir: std::env::temp_dir().join("raptor_test_v2"),
            ..Default::default()
        };

        let mut counter = DiskKmerCounterV2::new(config).unwrap();

        let seqs = vec![
            "ACGTACGTACGTACGT",
            "ACGTACGTACGTACGT", // Duplicate for counting
            "TGCATGCATGCATGCA",
        ];

        counter
            .distribute(seqs.iter().map(|s| s.as_bytes()))
            .unwrap();

        let (total, buckets, _) = counter.stats();
        assert!(total > 0);
        assert_eq!(buckets, 4);

        let counts = counter.count_all().unwrap();
        assert!(!counts.is_empty());

        // Verify we can decode k-mers
        for (&encoded, &count) in counts.iter().take(5) {
            let seq = decode_kmer(encoded, 11);
            assert_eq!(seq.len(), 11);
            assert!(count >= 1);
        }
    }

    #[test]
    fn disk_counter_matches_reference_with_ambiguous_bases_and_non_power_of_two_buckets() {
        let config = DiskCounterConfig {
            k: 5,
            num_buckets: 3,
            min_count: 1,
            temp_dir: std::env::temp_dir().join("raptor_test_v2_reference_match"),
            ..Default::default()
        };
        let mut counter = DiskKmerCounterV2::new(config).expect("counter should construct");

        let sequences = vec!["ACGTNACGTACGT", "ttttacgt", "NNNNN", "ACGTACGT", "acgtaa"];
        counter
            .distribute(sequences.iter().map(|s| s.as_bytes()))
            .expect("distribution should succeed");

        let observed = counter.count_all().expect("counting should succeed");

        let mut expected = AHashMap::new();
        for sequence in &sequences {
            for window in sequence.as_bytes().windows(5) {
                if let Some(kmer) = KmerU64::from_slice(window) {
                    let canonical = kmer.canonical().encoded;
                    *expected.entry(canonical).or_insert(0u32) += 1;
                }
            }
        }

        assert_eq!(observed, expected);
    }

    #[test]
    fn disk_counter_rejects_invalid_config() {
        let temp_base = std::env::temp_dir();

        let invalid_k = DiskCounterConfig {
            k: 0,
            temp_dir: temp_base.join("raptor_test_invalid_k"),
            ..Default::default()
        };
        assert!(DiskKmerCounterV2::new(invalid_k).is_err());

        let invalid_bucket_count = DiskCounterConfig {
            num_buckets: 0,
            temp_dir: temp_base.join("raptor_test_invalid_bucket_count"),
            ..Default::default()
        };
        assert!(DiskKmerCounterV2::new(invalid_bucket_count).is_err());

        let invalid_buffer = DiskCounterConfig {
            write_buffer_size: 0,
            temp_dir: temp_base.join("raptor_test_invalid_buffer"),
            ..Default::default()
        };
        assert!(DiskKmerCounterV2::new(invalid_buffer).is_err());
    }

    #[test]
    fn count_all_rejects_truncated_bucket_files() {
        let config = DiskCounterConfig {
            k: 11,
            num_buckets: 1,
            min_count: 1,
            temp_dir: std::env::temp_dir().join("raptor_test_v2_truncated_bucket"),
            ..Default::default()
        };
        let mut counter = DiskKmerCounterV2::new(config).expect("counter should construct");
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

    #[test]
    fn count_bucket_matches_reference_for_unsorted_duplicates() {
        let temp_dir = tempfile::tempdir().expect("temp dir");
        let path = temp_dir.path().join("bucket_00000.bin");
        let values = [
            9u64, 1, 9, 2, 1, 1, 9, 2, 2, 2, 5, 5, 5, 5, 7, 7, 8, 8, 8, 8,
        ];

        let mut bytes = Vec::with_capacity(values.len() * 8);
        for value in values {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        fs::write(&path, bytes).expect("bucket file should be written");

        let observed = DiskKmerCounterV2::count_bucket(&path, 2).expect("count should succeed");

        let mut expected = AHashMap::new();
        expected.insert(1, 3);
        expected.insert(2, 4);
        expected.insert(5, 4);
        expected.insert(7, 2);
        expected.insert(8, 4);
        expected.insert(9, 3);

        assert_eq!(observed, expected);
    }
}
