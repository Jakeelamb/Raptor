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
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

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
    total_kmers: u64,
}

impl DiskKmerCounterV2 {
    pub fn new(config: DiskCounterConfig) -> std::io::Result<Self> {
        fs::create_dir_all(&config.temp_dir)?;
        Ok(Self {
            config,
            bucket_paths: Vec::new(),
            bucket_counts: Vec::new(),
            total_kmers: 0,
        })
    }

    /// Pass 1: Distribute encoded k-mers to disk buckets
    pub fn distribute<I, S>(&mut self, sequences: I) -> std::io::Result<()>
    where
        I: Iterator<Item = S>,
        S: AsRef<[u8]>,
    {
        let k = self.config.k;
        let num_buckets = self.config.num_buckets;
        let mask = (num_buckets - 1) as u64;

        // Create buckets
        let mut buckets: Vec<KmerBucket> = (0..num_buckets)
            .map(|i| {
                let path = self.config.temp_dir.join(format!("bucket_{:05}.bin", i));
                KmerBucket::new(path, self.config.write_buffer_size)
            })
            .collect::<std::io::Result<Vec<_>>>()?;

        let mut total = 0u64;

        for seq in sequences {
            let bytes = seq.as_ref();
            if bytes.len() < k {
                continue;
            }

            // Extract k-mers using KmerU64 for proper encoding
            for i in 0..=bytes.len() - k {
                if let Some(kmer) = KmerU64::from_slice(&bytes[i..i + k]) {
                    // Use canonical form for consistent counting
                    let canonical = kmer.canonical();
                    let encoded = canonical.encoded;

                    // Bucket by lower bits of encoded value
                    let bucket_id = (encoded & mask) as usize;
                    buckets[bucket_id].add(encoded)?;
                    total += 1;
                }
            }
        }

        // Finalize buckets
        let results: Vec<(PathBuf, u64)> = buckets
            .into_iter()
            .map(|b| b.finalize())
            .collect::<std::io::Result<Vec<_>>>()?;

        self.bucket_paths = results.iter().map(|(p, _)| p.clone()).collect();
        self.bucket_counts = results.iter().map(|(_, c)| *c).collect();
        self.total_kmers = total;

        Ok(())
    }

    /// Pass 2: Count k-mers and return map of encoded -> count
    pub fn count_all(&self) -> std::io::Result<AHashMap<u64, u32>> {
        let min_count = self.config.min_count;

        // Process buckets in parallel
        let bucket_results: Vec<AHashMap<u64, u32>> = self.bucket_paths
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
        let num_kmers = file_size / 8;

        if num_kmers == 0 {
            return Ok(AHashMap::new());
        }

        // Read all k-mers
        let mut reader = BufReader::with_capacity(4 * 1024 * 1024, file);
        let mut kmers = Vec::with_capacity(num_kmers);
        let mut buf = [0u8; 8];

        while reader.read_exact(&mut buf).is_ok() {
            kmers.push(u64::from_le_bytes(buf));
        }

        // Sort for counting
        kmers.sort_unstable();

        // Count consecutive
        let mut counts = AHashMap::with_capacity(num_kmers / 5);
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

        if count >= min_count {
            counts.insert(current, count);
        }

        Ok(counts)
    }

    /// Get statistics
    pub fn stats(&self) -> (u64, usize, u64) {
        let disk_bytes: u64 = self.bucket_paths
            .iter()
            .filter_map(|p| fs::metadata(p).ok())
            .map(|m| m.len())
            .sum();
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

impl Drop for DiskKmerCounterV2 {
    fn drop(&mut self) {
        let _ = self.cleanup();
    }
}

/// Decode a u64-encoded k-mer back to a string
pub fn decode_kmer(encoded: u64, k: usize) -> String {
    KmerU64 { encoded, len: k as u8 }.decode()
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

    let mask = (1u64 << (k * 2)) - 1;
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
    let kmer = KmerU64 { encoded, len: k as u8 };
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
        let result = KmerU64 { encoded: extended, len: 4 };
        assert_eq!(result.decode(), "CGTA");
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

        counter.distribute(seqs.iter().map(|s| s.as_bytes())).unwrap();

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
}
