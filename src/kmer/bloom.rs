// src/kmer/bloom.rs
//! Bloom filter implementation for memory-efficient k-mer pre-filtering.
//!
//! Bloom filters allow probabilistic membership testing with no false negatives.
//! This is ideal for filtering singleton k-mers (sequencing errors) before
//! building the full count table, reducing memory by 30-50%.
//!
//! Two-pass k-mer counting:
//! 1. Pass 1: Insert all k-mers into Bloom filter
//! 2. Pass 2: Only count k-mers that appear in the Bloom filter twice
//!            (singleton k-mers from errors will only appear once)

use std::hash::{Hash, Hasher};
use ahash::AHasher;

/// A space-efficient probabilistic data structure for membership testing.
///
/// False positives are possible but false negatives are not.
/// Optimal for pre-filtering k-mers to remove singletons.
pub struct BloomFilter {
    /// Bit vector storing the filter state
    bits: Vec<u64>,
    /// Number of bits in the filter
    num_bits: usize,
    /// Number of hash functions to use
    num_hashes: usize,
    /// Seed for hash function variation
    seed: u64,
}

impl BloomFilter {
    /// Create a new Bloom filter with specified capacity and false positive rate.
    ///
    /// # Arguments
    /// * `expected_items` - Expected number of items to insert
    /// * `fp_rate` - Desired false positive rate (e.g., 0.01 for 1%)
    ///
    /// # Example
    /// ```ignore
    /// // Create filter for 10 million k-mers with 1% FP rate
    /// let bloom = BloomFilter::with_fp_rate(10_000_000, 0.01);
    /// ```
    pub fn with_fp_rate(expected_items: usize, fp_rate: f64) -> Self {
        // Optimal number of bits: m = -n * ln(p) / (ln(2)^2)
        let ln2_sq = std::f64::consts::LN_2 * std::f64::consts::LN_2;
        let num_bits = (-(expected_items as f64) * fp_rate.ln() / ln2_sq).ceil() as usize;

        // Optimal number of hash functions: k = (m/n) * ln(2)
        let num_hashes = ((num_bits as f64 / expected_items as f64) * std::f64::consts::LN_2)
            .ceil() as usize;

        Self::new(num_bits, num_hashes)
    }

    /// Create a new Bloom filter with explicit size parameters.
    ///
    /// # Arguments
    /// * `num_bits` - Number of bits in the filter
    /// * `num_hashes` - Number of hash functions to use
    pub fn new(num_bits: usize, num_hashes: usize) -> Self {
        // Round up to next multiple of 64 for u64 storage
        let num_bits = ((num_bits + 63) / 64) * 64;
        let num_u64s = num_bits / 64;

        Self {
            bits: vec![0u64; num_u64s],
            num_bits,
            num_hashes: num_hashes.max(1),
            seed: 0x517cc1b727220a95, // Random seed
        }
    }

    /// Create a Bloom filter with a specific memory budget.
    ///
    /// # Arguments
    /// * `bytes` - Maximum memory to use in bytes
    /// * `num_hashes` - Number of hash functions (typically 3-7)
    pub fn with_memory(bytes: usize, num_hashes: usize) -> Self {
        let num_bits = bytes * 8;
        Self::new(num_bits, num_hashes)
    }

    /// Insert an item into the Bloom filter.
    ///
    /// After insertion, `may_contain()` will return true for this item.
    #[inline]
    pub fn insert(&mut self, hash: u64) {
        for i in 0..self.num_hashes {
            let bit_idx = self.get_bit_index(hash, i);
            let word_idx = bit_idx / 64;
            let bit_offset = bit_idx % 64;
            self.bits[word_idx] |= 1u64 << bit_offset;
        }
    }

    /// Insert a hashable item into the Bloom filter.
    #[inline]
    pub fn insert_item<T: Hash>(&mut self, item: &T) {
        let hash = self.hash_item(item);
        self.insert(hash);
    }

    /// Check if an item may be in the filter.
    ///
    /// Returns `true` if the item might be present (could be false positive).
    /// Returns `false` if the item is definitely not present.
    #[inline]
    pub fn may_contain(&self, hash: u64) -> bool {
        for i in 0..self.num_hashes {
            let bit_idx = self.get_bit_index(hash, i);
            let word_idx = bit_idx / 64;
            let bit_offset = bit_idx % 64;
            if (self.bits[word_idx] & (1u64 << bit_offset)) == 0 {
                return false;
            }
        }
        true
    }

    /// Check if a hashable item may be in the filter.
    #[inline]
    pub fn may_contain_item<T: Hash>(&self, item: &T) -> bool {
        let hash = self.hash_item(item);
        self.may_contain(hash)
    }

    /// Get the bit index for hash function i.
    ///
    /// Uses double hashing: h(i) = h1 + i*h2 mod m
    #[inline]
    fn get_bit_index(&self, hash: u64, i: usize) -> usize {
        // Split the 64-bit hash into two 32-bit hashes
        let h1 = hash as u32 as u64;
        let h2 = (hash >> 32) as u32 as u64;

        // Double hashing formula
        let combined = h1.wrapping_add((i as u64).wrapping_mul(h2));
        (combined as usize) % self.num_bits
    }

    /// Hash an item using AHash.
    #[inline]
    fn hash_item<T: Hash>(&self, item: &T) -> u64 {
        let mut hasher = AHasher::default();
        item.hash(&mut hasher);
        hasher.finish() ^ self.seed
    }

    /// Get the number of bits in the filter.
    pub fn num_bits(&self) -> usize {
        self.num_bits
    }

    /// Get the number of hash functions.
    pub fn num_hashes(&self) -> usize {
        self.num_hashes
    }

    /// Get the memory usage in bytes.
    pub fn memory_bytes(&self) -> usize {
        self.bits.len() * 8
    }

    /// Estimate the current false positive rate.
    ///
    /// Based on the number of bits set.
    pub fn estimated_fp_rate(&self) -> f64 {
        let bits_set: usize = self.bits.iter().map(|w| w.count_ones() as usize).sum();
        let fill_ratio = bits_set as f64 / self.num_bits as f64;
        fill_ratio.powi(self.num_hashes as i32)
    }

    /// Clear all bits in the filter.
    pub fn clear(&mut self) {
        self.bits.fill(0);
    }

    /// Merge another Bloom filter into this one (OR operation).
    ///
    /// Both filters must have the same size.
    pub fn merge(&mut self, other: &BloomFilter) {
        assert_eq!(self.num_bits, other.num_bits, "Bloom filters must have same size");
        assert_eq!(self.num_hashes, other.num_hashes, "Bloom filters must have same num_hashes");

        for (a, b) in self.bits.iter_mut().zip(other.bits.iter()) {
            *a |= *b;
        }
    }
}

/// A counting Bloom filter that tracks approximate counts.
///
/// Uses 4 bits per counter, allowing counts up to 15.
/// Useful for identifying k-mers that appear multiple times.
pub struct CountingBloomFilter {
    /// 4-bit counters packed into u64s (16 counters per u64)
    counters: Vec<u64>,
    /// Number of counters
    num_counters: usize,
    /// Number of hash functions
    num_hashes: usize,
    /// Seed for hash variation
    seed: u64,
}

impl CountingBloomFilter {
    /// Create a new counting Bloom filter.
    ///
    /// # Arguments
    /// * `expected_items` - Expected number of items
    /// * `fp_rate` - Desired false positive rate
    pub fn with_fp_rate(expected_items: usize, fp_rate: f64) -> Self {
        let ln2_sq = std::f64::consts::LN_2 * std::f64::consts::LN_2;
        let num_counters = (-(expected_items as f64) * fp_rate.ln() / ln2_sq).ceil() as usize;
        let num_hashes = ((num_counters as f64 / expected_items as f64) * std::f64::consts::LN_2)
            .ceil() as usize;

        Self::new(num_counters, num_hashes)
    }

    /// Create with explicit parameters.
    pub fn new(num_counters: usize, num_hashes: usize) -> Self {
        // Round up to multiple of 16 (counters per u64)
        let num_counters = ((num_counters + 15) / 16) * 16;
        let num_u64s = num_counters / 16;

        Self {
            counters: vec![0u64; num_u64s],
            num_counters,
            num_hashes: num_hashes.max(1),
            seed: 0x517cc1b727220a95,
        }
    }

    /// Increment the count for a hash value.
    ///
    /// Returns the minimum count across all hash positions after increment.
    #[inline]
    pub fn insert(&mut self, hash: u64) -> u8 {
        let mut min_count = 15u8;

        for i in 0..self.num_hashes {
            let counter_idx = self.get_counter_index(hash, i);
            let word_idx = counter_idx / 16;
            let nibble_offset = (counter_idx % 16) * 4;

            let current = ((self.counters[word_idx] >> nibble_offset) & 0xF) as u8;
            if current < 15 {
                // Increment (saturating at 15)
                self.counters[word_idx] += 1u64 << nibble_offset;
                min_count = min_count.min(current + 1);
            } else {
                min_count = min_count.min(15);
            }
        }

        min_count
    }

    /// Get the estimated count for a hash value.
    ///
    /// Returns the minimum count across all hash positions.
    #[inline]
    pub fn get_count(&self, hash: u64) -> u8 {
        let mut min_count = 15u8;

        for i in 0..self.num_hashes {
            let counter_idx = self.get_counter_index(hash, i);
            let word_idx = counter_idx / 16;
            let nibble_offset = (counter_idx % 16) * 4;

            let count = ((self.counters[word_idx] >> nibble_offset) & 0xF) as u8;
            min_count = min_count.min(count);
        }

        min_count
    }

    /// Check if count is at least the given threshold.
    #[inline]
    pub fn count_at_least(&self, hash: u64, threshold: u8) -> bool {
        self.get_count(hash) >= threshold
    }

    /// Get counter index using double hashing.
    #[inline]
    fn get_counter_index(&self, hash: u64, i: usize) -> usize {
        let h1 = hash as u32 as u64;
        let h2 = (hash >> 32) as u32 as u64;
        let combined = h1.wrapping_add((i as u64).wrapping_mul(h2));
        (combined as usize) % self.num_counters
    }

    /// Get memory usage in bytes.
    pub fn memory_bytes(&self) -> usize {
        self.counters.len() * 8
    }

    /// Clear all counters.
    pub fn clear(&mut self) {
        self.counters.fill(0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bloom_filter_basic() {
        let mut bloom = BloomFilter::with_fp_rate(1000, 0.01);

        // Insert some items
        bloom.insert(12345);
        bloom.insert(67890);
        bloom.insert(11111);

        // Should find inserted items
        assert!(bloom.may_contain(12345));
        assert!(bloom.may_contain(67890));
        assert!(bloom.may_contain(11111));

        // Should (probably) not find items we didn't insert
        // Note: false positives are possible
        let mut false_positives = 0;
        for i in 0..1000 {
            if bloom.may_contain(i * 1000000 + 99999) {
                false_positives += 1;
            }
        }

        // With 1% FP rate, we expect roughly 10 false positives out of 1000
        assert!(false_positives < 50, "Too many false positives: {}", false_positives);
    }

    #[test]
    fn test_bloom_filter_hashable() {
        let mut bloom = BloomFilter::with_fp_rate(1000, 0.01);

        bloom.insert_item(&"ACGTACGT");
        bloom.insert_item(&"GCTAGCTA");

        assert!(bloom.may_contain_item(&"ACGTACGT"));
        assert!(bloom.may_contain_item(&"GCTAGCTA"));
    }

    #[test]
    fn test_counting_bloom_filter() {
        let mut cbf = CountingBloomFilter::with_fp_rate(1000, 0.01);

        // Insert same item multiple times
        let hash = 12345u64;
        let c1 = cbf.insert(hash);
        let c2 = cbf.insert(hash);
        let c3 = cbf.insert(hash);

        // Counts should be increasing
        assert!(c1 >= 1);
        assert!(c2 >= 2);
        assert!(c3 >= 3);

        // Check count - should be at least 3
        let count = cbf.get_count(hash);
        assert!(count >= 3, "Expected count >= 3, got {}", count);
        assert!(cbf.count_at_least(hash, 2));
        assert!(cbf.count_at_least(hash, 3));

        // Non-existent item should have low count (may have false positives)
        // But newly created filter should be mostly empty
    }

    #[test]
    fn test_counting_bloom_saturation() {
        let mut cbf = CountingBloomFilter::new(1000, 3);

        let hash = 42u64;

        // Insert 20 times, should saturate at 15
        for _ in 0..20 {
            cbf.insert(hash);
        }

        assert_eq!(cbf.get_count(hash), 15);
    }

    #[test]
    fn test_bloom_merge() {
        let mut bloom1 = BloomFilter::new(1000, 3);
        let mut bloom2 = BloomFilter::new(1000, 3);

        bloom1.insert(111);
        bloom1.insert(222);
        bloom2.insert(333);
        bloom2.insert(444);

        bloom1.merge(&bloom2);

        assert!(bloom1.may_contain(111));
        assert!(bloom1.may_contain(222));
        assert!(bloom1.may_contain(333));
        assert!(bloom1.may_contain(444));
    }

    #[test]
    fn test_memory_usage() {
        // 10 million items, 1% FP rate
        let bloom = BloomFilter::with_fp_rate(10_000_000, 0.01);

        // Should use roughly 12 MB (9.6 bits per item for 1% FP)
        let mb = bloom.memory_bytes() / (1024 * 1024);
        assert!(mb >= 10 && mb <= 20, "Unexpected memory usage: {} MB", mb);
    }
}
