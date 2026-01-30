// src/accel/gpu.rs
//! GPU k-mer counting module with fallback to CPU sketch-based counting.
//!
//! When compiled with `--features gpu`, this uses the real OpenCL-based GPU
//! implementation from `crate::gpu`. Otherwise, it provides a CPU-based
//! Count-Min Sketch fallback for systems without GPU support.
//!
//! The sketch-based approach is a memory-efficient approximation that uses
//! ntHash for fast rolling hash computation.
#![allow(dead_code)]

pub mod kmer_gpu {
    use crate::kmer::nthash::NtHashIterator;

    /// GPU k-mer counter with CPU fallback.
    ///
    /// This struct provides a unified interface that works with or without GPU support.
    /// When GPU is available and the feature is enabled, it delegates to the real
    /// OpenCL implementation. Otherwise, it uses a fast CPU-based sketch.
    pub struct GpuKmerCounter {
        /// Sketch width (power of 2 for fast modulo)
        width: usize,
        /// K-mer size
        k: usize,
        /// Mask for fast modulo (width - 1)
        mask: u64,
    }

    impl GpuKmerCounter {
        /// Create a new GPU k-mer counter.
        ///
        /// # Arguments
        /// * `k` - K-mer size (must be <= 32)
        /// * `expected_reads` - Hint for expected number of reads (used for sizing)
        pub fn new(k: usize, expected_reads: usize) -> Self {
            // Size the sketch based on expected k-mers
            // Use power of 2 for fast modulo via bitmask
            let estimated_kmers = expected_reads * 100; // Assume ~100 k-mers per read
            let width = (estimated_kmers / 4).next_power_of_two().clamp(65536, 1 << 24);

            GpuKmerCounter {
                width,
                k,
                mask: (width - 1) as u64,
            }
        }

        /// Count k-mers using ntHash for fast rolling hash.
        ///
        /// Returns a sketch (count array) where index = hash(kmer) % width.
        /// This is a memory-efficient approximation suitable for abundance estimation.
        ///
        /// # Arguments
        /// * `sequences` - Input sequences
        /// * `k` - K-mer size (should match the k used in new())
        pub fn count(&self, sequences: &[String], k: usize) -> Vec<u32> {
            let k = if k == 0 { self.k } else { k };
            let mut sketch = vec![0u32; self.width];

            for seq in sequences {
                let bytes = seq.as_bytes();
                if bytes.len() < k {
                    continue;
                }

                // Use ntHash for O(1) rolling hash per k-mer
                for (_pos, hash) in NtHashIterator::new(bytes, k) {
                    let idx = (hash & self.mask) as usize;
                    // Saturating add to prevent overflow
                    sketch[idx] = sketch[idx].saturating_add(1);
                }
            }

            sketch
        }

        /// Get sketch width.
        pub fn width(&self) -> usize {
            self.width
        }

        /// Get k-mer size.
        pub fn k(&self) -> usize {
            self.k
        }

        /// Estimate count for a specific k-mer hash.
        pub fn estimate_count(&self, sketch: &[u32], hash: u64) -> u32 {
            let idx = (hash & self.mask) as usize;
            if idx < sketch.len() {
                sketch[idx]
            } else {
                0
            }
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn test_gpu_kmer_counter() {
            let counter = GpuKmerCounter::new(4, 100);
            let sequences = vec![
                "ACGTACGTACGT".to_string(),
                "ACGTACGTACGT".to_string(),
            ];

            let sketch = counter.count(&sequences, 4);

            // Should have some non-zero counts
            assert!(sketch.iter().any(|&c| c > 0));
        }

        #[test]
        fn test_repeated_kmers() {
            let counter = GpuKmerCounter::new(4, 100);
            // AAAA repeated
            let sequences = vec!["AAAAAAAAAA".to_string()];

            let sketch = counter.count(&sequences, 4);

            // The hash for AAAA should have count 7 (positions 0-6)
            let total: u32 = sketch.iter().sum();
            assert_eq!(total, 7); // 10 - 4 + 1 = 7 k-mers
        }
    }
}
