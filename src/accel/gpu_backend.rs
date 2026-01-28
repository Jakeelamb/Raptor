use crate::accel::backend::{AdjacencyTable, ComputeBackend};
use std::collections::HashMap;

#[cfg(feature = "gpu")]
use crate::gpu::kmer_gpu::GpuKmerCounter;
#[cfg(feature = "gpu")]
use crate::gpu::overlap_gpu::GpuOverlapFinder;

/// GPU-accelerated compute backend using OpenCL
///
/// Automatically falls back to CPU for operations that fail or when GPU is unavailable
pub struct GpuBackend {
    #[cfg(feature = "gpu")]
    kmer_counter: Option<GpuKmerCounter>,
    #[cfg(feature = "gpu")]
    overlap_finder: Option<GpuOverlapFinder>,
    /// Whether to use the simple (sketch-based) k-mer counting kernel
    use_simple_counting: bool,
    /// Hash table size in bits (e.g., 24 for 16M entries)
    table_size_bits: usize,
    /// Maximum number of overlaps to return
    max_overlaps: usize,
}

impl GpuBackend {
    /// Create a new GPU backend
    ///
    /// # Arguments
    /// * `max_reads` - Maximum expected number of reads
    /// * `max_contigs` - Maximum expected number of contigs
    #[cfg(feature = "gpu")]
    pub fn new(max_reads: usize, max_contigs: usize) -> Result<Self, String> {
        // Default to 24-bit table (16M entries)
        let table_size_bits = 24;
        let max_overlaps = max_contigs * max_contigs / 10; // Expect ~10% overlap rate

        // Try to initialize GPU components
        let kmer_counter = match GpuKmerCounter::new(31, max_reads, table_size_bits) {
            Ok(counter) => {
                tracing::info!("GPU k-mer counter initialized");
                Some(counter)
            }
            Err(e) => {
                tracing::warn!("Failed to initialize GPU k-mer counter: {}", e);
                None
            }
        };

        let overlap_finder = match GpuOverlapFinder::new(max_contigs, max_overlaps) {
            Ok(finder) => {
                tracing::info!("GPU overlap finder initialized");
                Some(finder)
            }
            Err(e) => {
                tracing::warn!("Failed to initialize GPU overlap finder: {}", e);
                None
            }
        };

        if kmer_counter.is_none() && overlap_finder.is_none() {
            return Err("No GPU components could be initialized".to_string());
        }

        Ok(GpuBackend {
            kmer_counter,
            overlap_finder,
            use_simple_counting: false,
            table_size_bits,
            max_overlaps,
        })
    }

    #[cfg(not(feature = "gpu"))]
    pub fn new(_max_reads: usize, _max_contigs: usize) -> Result<Self, String> {
        Err("GPU support not compiled in. Build with --features gpu".to_string())
    }

    /// Configure to use simple sketch-based counting (faster but approximate)
    pub fn with_simple_counting(mut self) -> Self {
        self.use_simple_counting = true;
        self
    }

    /// Check if GPU is actually available for k-mer counting
    #[cfg(feature = "gpu")]
    pub fn has_kmer_gpu(&self) -> bool {
        self.kmer_counter.is_some()
    }

    #[cfg(not(feature = "gpu"))]
    pub fn has_kmer_gpu(&self) -> bool {
        false
    }

    /// Check if GPU is actually available for overlap detection
    #[cfg(feature = "gpu")]
    pub fn has_overlap_gpu(&self) -> bool {
        self.overlap_finder.is_some()
    }

    #[cfg(not(feature = "gpu"))]
    pub fn has_overlap_gpu(&self) -> bool {
        false
    }
}

impl ComputeBackend for GpuBackend {
    #[cfg(feature = "gpu")]
    fn count_kmers(&self, sequences: &[String], k: usize) -> HashMap<String, u32> {
        // Try GPU first
        if let Some(ref counter) = self.kmer_counter {
            // Need to reinitialize counter with correct k if it differs
            match GpuKmerCounter::new(k, sequences.len(), self.table_size_bits) {
                Ok(counter) => {
                    match counter.count(sequences) {
                        Ok(counts) => {
                            tracing::debug!("GPU k-mer counting successful: {} k-mers", counts.len());
                            return counts;
                        }
                        Err(e) => {
                            tracing::warn!("GPU k-mer counting failed, falling back to CPU: {}", e);
                        }
                    }
                }
                Err(e) => {
                    tracing::warn!("Failed to create GPU counter for k={}: {}", k, e);
                }
            }
        }

        // CPU fallback
        tracing::info!("Using CPU fallback for k-mer counting");
        crate::gpu::kmer_gpu::count_kmers_cpu(sequences, k)
    }

    #[cfg(not(feature = "gpu"))]
    fn count_kmers(&self, sequences: &[String], k: usize) -> HashMap<String, u32> {
        crate::gpu::kmer_gpu::count_kmers_cpu(sequences, k)
    }

    #[cfg(feature = "gpu")]
    fn find_overlaps(
        &self,
        contigs: &[String],
        min_overlap: usize,
        max_mismatch: usize,
    ) -> Vec<(usize, usize, usize)> {
        // Try GPU first
        if let Some(ref finder) = self.overlap_finder {
            match finder.find_overlaps(contigs, min_overlap, max_mismatch) {
                Ok(overlaps) => {
                    tracing::debug!("GPU overlap detection successful: {} overlaps", overlaps.len());
                    return overlaps;
                }
                Err(e) => {
                    tracing::warn!("GPU overlap detection failed, falling back to CPU: {}", e);
                }
            }
        }

        // CPU fallback
        tracing::info!("Using CPU fallback for overlap detection");
        crate::gpu::overlap_gpu::find_overlaps_cpu(contigs, min_overlap, max_mismatch)
    }

    #[cfg(not(feature = "gpu"))]
    fn find_overlaps(
        &self,
        contigs: &[String],
        min_overlap: usize,
        max_mismatch: usize,
    ) -> Vec<(usize, usize, usize)> {
        crate::gpu::overlap_gpu::find_overlaps_cpu(contigs, min_overlap, max_mismatch)
    }

    fn build_adjacency(&self, kmer_counts: &HashMap<String, u32>, k: usize) -> AdjacencyTable {
        // Adjacency building is currently CPU-only
        // (GPU would require more complex data structures)
        let mut adjacency = AdjacencyTable::new();
        let bases = ["A", "C", "G", "T"];

        for (kmer, &_count) in kmer_counts {
            if kmer.len() < k {
                continue;
            }

            // Try extensions on suffix (forward edges)
            let suffix = &kmer[1..];
            for base in &bases {
                let next = format!("{}{}", suffix, base);
                if let Some(&next_count) = kmer_counts.get(&next) {
                    adjacency.add_edge(kmer.clone(), next.clone(), next_count);
                }
            }
        }

        adjacency
    }

    fn name(&self) -> &'static str {
        #[cfg(feature = "gpu")]
        {
            "GPU (OpenCL)"
        }
        #[cfg(not(feature = "gpu"))]
        {
            "GPU (not available)"
        }
    }

    fn is_gpu(&self) -> bool {
        #[cfg(feature = "gpu")]
        {
            self.kmer_counter.is_some() || self.overlap_finder.is_some()
        }
        #[cfg(not(feature = "gpu"))]
        {
            false
        }
    }
}

/// Try to create a GPU backend, falling back to None if unavailable
#[cfg(feature = "gpu")]
pub fn try_create_gpu_backend(max_reads: usize, max_contigs: usize) -> Option<GpuBackend> {
    match GpuBackend::new(max_reads, max_contigs) {
        Ok(backend) => {
            if let Some(info) = crate::gpu::kmer_gpu::get_gpu_info() {
                tracing::info!("GPU backend created:\n{}", info);
            }
            Some(backend)
        }
        Err(e) => {
            tracing::warn!("Could not create GPU backend: {}", e);
            None
        }
    }
}

#[cfg(not(feature = "gpu"))]
pub fn try_create_gpu_backend(_max_reads: usize, _max_contigs: usize) -> Option<GpuBackend> {
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[cfg(feature = "gpu")]
    fn test_gpu_backend_creation() {
        // This test may fail if no GPU is available, which is OK
        let result = GpuBackend::new(1000, 100);
        // Just check that it doesn't panic
        match result {
            Ok(_) => println!("GPU backend created successfully"),
            Err(e) => println!("GPU backend creation failed (expected if no GPU): {}", e),
        }
    }

    #[test]
    fn test_try_create_returns_none_without_gpu() {
        // Without GPU feature or hardware, should return None
        #[cfg(not(feature = "gpu"))]
        {
            assert!(try_create_gpu_backend(1000, 100).is_none());
        }
    }
}
