//! Backend abstraction for compute operations
#![allow(dead_code)]

use std::collections::HashMap;
use ahash::AHashMap;

/// Adjacency table for k-mer graph traversal (String-based, legacy)
#[derive(Debug, Clone, Default)]
pub struct AdjacencyTable {
    /// Maps each k-mer to its neighbors with counts
    pub forward: HashMap<String, Vec<(String, u32)>>,
    /// Maps each k-mer to its predecessors with counts
    pub backward: HashMap<String, Vec<(String, u32)>>,
}

impl AdjacencyTable {
    pub fn new() -> Self {
        Self {
            forward: HashMap::new(),
            backward: HashMap::new(),
        }
    }

    pub fn add_edge(&mut self, from: String, to: String, count: u32) {
        self.forward.entry(from.clone()).or_default().push((to.clone(), count));
        self.backward.entry(to).or_default().push((from, count));
    }
}

/// High-performance adjacency table using u64-encoded k-mers.
/// Uses AHashMap for faster hashing of integer keys.
/// Memory: ~16 bytes per edge vs ~100+ bytes for String-based.
#[derive(Debug, Clone)]
pub struct AdjacencyTableU64 {
    /// Maps each k-mer (u64) to its successors with counts
    pub forward: AHashMap<u64, Vec<(u64, u32)>>,
    /// Maps each k-mer (u64) to its predecessors with counts
    pub backward: AHashMap<u64, Vec<(u64, u32)>>,
    /// K-mer size (needed for decoding)
    pub k: u8,
}

impl AdjacencyTableU64 {
    pub fn new(k: u8) -> Self {
        Self {
            forward: AHashMap::new(),
            backward: AHashMap::new(),
            k,
        }
    }

    pub fn with_capacity(k: u8, capacity: usize) -> Self {
        Self {
            forward: AHashMap::with_capacity(capacity),
            backward: AHashMap::with_capacity(capacity),
            k,
        }
    }

    #[inline]
    pub fn add_edge(&mut self, from: u64, to: u64, count: u32) {
        self.forward.entry(from).or_default().push((to, count));
        self.backward.entry(to).or_default().push((from, count));
    }

    /// Get forward neighbors of a k-mer
    #[inline]
    pub fn get_successors(&self, kmer: u64) -> Option<&Vec<(u64, u32)>> {
        self.forward.get(&kmer)
    }

    /// Get backward neighbors of a k-mer
    #[inline]
    pub fn get_predecessors(&self, kmer: u64) -> Option<&Vec<(u64, u32)>> {
        self.backward.get(&kmer)
    }
}

impl Default for AdjacencyTableU64 {
    fn default() -> Self {
        Self::new(31)
    }
}

/// Overlap between two contigs
#[derive(Debug, Clone)]
pub struct Overlap {
    pub from_idx: usize,
    pub to_idx: usize,
    pub overlap_len: usize,
    pub mismatches: usize,
}

/// Trait for compute backends (CPU or GPU)
///
/// This abstraction allows swapping between CPU (Rayon) and GPU (OpenCL)
/// implementations transparently.
pub trait ComputeBackend: Send + Sync {
    /// Count k-mers in the given sequences
    ///
    /// # Arguments
    /// * `sequences` - Input DNA sequences
    /// * `k` - K-mer size
    ///
    /// # Returns
    /// HashMap mapping canonical k-mers to their counts
    fn count_kmers(&self, sequences: &[String], k: usize) -> HashMap<String, u32>;

    /// Find overlaps between contig sequences
    ///
    /// # Arguments
    /// * `contigs` - Contig sequences to compare
    /// * `min_overlap` - Minimum overlap length
    /// * `max_mismatch` - Maximum allowed mismatches
    ///
    /// # Returns
    /// Vector of (from_idx, to_idx, overlap_len) tuples
    fn find_overlaps(&self, contigs: &[String], min_overlap: usize, max_mismatch: usize) -> Vec<(usize, usize, usize)>;

    /// Build k-mer adjacency table for graph assembly
    ///
    /// # Arguments
    /// * `kmer_counts` - K-mer counts from count_kmers()
    /// * `k` - K-mer size
    ///
    /// # Returns
    /// AdjacencyTable with forward/backward edges
    fn build_adjacency(&self, kmer_counts: &HashMap<String, u32>, k: usize) -> AdjacencyTable;

    /// Returns the name of this backend
    fn name(&self) -> &'static str;

    /// Returns whether this backend is GPU-accelerated
    fn is_gpu(&self) -> bool;
}

/// Threshold for preferring CPU over GPU (number of reads)
/// Below this threshold, GPU overhead may dominate
pub const GPU_THRESHOLD_READS: usize = 100_000;

/// Select the best backend based on data size and GPU availability
pub fn select_backend(
    prefer_gpu: bool,
    num_sequences: usize,
    #[cfg(feature = "gpu")] gpu_backend: Option<Box<dyn ComputeBackend>>,
    cpu_backend: Box<dyn ComputeBackend>,
) -> Box<dyn ComputeBackend> {
    #[cfg(feature = "gpu")]
    {
        if prefer_gpu && num_sequences >= GPU_THRESHOLD_READS {
            if let Some(gpu) = gpu_backend {
                tracing::info!("Using GPU backend for {} sequences", num_sequences);
                return gpu;
            }
            tracing::warn!("GPU requested but not available, falling back to CPU");
        } else if prefer_gpu && num_sequences < GPU_THRESHOLD_READS {
            tracing::info!(
                "Dataset too small ({} sequences < {} threshold), using CPU backend",
                num_sequences,
                GPU_THRESHOLD_READS
            );
        }
    }

    #[cfg(not(feature = "gpu"))]
    {
        if prefer_gpu {
            tracing::warn!("GPU support not compiled in, using CPU backend");
        }
        let _ = num_sequences; // suppress unused warning
    }

    cpu_backend
}
