//! Backend abstraction for compute operations
#![allow(dead_code)]

use ahash::AHashMap;
use std::collections::HashMap;

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
        upsert_edge_string(
            self.forward.entry(from.clone()).or_default(),
            to.clone(),
            count,
        );
        upsert_edge_string(self.backward.entry(to).or_default(), from, count);
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
        upsert_edge_u64(self.forward.entry(from).or_default(), to, count);
        upsert_edge_u64(self.backward.entry(to).or_default(), from, count);
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

#[inline]
fn upsert_edge_u64(edges: &mut Vec<(u64, u32)>, node: u64, count: u32) {
    match edges.binary_search_by_key(&node, |(n, _)| *n) {
        Ok(idx) => {
            edges[idx].1 = edges[idx].1.max(count);
        }
        Err(idx) => {
            edges.insert(idx, (node, count));
        }
    }
}

#[inline]
fn upsert_edge_string(edges: &mut Vec<(String, u32)>, node: String, count: u32) {
    match edges.binary_search_by(|(n, _)| n.as_str().cmp(node.as_str())) {
        Ok(idx) => {
            edges[idx].1 = edges[idx].1.max(count);
        }
        Err(idx) => {
            edges.insert(idx, (node, count));
        }
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
    fn find_overlaps(
        &self,
        contigs: &[String],
        min_overlap: usize,
        max_mismatch: usize,
    ) -> Vec<(usize, usize, usize)>;

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

#[cfg(test)]
mod tests {
    use super::{AdjacencyTable, AdjacencyTableU64};

    #[test]
    fn add_edge_u64_deduplicates_and_uses_max_weight() {
        let mut adjacency = AdjacencyTableU64::new(31);
        adjacency.add_edge(10, 20, 2);
        adjacency.add_edge(10, 20, 7);
        adjacency.add_edge(10, 20, 3);

        assert_eq!(adjacency.get_successors(10), Some(&vec![(20, 7)]));
        assert_eq!(adjacency.get_predecessors(20), Some(&vec![(10, 7)]));
    }

    #[test]
    fn add_edge_u64_has_stable_neighbor_order_independent_of_insert_order() {
        let mut adjacency_a = AdjacencyTableU64::new(31);
        adjacency_a.add_edge(1, 4, 1);
        adjacency_a.add_edge(1, 2, 1);
        adjacency_a.add_edge(1, 3, 1);

        let mut adjacency_b = AdjacencyTableU64::new(31);
        adjacency_b.add_edge(1, 2, 1);
        adjacency_b.add_edge(1, 3, 1);
        adjacency_b.add_edge(1, 4, 1);

        assert_eq!(
            adjacency_a.get_successors(1),
            Some(&vec![(2, 1), (3, 1), (4, 1)])
        );
        assert_eq!(adjacency_a.get_successors(1), adjacency_b.get_successors(1));
    }

    #[test]
    fn add_edge_string_deduplicates_and_uses_max_weight() {
        let mut adjacency = AdjacencyTable::new();
        adjacency.add_edge("AAA".to_string(), "AAT".to_string(), 4);
        adjacency.add_edge("AAA".to_string(), "AAT".to_string(), 9);
        adjacency.add_edge("AAA".to_string(), "AAT".to_string(), 1);

        assert_eq!(
            adjacency.forward.get("AAA"),
            Some(&vec![("AAT".to_string(), 9)])
        );
        assert_eq!(
            adjacency.backward.get("AAT"),
            Some(&vec![("AAA".to_string(), 9)])
        );
    }

    #[test]
    fn add_edge_string_has_stable_neighbor_order_independent_of_insert_order() {
        let mut adjacency_a = AdjacencyTable::new();
        adjacency_a.add_edge("AAA".to_string(), "AAT".to_string(), 1);
        adjacency_a.add_edge("AAA".to_string(), "AAC".to_string(), 1);
        adjacency_a.add_edge("AAA".to_string(), "AAG".to_string(), 1);

        let mut adjacency_b = AdjacencyTable::new();
        adjacency_b.add_edge("AAA".to_string(), "AAC".to_string(), 1);
        adjacency_b.add_edge("AAA".to_string(), "AAG".to_string(), 1);
        adjacency_b.add_edge("AAA".to_string(), "AAT".to_string(), 1);

        assert_eq!(
            adjacency_a.forward.get("AAA"),
            Some(&vec![
                ("AAC".to_string(), 1),
                ("AAG".to_string(), 1),
                ("AAT".to_string(), 1),
            ])
        );
        assert_eq!(adjacency_a.forward.get("AAA"), adjacency_b.forward.get("AAA"));
    }
}
