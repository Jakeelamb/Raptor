use crate::accel::backend::{AdjacencyTable, AdjacencyTableU64, ComputeBackend};
#[allow(deprecated)]
use crate::kmer::kmer::canonical_kmer;
use crate::kmer::kmer::{KmerU64, canonical_kmer_u64};
use ahash::AHashMap;
use rayon::prelude::*;
use std::collections::HashMap;

/// CPU-based compute backend using Rayon for parallelization
pub struct CpuBackend {
    /// Number of threads to use (0 = auto-detect)
    num_threads: usize,
}

impl CpuBackend {
    pub fn new() -> Self {
        Self { num_threads: 0 }
    }

    pub fn with_threads(num_threads: usize) -> Self {
        if num_threads > 0 {
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build_global()
                .ok(); // Ignore if already initialized
        }
        Self { num_threads }
    }
}

impl Default for CpuBackend {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuBackend {
    /// Count k-mers using u64 encoding for maximum performance.
    /// Uses sliding window encoding to avoid redundant work.
    /// Zero allocations in the inner loop.
    pub fn count_kmers_u64(&self, sequences: &[String], k: usize) -> AHashMap<u64, u32> {
        if k > 32 {
            panic!("k-mer size {} exceeds maximum of 32 for u64 encoding", k);
        }

        // Parallel k-mer counting with thread-local hashmaps
        let thread_counts: Vec<AHashMap<u64, u32>> = sequences
            .par_iter()
            .fold(
                || AHashMap::with_capacity(1024),
                |mut counts, seq| {
                    let bytes = seq.as_bytes();
                    if bytes.len() >= k {
                        // Initialize first k-mer
                        if let Some(mut kmer) = KmerU64::from_slice(&bytes[0..k]) {
                            let canonical = kmer.canonical();
                            *counts.entry(canonical.encoded).or_insert(0) += 1;

                            // Sliding window: extend and re-canonicalize
                            for i in k..bytes.len() {
                                if let Some(next) = kmer.extend(bytes[i]) {
                                    kmer = next;
                                    let canonical = kmer.canonical();
                                    *counts.entry(canonical.encoded).or_insert(0) += 1;
                                } else {
                                    // Invalid base encountered, restart
                                    if i + 1 >= k {
                                        if let Some(fresh) = KmerU64::from_slice(&bytes[i + 1 - k..i + 1]) {
                                            kmer = fresh;
                                            let canonical = kmer.canonical();
                                            *counts.entry(canonical.encoded).or_insert(0) += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    counts
                },
            )
            .collect();

        // Pre-size the merged map
        let total_estimate: usize = thread_counts.iter().map(|m| m.len()).sum();
        let mut merged = AHashMap::with_capacity(total_estimate);

        // Merge thread-local counts
        for local in thread_counts {
            for (kmer, count) in local {
                *merged.entry(kmer).or_insert(0) += count;
            }
        }

        merged
    }

    /// Build adjacency table using u64-encoded k-mers.
    /// Much faster than string-based version due to integer operations.
    pub fn build_adjacency_u64(&self, kmer_counts: &AHashMap<u64, u32>, k: usize) -> AdjacencyTableU64 {
        let mut adjacency = AdjacencyTableU64::with_capacity(k as u8, kmer_counts.len());

        // Mask for k-1 bases
        let suffix_mask: u64 = (1u64 << ((k - 1) * 2)) - 1;

        // Build adjacency by trying all possible extensions
        for (&kmer, &_count) in kmer_counts {
            // Get suffix (last k-1 bases)
            let suffix = kmer & suffix_mask;

            // Try all 4 possible extensions
            for base in 0u64..4 {
                let next = (suffix << 2) | base;
                let next_canonical = canonical_kmer_u64(next, k);

                if let Some(&next_count) = kmer_counts.get(&next_canonical) {
                    adjacency.add_edge(kmer, next_canonical, next_count);
                }
            }
        }

        adjacency
    }
}

impl ComputeBackend for CpuBackend {
    #[allow(deprecated)]
    fn count_kmers(&self, sequences: &[String], k: usize) -> HashMap<String, u32> {
        // Parallel k-mer counting with thread-local hashmaps
        let thread_counts: Vec<HashMap<String, u32>> = sequences
            .par_iter()
            .fold(
                || HashMap::new(),
                |mut counts, seq| {
                    if seq.len() >= k {
                        for i in 0..=seq.len() - k {
                            if let Some(kmer) = canonical_kmer(&seq[i..i + k]) {
                                *counts.entry(kmer).or_insert(0) += 1;
                            }
                        }
                    }
                    counts
                },
            )
            .collect();

        // Merge thread-local counts
        let mut merged = HashMap::new();
        for local in thread_counts {
            for (kmer, count) in local {
                *merged.entry(kmer).or_insert(0) += count;
            }
        }

        merged
    }

    /// Optimized overlap detection using suffix/prefix indexing.
    /// Complexity: O(n*k) instead of O(n^2) for n contigs.
    fn find_overlaps(
        &self,
        contigs: &[String],
        min_overlap: usize,
        max_mismatch: usize,
    ) -> Vec<(usize, usize, usize)> {
        if contigs.is_empty() {
            return Vec::new();
        }

        // Build suffix index: maps suffix -> list of (contig_idx, suffix_len)
        // This allows O(1) lookup of potential overlaps
        let mut suffix_index: AHashMap<String, Vec<(usize, usize)>> = AHashMap::new();

        // Index all suffixes of length >= min_overlap
        for (idx, contig) in contigs.iter().enumerate() {
            let len = contig.len();
            for suffix_len in min_overlap..=len.min(100) {  // Cap at 100 to avoid memory explosion
                let suffix = &contig[len - suffix_len..];
                suffix_index
                    .entry(suffix.to_string())
                    .or_default()
                    .push((idx, suffix_len));
            }
        }

        // Find overlaps by looking up prefixes in suffix index
        let results: Vec<(usize, usize, usize)> = contigs
            .par_iter()
            .enumerate()
            .flat_map(|(to_idx, contig)| {
                let mut overlaps = Vec::new();
                let len = contig.len();

                // Check all prefixes against suffix index
                for prefix_len in min_overlap..=len.min(100) {
                    let prefix = &contig[0..prefix_len];

                    if let Some(matches) = suffix_index.get(prefix) {
                        for &(from_idx, suffix_len) in matches {
                            if from_idx != to_idx && suffix_len == prefix_len {
                                // Verify the overlap with mismatch tolerance
                                let from = &contigs[from_idx];
                                let from_suffix = &from[from.len() - suffix_len..];

                                // Count mismatches
                                let mismatches = from_suffix
                                    .chars()
                                    .zip(prefix.chars())
                                    .filter(|(a, b)| a != b)
                                    .count();

                                if mismatches <= max_mismatch {
                                    overlaps.push((from_idx, to_idx, prefix_len));
                                }
                            }
                        }
                    }
                }

                overlaps
            })
            .collect();

        // Deduplicate (keep longest overlap per pair)
        let mut best_overlaps: AHashMap<(usize, usize), usize> = AHashMap::new();
        for (from, to, len) in results {
            let entry = best_overlaps.entry((from, to)).or_insert(0);
            if len > *entry {
                *entry = len;
            }
        }

        best_overlaps
            .into_iter()
            .map(|((from, to), len)| (from, to, len))
            .collect()
    }

    fn build_adjacency(&self, kmer_counts: &HashMap<String, u32>, k: usize) -> AdjacencyTable {
        let mut adjacency = AdjacencyTable::new();
        let bases = ["A", "C", "G", "T"];

        // Build adjacency by trying all possible extensions
        for (kmer, &count) in kmer_counts {
            if kmer.len() < k {
                continue;
            }

            // Try extensions on suffix (forward edges)
            let suffix = &kmer[1..];
            for base in &bases {
                let next = format!("{}{}", suffix, base);
                if kmer_counts.contains_key(&next) {
                    let next_count = *kmer_counts.get(&next).unwrap_or(&0);
                    adjacency.add_edge(kmer.clone(), next.clone(), next_count);
                }
            }
        }

        // Also check canonical forms for completeness
        for (kmer, &_count) in kmer_counts {
            if kmer.len() < k {
                continue;
            }

            let prefix = &kmer[..kmer.len() - 1];
            for base in &bases {
                let prev = format!("{}{}", base, prefix);
                if kmer_counts.contains_key(&prev) {
                    // This edge should already exist from the forward pass
                    // but ensure it's in the backward map
                    let prev_count = *kmer_counts.get(&prev).unwrap_or(&0);
                    if !adjacency.backward.contains_key(kmer) {
                        adjacency.backward.entry(kmer.clone()).or_default();
                    }
                    let entry = adjacency.backward.get_mut(kmer).unwrap();
                    if !entry.iter().any(|(k, _)| k == &prev) {
                        entry.push((prev.clone(), prev_count));
                    }
                }
            }
        }

        adjacency
    }

    fn name(&self) -> &'static str {
        "CPU (Rayon)"
    }

    fn is_gpu(&self) -> bool {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cpu_kmer_counting() {
        let backend = CpuBackend::new();
        let sequences = vec![
            "ACGTACGT".to_string(),
            "ACGTACGT".to_string(),
            "CGTACGTA".to_string(),
        ];

        let counts = backend.count_kmers(&sequences, 4);

        // ACGT appears at positions 0 and 4 in first two sequences
        assert!(counts.contains_key("ACGT"));
        assert!(*counts.get("ACGT").unwrap() >= 2);
    }

    #[test]
    fn test_cpu_overlap_detection() {
        let backend = CpuBackend::new();
        let contigs = vec![
            "ACGTACGTACGT".to_string(),
            "ACGTACGTAAAA".to_string(),  // Overlaps with first
        ];

        let overlaps = backend.find_overlaps(&contigs, 4, 2);

        // Should find some overlaps
        assert!(!overlaps.is_empty());
    }

    #[test]
    fn test_cpu_adjacency_building() {
        let backend = CpuBackend::new();
        let mut counts = HashMap::new();
        counts.insert("ACGT".to_string(), 10);
        counts.insert("CGTA".to_string(), 8);
        counts.insert("GTAC".to_string(), 5);

        let adj = backend.build_adjacency(&counts, 4);

        // ACGT -> CGTA should exist
        assert!(adj.forward.contains_key("ACGT"));
        let acgt_neighbors = adj.forward.get("ACGT").unwrap();
        assert!(acgt_neighbors.iter().any(|(k, _)| k == "CGTA"));
    }
}
