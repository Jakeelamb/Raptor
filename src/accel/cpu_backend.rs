use crate::accel::backend::{AdjacencyTable, ComputeBackend};
use crate::accel::simd::match_kmers_with_overlap;
use crate::kmer::kmer::canonical_kmer;
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

impl ComputeBackend for CpuBackend {
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

    fn find_overlaps(
        &self,
        contigs: &[String],
        min_overlap: usize,
        max_mismatch: usize,
    ) -> Vec<(usize, usize, usize)> {
        // Generate all pairs (excluding self-comparisons)
        let contig_pairs: Vec<(usize, usize)> = (0..contigs.len())
            .flat_map(|i| (0..contigs.len()).map(move |j| (i, j)))
            .filter(|&(i, j)| i != j)
            .collect();

        // Parallel overlap detection
        contig_pairs
            .par_iter()
            .filter_map(|&(i, j)| {
                let from = &contigs[i];
                let to = &contigs[j];

                match_kmers_with_overlap(from, to, min_overlap, max_mismatch).map(|(shift, _)| {
                    let overlap_len = from.len() - shift;
                    (i, j, overlap_len)
                })
            })
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
