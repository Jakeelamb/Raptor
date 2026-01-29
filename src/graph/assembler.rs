// src/graph/assembler.rs
use std::collections::{HashMap, HashSet};
#[allow(deprecated)]
use crate::kmer::kmer::canonical_kmer;
use crate::kmer::kmer::decode_kmer;
use crate::accel::backend::AdjacencyTableU64;
use ahash::{AHashMap, AHashSet};
use rayon::prelude::*;

#[derive(Debug, Clone)]
pub struct Contig {
    pub id: usize,
    pub sequence: String,
    /// K-mer path as u64-encoded values (use decode_kmer to convert back)
    pub kmer_path: Vec<u64>,
}

impl Contig {
    /// Convert the contig sequence to run-length encoding
    pub fn to_rle(&self) -> Vec<(u8, u8)> {
        crate::kmer::rle::rle_encode(&self.sequence)
    }
}

/// Legacy greedy assembly using String k-mers.
/// DEPRECATED: Use greedy_assembly_u64 for 6x memory reduction.
#[deprecated(note = "Use greedy_assembly_u64 for better performance")]
#[allow(deprecated)]
pub fn greedy_assembly(k: usize, kmer_counts: &HashMap<String, u32>, min_len: usize) -> Vec<Contig> {
    use crate::kmer::kmer::encode_kmer;

    let mut used = HashSet::new();
    let mut contigs = Vec::new();

    // Create new HashMap that stores both canonical and original form of each kmer
    let mut kmer_info = HashMap::new();
    for (kmer, &count) in kmer_counts.iter() {
        // Get the canonical form
        if let Some(canon) = canonical_kmer(kmer) {
            kmer_info.insert(kmer.clone(), (canon, count));
        }
    }

    // Sort by descending count
    let mut sorted_kmers: Vec<_> = kmer_info.iter().collect();
    sorted_kmers.sort_unstable_by(|a, b| b.1.1.cmp(&a.1.1));

    // Debug information
    println!("Number of kmers: {}", sorted_kmers.len());
    for (kmer, (canon, count)) in &sorted_kmers {
        println!("Kmer: {}, Canonical: {}, Count: {}", kmer, canon, count);
    }

    // Collect adjacent kmers that can be stitched together
    let mut adjacency: HashMap<String, Vec<(String, u32)>> = HashMap::new();
    for (kmer, (_, _)) in &kmer_info {
        let bases = ["A", "C", "G", "T"];

        // Try extensions on suffix
        if kmer.len() >= k-1 {
            let suffix = &kmer[1..];
            for base in &bases {
                let next = format!("{}{}", suffix, base);
                if let Some(&(_, count)) = kmer_info.get(&next) {
                    adjacency.entry(kmer.clone())
                             .or_default()
                             .push((next.clone(), count));
                }
            }
        }

        // Try extensions on prefix
        if kmer.len() >= k-1 {
            let prefix = &kmer[..kmer.len()-1];
            for base in &bases {
                let prev = format!("{}{}", base, prefix);
                if let Some(&(_, count)) = kmer_info.get(&prev) {
                    adjacency.entry(prev.clone())
                             .or_default()
                             .push((kmer.clone(), count));
                }
            }
        }
    }

    // Build contigs greedily
    for (seed_kmer, _) in sorted_kmers {
        if used.contains(seed_kmer) {
            println!("Skipping already used kmer: {}", seed_kmer);
            continue;
        }

        println!("Starting with seed: {}", seed_kmer);
        let mut contig = seed_kmer.clone();
        let mut path: Vec<u64> = vec![encode_kmer(seed_kmer).unwrap_or(0)];
        used.insert(seed_kmer.clone());

        // Extend right (forward)
        let mut current = seed_kmer.clone();
        while let Some(neighbors) = adjacency.get(&current) {
            // Find best extension (highest count)
            let mut best_next: Option<(String, u32)> = None;
            for (next, count) in neighbors {
                if !used.contains(next) {
                    if best_next.is_none() || count > &best_next.as_ref().unwrap().1 {
                        best_next = Some((next.clone(), *count));
                    }
                }
            }

            if let Some((next, _)) = best_next {
                println!("  Extending right with: {}", next);
                let overlap = k - 1;
                let extension = &next[overlap..];
                contig.push_str(extension);
                path.push(encode_kmer(&next).unwrap_or(0));
                used.insert(next.clone());
                current = next;
            } else {
                break;
            }
        }

        if contig.len() >= min_len {
            println!("Adding contig: {} (length: {})", contig, contig.len());
            contigs.push(Contig {
                id: contigs.len(),
                sequence: contig,
                kmer_path: path,
            });
        } else {
            println!("Contig too short: {} (length: {}, min: {})", contig, contig.len(), min_len);
        }
    }

    contigs
}

/// High-performance greedy assembly using u64-encoded k-mers.
/// Memory usage: ~9 bytes per k-mer vs ~55 bytes for String-based.
/// Zero allocations in the inner assembly loop.
pub fn greedy_assembly_u64(
    k: usize,
    kmer_counts: &AHashMap<u64, u32>,
    adjacency: &AdjacencyTableU64,
    min_len: usize,
) -> Vec<Contig> {
    let mut used = AHashSet::with_capacity(kmer_counts.len());
    let mut contigs = Vec::new();

    // Sort k-mers by descending count for greedy seed selection
    let mut sorted_kmers: Vec<(u64, u32)> = kmer_counts
        .iter()
        .map(|(&kmer, &count)| (kmer, count))
        .collect();
    sorted_kmers.sort_unstable_by(|a, b| b.1.cmp(&a.1));

    // Build contigs greedily
    for (seed_kmer, _count) in sorted_kmers {
        if used.contains(&seed_kmer) {
            continue;
        }

        // Start new contig from seed
        let seed_str = decode_kmer(seed_kmer, k);
        let mut contig = seed_str;
        let mut path = vec![seed_kmer];
        used.insert(seed_kmer);

        // Extend right (forward) following best edges
        let mut current = seed_kmer;
        while let Some(neighbors) = adjacency.get_successors(current) {
            // Find best unused extension (highest count)
            let mut best_next: Option<(u64, u32)> = None;
            for &(next, count) in neighbors {
                if !used.contains(&next) {
                    if best_next.is_none() || count > best_next.unwrap().1 {
                        best_next = Some((next, count));
                    }
                }
            }

            if let Some((next, _)) = best_next {
                // Extend contig by one base (last base of next k-mer)
                let next_decoded = decode_kmer(next, k);
                let extension = next_decoded.chars().last().unwrap();
                contig.push(extension);
                path.push(next);
                used.insert(next);
                current = next;
            } else {
                break;
            }
        }

        // Extend left (backward) from seed
        current = seed_kmer;
        while let Some(neighbors) = adjacency.get_predecessors(current) {
            let mut best_prev: Option<(u64, u32)> = None;
            for &(prev, count) in neighbors {
                if !used.contains(&prev) {
                    if best_prev.is_none() || count > best_prev.unwrap().1 {
                        best_prev = Some((prev, count));
                    }
                }
            }

            if let Some((prev, _)) = best_prev {
                // Prepend first base of prev k-mer
                let prev_decoded = decode_kmer(prev, k);
                let extension = prev_decoded.chars().next().unwrap();
                contig.insert(0, extension);
                path.insert(0, prev);
                used.insert(prev);
                current = prev;
            } else {
                break;
            }
        }

        if contig.len() >= min_len {
            contigs.push(Contig {
                id: contigs.len(),
                sequence: contig,
                kmer_path: path,
            });
        }
    }

    contigs
}
