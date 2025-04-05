// src/graph/assembler.rs
use std::collections::{HashMap, HashSet};
use crate::kmer::kmer::canonical_kmer;
use rayon::prelude::*;

#[derive(Debug, Clone)]
pub struct Contig {
    pub id: usize,
    pub sequence: String,
    pub kmer_path: Vec<String>,
}

impl Contig {
    /// Convert the contig sequence to run-length encoding
    pub fn to_rle(&self) -> Vec<(u8, u8)> {
        crate::kmer::rle::rle_encode(&self.sequence)
    }
}

pub fn greedy_assembly(k: usize, kmer_counts: &HashMap<String, u32>, min_len: usize) -> Vec<Contig> {
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
        let mut path = vec![seed_kmer.clone()];
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
                path.push(next.clone());
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
