// src/graph/assembler.rs
use crate::accel::backend::AdjacencyTableU64;
#[allow(deprecated)]
use crate::kmer::kmer::canonical_kmer;
use crate::kmer::kmer::decode_kmer;
use ahash::{AHashMap, AHashSet};
use std::collections::{HashMap, HashSet};

/// Base lookup table for decoding 2-bit encoded nucleotides
/// A=00, C=01, G=10, T=11
const BASES: [char; 4] = ['A', 'C', 'G', 'T'];

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
pub fn greedy_assembly(
    k: usize,
    kmer_counts: &HashMap<String, u32>,
    min_len: usize,
) -> Vec<Contig> {
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
    sorted_kmers.sort_unstable_by(|a, b| b.1 .1.cmp(&a.1 .1));

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
        if kmer.len() >= k - 1 {
            let suffix = &kmer[1..];
            for base in &bases {
                let next = format!("{}{}", suffix, base);
                if let Some(&(_, count)) = kmer_info.get(&next) {
                    adjacency
                        .entry(kmer.clone())
                        .or_default()
                        .push((next.clone(), count));
                }
            }
        }

        // Try extensions on prefix
        if kmer.len() >= k - 1 {
            let prefix = &kmer[..kmer.len() - 1];
            for base in &bases {
                let prev = format!("{}{}", base, prefix);
                if let Some(&(_, count)) = kmer_info.get(&prev) {
                    adjacency
                        .entry(prev.clone())
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
            println!(
                "Contig too short: {} (length: {}, min: {})",
                contig,
                contig.len(),
                min_len
            );
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
    sorted_kmers.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

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
                    match best_next {
                        None => best_next = Some((next, count)),
                        Some((best_kmer, best_count)) => {
                            if count > best_count || (count == best_count && next < best_kmer) {
                                best_next = Some((next, count));
                            }
                        }
                    }
                }
            }

            if let Some((next, _)) = best_next {
                // Extend contig by one base (last base of next k-mer)
                // OPTIMIZED: Extract last 2 bits directly instead of decoding entire k-mer
                let extension = BASES[(next & 0b11) as usize];
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
        let mut left_bases: Vec<char> = Vec::new();
        let mut left_path: Vec<u64> = Vec::new();
        while let Some(neighbors) = adjacency.get_predecessors(current) {
            let mut best_prev: Option<(u64, u32)> = None;
            for &(prev, count) in neighbors {
                if !used.contains(&prev) {
                    match best_prev {
                        None => best_prev = Some((prev, count)),
                        Some((best_kmer, best_count)) => {
                            if count > best_count || (count == best_count && prev < best_kmer) {
                                best_prev = Some((prev, count));
                            }
                        }
                    }
                }
            }

            if let Some((prev, _)) = best_prev {
                // Prepend first base of prev k-mer
                // OPTIMIZED: Extract first base by shifting right by (k-1)*2 bits
                let shift = (k - 1) * 2;
                let extension = BASES[((prev >> shift) & 0b11) as usize];
                left_bases.push(extension);
                left_path.push(prev);
                used.insert(prev);
                current = prev;
            } else {
                break;
            }
        }

        if !left_bases.is_empty() {
            left_bases.reverse();
            let mut prefixed_contig = String::with_capacity(left_bases.len() + contig.len());
            for base in left_bases {
                prefixed_contig.push(base);
            }
            prefixed_contig.push_str(&contig);
            contig = prefixed_contig;
        }

        if !left_path.is_empty() {
            left_path.reverse();
            left_path.extend(path);
            path = left_path;
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

/// Represents a bubble (divergent paths that reconverge) in the assembly graph
#[derive(Debug, Clone)]
pub struct Bubble {
    /// Starting k-mer where paths diverge
    pub start: u64,
    /// Ending k-mer where paths reconverge
    pub end: u64,
    /// The two alternative paths (as k-mer sequences)
    pub path1: Vec<u64>,
    pub path2: Vec<u64>,
    /// Coverage of each path
    pub coverage1: u32,
    pub coverage2: u32,
}

/// Remove tips (dead-end paths) from the assembly graph.
///
/// Tips are short paths with only one connection that likely represent
/// sequencing errors. Removing them reduces noise and improves assembly.
///
/// # Arguments
/// * `adjacency` - The k-mer adjacency graph (modified in place)
/// * `kmer_counts` - K-mer counts for coverage information
/// * `k` - K-mer size
/// * `max_tip_len` - Maximum length of a tip to remove (typically 2*k)
/// * `min_coverage` - Minimum coverage threshold; tips below this are removed
///
/// # Returns
/// Number of tips removed
pub fn remove_tips(
    adjacency: &mut AdjacencyTableU64,
    kmer_counts: &AHashMap<u64, u32>,
    k: usize,
    max_tip_len: usize,
    min_coverage: u32,
) -> usize {
    let mut to_remove = AHashSet::new();

    let mut kmers: Vec<u64> = kmer_counts.keys().copied().collect();
    kmers.sort_unstable();

    // Find k-mers that are dead-ends (in-degree=0 or out-degree=0)
    for kmer in kmers {
        let count = *kmer_counts.get(&kmer).unwrap_or(&0);

        // Skip high-coverage k-mers (probably real)
        if count >= min_coverage || to_remove.contains(&kmer) {
            continue;
        }

        let in_degree = adjacency
            .get_predecessors(kmer)
            .map(|v| v.len())
            .unwrap_or(0);
        let out_degree = adjacency.get_successors(kmer).map(|v| v.len()).unwrap_or(0);

        // Dead-end: only one direction connected
        let direction = if in_degree == 0 && out_degree > 0 {
            Some(TipDirection::Forward)
        } else if in_degree > 0 && out_degree == 0 {
            Some(TipDirection::Backward)
        } else {
            None
        };

        if let Some(direction) = direction {
            if let Some(path) = trace_tip_path(
                adjacency,
                kmer_counts,
                kmer,
                k,
                max_tip_len,
                min_coverage,
                direction,
            ) {
                to_remove.extend(path);
            }
        }
    }

    let mut ordered_to_remove: Vec<u64> = to_remove.into_iter().collect();
    ordered_to_remove.sort_unstable();

    // Remove the tips
    for &kmer in &ordered_to_remove {
        remove_kmer_from_graph(adjacency, kmer);
    }

    ordered_to_remove.len()
}

#[derive(Debug, Clone, Copy)]
enum TipDirection {
    Forward,
    Backward,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TipNeighborChoice {
    None,
    Multiple,
    Single(u64),
}

/// Trace a linear tip path from a dead-end.
///
/// Returns `None` if the path is too long, branches before reaching a junction,
/// or contains high-coverage k-mers that should be retained.
fn trace_tip_path(
    adjacency: &AdjacencyTableU64,
    kmer_counts: &AHashMap<u64, u32>,
    start: u64,
    _k: usize,
    max_len: usize,
    min_coverage: u32,
    direction: TipDirection,
) -> Option<Vec<u64>> {
    let mut current = start;
    let mut path = vec![start];
    let mut visited = AHashSet::new();
    visited.insert(start);

    loop {
        if path.len() > max_len {
            return None;
        }

        let neighbors = match direction {
            TipDirection::Forward => adjacency.get_successors(current),
            TipDirection::Backward => adjacency.get_predecessors(current),
        };

        match unique_unvisited_neighbor(neighbors, &visited) {
            TipNeighborChoice::Single(n) => {
                // Check if this node is a branch point (multiple connections)
                let in_deg = adjacency.get_predecessors(n).map(|v| v.len()).unwrap_or(0);
                let out_deg = adjacency.get_successors(n).map(|v| v.len()).unwrap_or(0);

                if in_deg > 1 || out_deg > 1 {
                    // Reached a junction, stop tracing
                    return Some(path);
                }

                if kmer_counts.get(&n).copied().unwrap_or(u32::MAX) >= min_coverage {
                    return None;
                }

                visited.insert(n);
                current = n;
                path.push(n);
            }
            TipNeighborChoice::None => return Some(path),
            TipNeighborChoice::Multiple => return None,
        }
    }
}

#[inline]
fn unique_unvisited_neighbor(
    neighbors: Option<&Vec<(u64, u32)>>,
    visited: &AHashSet<u64>,
) -> TipNeighborChoice {
    let Some(neighbors) = neighbors else {
        return TipNeighborChoice::None;
    };
    let mut candidate: Option<u64> = None;

    for &(next, _) in neighbors {
        if visited.contains(&next) {
            continue;
        }
        match candidate {
            None => candidate = Some(next),
            Some(_) => return TipNeighborChoice::Multiple,
        }
    }

    match candidate {
        Some(next) => TipNeighborChoice::Single(next),
        None => TipNeighborChoice::None,
    }
}

/// Detect bubbles in the assembly graph.
///
/// Bubbles are divergent paths that reconverge, often caused by SNPs,
/// small indels, or sequencing errors.
///
/// # Arguments
/// * `adjacency` - The k-mer adjacency graph
/// * `kmer_counts` - K-mer counts
/// * `k` - K-mer size
/// * `max_bubble_len` - Maximum bubble length to consider (typically 2*k)
///
/// # Returns
/// Vector of detected bubbles
pub fn detect_bubbles(
    adjacency: &AdjacencyTableU64,
    kmer_counts: &AHashMap<u64, u32>,
    k: usize,
    max_bubble_len: usize,
) -> Vec<Bubble> {
    let mut bubbles = Vec::new();
    let mut kmers: Vec<u64> = kmer_counts.keys().copied().collect();
    kmers.sort_unstable();

    // Find branching nodes (out-degree > 1)
    for kmer in kmers {
        if let Some(successors) = adjacency.get_successors(kmer) {
            if successors.len() >= 2 {
                // Try to find bubbles starting from this branch
                if let Some(bubble) = find_bubble_from_branch(
                    adjacency,
                    kmer_counts,
                    kmer,
                    successors,
                    k,
                    max_bubble_len,
                ) {
                    bubbles.push(bubble);
                }
            }
        }
    }

    bubbles
}

/// Try to find a bubble starting from a branching k-mer.
fn find_bubble_from_branch(
    adjacency: &AdjacencyTableU64,
    _kmer_counts: &AHashMap<u64, u32>,
    start: u64,
    successors: &[(u64, u32)],
    _k: usize,
    max_len: usize,
) -> Option<Bubble> {
    if successors.len() < 2 {
        return None;
    }

    // Take the two highest-coverage branches
    let mut sorted_succ: Vec<_> = successors.to_vec();
    sorted_succ.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    let (branch1, cov1) = sorted_succ[0];
    let (branch2, cov2) = sorted_succ[1];

    // Trace both paths looking for reconvergence
    let mut path1 = vec![branch1];
    let mut path2 = vec![branch2];
    let mut visited1 = AHashSet::new();
    let mut visited2 = AHashSet::new();
    visited1.insert(start);
    visited1.insert(branch1);
    visited2.insert(start);
    visited2.insert(branch2);

    let mut current1 = branch1;
    let mut current2 = branch2;

    for _ in 0..max_len {
        // Extend path 1
        if let Some(succ) = adjacency.get_successors(current1) {
            if let Some((next, _)) = best_unvisited_neighbor(succ, &visited1) {
                path1.push(next);
                visited1.insert(next);
                current1 = next;

                // Check if paths reconverge
                if visited2.contains(&next) {
                    return Some(Bubble {
                        start,
                        end: next,
                        path1,
                        path2: path2.into_iter().take_while(|k| *k != next).collect(),
                        coverage1: cov1,
                        coverage2: cov2,
                    });
                }
            }
        }

        // Extend path 2
        if let Some(succ) = adjacency.get_successors(current2) {
            if let Some((next, _)) = best_unvisited_neighbor(succ, &visited2) {
                path2.push(next);
                visited2.insert(next);
                current2 = next;

                // Check if paths reconverge
                if visited1.contains(&next) {
                    return Some(Bubble {
                        start,
                        end: next,
                        path1: path1.into_iter().take_while(|k| *k != next).collect(),
                        path2,
                        coverage1: cov1,
                        coverage2: cov2,
                    });
                }
            }
        }
    }

    None
}

/// Collapse a bubble by choosing the higher-coverage path.
///
/// # Arguments
/// * `adjacency` - The adjacency graph (modified in place)
/// * `bubble` - The bubble to collapse
///
/// # Returns
/// true if the bubble was collapsed, false otherwise
pub fn collapse_bubble(adjacency: &mut AdjacencyTableU64, bubble: &Bubble) -> bool {
    // Choose the path to remove (lower coverage)
    let path_to_remove = if bubble.coverage1 >= bubble.coverage2 {
        &bubble.path2
    } else {
        &bubble.path1
    };

    // Remove k-mers from the low-coverage path
    for kmer in path_to_remove {
        remove_kmer_from_graph(adjacency, *kmer);
    }

    true
}

#[inline]
fn best_unvisited_neighbor(
    neighbors: &[(u64, u32)],
    visited: &AHashSet<u64>,
) -> Option<(u64, u32)> {
    neighbors
        .iter()
        .copied()
        .filter(|(kmer, _)| !visited.contains(kmer))
        .max_by(|a, b| a.1.cmp(&b.1).then_with(|| b.0.cmp(&a.0)))
}

fn remove_kmer_from_graph(adjacency: &mut AdjacencyTableU64, kmer: u64) {
    if let Some(successors) = adjacency.forward.remove(&kmer) {
        for (next, _) in successors {
            if let Some(predecessors) = adjacency.backward.get_mut(&next) {
                predecessors.retain(|(pred, _)| *pred != kmer);
            }
        }
    }

    if let Some(predecessors) = adjacency.backward.remove(&kmer) {
        for (prev, _) in predecessors {
            if let Some(successors) = adjacency.forward.get_mut(&prev) {
                successors.retain(|(next, _)| *next != kmer);
            }
        }
    }
}

/// Clean up the assembly graph by removing tips and collapsing bubbles.
///
/// This is a high-level function that applies multiple rounds of cleanup.
///
/// # Arguments
/// * `adjacency` - The adjacency graph (modified in place)
/// * `kmer_counts` - K-mer counts
/// * `k` - K-mer size
/// * `min_coverage` - Minimum coverage for keeping a k-mer
///
/// # Returns
/// Tuple of (tips_removed, bubbles_collapsed)
pub fn cleanup_graph(
    adjacency: &mut AdjacencyTableU64,
    kmer_counts: &AHashMap<u64, u32>,
    k: usize,
    min_coverage: u32,
) -> (usize, usize) {
    let max_tip_len = 2 * k;
    let max_bubble_len = 2 * k;

    let mut total_tips = 0;
    let mut total_bubbles = 0;

    // Multiple rounds of cleanup
    for _round in 0..3 {
        // Remove tips
        let tips = remove_tips(adjacency, kmer_counts, k, max_tip_len, min_coverage);
        total_tips += tips;

        if tips == 0 {
            break;
        }

        // Detect and collapse bubbles
        let bubbles = detect_bubbles(adjacency, kmer_counts, k, max_bubble_len);
        for bubble in &bubbles {
            if collapse_bubble(adjacency, bubble) {
                total_bubbles += 1;
            }
        }
    }

    (total_tips, total_bubbles)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::kmer::encode_kmer;

    #[test]
    fn greedy_u64_uses_deterministic_seed_order_when_counts_tie() {
        let k = 3;
        let mut counts = AHashMap::new();
        counts.insert(encode_kmer("AAG").unwrap(), 10);
        counts.insert(encode_kmer("AAA").unwrap(), 10);
        counts.insert(encode_kmer("AAC").unwrap(), 10);

        let adjacency = AdjacencyTableU64::new(k as u8);
        let contigs = greedy_assembly_u64(k, &counts, &adjacency, k);
        let sequences: Vec<&str> = contigs.iter().map(|c| c.sequence.as_str()).collect();

        assert_eq!(sequences, vec!["AAA", "AAC", "AAG"]);
    }

    #[test]
    fn greedy_u64_tie_breaks_neighbors_by_kmer_value() {
        let k = 3;
        let mut counts = AHashMap::new();
        let aaa = encode_kmer("AAA").unwrap();
        let aag = encode_kmer("AAG").unwrap();
        let aac = encode_kmer("AAC").unwrap();
        counts.insert(aaa, 10);
        counts.insert(aag, 5);
        counts.insert(aac, 5);

        let mut adjacency = AdjacencyTableU64::new(k as u8);
        // Insert in reverse lexical order to ensure traversal does not depend on insertion order.
        adjacency.add_edge(aaa, aag, 5);
        adjacency.add_edge(aaa, aac, 5);

        let contigs = greedy_assembly_u64(k, &counts, &adjacency, k);
        assert_eq!(contigs.first().map(|c| c.sequence.as_str()), Some("AAAC"));
    }

    #[test]
    fn greedy_u64_extends_left_and_right_without_path_regression() {
        let k = 3;
        let tga = encode_kmer("TGA").unwrap();
        let gaa = encode_kmer("GAA").unwrap();
        let aaa = encode_kmer("AAA").unwrap();
        let aat = encode_kmer("AAT").unwrap();

        let mut counts = AHashMap::new();
        counts.insert(aaa, 10);
        counts.insert(gaa, 9);
        counts.insert(tga, 9);
        counts.insert(aat, 9);

        let mut adjacency = AdjacencyTableU64::new(k as u8);
        adjacency.add_edge(tga, gaa, 9);
        adjacency.add_edge(gaa, aaa, 10);
        adjacency.add_edge(aaa, aat, 9);

        let contigs = greedy_assembly_u64(k, &counts, &adjacency, k);
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0].sequence, "TGAAAT");
        assert_eq!(contigs[0].kmer_path, vec![tga, gaa, aaa, aat]);
    }

    #[test]
    fn remove_tips_removes_entire_linear_tip_path() {
        let k = 3;
        let taa = encode_kmer("TAA").unwrap();
        let aaa = encode_kmer("AAA").unwrap();
        let aat = encode_kmer("AAT").unwrap();
        let atc = encode_kmer("ATC").unwrap();
        let gtc = encode_kmer("GTC").unwrap();

        let mut counts = AHashMap::new();
        counts.insert(taa, 1);
        counts.insert(aaa, 1);
        counts.insert(aat, 1);
        counts.insert(atc, 10);
        counts.insert(gtc, 10);

        let mut adjacency = AdjacencyTableU64::new(k as u8);
        adjacency.add_edge(taa, aaa, 1);
        adjacency.add_edge(aaa, aat, 1);
        adjacency.add_edge(aat, atc, 1);
        adjacency.add_edge(gtc, atc, 10);

        let removed = remove_tips(&mut adjacency, &counts, k, 8, 3);
        assert_eq!(removed, 3);
        assert!(adjacency.get_successors(taa).is_none());
        assert!(adjacency.get_successors(aaa).is_none());
        assert!(adjacency.get_successors(aat).is_none());

        let predecessors = adjacency.get_predecessors(atc).cloned().unwrap_or_default();
        assert_eq!(predecessors, vec![(gtc, 10)]);
    }

    #[test]
    fn remove_tips_preserves_path_with_high_coverage_internal_node() {
        let k = 3;
        let taa = encode_kmer("TAA").unwrap();
        let aaa = encode_kmer("AAA").unwrap();
        let aat = encode_kmer("AAT").unwrap();
        let atc = encode_kmer("ATC").unwrap();
        let gtc = encode_kmer("GTC").unwrap();

        let mut counts = AHashMap::new();
        counts.insert(taa, 1);
        counts.insert(aaa, 12);
        counts.insert(aat, 1);
        counts.insert(atc, 10);
        counts.insert(gtc, 10);

        let mut adjacency = AdjacencyTableU64::new(k as u8);
        adjacency.add_edge(taa, aaa, 1);
        adjacency.add_edge(aaa, aat, 1);
        adjacency.add_edge(aat, atc, 1);
        adjacency.add_edge(gtc, atc, 10);

        let removed = remove_tips(&mut adjacency, &counts, k, 8, 3);
        assert_eq!(removed, 0);
        let successors = adjacency.get_successors(taa).cloned().unwrap_or_default();
        assert_eq!(successors, vec![(aaa, 1)]);
    }

    #[test]
    fn detect_bubbles_is_deterministic_with_equal_coverage_branches() {
        let k = 3;
        let aaa = encode_kmer("AAA").unwrap();
        let aac = encode_kmer("AAC").unwrap();
        let aag = encode_kmer("AAG").unwrap();
        let acc = encode_kmer("ACC").unwrap();

        let mut counts = AHashMap::new();
        counts.insert(aaa, 20);
        counts.insert(aac, 8);
        counts.insert(aag, 8);
        counts.insert(acc, 20);

        let mut adjacency_a = AdjacencyTableU64::new(k as u8);
        adjacency_a.add_edge(aaa, aag, 8);
        adjacency_a.add_edge(aaa, aac, 8);
        adjacency_a.add_edge(aac, acc, 8);
        adjacency_a.add_edge(aag, acc, 8);

        let mut adjacency_b = AdjacencyTableU64::new(k as u8);
        adjacency_b.add_edge(aaa, aac, 8);
        adjacency_b.add_edge(aaa, aag, 8);
        adjacency_b.add_edge(aag, acc, 8);
        adjacency_b.add_edge(aac, acc, 8);

        let bubbles_a = detect_bubbles(&adjacency_a, &counts, k, 8);
        let bubbles_b = detect_bubbles(&adjacency_b, &counts, k, 8);
        assert_eq!(bubbles_a.len(), 1);
        assert_eq!(bubbles_b.len(), 1);
        assert_eq!(bubbles_a[0].path1, bubbles_b[0].path1);
        assert_eq!(bubbles_a[0].path2, bubbles_b[0].path2);
        assert_eq!(bubbles_a[0].path1.first().copied(), Some(aac));
        assert_eq!(bubbles_a[0].path2.first().copied(), Some(aag));
    }
}
