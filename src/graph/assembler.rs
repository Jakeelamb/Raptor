// src/graph/assembler.rs
use crate::accel::backend::AdjacencyTableU64;
#[allow(deprecated)]
use crate::kmer::kmer::canonical_kmer;
use crate::kmer::kmer::decode_kmer;
use ahash::{AHashMap, AHashSet};
use std::collections::{HashMap, HashSet};

/// Base lookup table for decoding 2-bit encoded nucleotides
/// A=00, C=01, G=10, T=11
const BASE_BYTES: [u8; 4] = [b'A', b'C', b'G', b'T'];

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

    if k == 0 {
        return Vec::new();
    }

    let mut used = HashSet::new();
    let mut contigs = Vec::new();

    // Create new HashMap that stores both canonical and original form of each kmer
    let mut kmer_info = HashMap::new();
    for (kmer, &count) in kmer_counts.iter() {
        if kmer.len() != k {
            continue;
        }
        // Get the canonical form
        if let Some(canon) = canonical_kmer(kmer) {
            kmer_info.insert(kmer.clone(), (canon, count));
        }
    }
    if kmer_info.is_empty() {
        return contigs;
    }

    // Sort by descending count, then lexicographically by k-mer for deterministic ties.
    let mut sorted_kmers: Vec<_> = kmer_info.iter().collect();
    sorted_kmers.sort_unstable_by(|a, b| b.1 .1.cmp(&a.1 .1).then_with(|| a.0.cmp(b.0)));

    // Collect adjacent kmers that can be stitched together
    let mut adjacency: HashMap<String, Vec<(String, u32)>> = HashMap::new();
    for kmer in kmer_info.keys() {
        let suffix = &kmer[1..];
        for base in [b'A', b'C', b'G', b'T'] {
            let mut next = String::with_capacity(k);
            next.push_str(suffix);
            next.push(base as char);
            if let Some(&(_, count)) = kmer_info.get(&next) {
                adjacency
                    .entry(kmer.clone())
                    .or_default()
                    .push((next, count));
            }
        }
    }
    for neighbors in adjacency.values_mut() {
        neighbors.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        neighbors.dedup_by(|a, b| a.0 == b.0);
    }

    // Build contigs greedily
    for (seed_kmer, _) in sorted_kmers {
        if used.contains(seed_kmer) {
            continue;
        }

        let Some(seed_encoded) = encode_kmer(seed_kmer) else {
            continue;
        };
        let mut contig = seed_kmer.clone();
        let mut path: Vec<u64> = vec![seed_encoded];
        used.insert(seed_kmer.clone());

        // Extend right (forward)
        let mut current = seed_kmer.clone();
        while let Some(neighbors) = adjacency.get(&current) {
            if let Some((next, _)) = neighbors.iter().find(|(next, _)| !used.contains(next)) {
                let overlap = k - 1;
                let extension = &next[overlap..];
                contig.push_str(extension);
                if let Some(encoded) = encode_kmer(next) {
                    path.push(encoded);
                } else {
                    break;
                }
                used.insert(next.clone());
                current = next.clone();
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

/// High-performance greedy assembly using u64-encoded k-mers.
/// Memory usage: ~9 bytes per k-mer vs ~55 bytes for String-based.
/// Zero allocations in the inner assembly loop.
pub fn greedy_assembly_u64(
    k: usize,
    kmer_counts: &AHashMap<u64, u32>,
    adjacency: &AdjacencyTableU64,
    min_len: usize,
) -> Vec<Contig> {
    if !(1..=32).contains(&k) || k != adjacency.k as usize {
        return Vec::new();
    }

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
        let mut contig = decode_kmer(seed_kmer, k).into_bytes();
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
                let extension = BASE_BYTES[(next & 0b11) as usize];
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
        let mut left_bases: Vec<u8> = Vec::new();
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
                let extension = BASE_BYTES[((prev >> shift) & 0b11) as usize];
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
            let mut prefixed_contig = Vec::with_capacity(left_bases.len() + contig.len());
            prefixed_contig.extend_from_slice(&left_bases);
            prefixed_contig.extend_from_slice(&contig);
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
                sequence: String::from_utf8(contig)
                    .expect("contig bytes only contain valid DNA characters"),
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

#[inline]
fn sorted_kmers(kmer_counts: &AHashMap<u64, u32>) -> Vec<u64> {
    let mut kmers: Vec<u64> = kmer_counts.keys().copied().collect();
    kmers.sort_unstable();
    kmers
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
    let kmers = sorted_kmers(kmer_counts);
    remove_tips_with_kmer_order(adjacency, kmer_counts, &kmers, k, max_tip_len, min_coverage)
}

fn remove_tips_with_kmer_order(
    adjacency: &mut AdjacencyTableU64,
    kmer_counts: &AHashMap<u64, u32>,
    kmers: &[u64],
    k: usize,
    max_tip_len: usize,
    min_coverage: u32,
) -> usize {
    let mut to_remove = AHashSet::new();

    // Find k-mers that are dead-ends (in-degree=0 or out-degree=0)
    for &kmer in kmers {
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
    let kmers = sorted_kmers(kmer_counts);
    detect_bubbles_with_kmer_order(adjacency, kmer_counts, &kmers, k, max_bubble_len)
}

fn detect_bubbles_with_kmer_order(
    adjacency: &AdjacencyTableU64,
    kmer_counts: &AHashMap<u64, u32>,
    kmers: &[u64],
    k: usize,
    max_bubble_len: usize,
) -> Vec<Bubble> {
    let mut bubbles = Vec::new();

    // Find branching nodes (out-degree > 1)
    for &kmer in kmers {
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
        let mut progressed = false;

        // Extend path 1
        if let Some(succ) = adjacency.get_successors(current1) {
            if let Some((next, _)) = best_unvisited_neighbor(succ, &visited1) {
                progressed = true;
                path1.push(next);
                visited1.insert(next);
                current1 = next;

                // Check if paths reconverge
                if visited2.contains(&next) && is_valid_reconvergence(start, branch1, branch2, next)
                {
                    return Some(Bubble {
                        start,
                        end: next,
                        path1: path_without_end(&path1, next),
                        path2: path_without_end(&path2, next),
                        coverage1: cov1,
                        coverage2: cov2,
                    });
                }
            }
        }

        // Extend path 2
        if let Some(succ) = adjacency.get_successors(current2) {
            if let Some((next, _)) = best_unvisited_neighbor(succ, &visited2) {
                progressed = true;
                path2.push(next);
                visited2.insert(next);
                current2 = next;

                // Check if paths reconverge
                if visited1.contains(&next) && is_valid_reconvergence(start, branch1, branch2, next)
                {
                    return Some(Bubble {
                        start,
                        end: next,
                        path1: path_without_end(&path1, next),
                        path2: path_without_end(&path2, next),
                        coverage1: cov1,
                        coverage2: cov2,
                    });
                }
            }
        }

        if !progressed {
            break;
        }
    }

    None
}

#[inline]
fn is_valid_reconvergence(start: u64, branch1: u64, branch2: u64, node: u64) -> bool {
    node != start && node != branch1 && node != branch2
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

    remove_bubble_path_edges(adjacency, bubble.start, bubble.end, path_to_remove)
}

#[inline]
fn remove_bubble_path_edges(
    adjacency: &mut AdjacencyTableU64,
    start: u64,
    end: u64,
    path: &[u64],
) -> bool {
    let mut removed_any = false;

    if path.is_empty() {
        removed_any |= remove_directed_edge(adjacency, start, end);
        return removed_any;
    }

    let mut previous = start;
    for &node in path {
        removed_any |= remove_directed_edge(adjacency, previous, node);
        previous = node;
    }
    removed_any |= remove_directed_edge(adjacency, previous, end);

    // Remove now-isolated path nodes only; preserve shared nodes used elsewhere.
    for &node in path {
        prune_node_if_isolated(adjacency, node);
    }

    removed_any
}

#[inline]
fn remove_directed_edge(adjacency: &mut AdjacencyTableU64, from: u64, to: u64) -> bool {
    let mut removed = false;

    if let Some(successors) = adjacency.forward.get_mut(&from) {
        let before = successors.len();
        successors.retain(|(next, _)| *next != to);
        removed |= successors.len() != before;
    }
    if adjacency
        .forward
        .get(&from)
        .is_some_and(std::vec::Vec::is_empty)
    {
        adjacency.forward.remove(&from);
    }

    if let Some(predecessors) = adjacency.backward.get_mut(&to) {
        let before = predecessors.len();
        predecessors.retain(|(pred, _)| *pred != from);
        removed |= predecessors.len() != before;
    }
    if adjacency
        .backward
        .get(&to)
        .is_some_and(std::vec::Vec::is_empty)
    {
        adjacency.backward.remove(&to);
    }

    removed
}

#[inline]
fn prune_node_if_isolated(adjacency: &mut AdjacencyTableU64, node: u64) {
    let has_successors = adjacency
        .forward
        .get(&node)
        .is_some_and(|neighbors| !neighbors.is_empty());
    let has_predecessors = adjacency
        .backward
        .get(&node)
        .is_some_and(|neighbors| !neighbors.is_empty());

    if !has_successors && !has_predecessors {
        adjacency.forward.remove(&node);
        adjacency.backward.remove(&node);
    }
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

#[inline]
fn path_without_end(path: &[u64], end: u64) -> Vec<u64> {
    path.iter().copied().take_while(|k| *k != end).collect()
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
    let kmers = sorted_kmers(kmer_counts);

    // Multiple rounds of cleanup
    for _round in 0..3 {
        // Remove tips
        let tips = remove_tips_with_kmer_order(
            adjacency,
            kmer_counts,
            &kmers,
            k,
            max_tip_len,
            min_coverage,
        );
        total_tips += tips;

        // Detect and collapse bubbles
        let bubbles =
            detect_bubbles_with_kmer_order(adjacency, kmer_counts, &kmers, k, max_bubble_len);
        let mut collapsed_this_round = 0usize;
        for bubble in &bubbles {
            if collapse_bubble(adjacency, bubble) {
                total_bubbles += 1;
                collapsed_this_round += 1;
            }
        }

        // Stop early if no cleanup operation made progress.
        if tips == 0 && collapsed_this_round == 0 {
            break;
        }
    }

    (total_tips, total_bubbles)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::kmer::encode_kmer;
    use std::collections::HashMap;

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
    fn greedy_u64_rejects_invalid_k_and_k_mismatch_without_panicking() {
        let mut counts = AHashMap::new();
        counts.insert(encode_kmer("AAA").unwrap(), 10);

        let adjacency_k3 = AdjacencyTableU64::new(3);
        assert!(greedy_assembly_u64(0, &counts, &adjacency_k3, 1).is_empty());
        assert!(greedy_assembly_u64(33, &counts, &adjacency_k3, 1).is_empty());

        let adjacency_k4 = AdjacencyTableU64::new(4);
        assert!(greedy_assembly_u64(3, &counts, &adjacency_k4, 1).is_empty());
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
    #[allow(deprecated)]
    fn greedy_string_tie_breaks_seed_and_neighbor_order_deterministically() {
        let mut counts_a = HashMap::new();
        counts_a.insert("AAG".to_string(), 5);
        counts_a.insert("AAA".to_string(), 10);
        counts_a.insert("AAC".to_string(), 5);

        let mut counts_b = HashMap::new();
        counts_b.insert("AAC".to_string(), 5);
        counts_b.insert("AAA".to_string(), 10);
        counts_b.insert("AAG".to_string(), 5);

        let contigs_a = greedy_assembly(3, &counts_a, 3);
        let contigs_b = greedy_assembly(3, &counts_b, 3);

        let seqs_a: Vec<String> = contigs_a.into_iter().map(|c| c.sequence).collect();
        let seqs_b: Vec<String> = contigs_b.into_iter().map(|c| c.sequence).collect();

        assert_eq!(seqs_a, vec!["AAAC".to_string(), "AAG".to_string()]);
        assert_eq!(seqs_b, seqs_a);
    }

    #[test]
    #[allow(deprecated)]
    fn greedy_string_rejects_invalid_k_and_mismatched_kmers_without_panicking() {
        let mut counts = HashMap::new();
        counts.insert("AAA".to_string(), 10);
        counts.insert("AAC".to_string(), 9);
        counts.insert("AA".to_string(), 20);
        counts.insert("AA?".to_string(), 30);

        assert!(greedy_assembly(0, &counts, 1).is_empty());

        let k4 = greedy_assembly(4, &counts, 1);
        assert!(k4.is_empty());

        let k3 = greedy_assembly(3, &counts, 3);
        let seqs: Vec<String> = k3.into_iter().map(|c| c.sequence).collect();
        assert_eq!(seqs, vec!["AAAC".to_string()]);
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

    #[test]
    fn detect_bubbles_ignores_cycles_that_only_rejoin_the_branch_start() {
        let k = 3;
        let aaa = encode_kmer("AAA").unwrap();
        let aac = encode_kmer("AAC").unwrap();
        let aag = encode_kmer("AAG").unwrap();
        let agt = encode_kmer("AGT").unwrap();

        let mut counts = AHashMap::new();
        counts.insert(aaa, 20);
        counts.insert(aac, 10);
        counts.insert(aag, 10);
        counts.insert(agt, 9);

        let mut adjacency = AdjacencyTableU64::new(k as u8);
        adjacency.add_edge(aaa, aac, 10);
        adjacency.add_edge(aaa, aag, 10);
        adjacency.add_edge(aac, aaa, 10);
        adjacency.add_edge(aag, agt, 9);
        adjacency.add_edge(agt, aaa, 9);

        let bubbles = detect_bubbles(&adjacency, &counts, k, 8);
        assert!(bubbles.is_empty());
    }

    #[test]
    fn collapse_bubble_keeps_reconvergence_node() {
        let k = 3;
        let aaa = encode_kmer("AAA").unwrap();
        let aac = encode_kmer("AAC").unwrap();
        let aag = encode_kmer("AAG").unwrap();
        let acc = encode_kmer("ACC").unwrap();

        let mut counts = AHashMap::new();
        counts.insert(aaa, 12);
        counts.insert(aac, 8);
        counts.insert(aag, 8);
        counts.insert(acc, 10);

        let mut adjacency = AdjacencyTableU64::new(k as u8);
        adjacency.add_edge(aaa, aac, 8);
        adjacency.add_edge(aaa, aag, 8);
        adjacency.add_edge(aac, acc, 8);
        adjacency.add_edge(aag, acc, 8);

        let bubbles = detect_bubbles(&adjacency, &counts, k, 8);
        assert_eq!(bubbles.len(), 1);
        assert_eq!(bubbles[0].end, acc);

        let collapsed = collapse_bubble(&mut adjacency, &bubbles[0]);
        assert!(collapsed);

        // Reconvergence node must remain in the graph after collapsing one branch.
        assert!(adjacency.get_predecessors(acc).is_some());
        assert!(!adjacency.get_predecessors(acc).unwrap().is_empty());
    }

    #[test]
    fn collapse_bubble_preserves_external_edges_of_removed_path_nodes() {
        let k = 3;
        let taa = encode_kmer("TAA").unwrap();
        let aaa = encode_kmer("AAA").unwrap();
        let aac = encode_kmer("AAC").unwrap();
        let aag = encode_kmer("AAG").unwrap();
        let acc = encode_kmer("ACC").unwrap();

        let mut counts = AHashMap::new();
        counts.insert(taa, 7);
        counts.insert(aaa, 20);
        counts.insert(aac, 9);
        counts.insert(aag, 9);
        counts.insert(acc, 20);

        let mut adjacency = AdjacencyTableU64::new(k as u8);
        adjacency.add_edge(taa, aag, 7); // External support to the lower-priority branch node.
        adjacency.add_edge(aaa, aac, 9);
        adjacency.add_edge(aaa, aag, 9);
        adjacency.add_edge(aac, acc, 9);
        adjacency.add_edge(aag, acc, 9);

        let bubbles = detect_bubbles(&adjacency, &counts, k, 8);
        assert_eq!(bubbles.len(), 1);
        assert_eq!(bubbles[0].path1.first().copied(), Some(aac));
        assert_eq!(bubbles[0].path2.first().copied(), Some(aag));

        assert!(collapse_bubble(&mut adjacency, &bubbles[0]));
        assert!(!collapse_bubble(&mut adjacency, &bubbles[0]));

        // Bubble edge is removed.
        let start_successors = adjacency.get_successors(aaa).cloned().unwrap_or_default();
        assert_eq!(start_successors, vec![(aac, 9)]);

        // External edge into removed path node is preserved.
        let taa_successors = adjacency.get_successors(taa).cloned().unwrap_or_default();
        assert_eq!(taa_successors, vec![(aag, 7)]);
        assert!(adjacency.get_predecessors(aag).is_some());

        // Reconvergence is still reachable from surviving path.
        let acc_predecessors = adjacency.get_predecessors(acc).cloned().unwrap_or_default();
        assert_eq!(acc_predecessors, vec![(aac, 9)]);
    }

    #[test]
    fn cleanup_graph_collapses_bubbles_when_no_tips_are_removed() {
        let k = 3;
        let aaa = encode_kmer("AAA").unwrap();
        let aac = encode_kmer("AAC").unwrap();
        let aag = encode_kmer("AAG").unwrap();
        let acc = encode_kmer("ACC").unwrap();

        let mut counts = AHashMap::new();
        counts.insert(aaa, 20);
        counts.insert(aac, 10);
        counts.insert(aag, 10);
        counts.insert(acc, 20);

        let mut adjacency = AdjacencyTableU64::new(k as u8);
        adjacency.add_edge(aaa, aac, 10);
        adjacency.add_edge(aaa, aag, 10);
        adjacency.add_edge(aac, acc, 10);
        adjacency.add_edge(aag, acc, 10);

        let (tips_removed, bubbles_collapsed) = cleanup_graph(&mut adjacency, &counts, k, 3);
        assert_eq!(tips_removed, 0);
        assert_eq!(bubbles_collapsed, 1);

        let predecessors = adjacency.get_predecessors(acc).cloned().unwrap_or_default();
        assert_eq!(predecessors, vec![(aac, 10)]);
    }
}
