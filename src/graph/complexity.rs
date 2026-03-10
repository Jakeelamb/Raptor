use crate::eval::metrics::evaluate_lengths_in_place;
use petgraph::graphmap::DiGraphMap;
use petgraph::visit::EdgeRef;
use petgraph::Graph;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use tracing::info;

/// Statistics about paths in a graph
#[derive(Debug, Clone)]
pub struct PathStats {
    /// Total number of paths in the graph
    pub total_paths: usize,

    /// Average path length in segments
    pub average_length: f64,

    /// Median path length in segments
    pub median_length: f64,

    /// N10 of path lengths in segments
    pub path_n10: usize,

    /// N50 of path lengths in segments
    pub path_n50: usize,

    /// N90 of path lengths in segments
    pub path_n90: usize,

    /// N95 of path lengths in segments
    pub path_n95: usize,

    /// N99 of path lengths in segments
    pub path_n99: usize,

    /// L10 of path lengths in segments
    pub path_l10: usize,

    /// L50 of path lengths in segments
    pub path_l50: usize,

    /// L90 of path lengths in segments
    pub path_l90: usize,

    /// L95 of path lengths in segments
    pub path_l95: usize,

    /// L99 of path lengths in segments
    pub path_l99: usize,

    /// auN of path lengths in segments
    pub path_au_n: f64,

    /// Number of shared segments (nodes with multiple incoming/outgoing edges)
    pub branch_count: usize,

    /// Maximum depth of the graph
    pub max_depth: usize,

    /// Number of bubbles in the graph
    pub bubble_count: usize,

    /// Branchiness score (ratio of branches to total segments)
    pub branchiness: f64,
}

/// Compute statistics about paths in a GFA file
pub fn compute_path_stats(gfa_path: &str) -> Result<PathStats, std::io::Error> {
    info!("Computing complexity metrics for {}", gfa_path);

    // Read GFA file
    let file = File::open(gfa_path)?;
    let reader = BufReader::new(file);

    // Parse segments and links
    let mut segments = BTreeSet::new();
    let mut links = BTreeSet::new();
    let mut path_count = 0usize;
    let mut path_lengths = Vec::new();
    let mut node_path_count: HashMap<String, usize> = HashMap::new();
    let mut unique_nodes_in_path = HashSet::new();

    for line_result in reader.lines() {
        let line = line_result?;
        let mut fields = line.split('\t');
        let Some(record_type) = fields.next() else {
            continue;
        };

        match record_type {
            "S" => {
                if let Some(id) = fields.next() {
                    segments.insert(id.to_string());
                }
            }
            "L" => {
                if let (Some(from), Some(to)) = (fields.next(), fields.nth(1)) {
                    let from = from.to_string();
                    let to = to.to_string();
                    // Keep graph node accounting correct even when GFA omits explicit S records.
                    segments.insert(from.clone());
                    segments.insert(to.clone());
                    links.insert((from, to));
                }
            }
            "P" => {
                // P\tname\tsegment+,segment-\t...
                if let Some(segments_field) = fields.nth(1) {
                    add_path_record_from_p(
                        segments_field,
                        &mut path_count,
                        &mut path_lengths,
                        &mut node_path_count,
                        &mut unique_nodes_in_path,
                    );
                }
            }
            "W" => {
                // W\tsample\thap\tseqid\tstart\tend\t>seg1<seg2...
                if let Some(walk) = fields.nth(5) {
                    add_path_record_from_w(
                        walk,
                        &mut path_count,
                        &mut path_lengths,
                        &mut node_path_count,
                        &mut unique_nodes_in_path,
                    );
                }
            }
            _ => {}
        }
    }

    // Some GFAs (or hand-edited subsets) include paths without declaring all segments in S lines.
    // Preserve those referenced nodes so branch/depth/bubble metrics stay accurate.
    segments.extend(node_path_count.keys().cloned());

    // Build graph from links
    let mut graph = Graph::<String, ()>::new();
    let mut node_indices = HashMap::new();

    // Add all segments as nodes
    for seg in &segments {
        let idx = graph.add_node(seg.clone());
        node_indices.insert(seg.clone(), idx);
    }

    // Add all links as edges
    for (from, to) in links.iter() {
        if let (Some(&from_idx), Some(&to_idx)) = (node_indices.get(from), node_indices.get(to)) {
            graph.add_edge(from_idx, to_idx, ());
        }
    }

    // Calculate statistics

    // Find branch points (nodes with multiple incoming or outgoing edges)
    let mut branch_nodes = HashSet::new();

    // For each node, check if it has multiple incoming or outgoing edges
    for node_idx in graph.node_indices() {
        let in_count = graph
            .neighbors_directed(node_idx, petgraph::Direction::Incoming)
            .count();
        let out_count = graph
            .neighbors_directed(node_idx, petgraph::Direction::Outgoing)
            .count();

        if in_count > 1 || out_count > 1 {
            branch_nodes.insert(node_idx);
        }
    }

    // Nodes that appear in multiple paths are also branch points
    for (node, count) in node_path_count {
        if count > 1 {
            if let Some(&idx) = node_indices.get(&node) {
                branch_nodes.insert(idx);
            }
        }
    }

    // Calculate maximum depth of the graph using DFS
    let max_depth = calculate_max_depth(&graph);

    // Build a directed graph using indices
    // Use usize indices instead of String for node identifiers
    let mut digraph = DiGraphMap::<usize, ()>::new();

    // Add nodes using numeric indices
    let mut string_to_idx = HashMap::new();
    for (i, segment) in segments.iter().enumerate() {
        string_to_idx.insert(segment.clone(), i);
        digraph.add_node(i);
    }

    // Add edges using the numeric indices
    for (from, to) in links.iter() {
        if let (Some(&from_idx), Some(&to_idx)) = (string_to_idx.get(from), string_to_idx.get(to)) {
            digraph.add_edge(from_idx, to_idx, ());
        }
    }

    // Count bubbles (nodes with multiple paths that converge)
    let bubble_count = count_bubbles_simple(&digraph);

    // Calculate path-length distribution metrics.
    let path_length_stats = evaluate_lengths_in_place(&mut path_lengths);

    // Calculate branchiness
    let branchiness = if segments.is_empty() {
        0.0
    } else {
        branch_nodes.len() as f64 / segments.len() as f64
    };

    Ok(PathStats {
        total_paths: path_count,
        average_length: path_length_stats.avg_length,
        median_length: path_length_stats.median_length,
        path_n10: path_length_stats.n10,
        path_n50: path_length_stats.n50,
        path_n90: path_length_stats.n90,
        path_n95: path_length_stats.n95,
        path_n99: path_length_stats.n99,
        path_l10: path_length_stats.l10,
        path_l50: path_length_stats.l50,
        path_l90: path_length_stats.l90,
        path_l95: path_length_stats.l95,
        path_l99: path_length_stats.l99,
        path_au_n: path_length_stats.au_n,
        branch_count: branch_nodes.len(),
        max_depth,
        bubble_count,
        branchiness,
    })
}

#[inline]
fn add_path_record_from_p(
    segments_field: &str,
    path_count: &mut usize,
    path_lengths: &mut Vec<usize>,
    node_path_count: &mut HashMap<String, usize>,
    unique_nodes_in_path: &mut HashSet<String>,
) {
    unique_nodes_in_path.clear();
    let mut path_len = 0usize;

    for segment in segments_field.split(',') {
        let segment = segment.trim_end_matches(|c| c == '+' || c == '-');
        if segment.is_empty() {
            continue;
        }
        path_len += 1;
        unique_nodes_in_path.insert(segment.to_string());
    }

    finalize_path_record(
        path_len,
        unique_nodes_in_path,
        path_count,
        path_lengths,
        node_path_count,
    );
}

#[inline]
fn add_path_record_from_w(
    walk: &str,
    path_count: &mut usize,
    path_lengths: &mut Vec<usize>,
    node_path_count: &mut HashMap<String, usize>,
    unique_nodes_in_path: &mut HashSet<String>,
) {
    unique_nodes_in_path.clear();
    let mut path_len = 0usize;

    let bytes = walk.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        if bytes[i] != b'>' && bytes[i] != b'<' {
            i += 1;
            continue;
        }
        i += 1;
        let start = i;
        while i < bytes.len() && bytes[i] != b'>' && bytes[i] != b'<' {
            i += 1;
        }
        if start < i {
            path_len += 1;
            unique_nodes_in_path.insert(walk[start..i].to_string());
        }
    }

    finalize_path_record(
        path_len,
        unique_nodes_in_path,
        path_count,
        path_lengths,
        node_path_count,
    );
}

#[inline]
fn finalize_path_record(
    path_len: usize,
    unique_nodes_in_path: &mut HashSet<String>,
    path_count: &mut usize,
    path_lengths: &mut Vec<usize>,
    node_path_count: &mut HashMap<String, usize>,
) {
    if path_len == 0 {
        unique_nodes_in_path.clear();
        return;
    }

    *path_count += 1;
    path_lengths.push(path_len);
    for node in unique_nodes_in_path.drain() {
        *node_path_count.entry(node).or_insert(0) += 1;
    }
}

/// Calculate the maximum depth of the graph using DFS
fn calculate_max_depth<N, E>(graph: &Graph<N, E>) -> usize {
    if graph.node_count() == 0 {
        return 0;
    }

    // Condense SCCs first so cycles are handled deterministically.
    let sccs = petgraph::algo::kosaraju_scc(graph);
    if sccs.is_empty() {
        return 0;
    }

    let mut node_to_scc = HashMap::with_capacity(graph.node_count());
    let mut scc_weight = vec![0usize; sccs.len()];
    for (scc_idx, component) in sccs.iter().enumerate() {
        scc_weight[scc_idx] = component.len();
        for &node in component {
            node_to_scc.insert(node, scc_idx);
        }
    }

    let mut dag = DiGraphMap::<usize, ()>::new();
    for scc_idx in 0..sccs.len() {
        dag.add_node(scc_idx);
    }
    for edge in graph.edge_references() {
        let from_scc = node_to_scc[&edge.source()];
        let to_scc = node_to_scc[&edge.target()];
        if from_scc != to_scc {
            dag.add_edge(from_scc, to_scc, ());
        }
    }

    let topo = match petgraph::algo::toposort(&dag, None) {
        Ok(order) => order,
        Err(_) => return 0,
    };

    let mut max_nodes_to = scc_weight.clone();
    for scc in topo {
        let best_here = max_nodes_to[scc];
        for next in dag.neighbors_directed(scc, petgraph::Direction::Outgoing) {
            let candidate = best_here + scc_weight[next];
            if candidate > max_nodes_to[next] {
                max_nodes_to[next] = candidate;
            }
        }
    }

    max_nodes_to
        .into_iter()
        .max()
        .unwrap_or(0)
        .saturating_sub(1)
}

/// A simple bubble counting function
fn count_bubbles_simple<N, E>(graph: &DiGraphMap<N, E>) -> usize
where
    N: petgraph::graphmap::NodeTrait + std::hash::Hash + Eq + Copy,
{
    let mut bubble_count = 0;
    let mut reachable_from_primary = HashSet::new();
    let mut visited_secondary = HashSet::new();
    let mut stack = Vec::new();

    // For each node with multiple outgoing edges (potential bubble start)
    for node in graph.nodes() {
        let out_neighbors: Vec<_> = graph
            .neighbors_directed(node, petgraph::Direction::Outgoing)
            .collect();

        if out_neighbors.len() < 2 {
            continue;
        }

        // For each pair of alternative paths
        for i in 0..out_neighbors.len() {
            collect_reachable_outgoing(
                graph,
                out_neighbors[i],
                &mut reachable_from_primary,
                &mut stack,
            );

            for j in i + 1..out_neighbors.len() {
                if path_intersects_reachable(
                    graph,
                    out_neighbors[j],
                    &reachable_from_primary,
                    &mut visited_secondary,
                    &mut stack,
                ) {
                    bubble_count += 1;
                }
            }
        }
    }

    bubble_count
}

/// Collect all nodes reachable from `start` by following outgoing edges.
fn collect_reachable_outgoing<N, E>(
    graph: &DiGraphMap<N, E>,
    start: N,
    reachable: &mut HashSet<N>,
    stack: &mut Vec<N>,
) where
    N: petgraph::graphmap::NodeTrait + std::hash::Hash + Eq + Copy,
{
    reachable.clear();
    stack.clear();
    stack.push(start);

    while let Some(node) = stack.pop() {
        if !reachable.insert(node) {
            continue;
        }

        for neighbor in graph.neighbors_directed(node, petgraph::Direction::Outgoing) {
            if !reachable.contains(&neighbor) {
                stack.push(neighbor);
            }
        }
    }
}

/// Check whether the outgoing traversal from `start` intersects `reachable`.
fn path_intersects_reachable<N, E>(
    graph: &DiGraphMap<N, E>,
    start: N,
    reachable: &HashSet<N>,
    visited: &mut HashSet<N>,
    stack: &mut Vec<N>,
) -> bool
where
    N: petgraph::graphmap::NodeTrait + std::hash::Hash + Eq + Copy,
{
    visited.clear();
    stack.clear();
    stack.push(start);

    while let Some(node) = stack.pop() {
        if !visited.insert(node) {
            continue;
        }
        if reachable.contains(&node) {
            return true;
        }

        for neighbor in graph.neighbors_directed(node, petgraph::Direction::Outgoing) {
            if !visited.contains(&neighbor) {
                stack.push(neighbor);
            }
        }
    }

    false
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_compute_path_stats() {
        // Create a temporary GFA file for testing
        let mut temp_file = NamedTempFile::new().unwrap();

        // Write some basic GFA content to the temp file
        writeln!(temp_file, "H\tVN:Z:1.0").unwrap();
        writeln!(temp_file, "S\t1\tACGT").unwrap();
        writeln!(temp_file, "S\t2\tGCTA").unwrap();
        writeln!(temp_file, "S\t3\tTAGC").unwrap();
        writeln!(temp_file, "L\t1\t+\t2\t+\t3M").unwrap();
        writeln!(temp_file, "L\t2\t+\t3\t+\t3M").unwrap();
        writeln!(temp_file, "P\tpath1\t1+,2+,3+\t*").unwrap();
        writeln!(temp_file, "P\tpath2\t1+,3+\t*").unwrap();
        writeln!(temp_file, "P\tpath3\t2+,3+\t*").unwrap();

        // Compute path stats on the temporary file
        let stats = compute_path_stats(temp_file.path().to_str().unwrap()).unwrap();

        // Verify stats
        assert_eq!(stats.total_paths, 3);
        assert!((stats.average_length - (7.0 / 3.0)).abs() < 1e-12);
        assert_eq!(stats.median_length, 2.0);
        assert_eq!(stats.path_n10, 3);
        assert_eq!(stats.path_n50, 2);
        assert_eq!(stats.path_n90, 2);
        assert_eq!(stats.path_n95, 2);
        assert_eq!(stats.path_n99, 2);
        assert_eq!(stats.path_l10, 1);
        assert_eq!(stats.path_l50, 2);
        assert_eq!(stats.path_l90, 3);
        assert_eq!(stats.path_l95, 3);
        assert_eq!(stats.path_l99, 3);
        assert!((stats.path_au_n - (17.0 / 7.0)).abs() < 1e-12);
        // The branch count should be 3 because:
        // - Node 1 appears in two paths (path1, path2)
        // - Node 2 appears in two paths (path1, path3)
        // - Node 3 appears in all three paths
        assert_eq!(stats.branch_count, 3);
    }

    #[test]
    fn test_compute_path_stats_linear_graph_has_zero_branches() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "H\tVN:Z:1.0").unwrap();
        writeln!(temp_file, "S\t1\tAAAA").unwrap();
        writeln!(temp_file, "S\t2\tCCCC").unwrap();
        writeln!(temp_file, "S\t3\tGGGG").unwrap();
        writeln!(temp_file, "L\t1\t+\t2\t+\t3M").unwrap();
        writeln!(temp_file, "L\t2\t+\t3\t+\t3M").unwrap();
        writeln!(temp_file, "P\tlinear\t1+,2+,3+\t*").unwrap();

        let stats = compute_path_stats(temp_file.path().to_str().unwrap()).unwrap();
        assert_eq!(stats.total_paths, 1);
        assert_eq!(stats.average_length, 3.0);
        assert_eq!(stats.median_length, 3.0);
        assert_eq!(stats.path_n10, 3);
        assert_eq!(stats.path_n50, 3);
        assert_eq!(stats.path_n90, 3);
        assert_eq!(stats.path_n95, 3);
        assert_eq!(stats.path_n99, 3);
        assert_eq!(stats.path_l10, 1);
        assert_eq!(stats.path_l50, 1);
        assert_eq!(stats.path_l90, 1);
        assert_eq!(stats.path_l95, 1);
        assert_eq!(stats.path_l99, 1);
        assert_eq!(stats.path_au_n, 3.0);
        assert_eq!(stats.branch_count, 0);
        assert_eq!(stats.max_depth, 2);
        assert_eq!(stats.bubble_count, 0);
    }

    #[test]
    fn test_compute_path_stats_ignores_duplicate_links() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "H\tVN:Z:1.0").unwrap();
        writeln!(temp_file, "S\t1\tAAAA").unwrap();
        writeln!(temp_file, "S\t2\tCCCC").unwrap();
        writeln!(temp_file, "L\t1\t+\t2\t+\t3M").unwrap();
        writeln!(temp_file, "L\t1\t+\t2\t+\t3M").unwrap();
        writeln!(temp_file, "L\t1\t+\t2\t+\t3M").unwrap();
        writeln!(temp_file, "P\tonly\t1+,2+\t*").unwrap();

        let stats = compute_path_stats(temp_file.path().to_str().unwrap()).unwrap();
        assert_eq!(stats.total_paths, 1);
        assert_eq!(stats.average_length, 2.0);
        assert_eq!(stats.median_length, 2.0);
        assert_eq!(stats.path_n10, 2);
        assert_eq!(stats.path_n50, 2);
        assert_eq!(stats.path_n90, 2);
        assert_eq!(stats.path_n95, 2);
        assert_eq!(stats.path_n99, 2);
        assert_eq!(stats.path_l10, 1);
        assert_eq!(stats.path_l50, 1);
        assert_eq!(stats.path_l90, 1);
        assert_eq!(stats.path_l95, 1);
        assert_eq!(stats.path_l99, 1);
        assert_eq!(stats.path_au_n, 2.0);
        assert_eq!(stats.branch_count, 0);
        assert_eq!(stats.max_depth, 1);
        assert_eq!(stats.bubble_count, 0);
    }

    #[test]
    fn test_compute_path_stats_ignores_empty_path_records() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "H\tVN:Z:1.0").unwrap();
        writeln!(temp_file, "S\t1\tAAAA").unwrap();
        writeln!(temp_file, "S\t2\tCCCC").unwrap();
        writeln!(temp_file, "P\tempty\t\t*").unwrap();
        writeln!(temp_file, "P\tvalid\t1+,2+\t*").unwrap();

        let stats = compute_path_stats(temp_file.path().to_str().unwrap()).unwrap();
        assert_eq!(stats.total_paths, 1);
        assert_eq!(stats.average_length, 2.0);
        assert_eq!(stats.median_length, 2.0);
        assert_eq!(stats.path_n10, 2);
        assert_eq!(stats.path_n50, 2);
        assert_eq!(stats.path_n90, 2);
        assert_eq!(stats.path_n95, 2);
        assert_eq!(stats.path_n99, 2);
        assert_eq!(stats.path_l10, 1);
        assert_eq!(stats.path_l50, 1);
        assert_eq!(stats.path_l90, 1);
        assert_eq!(stats.path_l95, 1);
        assert_eq!(stats.path_l99, 1);
        assert_eq!(stats.path_au_n, 2.0);
    }

    #[test]
    fn test_compute_path_stats_counts_path_node_presence_not_repetitions() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "H\tVN:Z:1.0").unwrap();
        writeln!(temp_file, "S\t1\tAAAA").unwrap();
        writeln!(temp_file, "S\t2\tCCCC").unwrap();
        writeln!(temp_file, "S\t3\tGGGG").unwrap();
        writeln!(temp_file, "P\tpath1\t1+,1+,2+\t*").unwrap();
        writeln!(temp_file, "P\tpath2\t1+,3+\t*").unwrap();

        let stats = compute_path_stats(temp_file.path().to_str().unwrap()).unwrap();
        assert_eq!(stats.total_paths, 2);
        assert_eq!(stats.branch_count, 1);
    }

    #[test]
    fn test_compute_path_stats_supports_walk_records() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "H\tVN:Z:1.0").unwrap();
        writeln!(temp_file, "S\t1\tAAAA").unwrap();
        writeln!(temp_file, "S\t2\tCCCC").unwrap();
        writeln!(temp_file, "S\t3\tGGGG").unwrap();
        writeln!(temp_file, "L\t1\t+\t2\t+\t3M").unwrap();
        writeln!(temp_file, "L\t2\t+\t3\t+\t3M").unwrap();
        writeln!(temp_file, "W\tsample\t0\tchr1\t0\t12\t>1>2>3").unwrap();
        writeln!(temp_file, "W\tsample\t1\tchr1\t0\t8\t>1>3").unwrap();

        let stats = compute_path_stats(temp_file.path().to_str().unwrap()).unwrap();
        assert_eq!(stats.total_paths, 2);
        assert!((stats.average_length - 2.5).abs() < 1e-12);
        assert_eq!(stats.median_length, 2.5);
        assert_eq!(stats.path_n10, 3);
        assert_eq!(stats.path_n50, 3);
        assert_eq!(stats.path_n90, 2);
        assert_eq!(stats.path_n95, 2);
        assert_eq!(stats.path_n99, 2);
        assert_eq!(stats.path_l10, 1);
        assert_eq!(stats.path_l50, 1);
        assert_eq!(stats.path_l90, 2);
        assert_eq!(stats.path_l95, 2);
        assert_eq!(stats.path_l99, 2);
        assert!((stats.path_au_n - 2.6).abs() < 1e-12);
        assert_eq!(stats.branch_count, 2);
    }

    #[test]
    fn test_compute_path_stats_counts_both_p_and_w_records() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "H\tVN:Z:1.0").unwrap();
        writeln!(temp_file, "S\t1\tAAAA").unwrap();
        writeln!(temp_file, "S\t2\tCCCC").unwrap();
        writeln!(temp_file, "S\t3\tGGGG").unwrap();
        writeln!(temp_file, "P\tp1\t1+,2+,3+\t*").unwrap();
        writeln!(temp_file, "W\tsample\t0\tchr1\t0\t8\t>1>3").unwrap();

        let stats = compute_path_stats(temp_file.path().to_str().unwrap()).unwrap();
        assert_eq!(stats.total_paths, 2);
        assert_eq!(stats.path_n10, 3);
        assert_eq!(stats.path_n50, 3);
        assert_eq!(stats.path_n90, 2);
        assert_eq!(stats.path_n95, 2);
        assert_eq!(stats.path_n99, 2);
        assert_eq!(stats.path_l10, 1);
        assert_eq!(stats.path_l50, 1);
        assert_eq!(stats.path_l90, 2);
        assert_eq!(stats.path_l95, 2);
        assert_eq!(stats.path_l99, 2);
        assert!((stats.path_au_n - 2.6).abs() < 1e-12);
        assert_eq!(stats.branch_count, 2);
    }

    #[test]
    fn test_compute_path_stats_counts_shared_path_nodes_without_segment_records() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "H\tVN:Z:1.0").unwrap();
        writeln!(temp_file, "P\tp1\t1+,2+\t*").unwrap();
        writeln!(temp_file, "P\tp2\t1+,3+\t*").unwrap();

        let stats = compute_path_stats(temp_file.path().to_str().unwrap()).unwrap();
        assert_eq!(stats.total_paths, 2);
        assert_eq!(stats.path_n10, 2);
        assert_eq!(stats.path_n50, 2);
        assert_eq!(stats.path_l10, 1);
        assert_eq!(stats.path_l50, 1);
        assert_eq!(stats.branch_count, 1);
        assert!((stats.branchiness - (1.0 / 3.0)).abs() < 1e-12);
    }

    #[test]
    fn test_compute_path_stats_infers_link_nodes_without_segment_records() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "H\tVN:Z:1.0").unwrap();
        writeln!(temp_file, "L\t10\t+\t11\t+\t1M").unwrap();
        writeln!(temp_file, "L\t11\t+\t12\t+\t1M").unwrap();

        let stats = compute_path_stats(temp_file.path().to_str().unwrap()).unwrap();
        assert_eq!(stats.total_paths, 0);
        assert_eq!(stats.branch_count, 0);
        assert_eq!(stats.max_depth, 2);
        assert_eq!(stats.bubble_count, 0);
    }

    #[test]
    fn test_calculate_max_depth_handles_reconverging_paths() {
        let mut graph = Graph::<(), ()>::new();
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let n3 = graph.add_node(());
        let n4 = graph.add_node(());
        let n5 = graph.add_node(());

        graph.add_edge(n1, n2, ());
        graph.add_edge(n1, n3, ());
        graph.add_edge(n2, n4, ());
        graph.add_edge(n3, n4, ());
        graph.add_edge(n4, n5, ());

        assert_eq!(calculate_max_depth(&graph), 3);
    }

    #[test]
    fn test_calculate_max_depth_cycle_with_alternate_branch_is_deterministic() {
        let mut graph = Graph::<(), ()>::new();
        let n0 = graph.add_node(());
        let n1 = graph.add_node(());
        let n2 = graph.add_node(());
        let n3 = graph.add_node(());
        let n4 = graph.add_node(());

        graph.add_edge(n4, n0, ());
        graph.add_edge(n0, n1, ());
        graph.add_edge(n0, n2, ());
        graph.add_edge(n1, n3, ());
        graph.add_edge(n2, n3, ());
        graph.add_edge(n3, n1, ());

        assert_eq!(calculate_max_depth(&graph), 4);
    }

    #[test]
    fn test_count_bubbles_simple_detects_reconvergence_beyond_depth_ten() {
        let mut graph = DiGraphMap::<usize, ()>::new();
        // Bubble start with two outgoing choices.
        graph.add_edge(0, 1, ());
        graph.add_edge(0, 2, ());

        // First branch reaches node 100 after 11 hops.
        let mut first = 1usize;
        for next in 10usize..20usize {
            graph.add_edge(first, next, ());
            first = next;
        }
        graph.add_edge(first, 100, ());

        // Second branch reaches the same node 100 after 11 hops.
        let mut second = 2usize;
        for next in 200usize..210usize {
            graph.add_edge(second, next, ());
            second = next;
        }
        graph.add_edge(second, 100, ());

        assert_eq!(count_bubbles_simple(&graph), 1);
    }

    #[test]
    fn test_count_bubbles_simple_handles_cycles_without_recursion() {
        let mut graph = DiGraphMap::<usize, ()>::new();
        graph.add_edge(0, 1, ());
        graph.add_edge(0, 2, ());
        graph.add_edge(1, 3, ());
        graph.add_edge(3, 1, ());
        graph.add_edge(2, 4, ());
        graph.add_edge(4, 2, ());
        graph.add_edge(3, 5, ());
        graph.add_edge(4, 5, ());

        assert!(count_bubbles_simple(&graph) >= 1);
    }
}
