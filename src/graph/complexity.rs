use petgraph::graph::NodeIndex;
use petgraph::graphmap::DiGraphMap;
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
    let mut links = Vec::new();
    let mut paths = Vec::new();

    for line_result in reader.lines() {
        let line = line_result?;
        let parts: Vec<&str> = line.split('\t').collect();

        match parts.get(0) {
            Some(&"S") if parts.len() >= 3 => {
                // Segment line
                let id = parts[1].to_string();
                segments.insert(id);
            }
            Some(&"L") if parts.len() >= 5 => {
                // Link line
                let from = parts[1].to_string();
                let to = parts[3].to_string();
                links.push((from, to));
            }
            Some(&"P") if parts.len() >= 3 => {
                // Path line
                let path_segments = parts[2].to_string();
                let path_segments: Vec<String> = path_segments
                    .split(',')
                    .map(|s| s.trim_end_matches(|c| c == '+' || c == '-').to_string())
                    .collect();

                paths.push(path_segments);
            }
            _ => {}
        }
    }

    // Build graph from links
    let mut graph = Graph::<String, ()>::new();
    let mut node_indices = HashMap::new();

    // Add all segments as nodes
    for seg in &segments {
        let idx = graph.add_node(seg.clone());
        node_indices.insert(seg.clone(), idx);
    }

    // Add all links as edges
    for (from, to) in &links {
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

    // Also count nodes that appear in multiple different paths
    let mut node_path_count = HashMap::new();
    for path in &paths {
        let mut path_nodes = HashSet::new();
        for segment in path {
            path_nodes.insert(segment);
        }

        for node in path_nodes {
            *node_path_count.entry(node.clone()).or_insert(0) += 1;
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
    for (from, to) in &links {
        if let (Some(&from_idx), Some(&to_idx)) = (string_to_idx.get(from), string_to_idx.get(to)) {
            digraph.add_edge(from_idx, to_idx, ());
        }
    }

    // Count bubbles (nodes with multiple paths that converge)
    let bubble_count = count_bubbles_simple(&digraph);

    // Calculate average path length
    let avg_length = if paths.is_empty() {
        0.0
    } else {
        paths.iter().map(|p| p.len()).sum::<usize>() as f64 / paths.len() as f64
    };

    // Calculate branchiness
    let branchiness = if segments.is_empty() {
        0.0
    } else {
        branch_nodes.len() as f64 / segments.len() as f64
    };

    Ok(PathStats {
        total_paths: paths.len(),
        average_length: avg_length,
        branch_count: branch_nodes.len(),
        max_depth,
        bubble_count,
        branchiness,
    })
}

/// Calculate the maximum depth of the graph using DFS
fn calculate_max_depth<N, E>(graph: &Graph<N, E>) -> usize {
    if graph.node_count() == 0 {
        return 0;
    }

    // Find nodes with no incoming edges (start nodes)
    let start_nodes: Vec<NodeIndex> = graph
        .node_indices()
        .filter(|&n| {
            graph
                .neighbors_directed(n, petgraph::Direction::Incoming)
                .count()
                == 0
        })
        .collect();

    // If no start nodes, find the node with most outgoing edges
    let nodes_to_check: Vec<NodeIndex> = if start_nodes.is_empty() {
        graph.node_indices().collect::<Vec<_>>()
    } else {
        start_nodes
    };

    // Memoized DFS over outgoing edges, with cycle-safe recursion stack.
    let mut memo: HashMap<NodeIndex, usize> = HashMap::with_capacity(graph.node_count());
    let mut on_path: HashSet<NodeIndex> = HashSet::with_capacity(graph.node_count());
    nodes_to_check
        .into_iter()
        .map(|start| dfs_max_depth(graph, start, &mut memo, &mut on_path).saturating_sub(1))
        .max()
        .unwrap_or(0)
}

/// Recursive DFS to find maximum depth
fn dfs_max_depth<N, E>(
    graph: &Graph<N, E>,
    current: NodeIndex,
    memo: &mut HashMap<NodeIndex, usize>,
    on_path: &mut HashSet<NodeIndex>,
) -> usize {
    if let Some(&cached) = memo.get(&current) {
        return cached;
    }

    // Back-edge in a cycle: stop path growth at this point.
    if !on_path.insert(current) {
        return 0;
    }

    let mut max_nodes = 1usize;
    for neighbor in graph.neighbors_directed(current, petgraph::Direction::Outgoing) {
        let path_nodes = 1 + dfs_max_depth(graph, neighbor, memo, on_path);
        max_nodes = max_nodes.max(path_nodes);
    }

    on_path.remove(&current);
    memo.insert(current, max_nodes);
    max_nodes
}

/// A simple bubble counting function
fn count_bubbles_simple<N, E>(graph: &DiGraphMap<N, E>) -> usize
where
    N: petgraph::graphmap::NodeTrait + std::hash::Hash + Eq + Copy,
{
    let mut bubble_count = 0;

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
            for j in i + 1..out_neighbors.len() {
                let path1_start = out_neighbors[i];
                let path2_start = out_neighbors[j];

                // Check if these two paths have a common successor (convergence point)
                let has_common_successor =
                    paths_have_common_target(graph, path1_start, path2_start, 10);

                if has_common_successor {
                    bubble_count += 1;
                }
            }
        }
    }

    bubble_count
}

/// Check if two nodes have a common successor within a certain depth
fn paths_have_common_target<N, E>(
    graph: &DiGraphMap<N, E>,
    node1: N,
    node2: N,
    max_depth: usize,
) -> bool
where
    N: petgraph::graphmap::NodeTrait + std::hash::Hash + Eq + Copy,
{
    // Get all successors of node1 up to max_depth
    let mut visited1 = HashSet::new();
    let mut queue1 = vec![(node1, 1)];

    while let Some((node, depth)) = queue1.pop() {
        visited1.insert(node);

        if depth < max_depth {
            for neighbor in graph.neighbors_directed(node, petgraph::Direction::Outgoing) {
                if !visited1.contains(&neighbor) {
                    queue1.push((neighbor, depth + 1));
                }
            }
        }
    }

    // Get all successors of node2 up to max_depth
    let mut queue2 = vec![(node2, 1)];
    let mut visited2 = HashSet::new();

    while let Some((node, depth)) = queue2.pop() {
        // If this node is also a successor of node1, we found a common target
        if visited1.contains(&node) {
            return true;
        }

        visited2.insert(node);

        if depth < max_depth {
            for neighbor in graph.neighbors_directed(node, petgraph::Direction::Outgoing) {
                if !visited2.contains(&neighbor) {
                    queue2.push((neighbor, depth + 1));
                }
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
}
