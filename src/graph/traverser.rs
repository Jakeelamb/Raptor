use crate::graph::isoform_graph::IsoformGraph;
use crate::graph::isoform_traverse::TranscriptPath;
use petgraph::visit::EdgeRef;
use rayon::prelude::*;
use std::collections::HashSet;

/// Parallel version of enumerate_paths that scales better for complex graphs
pub fn enumerate_paths_parallel(
    graph: &IsoformGraph,
    start_nodes: &[usize],
    max_depth: usize,
) -> Vec<TranscriptPath> {
    if start_nodes.is_empty() {
        return Vec::new();
    }

    let mut ordered_starts = start_nodes.to_vec();
    ordered_starts.sort_unstable();
    ordered_starts.dedup();

    let mut grouped_paths: Vec<Vec<TranscriptPath>> = ordered_starts
        .par_iter()
        .map(|&start| enumerate_paths_from_start(graph, start, max_depth))
        .collect();

    let total_paths: usize = grouped_paths.iter().map(Vec::len).sum();
    let mut all_paths = Vec::with_capacity(total_paths);
    for group in grouped_paths.iter_mut() {
        all_paths.append(group);
    }
    all_paths
}

fn enumerate_paths_from_start(
    graph: &IsoformGraph,
    start: usize,
    max_depth: usize,
) -> Vec<TranscriptPath> {
    let mut visited = HashSet::new();
    visited.insert(start);
    let mut current_path = vec![start];
    let mut raw_paths = Vec::new();
    dfs_collect_paths(
        graph,
        &mut current_path,
        &mut visited,
        max_depth.max(1),
        &mut raw_paths,
    );

    let mut transcript_paths: Vec<TranscriptPath> = raw_paths
        .into_iter()
        .map(|nodes| {
            let length = nodes.len() * 25;
            TranscriptPath {
                confidence: calculate_path_confidence(graph, &nodes),
                nodes,
                length,
            }
        })
        .collect();
    transcript_paths.sort_unstable_by(|a, b| a.nodes.cmp(&b.nodes));
    transcript_paths
}

fn dfs_collect_paths(
    graph: &IsoformGraph,
    current_path: &mut Vec<usize>,
    visited: &mut HashSet<usize>,
    max_depth: usize,
    out: &mut Vec<Vec<usize>>,
) {
    if current_path.len() >= max_depth {
        out.push(current_path.clone());
        return;
    }

    let current = *current_path
        .last()
        .expect("path should always contain at least one node");
    let mut neighbors: Vec<usize> = graph
        .edges(current)
        .map(|edge| edge.target())
        .filter(|neighbor| !visited.contains(neighbor))
        .collect();
    neighbors.sort_unstable();
    neighbors.dedup();

    if neighbors.is_empty() {
        out.push(current_path.clone());
        return;
    }

    for neighbor in neighbors {
        visited.insert(neighbor);
        current_path.push(neighbor);
        dfs_collect_paths(graph, current_path, visited, max_depth, out);
        current_path.pop();
        visited.remove(&neighbor);
    }
}

/// Batch process paths in parallel for high performance
pub fn process_paths_parallel<F, T>(
    paths: &[TranscriptPath],
    batch_size: usize,
    processor: F,
) -> Vec<T>
where
    F: Fn(&TranscriptPath) -> T + Sync,
    T: Send,
{
    paths
        .par_chunks(batch_size.max(1))
        .flat_map(|chunk| chunk.iter().map(|path| processor(path)).collect::<Vec<_>>())
        .collect()
}

/// Filter a large number of paths in parallel by confidence threshold
pub fn filter_paths_parallel(paths: &[TranscriptPath], min_confidence: f32) -> Vec<TranscriptPath> {
    paths
        .par_iter()
        .filter(|path| path.confidence >= min_confidence)
        .cloned()
        .collect()
}

/// Calculate confidence for a path as the average of edge weights
fn calculate_path_confidence(graph: &IsoformGraph, path: &[usize]) -> f32 {
    if path.len() <= 1 {
        return 1.0; // Single node has perfect confidence
    }

    let mut sum_weights = 0.0;
    let mut count = 0;

    for i in 0..path.len() - 1 {
        if let Some(weight) = graph.edge_weight(path[i], path[i + 1]) {
            sum_weights += *weight;
            count += 1;
        }
    }

    if count == 0 {
        0.0
    } else {
        sum_weights / count as f32
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::graphmap::DiGraphMap;
    use rand::rngs::StdRng;
    use rand::seq::SliceRandom;
    use rand::SeedableRng;
    use rayon::ThreadPoolBuilder;

    fn summarize_paths(paths: &[TranscriptPath]) -> Vec<(Vec<usize>, u32)> {
        paths
            .iter()
            .map(|path| (path.nodes.clone(), path.confidence.to_bits()))
            .collect()
    }

    #[test]
    fn test_enumerate_paths_parallel() {
        // Create a simple test graph
        let mut graph = DiGraphMap::new();

        // Add nodes and edges
        graph.add_node(0);
        graph.add_node(1);
        graph.add_node(2);
        graph.add_node(3);

        graph.add_edge(0, 1, 0.9);
        graph.add_edge(1, 2, 0.8);
        graph.add_edge(1, 3, 0.6);

        // Enumerate paths
        let paths = enumerate_paths_parallel(&graph, &[0], 3);

        // Should find two paths: 0->1->2 and 0->1->3
        assert_eq!(paths.len(), 2);

        // Check if paths are as expected
        let path_strings: Vec<String> = paths.iter().map(|p| format!("{:?}", p.nodes)).collect();

        assert!(path_strings.contains(&"[0, 1, 2]".to_string()));
        assert!(path_strings.contains(&"[0, 1, 3]".to_string()));

        // Check confidence scores
        let path1 = paths.iter().find(|p| p.nodes == vec![0, 1, 2]).unwrap();
        let path2 = paths.iter().find(|p| p.nodes == vec![0, 1, 3]).unwrap();

        assert!((path1.confidence - 0.85).abs() < 0.01); // Average of 0.9 and 0.8
        assert!((path2.confidence - 0.75).abs() < 0.01); // Average of 0.9 and 0.6
    }

    #[test]
    fn test_filter_paths_parallel() {
        // Create test paths
        let paths = vec![
            TranscriptPath {
                nodes: vec![0, 1, 2],
                confidence: 0.9,
                length: 75,
            },
            TranscriptPath {
                nodes: vec![0, 1, 3],
                confidence: 0.7,
                length: 75,
            },
            TranscriptPath {
                nodes: vec![0, 2, 3],
                confidence: 0.3,
                length: 75,
            },
        ];

        // Filter with threshold 0.6
        let filtered = filter_paths_parallel(&paths, 0.6);

        // Should keep first two paths only
        assert_eq!(filtered.len(), 2);
        assert!(filtered.iter().any(|p| p.nodes == vec![0, 1, 2]));
        assert!(filtered.iter().any(|p| p.nodes == vec![0, 1, 3]));
    }

    #[test]
    fn enumerate_paths_parallel_is_stable_under_edge_order_start_permutation_and_thread_count() {
        let nodes: Vec<usize> = (0..7).collect();
        let edges = vec![
            (0, 1, 0.9),
            (0, 2, 0.8),
            (1, 3, 0.85),
            (2, 3, 0.75),
            (2, 4, 0.7),
            (3, 5, 0.95),
            (4, 5, 0.65),
            (4, 6, 0.6),
        ];
        let start_nodes = vec![4, 0, 2, 1];

        let mut graph_a = DiGraphMap::new();
        let mut graph_b = DiGraphMap::new();
        for &node in &nodes {
            graph_a.add_node(node);
            graph_b.add_node(node);
        }
        for &(u, v, w) in &edges {
            graph_a.add_edge(u, v, w);
        }
        for &(u, v, w) in edges.iter().rev() {
            graph_b.add_edge(u, v, w);
        }

        let baseline = ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .expect("single-thread pool should build")
            .install(|| summarize_paths(&enumerate_paths_parallel(&graph_a, &start_nodes, 8)));

        let multithread = ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("multi-thread pool should build")
            .install(|| summarize_paths(&enumerate_paths_parallel(&graph_b, &start_nodes, 8)));

        assert_eq!(multithread, baseline);

        let mut rng = StdRng::seed_from_u64(0x7A9E_11CE);
        for _ in 0..64 {
            let mut shuffled = start_nodes.clone();
            shuffled.shuffle(&mut rng);
            shuffled.extend_from_slice(&[0, 2, 2, 4]);

            let observed = ThreadPoolBuilder::new()
                .num_threads(3)
                .build()
                .expect("thread pool should build")
                .install(|| summarize_paths(&enumerate_paths_parallel(&graph_a, &shuffled, 8)));
            assert_eq!(observed, baseline);
        }
    }
}
