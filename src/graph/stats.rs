use crate::graph::stitch::Path;
use crate::graph::transcript::Transcript;
use petgraph::graphmap::DiGraphMap;
use std::collections::{HashMap, HashSet};

/// Compute the complexity metrics for the graph based on paths
pub struct GraphComplexityMetrics {
    pub total_segments: usize,
    pub total_paths: usize,
    pub shared_segments: usize,
    pub branchiness_percent: f32,
    pub average_path_length: f32,
    pub max_path_length: usize,
    pub min_path_length: usize,
}

fn compute_complexity_from_segments<T, F>(items: &[T], mut segments_of: F) -> GraphComplexityMetrics
where
    F: FnMut(&T) -> &[usize],
{
    let total_paths = items.len();
    let mut all_segments = HashSet::new();
    let mut segment_usage = HashMap::new();
    let mut path_lengths = Vec::with_capacity(items.len());
    let mut unique_segments_in_path = HashSet::new();

    for item in items {
        let segments = segments_of(item);
        path_lengths.push(segments.len());
        unique_segments_in_path.clear();
        for &segment in segments {
            all_segments.insert(segment);
            // Shared-segment metrics should track presence across paths, not
            // multiplicity within the same path.
            if unique_segments_in_path.insert(segment) {
                *segment_usage.entry(segment).or_insert(0) += 1;
            }
        }
    }

    let total_segments = all_segments.len();
    let shared_segments = segment_usage.values().filter(|&&count| count > 1).count();
    let branchiness_percent = if total_segments > 0 {
        (shared_segments as f32 / total_segments as f32) * 100.0
    } else {
        0.0
    };
    let average_path_length = if path_lengths.is_empty() {
        0.0
    } else {
        path_lengths.iter().sum::<usize>() as f32 / path_lengths.len() as f32
    };
    let max_path_length = path_lengths.iter().max().copied().unwrap_or(0);
    let min_path_length = path_lengths.iter().min().copied().unwrap_or(0);

    GraphComplexityMetrics {
        total_segments,
        total_paths,
        shared_segments,
        branchiness_percent,
        average_path_length,
        max_path_length,
        min_path_length,
    }
}

/// Count number of shared segments among multiple isoform paths
pub fn compute_branchiness(paths: &[Vec<String>]) -> usize {
    let mut seg_usage: HashMap<&str, usize> = HashMap::new();
    let mut unique_segments_in_path = HashSet::new();
    for path in paths {
        unique_segments_in_path.clear();
        for seg in path {
            if unique_segments_in_path.insert(seg.as_str()) {
                *seg_usage.entry(seg.as_str()).or_insert(0) += 1;
            }
        }
    }

    seg_usage.values().filter(|&&v| v > 1).count()
}

/// Compute path complexity metrics from a collection of paths
pub fn compute_path_complexity(paths: &[Path]) -> GraphComplexityMetrics {
    compute_complexity_from_segments(paths, |path| path.segments.as_slice())
}

/// Compute transcript-based complexity metrics
pub fn compute_transcript_complexity(transcripts: &[Transcript]) -> GraphComplexityMetrics {
    compute_complexity_from_segments(transcripts, |transcript| transcript.path.as_slice())
}

/// Format graph complexity metrics as a human-readable string
pub fn format_complexity_metrics(metrics: &GraphComplexityMetrics) -> String {
    let shared_percent = if metrics.total_segments > 0 {
        (metrics.shared_segments as f32 / metrics.total_segments as f32 * 100.0) as usize
    } else {
        0
    };

    format!(
        "Graph Complexity Analysis:\n\
         Total segments: {}\n\
         Total paths: {}\n\
         Shared segments: {} ({}%)\n\
         Graph branchiness: {:.1}%\n\
         Average path length: {:.2} segments\n\
         Path length range: {} to {} segments",
        metrics.total_segments,
        metrics.total_paths,
        metrics.shared_segments,
        shared_percent,
        metrics.branchiness_percent,
        metrics.average_path_length,
        metrics.min_path_length,
        metrics.max_path_length
    )
}

/// Compute in-degree and out-degree distributions for nodes in a graph
pub fn compute_degree_distribution<N, E>(
    graph: &DiGraphMap<N, E>,
) -> (HashMap<usize, usize>, HashMap<usize, usize>)
where
    N: petgraph::graphmap::NodeTrait,
{
    let mut in_degree_dist = HashMap::new();
    let mut out_degree_dist = HashMap::new();

    for node in graph.nodes() {
        let in_degree = graph
            .neighbors_directed(node, petgraph::Direction::Incoming)
            .count();
        let out_degree = graph
            .neighbors_directed(node, petgraph::Direction::Outgoing)
            .count();

        *in_degree_dist.entry(in_degree).or_insert(0) += 1;
        *out_degree_dist.entry(out_degree).or_insert(0) += 1;
    }

    (in_degree_dist, out_degree_dist)
}

/// Calculate bubble count (alternative paths between nodes)
pub fn count_bubbles<N, E>(graph: &DiGraphMap<N, E>) -> usize
where
    N: petgraph::graphmap::NodeTrait + std::hash::Hash + Eq + Copy,
{
    let mut bubble_count = 0;

    for source in graph.nodes() {
        let successors: Vec<_> = graph
            .neighbors_directed(source, petgraph::Direction::Outgoing)
            .collect();

        if successors.len() >= 2 {
            // Find all paths from each successor to see if they reconverge
            for i in 0..successors.len() {
                for j in i + 1..successors.len() {
                    let mut visited = HashSet::new();
                    let has_common_target =
                        find_common_successor(graph, successors[i], successors[j], &mut visited);

                    if has_common_target {
                        bubble_count += 1;
                    }
                }
            }
        }
    }

    bubble_count
}

/// Helper function to find if two nodes have a common successor
fn find_common_successor<N, E>(
    graph: &DiGraphMap<N, E>,
    node1: N,
    node2: N,
    visited: &mut HashSet<N>,
) -> bool
where
    N: petgraph::graphmap::NodeTrait + std::hash::Hash + Eq + Copy,
{
    // Get successors of both nodes
    let successors1: HashSet<_> = graph
        .neighbors_directed(node1, petgraph::Direction::Outgoing)
        .collect();
    let successors2: HashSet<_> = graph
        .neighbors_directed(node2, petgraph::Direction::Outgoing)
        .collect();

    // Check for direct common successor
    if !successors1.is_disjoint(&successors2) {
        return true;
    }

    // Mark both nodes as visited
    visited.insert(node1);
    visited.insert(node2);

    // Recursively check further successors (with cycle detection)
    for succ1 in successors1 {
        if visited.contains(&succ1) {
            continue;
        }

        for succ2 in &successors2 {
            if visited.contains(succ2) {
                continue;
            }

            if find_common_successor(graph, succ1, *succ2, visited) {
                return true;
            }
        }
    }

    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_branchiness() {
        let paths = vec![
            vec!["A".to_string(), "B".to_string(), "C".to_string()],
            vec!["A".to_string(), "D".to_string(), "E".to_string()],
            vec!["F".to_string(), "G".to_string(), "H".to_string()],
        ];

        let branchy = compute_branchiness(&paths);
        assert_eq!(branchy, 1); // Only segment "A" is shared
    }

    #[test]
    fn test_compute_branchiness_ignores_within_path_repetitions() {
        let paths = vec![
            vec!["A".to_string(), "A".to_string(), "B".to_string()],
            vec!["C".to_string(), "D".to_string()],
        ];

        let branchy = compute_branchiness(&paths);
        assert_eq!(branchy, 0);
    }

    #[test]
    fn test_compute_path_complexity() {
        let path1 = Path {
            id: 0,
            segments: vec![0, 1, 2],
            overlaps: vec![5, 5],
        };

        let path2 = Path {
            id: 1,
            segments: vec![0, 3, 4],
            overlaps: vec![5, 5],
        };

        let path3 = Path {
            id: 2,
            segments: vec![5, 6, 7],
            overlaps: vec![5, 5],
        };

        let paths = vec![path1, path2, path3];

        let metrics = compute_path_complexity(&paths);
        assert_eq!(metrics.total_paths, 3);
        assert_eq!(metrics.total_segments, 8); // segments 0-7
        assert_eq!(metrics.shared_segments, 1); // Only segment 0 is shared
        assert_eq!(metrics.average_path_length, 3.0);
    }

    #[test]
    fn test_compute_path_complexity_counts_shared_segments_by_path_presence() {
        let path1 = Path {
            id: 0,
            segments: vec![0, 0, 1],
            overlaps: vec![5, 5],
        };
        let path2 = Path {
            id: 1,
            segments: vec![2, 3],
            overlaps: vec![5],
        };

        let metrics = compute_path_complexity(&[path1, path2]);
        assert_eq!(metrics.total_segments, 4);
        assert_eq!(metrics.shared_segments, 0);
    }

    #[test]
    fn test_compute_transcript_complexity_counts_shared_segments_by_transcript_presence() {
        let transcripts = vec![
            Transcript::new(0, "AAAA".to_string(), vec![7, 7, 8], 0.9),
            Transcript::new(1, "CCCC".to_string(), vec![9, 10], 0.8),
        ];

        let metrics = compute_transcript_complexity(&transcripts);
        assert_eq!(metrics.total_segments, 4);
        assert_eq!(metrics.shared_segments, 0);
    }

    #[test]
    fn test_compute_transcript_complexity_matches_path_complexity_semantics() {
        let transcripts = vec![
            Transcript::new(0, "AAAA".to_string(), vec![0, 1, 2], 0.9),
            Transcript::new(1, "CCCC".to_string(), vec![0, 3, 4], 0.8),
            Transcript::new(2, "GGGG".to_string(), vec![5, 6, 7], 0.7),
        ];

        let metrics = compute_transcript_complexity(&transcripts);
        assert_eq!(metrics.total_paths, 3);
        assert_eq!(metrics.total_segments, 8);
        assert_eq!(metrics.shared_segments, 1);
        assert_eq!(metrics.average_path_length, 3.0);
        assert_eq!(metrics.min_path_length, 3);
        assert_eq!(metrics.max_path_length, 3);
    }

    #[test]
    fn test_format_complexity_metrics_handles_zero_segments() {
        let metrics = GraphComplexityMetrics {
            total_segments: 0,
            total_paths: 0,
            shared_segments: 0,
            branchiness_percent: 0.0,
            average_path_length: 0.0,
            max_path_length: 0,
            min_path_length: 0,
        };

        let rendered = format_complexity_metrics(&metrics);
        assert!(rendered.contains("Shared segments: 0 (0%)"));
    }

    #[test]
    fn test_compute_degree_distribution() {
        let mut graph = DiGraphMap::<i32, ()>::new();

        // Create a simple graph
        graph.add_edge(1, 2, ());
        graph.add_edge(1, 3, ());
        graph.add_edge(2, 4, ());
        graph.add_edge(3, 4, ());

        let (in_degree, out_degree) = compute_degree_distribution(&graph);

        // Node 1 has out-degree 2, nodes 2 and 3 have out-degree 1, node 4 has out-degree 0
        assert_eq!(*out_degree.get(&0).unwrap_or(&0), 1); // 1 node with out-degree 0
        assert_eq!(*out_degree.get(&1).unwrap_or(&0), 2); // 2 nodes with out-degree 1
        assert_eq!(*out_degree.get(&2).unwrap_or(&0), 1); // 1 node with out-degree 2

        // Node 1 has in-degree 0, nodes 2 and 3 have in-degree 1, node 4 has in-degree 2
        assert_eq!(*in_degree.get(&0).unwrap_or(&0), 1); // 1 node with in-degree 0
        assert_eq!(*in_degree.get(&1).unwrap_or(&0), 2); // 2 nodes with in-degree 1
        assert_eq!(*in_degree.get(&2).unwrap_or(&0), 1); // 1 node with in-degree 2
    }
}
