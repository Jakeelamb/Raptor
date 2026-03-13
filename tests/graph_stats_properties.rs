use petgraph::graphmap::DiGraphMap;
use proptest::prelude::*;
use raptor::graph::stats::{
    compute_branchiness, compute_path_complexity, compute_transcript_complexity, count_bubbles,
};
use raptor::graph::stitch::Path;
use raptor::graph::transcript::Transcript;
use std::collections::{HashMap, HashSet};

fn expected_shared_segments(paths: &[Vec<usize>]) -> usize {
    let mut usage = HashMap::new();
    let mut seen_in_path = HashSet::new();

    for path in paths {
        seen_in_path.clear();
        for &segment in path {
            if seen_in_path.insert(segment) {
                *usage.entry(segment).or_insert(0usize) += 1;
            }
        }
    }

    usage.values().filter(|&&count| count > 1).count()
}

fn expected_total_segments(paths: &[Vec<usize>]) -> usize {
    let mut all_segments = HashSet::new();
    for path in paths {
        for &segment in path {
            all_segments.insert(segment);
        }
    }
    all_segments.len()
}

fn build_graph(node_count: usize, edges: &[(usize, usize)]) -> DiGraphMap<usize, ()> {
    let mut graph = DiGraphMap::<usize, ()>::new();
    for node in 0..node_count {
        graph.add_node(node);
    }
    for &(from, to) in edges {
        if from < node_count && to < node_count && from != to {
            graph.add_edge(from, to, ());
        }
    }
    graph
}

fn reference_bubble_count(graph: &DiGraphMap<usize, ()>) -> usize {
    fn collect_reachable(graph: &DiGraphMap<usize, ()>, start: usize, out: &mut HashSet<usize>) {
        out.clear();
        let mut stack = vec![start];
        while let Some(node) = stack.pop() {
            if !out.insert(node) {
                continue;
            }
            for next in graph.neighbors_directed(node, petgraph::Direction::Outgoing) {
                if !out.contains(&next) {
                    stack.push(next);
                }
            }
        }
    }

    let mut bubbles = 0usize;
    let mut reachable_a = HashSet::new();
    let mut reachable_b = HashSet::new();

    for source in graph.nodes() {
        let successors: Vec<usize> = graph
            .neighbors_directed(source, petgraph::Direction::Outgoing)
            .collect();
        for i in 0..successors.len() {
            collect_reachable(graph, successors[i], &mut reachable_a);
            for &successor in successors.iter().skip(i + 1) {
                collect_reachable(graph, successor, &mut reachable_b);
                if !reachable_a.is_disjoint(&reachable_b) {
                    bubbles += 1;
                }
            }
        }
    }

    bubbles
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    #[test]
    fn compute_branchiness_matches_per_path_presence_even_with_repetitions(
        paths in prop::collection::vec(prop::collection::vec(0usize..24, 0..40), 1..24)
    ) {
        let as_strings: Vec<Vec<String>> = paths
            .iter()
            .map(|path| path.iter().map(|segment| format!("seg_{segment}")).collect())
            .collect();

        let observed = compute_branchiness(&as_strings);
        let expected = expected_shared_segments(&paths);
        prop_assert_eq!(observed, expected);
    }

    #[test]
    fn compute_path_complexity_shared_segments_follow_presence_semantics(
        paths in prop::collection::vec(prop::collection::vec(0usize..24, 0..40), 1..24)
    ) {
        let assembled_paths: Vec<Path> = paths
            .iter()
            .enumerate()
            .map(|(id, segments)| Path {
                id,
                segments: segments.clone(),
                overlaps: vec![0; segments.len().saturating_sub(1)],
            })
            .collect();

        let metrics = compute_path_complexity(&assembled_paths);
        prop_assert_eq!(metrics.total_paths, paths.len());
        prop_assert_eq!(metrics.total_segments, expected_total_segments(&paths));
        prop_assert_eq!(metrics.shared_segments, expected_shared_segments(&paths));
    }

    #[test]
    fn compute_transcript_complexity_shared_segments_follow_presence_semantics(
        paths in prop::collection::vec(prop::collection::vec(0usize..24, 0..40), 1..24)
    ) {
        let transcripts: Vec<Transcript> = paths
            .iter()
            .enumerate()
            .map(|(id, path)| Transcript::new(id, String::new(), path.clone(), 1.0))
            .collect();

        let metrics = compute_transcript_complexity(&transcripts);
        prop_assert_eq!(metrics.total_paths, paths.len());
        prop_assert_eq!(metrics.total_segments, expected_total_segments(&paths));
        prop_assert_eq!(metrics.shared_segments, expected_shared_segments(&paths));
    }

    #[test]
    fn count_bubbles_matches_reference_and_is_insertion_order_invariant(
        node_count in 1usize..=7,
        edge_items in prop::collection::vec((0usize..7, 0usize..7), 0..32),
        seed in any::<u64>(),
    ) {
        let mut edges = Vec::new();
        let mut seen = HashSet::new();
        for (from, to) in edge_items {
            if from < node_count && to < node_count && from != to && seen.insert((from, to)) {
                edges.push((from, to));
            }
        }

        let graph_a = build_graph(node_count, &edges);
        let expected = reference_bubble_count(&graph_a);
        let observed_a = count_bubbles(&graph_a);
        prop_assert_eq!(observed_a, expected);

        let mut edges_b = edges.clone();
        if !edges_b.is_empty() {
            let rotate_by = (seed as usize) % edges_b.len();
            edges_b.rotate_left(rotate_by);
            edges_b.reverse();
        }
        let graph_b = build_graph(node_count, &edges_b);
        let observed_b = count_bubbles(&graph_b);
        prop_assert_eq!(observed_b, expected);
    }
}
