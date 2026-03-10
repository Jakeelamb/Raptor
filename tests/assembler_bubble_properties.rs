use ahash::AHashMap;
use proptest::prelude::*;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use raptor::accel::backend::AdjacencyTableU64;
use raptor::graph::assembler::detect_bubbles;

fn add_linear_edges(edges: &mut Vec<(u64, u64, u32)>, nodes: &[u64], end: u64, cov: u32) {
    for pair in nodes.windows(2) {
        edges.push((pair[0], pair[1], cov));
    }
    if let Some(&last) = nodes.last() {
        edges.push((last, end, cov));
    }
}

fn make_counts(nodes: impl IntoIterator<Item = u64>, cov: u32) -> AHashMap<u64, u32> {
    let mut counts = AHashMap::new();
    for node in nodes {
        counts.insert(node, cov);
    }
    counts
}

fn build_graph(k: u8, edges: &[(u64, u64, u32)]) -> AdjacencyTableU64 {
    let mut adjacency = AdjacencyTableU64::new(k);
    for &(from, to, cov) in edges {
        adjacency.add_edge(from, to, cov);
    }
    adjacency
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(96))]

    #[test]
    fn detect_bubbles_is_invariant_to_edge_insertion_order(
        branch_a_len in 1usize..=6,
        branch_b_len in 1usize..=6,
        seed in any::<u64>(),
    ) {
        let start = 1u64;
        let branch_a_start = 2u64;
        let branch_b_start = 3u64;
        let end = 4u64;

        let mut next = 10u64;
        let mut branch_a_nodes = vec![branch_a_start];
        for _ in 1..branch_a_len {
            branch_a_nodes.push(next);
            next += 1;
        }

        let mut branch_b_nodes = vec![branch_b_start];
        for _ in 1..branch_b_len {
            branch_b_nodes.push(next);
            next += 1;
        }

        let mut edges = vec![(start, branch_a_start, 10), (start, branch_b_start, 10)];
        add_linear_edges(&mut edges, &branch_a_nodes, end, 9);
        add_linear_edges(&mut edges, &branch_b_nodes, end, 9);

        let mut all_nodes = vec![start, end];
        all_nodes.extend(branch_a_nodes.iter().copied());
        all_nodes.extend(branch_b_nodes.iter().copied());
        all_nodes.sort_unstable();
        all_nodes.dedup();
        let counts = make_counts(all_nodes, 20);

        let graph_a = build_graph(3, &edges);

        let mut shuffled = edges.clone();
        let mut rng = StdRng::seed_from_u64(seed);
        shuffled.shuffle(&mut rng);
        let graph_b = build_graph(3, &shuffled);

        let bubbles_a = detect_bubbles(&graph_a, &counts, 3, branch_a_len.max(branch_b_len) + 4);
        let bubbles_b = detect_bubbles(&graph_b, &counts, 3, branch_a_len.max(branch_b_len) + 4);

        prop_assert_eq!(bubbles_a.len(), 1);
        prop_assert_eq!(bubbles_b.len(), 1);
        prop_assert_eq!(bubbles_a[0].start, start);
        prop_assert_eq!(bubbles_a[0].end, end);
        prop_assert_eq!(&bubbles_a[0].path1, &bubbles_b[0].path1);
        prop_assert_eq!(&bubbles_a[0].path2, &bubbles_b[0].path2);
    }

    #[test]
    fn detect_bubbles_ignores_start_rejoin_cycles(
        branch_a_len in 1usize..=6,
        branch_b_len in 1usize..=6,
        seed in any::<u64>(),
    ) {
        let start = 1u64;
        let branch_a_start = 2u64;
        let branch_b_start = 3u64;

        let mut next = 10u64;
        let mut branch_a_nodes = vec![branch_a_start];
        for _ in 1..branch_a_len {
            branch_a_nodes.push(next);
            next += 1;
        }

        let mut branch_b_nodes = vec![branch_b_start];
        for _ in 1..branch_b_len {
            branch_b_nodes.push(next);
            next += 1;
        }

        let mut edges = vec![(start, branch_a_start, 10), (start, branch_b_start, 10)];

        for pair in branch_a_nodes.windows(2) {
            edges.push((pair[0], pair[1], 9));
        }
        edges.push((branch_a_nodes[branch_a_nodes.len() - 1], start, 9));

        for pair in branch_b_nodes.windows(2) {
            edges.push((pair[0], pair[1], 9));
        }
        edges.push((branch_b_nodes[branch_b_nodes.len() - 1], start, 9));

        let mut all_nodes = vec![start];
        all_nodes.extend(branch_a_nodes.iter().copied());
        all_nodes.extend(branch_b_nodes.iter().copied());
        all_nodes.sort_unstable();
        all_nodes.dedup();
        let counts = make_counts(all_nodes, 20);

        let mut shuffled = edges.clone();
        let mut rng = StdRng::seed_from_u64(seed);
        shuffled.shuffle(&mut rng);
        let graph = build_graph(3, &shuffled);

        let bubbles = detect_bubbles(&graph, &counts, 3, branch_a_len.max(branch_b_len) + 4);
        prop_assert!(bubbles.is_empty());
    }
}
