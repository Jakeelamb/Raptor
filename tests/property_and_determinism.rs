use ahash::AHashMap;
use petgraph::graphmap::DiGraphMap;
use proptest::prelude::*;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use raptor::accel::backend::AdjacencyTableU64;
use raptor::graph::assembler::{cleanup_graph, greedy_assembly_u64};
use raptor::graph::isoform_traverse::find_directed_paths;
use raptor::kmer::kmer::{encode_kmer, reverse_complement, KmerU64};

fn dna_string(max_len: usize) -> impl Strategy<Value = String> {
    prop::collection::vec(
        prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')],
        1..=max_len,
    )
    .prop_map(|chars| chars.into_iter().collect())
}

proptest! {
    #[test]
    fn kmer_u64_round_trip(seq in dna_string(32)) {
        let kmer = KmerU64::from_str(&seq).expect("generator emits only valid DNA");
        prop_assert_eq!(kmer.decode(), seq);
    }

    #[test]
    fn kmer_canonical_and_rc_invariants(seq in dna_string(32)) {
        let kmer = KmerU64::from_str(&seq).expect("generator emits only valid DNA");
        let rc = kmer.reverse_complement();

        prop_assert_eq!(rc.reverse_complement().encoded, kmer.encoded);
        prop_assert_eq!(rc.reverse_complement().len, kmer.len);

        let rc_str = reverse_complement(&seq);
        let expected = if seq <= rc_str {
            seq.clone()
        } else {
            rc_str
        };
        prop_assert_eq!(kmer.canonical().decode(), expected);
    }

    #[test]
    fn kmer_extend_matches_sliding_window(seq in dna_string(32), base in prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')]) {
        let kmer = KmerU64::from_str(&seq).expect("generator emits only valid DNA");
        let extended = kmer.extend(base as u8).expect("base is valid DNA");

        let expected = format!("{}{}", &seq[1..], base);
        prop_assert_eq!(extended.decode(), expected);
    }
}

#[test]
fn greedy_assembly_is_stable_under_randomized_insertion_order() {
    let k = 3;
    let kmers = vec![("AAA", 10), ("CAA", 10), ("GAA", 6), ("AAC", 5), ("AAG", 5)];
    let edges = vec![
        ("CAA", "AAA", 6),
        ("GAA", "AAA", 6),
        ("AAA", "AAG", 5),
        ("AAA", "AAC", 5),
    ];
    let expected = vec!["CAAAC".to_string(), "GAA".to_string(), "AAG".to_string()];

    let mut rng = StdRng::seed_from_u64(0x5EED_u64);

    for _ in 0..128 {
        let mut shuffled_kmers = kmers.clone();
        shuffled_kmers.shuffle(&mut rng);
        let mut counts = AHashMap::new();
        for (seq, count) in &shuffled_kmers {
            counts.insert(encode_kmer(seq).unwrap(), *count);
        }

        let mut shuffled_edges = edges.clone();
        shuffled_edges.shuffle(&mut rng);
        let mut adjacency = AdjacencyTableU64::new(k as u8);
        for (from, to, count) in &shuffled_edges {
            adjacency.add_edge(encode_kmer(from).unwrap(), encode_kmer(to).unwrap(), *count);
        }

        let contigs = greedy_assembly_u64(k, &counts, &adjacency, k);
        let observed: Vec<String> = contigs.into_iter().map(|c| c.sequence).collect();
        assert_eq!(observed, expected);
    }
}

#[test]
fn isoform_path_enumeration_is_stable_under_randomized_edge_order() {
    let nodes: Vec<usize> = (0..6).collect();
    let edges = vec![
        (0, 1, 0.9),
        (0, 2, 0.7),
        (1, 3, 0.8),
        (2, 3, 0.6),
        (3, 4, 0.95),
        (2, 5, 0.75),
    ];
    let start_nodes = vec![2, 0];
    let end_nodes = vec![4, 5];

    let mut baseline_graph = DiGraphMap::new();
    for &node in &nodes {
        baseline_graph.add_node(node);
    }
    for &(u, v, w) in &edges {
        baseline_graph.add_edge(u, v, w);
    }
    let baseline_paths: Vec<Vec<usize>> =
        find_directed_paths(&baseline_graph, &start_nodes, &end_nodes, 6)
            .into_iter()
            .map(|p| p.nodes)
            .collect();

    let mut rng = StdRng::seed_from_u64(0x150F0F1_u64);
    for _ in 0..128 {
        let mut shuffled_edges = edges.clone();
        shuffled_edges.shuffle(&mut rng);

        let mut graph = DiGraphMap::new();
        for &node in &nodes {
            graph.add_node(node);
        }
        for &(u, v, w) in &shuffled_edges {
            graph.add_edge(u, v, w);
        }

        let observed: Vec<Vec<usize>> = find_directed_paths(&graph, &start_nodes, &end_nodes, 6)
            .into_iter()
            .map(|p| p.nodes)
            .collect();
        assert_eq!(observed, baseline_paths);
    }
}

fn canonicalize_adjacency(adjacency: &AdjacencyTableU64) -> Vec<(u64, Vec<(u64, u32)>)> {
    let mut entries: Vec<(u64, Vec<(u64, u32)>)> = adjacency
        .forward
        .iter()
        .filter_map(|(&node, neighbors)| {
            let mut sorted_neighbors = neighbors.clone();
            sorted_neighbors.sort_unstable();
            if sorted_neighbors.is_empty() {
                None
            } else {
                Some((node, sorted_neighbors))
            }
        })
        .collect();
    entries.sort_unstable_by_key(|(node, _)| *node);
    entries
}

#[test]
fn cleanup_graph_is_stable_under_randomized_kmer_and_edge_insertion_order() {
    let k = 3;
    let kmers = vec![
        ("TAA", 1),
        ("AAA", 1),
        ("AAT", 1),
        ("ATC", 10),
        ("GTC", 10),
        ("TCA", 6),
        ("TCG", 6),
        ("CGA", 8),
    ];
    let edges = vec![
        ("TAA", "AAA", 1),
        ("AAA", "AAT", 1),
        ("AAT", "ATC", 1),
        ("GTC", "ATC", 10),
        ("ATC", "TCA", 6),
        ("ATC", "TCG", 6),
        ("TCA", "CGA", 6),
        ("TCG", "CGA", 6),
    ];

    let mut baseline_counts = AHashMap::new();
    for (seq, count) in &kmers {
        baseline_counts.insert(encode_kmer(seq).unwrap(), *count);
    }
    let mut baseline_adjacency = AdjacencyTableU64::new(k as u8);
    for (from, to, count) in &edges {
        baseline_adjacency.add_edge(encode_kmer(from).unwrap(), encode_kmer(to).unwrap(), *count);
    }
    let baseline_summary = cleanup_graph(&mut baseline_adjacency, &baseline_counts, k, 3);
    let baseline_graph = canonicalize_adjacency(&baseline_adjacency);

    let mut rng = StdRng::seed_from_u64(0xC1EA_D00D_u64);
    for _ in 0..128 {
        let mut shuffled_kmers = kmers.clone();
        shuffled_kmers.shuffle(&mut rng);
        let mut counts = AHashMap::new();
        for (seq, count) in &shuffled_kmers {
            counts.insert(encode_kmer(seq).unwrap(), *count);
        }

        let mut shuffled_edges = edges.clone();
        shuffled_edges.shuffle(&mut rng);
        let mut adjacency = AdjacencyTableU64::new(k as u8);
        for (from, to, count) in &shuffled_edges {
            adjacency.add_edge(encode_kmer(from).unwrap(), encode_kmer(to).unwrap(), *count);
        }

        let observed_summary = cleanup_graph(&mut adjacency, &counts, k, 3);
        let observed_graph = canonicalize_adjacency(&adjacency);
        assert_eq!(observed_summary, baseline_summary);
        assert_eq!(observed_graph, baseline_graph);
    }
}

#[test]
fn cleanup_graph_collapses_bubble_without_tip_removal_under_randomized_edge_order() {
    let k = 3;
    let kmers = vec![("AAA", 20), ("AAC", 10), ("AAG", 10), ("ACC", 20)];
    let edges = vec![
        ("AAA", "AAC", 10),
        ("AAA", "AAG", 10),
        ("AAC", "ACC", 10),
        ("AAG", "ACC", 10),
    ];

    let mut expected_counts = AHashMap::new();
    for (seq, count) in &kmers {
        expected_counts.insert(encode_kmer(seq).unwrap(), *count);
    }
    let aaa = encode_kmer("AAA").unwrap();
    let aac = encode_kmer("AAC").unwrap();
    let acc = encode_kmer("ACC").unwrap();

    let mut rng = StdRng::seed_from_u64(0x00B0_BB1E_u64);
    for _ in 0..128 {
        let mut shuffled_edges = edges.clone();
        shuffled_edges.shuffle(&mut rng);
        let mut adjacency = AdjacencyTableU64::new(k as u8);
        for (from, to, count) in &shuffled_edges {
            adjacency.add_edge(encode_kmer(from).unwrap(), encode_kmer(to).unwrap(), *count);
        }

        let summary = cleanup_graph(&mut adjacency, &expected_counts, k, 3);
        assert_eq!(summary, (0, 1));

        let mut observed_graph = canonicalize_adjacency(&adjacency);
        observed_graph.sort_unstable_by_key(|(node, _)| *node);
        assert_eq!(
            observed_graph,
            vec![(aaa, vec![(aac, 10)]), (aac, vec![(acc, 10)])]
        );
    }
}
