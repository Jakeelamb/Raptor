use ahash::AHashMap;
use proptest::prelude::*;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use raptor::accel::backend::AdjacencyTableU64;
use raptor::graph::assembler::greedy_assembly_u64;
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
