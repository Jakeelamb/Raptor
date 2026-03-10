use proptest::prelude::*;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use raptor::kmer::variable_k::{optimal_k, select_best_k};
use std::collections::HashMap;

fn dna_string_strategy(min_len: usize, max_len: usize) -> impl Strategy<Value = String> {
    prop::collection::vec(
        prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')],
        min_len..=max_len,
    )
    .prop_map(|chars| chars.into_iter().collect::<String>())
}

proptest! {
    #[test]
    fn select_best_k_is_invariant_to_histogram_insertion_order(
        entries in prop::collection::btree_map(1usize..=32, 0usize..20_000, 1..16),
        seed in any::<u64>()
    ) {
        let mut baseline = HashMap::new();
        let mut insertion_order = Vec::with_capacity(entries.len());
        for (k, count) in entries {
            baseline.insert(k, count);
            insertion_order.push((k, count));
        }
        let expected = select_best_k(&baseline);

        let mut shuffled = insertion_order;
        let mut rng = StdRng::seed_from_u64(seed);
        shuffled.shuffle(&mut rng);

        let mut observed_map = HashMap::new();
        for (k, count) in shuffled {
            observed_map.insert(k, count);
        }

        let observed = select_best_k(&observed_map);
        prop_assert_eq!(observed, expected);
    }

    #[test]
    fn optimal_k_is_always_within_nthash_bounds(
        sequences in prop::collection::vec(dna_string_strategy(1, 96), 1..32),
        max_k in 1usize..=64
    ) {
        let k = optimal_k(&sequences, max_k);
        let upper = max_k.min(32);

        prop_assert!(k >= 1);
        prop_assert!(k <= upper);
    }
}
