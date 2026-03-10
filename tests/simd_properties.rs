use proptest::prelude::*;
use raptor::accel::simd::match_kmers_with_overlap;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    #[test]
    fn match_kmers_with_overlap_returns_none_when_min_overlap_exceeds_shorter_input(
        query in "[ACGTN]{0,64}",
        target in "[ACGTN]{0,64}",
        extra in 1usize..8
    ) {
        let min_overlap = query.len().min(target.len()) + extra;
        prop_assert_eq!(
            match_kmers_with_overlap(&query, &target, min_overlap, 0),
            None
        );
    }
}
