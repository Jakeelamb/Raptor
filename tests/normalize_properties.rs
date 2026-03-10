use proptest::prelude::*;
use raptor::kmer::normalize::estimate_read_abundance;
use raptor::kmer::nthash::NtHashIterator;

fn reference_abundance(seq: &str, k: usize, sketch: &[u32]) -> (u32, u32) {
    if k == 0 || seq.len() < k || sketch.is_empty() {
        return (0, 0);
    }

    let sketch_len = sketch.len() as u64;
    let mut values: Vec<u32> = NtHashIterator::new(seq.as_bytes(), k)
        .map(|(_, hash)| sketch[(hash % sketch_len) as usize])
        .collect();

    if values.is_empty() {
        return (0, 0);
    }

    let min = values.iter().copied().min().unwrap_or(0);
    let mid = values.len() / 2;
    values.select_nth_unstable(mid);
    (values[mid], min)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]
    #[test]
    fn estimate_read_abundance_matches_reference_across_sketch_layouts(
        seq in "[ACGTN]{0,96}",
        k in 0usize..20,
        sketch in prop::collection::vec(0u32..10_000, 0..128)
    ) {
        let observed = estimate_read_abundance(&seq, k, &sketch);
        let (expected_median, expected_min) = reference_abundance(&seq, k, &sketch);

        prop_assert_eq!(observed.median, expected_median);
        prop_assert_eq!(observed.min, expected_min);

        for _ in 0..4 {
            let repeat = estimate_read_abundance(&seq, k, &sketch);
            prop_assert_eq!(repeat.median, observed.median);
            prop_assert_eq!(repeat.min, observed.min);
        }
    }
}
