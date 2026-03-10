use proptest::prelude::*;
use raptor::pipeline::large_genome_assembler::InsertSizeStats;

fn expected_midpoint_median(values: &[usize]) -> usize {
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 1 {
        sorted[mid]
    } else {
        (sorted[mid - 1] + sorted[mid]) / 2
    }
}

proptest! {
    #[test]
    fn insert_size_stats_reports_exact_median_midpoint(samples in prop::collection::vec(0usize..100_000, 1..256)) {
        let stats = InsertSizeStats::from_samples(&samples);
        let expected = expected_midpoint_median(&samples);
        prop_assert_eq!(stats.median, expected);
    }

    #[test]
    fn insert_size_stats_is_permutation_invariant(samples in prop::collection::vec(0usize..100_000, 1..256)) {
        let baseline = InsertSizeStats::from_samples(&samples);

        let mut permuted = samples.clone();
        if permuted.len() > 1 {
            let shift = permuted.len() / 3;
            permuted.rotate_left(shift.max(1));
        }
        let observed = InsertSizeStats::from_samples(&permuted);

        prop_assert_eq!(observed.min, baseline.min);
        prop_assert_eq!(observed.max, baseline.max);
        prop_assert_eq!(observed.median, baseline.median);
        prop_assert!((observed.mean - baseline.mean).abs() < 1e-9);
        prop_assert!((observed.std_dev - baseline.std_dev).abs() < 1e-9);
    }
}
