use proptest::prelude::*;
use raptor::graph::transcript::Transcript;
use raptor::quant::tpm::compute_tpm;

fn make_transcripts(lengths: &[usize]) -> Vec<Transcript> {
    lengths
        .iter()
        .enumerate()
        .map(|(idx, &len)| Transcript {
            id: idx + 1,
            sequence: "A".repeat(len),
            path: vec![idx + 1],
            confidence: 1.0,
            length: len,
            strand: '+',
            tpm: None,
            splicing: "linear".to_string(),
        })
        .collect()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn compute_tpm_normalization_and_scaling_invariants_hold(
        (lengths, counts, scale) in (1usize..32).prop_flat_map(|n| {
            (
                prop::collection::vec(1usize..20_000, n),
                prop::collection::vec(0usize..50_000, n),
                1usize..128,
            )
        })
    ) {
        let transcripts = make_transcripts(&lengths);
        let base = compute_tpm(&counts, &transcripts);

        prop_assert_eq!(base.len(), counts.len());
        prop_assert!(base.iter().all(|v| v.is_finite() && *v >= 0.0));

        if counts.iter().any(|&c| c > 0) {
            let sum = base.iter().sum::<f64>();
            prop_assert!((sum - 1_000_000.0).abs() < 1e-6);

            let scaled_counts: Vec<usize> = counts
                .iter()
                .map(|&count| count.saturating_mul(scale))
                .collect();
            let scaled = compute_tpm(&scaled_counts, &transcripts);
            prop_assert_eq!(scaled.len(), base.len());

            for (lhs, rhs) in base.iter().zip(scaled.iter()) {
                prop_assert!((lhs - rhs).abs() < 1e-6);
            }
        } else {
            prop_assert!(base.iter().all(|&v| v == 0.0));
        }
    }
}
