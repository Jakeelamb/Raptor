use proptest::prelude::*;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use raptor::eval::metrics::evaluate_lengths;
use raptor::graph::transcript::{calculate_transcript_stats, Transcript};

fn transcripts_from_lengths(lengths: &[usize]) -> Vec<Transcript> {
    lengths
        .iter()
        .enumerate()
        .map(|(id, &len)| Transcript {
            id,
            sequence: "A".repeat(len),
            path: vec![id],
            confidence: 1.0,
            length: len,
            strand: '+',
            tpm: None,
            splicing: "linear".to_string(),
        })
        .collect()
}

proptest! {
    #[test]
    fn evaluate_lengths_reports_consistent_invariants(
        lengths in prop::collection::vec(1usize..10_000, 1..64)
    ) {
        let stats = evaluate_lengths(&lengths);
        let total: usize = lengths.iter().sum();

        prop_assert_eq!(stats.total, lengths.len());
        prop_assert_eq!(stats.total_bases, total);
        prop_assert!(stats.n50 >= stats.n75);
        prop_assert!(stats.n75 >= stats.n90);
        prop_assert!(stats.n90 >= stats.n95);
        prop_assert!(stats.l50 <= stats.l90);
        prop_assert!(stats.l90 <= stats.l95);
        prop_assert!(stats.longest >= stats.n50);
        prop_assert!(stats.au_n + 1e-9 >= stats.avg_length);
        prop_assert!(stats.au_n <= stats.longest as f64 + 1e-9);
    }

    #[test]
    fn evaluate_lengths_is_permutation_invariant(
        lengths in prop::collection::vec(1usize..2_000, 1..64),
        seed in any::<u64>()
    ) {
        let baseline = evaluate_lengths(&lengths);

        let mut shuffled = lengths.clone();
        let mut rng = StdRng::seed_from_u64(seed);
        shuffled.shuffle(&mut rng);
        let observed = evaluate_lengths(&shuffled);

        prop_assert_eq!(observed.total, baseline.total);
        prop_assert_eq!(observed.total_bases, baseline.total_bases);
        prop_assert_eq!(observed.n50, baseline.n50);
        prop_assert_eq!(observed.n75, baseline.n75);
        prop_assert_eq!(observed.n90, baseline.n90);
        prop_assert_eq!(observed.n95, baseline.n95);
        prop_assert_eq!(observed.l50, baseline.l50);
        prop_assert_eq!(observed.l90, baseline.l90);
        prop_assert_eq!(observed.l95, baseline.l95);
        prop_assert_eq!(observed.longest, baseline.longest);
        prop_assert!((observed.avg_length - baseline.avg_length).abs() < 1e-12);
        prop_assert!((observed.au_n - baseline.au_n).abs() < 1e-12);
    }

    #[test]
    fn transcript_stats_length_metrics_match_evaluate_lengths_and_are_permutation_invariant(
        lengths in prop::collection::vec(1usize..2_000, 1..64),
        seed in any::<u64>()
    ) {
        let expected = evaluate_lengths(&lengths);
        let transcripts = transcripts_from_lengths(&lengths);
        let observed = calculate_transcript_stats(&transcripts);

        prop_assert_eq!(observed.get("count").copied().unwrap_or(-1.0) as usize, lengths.len());
        prop_assert_eq!(observed.get("total_length").copied().unwrap_or(-1.0) as usize, expected.total_bases);
        prop_assert_eq!(observed.get("mean_length").copied().unwrap_or(-1.0), expected.avg_length);
        prop_assert_eq!(observed.get("n50").copied().unwrap_or(-1.0) as usize, expected.n50);
        prop_assert_eq!(observed.get("n75").copied().unwrap_or(-1.0) as usize, expected.n75);
        prop_assert_eq!(observed.get("n90").copied().unwrap_or(-1.0) as usize, expected.n90);
        prop_assert_eq!(observed.get("n95").copied().unwrap_or(-1.0) as usize, expected.n95);
        prop_assert_eq!(observed.get("l50").copied().unwrap_or(-1.0) as usize, expected.l50);
        prop_assert_eq!(observed.get("l90").copied().unwrap_or(-1.0) as usize, expected.l90);
        prop_assert_eq!(observed.get("l95").copied().unwrap_or(-1.0) as usize, expected.l95);
        prop_assert!((observed.get("au_n").copied().unwrap_or(-1.0) - expected.au_n).abs() < 1e-12);

        let mut shuffled = lengths.clone();
        let mut rng = StdRng::seed_from_u64(seed);
        shuffled.shuffle(&mut rng);
        let shuffled_transcripts = transcripts_from_lengths(&shuffled);
        let shuffled_stats = calculate_transcript_stats(&shuffled_transcripts);

        for key in ["count", "total_length", "mean_length", "n50", "n75", "n90", "n95", "l50", "l90", "l95", "au_n"] {
            let lhs = observed.get(key).copied().unwrap_or(f64::NAN);
            let rhs = shuffled_stats.get(key).copied().unwrap_or(f64::NAN);
            prop_assert!((lhs - rhs).abs() < 1e-12, "metric `{}` differed after permutation: {} vs {}", key, lhs, rhs);
        }
    }
}
