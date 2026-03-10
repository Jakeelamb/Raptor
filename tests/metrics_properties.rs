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
        prop_assert!(stats.n25 >= stats.n50);
        prop_assert!(stats.n50 >= stats.n75);
        prop_assert!(stats.n75 >= stats.n90);
        prop_assert!(stats.n90 >= stats.n95);
        prop_assert!(stats.n95 >= stats.n99);
        prop_assert!(stats.l25 <= stats.l50);
        prop_assert!(stats.l50 <= stats.l75);
        prop_assert!(stats.l75 <= stats.l90);
        prop_assert!(stats.l50 <= stats.l90);
        prop_assert!(stats.l90 <= stats.l95);
        prop_assert!(stats.l95 <= stats.l99);
        prop_assert!(stats.longest >= stats.n50);
        prop_assert!(stats.median_length <= stats.longest as f64 + 1e-9);
        prop_assert!(stats.au_n + 1e-9 >= stats.avg_length);
        prop_assert!(stats.au_n <= stats.longest as f64 + 1e-9);
        prop_assert!(stats.contigs_ge_50kb <= stats.contigs_ge_10kb);
        prop_assert!(stats.contigs_ge_100kb <= stats.contigs_ge_50kb);
        prop_assert!(stats.contigs_ge_10kb <= stats.contigs_ge_1kb);
        prop_assert!(stats.contigs_ge_1kb <= stats.total);
        prop_assert!(stats.contigs_ge_100kb_frac <= stats.contigs_ge_50kb_frac + 1e-12);
        prop_assert!(stats.contigs_ge_50kb_frac <= stats.contigs_ge_10kb_frac + 1e-12);
        prop_assert!(stats.contigs_ge_10kb_frac <= stats.contigs_ge_1kb_frac + 1e-12);
        prop_assert!(stats.contigs_ge_1kb_frac <= 1.0 + 1e-12);
        prop_assert!(stats.bases_ge_50kb <= stats.bases_ge_10kb);
        prop_assert!(stats.bases_ge_100kb <= stats.bases_ge_50kb);
        prop_assert!(stats.bases_ge_10kb <= stats.bases_ge_1kb);
        prop_assert!(stats.bases_ge_1kb <= stats.total_bases);
        prop_assert!(stats.bases_ge_100kb_frac <= stats.bases_ge_50kb_frac + 1e-12);
        prop_assert!(stats.bases_ge_50kb_frac <= stats.bases_ge_10kb_frac + 1e-12);
        prop_assert!(stats.bases_ge_10kb_frac <= stats.bases_ge_1kb_frac + 1e-12);
        prop_assert!(stats.bases_ge_1kb_frac <= 1.0 + 1e-12);
        prop_assert!((stats.contigs_ge_1kb_frac - (stats.contigs_ge_1kb as f64 / stats.total as f64)).abs() < 1e-12);
        if stats.total_bases == 0 {
            prop_assert_eq!(stats.bases_ge_1kb_frac, 0.0);
            prop_assert_eq!(stats.bases_ge_10kb_frac, 0.0);
            prop_assert_eq!(stats.bases_ge_50kb_frac, 0.0);
            prop_assert_eq!(stats.bases_ge_100kb_frac, 0.0);
        } else {
            prop_assert!((stats.bases_ge_1kb_frac - (stats.bases_ge_1kb as f64 / stats.total_bases as f64)).abs() < 1e-12);
        }
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
        prop_assert_eq!(observed.n25, baseline.n25);
        prop_assert_eq!(observed.n50, baseline.n50);
        prop_assert_eq!(observed.n75, baseline.n75);
        prop_assert_eq!(observed.n90, baseline.n90);
        prop_assert_eq!(observed.n95, baseline.n95);
        prop_assert_eq!(observed.n99, baseline.n99);
        prop_assert_eq!(observed.l25, baseline.l25);
        prop_assert_eq!(observed.l50, baseline.l50);
        prop_assert_eq!(observed.l75, baseline.l75);
        prop_assert_eq!(observed.l90, baseline.l90);
        prop_assert_eq!(observed.l95, baseline.l95);
        prop_assert_eq!(observed.l99, baseline.l99);
        prop_assert_eq!(observed.longest, baseline.longest);
        prop_assert!((observed.median_length - baseline.median_length).abs() < 1e-12);
        prop_assert_eq!(observed.contigs_ge_1kb, baseline.contigs_ge_1kb);
        prop_assert_eq!(observed.contigs_ge_10kb, baseline.contigs_ge_10kb);
        prop_assert_eq!(observed.contigs_ge_50kb, baseline.contigs_ge_50kb);
        prop_assert_eq!(observed.contigs_ge_100kb, baseline.contigs_ge_100kb);
        prop_assert_eq!(observed.bases_ge_1kb, baseline.bases_ge_1kb);
        prop_assert_eq!(observed.bases_ge_10kb, baseline.bases_ge_10kb);
        prop_assert_eq!(observed.bases_ge_50kb, baseline.bases_ge_50kb);
        prop_assert_eq!(observed.bases_ge_100kb, baseline.bases_ge_100kb);
        prop_assert!((observed.contigs_ge_1kb_frac - baseline.contigs_ge_1kb_frac).abs() < 1e-12);
        prop_assert!((observed.contigs_ge_10kb_frac - baseline.contigs_ge_10kb_frac).abs() < 1e-12);
        prop_assert!((observed.contigs_ge_50kb_frac - baseline.contigs_ge_50kb_frac).abs() < 1e-12);
        prop_assert!((observed.contigs_ge_100kb_frac - baseline.contigs_ge_100kb_frac).abs() < 1e-12);
        prop_assert!((observed.bases_ge_1kb_frac - baseline.bases_ge_1kb_frac).abs() < 1e-12);
        prop_assert!((observed.bases_ge_10kb_frac - baseline.bases_ge_10kb_frac).abs() < 1e-12);
        prop_assert!((observed.bases_ge_50kb_frac - baseline.bases_ge_50kb_frac).abs() < 1e-12);
        prop_assert!((observed.bases_ge_100kb_frac - baseline.bases_ge_100kb_frac).abs() < 1e-12);
        prop_assert!((observed.avg_length - baseline.avg_length).abs() < 1e-12);
        prop_assert!((observed.au_n - baseline.au_n).abs() < 1e-12);
    }

    #[test]
    fn evaluate_lengths_handles_zero_length_contigs_consistently(
        lengths in prop::collection::vec(0usize..2_000, 1..64),
        seed in any::<u64>()
    ) {
        let baseline = evaluate_lengths(&lengths);
        let total: usize = lengths.iter().sum();

        prop_assert_eq!(baseline.total, lengths.len());
        prop_assert_eq!(baseline.total_bases, total);
        if total == 0 {
            prop_assert_eq!(baseline.n25, 0);
            prop_assert_eq!(baseline.n50, 0);
            prop_assert_eq!(baseline.n75, 0);
            prop_assert_eq!(baseline.n90, 0);
            prop_assert_eq!(baseline.n95, 0);
            prop_assert_eq!(baseline.n99, 0);
            prop_assert_eq!(baseline.l25, 0);
            prop_assert_eq!(baseline.l50, 0);
            prop_assert_eq!(baseline.l75, 0);
            prop_assert_eq!(baseline.l90, 0);
            prop_assert_eq!(baseline.l95, 0);
            prop_assert_eq!(baseline.l99, 0);
            prop_assert_eq!(baseline.longest, 0);
            prop_assert_eq!(baseline.au_n, 0.0);
            prop_assert_eq!(baseline.median_length, 0.0);
        }

        let mut shuffled = lengths.clone();
        let mut rng = StdRng::seed_from_u64(seed);
        shuffled.shuffle(&mut rng);
        let observed = evaluate_lengths(&shuffled);

        prop_assert_eq!(observed.total, baseline.total);
        prop_assert_eq!(observed.total_bases, baseline.total_bases);
        prop_assert_eq!(observed.n25, baseline.n25);
        prop_assert_eq!(observed.n50, baseline.n50);
        prop_assert_eq!(observed.n75, baseline.n75);
        prop_assert_eq!(observed.n90, baseline.n90);
        prop_assert_eq!(observed.n95, baseline.n95);
        prop_assert_eq!(observed.n99, baseline.n99);
        prop_assert_eq!(observed.l25, baseline.l25);
        prop_assert_eq!(observed.l50, baseline.l50);
        prop_assert_eq!(observed.l75, baseline.l75);
        prop_assert_eq!(observed.l90, baseline.l90);
        prop_assert_eq!(observed.l95, baseline.l95);
        prop_assert_eq!(observed.l99, baseline.l99);
        prop_assert_eq!(observed.longest, baseline.longest);
        prop_assert!((observed.median_length - baseline.median_length).abs() < 1e-12);
        prop_assert_eq!(observed.contigs_ge_1kb, baseline.contigs_ge_1kb);
        prop_assert_eq!(observed.contigs_ge_10kb, baseline.contigs_ge_10kb);
        prop_assert_eq!(observed.contigs_ge_50kb, baseline.contigs_ge_50kb);
        prop_assert_eq!(observed.contigs_ge_100kb, baseline.contigs_ge_100kb);
        prop_assert_eq!(observed.bases_ge_1kb, baseline.bases_ge_1kb);
        prop_assert_eq!(observed.bases_ge_10kb, baseline.bases_ge_10kb);
        prop_assert_eq!(observed.bases_ge_50kb, baseline.bases_ge_50kb);
        prop_assert_eq!(observed.bases_ge_100kb, baseline.bases_ge_100kb);
        prop_assert!((observed.contigs_ge_1kb_frac - baseline.contigs_ge_1kb_frac).abs() < 1e-12);
        prop_assert!((observed.contigs_ge_10kb_frac - baseline.contigs_ge_10kb_frac).abs() < 1e-12);
        prop_assert!((observed.contigs_ge_50kb_frac - baseline.contigs_ge_50kb_frac).abs() < 1e-12);
        prop_assert!((observed.contigs_ge_100kb_frac - baseline.contigs_ge_100kb_frac).abs() < 1e-12);
        prop_assert!((observed.bases_ge_1kb_frac - baseline.bases_ge_1kb_frac).abs() < 1e-12);
        prop_assert!((observed.bases_ge_10kb_frac - baseline.bases_ge_10kb_frac).abs() < 1e-12);
        prop_assert!((observed.bases_ge_50kb_frac - baseline.bases_ge_50kb_frac).abs() < 1e-12);
        prop_assert!((observed.bases_ge_100kb_frac - baseline.bases_ge_100kb_frac).abs() < 1e-12);
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
        prop_assert_eq!(observed.get("ungapped_total_length").copied().unwrap_or(-1.0) as usize, expected.total_bases);
        prop_assert_eq!(observed.get("mean_length").copied().unwrap_or(-1.0), expected.avg_length);
        prop_assert_eq!(observed.get("median_length").copied().unwrap_or(-1.0), expected.median_length);
        prop_assert_eq!(observed.get("n25").copied().unwrap_or(-1.0) as usize, expected.n25);
        prop_assert_eq!(observed.get("n50").copied().unwrap_or(-1.0) as usize, expected.n50);
        prop_assert_eq!(observed.get("n75").copied().unwrap_or(-1.0) as usize, expected.n75);
        prop_assert_eq!(observed.get("n90").copied().unwrap_or(-1.0) as usize, expected.n90);
        prop_assert_eq!(observed.get("n95").copied().unwrap_or(-1.0) as usize, expected.n95);
        prop_assert_eq!(observed.get("n99").copied().unwrap_or(-1.0) as usize, expected.n99);
        prop_assert_eq!(observed.get("l25").copied().unwrap_or(-1.0) as usize, expected.l25);
        prop_assert_eq!(observed.get("l50").copied().unwrap_or(-1.0) as usize, expected.l50);
        prop_assert_eq!(observed.get("l75").copied().unwrap_or(-1.0) as usize, expected.l75);
        prop_assert_eq!(observed.get("l90").copied().unwrap_or(-1.0) as usize, expected.l90);
        prop_assert_eq!(observed.get("l95").copied().unwrap_or(-1.0) as usize, expected.l95);
        prop_assert_eq!(observed.get("l99").copied().unwrap_or(-1.0) as usize, expected.l99);
        prop_assert!((observed.get("au_n").copied().unwrap_or(-1.0) - expected.au_n).abs() < 1e-12);
        prop_assert_eq!(observed.get("ungapped_n50").copied().unwrap_or(-1.0) as usize, expected.n50);
        prop_assert!((observed.get("ungapped_au_n").copied().unwrap_or(-1.0) - expected.au_n).abs() < 1e-12);
        prop_assert_eq!(observed.get("contigs_ge_1kb").copied().unwrap_or(-1.0) as usize, expected.contigs_ge_1kb);
        prop_assert_eq!(observed.get("contigs_ge_10kb").copied().unwrap_or(-1.0) as usize, expected.contigs_ge_10kb);
        prop_assert_eq!(observed.get("contigs_ge_50kb").copied().unwrap_or(-1.0) as usize, expected.contigs_ge_50kb);
        prop_assert_eq!(observed.get("contigs_ge_100kb").copied().unwrap_or(-1.0) as usize, expected.contigs_ge_100kb);
        prop_assert_eq!(observed.get("bases_ge_1kb").copied().unwrap_or(-1.0) as usize, expected.bases_ge_1kb);
        prop_assert_eq!(observed.get("bases_ge_10kb").copied().unwrap_or(-1.0) as usize, expected.bases_ge_10kb);
        prop_assert_eq!(observed.get("bases_ge_50kb").copied().unwrap_or(-1.0) as usize, expected.bases_ge_50kb);
        prop_assert_eq!(observed.get("bases_ge_100kb").copied().unwrap_or(-1.0) as usize, expected.bases_ge_100kb);
        prop_assert!((observed.get("contigs_ge_1kb_frac").copied().unwrap_or(-1.0) - expected.contigs_ge_1kb_frac).abs() < 1e-12);
        prop_assert!((observed.get("contigs_ge_10kb_frac").copied().unwrap_or(-1.0) - expected.contigs_ge_10kb_frac).abs() < 1e-12);
        prop_assert!((observed.get("contigs_ge_50kb_frac").copied().unwrap_or(-1.0) - expected.contigs_ge_50kb_frac).abs() < 1e-12);
        prop_assert!((observed.get("contigs_ge_100kb_frac").copied().unwrap_or(-1.0) - expected.contigs_ge_100kb_frac).abs() < 1e-12);
        prop_assert!((observed.get("bases_ge_1kb_frac").copied().unwrap_or(-1.0) - expected.bases_ge_1kb_frac).abs() < 1e-12);
        prop_assert!((observed.get("bases_ge_10kb_frac").copied().unwrap_or(-1.0) - expected.bases_ge_10kb_frac).abs() < 1e-12);
        prop_assert!((observed.get("bases_ge_50kb_frac").copied().unwrap_or(-1.0) - expected.bases_ge_50kb_frac).abs() < 1e-12);
        prop_assert!((observed.get("bases_ge_100kb_frac").copied().unwrap_or(-1.0) - expected.bases_ge_100kb_frac).abs() < 1e-12);
        prop_assert_eq!(observed.get("acgt_bases").copied().unwrap_or(-1.0) as usize, expected.total_bases);
        prop_assert_eq!(observed.get("n_bases").copied().unwrap_or(-1.0), 0.0);
        prop_assert_eq!(observed.get("ambiguous_bases").copied().unwrap_or(-1.0), 0.0);
        prop_assert_eq!(observed.get("gc_content").copied().unwrap_or(-1.0), 0.0);
        prop_assert_eq!(observed.get("gc_content_acgt").copied().unwrap_or(-1.0), 0.0);
        prop_assert_eq!(observed.get("n_content").copied().unwrap_or(-1.0), 0.0);
        prop_assert_eq!(observed.get("ambiguous_content").copied().unwrap_or(-1.0), 0.0);
        prop_assert_eq!(observed.get("length_field_mismatch_count").copied().unwrap_or(-1.0), 0.0);
        let expected_total_rle_runs: usize = lengths
            .iter()
            .map(|len| len.div_ceil(u8::MAX as usize))
            .sum();
        prop_assert_eq!(
            observed.get("total_rle_runs").copied().unwrap_or(-1.0) as usize,
            expected_total_rle_runs
        );
        prop_assert!(
            (observed
                .get("length_weighted_rle_ratio")
                .copied()
                .unwrap_or(-1.0)
                - (expected_total_rle_runs as f64 / expected.total_bases as f64))
                .abs()
                < 1e-12
        );

        let mut shuffled = lengths.clone();
        let mut rng = StdRng::seed_from_u64(seed);
        shuffled.shuffle(&mut rng);
        let shuffled_transcripts = transcripts_from_lengths(&shuffled);
        let shuffled_stats = calculate_transcript_stats(&shuffled_transcripts);

        for key in [
            "count",
            "total_length",
            "ungapped_total_length",
            "mean_length",
            "median_length",
            "n25",
            "n50",
            "n75",
            "n90",
            "n95",
            "n99",
            "l25",
            "l50",
            "l75",
            "l90",
            "l95",
            "l99",
            "au_n",
            "ungapped_n50",
            "ungapped_au_n",
            "contigs_ge_1kb",
            "contigs_ge_10kb",
            "contigs_ge_50kb",
            "contigs_ge_100kb",
            "bases_ge_1kb",
            "bases_ge_10kb",
            "bases_ge_50kb",
            "bases_ge_100kb",
            "contigs_ge_1kb_frac",
            "contigs_ge_10kb_frac",
            "contigs_ge_50kb_frac",
            "contigs_ge_100kb_frac",
            "bases_ge_1kb_frac",
            "bases_ge_10kb_frac",
            "bases_ge_50kb_frac",
            "bases_ge_100kb_frac",
            "gc_content",
            "gc_content_acgt",
            "n_content",
            "ambiguous_content",
            "acgt_bases",
            "n_bases",
            "ambiguous_bases",
            "length_field_mismatch_count",
            "total_rle_runs",
            "mean_rle_ratio",
            "length_weighted_rle_ratio",
        ] {
            let lhs = observed.get(key).copied().unwrap_or(f64::NAN);
            let rhs = shuffled_stats.get(key).copied().unwrap_or(f64::NAN);
            prop_assert!((lhs - rhs).abs() < 1e-12, "metric `{}` differed after permutation: {} vs {}", key, lhs, rhs);
        }
    }
}
