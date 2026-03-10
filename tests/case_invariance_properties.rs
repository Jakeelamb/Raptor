use proptest::prelude::*;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use raptor::graph::transcript::{calculate_transcript_stats, Transcript, TRANSCRIPT_STATS_KEYS};
use raptor::stats::calculate_stats;
use std::io::Write;
use tempfile::NamedTempFile;

fn randomized_case(sequence: &str, seed: u64) -> String {
    let mut rng = StdRng::seed_from_u64(seed);
    sequence
        .bytes()
        .map(|base| {
            if base.is_ascii_alphabetic() {
                if rng.gen_bool(0.5) {
                    base.to_ascii_uppercase() as char
                } else {
                    base.to_ascii_lowercase() as char
                }
            } else {
                base as char
            }
        })
        .collect()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    #[test]
    fn fasta_stats_are_case_invariant(
        contigs in prop::collection::vec("[ACGTNRYacgtnry]{0,80}", 1..12),
        seed in any::<u64>(),
    ) {
        let mut baseline = NamedTempFile::new().unwrap();
        let mut perturbed = NamedTempFile::new().unwrap();
        let mut case_rng = StdRng::seed_from_u64(seed);

        for (idx, seq) in contigs.iter().enumerate() {
            writeln!(baseline, ">contig_{idx}").unwrap();
            writeln!(baseline, "{seq}").unwrap();

            writeln!(perturbed, ">contig_{idx}").unwrap();
            let perturbed_seq = randomized_case(seq, case_rng.gen());
            writeln!(perturbed, "{perturbed_seq}").unwrap();
        }

        let a = calculate_stats(baseline.path().to_str().unwrap()).unwrap();
        let b = calculate_stats(perturbed.path().to_str().unwrap()).unwrap();

        prop_assert_eq!(a.total_contigs, b.total_contigs);
        prop_assert_eq!(a.total_length, b.total_length);
        prop_assert_eq!(a.ungapped_total_length, b.ungapped_total_length);
        prop_assert_eq!(a.gc_bases, b.gc_bases);
        prop_assert_eq!(a.acgt_bases, b.acgt_bases);
        prop_assert_eq!(a.n_bases, b.n_bases);
        prop_assert_eq!(a.ambiguous_bases, b.ambiguous_bases);
        prop_assert_eq!(a.total_rle_runs, b.total_rle_runs);
        prop_assert_eq!(a.n_run_count, b.n_run_count);
        prop_assert_eq!(a.max_n_run, b.max_n_run);
        prop_assert!((a.mean_rle_ratio - b.mean_rle_ratio).abs() < 1e-12);
        prop_assert!((a.length_weighted_rle_ratio - b.length_weighted_rle_ratio).abs() < 1e-12);
        prop_assert!((a.mean_n_run_length - b.mean_n_run_length).abs() < 1e-12);
        prop_assert!((a.gc_content - b.gc_content).abs() < 1e-12);
        prop_assert!((a.n_content - b.n_content).abs() < 1e-12);
        prop_assert!((a.ambiguous_content - b.ambiguous_content).abs() < 1e-12);
    }

    #[test]
    fn transcript_stats_are_case_invariant(
        sequences in prop::collection::vec("[ACGTNRYacgtnry]{1,120}", 1..24),
        seed in any::<u64>(),
    ) {
        let baseline: Vec<Transcript> = sequences
            .iter()
            .enumerate()
            .map(|(id, seq)| Transcript::new(id, seq.clone(), vec![id], 0.5))
            .collect();

        let mut case_rng = StdRng::seed_from_u64(seed);
        let perturbed: Vec<Transcript> = sequences
            .iter()
            .enumerate()
            .map(|(id, seq)| {
                let seq = randomized_case(seq, case_rng.gen());
                Transcript::new(id, seq, vec![id], 0.5)
            })
            .collect();

        let a = calculate_transcript_stats(&baseline);
        let b = calculate_transcript_stats(&perturbed);

        for &key in TRANSCRIPT_STATS_KEYS {
            let lhs = a.get(key).copied().unwrap_or(f64::NAN);
            let rhs = b.get(key).copied().unwrap_or(f64::NAN);
            prop_assert!((lhs - rhs).abs() < 1e-12, "metric mismatch for key `{key}`: {lhs} vs {rhs}");
        }
    }
}
