use proptest::prelude::*;
use raptor::graph::stats::{
    compute_branchiness, compute_path_complexity, compute_transcript_complexity,
};
use raptor::graph::stitch::Path;
use raptor::graph::transcript::Transcript;
use std::collections::{HashMap, HashSet};

fn expected_shared_segments(paths: &[Vec<usize>]) -> usize {
    let mut usage = HashMap::new();
    let mut seen_in_path = HashSet::new();

    for path in paths {
        seen_in_path.clear();
        for &segment in path {
            if seen_in_path.insert(segment) {
                *usage.entry(segment).or_insert(0usize) += 1;
            }
        }
    }

    usage.values().filter(|&&count| count > 1).count()
}

fn expected_total_segments(paths: &[Vec<usize>]) -> usize {
    let mut all_segments = HashSet::new();
    for path in paths {
        for &segment in path {
            all_segments.insert(segment);
        }
    }
    all_segments.len()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    #[test]
    fn compute_branchiness_matches_per_path_presence_even_with_repetitions(
        paths in prop::collection::vec(prop::collection::vec(0usize..24, 0..40), 1..24)
    ) {
        let as_strings: Vec<Vec<String>> = paths
            .iter()
            .map(|path| path.iter().map(|segment| format!("seg_{segment}")).collect())
            .collect();

        let observed = compute_branchiness(&as_strings);
        let expected = expected_shared_segments(&paths);
        prop_assert_eq!(observed, expected);
    }

    #[test]
    fn compute_path_complexity_shared_segments_follow_presence_semantics(
        paths in prop::collection::vec(prop::collection::vec(0usize..24, 0..40), 1..24)
    ) {
        let assembled_paths: Vec<Path> = paths
            .iter()
            .enumerate()
            .map(|(id, segments)| Path {
                id,
                segments: segments.clone(),
                overlaps: vec![0; segments.len().saturating_sub(1)],
            })
            .collect();

        let metrics = compute_path_complexity(&assembled_paths);
        prop_assert_eq!(metrics.total_paths, paths.len());
        prop_assert_eq!(metrics.total_segments, expected_total_segments(&paths));
        prop_assert_eq!(metrics.shared_segments, expected_shared_segments(&paths));
    }

    #[test]
    fn compute_transcript_complexity_shared_segments_follow_presence_semantics(
        paths in prop::collection::vec(prop::collection::vec(0usize..24, 0..40), 1..24)
    ) {
        let transcripts: Vec<Transcript> = paths
            .iter()
            .enumerate()
            .map(|(id, path)| Transcript::new(id, String::new(), path.clone(), 1.0))
            .collect();

        let metrics = compute_transcript_complexity(&transcripts);
        prop_assert_eq!(metrics.total_paths, paths.len());
        prop_assert_eq!(metrics.total_segments, expected_total_segments(&paths));
        prop_assert_eq!(metrics.shared_segments, expected_shared_segments(&paths));
    }
}
