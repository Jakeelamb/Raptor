use proptest::prelude::*;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use raptor::graph::navigation::{extract_path_metadata, parse_gfa_paths, parse_gfa_paths_reader};
use std::collections::HashMap;
use std::io::Cursor;

fn build_path_lines(paths: &[Vec<(u8, bool)>]) -> Vec<String> {
    paths
        .iter()
        .enumerate()
        .map(|(idx, path)| {
            let segments = path
                .iter()
                .map(|(seg, forward)| format!("{}{}", seg, if *forward { '+' } else { '-' }))
                .collect::<Vec<_>>()
                .join(",");
            format!("P\tpath_{}\t{}\t*", idx, segments)
        })
        .collect()
}

fn metadata_snapshot(
    metadata: &[raptor::graph::navigation::PathMetadata],
) -> Vec<(String, usize, usize, bool)> {
    metadata
        .iter()
        .map(|m| {
            (
                m.id.clone(),
                m.segment_count,
                m.unique_segment_count,
                m.has_inversions,
            )
        })
        .collect()
}

proptest! {
    #[test]
    fn parse_gfa_paths_is_permutation_invariant_for_unique_path_ids(
        path_segments in prop::collection::vec(
            prop::collection::vec((0u8..16, any::<bool>()), 1..8),
            1..32
        ),
        seed in any::<u64>()
    ) {
        let lines = build_path_lines(&path_segments);
        let baseline = parse_gfa_paths(&lines);

        let mut shuffled = lines.clone();
        let mut rng = StdRng::seed_from_u64(seed);
        shuffled.shuffle(&mut rng);
        let observed = parse_gfa_paths(&shuffled);

        prop_assert_eq!(&observed, &baseline);

        let reader_observed = parse_gfa_paths_reader(Cursor::new(lines.join("\n")))
            .expect("reader parser should succeed for generated valid GFA P-lines");
        prop_assert_eq!(&reader_observed, &baseline);
    }

    #[test]
    fn extract_path_metadata_is_sorted_and_order_independent(
        path_segments in prop::collection::vec(
            prop::collection::vec((0u8..16, any::<bool>()), 1..8),
            1..32
        ),
        seed in any::<u64>()
    ) {
        let lines = build_path_lines(&path_segments);
        let baseline_map = parse_gfa_paths(&lines);
        let baseline_metadata = extract_path_metadata(&baseline_map);
        let baseline_snapshot = metadata_snapshot(&baseline_metadata);

        let ids: Vec<&str> = baseline_metadata.iter().map(|m| m.id.as_str()).collect();
        let mut expected_ids: Vec<String> = (0..path_segments.len())
            .map(|idx| format!("path_{}", idx))
            .collect();
        expected_ids.sort_unstable();
        prop_assert_eq!(
            ids,
            expected_ids.iter().map(String::as_str).collect::<Vec<_>>()
        );

        let mut entries: Vec<(String, Vec<(String, char)>)> = baseline_map.into_iter().collect();
        let mut rng = StdRng::seed_from_u64(seed);
        entries.shuffle(&mut rng);
        let mut shuffled_map = HashMap::new();
        for (id, segments) in entries {
            shuffled_map.insert(id, segments);
        }

        let shuffled_snapshot = metadata_snapshot(&extract_path_metadata(&shuffled_map));
        prop_assert_eq!(shuffled_snapshot, baseline_snapshot);
    }
}
