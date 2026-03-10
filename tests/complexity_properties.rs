use proptest::prelude::*;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use raptor::graph::complexity::compute_path_stats;
use std::collections::BTreeSet;
use std::io::Write;
use tempfile::NamedTempFile;

fn write_gfa(
    segments: &[String],
    links: &[(usize, usize)],
    paths: &[Vec<usize>],
    seed: u64,
) -> NamedTempFile {
    let mut records = Vec::new();

    for segment in segments {
        records.push(format!("S\t{}\tAAAA", segment));
    }

    for &(from, to) in links {
        records.push(format!("L\t{}\t+\t{}\t+\t0M", segments[from], segments[to]));
    }

    for (idx, path) in paths.iter().enumerate() {
        let segment_list = path
            .iter()
            .map(|&node| format!("{}+", segments[node]))
            .collect::<Vec<_>>()
            .join(",");
        records.push(format!("P\tpath_{}\t{}\t*", idx, segment_list));
    }

    let mut rng = StdRng::seed_from_u64(seed);
    records.shuffle(&mut rng);

    let mut file = NamedTempFile::new().expect("create temp file");
    writeln!(file, "H\tVN:Z:1.0").expect("write header");
    for record in records {
        writeln!(file, "{}", record).expect("write record");
    }
    file
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn compute_path_stats_is_invariant_to_gfa_record_order(
        segment_count in 1usize..=6,
        link_items in prop::collection::vec((0usize..6, 0usize..6), 0..16),
        path_items in prop::collection::vec(prop::collection::vec(0usize..6, 1..6), 0..8),
        seed in any::<u64>()
    ) {
        let segments: Vec<String> = (0..segment_count).map(|i| (i + 1).to_string()).collect();

        let mut unique_links = BTreeSet::new();
        for (from, to) in link_items {
            if from < segment_count && to < segment_count && from != to {
                unique_links.insert((from, to));
            }
        }
        let links: Vec<(usize, usize)> = unique_links.into_iter().collect();

        let mut paths = Vec::new();
        for path in path_items {
            let filtered: Vec<usize> = path
                .into_iter()
                .filter(|&node| node < segment_count)
                .collect();
            if !filtered.is_empty() {
                paths.push(filtered);
            }
        }
        if paths.is_empty() {
            paths.push(vec![0]);
        }

        let file_a = write_gfa(&segments, &links, &paths, seed);
        let file_b = write_gfa(&segments, &links, &paths, seed ^ 0xD3E7_9B1A_5C4F_2081);

        let stats_a = compute_path_stats(file_a.path().to_str().expect("utf8 path"))
            .expect("compute path stats for order A");
        let stats_b = compute_path_stats(file_b.path().to_str().expect("utf8 path"))
            .expect("compute path stats for order B");

        prop_assert_eq!(stats_a.total_paths, stats_b.total_paths);
        prop_assert_eq!(stats_a.branch_count, stats_b.branch_count);
        prop_assert_eq!(stats_a.max_depth, stats_b.max_depth);
        prop_assert_eq!(stats_a.bubble_count, stats_b.bubble_count);
        prop_assert!((stats_a.average_length - stats_b.average_length).abs() < 1e-12);
        prop_assert!((stats_a.branchiness - stats_b.branchiness).abs() < 1e-12);
    }
}
