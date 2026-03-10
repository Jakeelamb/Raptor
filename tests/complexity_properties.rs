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

fn write_gfa_with_walks(
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
        let walk = path
            .iter()
            .map(|&node| format!(">{}", segments[node]))
            .collect::<Vec<_>>()
            .join("");
        records.push(format!(
            "W\tsample\t{}\tchr1\t0\t{}\t{}",
            idx,
            path.len(),
            walk
        ));
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
        prop_assert!((stats_a.median_length - stats_b.median_length).abs() < 1e-12);
        prop_assert_eq!(stats_a.path_n50, stats_b.path_n50);
        prop_assert_eq!(stats_a.path_n90, stats_b.path_n90);
        prop_assert_eq!(stats_a.path_n95, stats_b.path_n95);
        prop_assert_eq!(stats_a.path_n99, stats_b.path_n99);
        prop_assert_eq!(stats_a.path_l50, stats_b.path_l50);
        prop_assert_eq!(stats_a.path_l90, stats_b.path_l90);
        prop_assert_eq!(stats_a.path_l95, stats_b.path_l95);
        prop_assert_eq!(stats_a.path_l99, stats_b.path_l99);
        prop_assert!((stats_a.path_au_n - stats_b.path_au_n).abs() < 1e-12);
        prop_assert!((stats_a.branchiness - stats_b.branchiness).abs() < 1e-12);
    }

    #[test]
    fn compute_path_stats_is_invariant_to_duplicate_links(
        segment_count in 1usize..=6,
        link_items in prop::collection::vec((0usize..6, 0usize..6), 0..16),
        path_items in prop::collection::vec(prop::collection::vec(0usize..6, 1..6), 0..8),
        seed in any::<u64>()
    ) {
        let segments: Vec<String> = (0..segment_count).map(|i| (i + 1).to_string()).collect();

        let mut links_with_duplicates = Vec::new();
        let mut unique_links = BTreeSet::new();
        for (from, to) in link_items {
            if from < segment_count && to < segment_count && from != to {
                links_with_duplicates.push((from, to));
                if ((from + to) & 1) == 0 {
                    links_with_duplicates.push((from, to));
                }
                unique_links.insert((from, to));
            }
        }
        let dedup_links: Vec<(usize, usize)> = unique_links.into_iter().collect();

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

        let file_with_duplicates = write_gfa(&segments, &links_with_duplicates, &paths, seed);
        let file_deduped = write_gfa(&segments, &dedup_links, &paths, seed ^ 0xA73C_91E4_6B2F_D0C8);

        let stats_with_duplicates = compute_path_stats(
            file_with_duplicates.path().to_str().expect("utf8 path")
        ).expect("compute stats with duplicate links");
        let stats_deduped = compute_path_stats(
            file_deduped.path().to_str().expect("utf8 path")
        ).expect("compute stats with deduplicated links");

        prop_assert_eq!(stats_with_duplicates.total_paths, stats_deduped.total_paths);
        prop_assert_eq!(stats_with_duplicates.branch_count, stats_deduped.branch_count);
        prop_assert_eq!(stats_with_duplicates.max_depth, stats_deduped.max_depth);
        prop_assert_eq!(stats_with_duplicates.bubble_count, stats_deduped.bubble_count);
        prop_assert!((stats_with_duplicates.average_length - stats_deduped.average_length).abs() < 1e-12);
        prop_assert!((stats_with_duplicates.median_length - stats_deduped.median_length).abs() < 1e-12);
        prop_assert_eq!(stats_with_duplicates.path_n50, stats_deduped.path_n50);
        prop_assert_eq!(stats_with_duplicates.path_n90, stats_deduped.path_n90);
        prop_assert_eq!(stats_with_duplicates.path_n95, stats_deduped.path_n95);
        prop_assert_eq!(stats_with_duplicates.path_n99, stats_deduped.path_n99);
        prop_assert_eq!(stats_with_duplicates.path_l50, stats_deduped.path_l50);
        prop_assert_eq!(stats_with_duplicates.path_l90, stats_deduped.path_l90);
        prop_assert_eq!(stats_with_duplicates.path_l95, stats_deduped.path_l95);
        prop_assert_eq!(stats_with_duplicates.path_l99, stats_deduped.path_l99);
        prop_assert!((stats_with_duplicates.path_au_n - stats_deduped.path_au_n).abs() < 1e-12);
        prop_assert!((stats_with_duplicates.branchiness - stats_deduped.branchiness).abs() < 1e-12);
    }

    #[test]
    fn compute_path_stats_matches_equivalent_p_and_w_path_records(
        segment_count in 1usize..=6,
        link_items in prop::collection::vec((0usize..6, 0usize..6), 0..16),
        path_items in prop::collection::vec(prop::collection::vec(0usize..6, 1..6), 1..8),
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

        let file_p = write_gfa(&segments, &links, &paths, seed);
        let file_w = write_gfa_with_walks(&segments, &links, &paths, seed ^ 0x9A5C_72D4_18EF_3341);

        let stats_p = compute_path_stats(file_p.path().to_str().expect("utf8 path"))
            .expect("compute stats from P records");
        let stats_w = compute_path_stats(file_w.path().to_str().expect("utf8 path"))
            .expect("compute stats from W records");

        prop_assert_eq!(stats_p.total_paths, stats_w.total_paths);
        prop_assert_eq!(stats_p.branch_count, stats_w.branch_count);
        prop_assert_eq!(stats_p.max_depth, stats_w.max_depth);
        prop_assert_eq!(stats_p.bubble_count, stats_w.bubble_count);
        prop_assert!((stats_p.average_length - stats_w.average_length).abs() < 1e-12);
        prop_assert!((stats_p.median_length - stats_w.median_length).abs() < 1e-12);
        prop_assert_eq!(stats_p.path_n50, stats_w.path_n50);
        prop_assert_eq!(stats_p.path_n90, stats_w.path_n90);
        prop_assert_eq!(stats_p.path_n95, stats_w.path_n95);
        prop_assert_eq!(stats_p.path_n99, stats_w.path_n99);
        prop_assert_eq!(stats_p.path_l50, stats_w.path_l50);
        prop_assert_eq!(stats_p.path_l90, stats_w.path_l90);
        prop_assert_eq!(stats_p.path_l95, stats_w.path_l95);
        prop_assert_eq!(stats_p.path_l99, stats_w.path_l99);
        prop_assert!((stats_p.path_au_n - stats_w.path_au_n).abs() < 1e-12);
        prop_assert!((stats_p.branchiness - stats_w.branchiness).abs() < 1e-12);
    }

    #[test]
    fn compute_path_stats_counts_long_reconverging_bubbles(
        depth in 11usize..40,
        seed in any::<u64>()
    ) {
        let reconverge = (2 * depth) + 1;
        let segment_count = reconverge + 1;
        let segments: Vec<String> = (0..segment_count).map(|i| i.to_string()).collect();

        let mut links = Vec::new();
        // Branch A: 0 -> 1 -> ... -> depth -> reconverge
        links.push((0, 1));
        for node in 1..depth {
            links.push((node, node + 1));
        }
        links.push((depth, reconverge));

        // Branch B: 0 -> (depth + 1) -> ... -> (2 * depth) -> reconverge
        let second_start = depth + 1;
        links.push((0, second_start));
        for node in second_start..(2 * depth) {
            links.push((node, node + 1));
        }
        links.push((2 * depth, reconverge));

        let path_a: Vec<usize> = std::iter::once(0)
            .chain(1..=depth)
            .chain(std::iter::once(reconverge))
            .collect();
        let path_b: Vec<usize> = std::iter::once(0)
            .chain(second_start..=(2 * depth))
            .chain(std::iter::once(reconverge))
            .collect();

        let file = write_gfa(&segments, &links, &[path_a, path_b], seed);
        let stats = compute_path_stats(file.path().to_str().expect("utf8 path"))
            .expect("compute path stats for long reconverging bubble");

        prop_assert_eq!(stats.bubble_count, 1);
    }
}
