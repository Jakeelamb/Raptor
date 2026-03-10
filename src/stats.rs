use crate::eval::metrics::{evaluate_lengths_in_place, BaseComposition};
use crate::io::fasta::try_open_fasta;
use serde::Serialize;
use std::io::BufRead;

#[derive(Serialize)]
pub struct Stats {
    pub total_contigs: usize,
    pub total_length: usize,
    pub average_length: f64,
    pub median_length: f64,
    pub gc_bases: usize,
    pub acgt_bases: usize,
    pub n_bases: usize,
    pub ambiguous_bases: usize,
    pub gc_content: f64,
    pub n_content: f64,
    pub ambiguous_content: f64,
    pub n25: usize,
    pub n50: usize,
    pub n75: usize,
    pub n90: usize,
    pub n95: usize,
    pub n99: usize,
    pub l25: usize,
    pub l50: usize,
    pub l75: usize,
    pub l90: usize,
    pub l95: usize,
    pub l99: usize,
    pub au_n: f64,
    pub longest_contig: usize,
    pub contigs_ge_1kb: usize,
    pub contigs_ge_10kb: usize,
    pub contigs_ge_50kb: usize,
    pub contigs_ge_100kb: usize,
    pub bases_ge_1kb: usize,
    pub bases_ge_10kb: usize,
    pub bases_ge_50kb: usize,
    pub bases_ge_100kb: usize,
    pub contigs_ge_1kb_frac: f64,
    pub contigs_ge_10kb_frac: f64,
    pub contigs_ge_50kb_frac: f64,
    pub contigs_ge_100kb_frac: f64,
    pub bases_ge_1kb_frac: f64,
    pub bases_ge_10kb_frac: f64,
    pub bases_ge_50kb_frac: f64,
    pub bases_ge_100kb_frac: f64,
    // Graph-related stats
    pub path_count: Option<usize>,
    pub avg_path_length: Option<f64>,
    pub branch_count: Option<usize>,
    pub graph_max_depth: Option<usize>,
    pub graph_bubble_count: Option<usize>,
    pub graph_branchiness: Option<f64>,
    pub graph_path_median_length: Option<f64>,
    pub graph_path_n50: Option<usize>,
    pub graph_path_n90: Option<usize>,
    pub graph_path_au_n: Option<f64>,
}

pub fn calculate_stats(path: &str) -> std::io::Result<Stats> {
    let reader = try_open_fasta(path)?;
    let mut lengths = vec![];
    let mut in_sequence = false;
    let mut current_len = 0usize;
    let mut composition = BaseComposition::default();

    for line_result in reader.lines() {
        let line = line_result?;
        if line.starts_with('>') {
            // If we were in a sequence, add its final length.
            if in_sequence {
                lengths.push(current_len);
            }
            in_sequence = true;
            current_len = 0;
        } else if in_sequence {
            let seq = line.trim().as_bytes();
            current_len += seq.len();
            composition.add_sequence(seq);
        }
    }

    // Add the last sequence if there is one
    if in_sequence {
        lengths.push(current_len);
    }

    let length_stats = evaluate_lengths_in_place(&mut lengths);
    let gc_content = composition.gc_content();
    let n_content = composition.n_content(length_stats.total_bases);
    let ambiguous_content = composition.ambiguous_content(length_stats.total_bases);

    Ok(Stats {
        total_contigs: length_stats.total,
        total_length: length_stats.total_bases,
        average_length: length_stats.avg_length,
        median_length: length_stats.median_length,
        gc_bases: composition.gc_bases,
        acgt_bases: composition.acgt_bases,
        n_bases: composition.n_bases,
        ambiguous_bases: composition.ambiguous_bases,
        gc_content,
        n_content,
        ambiguous_content,
        n25: length_stats.n25,
        n50: length_stats.n50,
        n75: length_stats.n75,
        n90: length_stats.n90,
        n95: length_stats.n95,
        n99: length_stats.n99,
        l25: length_stats.l25,
        l50: length_stats.l50,
        l75: length_stats.l75,
        l90: length_stats.l90,
        l95: length_stats.l95,
        l99: length_stats.l99,
        au_n: length_stats.au_n,
        longest_contig: length_stats.longest,
        contigs_ge_1kb: length_stats.contigs_ge_1kb,
        contigs_ge_10kb: length_stats.contigs_ge_10kb,
        contigs_ge_50kb: length_stats.contigs_ge_50kb,
        contigs_ge_100kb: length_stats.contigs_ge_100kb,
        bases_ge_1kb: length_stats.bases_ge_1kb,
        bases_ge_10kb: length_stats.bases_ge_10kb,
        bases_ge_50kb: length_stats.bases_ge_50kb,
        bases_ge_100kb: length_stats.bases_ge_100kb,
        contigs_ge_1kb_frac: length_stats.contigs_ge_1kb_frac,
        contigs_ge_10kb_frac: length_stats.contigs_ge_10kb_frac,
        contigs_ge_50kb_frac: length_stats.contigs_ge_50kb_frac,
        contigs_ge_100kb_frac: length_stats.contigs_ge_100kb_frac,
        bases_ge_1kb_frac: length_stats.bases_ge_1kb_frac,
        bases_ge_10kb_frac: length_stats.bases_ge_10kb_frac,
        bases_ge_50kb_frac: length_stats.bases_ge_50kb_frac,
        bases_ge_100kb_frac: length_stats.bases_ge_100kb_frac,
        path_count: None,
        avg_path_length: None,
        branch_count: None,
        graph_max_depth: None,
        graph_bubble_count: None,
        graph_branchiness: None,
        graph_path_median_length: None,
        graph_path_n50: None,
        graph_path_n90: None,
        graph_path_au_n: None,
    })
}

/// Generate and display graph complexity stats from a GFA file
pub fn calculate_graph_stats(path: &str) -> Option<crate::graph::complexity::PathStats> {
    // Check if the file exists and has a GFA extension
    if !path.to_lowercase().ends_with(".gfa") {
        return None;
    }

    // Use the new compute_path_stats function that takes a GFA file path
    match crate::graph::complexity::compute_path_stats(path) {
        Ok(stats) => Some(stats),
        Err(e) => {
            eprintln!("Error computing graph stats: {}", e);
            None
        }
    }
}

/// Update Stats object with graph complexity information
pub fn update_with_graph_stats(
    stats: &mut Stats,
    graph_stats: &crate::graph::complexity::PathStats,
) {
    stats.path_count = Some(graph_stats.total_paths);
    stats.avg_path_length = Some(graph_stats.average_length);
    stats.branch_count = Some(graph_stats.branch_count);
    stats.graph_max_depth = Some(graph_stats.max_depth);
    stats.graph_bubble_count = Some(graph_stats.bubble_count);
    stats.graph_branchiness = Some(graph_stats.branchiness);
    stats.graph_path_median_length = Some(graph_stats.median_length);
    stats.graph_path_n50 = Some(graph_stats.path_n50);
    stats.graph_path_n90 = Some(graph_stats.path_n90);
    stats.graph_path_au_n = Some(graph_stats.path_au_n);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::{NamedTempFile, TempDir};

    #[test]
    fn test_calculate_stats() {
        // Create a temporary FASTA file
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">contig_1").unwrap();
        writeln!(file, "ATCGATCGATCGATCGATCG").unwrap(); // 20 bp
        writeln!(file, ">contig_2").unwrap();
        writeln!(file, "GCTAGCTAGCTAGCTAGCTAGCTA").unwrap(); // 24 bp
        writeln!(file, ">contig_3").unwrap();
        writeln!(file, "ATCG").unwrap(); // 4 bp

        let stats = calculate_stats(file.path().to_str().unwrap()).unwrap();

        assert_eq!(stats.total_contigs, 3);
        assert_eq!(stats.total_length, 48);
        assert_eq!(stats.average_length, 16.0);
        assert_eq!(stats.median_length, 20.0);
        assert_eq!(stats.gc_bases, 24);
        assert_eq!(stats.acgt_bases, 48);
        assert_eq!(stats.n_bases, 0);
        assert_eq!(stats.ambiguous_bases, 0);
        assert!((stats.gc_content - 0.5).abs() < 1e-12);
        assert_eq!(stats.n_content, 0.0);
        assert_eq!(stats.ambiguous_content, 0.0);
        assert_eq!(stats.n25, 24);
        assert_eq!(stats.n50, 24);
        assert_eq!(stats.n75, 20);
        assert_eq!(stats.n90, 20);
        assert_eq!(stats.n95, 4);
        assert_eq!(stats.n99, 4);
        assert_eq!(stats.l25, 1);
        assert_eq!(stats.l50, 1);
        assert_eq!(stats.l75, 2);
        assert_eq!(stats.l90, 2);
        assert_eq!(stats.l95, 3);
        assert_eq!(stats.l99, 3);
        assert!((stats.au_n - 20.6666666667).abs() < 1e-6);
        assert_eq!(stats.longest_contig, 24);
        assert_eq!(stats.contigs_ge_1kb, 0);
        assert_eq!(stats.contigs_ge_10kb, 0);
        assert_eq!(stats.contigs_ge_50kb, 0);
        assert_eq!(stats.contigs_ge_100kb, 0);
        assert_eq!(stats.bases_ge_1kb, 0);
        assert_eq!(stats.bases_ge_10kb, 0);
        assert_eq!(stats.bases_ge_50kb, 0);
        assert_eq!(stats.bases_ge_100kb, 0);
        assert_eq!(stats.contigs_ge_1kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_10kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_50kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_100kb_frac, 0.0);
        assert_eq!(stats.bases_ge_1kb_frac, 0.0);
        assert_eq!(stats.bases_ge_10kb_frac, 0.0);
        assert_eq!(stats.bases_ge_50kb_frac, 0.0);
        assert_eq!(stats.bases_ge_100kb_frac, 0.0);
    }

    #[test]
    fn test_calculate_stats_multiline_sequences() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">contig_1").unwrap();
        writeln!(file, "ATCGATCG").unwrap();
        writeln!(file, "ATCG").unwrap(); // 12 bp total
        writeln!(file, ">contig_2").unwrap();
        writeln!(file, "GCTA").unwrap(); // 4 bp

        let stats = calculate_stats(file.path().to_str().unwrap()).unwrap();
        assert_eq!(stats.total_contigs, 2);
        assert_eq!(stats.total_length, 16);
        assert_eq!(stats.median_length, 8.0);
        assert_eq!(stats.gc_bases, 8);
        assert_eq!(stats.acgt_bases, 16);
        assert_eq!(stats.n_bases, 0);
        assert_eq!(stats.ambiguous_bases, 0);
        assert!((stats.gc_content - 0.5).abs() < 1e-12);
        assert_eq!(stats.n_content, 0.0);
        assert_eq!(stats.ambiguous_content, 0.0);
        assert_eq!(stats.n25, 12);
        assert_eq!(stats.n50, 12);
        assert_eq!(stats.n75, 12);
        assert_eq!(stats.n90, 4);
        assert_eq!(stats.n95, 4);
        assert_eq!(stats.n99, 4);
        assert_eq!(stats.l25, 1);
        assert_eq!(stats.l50, 1);
        assert_eq!(stats.l75, 1);
        assert_eq!(stats.l90, 2);
        assert_eq!(stats.l95, 2);
        assert_eq!(stats.l99, 2);
        assert!((stats.au_n - 10.0).abs() < 1e-6);
        assert_eq!(stats.longest_contig, 12);
        assert_eq!(stats.contigs_ge_1kb, 0);
        assert_eq!(stats.contigs_ge_10kb, 0);
        assert_eq!(stats.contigs_ge_50kb, 0);
        assert_eq!(stats.contigs_ge_100kb, 0);
        assert_eq!(stats.bases_ge_1kb, 0);
        assert_eq!(stats.bases_ge_10kb, 0);
        assert_eq!(stats.bases_ge_50kb, 0);
        assert_eq!(stats.bases_ge_100kb, 0);
        assert_eq!(stats.contigs_ge_1kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_10kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_50kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_100kb_frac, 0.0);
        assert_eq!(stats.bases_ge_1kb_frac, 0.0);
        assert_eq!(stats.bases_ge_10kb_frac, 0.0);
        assert_eq!(stats.bases_ge_50kb_frac, 0.0);
        assert_eq!(stats.bases_ge_100kb_frac, 0.0);
    }

    #[test]
    fn test_calculate_stats_reports_length_bucket_metrics() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">contig_1").unwrap();
        writeln!(file, "{}", "A".repeat(50_000)).unwrap();
        writeln!(file, ">contig_2").unwrap();
        writeln!(file, "{}", "C".repeat(10_000)).unwrap();
        writeln!(file, ">contig_3").unwrap();
        writeln!(file, "{}", "G".repeat(1_000)).unwrap();
        writeln!(file, ">contig_4").unwrap();
        writeln!(file, "{}", "T".repeat(999)).unwrap();

        let stats = calculate_stats(file.path().to_str().unwrap()).unwrap();
        assert_eq!(stats.contigs_ge_1kb, 3);
        assert_eq!(stats.contigs_ge_10kb, 2);
        assert_eq!(stats.contigs_ge_50kb, 1);
        assert_eq!(stats.contigs_ge_100kb, 0);
        assert_eq!(stats.bases_ge_1kb, 61_000);
        assert_eq!(stats.bases_ge_10kb, 60_000);
        assert_eq!(stats.bases_ge_50kb, 50_000);
        assert_eq!(stats.bases_ge_100kb, 0);
        assert!((stats.contigs_ge_1kb_frac - 0.75).abs() < 1e-12);
        assert!((stats.contigs_ge_10kb_frac - 0.5).abs() < 1e-12);
        assert!((stats.contigs_ge_50kb_frac - 0.25).abs() < 1e-12);
        assert_eq!(stats.contigs_ge_100kb_frac, 0.0);
        assert!((stats.bases_ge_1kb_frac - (61_000.0 / 61_999.0)).abs() < 1e-12);
        assert!((stats.bases_ge_10kb_frac - (60_000.0 / 61_999.0)).abs() < 1e-12);
        assert!((stats.bases_ge_50kb_frac - (50_000.0 / 61_999.0)).abs() < 1e-12);
        assert_eq!(stats.bases_ge_100kb_frac, 0.0);
    }

    #[test]
    fn test_calculate_stats_reports_100kb_bucket_metrics() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">contig_1").unwrap();
        writeln!(file, "{}", "A".repeat(100_000)).unwrap();
        writeln!(file, ">contig_2").unwrap();
        writeln!(file, "{}", "C".repeat(99_999)).unwrap();

        let stats = calculate_stats(file.path().to_str().unwrap()).unwrap();
        assert_eq!(stats.contigs_ge_100kb, 1);
        assert_eq!(stats.bases_ge_100kb, 100_000);
        assert!((stats.contigs_ge_100kb_frac - 0.5).abs() < 1e-12);
        assert!((stats.bases_ge_100kb_frac - (100_000.0 / 199_999.0)).abs() < 1e-12);
    }

    #[test]
    fn test_update_with_graph_stats_populates_extended_metrics() {
        let mut stats = Stats {
            total_contigs: 0,
            total_length: 0,
            average_length: 0.0,
            median_length: 0.0,
            gc_bases: 0,
            acgt_bases: 0,
            n_bases: 0,
            ambiguous_bases: 0,
            gc_content: 0.0,
            n_content: 0.0,
            ambiguous_content: 0.0,
            n25: 0,
            n50: 0,
            n75: 0,
            n90: 0,
            n95: 0,
            n99: 0,
            l25: 0,
            l50: 0,
            l75: 0,
            l90: 0,
            l95: 0,
            l99: 0,
            au_n: 0.0,
            longest_contig: 0,
            contigs_ge_1kb: 0,
            contigs_ge_10kb: 0,
            contigs_ge_50kb: 0,
            contigs_ge_100kb: 0,
            bases_ge_1kb: 0,
            bases_ge_10kb: 0,
            bases_ge_50kb: 0,
            bases_ge_100kb: 0,
            contigs_ge_1kb_frac: 0.0,
            contigs_ge_10kb_frac: 0.0,
            contigs_ge_50kb_frac: 0.0,
            contigs_ge_100kb_frac: 0.0,
            bases_ge_1kb_frac: 0.0,
            bases_ge_10kb_frac: 0.0,
            bases_ge_50kb_frac: 0.0,
            bases_ge_100kb_frac: 0.0,
            path_count: None,
            avg_path_length: None,
            branch_count: None,
            graph_max_depth: None,
            graph_bubble_count: None,
            graph_branchiness: None,
            graph_path_median_length: None,
            graph_path_n50: None,
            graph_path_n90: None,
            graph_path_au_n: None,
        };
        let graph_stats = crate::graph::complexity::PathStats {
            total_paths: 7,
            average_length: 3.5,
            median_length: 3.0,
            path_n50: 4,
            path_n90: 2,
            path_au_n: 3.2,
            branch_count: 2,
            max_depth: 5,
            bubble_count: 1,
            branchiness: 0.4,
        };

        update_with_graph_stats(&mut stats, &graph_stats);

        assert_eq!(stats.path_count, Some(7));
        assert_eq!(stats.avg_path_length, Some(3.5));
        assert_eq!(stats.branch_count, Some(2));
        assert_eq!(stats.graph_max_depth, Some(5));
        assert_eq!(stats.graph_bubble_count, Some(1));
        assert_eq!(stats.graph_branchiness, Some(0.4));
        assert_eq!(stats.graph_path_median_length, Some(3.0));
        assert_eq!(stats.graph_path_n50, Some(4));
        assert_eq!(stats.graph_path_n90, Some(2));
        assert_eq!(stats.graph_path_au_n, Some(3.2));
    }

    #[test]
    fn test_calculate_stats_reports_base_composition_metrics() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">contig_1").unwrap();
        writeln!(file, "GCGCNNNN").unwrap(); // 4 GC, 4 N
        writeln!(file, ">contig_2").unwrap();
        writeln!(file, "atuy").unwrap(); // 1 A, 1 T, 1 U, 1 ambiguous

        let stats = calculate_stats(file.path().to_str().unwrap()).unwrap();
        assert_eq!(stats.total_length, 12);
        assert_eq!(stats.gc_bases, 4);
        assert_eq!(stats.acgt_bases, 7);
        assert_eq!(stats.n_bases, 4);
        assert_eq!(stats.ambiguous_bases, 1);
        assert!((stats.gc_content - (4.0 / 7.0)).abs() < 1e-12);
        assert!((stats.n_content - (4.0 / 12.0)).abs() < 1e-12);
        assert!((stats.ambiguous_content - (1.0 / 12.0)).abs() < 1e-12);
    }

    #[test]
    fn test_calculate_stats_counts_zero_length_contigs() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">contig_1").unwrap();
        writeln!(file, "AT").unwrap();
        writeln!(file, ">contig_2").unwrap(); // empty contig
        writeln!(file, ">contig_3").unwrap();
        writeln!(file, "G").unwrap();
        writeln!(file, ">contig_4").unwrap(); // empty contig at EOF

        let stats = calculate_stats(file.path().to_str().unwrap()).unwrap();
        assert_eq!(stats.total_contigs, 4);
        assert_eq!(stats.total_length, 3);
        assert!((stats.average_length - 0.75).abs() < 1e-12);
        assert_eq!(stats.median_length, 0.5);
        assert_eq!(stats.n25, 2);
        assert_eq!(stats.n50, 2);
        assert_eq!(stats.n75, 1);
        assert_eq!(stats.n90, 1);
        assert_eq!(stats.n95, 1);
        assert_eq!(stats.n99, 1);
        assert_eq!(stats.l25, 1);
        assert_eq!(stats.l50, 1);
        assert_eq!(stats.l75, 2);
        assert_eq!(stats.l90, 2);
        assert_eq!(stats.l95, 2);
        assert_eq!(stats.l99, 2);
    }

    #[test]
    fn test_calculate_stats_returns_not_found_for_missing_fasta() {
        let temp_dir = TempDir::new().unwrap();
        let missing = temp_dir.path().join("missing.fasta");

        match calculate_stats(missing.to_str().unwrap()) {
            Ok(_) => panic!("expected not found error for missing FASTA"),
            Err(err) => assert_eq!(err.kind(), std::io::ErrorKind::NotFound),
        }
    }
}
