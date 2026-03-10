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
    pub contigs_with_n: usize,
    pub contigs_with_ambiguous: usize,
    pub contigs_all_acgt: usize,
    pub gc_content: f64,
    pub n_content: f64,
    pub ambiguous_content: f64,
    pub contigs_with_n_frac: f64,
    pub contigs_with_ambiguous_frac: f64,
    pub contigs_all_acgt_frac: f64,
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

#[derive(Debug, Default, Clone, Copy)]
struct ContigQualitySummary {
    with_n: usize,
    with_ambiguous: usize,
    all_acgt: usize,
}

#[inline]
fn update_contig_quality(summary: &mut ContigQualitySummary, composition: BaseComposition) {
    if composition.n_bases > 0 {
        summary.with_n += 1;
    }
    if composition.ambiguous_bases > 0 {
        summary.with_ambiguous += 1;
    }
    if composition.n_bases == 0 && composition.ambiguous_bases == 0 {
        summary.all_acgt += 1;
    }
}

pub fn calculate_stats(path: &str) -> std::io::Result<Stats> {
    let mut reader = try_open_fasta(path)?;
    let mut lengths = vec![];
    let mut in_sequence = false;
    let mut current_len = 0usize;
    let mut composition = BaseComposition::default();
    let mut current_contig_composition = BaseComposition::default();
    let mut contig_quality = ContigQualitySummary::default();
    let mut line = String::new();

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }

        let trimmed_line = line.trim_end_matches(['\n', '\r']);
        if trimmed_line.starts_with('>') {
            // If we were in a sequence, add its final length.
            if in_sequence {
                lengths.push(current_len);
                update_contig_quality(&mut contig_quality, current_contig_composition);
            }
            in_sequence = true;
            current_len = 0;
            current_contig_composition = BaseComposition::default();
        } else if in_sequence {
            let seq = trimmed_line.trim().as_bytes();
            current_len = current_len.checked_add(seq.len()).ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("contig length overflow while reading {}", path),
                )
            })?;
            composition.add_sequence(seq);
            current_contig_composition.add_sequence(seq);
        }
    }

    // Add the last sequence if there is one
    if in_sequence {
        lengths.push(current_len);
        update_contig_quality(&mut contig_quality, current_contig_composition);
    }

    let length_stats = evaluate_lengths_in_place(&mut lengths);
    let gc_content = composition.gc_content();
    let n_content = composition.n_content(length_stats.total_bases);
    let ambiguous_content = composition.ambiguous_content(length_stats.total_bases);
    let total_contigs = length_stats.total as f64;
    let (contigs_with_n_frac, contigs_with_ambiguous_frac, contigs_all_acgt_frac) =
        if total_contigs > 0.0 {
            (
                contig_quality.with_n as f64 / total_contigs,
                contig_quality.with_ambiguous as f64 / total_contigs,
                contig_quality.all_acgt as f64 / total_contigs,
            )
        } else {
            (0.0, 0.0, 0.0)
        };

    Ok(Stats {
        total_contigs: length_stats.total,
        total_length: length_stats.total_bases,
        average_length: length_stats.avg_length,
        median_length: length_stats.median_length,
        gc_bases: composition.gc_bases,
        acgt_bases: composition.acgt_bases,
        n_bases: composition.n_bases,
        ambiguous_bases: composition.ambiguous_bases,
        contigs_with_n: contig_quality.with_n,
        contigs_with_ambiguous: contig_quality.with_ambiguous,
        contigs_all_acgt: contig_quality.all_acgt,
        gc_content,
        n_content,
        ambiguous_content,
        contigs_with_n_frac,
        contigs_with_ambiguous_frac,
        contigs_all_acgt_frac,
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
    use proptest::prelude::*;
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
        assert_eq!(stats.contigs_with_n, 0);
        assert_eq!(stats.contigs_with_ambiguous, 0);
        assert_eq!(stats.contigs_all_acgt, 3);
        assert!((stats.gc_content - 0.5).abs() < 1e-12);
        assert_eq!(stats.n_content, 0.0);
        assert_eq!(stats.ambiguous_content, 0.0);
        assert_eq!(stats.contigs_with_n_frac, 0.0);
        assert_eq!(stats.contigs_with_ambiguous_frac, 0.0);
        assert_eq!(stats.contigs_all_acgt_frac, 1.0);
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
        assert_eq!(stats.contigs_with_n, 0);
        assert_eq!(stats.contigs_with_ambiguous, 0);
        assert_eq!(stats.contigs_all_acgt, 2);
        assert!((stats.gc_content - 0.5).abs() < 1e-12);
        assert_eq!(stats.n_content, 0.0);
        assert_eq!(stats.ambiguous_content, 0.0);
        assert_eq!(stats.contigs_with_n_frac, 0.0);
        assert_eq!(stats.contigs_with_ambiguous_frac, 0.0);
        assert_eq!(stats.contigs_all_acgt_frac, 1.0);
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
            contigs_with_n: 0,
            contigs_with_ambiguous: 0,
            contigs_all_acgt: 0,
            gc_content: 0.0,
            n_content: 0.0,
            ambiguous_content: 0.0,
            contigs_with_n_frac: 0.0,
            contigs_with_ambiguous_frac: 0.0,
            contigs_all_acgt_frac: 0.0,
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
        assert_eq!(stats.contigs_with_n, 1);
        assert_eq!(stats.contigs_with_ambiguous, 1);
        assert_eq!(stats.contigs_all_acgt, 0);
        assert!((stats.gc_content - (4.0 / 7.0)).abs() < 1e-12);
        assert!((stats.n_content - (4.0 / 12.0)).abs() < 1e-12);
        assert!((stats.ambiguous_content - (1.0 / 12.0)).abs() < 1e-12);
        assert!((stats.contigs_with_n_frac - 0.5).abs() < 1e-12);
        assert!((stats.contigs_with_ambiguous_frac - 0.5).abs() < 1e-12);
        assert_eq!(stats.contigs_all_acgt_frac, 0.0);
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
        assert_eq!(stats.contigs_with_n, 0);
        assert_eq!(stats.contigs_with_ambiguous, 0);
        assert_eq!(stats.contigs_all_acgt, 4);
        assert_eq!(stats.contigs_with_n_frac, 0.0);
        assert_eq!(stats.contigs_with_ambiguous_frac, 0.0);
        assert_eq!(stats.contigs_all_acgt_frac, 1.0);
    }

    #[test]
    fn test_calculate_stats_reports_contig_ambiguity_metrics() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">contig_1").unwrap();
        writeln!(file, "ACGT").unwrap();
        writeln!(file, ">contig_2").unwrap();
        writeln!(file, "NN").unwrap();
        writeln!(file, ">contig_3").unwrap();
        writeln!(file, "ARY").unwrap();

        let stats = calculate_stats(file.path().to_str().unwrap()).unwrap();
        assert_eq!(stats.total_contigs, 3);
        assert_eq!(stats.contigs_with_n, 1);
        assert_eq!(stats.contigs_with_ambiguous, 1);
        assert_eq!(stats.contigs_all_acgt, 1);
        assert!((stats.contigs_with_n_frac - (1.0 / 3.0)).abs() < 1e-12);
        assert!((stats.contigs_with_ambiguous_frac - (1.0 / 3.0)).abs() < 1e-12);
        assert!((stats.contigs_all_acgt_frac - (1.0 / 3.0)).abs() < 1e-12);
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

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(64))]

        #[test]
        fn calculate_stats_is_invariant_to_sequence_line_wrapping(
            contigs in prop::collection::vec("[ACGTNacgtnYRyr]{0,80}", 1..12)
        ) {
            let mut single_line = NamedTempFile::new().unwrap();
            let mut wrapped = NamedTempFile::new().unwrap();

            for (idx, sequence) in contigs.iter().enumerate() {
                writeln!(single_line, ">contig_{}", idx).unwrap();
                writeln!(single_line, "{}", sequence).unwrap();

                writeln!(wrapped, ">contig_{}", idx).unwrap();
                let width = (idx % 7) + 1;
                if sequence.is_empty() {
                    writeln!(wrapped).unwrap();
                } else {
                    for chunk in sequence.as_bytes().chunks(width) {
                        writeln!(wrapped, "{}", std::str::from_utf8(chunk).unwrap()).unwrap();
                    }
                }
            }

            let stats_a = calculate_stats(single_line.path().to_str().unwrap()).unwrap();
            let stats_b = calculate_stats(wrapped.path().to_str().unwrap()).unwrap();

            prop_assert_eq!(stats_a.total_contigs, stats_b.total_contigs);
            prop_assert_eq!(stats_a.total_length, stats_b.total_length);
            prop_assert_eq!(stats_a.gc_bases, stats_b.gc_bases);
            prop_assert_eq!(stats_a.acgt_bases, stats_b.acgt_bases);
            prop_assert_eq!(stats_a.n_bases, stats_b.n_bases);
            prop_assert_eq!(stats_a.ambiguous_bases, stats_b.ambiguous_bases);
            prop_assert_eq!(stats_a.contigs_with_n, stats_b.contigs_with_n);
            prop_assert_eq!(stats_a.contigs_with_ambiguous, stats_b.contigs_with_ambiguous);
            prop_assert_eq!(stats_a.contigs_all_acgt, stats_b.contigs_all_acgt);
            prop_assert_eq!(stats_a.n25, stats_b.n25);
            prop_assert_eq!(stats_a.n50, stats_b.n50);
            prop_assert_eq!(stats_a.n75, stats_b.n75);
            prop_assert_eq!(stats_a.n90, stats_b.n90);
            prop_assert_eq!(stats_a.n95, stats_b.n95);
            prop_assert_eq!(stats_a.n99, stats_b.n99);
            prop_assert_eq!(stats_a.l25, stats_b.l25);
            prop_assert_eq!(stats_a.l50, stats_b.l50);
            prop_assert_eq!(stats_a.l75, stats_b.l75);
            prop_assert_eq!(stats_a.l90, stats_b.l90);
            prop_assert_eq!(stats_a.l95, stats_b.l95);
            prop_assert_eq!(stats_a.l99, stats_b.l99);
            prop_assert!((stats_a.average_length - stats_b.average_length).abs() < 1e-12);
            prop_assert!((stats_a.median_length - stats_b.median_length).abs() < 1e-12);
            prop_assert!((stats_a.gc_content - stats_b.gc_content).abs() < 1e-12);
            prop_assert!((stats_a.n_content - stats_b.n_content).abs() < 1e-12);
            prop_assert!((stats_a.ambiguous_content - stats_b.ambiguous_content).abs() < 1e-12);
            prop_assert!((stats_a.contigs_with_n_frac - stats_b.contigs_with_n_frac).abs() < 1e-12);
            prop_assert!((stats_a.contigs_with_ambiguous_frac - stats_b.contigs_with_ambiguous_frac).abs() < 1e-12);
            prop_assert!((stats_a.contigs_all_acgt_frac - stats_b.contigs_all_acgt_frac).abs() < 1e-12);
            prop_assert!((stats_a.au_n - stats_b.au_n).abs() < 1e-12);
        }
    }
}
