use crate::eval::metrics::{evaluate_lengths_in_place, normalized_run_base, BaseComposition};
use crate::io::fasta::try_open_fasta;
use serde::Serialize;
use std::io::BufRead;

const TSV_HEADER: &str = "contigs\ttotal_len\tavg_len\tmedian_len\tgc_bases\tacgt_bases\tn_bases\tambiguous_bases\tmean_rle_ratio\tlength_weighted_rle_ratio\ttotal_rle_runs\tgc_content\tn_content\tambiguous_content\tn10\tn25\tn50\tn75\tn90\tn95\tn99\tl10\tl25\tl50\tl75\tl90\tl95\tl99\taun\teffective_contig_count\tungapped_effective_contig_count\tungapped_total_len\tungapped_n50\tungapped_aun\tlongest\tcontigs_ge_1kb\tcontigs_ge_10kb\tcontigs_ge_50kb\tcontigs_ge_100kb\tcontigs_ge_1mb\tbases_ge_1kb\tbases_ge_10kb\tbases_ge_50kb\tbases_ge_100kb\tbases_ge_1mb\tcontigs_ge_1kb_frac\tcontigs_ge_10kb_frac\tcontigs_ge_50kb_frac\tcontigs_ge_100kb_frac\tcontigs_ge_1mb_frac\tbases_ge_1kb_frac\tbases_ge_10kb_frac\tbases_ge_50kb_frac\tbases_ge_100kb_frac\tbases_ge_1mb_frac\tn_run_count\tmax_n_run\tmean_n_run_length\tn_runs_per_100kb\tn_bases_per_100kb\tambiguous_bases_per_100kb\tcontigs_with_n\tcontigs_with_ambiguous\tcontigs_all_acgt\tcontigs_with_n_frac\tcontigs_with_ambiguous_frac\tcontigs_all_acgt_frac";

#[derive(Serialize)]
pub struct Stats {
    pub total_contigs: usize,
    pub total_length: usize,
    pub ungapped_total_length: usize,
    pub average_length: f64,
    pub median_length: f64,
    pub n10: usize,
    pub gc_bases: usize,
    pub acgt_bases: usize,
    pub n_bases: usize,
    pub ambiguous_bases: usize,
    pub mean_rle_ratio: f64,
    pub length_weighted_rle_ratio: f64,
    pub total_rle_runs: usize,
    pub n_run_count: usize,
    pub max_n_run: usize,
    pub mean_n_run_length: f64,
    pub n_runs_per_100kb: f64,
    pub n_bases_per_100kb: f64,
    pub ambiguous_bases_per_100kb: f64,
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
    pub ungapped_n50: usize,
    pub ungapped_au_n: f64,
    pub ungapped_effective_contig_count: f64,
    pub l10: usize,
    pub l25: usize,
    pub l50: usize,
    pub l75: usize,
    pub l90: usize,
    pub l95: usize,
    pub l99: usize,
    pub au_n: f64,
    pub effective_contig_count: f64,
    pub longest_contig: usize,
    pub contigs_ge_1kb: usize,
    pub contigs_ge_10kb: usize,
    pub contigs_ge_50kb: usize,
    pub contigs_ge_100kb: usize,
    pub contigs_ge_1mb: usize,
    pub bases_ge_1kb: usize,
    pub bases_ge_10kb: usize,
    pub bases_ge_50kb: usize,
    pub bases_ge_100kb: usize,
    pub bases_ge_1mb: usize,
    pub contigs_ge_1kb_frac: f64,
    pub contigs_ge_10kb_frac: f64,
    pub contigs_ge_50kb_frac: f64,
    pub contigs_ge_100kb_frac: f64,
    pub contigs_ge_1mb_frac: f64,
    pub bases_ge_1kb_frac: f64,
    pub bases_ge_10kb_frac: f64,
    pub bases_ge_50kb_frac: f64,
    pub bases_ge_100kb_frac: f64,
    pub bases_ge_1mb_frac: f64,
    // Graph-related stats
    pub path_count: Option<usize>,
    pub avg_path_length: Option<f64>,
    pub branch_count: Option<usize>,
    pub graph_max_depth: Option<usize>,
    pub graph_bubble_count: Option<usize>,
    pub graph_branchiness: Option<f64>,
    pub graph_path_median_length: Option<f64>,
    pub graph_path_n10: Option<usize>,
    pub graph_path_n25: Option<usize>,
    pub graph_path_n50: Option<usize>,
    pub graph_path_n75: Option<usize>,
    pub graph_path_n90: Option<usize>,
    pub graph_path_n95: Option<usize>,
    pub graph_path_n99: Option<usize>,
    pub graph_path_l10: Option<usize>,
    pub graph_path_l25: Option<usize>,
    pub graph_path_l50: Option<usize>,
    pub graph_path_l75: Option<usize>,
    pub graph_path_l90: Option<usize>,
    pub graph_path_l95: Option<usize>,
    pub graph_path_l99: Option<usize>,
    pub graph_path_au_n: Option<f64>,
    pub graph_path_effective_count: Option<f64>,
}

#[derive(Debug, Default, Clone, Copy)]
struct ContigQualitySummary {
    with_n: usize,
    with_ambiguous: usize,
    all_acgt: usize,
}

#[derive(Debug, Default, Clone, Copy)]
struct NRunSummary {
    run_count: usize,
    max_run: usize,
    active_run: usize,
}

impl NRunSummary {
    #[inline]
    fn finish_active_run(&mut self) {
        if self.active_run > 0 {
            self.run_count += 1;
            self.max_run = self.max_run.max(self.active_run);
            self.active_run = 0;
        }
    }

    #[inline]
    fn add_base(&mut self, base: u8) {
        if matches!(base, b'N' | b'n') {
            self.active_run += 1;
        } else {
            self.finish_active_run();
        }
    }

    #[inline]
    fn finish_contig(&mut self) {
        self.finish_active_run();
    }
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

#[inline]
fn per_100kb(count: usize, total_bases: usize) -> f64 {
    if total_bases == 0 {
        0.0
    } else {
        count as f64 * 100_000.0 / total_bases as f64
    }
}

#[inline]
pub fn tsv_header() -> &'static str {
    TSV_HEADER
}

pub fn tsv_row(stats: &Stats) -> String {
    use std::fmt::Write as _;

    let mut row = String::with_capacity(1024);
    macro_rules! push_fmt {
        ($($arg:tt)*) => {{
            if !row.is_empty() {
                row.push('\t');
            }
            write!(&mut row, $($arg)*).expect("writing to String cannot fail");
        }};
    }

    push_fmt!("{}", stats.total_contigs);
    push_fmt!("{}", stats.total_length);
    push_fmt!("{:.2}", stats.average_length);
    push_fmt!("{:.2}", stats.median_length);
    push_fmt!("{}", stats.gc_bases);
    push_fmt!("{}", stats.acgt_bases);
    push_fmt!("{}", stats.n_bases);
    push_fmt!("{}", stats.ambiguous_bases);
    push_fmt!("{:.6}", stats.mean_rle_ratio);
    push_fmt!("{:.6}", stats.length_weighted_rle_ratio);
    push_fmt!("{}", stats.total_rle_runs);
    push_fmt!("{:.6}", stats.gc_content);
    push_fmt!("{:.6}", stats.n_content);
    push_fmt!("{:.6}", stats.ambiguous_content);
    push_fmt!("{}", stats.n10);
    push_fmt!("{}", stats.n25);
    push_fmt!("{}", stats.n50);
    push_fmt!("{}", stats.n75);
    push_fmt!("{}", stats.n90);
    push_fmt!("{}", stats.n95);
    push_fmt!("{}", stats.n99);
    push_fmt!("{}", stats.l10);
    push_fmt!("{}", stats.l25);
    push_fmt!("{}", stats.l50);
    push_fmt!("{}", stats.l75);
    push_fmt!("{}", stats.l90);
    push_fmt!("{}", stats.l95);
    push_fmt!("{}", stats.l99);
    push_fmt!("{:.2}", stats.au_n);
    push_fmt!("{:.6}", stats.effective_contig_count);
    push_fmt!("{:.6}", stats.ungapped_effective_contig_count);
    push_fmt!("{}", stats.ungapped_total_length);
    push_fmt!("{}", stats.ungapped_n50);
    push_fmt!("{:.2}", stats.ungapped_au_n);
    push_fmt!("{}", stats.longest_contig);
    push_fmt!("{}", stats.contigs_ge_1kb);
    push_fmt!("{}", stats.contigs_ge_10kb);
    push_fmt!("{}", stats.contigs_ge_50kb);
    push_fmt!("{}", stats.contigs_ge_100kb);
    push_fmt!("{}", stats.contigs_ge_1mb);
    push_fmt!("{}", stats.bases_ge_1kb);
    push_fmt!("{}", stats.bases_ge_10kb);
    push_fmt!("{}", stats.bases_ge_50kb);
    push_fmt!("{}", stats.bases_ge_100kb);
    push_fmt!("{}", stats.bases_ge_1mb);
    push_fmt!("{:.6}", stats.contigs_ge_1kb_frac);
    push_fmt!("{:.6}", stats.contigs_ge_10kb_frac);
    push_fmt!("{:.6}", stats.contigs_ge_50kb_frac);
    push_fmt!("{:.6}", stats.contigs_ge_100kb_frac);
    push_fmt!("{:.6}", stats.contigs_ge_1mb_frac);
    push_fmt!("{:.6}", stats.bases_ge_1kb_frac);
    push_fmt!("{:.6}", stats.bases_ge_10kb_frac);
    push_fmt!("{:.6}", stats.bases_ge_50kb_frac);
    push_fmt!("{:.6}", stats.bases_ge_100kb_frac);
    push_fmt!("{:.6}", stats.bases_ge_1mb_frac);
    push_fmt!("{}", stats.n_run_count);
    push_fmt!("{}", stats.max_n_run);
    push_fmt!("{:.6}", stats.mean_n_run_length);
    push_fmt!("{:.6}", stats.n_runs_per_100kb);
    push_fmt!("{:.6}", stats.n_bases_per_100kb);
    push_fmt!("{:.6}", stats.ambiguous_bases_per_100kb);
    push_fmt!("{}", stats.contigs_with_n);
    push_fmt!("{}", stats.contigs_with_ambiguous);
    push_fmt!("{}", stats.contigs_all_acgt);
    push_fmt!("{:.6}", stats.contigs_with_n_frac);
    push_fmt!("{:.6}", stats.contigs_with_ambiguous_frac);
    push_fmt!("{:.6}", stats.contigs_all_acgt_frac);

    row
}

pub fn calculate_stats(path: &str) -> std::io::Result<Stats> {
    let mut reader = try_open_fasta(path)?;
    let mut lengths = vec![];
    let mut ungapped_lengths = vec![];
    let mut in_sequence = false;
    let mut current_len = 0usize;
    let mut current_ungapped_len = 0usize;
    let mut composition = BaseComposition::default();
    let mut current_contig_composition = BaseComposition::default();
    let mut contig_quality = ContigQualitySummary::default();
    let mut n_runs = NRunSummary::default();
    let mut total_rle_ratio = 0.0f64;
    let mut total_rle_runs = 0usize;
    let mut current_rle_len = 0usize;
    let mut current_rle_last_base: Option<u8> = None;
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
                ungapped_lengths.push(current_ungapped_len);
                update_contig_quality(&mut contig_quality, current_contig_composition);
                n_runs.finish_contig();
                total_rle_ratio += if current_len == 0 {
                    1.0
                } else {
                    current_rle_len as f64 / current_len as f64
                };
                total_rle_runs = total_rle_runs.checked_add(current_rle_len).ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("total RLE run count overflow while reading {}", path),
                    )
                })?;
            }
            in_sequence = true;
            current_len = 0;
            current_ungapped_len = 0;
            current_contig_composition = BaseComposition::default();
            current_rle_len = 0;
            current_rle_last_base = None;
        } else if in_sequence {
            let seq = trimmed_line.trim().as_bytes();
            current_len = current_len.checked_add(seq.len()).ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("contig length overflow while reading {}", path),
                )
            })?;
            let mut ungapped_bases_in_line = 0usize;
            for &base in seq {
                match base {
                    b'A' | b'a' | b'T' | b't' | b'U' | b'u' => {
                        composition.acgt_bases += 1;
                        current_contig_composition.acgt_bases += 1;
                        ungapped_bases_in_line += 1;
                    }
                    b'G' | b'g' | b'C' | b'c' => {
                        composition.gc_bases += 1;
                        composition.acgt_bases += 1;
                        current_contig_composition.gc_bases += 1;
                        current_contig_composition.acgt_bases += 1;
                        ungapped_bases_in_line += 1;
                    }
                    b'N' | b'n' => {
                        composition.n_bases += 1;
                        current_contig_composition.n_bases += 1;
                    }
                    _ => {
                        composition.ambiguous_bases += 1;
                        current_contig_composition.ambiguous_bases += 1;
                        ungapped_bases_in_line += 1;
                    }
                }
                n_runs.add_base(base);

                let run_base = normalized_run_base(base);
                if current_rle_last_base != Some(run_base) {
                    current_rle_len = current_rle_len.checked_add(1).ok_or_else(|| {
                        std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            format!("RLE run count overflow while reading {}", path),
                        )
                    })?;
                    current_rle_last_base = Some(run_base);
                }
            }
            current_ungapped_len = current_ungapped_len
                .checked_add(ungapped_bases_in_line)
                .ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("ungapped contig length overflow while reading {}", path),
                    )
                })?;
        } else if !trimmed_line.trim().is_empty() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "encountered sequence data before first FASTA header while reading {}",
                    path
                ),
            ));
        }
    }

    // Add the last sequence if there is one
    if in_sequence {
        lengths.push(current_len);
        ungapped_lengths.push(current_ungapped_len);
        update_contig_quality(&mut contig_quality, current_contig_composition);
        n_runs.finish_contig();
        total_rle_ratio += if current_len == 0 {
            1.0
        } else {
            current_rle_len as f64 / current_len as f64
        };
        total_rle_runs = total_rle_runs.checked_add(current_rle_len).ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("total RLE run count overflow while reading {}", path),
            )
        })?;
    }

    let length_stats = evaluate_lengths_in_place(&mut lengths);
    let ungapped_length_stats = evaluate_lengths_in_place(&mut ungapped_lengths);
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
    let mean_rle_ratio = if total_contigs > 0.0 {
        total_rle_ratio / total_contigs
    } else {
        0.0
    };
    let length_weighted_rle_ratio = if length_stats.total_bases > 0 {
        total_rle_runs as f64 / length_stats.total_bases as f64
    } else {
        0.0
    };
    let mean_n_run_length = if n_runs.run_count > 0 {
        composition.n_bases as f64 / n_runs.run_count as f64
    } else {
        0.0
    };

    Ok(Stats {
        total_contigs: length_stats.total,
        total_length: length_stats.total_bases,
        ungapped_total_length: ungapped_length_stats.total_bases,
        average_length: length_stats.avg_length,
        median_length: length_stats.median_length,
        n10: length_stats.n10,
        gc_bases: composition.gc_bases,
        acgt_bases: composition.acgt_bases,
        n_bases: composition.n_bases,
        ambiguous_bases: composition.ambiguous_bases,
        mean_rle_ratio,
        length_weighted_rle_ratio,
        total_rle_runs,
        n_run_count: n_runs.run_count,
        max_n_run: n_runs.max_run,
        mean_n_run_length,
        n_runs_per_100kb: per_100kb(n_runs.run_count, length_stats.total_bases),
        n_bases_per_100kb: per_100kb(composition.n_bases, length_stats.total_bases),
        ambiguous_bases_per_100kb: per_100kb(composition.ambiguous_bases, length_stats.total_bases),
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
        ungapped_n50: ungapped_length_stats.n50,
        ungapped_au_n: ungapped_length_stats.au_n,
        ungapped_effective_contig_count: ungapped_length_stats.effective_count,
        l10: length_stats.l10,
        l25: length_stats.l25,
        l50: length_stats.l50,
        l75: length_stats.l75,
        l90: length_stats.l90,
        l95: length_stats.l95,
        l99: length_stats.l99,
        au_n: length_stats.au_n,
        effective_contig_count: length_stats.effective_count,
        longest_contig: length_stats.longest,
        contigs_ge_1kb: length_stats.contigs_ge_1kb,
        contigs_ge_10kb: length_stats.contigs_ge_10kb,
        contigs_ge_50kb: length_stats.contigs_ge_50kb,
        contigs_ge_100kb: length_stats.contigs_ge_100kb,
        contigs_ge_1mb: length_stats.contigs_ge_1mb,
        bases_ge_1kb: length_stats.bases_ge_1kb,
        bases_ge_10kb: length_stats.bases_ge_10kb,
        bases_ge_50kb: length_stats.bases_ge_50kb,
        bases_ge_100kb: length_stats.bases_ge_100kb,
        bases_ge_1mb: length_stats.bases_ge_1mb,
        contigs_ge_1kb_frac: length_stats.contigs_ge_1kb_frac,
        contigs_ge_10kb_frac: length_stats.contigs_ge_10kb_frac,
        contigs_ge_50kb_frac: length_stats.contigs_ge_50kb_frac,
        contigs_ge_100kb_frac: length_stats.contigs_ge_100kb_frac,
        contigs_ge_1mb_frac: length_stats.contigs_ge_1mb_frac,
        bases_ge_1kb_frac: length_stats.bases_ge_1kb_frac,
        bases_ge_10kb_frac: length_stats.bases_ge_10kb_frac,
        bases_ge_50kb_frac: length_stats.bases_ge_50kb_frac,
        bases_ge_100kb_frac: length_stats.bases_ge_100kb_frac,
        bases_ge_1mb_frac: length_stats.bases_ge_1mb_frac,
        path_count: None,
        avg_path_length: None,
        branch_count: None,
        graph_max_depth: None,
        graph_bubble_count: None,
        graph_branchiness: None,
        graph_path_median_length: None,
        graph_path_n10: None,
        graph_path_n25: None,
        graph_path_n50: None,
        graph_path_n75: None,
        graph_path_n90: None,
        graph_path_n95: None,
        graph_path_n99: None,
        graph_path_l10: None,
        graph_path_l25: None,
        graph_path_l50: None,
        graph_path_l75: None,
        graph_path_l90: None,
        graph_path_l95: None,
        graph_path_l99: None,
        graph_path_au_n: None,
        graph_path_effective_count: None,
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
    stats.graph_path_n10 = Some(graph_stats.path_n10);
    stats.graph_path_n25 = Some(graph_stats.path_n25);
    stats.graph_path_n50 = Some(graph_stats.path_n50);
    stats.graph_path_n75 = Some(graph_stats.path_n75);
    stats.graph_path_n90 = Some(graph_stats.path_n90);
    stats.graph_path_n95 = Some(graph_stats.path_n95);
    stats.graph_path_n99 = Some(graph_stats.path_n99);
    stats.graph_path_l10 = Some(graph_stats.path_l10);
    stats.graph_path_l25 = Some(graph_stats.path_l25);
    stats.graph_path_l50 = Some(graph_stats.path_l50);
    stats.graph_path_l75 = Some(graph_stats.path_l75);
    stats.graph_path_l90 = Some(graph_stats.path_l90);
    stats.graph_path_l95 = Some(graph_stats.path_l95);
    stats.graph_path_l99 = Some(graph_stats.path_l99);
    stats.graph_path_au_n = Some(graph_stats.path_au_n);
    stats.graph_path_effective_count = Some(graph_stats.path_effective_count);
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use std::collections::HashMap;
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
        assert_eq!(stats.ungapped_total_length, 48);
        assert_eq!(stats.average_length, 16.0);
        assert_eq!(stats.median_length, 20.0);
        assert_eq!(stats.gc_bases, 24);
        assert_eq!(stats.acgt_bases, 48);
        assert_eq!(stats.n_bases, 0);
        assert_eq!(stats.ambiguous_bases, 0);
        assert_eq!(stats.mean_rle_ratio, 1.0);
        assert_eq!(stats.length_weighted_rle_ratio, 1.0);
        assert_eq!(stats.total_rle_runs, 48);
        assert_eq!(stats.n_run_count, 0);
        assert_eq!(stats.max_n_run, 0);
        assert_eq!(stats.mean_n_run_length, 0.0);
        assert_eq!(stats.n_runs_per_100kb, 0.0);
        assert_eq!(stats.n_bases_per_100kb, 0.0);
        assert_eq!(stats.ambiguous_bases_per_100kb, 0.0);
        assert_eq!(stats.contigs_with_n, 0);
        assert_eq!(stats.contigs_with_ambiguous, 0);
        assert_eq!(stats.contigs_all_acgt, 3);
        assert!((stats.gc_content - 0.5).abs() < 1e-12);
        assert_eq!(stats.n_content, 0.0);
        assert_eq!(stats.ambiguous_content, 0.0);
        assert_eq!(stats.contigs_with_n_frac, 0.0);
        assert_eq!(stats.contigs_with_ambiguous_frac, 0.0);
        assert_eq!(stats.contigs_all_acgt_frac, 1.0);
        assert_eq!(stats.n10, 24);
        assert_eq!(stats.n25, 24);
        assert_eq!(stats.n50, 24);
        assert_eq!(stats.n75, 20);
        assert_eq!(stats.n90, 20);
        assert_eq!(stats.n95, 4);
        assert_eq!(stats.n99, 4);
        assert_eq!(stats.ungapped_n50, 24);
        assert!((stats.ungapped_au_n - 20.6666666667).abs() < 1e-6);
        assert!((stats.ungapped_effective_contig_count - (48.0 / 20.6666666667)).abs() < 1e-6);
        assert_eq!(stats.l10, 1);
        assert_eq!(stats.l25, 1);
        assert_eq!(stats.l50, 1);
        assert_eq!(stats.l75, 2);
        assert_eq!(stats.l90, 2);
        assert_eq!(stats.l95, 3);
        assert_eq!(stats.l99, 3);
        assert!((stats.au_n - 20.6666666667).abs() < 1e-6);
        assert!((stats.effective_contig_count - (48.0 / 20.6666666667)).abs() < 1e-6);
        assert_eq!(stats.longest_contig, 24);
        assert_eq!(stats.contigs_ge_1kb, 0);
        assert_eq!(stats.contigs_ge_10kb, 0);
        assert_eq!(stats.contigs_ge_50kb, 0);
        assert_eq!(stats.contigs_ge_100kb, 0);
        assert_eq!(stats.contigs_ge_1mb, 0);
        assert_eq!(stats.bases_ge_1kb, 0);
        assert_eq!(stats.bases_ge_10kb, 0);
        assert_eq!(stats.bases_ge_50kb, 0);
        assert_eq!(stats.bases_ge_100kb, 0);
        assert_eq!(stats.bases_ge_1mb, 0);
        assert_eq!(stats.contigs_ge_1kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_10kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_50kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_100kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_1mb_frac, 0.0);
        assert_eq!(stats.bases_ge_1kb_frac, 0.0);
        assert_eq!(stats.bases_ge_10kb_frac, 0.0);
        assert_eq!(stats.bases_ge_50kb_frac, 0.0);
        assert_eq!(stats.bases_ge_100kb_frac, 0.0);
        assert_eq!(stats.bases_ge_1mb_frac, 0.0);
    }

    #[test]
    fn stats_tsv_row_matches_header_and_l_metrics() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">contig_1").unwrap();
        writeln!(file, "ATCGATCGATCGATCGATCG").unwrap(); // 20 bp
        writeln!(file, ">contig_2").unwrap();
        writeln!(file, "GCTAGCTAGCTAGCTAGCTAGCTA").unwrap(); // 24 bp
        writeln!(file, ">contig_3").unwrap();
        writeln!(file, "ATCG").unwrap(); // 4 bp

        let stats = calculate_stats(file.path().to_str().unwrap()).unwrap();
        let header_cols: Vec<&str> = tsv_header().split('\t').collect();
        let row = tsv_row(&stats);
        let row_cols: Vec<&str> = row.split('\t').collect();

        assert_eq!(header_cols.len(), row_cols.len());

        let row_by_name: HashMap<&str, &str> = header_cols
            .iter()
            .copied()
            .zip(row_cols.iter().copied())
            .collect();

        assert_eq!(row_by_name.get("l10").copied(), Some("1"));
        assert_eq!(row_by_name.get("l25").copied(), Some("1"));
        assert_eq!(row_by_name.get("l50").copied(), Some("1"));
        assert_eq!(row_by_name.get("l75").copied(), Some("2"));
        assert_eq!(row_by_name.get("l90").copied(), Some("2"));
        assert_eq!(row_by_name.get("l95").copied(), Some("3"));
        assert_eq!(row_by_name.get("l99").copied(), Some("3"));
        assert_eq!(row_by_name.get("ungapped_total_len").copied(), Some("48"));
        assert_eq!(row_by_name.get("ungapped_n50").copied(), Some("24"));
        assert_eq!(row_by_name.get("ungapped_aun").copied(), Some("20.67"));
    }

    #[test]
    fn stats_tsv_row_reports_ungapped_metrics_with_gaps() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">contig_1").unwrap();
        writeln!(file, "AAANNN").unwrap(); // len 6, ungapped len 3
        writeln!(file, ">contig_2").unwrap();
        writeln!(file, "GGGG").unwrap(); // len 4, ungapped len 4

        let stats = calculate_stats(file.path().to_str().unwrap()).unwrap();
        let header_cols: Vec<&str> = tsv_header().split('\t').collect();
        let row = tsv_row(&stats);
        let row_cols: Vec<&str> = row.split('\t').collect();
        let row_by_name: HashMap<&str, &str> = header_cols
            .iter()
            .copied()
            .zip(row_cols.iter().copied())
            .collect();

        assert_eq!(row_by_name.get("total_len").copied(), Some("10"));
        assert_eq!(row_by_name.get("ungapped_total_len").copied(), Some("7"));
        assert_eq!(row_by_name.get("n50").copied(), Some("6"));
        assert_eq!(row_by_name.get("ungapped_n50").copied(), Some("4"));
        assert_eq!(row_by_name.get("ungapped_aun").copied(), Some("3.57"));
    }

    #[test]
    fn stats_tsv_row_reports_n_run_density_metrics() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">contig_1").unwrap();
        writeln!(file, "AANN").unwrap();
        writeln!(file, ">contig_2").unwrap();
        writeln!(file, "NNRA").unwrap();

        let stats = calculate_stats(file.path().to_str().unwrap()).unwrap();
        let header_cols: Vec<&str> = tsv_header().split('\t').collect();
        let row = tsv_row(&stats);
        let row_cols: Vec<&str> = row.split('\t').collect();
        let row_by_name: HashMap<&str, &str> = header_cols
            .iter()
            .copied()
            .zip(row_cols.iter().copied())
            .collect();

        assert_eq!(row_by_name.get("n_run_count").copied(), Some("2"));
        assert_eq!(row_by_name.get("max_n_run").copied(), Some("2"));
        assert_eq!(
            row_by_name.get("mean_n_run_length").copied(),
            Some("2.000000")
        );
        assert_eq!(
            row_by_name.get("n_runs_per_100kb").copied(),
            Some("25000.000000")
        );
        assert_eq!(
            row_by_name.get("n_bases_per_100kb").copied(),
            Some("50000.000000")
        );
        assert_eq!(
            row_by_name.get("ambiguous_bases_per_100kb").copied(),
            Some("12500.000000")
        );
    }

    proptest! {
        #[test]
        fn stats_tsv_row_column_count_stays_in_sync_for_random_inputs(
            seqs in prop::collection::vec("[ACGTNacgtnRYSWKMBDHVryswkmbdhv]{0,64}", 1..32)
        ) {
            let mut file = NamedTempFile::new().expect("create temp fasta");
            for (idx, seq) in seqs.iter().enumerate() {
                writeln!(file, ">contig_{idx}").expect("write header");
                writeln!(file, "{seq}").expect("write sequence");
            }

            let stats = calculate_stats(file.path().to_str().expect("utf8 path"))
                .expect("calculate stats");

            let header_cols = tsv_header().split('\t').count();
            let row_cols = tsv_row(&stats).split('\t').count();
            prop_assert_eq!(row_cols, header_cols);
        }
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
        assert_eq!(stats.ungapped_total_length, 16);
        assert_eq!(stats.median_length, 8.0);
        assert_eq!(stats.gc_bases, 8);
        assert_eq!(stats.acgt_bases, 16);
        assert_eq!(stats.n_bases, 0);
        assert_eq!(stats.ambiguous_bases, 0);
        assert_eq!(stats.mean_rle_ratio, 1.0);
        assert_eq!(stats.length_weighted_rle_ratio, 1.0);
        assert_eq!(stats.total_rle_runs, 16);
        assert_eq!(stats.n_run_count, 0);
        assert_eq!(stats.max_n_run, 0);
        assert_eq!(stats.mean_n_run_length, 0.0);
        assert_eq!(stats.n_runs_per_100kb, 0.0);
        assert_eq!(stats.n_bases_per_100kb, 0.0);
        assert_eq!(stats.ambiguous_bases_per_100kb, 0.0);
        assert_eq!(stats.contigs_with_n, 0);
        assert_eq!(stats.contigs_with_ambiguous, 0);
        assert_eq!(stats.contigs_all_acgt, 2);
        assert!((stats.gc_content - 0.5).abs() < 1e-12);
        assert_eq!(stats.n_content, 0.0);
        assert_eq!(stats.ambiguous_content, 0.0);
        assert_eq!(stats.contigs_with_n_frac, 0.0);
        assert_eq!(stats.contigs_with_ambiguous_frac, 0.0);
        assert_eq!(stats.contigs_all_acgt_frac, 1.0);
        assert_eq!(stats.n10, 12);
        assert_eq!(stats.n25, 12);
        assert_eq!(stats.n50, 12);
        assert_eq!(stats.n75, 12);
        assert_eq!(stats.n90, 4);
        assert_eq!(stats.n95, 4);
        assert_eq!(stats.n99, 4);
        assert_eq!(stats.ungapped_n50, 12);
        assert!((stats.ungapped_au_n - 10.0).abs() < 1e-6);
        assert!((stats.ungapped_effective_contig_count - 1.6).abs() < 1e-12);
        assert_eq!(stats.l10, 1);
        assert_eq!(stats.l25, 1);
        assert_eq!(stats.l50, 1);
        assert_eq!(stats.l75, 1);
        assert_eq!(stats.l90, 2);
        assert_eq!(stats.l95, 2);
        assert_eq!(stats.l99, 2);
        assert!((stats.au_n - 10.0).abs() < 1e-6);
        assert!((stats.effective_contig_count - 1.6).abs() < 1e-12);
        assert_eq!(stats.longest_contig, 12);
        assert_eq!(stats.contigs_ge_1kb, 0);
        assert_eq!(stats.contigs_ge_10kb, 0);
        assert_eq!(stats.contigs_ge_50kb, 0);
        assert_eq!(stats.contigs_ge_100kb, 0);
        assert_eq!(stats.contigs_ge_1mb, 0);
        assert_eq!(stats.bases_ge_1kb, 0);
        assert_eq!(stats.bases_ge_10kb, 0);
        assert_eq!(stats.bases_ge_50kb, 0);
        assert_eq!(stats.bases_ge_100kb, 0);
        assert_eq!(stats.bases_ge_1mb, 0);
        assert_eq!(stats.contigs_ge_1kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_10kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_50kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_100kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_1mb_frac, 0.0);
        assert_eq!(stats.bases_ge_1kb_frac, 0.0);
        assert_eq!(stats.bases_ge_10kb_frac, 0.0);
        assert_eq!(stats.bases_ge_50kb_frac, 0.0);
        assert_eq!(stats.bases_ge_100kb_frac, 0.0);
        assert_eq!(stats.bases_ge_1mb_frac, 0.0);
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
        assert_eq!(stats.contigs_ge_1mb, 0);
        assert_eq!(stats.bases_ge_1kb, 61_000);
        assert_eq!(stats.bases_ge_10kb, 60_000);
        assert_eq!(stats.bases_ge_50kb, 50_000);
        assert_eq!(stats.bases_ge_100kb, 0);
        assert_eq!(stats.bases_ge_1mb, 0);
        assert!((stats.contigs_ge_1kb_frac - 0.75).abs() < 1e-12);
        assert!((stats.contigs_ge_10kb_frac - 0.5).abs() < 1e-12);
        assert!((stats.contigs_ge_50kb_frac - 0.25).abs() < 1e-12);
        assert_eq!(stats.contigs_ge_100kb_frac, 0.0);
        assert_eq!(stats.contigs_ge_1mb_frac, 0.0);
        assert!((stats.bases_ge_1kb_frac - (61_000.0 / 61_999.0)).abs() < 1e-12);
        assert!((stats.bases_ge_10kb_frac - (60_000.0 / 61_999.0)).abs() < 1e-12);
        assert!((stats.bases_ge_50kb_frac - (50_000.0 / 61_999.0)).abs() < 1e-12);
        assert_eq!(stats.bases_ge_100kb_frac, 0.0);
        assert_eq!(stats.bases_ge_1mb_frac, 0.0);
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
        assert_eq!(stats.contigs_ge_1mb, 0);
        assert_eq!(stats.bases_ge_1mb, 0);
        assert_eq!(stats.contigs_ge_1mb_frac, 0.0);
        assert_eq!(stats.bases_ge_1mb_frac, 0.0);
    }

    #[test]
    fn test_calculate_stats_reports_1mb_bucket_metrics() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">contig_1").unwrap();
        writeln!(file, "{}", "A".repeat(1_500_000)).unwrap();
        writeln!(file, ">contig_2").unwrap();
        writeln!(file, "{}", "C".repeat(900_000)).unwrap();
        writeln!(file, ">contig_3").unwrap();
        writeln!(file, "{}", "G".repeat(100_000)).unwrap();

        let stats = calculate_stats(file.path().to_str().unwrap()).unwrap();
        assert_eq!(stats.contigs_ge_1mb, 1);
        assert_eq!(stats.bases_ge_1mb, 1_500_000);
        assert!((stats.contigs_ge_1mb_frac - (1.0 / 3.0)).abs() < 1e-12);
        assert!((stats.bases_ge_1mb_frac - (1_500_000.0 / 2_500_000.0)).abs() < 1e-12);
    }

    #[test]
    fn test_update_with_graph_stats_populates_extended_metrics() {
        let mut stats = Stats {
            total_contigs: 0,
            total_length: 0,
            ungapped_total_length: 0,
            average_length: 0.0,
            median_length: 0.0,
            n10: 0,
            gc_bases: 0,
            acgt_bases: 0,
            n_bases: 0,
            ambiguous_bases: 0,
            mean_rle_ratio: 0.0,
            length_weighted_rle_ratio: 0.0,
            total_rle_runs: 0,
            n_run_count: 0,
            max_n_run: 0,
            mean_n_run_length: 0.0,
            n_runs_per_100kb: 0.0,
            n_bases_per_100kb: 0.0,
            ambiguous_bases_per_100kb: 0.0,
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
            ungapped_n50: 0,
            ungapped_au_n: 0.0,
            ungapped_effective_contig_count: 0.0,
            l10: 0,
            l25: 0,
            l50: 0,
            l75: 0,
            l90: 0,
            l95: 0,
            l99: 0,
            au_n: 0.0,
            effective_contig_count: 0.0,
            longest_contig: 0,
            contigs_ge_1kb: 0,
            contigs_ge_10kb: 0,
            contigs_ge_50kb: 0,
            contigs_ge_100kb: 0,
            contigs_ge_1mb: 0,
            bases_ge_1kb: 0,
            bases_ge_10kb: 0,
            bases_ge_50kb: 0,
            bases_ge_100kb: 0,
            bases_ge_1mb: 0,
            contigs_ge_1kb_frac: 0.0,
            contigs_ge_10kb_frac: 0.0,
            contigs_ge_50kb_frac: 0.0,
            contigs_ge_100kb_frac: 0.0,
            contigs_ge_1mb_frac: 0.0,
            bases_ge_1kb_frac: 0.0,
            bases_ge_10kb_frac: 0.0,
            bases_ge_50kb_frac: 0.0,
            bases_ge_100kb_frac: 0.0,
            bases_ge_1mb_frac: 0.0,
            path_count: None,
            avg_path_length: None,
            branch_count: None,
            graph_max_depth: None,
            graph_bubble_count: None,
            graph_branchiness: None,
            graph_path_median_length: None,
            graph_path_n10: None,
            graph_path_n25: None,
            graph_path_n50: None,
            graph_path_n75: None,
            graph_path_n90: None,
            graph_path_n95: None,
            graph_path_n99: None,
            graph_path_l10: None,
            graph_path_l25: None,
            graph_path_l50: None,
            graph_path_l75: None,
            graph_path_l90: None,
            graph_path_l95: None,
            graph_path_l99: None,
            graph_path_au_n: None,
            graph_path_effective_count: None,
        };
        let graph_stats = crate::graph::complexity::PathStats {
            total_paths: 7,
            average_length: 3.5,
            median_length: 3.0,
            path_n10: 5,
            path_n25: 4,
            path_n50: 4,
            path_n75: 3,
            path_n90: 2,
            path_n95: 2,
            path_n99: 2,
            path_l10: 1,
            path_l25: 2,
            path_l50: 3,
            path_l75: 4,
            path_l90: 6,
            path_l95: 7,
            path_l99: 7,
            path_au_n: 3.2,
            path_effective_count: 2.1875,
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
        assert_eq!(stats.graph_path_n10, Some(5));
        assert_eq!(stats.graph_path_n25, Some(4));
        assert_eq!(stats.graph_path_n50, Some(4));
        assert_eq!(stats.graph_path_n75, Some(3));
        assert_eq!(stats.graph_path_n90, Some(2));
        assert_eq!(stats.graph_path_n95, Some(2));
        assert_eq!(stats.graph_path_n99, Some(2));
        assert_eq!(stats.graph_path_l10, Some(1));
        assert_eq!(stats.graph_path_l25, Some(2));
        assert_eq!(stats.graph_path_l50, Some(3));
        assert_eq!(stats.graph_path_l75, Some(4));
        assert_eq!(stats.graph_path_l90, Some(6));
        assert_eq!(stats.graph_path_l95, Some(7));
        assert_eq!(stats.graph_path_l99, Some(7));
        assert_eq!(stats.graph_path_au_n, Some(3.2));
        assert_eq!(stats.graph_path_effective_count, Some(2.1875));
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
        assert_eq!(stats.ungapped_total_length, 8);
        assert_eq!(stats.gc_bases, 4);
        assert_eq!(stats.acgt_bases, 7);
        assert_eq!(stats.n_bases, 4);
        assert_eq!(stats.ambiguous_bases, 1);
        assert!((stats.mean_rle_ratio - 0.8125).abs() < 1e-12);
        assert!((stats.length_weighted_rle_ratio - 0.75).abs() < 1e-12);
        assert_eq!(stats.total_rle_runs, 9);
        assert_eq!(stats.n_run_count, 1);
        assert_eq!(stats.max_n_run, 4);
        assert!((stats.mean_n_run_length - 4.0).abs() < 1e-12);
        assert!((stats.n_runs_per_100kb - (1.0 * 100_000.0 / 12.0)).abs() < 1e-12);
        assert!((stats.n_bases_per_100kb - (4.0 * 100_000.0 / 12.0)).abs() < 1e-12);
        assert!((stats.ambiguous_bases_per_100kb - (1.0 * 100_000.0 / 12.0)).abs() < 1e-12);
        assert_eq!(stats.contigs_with_n, 1);
        assert_eq!(stats.contigs_with_ambiguous, 1);
        assert_eq!(stats.contigs_all_acgt, 0);
        assert_eq!(stats.ungapped_n50, 4);
        assert!((stats.ungapped_au_n - 4.0).abs() < 1e-12);
        assert!((stats.ungapped_effective_contig_count - 2.0).abs() < 1e-12);
        assert!((stats.gc_content - (4.0 / 7.0)).abs() < 1e-12);
        assert!((stats.n_content - (4.0 / 12.0)).abs() < 1e-12);
        assert!((stats.ambiguous_content - (1.0 / 12.0)).abs() < 1e-12);
        assert!((stats.contigs_with_n_frac - 0.5).abs() < 1e-12);
        assert!((stats.contigs_with_ambiguous_frac - 0.5).abs() < 1e-12);
        assert_eq!(stats.contigs_all_acgt_frac, 0.0);
        assert!((stats.effective_contig_count - (12.0 / 6.6666666667)).abs() < 1e-6);
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
        assert_eq!(stats.ungapped_total_length, 3);
        assert!((stats.average_length - 0.75).abs() < 1e-12);
        assert_eq!(stats.median_length, 0.5);
        assert_eq!(stats.n10, 2);
        assert_eq!(stats.n25, 2);
        assert_eq!(stats.n50, 2);
        assert_eq!(stats.n75, 1);
        assert_eq!(stats.n90, 1);
        assert_eq!(stats.n95, 1);
        assert_eq!(stats.n99, 1);
        assert_eq!(stats.ungapped_n50, 2);
        assert_eq!(stats.l10, 1);
        assert_eq!(stats.l25, 1);
        assert_eq!(stats.l50, 1);
        assert_eq!(stats.l75, 2);
        assert_eq!(stats.l90, 2);
        assert_eq!(stats.l95, 2);
        assert_eq!(stats.l99, 2);
        assert_eq!(stats.contigs_with_n, 0);
        assert_eq!(stats.contigs_with_ambiguous, 0);
        assert_eq!(stats.contigs_all_acgt, 4);
        assert_eq!(stats.mean_rle_ratio, 1.0);
        assert_eq!(stats.length_weighted_rle_ratio, 1.0);
        assert_eq!(stats.total_rle_runs, 3);
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
    fn test_calculate_stats_reports_n_run_metrics_across_wrapped_lines_and_contigs() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">contig_1").unwrap();
        writeln!(file, "AANN").unwrap();
        writeln!(file, "NNCC").unwrap(); // same run continues across line boundary (run len 4)
        writeln!(file, ">contig_2").unwrap();
        writeln!(file, "NN").unwrap(); // second run (len 2)
        writeln!(file, "A").unwrap();
        writeln!(file, "NNN").unwrap(); // third run (len 3)
        writeln!(file, ">contig_3").unwrap();
        writeln!(file, "ACGT").unwrap();

        let stats = calculate_stats(file.path().to_str().unwrap()).unwrap();
        assert_eq!(stats.n_run_count, 3);
        assert_eq!(stats.max_n_run, 4);
        assert!((stats.mean_n_run_length - 3.0).abs() < 1e-12);
        assert!((stats.n_runs_per_100kb - (3.0 * 100_000.0 / 18.0)).abs() < 1e-12);
        assert!((stats.n_bases_per_100kb - (9.0 * 100_000.0 / 18.0)).abs() < 1e-12);
    }

    #[test]
    fn test_calculate_stats_reports_mean_rle_ratio_across_wrapped_lines() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">contig_1").unwrap();
        writeln!(file, "AAA").unwrap();
        writeln!(file, "AAA").unwrap(); // one run across wrapped lines => 1/6
        writeln!(file, ">contig_2").unwrap();
        writeln!(file, "ATAT").unwrap(); // four runs => 4/4

        let stats = calculate_stats(file.path().to_str().unwrap()).unwrap();
        assert!((stats.mean_rle_ratio - ((1.0 / 6.0 + 1.0) / 2.0)).abs() < 1e-12);
        assert!((stats.length_weighted_rle_ratio - 0.5).abs() < 1e-12);
        assert_eq!(stats.total_rle_runs, 5);
    }

    #[test]
    fn test_calculate_stats_rle_metrics_are_case_insensitive() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">contig_1").unwrap();
        writeln!(file, "AaAA").unwrap();
        writeln!(file, ">contig_2").unwrap();
        writeln!(file, "tT").unwrap();

        let stats = calculate_stats(file.path().to_str().unwrap()).unwrap();
        assert_eq!(stats.total_rle_runs, 2);
        assert!((stats.mean_rle_ratio - ((1.0 / 4.0 + 1.0 / 2.0) / 2.0)).abs() < 1e-12);
        assert!((stats.length_weighted_rle_ratio - (2.0 / 6.0)).abs() < 1e-12);
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

    #[test]
    fn test_calculate_stats_rejects_sequence_data_before_first_header() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "ACGT").unwrap();
        writeln!(file, ">contig_1").unwrap();
        writeln!(file, "TTTT").unwrap();

        match calculate_stats(file.path().to_str().unwrap()) {
            Ok(_) => panic!("expected invalid FASTA error"),
            Err(err) => {
                assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
                assert!(err
                    .to_string()
                    .contains("sequence data before first FASTA header"));
            }
        }
    }

    #[test]
    fn test_calculate_stats_allows_leading_blank_lines_before_first_header() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file).unwrap();
        writeln!(file, "   ").unwrap();
        writeln!(file, ">contig_1").unwrap();
        writeln!(file, "ACGT").unwrap();

        let stats = calculate_stats(file.path().to_str().unwrap()).unwrap();
        assert_eq!(stats.total_contigs, 1);
        assert_eq!(stats.total_length, 4);
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
            prop_assert_eq!(stats_a.ungapped_total_length, stats_b.ungapped_total_length);
            prop_assert_eq!(stats_a.gc_bases, stats_b.gc_bases);
            prop_assert_eq!(stats_a.acgt_bases, stats_b.acgt_bases);
            prop_assert_eq!(stats_a.n_bases, stats_b.n_bases);
            prop_assert_eq!(stats_a.ambiguous_bases, stats_b.ambiguous_bases);
            prop_assert!((stats_a.mean_rle_ratio - stats_b.mean_rle_ratio).abs() < 1e-12);
            prop_assert!((stats_a.length_weighted_rle_ratio - stats_b.length_weighted_rle_ratio).abs() < 1e-12);
            prop_assert_eq!(stats_a.total_rle_runs, stats_b.total_rle_runs);
            prop_assert_eq!(stats_a.n_run_count, stats_b.n_run_count);
            prop_assert_eq!(stats_a.max_n_run, stats_b.max_n_run);
            prop_assert!((stats_a.mean_n_run_length - stats_b.mean_n_run_length).abs() < 1e-12);
            prop_assert!((stats_a.n_runs_per_100kb - stats_b.n_runs_per_100kb).abs() < 1e-12);
            prop_assert!((stats_a.n_bases_per_100kb - stats_b.n_bases_per_100kb).abs() < 1e-12);
            prop_assert!((stats_a.ambiguous_bases_per_100kb - stats_b.ambiguous_bases_per_100kb).abs() < 1e-12);
            prop_assert_eq!(stats_a.contigs_with_n, stats_b.contigs_with_n);
            prop_assert_eq!(stats_a.contigs_with_ambiguous, stats_b.contigs_with_ambiguous);
            prop_assert_eq!(stats_a.contigs_all_acgt, stats_b.contigs_all_acgt);
            prop_assert_eq!(stats_a.n10, stats_b.n10);
            prop_assert_eq!(stats_a.n25, stats_b.n25);
            prop_assert_eq!(stats_a.n50, stats_b.n50);
            prop_assert_eq!(stats_a.n75, stats_b.n75);
            prop_assert_eq!(stats_a.n90, stats_b.n90);
            prop_assert_eq!(stats_a.n95, stats_b.n95);
            prop_assert_eq!(stats_a.n99, stats_b.n99);
            prop_assert_eq!(stats_a.ungapped_n50, stats_b.ungapped_n50);
            prop_assert!((stats_a.ungapped_au_n - stats_b.ungapped_au_n).abs() < 1e-12);
            prop_assert!(
                (stats_a.ungapped_effective_contig_count - stats_b.ungapped_effective_contig_count)
                    .abs()
                    < 1e-12
            );
            prop_assert_eq!(stats_a.l10, stats_b.l10);
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
            prop_assert!(
                (stats_a.effective_contig_count - stats_b.effective_contig_count).abs() < 1e-12
            );
        }
    }
}
