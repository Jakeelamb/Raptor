use crate::accel::CpuBackend;
use crate::eval::metrics::{evaluate_lengths_in_place, normalized_run_base, BaseComposition};
use crate::graph::assembler::{greedy_assembly_u64, Contig};
use crate::graph::overlap::find_overlaps;
use crate::graph::stitch::OverlapGraphBuilder;
use crate::io::fasta::FastaWriter;
use crate::io::fastq::{
    stream_fastq_records_checked, stream_paired_fastq_records_checked, try_open_fastq, FastqRecord,
};
use crate::io::gfa::GfaWriter;
use crate::io::gfa2::Gfa2Writer;
use crate::kmer::variable_k::{kmer_coverage_histogram, optimal_k, select_best_k};
use ahash::AHashMap;
use serde::Serialize;
use std::fs;
use std::io::{self, Write};
use tracing::{info, warn};

const READ_OVERLAP_RESCUE_MAX_READS: usize = 200_000;

#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
struct AssemblyQualitySummary {
    total_contigs: usize,
    total_bases: usize,
    ungapped_total_bases: usize,
    avg_length: f64,
    median_length: f64,
    n10: usize,
    n25: usize,
    n50: usize,
    n75: usize,
    n90: usize,
    n95: usize,
    n99: usize,
    l10: usize,
    l25: usize,
    l50: usize,
    l75: usize,
    l90: usize,
    l95: usize,
    l99: usize,
    au_n: f64,
    effective_contig_count: f64,
    ungapped_n50: usize,
    ungapped_n90: usize,
    ungapped_n95: usize,
    ungapped_n99: usize,
    ungapped_au_n: f64,
    ungapped_effective_contig_count: f64,
    gap_bases: usize,
    gap_bases_frac: f64,
    longest: usize,
    longest_frac: f64,
    gc_bases: usize,
    acgt_bases: usize,
    n_bases: usize,
    ambiguous_bases: usize,
    contigs_with_n: usize,
    contigs_with_ambiguous: usize,
    contigs_all_acgt: usize,
    distinct_sequences: usize,
    duplicate_sequences: usize,
    duplicate_bases: usize,
    canonical_distinct_sequences: usize,
    canonical_duplicate_sequences: usize,
    canonical_duplicate_bases: usize,
    contigs_with_n_frac: f64,
    contigs_with_ambiguous_frac: f64,
    contigs_all_acgt_frac: f64,
    distinct_sequences_frac: f64,
    duplicate_sequences_frac: f64,
    duplicate_bases_frac: f64,
    canonical_distinct_sequences_frac: f64,
    canonical_duplicate_sequences_frac: f64,
    canonical_duplicate_bases_frac: f64,
    mean_rle_ratio: f64,
    length_weighted_rle_ratio: f64,
    total_rle_runs: usize,
    gc_content: f64,
    n_content: f64,
    ambiguous_content: f64,
    n_runs: usize,
    longest_n_run: usize,
    mean_n_run_length: f64,
    n_runs_per_100kb: f64,
    n_bases_per_100kb: f64,
    ambiguous_bases_per_100kb: f64,
    contigs_ge_1kb: usize,
    contigs_ge_10kb: usize,
    contigs_ge_50kb: usize,
    contigs_ge_100kb: usize,
    contigs_ge_1mb: usize,
    bases_ge_1kb: usize,
    bases_ge_10kb: usize,
    bases_ge_50kb: usize,
    bases_ge_100kb: usize,
    bases_ge_1mb: usize,
    contigs_ge_1kb_frac: f64,
    contigs_ge_10kb_frac: f64,
    contigs_ge_50kb_frac: f64,
    contigs_ge_100kb_frac: f64,
    contigs_ge_1mb_frac: f64,
    bases_ge_1kb_frac: f64,
    bases_ge_10kb_frac: f64,
    bases_ge_50kb_frac: f64,
    bases_ge_100kb_frac: f64,
    bases_ge_1mb_frac: f64,
}

#[derive(Debug, Clone, Copy, Default)]
struct NRunSummary {
    count: usize,
    longest: usize,
}

#[derive(Debug, Clone, Copy, Default)]
struct ContigQualitySummary {
    with_n: usize,
    with_ambiguous: usize,
    all_acgt: usize,
}

#[derive(Debug, Clone, Copy, Default)]
struct SequenceAnalysis {
    ungapped_len: usize,
    rle_runs: usize,
    has_n: bool,
    has_ambiguous: bool,
}

#[derive(Debug, Clone, Copy, Default)]
struct DuplicateSummary {
    distinct_sequences: usize,
    duplicate_sequences: usize,
    duplicate_bases: usize,
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
fn analyze_sequence(
    sequence: &[u8],
    composition: &mut BaseComposition,
    n_runs: &mut NRunSummary,
) -> SequenceAnalysis {
    let mut ungapped_len = sequence.len();
    let mut rle_runs = 0usize;
    let mut previous_run_base = None;
    let mut run_len = 0usize;
    let mut has_n = false;
    let mut has_ambiguous = false;

    for &base in sequence {
        let run_base = normalized_run_base(base);
        if previous_run_base != Some(run_base) {
            previous_run_base = Some(run_base);
            rle_runs = rle_runs.saturating_add(1);
        }

        match base {
            b'A' | b'a' | b'T' | b't' | b'U' | b'u' => {
                composition.acgt_bases += 1;
                if run_len > 0 {
                    n_runs.count += 1;
                    n_runs.longest = n_runs.longest.max(run_len);
                    run_len = 0;
                }
            }
            b'G' | b'g' | b'C' | b'c' => {
                composition.gc_bases += 1;
                composition.acgt_bases += 1;
                if run_len > 0 {
                    n_runs.count += 1;
                    n_runs.longest = n_runs.longest.max(run_len);
                    run_len = 0;
                }
            }
            b'N' | b'n' => {
                composition.n_bases += 1;
                ungapped_len = ungapped_len.saturating_sub(1);
                run_len += 1;
                has_n = true;
            }
            _ => {
                composition.ambiguous_bases += 1;
                has_ambiguous = true;
                if run_len > 0 {
                    n_runs.count += 1;
                    n_runs.longest = n_runs.longest.max(run_len);
                    run_len = 0;
                }
            }
        }
    }

    if run_len > 0 {
        n_runs.count += 1;
        n_runs.longest = n_runs.longest.max(run_len);
    }

    SequenceAnalysis {
        ungapped_len,
        rle_runs,
        has_n,
        has_ambiguous,
    }
}

#[inline]
fn normalize_canonical_base(base: u8) -> u8 {
    match base {
        b'U' | b'u' => b'T',
        _ => base.to_ascii_uppercase(),
    }
}

#[inline]
fn complement_iupac_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'N' => b'N',
        _ => base,
    }
}

fn reverse_complement_dna(sequence: &str) -> String {
    let mut reverse = String::with_capacity(sequence.len());
    for base in sequence.bytes().rev() {
        reverse.push(match base {
            b'A' | b'a' => 'T',
            b'C' | b'c' => 'G',
            b'G' | b'g' => 'C',
            b'T' | b't' | b'U' | b'u' => 'A',
            b'R' | b'r' => 'Y',
            b'Y' | b'y' => 'R',
            b'S' | b's' => 'S',
            b'W' | b'w' => 'W',
            b'K' | b'k' => 'M',
            b'M' | b'm' => 'K',
            b'B' | b'b' => 'V',
            b'V' | b'v' => 'B',
            b'D' | b'd' => 'H',
            b'H' | b'h' => 'D',
            b'N' | b'n' => 'N',
            _ => 'N',
        });
    }
    reverse
}

#[inline]
fn canonical_sequence_key(sequence: &str) -> Vec<u8> {
    let bytes = sequence.as_bytes();
    let mut forward = Vec::with_capacity(bytes.len());
    let mut reverse_complement = Vec::with_capacity(bytes.len());

    for &base in bytes {
        forward.push(normalize_canonical_base(base));
    }

    for &base in bytes.iter().rev() {
        let normalized = normalize_canonical_base(base);
        reverse_complement.push(complement_iupac_base(normalized));
    }

    if reverse_complement < forward {
        reverse_complement
    } else {
        forward
    }
}

#[inline]
fn summarize_duplicate_counts<K, F>(
    counts: &std::collections::HashMap<K, usize>,
    mut key_len: F,
) -> DuplicateSummary
where
    K: Eq + std::hash::Hash,
    F: FnMut(&K) -> usize,
{
    let mut summary = DuplicateSummary::default();
    for (key, count) in counts {
        summary.distinct_sequences += 1;
        let duplicates = count.saturating_sub(1);
        summary.duplicate_sequences = summary.duplicate_sequences.saturating_add(duplicates);
        summary.duplicate_bases = summary
            .duplicate_bases
            .saturating_add(key_len(key).saturating_mul(duplicates));
    }
    summary
}

#[inline]
fn summarize_assembly_quality(contigs: &[Contig]) -> AssemblyQualitySummary {
    let mut lengths = Vec::with_capacity(contigs.len());
    let mut ungapped_lengths = Vec::with_capacity(contigs.len());
    let mut composition = BaseComposition::default();
    let mut n_runs = NRunSummary::default();
    let mut contig_quality = ContigQualitySummary::default();
    let mut sequence_counts: std::collections::HashMap<&str, usize> =
        std::collections::HashMap::with_capacity(contigs.len());
    let mut canonical_sequence_counts: std::collections::HashMap<Vec<u8>, usize> =
        std::collections::HashMap::with_capacity(contigs.len());
    let mut total_rle_ratio_scaled = 0u128;
    let mut total_rle_runs = 0usize;

    for contig in contigs {
        let contig_len = contig.sequence.len();
        *sequence_counts.entry(contig.sequence.as_str()).or_insert(0) += 1;
        *canonical_sequence_counts
            .entry(canonical_sequence_key(&contig.sequence))
            .or_insert(0) += 1;
        lengths.push(contig_len);
        let sequence = contig.sequence.as_bytes();
        let analysis = analyze_sequence(sequence, &mut composition, &mut n_runs);
        ungapped_lengths.push(analysis.ungapped_len);
        total_rle_runs = total_rle_runs.saturating_add(analysis.rle_runs);
        if analysis.has_n {
            contig_quality.with_n += 1;
        }
        if analysis.has_ambiguous {
            contig_quality.with_ambiguous += 1;
        }
        if !analysis.has_n && !analysis.has_ambiguous {
            contig_quality.all_acgt += 1;
        }
        if contig_len == 0 {
            total_rle_ratio_scaled = total_rle_ratio_scaled.saturating_add(RLE_RATIO_SCALE);
        } else {
            let scaled = ((analysis.rle_runs as u128)
                .saturating_mul(RLE_RATIO_SCALE)
                .saturating_add((contig_len as u128) / 2))
                / contig_len as u128;
            total_rle_ratio_scaled = total_rle_ratio_scaled.saturating_add(scaled);
        }
    }

    let length_stats = evaluate_lengths_in_place(&mut lengths);
    let ungapped_length_stats = evaluate_lengths_in_place(&mut ungapped_lengths);
    let longest_frac = if length_stats.total_bases > 0 {
        length_stats.longest as f64 / length_stats.total_bases as f64
    } else {
        0.0
    };
    let gap_bases = length_stats
        .total_bases
        .saturating_sub(ungapped_length_stats.total_bases);
    let gap_bases_frac = if length_stats.total_bases > 0 {
        gap_bases as f64 / length_stats.total_bases as f64
    } else {
        0.0
    };
    let mean_n_run_length = if n_runs.count > 0 {
        composition.n_bases as f64 / n_runs.count as f64
    } else {
        0.0
    };
    let mean_rle_ratio = if length_stats.total > 0 {
        total_rle_ratio_scaled as f64 / (length_stats.total as f64 * RLE_RATIO_SCALE as f64)
    } else {
        0.0
    };
    let length_weighted_rle_ratio = if length_stats.total_bases > 0 {
        total_rle_runs as f64 / length_stats.total_bases as f64
    } else {
        0.0
    };
    let duplicate_summary = summarize_duplicate_counts(&sequence_counts, |sequence| sequence.len());
    let canonical_duplicate_summary =
        summarize_duplicate_counts(&canonical_sequence_counts, |sequence| sequence.len());
    let total_contigs_f = length_stats.total as f64;
    let (contigs_with_n_frac, contigs_with_ambiguous_frac, contigs_all_acgt_frac) =
        if total_contigs_f > 0.0 {
            (
                contig_quality.with_n as f64 / total_contigs_f,
                contig_quality.with_ambiguous as f64 / total_contigs_f,
                contig_quality.all_acgt as f64 / total_contigs_f,
            )
        } else {
            (0.0, 0.0, 0.0)
        };
    let (distinct_sequences_frac, duplicate_sequences_frac) = if total_contigs_f > 0.0 {
        (
            duplicate_summary.distinct_sequences as f64 / total_contigs_f,
            duplicate_summary.duplicate_sequences as f64 / total_contigs_f,
        )
    } else {
        (0.0, 0.0)
    };
    let (canonical_distinct_sequences_frac, canonical_duplicate_sequences_frac) =
        if total_contigs_f > 0.0 {
            (
                canonical_duplicate_summary.distinct_sequences as f64 / total_contigs_f,
                canonical_duplicate_summary.duplicate_sequences as f64 / total_contigs_f,
            )
        } else {
            (0.0, 0.0)
        };
    let duplicate_bases_frac = if length_stats.total_bases > 0 {
        duplicate_summary.duplicate_bases as f64 / length_stats.total_bases as f64
    } else {
        0.0
    };
    let canonical_duplicate_bases_frac = if length_stats.total_bases > 0 {
        canonical_duplicate_summary.duplicate_bases as f64 / length_stats.total_bases as f64
    } else {
        0.0
    };
    AssemblyQualitySummary {
        total_contigs: length_stats.total,
        total_bases: length_stats.total_bases,
        ungapped_total_bases: ungapped_length_stats.total_bases,
        avg_length: length_stats.avg_length,
        median_length: length_stats.median_length,
        n10: length_stats.n10,
        n25: length_stats.n25,
        n50: length_stats.n50,
        n75: length_stats.n75,
        n90: length_stats.n90,
        n95: length_stats.n95,
        n99: length_stats.n99,
        l10: length_stats.l10,
        l25: length_stats.l25,
        l50: length_stats.l50,
        l75: length_stats.l75,
        l90: length_stats.l90,
        l95: length_stats.l95,
        l99: length_stats.l99,
        au_n: length_stats.au_n,
        effective_contig_count: length_stats.effective_count,
        ungapped_n50: ungapped_length_stats.n50,
        ungapped_n90: ungapped_length_stats.n90,
        ungapped_n95: ungapped_length_stats.n95,
        ungapped_n99: ungapped_length_stats.n99,
        ungapped_au_n: ungapped_length_stats.au_n,
        ungapped_effective_contig_count: ungapped_length_stats.effective_count,
        gap_bases,
        gap_bases_frac,
        longest: length_stats.longest,
        longest_frac,
        gc_bases: composition.gc_bases,
        acgt_bases: composition.acgt_bases,
        n_bases: composition.n_bases,
        ambiguous_bases: composition.ambiguous_bases,
        contigs_with_n: contig_quality.with_n,
        contigs_with_ambiguous: contig_quality.with_ambiguous,
        contigs_all_acgt: contig_quality.all_acgt,
        distinct_sequences: duplicate_summary.distinct_sequences,
        duplicate_sequences: duplicate_summary.duplicate_sequences,
        duplicate_bases: duplicate_summary.duplicate_bases,
        canonical_distinct_sequences: canonical_duplicate_summary.distinct_sequences,
        canonical_duplicate_sequences: canonical_duplicate_summary.duplicate_sequences,
        canonical_duplicate_bases: canonical_duplicate_summary.duplicate_bases,
        contigs_with_n_frac,
        contigs_with_ambiguous_frac,
        contigs_all_acgt_frac,
        distinct_sequences_frac,
        duplicate_sequences_frac,
        duplicate_bases_frac,
        canonical_distinct_sequences_frac,
        canonical_duplicate_sequences_frac,
        canonical_duplicate_bases_frac,
        mean_rle_ratio,
        length_weighted_rle_ratio,
        total_rle_runs,
        gc_content: composition.gc_content(),
        n_content: composition.n_content(length_stats.total_bases),
        ambiguous_content: composition.ambiguous_content(length_stats.total_bases),
        n_runs: n_runs.count,
        longest_n_run: n_runs.longest,
        mean_n_run_length,
        n_runs_per_100kb: per_100kb(n_runs.count, length_stats.total_bases),
        n_bases_per_100kb: per_100kb(composition.n_bases, length_stats.total_bases),
        ambiguous_bases_per_100kb: per_100kb(composition.ambiguous_bases, length_stats.total_bases),
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
    }
}

fn write_assembly_quality_reports(
    output_path: &str,
    quality: AssemblyQualitySummary,
) -> io::Result<(String, String)> {
    let json_path = get_output_filename(output_path, "assembly_metrics.json");
    let tsv_path = get_output_filename(output_path, "assembly_metrics.tsv");

    let json_file = fs::File::create(&json_path)?;
    serde_json::to_writer_pretty(json_file, &quality).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("failed to serialize assembly quality metrics to JSON: {err}"),
        )
    })?;

    let mut tsv_file = fs::File::create(&tsv_path)?;
    writeln!(tsv_file, "metric\tvalue")?;

    let rows = [
        ("total_contigs", quality.total_contigs.to_string()),
        ("total_bases", quality.total_bases.to_string()),
        (
            "ungapped_total_bases",
            quality.ungapped_total_bases.to_string(),
        ),
        ("avg_length", format!("{:.12}", quality.avg_length)),
        ("median_length", format!("{:.12}", quality.median_length)),
        ("n10", quality.n10.to_string()),
        ("n25", quality.n25.to_string()),
        ("n50", quality.n50.to_string()),
        ("n75", quality.n75.to_string()),
        ("n90", quality.n90.to_string()),
        ("n95", quality.n95.to_string()),
        ("n99", quality.n99.to_string()),
        ("l10", quality.l10.to_string()),
        ("l25", quality.l25.to_string()),
        ("l50", quality.l50.to_string()),
        ("l75", quality.l75.to_string()),
        ("l90", quality.l90.to_string()),
        ("l95", quality.l95.to_string()),
        ("l99", quality.l99.to_string()),
        ("au_n", format!("{:.12}", quality.au_n)),
        (
            "effective_contig_count",
            format!("{:.12}", quality.effective_contig_count),
        ),
        ("ungapped_n50", quality.ungapped_n50.to_string()),
        ("ungapped_n90", quality.ungapped_n90.to_string()),
        ("ungapped_n95", quality.ungapped_n95.to_string()),
        ("ungapped_n99", quality.ungapped_n99.to_string()),
        ("ungapped_au_n", format!("{:.12}", quality.ungapped_au_n)),
        (
            "ungapped_effective_contig_count",
            format!("{:.12}", quality.ungapped_effective_contig_count),
        ),
        ("gap_bases", quality.gap_bases.to_string()),
        ("gap_bases_frac", format!("{:.12}", quality.gap_bases_frac)),
        ("longest", quality.longest.to_string()),
        ("longest_frac", format!("{:.12}", quality.longest_frac)),
        ("gc_bases", quality.gc_bases.to_string()),
        ("acgt_bases", quality.acgt_bases.to_string()),
        ("n_bases", quality.n_bases.to_string()),
        ("ambiguous_bases", quality.ambiguous_bases.to_string()),
        ("contigs_with_n", quality.contigs_with_n.to_string()),
        (
            "contigs_with_ambiguous",
            quality.contigs_with_ambiguous.to_string(),
        ),
        ("contigs_all_acgt", quality.contigs_all_acgt.to_string()),
        ("distinct_sequences", quality.distinct_sequences.to_string()),
        (
            "duplicate_sequences",
            quality.duplicate_sequences.to_string(),
        ),
        ("duplicate_bases", quality.duplicate_bases.to_string()),
        (
            "canonical_distinct_sequences",
            quality.canonical_distinct_sequences.to_string(),
        ),
        (
            "canonical_duplicate_sequences",
            quality.canonical_duplicate_sequences.to_string(),
        ),
        (
            "canonical_duplicate_bases",
            quality.canonical_duplicate_bases.to_string(),
        ),
        (
            "contigs_with_n_frac",
            format!("{:.12}", quality.contigs_with_n_frac),
        ),
        (
            "contigs_with_ambiguous_frac",
            format!("{:.12}", quality.contigs_with_ambiguous_frac),
        ),
        (
            "contigs_all_acgt_frac",
            format!("{:.12}", quality.contigs_all_acgt_frac),
        ),
        (
            "distinct_sequences_frac",
            format!("{:.12}", quality.distinct_sequences_frac),
        ),
        (
            "duplicate_sequences_frac",
            format!("{:.12}", quality.duplicate_sequences_frac),
        ),
        (
            "duplicate_bases_frac",
            format!("{:.12}", quality.duplicate_bases_frac),
        ),
        (
            "canonical_distinct_sequences_frac",
            format!("{:.12}", quality.canonical_distinct_sequences_frac),
        ),
        (
            "canonical_duplicate_sequences_frac",
            format!("{:.12}", quality.canonical_duplicate_sequences_frac),
        ),
        (
            "canonical_duplicate_bases_frac",
            format!("{:.12}", quality.canonical_duplicate_bases_frac),
        ),
        ("mean_rle_ratio", format!("{:.12}", quality.mean_rle_ratio)),
        (
            "length_weighted_rle_ratio",
            format!("{:.12}", quality.length_weighted_rle_ratio),
        ),
        ("total_rle_runs", quality.total_rle_runs.to_string()),
        ("gc_content", format!("{:.12}", quality.gc_content)),
        ("n_content", format!("{:.12}", quality.n_content)),
        (
            "ambiguous_content",
            format!("{:.12}", quality.ambiguous_content),
        ),
        ("n_runs", quality.n_runs.to_string()),
        ("longest_n_run", quality.longest_n_run.to_string()),
        (
            "mean_n_run_length",
            format!("{:.12}", quality.mean_n_run_length),
        ),
        (
            "n_runs_per_100kb",
            format!("{:.12}", quality.n_runs_per_100kb),
        ),
        (
            "n_bases_per_100kb",
            format!("{:.12}", quality.n_bases_per_100kb),
        ),
        (
            "ambiguous_bases_per_100kb",
            format!("{:.12}", quality.ambiguous_bases_per_100kb),
        ),
        ("contigs_ge_1kb", quality.contigs_ge_1kb.to_string()),
        ("contigs_ge_10kb", quality.contigs_ge_10kb.to_string()),
        ("contigs_ge_50kb", quality.contigs_ge_50kb.to_string()),
        ("contigs_ge_100kb", quality.contigs_ge_100kb.to_string()),
        ("contigs_ge_1mb", quality.contigs_ge_1mb.to_string()),
        ("bases_ge_1kb", quality.bases_ge_1kb.to_string()),
        ("bases_ge_10kb", quality.bases_ge_10kb.to_string()),
        ("bases_ge_50kb", quality.bases_ge_50kb.to_string()),
        ("bases_ge_100kb", quality.bases_ge_100kb.to_string()),
        ("bases_ge_1mb", quality.bases_ge_1mb.to_string()),
        (
            "contigs_ge_1kb_frac",
            format!("{:.12}", quality.contigs_ge_1kb_frac),
        ),
        (
            "contigs_ge_10kb_frac",
            format!("{:.12}", quality.contigs_ge_10kb_frac),
        ),
        (
            "contigs_ge_50kb_frac",
            format!("{:.12}", quality.contigs_ge_50kb_frac),
        ),
        (
            "contigs_ge_100kb_frac",
            format!("{:.12}", quality.contigs_ge_100kb_frac),
        ),
        (
            "contigs_ge_1mb_frac",
            format!("{:.12}", quality.contigs_ge_1mb_frac),
        ),
        (
            "bases_ge_1kb_frac",
            format!("{:.12}", quality.bases_ge_1kb_frac),
        ),
        (
            "bases_ge_10kb_frac",
            format!("{:.12}", quality.bases_ge_10kb_frac),
        ),
        (
            "bases_ge_50kb_frac",
            format!("{:.12}", quality.bases_ge_50kb_frac),
        ),
        (
            "bases_ge_100kb_frac",
            format!("{:.12}", quality.bases_ge_100kb_frac),
        ),
        (
            "bases_ge_1mb_frac",
            format!("{:.12}", quality.bases_ge_1mb_frac),
        ),
    ];

    for (metric, value) in rows {
        writeln!(tsv_file, "{}\t{}", metric, value)?;
    }

    Ok((json_path, tsv_path))
}

#[inline]
fn should_try_read_overlap_rescue(contigs: &[Contig], read_len: usize) -> bool {
    if read_len == 0 {
        return false;
    }
    let longest = contigs
        .iter()
        .map(|contig| contig.sequence.len())
        .max()
        .unwrap_or(0);
    longest < read_len
}

fn assemble_read_overlap_contigs(sequences: &[String], min_len: usize) -> Vec<Contig> {
    if sequences.is_empty() {
        return Vec::new();
    }

    let mut reads: Vec<&str> = sequences.iter().map(String::as_str).collect();
    reads.sort_unstable();
    reads.dedup();

    let min_read_len = reads.iter().map(|seq| seq.len()).min().unwrap_or(0);
    if min_read_len < 2 {
        return Vec::new();
    }

    let min_overlap = (min_read_len / 3).max(8).min(min_read_len - 1);
    let prefix_index = build_read_prefix_index(&reads, min_overlap, min_read_len);
    let mut successors: Vec<Vec<(usize, usize)>> = vec![Vec::new(); reads.len()];
    let mut predecessors: Vec<Vec<(usize, usize)>> = vec![Vec::new(); reads.len()];

    for (left_idx, left) in reads.iter().enumerate() {
        for overlap in (min_overlap..=left.len().min(min_read_len)).rev() {
            let suffix = &left[left.len() - overlap..];
            let Some(candidates) = prefix_index.get(suffix) else {
                continue;
            };
            for &right_idx in candidates {
                if left_idx == right_idx {
                    continue;
                }
                successors[left_idx].push((right_idx, overlap));
                predecessors[right_idx].push((left_idx, overlap));
            }
        }
    }

    for edges in &mut successors {
        edges.sort_unstable_by(|a, b| {
            b.1.cmp(&a.1)
                .then_with(|| read_tie_break_key(reads[b.0]).cmp(&read_tie_break_key(reads[a.0])))
                .then_with(|| a.0.cmp(&b.0))
        });
    }
    for edges in &mut predecessors {
        edges.sort_unstable_by(|a, b| {
            b.1.cmp(&a.1)
                .then_with(|| read_tie_break_key(reads[b.0]).cmp(&read_tie_break_key(reads[a.0])))
                .then_with(|| a.0.cmp(&b.0))
        });
    }

    let mut stitched = Vec::with_capacity(reads.len());
    for seed_idx in 0..reads.len() {
        stitched.push(stitch_read_overlap_path(
            seed_idx,
            &reads,
            &successors,
            &predecessors,
        ));
    }
    remove_contained_sequences(&mut stitched);

    let mut contigs = Vec::new();
    for sequence in stitched {
        if sequence.len() >= min_len {
            contigs.push(Contig {
                id: contigs.len(),
                sequence,
                kmer_path: Vec::new(),
            });
        }
    }
    contigs
}

fn build_read_prefix_index<'a>(
    reads: &[&'a str],
    min_overlap: usize,
    max_overlap: usize,
) -> AHashMap<&'a str, Vec<usize>> {
    let mut index: AHashMap<&'a str, Vec<usize>> = AHashMap::new();
    for (idx, read) in reads.iter().enumerate() {
        let upper = read.len().min(max_overlap);
        for overlap in min_overlap..=upper {
            index.entry(&read[..overlap]).or_default().push(idx);
        }
    }
    index
}

#[inline]
fn read_tie_break_key(read: &str) -> (usize, &str) {
    (read.len(), read)
}

fn exact_suffix_prefix_overlap(left: &str, right: &str, min_overlap: usize) -> Option<usize> {
    let max_overlap = left.len().min(right.len());
    if max_overlap < min_overlap {
        return None;
    }

    (min_overlap..=max_overlap)
        .rev()
        .find(|&overlap| left[left.len() - overlap..] == right[..overlap])
}

fn stitch_read_overlap_path(
    seed_idx: usize,
    reads: &[&str],
    successors: &[Vec<(usize, usize)>],
    predecessors: &[Vec<(usize, usize)>],
) -> String {
    let mut used = vec![false; reads.len()];
    used[seed_idx] = true;

    let mut contig = reads[seed_idx].to_string();

    let mut current = seed_idx;
    while let Some((prev, overlap)) = predecessors[current]
        .iter()
        .copied()
        .find(|(candidate, _)| !used[*candidate])
    {
        let prefix = &reads[prev][..reads[prev].len() - overlap];
        let mut extended = String::with_capacity(prefix.len() + contig.len());
        extended.push_str(prefix);
        extended.push_str(&contig);
        contig = extended;
        used[prev] = true;
        current = prev;
    }

    current = seed_idx;
    while let Some((next, overlap)) = successors[current]
        .iter()
        .copied()
        .find(|(candidate, _)| !used[*candidate])
    {
        contig.push_str(&reads[next][overlap..]);
        used[next] = true;
        current = next;
    }

    contig
}

fn remove_contained_sequences(sequences: &mut Vec<String>) {
    sequences.sort_unstable_by(|a, b| b.len().cmp(&a.len()).then_with(|| a.cmp(b)));
    let mut kept: Vec<String> = Vec::with_capacity(sequences.len());
    'candidate: for sequence in sequences.drain(..) {
        for existing in &kept {
            if existing.contains(&sequence) {
                continue 'candidate;
            }
        }
        kept.push(sequence);
    }
    *sequences = kept;
}

fn maybe_rescue_fragmented_contigs_with_read_overlaps(
    contigs: Vec<Contig>,
    sequences: &[String],
    min_len: usize,
) -> Vec<Contig> {
    if sequences.len() > READ_OVERLAP_RESCUE_MAX_READS {
        return contigs;
    }

    let min_read_len = sequences.iter().map(|seq| seq.len()).min().unwrap_or(0);
    if !should_try_read_overlap_rescue(&contigs, min_read_len) {
        return contigs;
    }

    let rescued = assemble_read_overlap_contigs(sequences, min_len);
    let original_longest = contigs
        .iter()
        .map(|contig| contig.sequence.len())
        .max()
        .unwrap_or(0);
    let rescued_longest = rescued
        .iter()
        .map(|contig| contig.sequence.len())
        .max()
        .unwrap_or(0);

    if rescued_longest > original_longest {
        info!(
            "Read-overlap rescue improved longest contig from {} bp to {} bp ({} -> {} contigs)",
            original_longest,
            rescued_longest,
            contigs.len(),
            rescued.len()
        );
        rescued
    } else {
        contigs
    }
}

fn split_contigs_at_paired_start_gaps(
    contigs: Vec<Contig>,
    paired_sequences: &[(String, String)],
    min_len: usize,
) -> Vec<Contig> {
    let max_read_len = paired_sequences
        .iter()
        .map(|(read1, read2)| read1.len().max(read2.len()))
        .max()
        .unwrap_or(0);
    let mut split_contigs = Vec::with_capacity(contigs.len());
    let mut split_count = 0usize;
    for contig in contigs {
        if !should_try_paired_start_gap_split(&contig.sequence, max_read_len) {
            split_contigs.push(contig);
            continue;
        }
        let pieces = paired_start_gap_splits(&contig.sequence, paired_sequences, min_len);
        if pieces.len() <= 1 {
            split_contigs.push(contig);
            continue;
        }
        split_count += pieces.len() - 1;
        for sequence in pieces {
            split_contigs.push(Contig {
                id: split_contigs.len(),
                sequence,
                kmer_path: Vec::new(),
            });
        }
    }
    if split_count > 0 {
        info!(
            "Split {} paired-supported compact-overlap contig segment(s) by read-pair start gaps",
            split_count
        );
    }
    split_contigs
}

fn should_try_paired_start_gap_split(contig: &str, max_read_len: usize) -> bool {
    max_read_len > 0 && contig.len() <= max_read_len.saturating_mul(5)
}

fn count_kmers_u64_filtered_with_optional_gpu(
    cpu_backend: &CpuBackend,
    sequences: &[String],
    k: usize,
    min_kmer_count: u32,
    use_gpu: bool,
) -> AHashMap<u64, u32> {
    if use_gpu {
        #[cfg(feature = "gpu")]
        {
            match count_kmers_u64_filtered_gpu(sequences, k, min_kmer_count) {
                Ok(counts) => {
                    info!(
                        "GPU k-mer counting produced {} filtered canonical k-mers",
                        counts.len()
                    );
                    return counts;
                }
                Err(err) => {
                    warn!("GPU k-mer counting failed, falling back to CPU: {}", err);
                }
            }
        }
        #[cfg(not(feature = "gpu"))]
        {
            warn!("GPU requested but binary was not built with --features gpu; using CPU k-mer counting");
        }
    }

    info!("Using CPU k-mer counting and CPU graph build");
    cpu_backend.count_kmers_u64_filtered(sequences, k, min_kmer_count)
}

#[cfg(feature = "gpu")]
fn count_kmers_u64_filtered_gpu(
    sequences: &[String],
    k: usize,
    min_kmer_count: u32,
) -> Result<AHashMap<u64, u32>, String> {
    use crate::gpu::kmer_gpu::GpuKmerCounter;
    use crate::kmer::kmer::KmerU64;

    if !(1..=32).contains(&k) || sequences.is_empty() {
        return Ok(AHashMap::new());
    }

    let counter = GpuKmerCounter::new(k, sequences.len(), 24)?;
    let string_counts = counter.count(sequences)?;
    let mut counts = AHashMap::with_capacity(string_counts.len());
    for (kmer, count) in string_counts {
        if count < min_kmer_count {
            continue;
        }
        let Some(encoded) = KmerU64::from_slice(kmer.as_bytes()) else {
            continue;
        };
        counts.insert(encoded.canonical().encoded, count);
    }
    Ok(counts)
}

fn paired_start_gap_splits(
    contig: &str,
    paired_sequences: &[(String, String)],
    min_len: usize,
) -> Vec<String> {
    let mut spans = Vec::new();
    let mut max_read_len = 0usize;
    for (read1, read2) in paired_sequences {
        max_read_len = max_read_len.max(read1.len()).max(read2.len());
        let Some(read1_start) = contig.find(read1) else {
            continue;
        };
        let Some(read2_start) = contig.find(read2) else {
            continue;
        };
        let start = read1_start.min(read2_start);
        let end = (read1_start + read1.len()).max(read2_start + read2.len());
        spans.push((start, end));
    }
    if spans.len() < 4 || max_read_len == 0 {
        return vec![contig.to_string()];
    }
    spans.sort_unstable();
    spans.dedup();

    let gap_threshold = ((max_read_len * 2) / 3).max(1);
    let mut groups: Vec<Vec<(usize, usize)>> = Vec::new();
    for span in spans {
        if let Some(group) = groups.last_mut() {
            let previous_start = group.last().map(|(start, _)| *start).unwrap_or(span.0);
            if span.0.saturating_sub(previous_start) <= gap_threshold {
                group.push(span);
                continue;
            }
        }
        groups.push(vec![span]);
    }

    if groups.len() <= 1 {
        return vec![contig.to_string()];
    }

    let mut pieces = Vec::with_capacity(groups.len());
    for group in groups {
        let start = group.iter().map(|(start, _)| *start).min().unwrap_or(0);
        let end = group
            .iter()
            .map(|(_, end)| *end)
            .max()
            .unwrap_or(start)
            .min(contig.len());
        if end.saturating_sub(start) < min_len {
            return vec![contig.to_string()];
        }
        pieces.push(contig[start..end].to_string());
    }

    if pieces.iter().any(|piece| piece.len() == contig.len()) {
        vec![contig.to_string()]
    } else {
        pieces
    }
}

#[inline]
fn sequence_only_record(record: FastqRecord) -> FastqRecord {
    FastqRecord {
        header: String::new(),
        sequence: record.sequence,
        plus: String::new(),
        quality: String::new(),
    }
}

const ESTIMATED_RECORDS_PER_MB: usize = 5_000;
const DEFAULT_SEQUENCE_CAPACITY: usize = 100_000;
const MIN_SEQUENCE_CAPACITY: usize = 10_000;
const MAX_SEQUENCE_CAPACITY: usize = 2_000_000;
const RLE_RATIO_SCALE: u128 = 1_000_000_000_000;

#[inline]
fn estimate_sequence_capacity(file_size_bytes: Option<u64>) -> usize {
    let estimated = file_size_bytes
        .map(|bytes| {
            let mb = (bytes / (1024 * 1024)).min(usize::MAX as u64) as usize;
            mb.saturating_mul(ESTIMATED_RECORDS_PER_MB)
        })
        .unwrap_or(DEFAULT_SEQUENCE_CAPACITY);

    estimated.clamp(MIN_SEQUENCE_CAPACITY, MAX_SEQUENCE_CAPACITY)
}

pub fn assemble_reads(
    input_path: &str,
    input2_path: Option<&str>,
    output_path: &str,
    min_len: usize,
    _output_gfa: bool,
    _output_gfa2: bool,
    _adaptive_k: bool,
    _use_rle: bool,
    collapse_repeats: bool,
    min_repeat_len: usize,
    polish: bool,
    polish_window: usize,
    streaming: bool,
    export_metadata: bool,
    json_metadata: Option<String>,
    tsv_metadata: Option<String>,
    isoforms: bool,
    gtf_path: Option<String>,
    gff3_path: Option<String>,
    max_path_depth: usize,
    min_confidence: f64,
    compute_tpm: bool,
    polish_isoforms: bool,
    samples_path: Option<String>,
    min_tpm: f64,
    long_reads: Option<String>,
    counts_matrix: bool,
) -> io::Result<()> {
    // Default to CPU backend for backwards compatibility
    assemble_reads_with_gpu(
        input_path,
        input2_path,
        output_path,
        min_len,
        _output_gfa,
        _output_gfa2,
        _adaptive_k,
        _use_rle,
        collapse_repeats,
        min_repeat_len,
        polish,
        polish_window,
        streaming,
        export_metadata,
        json_metadata,
        tsv_metadata,
        isoforms,
        gtf_path,
        gff3_path,
        max_path_depth,
        min_confidence,
        compute_tpm,
        polish_isoforms,
        samples_path,
        min_tpm,
        long_reads,
        counts_matrix,
        false, // use_gpu = false by default
    )
}

/// Assemble reads with optional GPU acceleration
pub fn assemble_reads_with_gpu(
    input_path: &str,
    input2_path: Option<&str>,
    output_path: &str,
    min_len: usize,
    _output_gfa: bool,
    _output_gfa2: bool,
    _adaptive_k: bool,
    _use_rle: bool,
    collapse_repeats: bool,
    min_repeat_len: usize,
    polish: bool,
    polish_window: usize,
    _streaming: bool,
    export_metadata: bool,
    json_metadata: Option<String>,
    tsv_metadata: Option<String>,
    isoforms: bool,
    gtf_path: Option<String>,
    gff3_path: Option<String>,
    max_path_depth: usize,
    min_confidence: f64,
    compute_tpm: bool,
    polish_isoforms: bool,
    samples_path: Option<String>,
    min_tpm: f64,
    long_reads: Option<String>,
    counts_matrix: bool,
    use_gpu: bool,
) -> io::Result<()> {
    info!("Starting assembly from: {}", input_path);
    if let Some(input2_path) = input2_path {
        info!("Using paired-end mate input: {}", input2_path);
    }

    // Determine k-mer size - either adaptive or fixed optimal
    let max_k = 41;
    let mut k: usize;

    let mut records: Vec<FastqRecord> = Vec::new();
    let need_sequence_only_records = polish || (isoforms && (polish_isoforms || compute_tpm));

    // Use streaming mode for memory efficiency with large files
    let cpu_backend = CpuBackend::new();

    // Stream sequences and count k-mers in chunks for memory efficiency
    info!("Streaming FASTQ records for k-mer counting...");
    let reader = try_open_fastq(input_path)?;

    // Pre-size sequences vector based on file size, with hard caps to avoid
    // pathological over-allocation for very large inputs.
    let estimated_count =
        estimate_sequence_capacity(std::fs::metadata(input_path).ok().map(|m| m.len()));

    let mut sequences: Vec<String> = Vec::with_capacity(estimated_count);
    let mut paired_sequences: Vec<(String, String)> = Vec::new();
    if let Some(input2_path) = input2_path {
        let reader2 = try_open_fastq(input2_path)?;
        for pair in stream_paired_fastq_records_checked(reader, reader2) {
            let (r1, r2) = pair?;
            let mate_sequence = reverse_complement_dna(&r2.sequence);
            paired_sequences.push((r1.sequence.clone(), mate_sequence.clone()));
            sequences.push(r1.sequence);
            sequences.push(mate_sequence);
        }
    } else {
        // First pass: collect sequences for k optimization (sample if large)
        for record in stream_fastq_records_checked(reader) {
            sequences.push(record?.sequence);
        }
    }

    let num_sequences = sequences.len();
    info!("Loaded {} sequences", num_sequences);
    let can_try_read_overlap_rescue = num_sequences <= READ_OVERLAP_RESCUE_MAX_READS;

    if use_gpu {
        info!("GPU requested for k-mer counting; graph build remains CPU");
    }

    // Determine k-mer size from a bounded prefix sample without cloning.
    let sample_size = sequences.len().min(10_000);
    let sample = &sequences[..sample_size];

    if _adaptive_k {
        let hist = kmer_coverage_histogram(sample, max_k);
        k = select_best_k(&hist);
        info!("Using adaptive k-mer: {}", k);
    } else {
        k = optimal_k(sample, max_k);
        info!("Using optimal k-mer size: {}", k);
    }

    // Ensure k <= 32 for u64 encoding
    if k > 32 {
        info!("Capping k-mer size at 32 for u64 encoding (was {})", k);
        k = 32;
    }

    // Count k-mers using optimized u64 path with Bloom filter pre-filtering
    // This filters singleton k-mers (sequencing errors) for 30-50% memory reduction
    info!(
        "Counting k-mers with k={} using optimized u64 encoding with Bloom filter",
        k
    );
    let min_kmer_count = 2; // Filter k-mers appearing only once
    let kmer_counts_u64 = count_kmers_u64_filtered_with_optional_gpu(
        &cpu_backend,
        &sequences,
        k,
        min_kmer_count,
        use_gpu,
    );
    info!(
        "Found {} unique k-mers (after filtering singletons)",
        kmer_counts_u64.len()
    );

    if need_sequence_only_records {
        records = sequences
            .iter()
            .map(|sequence| {
                sequence_only_record(FastqRecord {
                    header: String::new(),
                    sequence: sequence.clone(),
                    plus: String::new(),
                    quality: String::new(),
                })
            })
            .collect();
        info!(
            "Reused {} input sequences for polishing/TPM (no second FASTQ pass)",
            records.len()
        );
    } else if !can_try_read_overlap_rescue {
        // Sequence strings are no longer needed after k-mer counting.
        // Releasing this buffer early reduces peak RSS during graph cleanup/polishing.
        let released_records = sequences.len();
        sequences.clear();
        sequences.shrink_to_fit();
        info!(
            "Released {} input sequences from memory after k-mer counting",
            released_records
        );
    }

    // Build adjacency table for assembly
    info!("Building adjacency table...");
    let mut adjacency = cpu_backend.build_adjacency_u64(&kmer_counts_u64, k);
    info!(
        "Adjacency table built with {} forward edges",
        adjacency.forward.len()
    );

    // Clean up graph: remove tips and collapse bubbles
    // This reduces noise from sequencing errors and improves assembly quality
    {
        use crate::graph::assembler::cleanup_graph;
        info!("Cleaning up assembly graph (removing tips and bubbles)...");
        let (tips_removed, bubbles_collapsed) =
            cleanup_graph(&mut adjacency, &kmer_counts_u64, k, min_kmer_count);
        info!(
            "Graph cleanup: removed {} tips, collapsed {} bubbles",
            tips_removed, bubbles_collapsed
        );
    }

    // Perform greedy assembly using u64 k-mers
    info!("Assembling contigs with minimum length: {}", min_len);
    let mut contigs = greedy_assembly_u64(k, &kmer_counts_u64, &adjacency, min_len);
    if can_try_read_overlap_rescue {
        contigs = maybe_rescue_fragmented_contigs_with_read_overlaps(contigs, &sequences, min_len);
    }

    if !need_sequence_only_records && !sequences.is_empty() {
        let released_records = sequences.len();
        sequences.clear();
        sequences.shrink_to_fit();
        info!(
            "Released {} input sequences from memory after contig construction",
            released_records
        );
    }

    // Collapse repeats if requested
    if collapse_repeats {
        use crate::graph::simplify::collapse_repeats;
        let before_count = contigs.len();
        contigs = collapse_repeats(contigs, min_repeat_len);
        info!(
            "Collapsed repeats: {} -> {} contigs (removed {})",
            before_count,
            contigs.len(),
            before_count - contigs.len()
        );
    }

    if !paired_sequences.is_empty() {
        contigs = split_contigs_at_paired_start_gaps(contigs, &paired_sequences, min_len);
    }

    // Polish contigs if requested (using parallel implementation for 4-8x speedup)
    if polish {
        use crate::graph::polish::polish_contig_parallel;
        info!(
            "Polishing contigs using aligned reads (window size: {}, parallel)",
            polish_window
        );

        let chunk_size = 10000; // Process 10kb chunks in parallel
        for contig in &mut contigs {
            contig.sequence =
                polish_contig_parallel(&contig.sequence, &records, polish_window, chunk_size);
        }

        info!("Completed contig polishing");
    }

    canonicalize_contig_output_order(&mut contigs);

    let quality = summarize_assembly_quality(&contigs);
    info!(
        "Contig statistics: {} contigs, {} bp total, Mean/Median: {:.1}/{:.1} bp, N10/N25/N50/N75/N90/N95/N99: {}/{}/{}/{}/{}/{}/{} bp, L10/L25/L50/L75/L90/L95/L99: {}/{}/{}/{}/{}/{}/{}, auN/effective count: {:.1}/{:.2}, Ungapped bases/gaps/N50/N90/N95/N99/auN/effective count: {}/{}/{}/{}/{}/{}/{:.1}/{:.2}, Longest: {} bp ({:.2}%), GC/N/Ambiguous bases: {}/{}/{} (fractions {:.2}%/{:.2}%/{:.2}%), Contigs with N/Ambiguous/All-ACGT: {}/{}/{} ({:.2}%/{:.2}%/{:.2}%), RLE mean/weighted runs ratio: {:.4}/{:.4} ({} runs), N-runs: count {}, max {}, mean {:.1} bp ({:.1} per 100kb), N/Ambiguous bases per 100kb: {:.1}/{:.1}, >=1kb/10kb/50kb/100kb/1mb contigs: {}/{}/{}/{}/{} ({:.1}%/{:.1}%/{:.1}%/{:.1}%/{:.1}%), span: {}/{}/{}/{}/{} bp ({:.1}%/{:.1}%/{:.1}%/{:.1}%/{:.1}%)",
        quality.total_contigs,
        quality.total_bases,
        quality.avg_length,
        quality.median_length,
        quality.n10,
        quality.n25,
        quality.n50,
        quality.n75,
        quality.n90,
        quality.n95,
        quality.n99,
        quality.l10,
        quality.l25,
        quality.l50,
        quality.l75,
        quality.l90,
        quality.l95,
        quality.l99,
        quality.au_n,
        quality.effective_contig_count,
        quality.ungapped_total_bases,
        quality.gap_bases,
        quality.ungapped_n50,
        quality.ungapped_n90,
        quality.ungapped_n95,
        quality.ungapped_n99,
        quality.ungapped_au_n,
        quality.ungapped_effective_contig_count,
        quality.longest,
        quality.longest_frac * 100.0,
        quality.gc_bases,
        quality.n_bases,
        quality.ambiguous_bases,
        quality.gc_content * 100.0,
        quality.n_content * 100.0,
        quality.ambiguous_content * 100.0,
        quality.contigs_with_n,
        quality.contigs_with_ambiguous,
        quality.contigs_all_acgt,
        quality.contigs_with_n_frac * 100.0,
        quality.contigs_with_ambiguous_frac * 100.0,
        quality.contigs_all_acgt_frac * 100.0,
        quality.mean_rle_ratio,
        quality.length_weighted_rle_ratio,
        quality.total_rle_runs,
        quality.n_runs,
        quality.longest_n_run,
        quality.mean_n_run_length,
        quality.n_runs_per_100kb,
        quality.n_bases_per_100kb,
        quality.ambiguous_bases_per_100kb,
        quality.contigs_ge_1kb,
        quality.contigs_ge_10kb,
        quality.contigs_ge_50kb,
        quality.contigs_ge_100kb,
        quality.contigs_ge_1mb,
        quality.contigs_ge_1kb_frac * 100.0,
        quality.contigs_ge_10kb_frac * 100.0,
        quality.contigs_ge_50kb_frac * 100.0,
        quality.contigs_ge_100kb_frac * 100.0,
        quality.contigs_ge_1mb_frac * 100.0,
        quality.bases_ge_1kb,
        quality.bases_ge_10kb,
        quality.bases_ge_50kb,
        quality.bases_ge_100kb,
        quality.bases_ge_1mb,
        quality.bases_ge_1kb_frac * 100.0,
        quality.bases_ge_10kb_frac * 100.0,
        quality.bases_ge_50kb_frac * 100.0,
        quality.bases_ge_100kb_frac * 100.0,
        quality.bases_ge_1mb_frac * 100.0
    );

    let (quality_json_path, quality_tsv_path) =
        write_assembly_quality_reports(output_path, quality)?;
    info!(
        "Assembly quality reports written to {} and {}",
        quality_json_path, quality_tsv_path
    );

    // Write FASTA output
    let mut writer = FastaWriter::try_new(output_path)?;
    for (i, contig) in contigs.iter().enumerate() {
        writer.write_contig(contig, i + 1)?;

        // Write RLE version if requested
        if _use_rle {
            writer.write_rle_contig(contig, i + 1)?;
        }
    }

    info!(
        "Assembly complete: {} contigs written to {}",
        contigs.len(),
        output_path
    );

    // Export metadata if requested
    if export_metadata || json_metadata.is_some() || tsv_metadata.is_some() {
        use crate::io::metadata::generate_metadata;
        info!("Generating contig metadata");
        let meta = generate_metadata(&contigs);

        // Handle standard metadata JSON export
        if export_metadata {
            info!("Exporting contig metadata to JSON");
            let json = serde_json::to_string_pretty(&meta).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("failed to serialize contig metadata to JSON: {err}"),
                )
            })?;
            let meta_path = format!("{}.contig_meta.json", output_path);
            fs::write(&meta_path, json)?;
            info!("Metadata written to {}", meta_path);
        }

        // Handle custom JSON metadata path
        if let Some(path) = &json_metadata {
            info!("Writing JSON metadata to custom path: {}", path);
            let json = serde_json::to_string_pretty(&meta).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("failed to serialize contig metadata for {}: {err}", path),
                )
            })?;
            fs::write(path, json)?;
        }

        // Handle custom TSV metadata path
        if let Some(path) = &tsv_metadata {
            info!("Writing TSV metadata to custom path: {}", path);
            let mut file = std::fs::File::create(path)?;
            writeln!(file, "contig_id\tlength\trle_compression\tgc_content")?;
            for m in &meta {
                writeln!(
                    file,
                    "{}\t{}\t{:.4}\t{:.4}",
                    m.id, m.length, m.rle_compression, m.gc_content
                )?;
            }
        }
    }

    // If GFA or GFA2 output is requested, find overlaps between contigs
    if (_output_gfa || _output_gfa2 || isoforms) && !contigs.is_empty() {
        // Borrow contig sequences to avoid an extra full-sequence clone pass.
        let contig_seqs: Vec<&str> = contigs.iter().map(|c| c.sequence.as_str()).collect();

        // Build overlap graph
        let min_overlap = (k / 2).max(15); // Use at least half of k but minimum 15bp
        let max_mismatches = 3; // Allow up to 3 mismatches in the overlap

        info!("Finding overlaps (single pass)");

        let builder = OverlapGraphBuilder::new(min_overlap, max_mismatches, max_mismatches);
        let graph = builder.build_overlap_graph(&contig_seqs);
        let (_, paths) = builder.stitch_contigs(&graph);

        // Use canonical overlap detection path to build graph links.
        let links = find_overlaps(&contigs, min_overlap, max_mismatches);

        info!("Found {} overlaps", links.len());

        // Process isoforms if requested
        if isoforms {
            info!("Performing isoform inference from assembly graph");

            // Derive per-contig expression support from assembled k-mer paths.
            let kmer_counts_converted = derive_contig_expression_map(&contigs, &kmer_counts_u64);

            let mut transcripts = crate::pipeline::isoform_processor::process_isoforms(
                &contigs,
                &links,
                &kmer_counts_converted,
                max_path_depth,
                min_confidence,
            )
            .unwrap_or_else(|e| {
                warn!("Error processing isoforms: {}", e);
                Vec::new()
            });

            info!("Generated {} transcript isoforms", transcripts.len());

            // Export counts matrix
            if counts_matrix {
                info!("Writing isoform counts matrix");
                let counts_matrix_path = format!("{}_isoform.counts.matrix", output_path);
                if let Err(e) = crate::quant::matrix::write_isoform_counts_matrix(
                    &transcripts,
                    &counts_matrix_path,
                ) {
                    warn!("Failed to write counts matrix: {}", e);
                } else {
                    info!("Counts matrix written to: {}", counts_matrix_path);
                }
            }

            // Export GTF if requested
            if let Some(gtf_path) = &gtf_path {
                info!("Writing isoform GTF to: {}", gtf_path);
                use crate::io::gtf::write_gtf;
                if let Err(e) = write_gtf(&transcripts, gtf_path) {
                    warn!("Failed to write GTF file: {}", e);
                } else {
                    info!(
                        "GTF output complete: {} transcripts written",
                        transcripts.len()
                    );
                }
            }

            // Export GFF3 if requested
            if let Some(gff3_path) = &gff3_path {
                info!("Writing isoform GFF3 to: {}", gff3_path);
                use crate::io::gff3::write_gff3;
                if let Err(e) = write_gff3(&transcripts, gff3_path) {
                    warn!("Failed to write GFF3 file: {}", e);
                } else {
                    info!(
                        "GFF3 output complete: {} transcripts written",
                        transcripts.len()
                    );
                }
            }

            // Polish isoform sequences if requested
            if polish_isoforms {
                info!("Polishing isoform sequences with aligned reads");
                use crate::polish::align::polish_sequence;
                for t in &mut transcripts {
                    t.sequence = polish_sequence(&t.sequence, &records, 25);
                }
                info!("Isoform polishing complete");
            }

            // Compute TPM values if requested
            if compute_tpm {
                info!("Computing TPM expression values for transcripts");
                use crate::quant::tpm::{compute_tpm, count_reads, filter_by_tpm, write_tpm_table};

                let counts = count_reads(&transcripts, &records);
                let tpms = compute_tpm(&counts, &transcripts);

                // Filter transcripts by TPM if requested
                if min_tpm > 0.0 {
                    info!("Filtering transcripts with TPM < {}", min_tpm);
                    let (filtered_transcripts, filtered_tpms) =
                        filter_by_tpm(&transcripts, &tpms, min_tpm);
                    let filtered_count = transcripts.len() - filtered_transcripts.len();
                    info!(
                        "Filtered out {} transcripts with low expression",
                        filtered_count
                    );

                    // Replace transcripts with filtered set
                    transcripts = filtered_transcripts;

                    // Update TPMs to match filtered transcripts
                    let tpm_path = format!("{}_isoform.tpm.tsv", output_path);
                    write_tpm_table(&transcripts, &filtered_tpms, &tpm_path)?;
                    info!("TPM values written to {}", tpm_path);

                    // Write transcript metrics to JSON if requested
                    if let Some(json_path) = &json_metadata {
                        info!("Writing transcript metrics to JSON: {}", json_path);
                        use crate::io::metadata::write_transcript_metrics;
                        if let Err(e) =
                            write_transcript_metrics(&transcripts, &filtered_tpms, json_path)
                        {
                            warn!("Failed to write transcript metrics to JSON: {}", e);
                        } else {
                            info!("Transcript metrics written to {}", json_path);
                        }
                    }
                } else {
                    let tpm_path = format!("{}_isoform.tpm.tsv", output_path);
                    write_tpm_table(&transcripts, &tpms, &tpm_path)?;
                    info!("TPM values written to {}", tpm_path);

                    // Write transcript metrics to JSON if requested
                    if let Some(json_path) = &json_metadata {
                        info!("Writing transcript metrics to JSON: {}", json_path);
                        use crate::io::metadata::write_transcript_metrics;
                        if let Err(e) = write_transcript_metrics(&transcripts, &tpms, json_path) {
                            warn!("Failed to write transcript metrics to JSON: {}", e);
                        } else {
                            info!("Transcript metrics written to {}", json_path);
                        }
                    }
                }

                // Process multiple samples if provided
                if let Some(sample_file) = &samples_path {
                    use crate::io::sam::parse_sam_transcript_hits;
                    use crate::quant::matrix::write_counts_matrix;
                    use std::collections::HashMap;

                    info!("Processing multi-sample data from {}", sample_file);
                    let sample_content = std::fs::read_to_string(sample_file)?;

                    let mut sample_tpms: HashMap<String, Vec<f64>> = HashMap::new();

                    // Process each sample
                    for line in sample_content.lines() {
                        if line.trim().is_empty() || line.starts_with('#') {
                            continue; // Skip empty lines and comments
                        }

                        let parts: Vec<&str> = line.split(',').collect();
                        if parts.len() < 2 {
                            warn!("Invalid sample line: {}", line);
                            continue;
                        }

                        let sample_name = parts[0].trim().to_string();
                        let sam_path = parts[1].trim();

                        info!(
                            "Processing sample: {} from alignment {}",
                            sample_name, sam_path
                        );

                        // Parse SAM file for transcript hits
                        let hits = match parse_sam_transcript_hits(sam_path) {
                            Ok(h) => h,
                            Err(e) => {
                                warn!("Failed to parse SAM file {}: {}", sam_path, e);
                                continue;
                            }
                        };

                        // Count hits per transcript
                        let sam_counts: Vec<usize> = transcripts
                            .iter()
                            .map(|t| *hits.get(&format!("transcript_{}", t.id)).unwrap_or(&0))
                            .collect();

                        // Compute TPM values
                        let sam_tpms = compute_tpm(&sam_counts, &transcripts);
                        sample_tpms.insert(sample_name, sam_tpms);
                    }

                    if !sample_tpms.is_empty() {
                        // Write the matrix output
                        let matrix_path = format!("{}_isoform.counts.matrix", output_path);
                        if let Err(e) =
                            write_counts_matrix(&sample_tpms, &transcripts, &matrix_path)
                        {
                            warn!("Failed to write multi-sample counts matrix: {}", e);
                        } else {
                            info!("Multi-sample counts matrix written to {}", matrix_path);
                        }
                    } else {
                        warn!("No valid samples found in {}", sample_file);
                    }
                }
            }

            // Apply transcript polishing with long reads if requested
            if let Some(polish_sam_path) = &long_reads {
                info!(
                    "Polishing transcript sequences with long read alignments from {}",
                    polish_sam_path
                );
                use crate::polish::longread::polish_transcripts;

                match polish_transcripts(&mut transcripts, polish_sam_path) {
                    Ok(count) => {
                        info!("Successfully polished {} transcripts", count);
                    }
                    Err(e) => {
                        warn!("Error during transcript polishing: {}", e);
                    }
                }
            }

            // Display transcript evaluation metrics
            let mut transcript_lengths: Vec<usize> =
                transcripts.iter().map(|t| t.sequence.len()).collect();
            let stats = evaluate_lengths_in_place(&mut transcript_lengths);
            info!(
        "Transcript statistics: {} transcripts, {} bp total, Avg: {:.1} bp, N10/N25/N50/N75/N90/N95/N99: {}/{}/{}/{}/{}/{}/{} bp, L10/L25/L50/L75/L90/L95/L99: {}/{}/{}/{}/{}/{}/{}, auN: {:.1}",
                stats.total,
                stats.total_bases,
                stats.avg_length,
                stats.n10,
                stats.n25,
                stats.n50,
                stats.n75,
                stats.n90,
                stats.n95,
                stats.n99,
                stats.l10,
                stats.l25,
                stats.l50,
                stats.l75,
                stats.l90,
                stats.l95,
                stats.l99,
                stats.au_n
            );
        }
        // Output GFA if requested
        if _output_gfa {
            let gfa_path = get_output_filename(output_path, "gfa");
            info!("Writing GFA to: {}", gfa_path);

            // Write GFA output
            let mut gfa_writer = GfaWriter::try_new(&gfa_path)?;

            if _use_rle {
                gfa_writer.write_rle_segments(&contigs)?;
            } else {
                gfa_writer.write_segments(&contigs)?;
            }

            gfa_writer.write_links(&links)?;
            gfa_writer.write_assembly_paths(&paths)?;

            info!(
                "GFA output complete: {} segments, {} links, {} paths written",
                contigs.len(),
                links.len(),
                paths.len()
            );
        }

        // Output GFA2 if requested
        if _output_gfa2 {
            let gfa2_path = get_output_filename(output_path, "gfa2");
            info!("Writing GFA2 to: {}", gfa2_path);

            // Write GFA2 output
            let mut gfa2_writer = Gfa2Writer::try_new(&gfa2_path)?;
            gfa2_writer.write_segments(&contigs)?;
            gfa2_writer.write_links(&links)?;
            gfa2_writer.write_paths(&paths)?;

            info!(
                "GFA2 output complete: {} segments, {} links, {} paths written",
                contigs.len(),
                links.len(),
                paths.len()
            );
        }
    }

    Ok(())
}

// Helper function to generate output filenames
fn get_output_filename(output_path: &str, extension: &str) -> String {
    let path = if output_path.ends_with(".gz") {
        output_path.strip_suffix(".gz").unwrap_or(output_path)
    } else {
        output_path
    };

    if path.ends_with(".fasta") || path.ends_with(".fa") {
        format!(
            "{}.{}",
            &path[..path.rfind('.').unwrap_or(path.len())],
            extension
        )
    } else {
        format!("{}.{}", path, extension)
    }
}

fn derive_contig_expression_map(
    contigs: &[crate::graph::assembler::Contig],
    kmer_counts: &ahash::AHashMap<u64, u32>,
) -> std::collections::HashMap<usize, usize> {
    let mut expression_by_contig = std::collections::HashMap::with_capacity(contigs.len());

    for contig in contigs {
        let mut support_sum = 0u128;
        let mut support_obs = 0u128;

        for kmer in &contig.kmer_path {
            if let Some(&count) = kmer_counts.get(kmer) {
                support_sum += count as u128;
                support_obs += 1;
            }
        }

        if support_obs > 0 {
            // Rounded mean support keeps deterministic integer weights for graph scoring.
            let mean_support = ((support_sum + (support_obs / 2)) / support_obs) as usize;
            expression_by_contig.insert(contig.id, mean_support.max(1));
        }
    }

    expression_by_contig
}

#[inline]
fn canonicalize_contig_output_order(contigs: &mut [Contig]) {
    contigs.sort_unstable_by(|a, b| {
        b.sequence
            .len()
            .cmp(&a.sequence.len())
            .then_with(|| a.sequence.cmp(&b.sequence))
            .then_with(|| a.kmer_path.cmp(&b.kmer_path))
            .then_with(|| a.id.cmp(&b.id))
    });

    for (new_id, contig) in contigs.iter_mut().enumerate() {
        contig.id = new_id;
    }
}

#[cfg(test)]
mod tests {
    use super::{
        assemble_read_overlap_contigs, assemble_reads_with_gpu, build_read_prefix_index,
        canonical_sequence_key, canonicalize_contig_output_order, derive_contig_expression_map,
        estimate_sequence_capacity, maybe_rescue_fragmented_contigs_with_read_overlaps,
        sequence_only_record, should_try_paired_start_gap_split,
        split_contigs_at_paired_start_gaps, summarize_assembly_quality,
        write_assembly_quality_reports, AssemblyQualitySummary,
    };
    use crate::graph::assembler::Contig;
    use crate::io::fastq::FastqRecord;
    use ahash::AHashMap;
    use proptest::prelude::*;
    use rand::rngs::StdRng;
    use rand::seq::SliceRandom;
    use rand::SeedableRng;
    use std::io::{self, Write};
    use tempfile::{NamedTempFile, TempDir};

    #[inline]
    fn is_n(base: u8) -> bool {
        matches!(base, b'N' | b'n')
    }

    #[inline]
    fn is_acgt_or_u(base: u8) -> bool {
        matches!(
            base,
            b'A' | b'a' | b'T' | b't' | b'U' | b'u' | b'G' | b'g' | b'C' | b'c'
        )
    }

    #[inline]
    fn is_ambiguous(base: u8) -> bool {
        !is_n(base) && !is_acgt_or_u(base)
    }

    #[inline]
    fn count_rle_runs(sequence: &str) -> usize {
        let mut bytes = sequence.bytes();
        let Some(first) = bytes.next() else {
            return 0;
        };

        let mut runs = 1usize;
        let mut prev = first.to_ascii_uppercase();
        for base in bytes {
            let current = base.to_ascii_uppercase();
            if current != prev {
                runs += 1;
                prev = current;
            }
        }
        runs
    }

    #[inline]
    fn n_run_stats(sequence: &str) -> (usize, usize) {
        let mut count = 0usize;
        let mut longest = 0usize;
        let mut active = 0usize;
        for base in sequence.bytes() {
            if is_n(base) {
                active += 1;
            } else if active > 0 {
                count += 1;
                longest = longest.max(active);
                active = 0;
            }
        }
        if active > 0 {
            count += 1;
            longest = longest.max(active);
        }
        (count, longest)
    }

    #[test]
    fn derive_contig_expression_map_averages_observed_kmer_support() {
        let contigs = vec![
            Contig {
                id: 7,
                sequence: "AAAC".to_string(),
                kmer_path: vec![11, 12, 13],
            },
            Contig {
                id: 8,
                sequence: "GGG".to_string(),
                kmer_path: vec![99],
            },
            Contig {
                id: 9,
                sequence: "TTT".to_string(),
                kmer_path: vec![],
            },
        ];

        let mut kmer_counts = AHashMap::new();
        kmer_counts.insert(11, 3);
        kmer_counts.insert(12, 5);
        kmer_counts.insert(13, 4);

        let expression = derive_contig_expression_map(&contigs, &kmer_counts);

        assert_eq!(expression.get(&7), Some(&4));
        assert!(!expression.contains_key(&8));
        assert!(!expression.contains_key(&9));
    }

    #[test]
    fn canonicalize_contig_output_order_sorts_and_reindexes_deterministically() {
        let mut contigs = vec![
            Contig {
                id: 9,
                sequence: "GGGG".to_string(),
                kmer_path: vec![3, 4],
            },
            Contig {
                id: 4,
                sequence: "AAAA".to_string(),
                kmer_path: vec![1, 2],
            },
            Contig {
                id: 7,
                sequence: "AAA".to_string(),
                kmer_path: vec![8],
            },
            Contig {
                id: 2,
                sequence: "AAAA".to_string(),
                kmer_path: vec![1, 3],
            },
        ];

        canonicalize_contig_output_order(&mut contigs);

        let fingerprints: Vec<(usize, &str, &[u64])> = contigs
            .iter()
            .map(|contig| {
                (
                    contig.id,
                    contig.sequence.as_str(),
                    contig.kmer_path.as_slice(),
                )
            })
            .collect();
        assert_eq!(
            fingerprints,
            vec![
                (0, "AAAA", &[1, 2][..]),
                (1, "AAAA", &[1, 3][..]),
                (2, "GGGG", &[3, 4][..]),
                (3, "AAA", &[8][..]),
            ]
        );
    }

    #[test]
    fn paired_start_gap_split_breaks_compact_overlap_fusion() {
        fn pseudo_dna(seed: u64, len: usize) -> String {
            let mut state = seed;
            let mut sequence = String::with_capacity(len);
            for _ in 0..len {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                let base = match (state >> 32) & 3 {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    _ => 'T',
                };
                sequence.push(base);
            }
            sequence
        }

        let left_unique = pseudo_dna(11, 120);
        let shared = pseudo_dna(29, 96);
        let right_unique = pseudo_dna(47, 120);
        let fused = format!("{left_unique}{shared}{right_unique}");
        let paired_sequences = vec![
            (fused[0..75].to_string(), fused[85..160].to_string()),
            (fused[48..123].to_string(), fused[133..208].to_string()),
            (fused[56..131].to_string(), fused[141..216].to_string()),
            (fused[120..195].to_string(), fused[205..280].to_string()),
            (fused[168..243].to_string(), fused[253..328].to_string()),
            (fused[176..251].to_string(), fused[261..336].to_string()),
        ];
        let contigs = vec![Contig {
            id: 0,
            sequence: fused,
            kmer_path: Vec::new(),
        }];

        let split = split_contigs_at_paired_start_gaps(contigs, &paired_sequences, 25);
        let lengths: Vec<usize> = split.iter().map(|contig| contig.sequence.len()).collect();

        assert_eq!(lengths, vec![216, 216]);
        assert_eq!(split[0].sequence, format!("{left_unique}{shared}"));
        assert_eq!(split[1].sequence, format!("{shared}{right_unique}"));
    }

    #[test]
    fn paired_start_gap_split_skips_long_public_scale_contigs() {
        assert!(should_try_paired_start_gap_split(&"A".repeat(380), 76));
        assert!(!should_try_paired_start_gap_split(&"A".repeat(381), 76));
    }

    #[test]
    fn sequence_only_record_drops_non_sequence_fields() {
        let record = FastqRecord {
            header: "@read_1".to_string(),
            sequence: "ACGT".to_string(),
            plus: "+".to_string(),
            quality: "IIII".to_string(),
        };

        let compact = sequence_only_record(record);
        assert!(compact.header.is_empty());
        assert_eq!(compact.sequence, "ACGT");
        assert!(compact.plus.is_empty());
        assert!(compact.quality.is_empty());
    }

    #[test]
    fn read_overlap_rescue_builds_longer_contigs_from_overlapping_reads() {
        let reads = vec![
            "ACGTACGTACGTACGTACGT".to_string(),
            "ACGTACGTACGTGGGGGGGG".to_string(),
            "ACGTGGGGGGGGTTTTTTTT".to_string(),
            "GGGGGGGGTTTTTTTTCCCC".to_string(),
        ];

        let rescued = assemble_read_overlap_contigs(&reads, 10);
        let longest = rescued
            .iter()
            .map(|contig| contig.sequence.len())
            .max()
            .unwrap_or(0);

        assert!(
            longest > reads[0].len(),
            "overlap rescue should extend beyond a single read, got longest={longest}"
        );
    }

    #[test]
    fn read_prefix_index_finds_exact_suffix_prefix_candidates() {
        let reads = vec!["AAACCC", "CCCGGG", "GGGTTT"];
        let index = build_read_prefix_index(&reads, 3, 6);
        assert_eq!(index.get("CCC").cloned().unwrap_or_default(), vec![1]);
        assert_eq!(index.get("GGG").cloned().unwrap_or_default(), vec![2]);
        assert!(index.get("TTT").is_none());
    }

    #[test]
    fn fragmented_kmer_contigs_are_replaced_by_better_read_overlap_rescue() {
        let fragmented = vec![Contig {
            id: 0,
            sequence: "ACGTACGT".to_string(),
            kmer_path: vec![],
        }];
        let reads = vec![
            "ACGTACGTACGTACGTACGT".to_string(),
            "ACGTACGTACGTGGGGGGGG".to_string(),
            "ACGTGGGGGGGGTTTTTTTT".to_string(),
            "GGGGGGGGTTTTTTTTCCCC".to_string(),
        ];

        let rescued = maybe_rescue_fragmented_contigs_with_read_overlaps(fragmented, &reads, 8);
        let longest = rescued
            .iter()
            .map(|contig| contig.sequence.len())
            .max()
            .unwrap_or(0);

        assert!(longest > reads[0].len());
    }

    #[test]
    fn estimate_sequence_capacity_applies_floor_default_and_cap() {
        // Tiny files keep a floor to avoid churn from repeated reallocations.
        assert_eq!(estimate_sequence_capacity(Some(128 * 1024)), 10_000);

        // Missing metadata falls back to a stable default.
        assert_eq!(estimate_sequence_capacity(None), 100_000);

        // Very large files are capped to avoid pathological pre-allocation.
        let two_tb = 2_u64 * 1024 * 1024 * 1024 * 1024;
        assert_eq!(estimate_sequence_capacity(Some(two_tb)), 2_000_000);
    }

    #[test]
    fn summarize_assembly_quality_reports_base_composition_and_length_metrics() {
        let contigs = vec![
            Contig {
                id: 0,
                sequence: "GGCC".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 1,
                sequence: "AANT".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 2,
                sequence: "ARYT".to_string(),
                kmer_path: vec![],
            },
        ];

        let summary = summarize_assembly_quality(&contigs);
        assert_eq!(summary.total_contigs, 3);
        assert_eq!(summary.total_bases, 12);
        assert_eq!(summary.ungapped_total_bases, 11);
        assert_eq!(summary.median_length, 4.0);
        assert_eq!(summary.n10, 4);
        assert_eq!(summary.n25, 4);
        assert_eq!(summary.n50, 4);
        assert_eq!(summary.n75, 4);
        assert_eq!(summary.n90, 4);
        assert_eq!(summary.n95, 4);
        assert_eq!(summary.n99, 4);
        assert_eq!(summary.l10, 1);
        assert_eq!(summary.l25, 1);
        assert_eq!(summary.l50, 2);
        assert_eq!(summary.l75, 3);
        assert_eq!(summary.l90, 3);
        assert_eq!(summary.l95, 3);
        assert_eq!(summary.l99, 3);
        assert_eq!(summary.longest, 4);
        assert!((summary.longest_frac - (1.0 / 3.0)).abs() < 1e-12);
        assert!((summary.avg_length - 4.0).abs() < 1e-12);
        assert!((summary.au_n - 4.0).abs() < 1e-12);
        assert!((summary.effective_contig_count - 3.0).abs() < 1e-12);
        assert_eq!(summary.ungapped_n50, 4);
        assert_eq!(summary.ungapped_n90, 3);
        assert_eq!(summary.ungapped_n95, 3);
        assert_eq!(summary.ungapped_n99, 3);
        assert!((summary.ungapped_au_n - (41.0 / 11.0)).abs() < 1e-12);
        assert!((summary.ungapped_effective_contig_count - (11.0 / (41.0 / 11.0))).abs() < 1e-12);
        assert_eq!(summary.gap_bases, 1);
        assert!((summary.gap_bases_frac - (1.0 / 12.0)).abs() < 1e-12);
        assert_eq!(summary.gc_bases, 4);
        assert_eq!(summary.acgt_bases, 9);
        assert_eq!(summary.n_bases, 1);
        assert_eq!(summary.ambiguous_bases, 2);
        assert_eq!(summary.contigs_with_n, 1);
        assert_eq!(summary.contigs_with_ambiguous, 1);
        assert_eq!(summary.contigs_all_acgt, 1);
        assert_eq!(summary.distinct_sequences, 3);
        assert_eq!(summary.duplicate_sequences, 0);
        assert_eq!(summary.duplicate_bases, 0);
        assert_eq!(summary.canonical_distinct_sequences, 3);
        assert_eq!(summary.canonical_duplicate_sequences, 0);
        assert_eq!(summary.canonical_duplicate_bases, 0);
        assert!((summary.contigs_with_n_frac - (1.0 / 3.0)).abs() < 1e-12);
        assert!((summary.contigs_with_ambiguous_frac - (1.0 / 3.0)).abs() < 1e-12);
        assert!((summary.contigs_all_acgt_frac - (1.0 / 3.0)).abs() < 1e-12);
        assert_eq!(summary.distinct_sequences_frac, 1.0);
        assert_eq!(summary.duplicate_sequences_frac, 0.0);
        assert_eq!(summary.duplicate_bases_frac, 0.0);
        assert_eq!(summary.canonical_distinct_sequences_frac, 1.0);
        assert_eq!(summary.canonical_duplicate_sequences_frac, 0.0);
        assert_eq!(summary.canonical_duplicate_bases_frac, 0.0);
        assert!((summary.mean_rle_ratio - 0.75).abs() < 1e-12);
        assert!((summary.length_weighted_rle_ratio - 0.75).abs() < 1e-12);
        assert_eq!(summary.total_rle_runs, 9);
        assert!((summary.gc_content - (4.0 / 9.0)).abs() < 1e-12);
        assert!((summary.n_content - (1.0 / 12.0)).abs() < 1e-12);
        assert!((summary.ambiguous_content - (2.0 / 12.0)).abs() < 1e-12);
        assert_eq!(summary.n_runs, 1);
        assert_eq!(summary.longest_n_run, 1);
        assert!((summary.mean_n_run_length - 1.0).abs() < 1e-12);
        assert!((summary.n_runs_per_100kb - (1.0 * 100_000.0 / 12.0)).abs() < 1e-12);
        assert!((summary.n_bases_per_100kb - (1.0 * 100_000.0 / 12.0)).abs() < 1e-12);
        assert!((summary.ambiguous_bases_per_100kb - (2.0 * 100_000.0 / 12.0)).abs() < 1e-12);
        assert_eq!(summary.contigs_ge_1kb, 0);
        assert_eq!(summary.contigs_ge_1kb_frac, 0.0);
        assert_eq!(summary.bases_ge_1kb_frac, 0.0);
    }

    #[test]
    fn summarize_assembly_quality_reports_1kb_bucket_fractions() {
        let contigs = vec![
            Contig {
                id: 10,
                sequence: "A".repeat(1_200),
                kmer_path: vec![],
            },
            Contig {
                id: 11,
                sequence: "T".repeat(800),
                kmer_path: vec![],
            },
        ];

        let summary = summarize_assembly_quality(&contigs);
        assert_eq!(summary.contigs_ge_1kb, 1);
        assert_eq!(summary.contigs_ge_10kb, 0);
        assert_eq!(summary.contigs_ge_50kb, 0);
        assert_eq!(summary.contigs_ge_100kb, 0);
        assert_eq!(summary.contigs_ge_1mb, 0);
        assert_eq!(summary.bases_ge_1kb, 1_200);
        assert_eq!(summary.bases_ge_10kb, 0);
        assert_eq!(summary.bases_ge_50kb, 0);
        assert_eq!(summary.bases_ge_100kb, 0);
        assert_eq!(summary.bases_ge_1mb, 0);
        assert!((summary.contigs_ge_1kb_frac - 0.5).abs() < 1e-12);
        assert_eq!(summary.contigs_ge_10kb_frac, 0.0);
        assert_eq!(summary.contigs_ge_50kb_frac, 0.0);
        assert_eq!(summary.contigs_ge_100kb_frac, 0.0);
        assert_eq!(summary.contigs_ge_1mb_frac, 0.0);
        assert!((summary.bases_ge_1kb_frac - 0.6).abs() < 1e-12);
        assert_eq!(summary.bases_ge_10kb_frac, 0.0);
        assert_eq!(summary.bases_ge_50kb_frac, 0.0);
        assert_eq!(summary.bases_ge_100kb_frac, 0.0);
        assert_eq!(summary.bases_ge_1mb_frac, 0.0);
    }

    #[test]
    fn summarize_assembly_quality_reports_multi_threshold_contig_and_span_metrics() {
        let contigs = vec![
            Contig {
                id: 0,
                sequence: "A".repeat(100_000),
                kmer_path: vec![],
            },
            Contig {
                id: 1,
                sequence: "C".repeat(50_000),
                kmer_path: vec![],
            },
            Contig {
                id: 2,
                sequence: "G".repeat(10_000),
                kmer_path: vec![],
            },
            Contig {
                id: 3,
                sequence: "T".repeat(1_000),
                kmer_path: vec![],
            },
            Contig {
                id: 4,
                sequence: "N".repeat(999),
                kmer_path: vec![],
            },
        ];

        let summary = summarize_assembly_quality(&contigs);
        assert_eq!(summary.total_bases, 161_999);
        assert_eq!(summary.ungapped_total_bases, 161_000);
        assert_eq!(summary.median_length, 10_000.0);
        assert_eq!(summary.n10, 100_000);
        assert_eq!(summary.n50, 100_000);
        assert_eq!(summary.n90, 50_000);
        assert_eq!(summary.n95, 10_000);
        assert_eq!(summary.n99, 1_000);
        assert_eq!(summary.l10, 1);
        assert_eq!(summary.l25, 1);
        assert_eq!(summary.l50, 1);
        assert_eq!(summary.l75, 2);
        assert_eq!(summary.l90, 2);
        assert_eq!(summary.l95, 3);
        assert_eq!(summary.l99, 4);
        assert_eq!(summary.contigs_ge_1kb, 4);
        assert_eq!(summary.contigs_ge_10kb, 3);
        assert_eq!(summary.contigs_ge_50kb, 2);
        assert_eq!(summary.contigs_ge_100kb, 1);
        assert_eq!(summary.contigs_ge_1mb, 0);
        assert_eq!(summary.bases_ge_1kb, 161_000);
        assert_eq!(summary.bases_ge_10kb, 160_000);
        assert_eq!(summary.bases_ge_50kb, 150_000);
        assert_eq!(summary.bases_ge_100kb, 100_000);
        assert_eq!(summary.bases_ge_1mb, 0);
        assert!((summary.longest_frac - (100_000.0 / 161_999.0)).abs() < 1e-12);
        assert!((summary.contigs_ge_1kb_frac - 0.8).abs() < 1e-12);
        assert!((summary.contigs_ge_10kb_frac - 0.6).abs() < 1e-12);
        assert!((summary.contigs_ge_50kb_frac - 0.4).abs() < 1e-12);
        assert!((summary.contigs_ge_100kb_frac - 0.2).abs() < 1e-12);
        assert_eq!(summary.contigs_ge_1mb_frac, 0.0);
        assert!((summary.bases_ge_1kb_frac - (161_000.0 / 161_999.0)).abs() < 1e-12);
        assert!((summary.bases_ge_10kb_frac - (160_000.0 / 161_999.0)).abs() < 1e-12);
        assert!((summary.bases_ge_50kb_frac - (150_000.0 / 161_999.0)).abs() < 1e-12);
        assert!((summary.bases_ge_100kb_frac - (100_000.0 / 161_999.0)).abs() < 1e-12);
        assert_eq!(summary.bases_ge_1mb_frac, 0.0);
        assert_eq!(summary.contigs_with_n, 1);
        assert_eq!(summary.contigs_with_ambiguous, 0);
        assert_eq!(summary.contigs_all_acgt, 4);
        assert!((summary.contigs_with_n_frac - 0.2).abs() < 1e-12);
        assert_eq!(summary.contigs_with_ambiguous_frac, 0.0);
        assert!((summary.contigs_all_acgt_frac - 0.8).abs() < 1e-12);
        assert_eq!(summary.n_runs, 1);
        assert_eq!(summary.longest_n_run, 999);
        assert!((summary.mean_n_run_length - 999.0).abs() < 1e-12);
        assert!((summary.mean_rle_ratio - 0.0004262002002002).abs() < 1e-12);
        assert!((summary.length_weighted_rle_ratio - (5.0 / 161_999.0)).abs() < 1e-12);
        assert_eq!(summary.total_rle_runs, 5);
        assert_eq!(summary.ungapped_n50, 100_000);
        assert_eq!(summary.ungapped_n90, 50_000);
        assert_eq!(summary.ungapped_n95, 10_000);
        assert_eq!(summary.ungapped_n99, 10_000);
        assert!((summary.ungapped_au_n - (12_601_000_000.0 / 161_000.0)).abs() < 1e-9);
        assert_eq!(summary.gap_bases, 999);
        assert!((summary.gap_bases_frac - (999.0 / 161_999.0)).abs() < 1e-12);
        assert!((summary.effective_contig_count - (161_999.0 / summary.au_n)).abs() < 1e-9);
    }

    #[test]
    fn summarize_assembly_quality_reports_n_run_metrics() {
        let contigs = vec![
            Contig {
                id: 0,
                sequence: "AANNNCC".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 1,
                sequence: "NNAAAN".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 2,
                sequence: "CGT".to_string(),
                kmer_path: vec![],
            },
        ];

        let summary = summarize_assembly_quality(&contigs);
        assert_eq!(summary.n_bases, 6);
        assert_eq!(summary.ungapped_total_bases, 10);
        assert_eq!(summary.contigs_with_n, 2);
        assert_eq!(summary.contigs_with_ambiguous, 0);
        assert_eq!(summary.contigs_all_acgt, 1);
        assert!((summary.contigs_with_n_frac - (2.0 / 3.0)).abs() < 1e-12);
        assert_eq!(summary.contigs_with_ambiguous_frac, 0.0);
        assert!((summary.contigs_all_acgt_frac - (1.0 / 3.0)).abs() < 1e-12);
        assert_eq!(summary.n_runs, 3);
        assert_eq!(summary.longest_n_run, 3);
        assert!((summary.mean_n_run_length - 2.0).abs() < 1e-12);
        assert!((summary.mean_rle_ratio - ((3.0 / 7.0 + 3.0 / 6.0 + 1.0) / 3.0)).abs() < 1e-12);
        assert!((summary.length_weighted_rle_ratio - (9.0 / 16.0)).abs() < 1e-12);
        assert_eq!(summary.total_rle_runs, 9);
        assert!((summary.n_runs_per_100kb - (3.0 * 100_000.0 / 16.0)).abs() < 1e-12);
        assert_eq!(summary.ungapped_n50, 3);
        assert_eq!(summary.ungapped_n90, 3);
        assert_eq!(summary.ungapped_n95, 3);
        assert_eq!(summary.ungapped_n99, 3);
        assert_eq!(summary.gap_bases, 6);
        assert!((summary.gap_bases_frac - (6.0 / 16.0)).abs() < 1e-12);
    }

    #[test]
    fn summarize_assembly_quality_rle_metrics_are_case_insensitive() {
        let contigs = vec![
            Contig {
                id: 0,
                sequence: "AaAA".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 1,
                sequence: "tT".to_string(),
                kmer_path: vec![],
            },
        ];

        let summary = summarize_assembly_quality(&contigs);
        assert_eq!(summary.total_rle_runs, 2);
        assert!((summary.mean_rle_ratio - ((1.0 / 4.0 + 1.0 / 2.0) / 2.0)).abs() < 1e-12);
        assert!((summary.length_weighted_rle_ratio - (2.0 / 6.0)).abs() < 1e-12);
    }

    #[test]
    fn summarize_assembly_quality_reports_duplicate_sequence_metrics() {
        let contigs = vec![
            Contig {
                id: 0,
                sequence: "ACGT".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 1,
                sequence: "TT".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 2,
                sequence: "ACGT".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 3,
                sequence: "ACGT".to_string(),
                kmer_path: vec![],
            },
        ];

        let summary = summarize_assembly_quality(&contigs);
        assert_eq!(summary.total_contigs, 4);
        assert_eq!(summary.total_bases, 14);
        assert_eq!(summary.distinct_sequences, 2);
        assert_eq!(summary.duplicate_sequences, 2);
        assert_eq!(summary.duplicate_bases, 8);
        assert_eq!(summary.canonical_distinct_sequences, 2);
        assert_eq!(summary.canonical_duplicate_sequences, 2);
        assert_eq!(summary.canonical_duplicate_bases, 8);
        assert!((summary.distinct_sequences_frac - 0.5).abs() < 1e-12);
        assert!((summary.duplicate_sequences_frac - 0.5).abs() < 1e-12);
        assert!((summary.duplicate_bases_frac - (8.0 / 14.0)).abs() < 1e-12);
        assert!((summary.canonical_distinct_sequences_frac - 0.5).abs() < 1e-12);
        assert!((summary.canonical_duplicate_sequences_frac - 0.5).abs() < 1e-12);
        assert!((summary.canonical_duplicate_bases_frac - (8.0 / 14.0)).abs() < 1e-12);
    }

    #[test]
    fn summarize_assembly_quality_reports_canonical_duplicate_sequence_metrics() {
        let contigs = vec![
            Contig {
                id: 0,
                sequence: "AAGT".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 1,
                sequence: "actt".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 2,
                sequence: "AaGt".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 3,
                sequence: "CCCC".to_string(),
                kmer_path: vec![],
            },
        ];

        let summary = summarize_assembly_quality(&contigs);
        assert_eq!(summary.total_contigs, 4);
        assert_eq!(summary.total_bases, 16);
        // Exact-string duplicate metrics remain strict.
        assert_eq!(summary.distinct_sequences, 4);
        assert_eq!(summary.duplicate_sequences, 0);
        assert_eq!(summary.duplicate_bases, 0);
        assert_eq!(summary.distinct_sequences_frac, 1.0);
        assert_eq!(summary.duplicate_sequences_frac, 0.0);
        assert_eq!(summary.duplicate_bases_frac, 0.0);
        // Canonical duplicate metrics collapse case and reverse-complement variants.
        assert_eq!(summary.canonical_distinct_sequences, 2);
        assert_eq!(summary.canonical_duplicate_sequences, 2);
        assert_eq!(summary.canonical_duplicate_bases, 8);
        assert!((summary.canonical_distinct_sequences_frac - 0.5).abs() < 1e-12);
        assert!((summary.canonical_duplicate_sequences_frac - 0.5).abs() < 1e-12);
        assert!((summary.canonical_duplicate_bases_frac - 0.5).abs() < 1e-12);
    }

    #[test]
    fn summarize_assembly_quality_is_invariant_to_contig_order() {
        let base_contigs = vec![
            Contig {
                id: 0,
                sequence: "AANNNCC".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 1,
                sequence: "NNAAAN".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 2,
                sequence: "CGTARY".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 3,
                sequence: "GGCC".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 4,
                sequence: "N".repeat(1000),
                kmer_path: vec![],
            },
        ];

        let expected = summarize_assembly_quality(&base_contigs);
        let mut permuted = base_contigs.clone();
        let mut rng = StdRng::seed_from_u64(0xA55E_4B1E);

        for _ in 0..128 {
            permuted.shuffle(&mut rng);
            let observed = summarize_assembly_quality(&permuted);
            assert_eq!(observed, expected);
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(64))]

        #[test]
        fn summarize_assembly_quality_satisfies_internal_invariants(
            sequences in prop::collection::vec("[ACGTNRYacgtnry]{0,120}", 1..24)
        ) {
            let contigs: Vec<Contig> = sequences
                .iter()
                .enumerate()
                .map(|(id, sequence)| Contig {
                    id,
                    sequence: sequence.clone(),
                    kmer_path: vec![],
                })
                .collect();

            let summary = summarize_assembly_quality(&contigs);
            let lengths: Vec<usize> = sequences.iter().map(String::len).collect();
            let total_bases: usize = lengths.iter().sum();
            let n_bases: usize = sequences
                .iter()
                .map(|sequence| sequence.bytes().filter(|&base| is_n(base)).count())
                .sum();
            let ambiguous_bases: usize = sequences
                .iter()
                .map(|sequence| sequence.bytes().filter(|&base| is_ambiguous(base)).count())
                .sum();
            let ungapped_total_bases = total_bases.saturating_sub(n_bases);
            let gc_bases: usize = sequences
                .iter()
                .map(|sequence| {
                    sequence
                        .bytes()
                        .filter(|&base| matches!(base, b'G' | b'g' | b'C' | b'c'))
                        .count()
                })
                .sum();
            let acgt_bases: usize = sequences
                .iter()
                .map(|sequence| {
                    sequence
                        .bytes()
                        .filter(|&base| matches!(
                            base,
                            b'A' | b'a' | b'T' | b't' | b'U' | b'u' | b'G' | b'g' | b'C' | b'c'
                        ))
                        .count()
                })
                .sum();
            let contigs_with_n = sequences
                .iter()
                .filter(|sequence| sequence.bytes().any(is_n))
                .count();
            let contigs_with_ambiguous = sequences
                .iter()
                .filter(|sequence| sequence.bytes().any(is_ambiguous))
                .count();
            let contigs_all_acgt = sequences
                .iter()
                .filter(|sequence| sequence.bytes().all(|base| is_acgt_or_u(base)))
                .count();
            let mut sequence_counts = std::collections::HashMap::new();
            let mut canonical_sequence_counts = std::collections::HashMap::new();
            for sequence in &sequences {
                *sequence_counts.entry(sequence.as_str()).or_insert(0usize) += 1;
                *canonical_sequence_counts
                    .entry(canonical_sequence_key(sequence))
                    .or_insert(0usize) += 1;
            }
            let distinct_sequences = sequence_counts.len();
            let duplicate_sequences: usize = sequence_counts
                .values()
                .map(|count| count.saturating_sub(1))
                .sum();
            let duplicate_bases: usize = sequence_counts
                .iter()
                .map(|(sequence, count)| sequence.len().saturating_mul(count.saturating_sub(1)))
                .sum();
            let canonical_distinct_sequences = canonical_sequence_counts.len();
            let canonical_duplicate_sequences: usize = canonical_sequence_counts
                .values()
                .map(|count| count.saturating_sub(1))
                .sum();
            let canonical_duplicate_bases: usize = canonical_sequence_counts
                .iter()
                .map(|(sequence, count)| sequence.len().saturating_mul(count.saturating_sub(1)))
                .sum();
            let total_rle_runs: usize = sequences.iter().map(|sequence| count_rle_runs(sequence)).sum();
            let (n_runs, longest_n_run) = sequences
                .iter()
                .fold((0usize, 0usize), |(total_runs, max_run), sequence| {
                    let (runs, longest) = n_run_stats(sequence);
                    (total_runs + runs, max_run.max(longest))
                });

            prop_assert_eq!(summary.total_contigs, sequences.len());
            prop_assert_eq!(summary.total_bases, total_bases);
            prop_assert_eq!(summary.ungapped_total_bases, ungapped_total_bases);
            prop_assert_eq!(summary.gap_bases, total_bases.saturating_sub(ungapped_total_bases));
            prop_assert_eq!(summary.gc_bases, gc_bases);
            prop_assert_eq!(summary.acgt_bases, acgt_bases);
            prop_assert_eq!(summary.n_bases, n_bases);
            prop_assert_eq!(summary.ambiguous_bases, ambiguous_bases);
            prop_assert_eq!(summary.contigs_with_n, contigs_with_n);
            prop_assert_eq!(summary.contigs_with_ambiguous, contigs_with_ambiguous);
            prop_assert_eq!(summary.contigs_all_acgt, contigs_all_acgt);
            prop_assert_eq!(summary.distinct_sequences, distinct_sequences);
            prop_assert_eq!(summary.duplicate_sequences, duplicate_sequences);
            prop_assert_eq!(summary.duplicate_bases, duplicate_bases);
            prop_assert_eq!(summary.canonical_distinct_sequences, canonical_distinct_sequences);
            prop_assert_eq!(summary.canonical_duplicate_sequences, canonical_duplicate_sequences);
            prop_assert_eq!(summary.canonical_duplicate_bases, canonical_duplicate_bases);
            prop_assert_eq!(summary.total_rle_runs, total_rle_runs);
            prop_assert_eq!(summary.n_runs, n_runs);
            prop_assert_eq!(summary.longest_n_run, longest_n_run);
            prop_assert!(summary.ungapped_total_bases <= summary.total_bases);
            prop_assert_eq!(summary.distinct_sequences + summary.duplicate_sequences, summary.total_contigs);
            prop_assert_eq!(
                summary.canonical_distinct_sequences + summary.canonical_duplicate_sequences,
                summary.total_contigs
            );
            prop_assert!(summary.duplicate_bases <= summary.total_bases);
            prop_assert!(summary.canonical_duplicate_bases <= summary.total_bases);
            prop_assert!(summary.canonical_duplicate_sequences >= summary.duplicate_sequences);
            prop_assert!(summary.canonical_distinct_sequences <= summary.distinct_sequences);

            prop_assert!(summary.contigs_ge_1mb <= summary.contigs_ge_100kb);
            prop_assert!(summary.contigs_ge_100kb <= summary.contigs_ge_50kb);
            prop_assert!(summary.contigs_ge_50kb <= summary.contigs_ge_10kb);
            prop_assert!(summary.contigs_ge_10kb <= summary.contigs_ge_1kb);
            prop_assert!(summary.contigs_ge_1kb <= summary.total_contigs);

            prop_assert!(summary.bases_ge_1mb <= summary.bases_ge_100kb);
            prop_assert!(summary.bases_ge_100kb <= summary.bases_ge_50kb);
            prop_assert!(summary.bases_ge_50kb <= summary.bases_ge_10kb);
            prop_assert!(summary.bases_ge_10kb <= summary.bases_ge_1kb);
            prop_assert!(summary.bases_ge_1kb <= summary.total_bases);
            prop_assert!(summary.ungapped_n50 >= summary.ungapped_n90);
            prop_assert!(summary.ungapped_n90 >= summary.ungapped_n95);
            prop_assert!(summary.ungapped_n95 >= summary.ungapped_n99);

            for fraction in [
                summary.contigs_with_n_frac,
                summary.contigs_with_ambiguous_frac,
                summary.contigs_all_acgt_frac,
                summary.distinct_sequences_frac,
                summary.duplicate_sequences_frac,
                summary.canonical_distinct_sequences_frac,
                summary.canonical_duplicate_sequences_frac,
                summary.contigs_ge_1kb_frac,
                summary.contigs_ge_10kb_frac,
                summary.contigs_ge_50kb_frac,
                summary.contigs_ge_100kb_frac,
                summary.contigs_ge_1mb_frac,
                summary.duplicate_bases_frac,
                summary.canonical_duplicate_bases_frac,
                summary.bases_ge_1kb_frac,
                summary.bases_ge_10kb_frac,
                summary.bases_ge_50kb_frac,
                summary.bases_ge_100kb_frac,
                summary.bases_ge_1mb_frac,
                summary.gap_bases_frac,
                summary.gc_content,
                summary.n_content,
                summary.ambiguous_content,
            ] {
                prop_assert!((0.0..=1.0).contains(&fraction));
            }

            if summary.total_bases > 0 {
                let total_bases_f = summary.total_bases as f64;
                prop_assert!((summary.n_content - summary.n_bases as f64 / total_bases_f).abs() < 1e-12);
                prop_assert!((summary.ambiguous_content - summary.ambiguous_bases as f64 / total_bases_f).abs() < 1e-12);
                prop_assert!((summary.gap_bases_frac - summary.gap_bases as f64 / total_bases_f).abs() < 1e-12);
                prop_assert!((summary.duplicate_bases_frac - summary.duplicate_bases as f64 / total_bases_f).abs() < 1e-12);
                prop_assert!((
                    summary.canonical_duplicate_bases_frac
                        - summary.canonical_duplicate_bases as f64 / total_bases_f
                ).abs() < 1e-12);
                prop_assert!((summary.length_weighted_rle_ratio - summary.total_rle_runs as f64 / total_bases_f).abs() < 1e-12);
                prop_assert!((summary.n_runs_per_100kb - summary.n_runs as f64 * 100_000.0 / total_bases_f).abs() < 1e-12);
                prop_assert!((summary.n_bases_per_100kb - summary.n_bases as f64 * 100_000.0 / total_bases_f).abs() < 1e-12);
                prop_assert!((summary.ambiguous_bases_per_100kb - summary.ambiguous_bases as f64 * 100_000.0 / total_bases_f).abs() < 1e-12);
            } else {
                prop_assert_eq!(summary.n_content, 0.0);
                prop_assert_eq!(summary.ambiguous_content, 0.0);
                prop_assert_eq!(summary.gap_bases_frac, 0.0);
                prop_assert_eq!(summary.duplicate_bases_frac, 0.0);
                prop_assert_eq!(summary.canonical_duplicate_bases_frac, 0.0);
                prop_assert_eq!(summary.length_weighted_rle_ratio, 0.0);
                prop_assert_eq!(summary.n_runs_per_100kb, 0.0);
                prop_assert_eq!(summary.n_bases_per_100kb, 0.0);
                prop_assert_eq!(summary.ambiguous_bases_per_100kb, 0.0);
            }
        }
    }

    #[test]
    fn assemble_reads_with_gpu_rejects_truncated_fastq_input() {
        let mut input = NamedTempFile::new().expect("create input fastq");
        writeln!(input, "@read_1").expect("write header");
        writeln!(input, "ACGTACGT").expect("write sequence");
        writeln!(input, "+").expect("write plus line");
        input.flush().expect("flush input");

        let output = NamedTempFile::new().expect("create output path");
        let err = assemble_reads_with_gpu(
            input.path().to_str().expect("utf8 input path"),
            None,
            output.path().to_str().expect("utf8 output path"),
            1,     // min_len
            false, // output_gfa
            false, // output_gfa2
            false, // adaptive_k
            false, // use_rle
            false, // collapse_repeats
            0,     // min_repeat_len
            false, // polish
            21,    // polish_window
            false, // streaming
            false, // export_metadata
            None,  // json_metadata
            None,  // tsv_metadata
            false, // isoforms
            None,  // gtf_path
            None,  // gff3_path
            100,   // max_path_depth
            0.0,   // min_confidence
            false, // compute_tpm
            false, // polish_isoforms
            None,  // samples_path
            0.0,   // min_tpm
            None,  // long_reads
            false, // counts_matrix
            false, // use_gpu
        )
        .expect_err("truncated FASTQ must return an error");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn assemble_reads_with_gpu_rejects_mismatched_paired_fastq_inputs() {
        let mut r1 = NamedTempFile::new().expect("create r1 fastq");
        writeln!(r1, "@read_1/1").expect("write header");
        writeln!(r1, "ACGTACGT").expect("write sequence");
        writeln!(r1, "+").expect("write plus line");
        writeln!(r1, "IIIIIIII").expect("write quality");
        writeln!(r1, "@read_2/1").expect("write header");
        writeln!(r1, "TGCATGCA").expect("write sequence");
        writeln!(r1, "+").expect("write plus line");
        writeln!(r1, "IIIIIIII").expect("write quality");
        r1.flush().expect("flush r1");

        let mut r2 = NamedTempFile::new().expect("create r2 fastq");
        writeln!(r2, "@read_1/2").expect("write header");
        writeln!(r2, "ACGTACGT").expect("write sequence");
        writeln!(r2, "+").expect("write plus line");
        writeln!(r2, "IIIIIIII").expect("write quality");
        r2.flush().expect("flush r2");

        let output = NamedTempFile::new().expect("create output path");
        let err = assemble_reads_with_gpu(
            r1.path().to_str().expect("utf8 r1 path"),
            Some(r2.path().to_str().expect("utf8 r2 path")),
            output.path().to_str().expect("utf8 output path"),
            1,     // min_len
            false, // output_gfa
            false, // output_gfa2
            false, // adaptive_k
            false, // use_rle
            false, // collapse_repeats
            0,     // min_repeat_len
            false, // polish
            21,    // polish_window
            false, // streaming
            false, // export_metadata
            None,  // json_metadata
            None,  // tsv_metadata
            false, // isoforms
            None,  // gtf_path
            None,  // gff3_path
            100,   // max_path_depth
            0.0,   // min_confidence
            false, // compute_tpm
            false, // polish_isoforms
            None,  // samples_path
            0.0,   // min_tpm
            None,  // long_reads
            false, // counts_matrix
            false, // use_gpu
        )
        .expect_err("mismatched paired FASTQ inputs must return an error");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn assemble_reads_with_gpu_rejects_directory_output_path() {
        let mut input = NamedTempFile::new().expect("create input fastq");
        writeln!(input, "@read_1").expect("write header");
        writeln!(input, "ACGTACGT").expect("write sequence");
        writeln!(input, "+").expect("write plus line");
        writeln!(input, "IIIIIIII").expect("write quality");
        input.flush().expect("flush input");

        let temp_dir = TempDir::new().expect("create temp dir");
        let output_dir = temp_dir.path().join("out_dir");
        std::fs::create_dir_all(&output_dir).expect("create output directory");

        let err = assemble_reads_with_gpu(
            input.path().to_str().expect("utf8 input path"),
            None,
            output_dir.to_str().expect("utf8 output path"),
            1,     // min_len
            false, // output_gfa
            false, // output_gfa2
            false, // adaptive_k
            false, // use_rle
            false, // collapse_repeats
            0,     // min_repeat_len
            false, // polish
            21,    // polish_window
            false, // streaming
            false, // export_metadata
            None,  // json_metadata
            None,  // tsv_metadata
            false, // isoforms
            None,  // gtf_path
            None,  // gff3_path
            100,   // max_path_depth
            0.0,   // min_confidence
            false, // compute_tpm
            false, // polish_isoforms
            None,  // samples_path
            0.0,   // min_tpm
            None,  // long_reads
            false, // counts_matrix
            false, // use_gpu
        )
        .expect_err("directory output path must return io::Error");
        assert_eq!(err.kind(), io::ErrorKind::IsADirectory);
    }

    #[test]
    fn write_assembly_quality_reports_emits_stable_json_and_tsv_sidecars() {
        let temp_dir = TempDir::new().expect("create temp dir");
        let output_base = temp_dir.path().join("assembled.fasta");
        let summary = AssemblyQualitySummary {
            total_contigs: 5,
            total_bases: 1234,
            ungapped_total_bases: 1200,
            avg_length: 246.8,
            median_length: 210.0,
            n10: 320,
            n25: 300,
            n50: 250,
            n75: 200,
            n90: 150,
            n95: 140,
            n99: 100,
            l10: 1,
            l25: 1,
            l50: 2,
            l75: 3,
            l90: 4,
            l95: 5,
            l99: 5,
            au_n: 260.5,
            effective_contig_count: 1234.0 / 260.5,
            ungapped_n50: 240,
            ungapped_n90: 180,
            ungapped_n95: 170,
            ungapped_n99: 150,
            ungapped_au_n: 255.25,
            ungapped_effective_contig_count: 1200.0 / 255.25,
            gap_bases: 34,
            gap_bases_frac: 34.0 / 1234.0,
            longest: 420,
            longest_frac: 0.340356564,
            gc_bases: 580,
            acgt_bases: 1100,
            n_bases: 90,
            ambiguous_bases: 44,
            contigs_with_n: 2,
            contigs_with_ambiguous: 1,
            contigs_all_acgt: 2,
            distinct_sequences: 4,
            duplicate_sequences: 1,
            duplicate_bases: 120,
            canonical_distinct_sequences: 3,
            canonical_duplicate_sequences: 2,
            canonical_duplicate_bases: 240,
            contigs_with_n_frac: 0.4,
            contigs_with_ambiguous_frac: 0.2,
            contigs_all_acgt_frac: 0.4,
            distinct_sequences_frac: 0.8,
            duplicate_sequences_frac: 0.2,
            duplicate_bases_frac: 120.0 / 1234.0,
            canonical_distinct_sequences_frac: 0.6,
            canonical_duplicate_sequences_frac: 0.4,
            canonical_duplicate_bases_frac: 240.0 / 1234.0,
            mean_rle_ratio: 0.75,
            length_weighted_rle_ratio: 0.8,
            total_rle_runs: 987,
            gc_content: 0.5,
            n_content: 0.1,
            ambiguous_content: 0.02,
            n_runs: 4,
            longest_n_run: 21,
            mean_n_run_length: 22.5,
            n_runs_per_100kb: 4.0 * 100_000.0 / 1234.0,
            n_bases_per_100kb: 90.0 * 100_000.0 / 1234.0,
            ambiguous_bases_per_100kb: 44.0 * 100_000.0 / 1234.0,
            contigs_ge_1kb: 1,
            contigs_ge_10kb: 0,
            contigs_ge_50kb: 0,
            contigs_ge_100kb: 0,
            contigs_ge_1mb: 0,
            bases_ge_1kb: 1000,
            bases_ge_10kb: 0,
            bases_ge_50kb: 0,
            bases_ge_100kb: 0,
            bases_ge_1mb: 0,
            contigs_ge_1kb_frac: 0.2,
            contigs_ge_10kb_frac: 0.0,
            contigs_ge_50kb_frac: 0.0,
            contigs_ge_100kb_frac: 0.0,
            contigs_ge_1mb_frac: 0.0,
            bases_ge_1kb_frac: 0.81,
            bases_ge_10kb_frac: 0.0,
            bases_ge_50kb_frac: 0.0,
            bases_ge_100kb_frac: 0.0,
            bases_ge_1mb_frac: 0.0,
        };

        let (json_path, tsv_path) = write_assembly_quality_reports(
            output_base.to_str().expect("utf8 output path"),
            summary,
        )
        .expect("write quality reports");

        assert!(json_path.ends_with("assembled.assembly_metrics.json"));
        assert!(tsv_path.ends_with("assembled.assembly_metrics.tsv"));

        let json = std::fs::read_to_string(&json_path).expect("read json report");
        let parsed: serde_json::Value = serde_json::from_str(&json).expect("parse json");
        assert_eq!(parsed["total_contigs"], 5);
        assert_eq!(parsed["n10"], 320);
        assert_eq!(parsed["n50"], 250);
        assert_eq!(parsed["ungapped_total_bases"], 1200);
        assert_eq!(parsed["ungapped_n50"], 240);
        assert_eq!(parsed["ungapped_n90"], 180);
        assert_eq!(parsed["ungapped_n95"], 170);
        assert_eq!(parsed["ungapped_n99"], 150);
        assert_eq!(parsed["ungapped_au_n"], 255.25);
        assert_eq!(parsed["effective_contig_count"], 1234.0 / 260.5);
        assert_eq!(parsed["ungapped_effective_contig_count"], 1200.0 / 255.25);
        assert_eq!(parsed["gap_bases"], 34);
        assert_eq!(parsed["gap_bases_frac"], 34.0 / 1234.0);
        assert_eq!(parsed["l10"], 1);
        assert_eq!(parsed["l75"], 3);
        assert_eq!(parsed["acgt_bases"], 1100);
        assert_eq!(parsed["gc_content"], 0.5);
        assert_eq!(parsed["longest_frac"], 0.340356564);
        assert_eq!(parsed["mean_rle_ratio"], 0.75);
        assert_eq!(parsed["length_weighted_rle_ratio"], 0.8);
        assert_eq!(parsed["total_rle_runs"], 987);
        assert_eq!(parsed["n_runs"], 4);
        assert_eq!(parsed["longest_n_run"], 21);
        assert_eq!(parsed["mean_n_run_length"], 22.5);
        assert_eq!(parsed["n_runs_per_100kb"], 4.0 * 100_000.0 / 1234.0);
        assert_eq!(parsed["n_bases_per_100kb"], 90.0 * 100_000.0 / 1234.0);
        assert_eq!(
            parsed["ambiguous_bases_per_100kb"],
            44.0 * 100_000.0 / 1234.0
        );
        assert_eq!(parsed["contigs_with_n"], 2);
        assert_eq!(parsed["contigs_with_ambiguous"], 1);
        assert_eq!(parsed["contigs_all_acgt"], 2);
        assert_eq!(parsed["distinct_sequences"], 4);
        assert_eq!(parsed["duplicate_sequences"], 1);
        assert_eq!(parsed["duplicate_bases"], 120);
        assert_eq!(parsed["canonical_distinct_sequences"], 3);
        assert_eq!(parsed["canonical_duplicate_sequences"], 2);
        assert_eq!(parsed["canonical_duplicate_bases"], 240);
        assert_eq!(parsed["contigs_with_n_frac"], 0.4);
        assert_eq!(parsed["contigs_with_ambiguous_frac"], 0.2);
        assert_eq!(parsed["contigs_all_acgt_frac"], 0.4);
        assert_eq!(parsed["distinct_sequences_frac"], 0.8);
        assert_eq!(parsed["duplicate_sequences_frac"], 0.2);
        assert_eq!(parsed["duplicate_bases_frac"], 120.0 / 1234.0);
        assert_eq!(parsed["canonical_distinct_sequences_frac"], 0.6);
        assert_eq!(parsed["canonical_duplicate_sequences_frac"], 0.4);
        assert_eq!(parsed["canonical_duplicate_bases_frac"], 240.0 / 1234.0);

        let tsv = std::fs::read_to_string(&tsv_path).expect("read tsv report");
        let mut lines = tsv.lines();
        assert_eq!(lines.next(), Some("metric\tvalue"));
        assert!(tsv.contains("n10\t320"));
        assert!(tsv.contains("n50\t250"));
        assert!(tsv.contains("ungapped_total_bases\t1200"));
        assert!(tsv.contains("ungapped_n50\t240"));
        assert!(tsv.contains("ungapped_n90\t180"));
        assert!(tsv.contains("ungapped_n95\t170"));
        assert!(tsv.contains("ungapped_n99\t150"));
        assert!(tsv.contains("ungapped_au_n\t255.250000000000"));
        assert!(tsv.contains("effective_contig_count\t4.737044145873"));
        assert!(tsv.contains("ungapped_effective_contig_count\t4.701273261508"));
        assert!(tsv.contains("gap_bases\t34"));
        assert!(tsv.contains("gap_bases_frac\t0.027552674230"));
        assert!(tsv.contains("l10\t1"));
        assert!(tsv.contains("l75\t3"));
        assert!(tsv.contains("acgt_bases\t1100"));
        assert!(tsv.contains("longest_frac\t0.340356564000"));
        assert!(tsv.contains("gc_content\t0.500000000000"));
        assert!(tsv.contains("mean_rle_ratio\t0.750000000000"));
        assert!(tsv.contains("length_weighted_rle_ratio\t0.800000000000"));
        assert!(tsv.contains("total_rle_runs\t987"));
        assert!(tsv.contains("n_runs\t4"));
        assert!(tsv.contains("longest_n_run\t21"));
        assert!(tsv.contains("mean_n_run_length\t22.500000000000"));
        assert!(tsv.contains("n_runs_per_100kb\t324.149108589951"));
        assert!(tsv.contains("n_bases_per_100kb\t7293.354943273906"));
        assert!(tsv.contains("ambiguous_bases_per_100kb\t3565.640194489465"));
        assert!(tsv.contains("contigs_with_n\t2"));
        assert!(tsv.contains("contigs_with_ambiguous\t1"));
        assert!(tsv.contains("contigs_all_acgt\t2"));
        assert!(tsv.contains("distinct_sequences\t4"));
        assert!(tsv.contains("duplicate_sequences\t1"));
        assert!(tsv.contains("duplicate_bases\t120"));
        assert!(tsv.contains("canonical_distinct_sequences\t3"));
        assert!(tsv.contains("canonical_duplicate_sequences\t2"));
        assert!(tsv.contains("canonical_duplicate_bases\t240"));
        assert!(tsv.contains("contigs_with_n_frac\t0.400000000000"));
        assert!(tsv.contains("contigs_with_ambiguous_frac\t0.200000000000"));
        assert!(tsv.contains("contigs_all_acgt_frac\t0.400000000000"));
        assert!(tsv.contains("distinct_sequences_frac\t0.800000000000"));
        assert!(tsv.contains("duplicate_sequences_frac\t0.200000000000"));
        assert!(tsv.contains("duplicate_bases_frac\t0.097244732577"));
        assert!(tsv.contains("canonical_distinct_sequences_frac\t0.600000000000"));
        assert!(tsv.contains("canonical_duplicate_sequences_frac\t0.400000000000"));
        assert!(tsv.contains("canonical_duplicate_bases_frac\t0.194489465154"));
        assert!(tsv.contains("bases_ge_1kb_frac\t0.810000000000"));
        assert!(tsv.contains("contigs_ge_1mb\t0"));
        assert!(tsv.contains("bases_ge_1mb\t0"));
    }
}
