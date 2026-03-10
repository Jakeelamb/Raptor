use crate::accel::CpuBackend;
use crate::eval::metrics::{evaluate_lengths_in_place, normalized_run_base, BaseComposition};
use crate::graph::assembler::{greedy_assembly_u64, Contig};
use crate::graph::overlap::find_overlaps;
use crate::graph::stitch::OverlapGraphBuilder;
use crate::io::fasta::FastaWriter;
use crate::io::fastq::{stream_fastq_records_checked, try_open_fastq, FastqRecord};
use crate::io::gfa::GfaWriter;
use crate::io::gfa2::Gfa2Writer;
use crate::kmer::variable_k::{kmer_coverage_histogram, optimal_k, select_best_k};
use serde::Serialize;
use std::fs;
use std::io::{self, Write};
use tracing::{info, warn};

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
    ungapped_au_n: f64,
    ungapped_effective_contig_count: f64,
    longest: usize,
    longest_frac: f64,
    gc_bases: usize,
    acgt_bases: usize,
    n_bases: usize,
    ambiguous_bases: usize,
    contigs_with_n: usize,
    contigs_with_ambiguous: usize,
    contigs_all_acgt: usize,
    contigs_with_n_frac: f64,
    contigs_with_ambiguous_frac: f64,
    contigs_all_acgt_frac: f64,
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
fn summarize_assembly_quality(contigs: &[Contig]) -> AssemblyQualitySummary {
    let mut lengths = Vec::with_capacity(contigs.len());
    let mut ungapped_lengths = Vec::with_capacity(contigs.len());
    let mut composition = BaseComposition::default();
    let mut n_runs = NRunSummary::default();
    let mut contig_quality = ContigQualitySummary::default();
    let mut total_rle_ratio_scaled = 0u128;
    let mut total_rle_runs = 0usize;

    for contig in contigs {
        let contig_len = contig.sequence.len();
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
        ungapped_au_n: ungapped_length_stats.au_n,
        ungapped_effective_contig_count: ungapped_length_stats.effective_count,
        longest: length_stats.longest,
        longest_frac,
        gc_bases: composition.gc_bases,
        acgt_bases: composition.acgt_bases,
        n_bases: composition.n_bases,
        ambiguous_bases: composition.ambiguous_bases,
        contigs_with_n: contig_quality.with_n,
        contigs_with_ambiguous: contig_quality.with_ambiguous,
        contigs_all_acgt: contig_quality.all_acgt,
        contigs_with_n_frac,
        contigs_with_ambiguous_frac,
        contigs_all_acgt_frac,
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
        ("ungapped_au_n", format!("{:.12}", quality.ungapped_au_n)),
        (
            "ungapped_effective_contig_count",
            format!("{:.12}", quality.ungapped_effective_contig_count),
        ),
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

    // Determine k-mer size - either adaptive or fixed optimal
    let max_k = 41;
    let mut k: usize;

    let mut records: Vec<FastqRecord> = Vec::new();

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
    // First pass: collect sequences for k optimization (sample if large)
    for record in stream_fastq_records_checked(reader) {
        sequences.push(record?.sequence);
    }

    let num_sequences = sequences.len();
    info!("Loaded {} sequences", num_sequences);

    if use_gpu {
        info!("GPU requested; current assembly path uses CPU k-mer/graph kernels");
    }
    info!("Using CPU backend for k-mer counting and graph build");

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
    let kmer_counts_u64 = cpu_backend.count_kmers_u64_filtered(&sequences, k, min_kmer_count);
    info!(
        "Found {} unique k-mers (after filtering singletons)",
        kmer_counts_u64.len()
    );

    // Sequence strings are no longer needed after k-mer counting.
    // Releasing this buffer early reduces peak RSS during graph cleanup/polishing.
    let released_records = sequences.len();
    sequences.clear();
    sequences.shrink_to_fit();
    info!(
        "Released {} input sequences from memory after k-mer counting",
        released_records
    );

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

    // Keep records for polishing if needed
    if polish || (isoforms && (polish_isoforms || compute_tpm)) {
        let reader = try_open_fastq(input_path)?;
        records = stream_fastq_records_checked(reader)
            .map(|record| record.map(sequence_only_record))
            .collect::<io::Result<Vec<_>>>()?;
        info!(
            "Loaded {} reads in sequence-only form for polishing/TPM",
            records.len()
        );
    }

    // Perform greedy assembly using u64 k-mers
    info!("Assembling contigs with minimum length: {}", min_len);
    let mut contigs = greedy_assembly_u64(k, &kmer_counts_u64, &adjacency, min_len);

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
        "Contig statistics: {} contigs, {} bp total, Mean/Median: {:.1}/{:.1} bp, N10/N25/N50/N75/N90/N95/N99: {}/{}/{}/{}/{}/{}/{} bp, L10/L25/L50/L75/L90/L95/L99: {}/{}/{}/{}/{}/{}/{}, auN/effective count: {:.1}/{:.2}, Ungapped bases/N50/auN/effective count: {}/{}/{:.1}/{:.2}, Longest: {} bp ({:.2}%), GC/N/Ambiguous bases: {}/{}/{} (fractions {:.2}%/{:.2}%/{:.2}%), Contigs with N/Ambiguous/All-ACGT: {}/{}/{} ({:.2}%/{:.2}%/{:.2}%), RLE mean/weighted runs ratio: {:.4}/{:.4} ({} runs), N-runs: count {}, max {}, mean {:.1} bp ({:.1} per 100kb), N/Ambiguous bases per 100kb: {:.1}/{:.1}, >=1kb/10kb/50kb/100kb/1mb contigs: {}/{}/{}/{}/{} ({:.1}%/{:.1}%/{:.1}%/{:.1}%/{:.1}%), span: {}/{}/{}/{}/{} bp ({:.1}%/{:.1}%/{:.1}%/{:.1}%/{:.1}%)",
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
        quality.ungapped_n50,
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
    let mut writer = FastaWriter::new(output_path);
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
            let mut gfa_writer = GfaWriter::new(&gfa_path);

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
            let mut gfa2_writer = Gfa2Writer::new(&gfa2_path);
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
        assemble_reads_with_gpu, canonicalize_contig_output_order, derive_contig_expression_map,
        estimate_sequence_capacity, sequence_only_record, summarize_assembly_quality,
        write_assembly_quality_reports, AssemblyQualitySummary,
    };
    use crate::graph::assembler::Contig;
    use crate::io::fastq::FastqRecord;
    use ahash::AHashMap;
    use rand::rngs::StdRng;
    use rand::seq::SliceRandom;
    use rand::SeedableRng;
    use std::io::{self, Write};
    use tempfile::{NamedTempFile, TempDir};

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
        assert!((summary.ungapped_au_n - (41.0 / 11.0)).abs() < 1e-12);
        assert!((summary.ungapped_effective_contig_count - (11.0 / (41.0 / 11.0))).abs() < 1e-12);
        assert_eq!(summary.gc_bases, 4);
        assert_eq!(summary.acgt_bases, 9);
        assert_eq!(summary.n_bases, 1);
        assert_eq!(summary.ambiguous_bases, 2);
        assert_eq!(summary.contigs_with_n, 1);
        assert_eq!(summary.contigs_with_ambiguous, 1);
        assert_eq!(summary.contigs_all_acgt, 1);
        assert!((summary.contigs_with_n_frac - (1.0 / 3.0)).abs() < 1e-12);
        assert!((summary.contigs_with_ambiguous_frac - (1.0 / 3.0)).abs() < 1e-12);
        assert!((summary.contigs_all_acgt_frac - (1.0 / 3.0)).abs() < 1e-12);
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
        assert!((summary.ungapped_au_n - (12_601_000_000.0 / 161_000.0)).abs() < 1e-9);
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
            ungapped_au_n: 255.25,
            ungapped_effective_contig_count: 1200.0 / 255.25,
            longest: 420,
            longest_frac: 0.340356564,
            gc_bases: 580,
            acgt_bases: 1100,
            n_bases: 90,
            ambiguous_bases: 44,
            contigs_with_n: 2,
            contigs_with_ambiguous: 1,
            contigs_all_acgt: 2,
            contigs_with_n_frac: 0.4,
            contigs_with_ambiguous_frac: 0.2,
            contigs_all_acgt_frac: 0.4,
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
        assert_eq!(parsed["ungapped_au_n"], 255.25);
        assert_eq!(parsed["effective_contig_count"], 1234.0 / 260.5);
        assert_eq!(parsed["ungapped_effective_contig_count"], 1200.0 / 255.25);
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
        assert_eq!(parsed["contigs_with_n_frac"], 0.4);
        assert_eq!(parsed["contigs_with_ambiguous_frac"], 0.2);
        assert_eq!(parsed["contigs_all_acgt_frac"], 0.4);

        let tsv = std::fs::read_to_string(&tsv_path).expect("read tsv report");
        let mut lines = tsv.lines();
        assert_eq!(lines.next(), Some("metric\tvalue"));
        assert!(tsv.contains("n10\t320"));
        assert!(tsv.contains("n50\t250"));
        assert!(tsv.contains("ungapped_total_bases\t1200"));
        assert!(tsv.contains("ungapped_n50\t240"));
        assert!(tsv.contains("ungapped_au_n\t255.250000000000"));
        assert!(tsv.contains("effective_contig_count\t4.737044145873"));
        assert!(tsv.contains("ungapped_effective_contig_count\t4.701273261508"));
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
        assert!(tsv.contains("contigs_with_n_frac\t0.400000000000"));
        assert!(tsv.contains("contigs_with_ambiguous_frac\t0.200000000000"));
        assert!(tsv.contains("contigs_all_acgt_frac\t0.400000000000"));
        assert!(tsv.contains("bases_ge_1kb_frac\t0.810000000000"));
        assert!(tsv.contains("contigs_ge_1mb\t0"));
        assert!(tsv.contains("bases_ge_1mb\t0"));
    }
}
