use crate::accel::CpuBackend;
use crate::eval::metrics::{evaluate_lengths_in_place, BaseComposition};
use crate::graph::assembler::{greedy_assembly_u64, Contig};
use crate::graph::overlap::find_overlaps;
use crate::graph::stitch::OverlapGraphBuilder;
use crate::io::fasta::FastaWriter;
use crate::io::fastq::{open_fastq, stream_fastq_records, FastqRecord};
use crate::io::gfa::GfaWriter;
use crate::io::gfa2::Gfa2Writer;
use crate::kmer::variable_k::{kmer_coverage_histogram, optimal_k, select_best_k};
use std::fs;
use std::io::{self, Write};
use tracing::{info, warn};

#[derive(Debug, Clone, Copy)]
struct AssemblyQualitySummary {
    total_contigs: usize,
    total_bases: usize,
    avg_length: f64,
    n50: usize,
    n90: usize,
    n95: usize,
    au_n: f64,
    longest: usize,
    gc_content: f64,
    n_content: f64,
    ambiguous_content: f64,
    contigs_ge_1kb: usize,
    contigs_ge_1kb_frac: f64,
    bases_ge_1kb_frac: f64,
}

#[inline]
fn summarize_assembly_quality(contigs: &[Contig]) -> AssemblyQualitySummary {
    let mut lengths = Vec::with_capacity(contigs.len());
    let mut composition = BaseComposition::default();

    for contig in contigs {
        lengths.push(contig.sequence.len());
        composition.add_sequence(contig.sequence.as_bytes());
    }

    let length_stats = evaluate_lengths_in_place(&mut lengths);
    AssemblyQualitySummary {
        total_contigs: length_stats.total,
        total_bases: length_stats.total_bases,
        avg_length: length_stats.avg_length,
        n50: length_stats.n50,
        n90: length_stats.n90,
        n95: length_stats.n95,
        au_n: length_stats.au_n,
        longest: length_stats.longest,
        gc_content: composition.gc_content(),
        n_content: composition.n_content(length_stats.total_bases),
        ambiguous_content: composition.ambiguous_content(length_stats.total_bases),
        contigs_ge_1kb: length_stats.contigs_ge_1kb,
        contigs_ge_1kb_frac: length_stats.contigs_ge_1kb_frac,
        bases_ge_1kb_frac: length_stats.bases_ge_1kb_frac,
    }
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
    let reader = open_fastq(input_path);

    // Pre-size sequences vector based on estimated count
    // Estimate: assume average 4 lines per record, ~200 bytes per line for typical FASTQ
    // This reduces reallocation overhead by 10-15%
    const ESTIMATED_RECORDS_PER_MB: usize = 5000;
    let estimated_count = std::fs::metadata(input_path)
        .map(|m| (m.len() as usize / (1024 * 1024)) * ESTIMATED_RECORDS_PER_MB)
        .unwrap_or(100_000)
        .max(10_000); // At least 10k capacity

    let mut sequences: Vec<String> = Vec::with_capacity(estimated_count);
    // First pass: collect sequences for k optimization (sample if large)
    for record in stream_fastq_records(reader) {
        sequences.push(record.sequence);
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
        let reader = open_fastq(input_path);
        records = stream_fastq_records(reader).collect();
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

    let quality = summarize_assembly_quality(&contigs);
    info!(
        "Contig statistics: {} contigs, {} bp total, Avg: {:.1} bp, N50/N90/N95: {}/{}/{} bp, auN: {:.1}, Longest: {} bp, GC: {:.2}%, N: {:.2}%, Ambiguous: {:.2}%, >=1kb: {} ({:.1}% contigs, {:.1}% bases)",
        quality.total_contigs,
        quality.total_bases,
        quality.avg_length,
        quality.n50,
        quality.n90,
        quality.n95,
        quality.au_n,
        quality.longest,
        quality.gc_content * 100.0,
        quality.n_content * 100.0,
        quality.ambiguous_content * 100.0,
        quality.contigs_ge_1kb,
        quality.contigs_ge_1kb_frac * 100.0,
        quality.bases_ge_1kb_frac * 100.0
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
        // Extract raw sequences for overlap finding
        let contig_seqs: Vec<String> = contigs.iter().map(|c| c.sequence.clone()).collect();

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

            // Default values for new parameters
            let min_tpm_value = None; // No TPM filtering by default
            let strand_aware_value = false; // Non-strand-aware by default
            let bam_path_value = None; // No BAM file for quantification
            let long_reads_value = None; // No long reads for polishing

            // Derive per-contig expression support from assembled k-mer paths.
            let kmer_counts_converted = derive_contig_expression_map(&contigs, &kmer_counts_u64);

            // Create a function pointer with the expected signature
            let get_output_filename_fn: fn(&str, Option<&str>) -> String =
                |base: &str, _ext: Option<&str>| get_output_filename(base, "");

            let mut transcripts = crate::pipeline::isoform_processor::process_isoforms(
                &contigs,
                &links,
                &kmer_counts_converted,
                output_path,
                gtf_path.as_deref(),
                max_path_depth,
                min_confidence,
                min_tpm_value,
                strand_aware_value,
                bam_path_value,
                long_reads_value,
                get_output_filename_fn,
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
                "Transcript statistics: {} transcripts, {} bp total, Avg: {:.1} bp, N25/N50/N75/N90/N95/N99: {}/{}/{}/{}/{}/{} bp, L25/L50/L75/L90/L95/L99: {}/{}/{}/{}/{}/{}, auN: {:.1}",
                stats.total,
                stats.total_bases,
                stats.avg_length,
                stats.n25,
                stats.n50,
                stats.n75,
                stats.n90,
                stats.n95,
                stats.n99,
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

#[cfg(test)]
mod tests {
    use super::{derive_contig_expression_map, summarize_assembly_quality};
    use crate::graph::assembler::Contig;
    use ahash::AHashMap;

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
        assert_eq!(summary.n50, 4);
        assert_eq!(summary.n90, 4);
        assert_eq!(summary.n95, 4);
        assert_eq!(summary.longest, 4);
        assert!((summary.avg_length - 4.0).abs() < 1e-12);
        assert!((summary.au_n - 4.0).abs() < 1e-12);
        assert!((summary.gc_content - (4.0 / 9.0)).abs() < 1e-12);
        assert!((summary.n_content - (1.0 / 12.0)).abs() < 1e-12);
        assert!((summary.ambiguous_content - (2.0 / 12.0)).abs() < 1e-12);
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
        assert!((summary.contigs_ge_1kb_frac - 0.5).abs() < 1e-12);
        assert!((summary.bases_ge_1kb_frac - 0.6).abs() < 1e-12);
    }
}
