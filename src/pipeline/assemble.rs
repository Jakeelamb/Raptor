use crate::io::fastq::{open_fastq, stream_fastq_records, FastqRecord};
use crate::io::fasta::FastaWriter;
use crate::io::gfa::GfaWriter;
use crate::io::gfa2::Gfa2Writer;
use crate::kmer::variable_k::{optimal_k, kmer_coverage_histogram, select_best_k};
use crate::graph::assembler::greedy_assembly_u64;
use crate::graph::overlap::find_overlaps;
use crate::graph::stitch::OverlapGraphBuilder;
use crate::accel::{create_backend, CpuBackend};
use std::collections::HashMap;
use tracing::{info, warn};
use std::fs;
use std::io::Write;
use crate::eval::metrics::evaluate_lengths;

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
) {
    // Default to CPU backend for backwards compatibility
    assemble_reads_with_gpu(
        input_path, output_path, min_len, _output_gfa, _output_gfa2,
        _adaptive_k, _use_rle, collapse_repeats, min_repeat_len,
        polish, polish_window, streaming, export_metadata,
        json_metadata, tsv_metadata, isoforms, gtf_path, gff3_path,
        max_path_depth, min_confidence, compute_tpm, polish_isoforms,
        samples_path, min_tpm, long_reads, counts_matrix,
        false  // use_gpu = false by default
    );
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
) {
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
    let mut num_sequences = 0;

    // Chunk size for streaming processing (adjust based on available memory)
    const CHUNK_SIZE: usize = 100_000;

    // First pass: collect sequences for k optimization (sample if large)
    for record in stream_fastq_records(reader) {
        sequences.push(record.sequence);
        num_sequences += 1;
    }

    info!("Loaded {} sequences", num_sequences);

    // Create the appropriate compute backend
    let backend = create_backend(use_gpu, num_sequences, num_sequences / 10 + 100);
    info!("Using compute backend: {}", backend.name());

    // Determine k-mer size from sample
    let sample_size = sequences.len().min(10_000);
    let sample: Vec<String> = sequences.iter().take(sample_size).cloned().collect();

    if _adaptive_k {
        let hist = kmer_coverage_histogram(&sample, max_k);
        k = select_best_k(&hist);
        info!("Using adaptive k-mer: {}", k);
    } else {
        k = optimal_k(&sample, max_k);
        info!("Using optimal k-mer size: {}", k);
    }

    // Ensure k <= 32 for u64 encoding
    if k > 32 {
        info!("Capping k-mer size at 32 for u64 encoding (was {})", k);
        k = 32;
    }

    // Count k-mers using optimized u64 path with Bloom filter pre-filtering
    // This filters singleton k-mers (sequencing errors) for 30-50% memory reduction
    info!("Counting k-mers with k={} using optimized u64 encoding with Bloom filter", k);
    let min_kmer_count = 2; // Filter k-mers appearing only once
    let kmer_counts_u64 = cpu_backend.count_kmers_u64_filtered(&sequences, k, min_kmer_count);
    info!("Found {} unique k-mers (after filtering singletons)", kmer_counts_u64.len());

    // Build adjacency table for assembly
    info!("Building adjacency table...");
    let mut adjacency = cpu_backend.build_adjacency_u64(&kmer_counts_u64, k);
    info!("Adjacency table built with {} forward edges", adjacency.forward.len());

    // Clean up graph: remove tips and collapse bubbles
    // This reduces noise from sequencing errors and improves assembly quality
    {
        use crate::graph::assembler::cleanup_graph;
        info!("Cleaning up assembly graph (removing tips and bubbles)...");
        let (tips_removed, bubbles_collapsed) = cleanup_graph(&mut adjacency, &kmer_counts_u64, k, min_kmer_count);
        info!("Graph cleanup: removed {} tips, collapsed {} bubbles", tips_removed, bubbles_collapsed);
    }

    // Keep records for polishing if needed
    if polish || (isoforms && (polish_isoforms || compute_tpm)) {
        let reader = open_fastq(input_path);
        records = stream_fastq_records(reader).collect();
    }

    // Perform greedy assembly using u64 k-mers
    info!("Assembling contigs with minimum length: {}", min_len);
    let mut contigs = greedy_assembly_u64(k, &kmer_counts_u64, &adjacency, min_len);

    // Legacy path for GPU backend (uses String k-mers)
    // Convert u64 counts to String counts for compatibility with existing GPU code
    let kmer_counts: HashMap<String, u32> = kmer_counts_u64
        .iter()
        .map(|(&kmer, &count)| (crate::kmer::kmer::decode_kmer(kmer, k), count))
        .collect();

    // Collapse repeats if requested
    if collapse_repeats {
        use crate::graph::simplify::collapse_repeats;
        let before_count = contigs.len();
        contigs = collapse_repeats(contigs, min_repeat_len);
        info!("Collapsed repeats: {} -> {} contigs (removed {})",
              before_count, contigs.len(), before_count - contigs.len());
    }

    // Polish contigs if requested (using parallel implementation for 4-8x speedup)
    if polish {
        use crate::graph::polish::polish_contig_parallel;
        info!("Polishing contigs using aligned reads (window size: {}, parallel)", polish_window);

        let chunk_size = 10000; // Process 10kb chunks in parallel
        for contig in &mut contigs {
            contig.sequence = polish_contig_parallel(&contig.sequence, &records, polish_window, chunk_size);
        }

        info!("Completed contig polishing");
    }

    // Write FASTA output
    let mut writer = FastaWriter::new(output_path);
    for (i, contig) in contigs.iter().enumerate() {
        writer.write_contig(contig, i + 1).expect("Failed to write contig");

        // Write RLE version if requested
        if _use_rle {
            writer.write_rle_contig(contig, i + 1).expect("Failed to write RLE contig");
        }
    }

    info!("Assembly complete: {} contigs written to {}", contigs.len(), output_path);

    // Export metadata if requested
    if export_metadata || json_metadata.is_some() || tsv_metadata.is_some() {
        use crate::io::metadata::generate_metadata;
        info!("Generating contig metadata");
        let meta = generate_metadata(&contigs);

        // Handle standard metadata JSON export
        if export_metadata {
            info!("Exporting contig metadata to JSON");
            let json = serde_json::to_string_pretty(&meta).unwrap();
            let meta_path = format!("{}.contig_meta.json", output_path);
            fs::write(&meta_path, json).expect("Failed to write metadata");
            info!("Metadata written to {}", meta_path);
        }

        // Handle custom JSON metadata path
        if let Some(path) = &json_metadata {
            info!("Writing JSON metadata to custom path: {}", path);
            let json = serde_json::to_string_pretty(&meta).unwrap();
            fs::write(path, json).expect("Failed to write JSON metadata");
        }

        // Handle custom TSV metadata path
        if let Some(path) = &tsv_metadata {
            info!("Writing TSV metadata to custom path: {}", path);
            let mut file = std::fs::File::create(path).expect("Failed to create TSV metadata file");
            writeln!(file, "contig_id\tlength\trle_compression\tgc_content").unwrap();
            for m in &meta {
                writeln!(file, "{}\t{}\t{:.4}\t{:.4}",
                    m.id, m.length, m.rle_compression, m.gc_content).unwrap();
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

        // Use backend for overlap detection (currently unused, kept for future integration)
        info!("Finding overlaps using {}", backend.name());
        let _backend_overlaps = backend.find_overlaps(&contig_seqs, min_overlap, max_mismatches);

        let builder = OverlapGraphBuilder::new(min_overlap, max_mismatches, max_mismatches);
        let graph = builder.build_overlap_graph(&contig_seqs);
        let (_, paths) = builder.stitch_contigs(&graph);

        // Convert backend overlaps to the link format expected by GFA writers
        // Also use the original find_overlaps for compatibility with existing code
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

            // Convert kmer_counts to the expected HashMap<usize, usize> type
            let kmer_counts_converted: HashMap<usize, usize> = kmer_counts
                .iter()
                .filter_map(|(k, v)| {
                    if let Ok(id) = k.parse::<usize>() {
                        Some((id, *v as usize))
                    } else {
                        None
                    }
                })
                .collect();

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
                get_output_filename_fn
            ).unwrap_or_else(|e| {
                warn!("Error processing isoforms: {}", e);
                Vec::new()
            });

            info!("Generated {} transcript isoforms", transcripts.len());

            // Export counts matrix
            if counts_matrix {
                info!("Writing isoform counts matrix");
                let counts_matrix_path = format!("{}_isoform.counts.matrix", output_path);
                if let Err(e) = crate::quant::matrix::write_isoform_counts_matrix(&transcripts, &counts_matrix_path) {
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
                    info!("GTF output complete: {} transcripts written", transcripts.len());
                }
            }

            // Export GFF3 if requested
            if let Some(gff3_path) = &gff3_path {
                info!("Writing isoform GFF3 to: {}", gff3_path);
                use crate::io::gff3::write_gff3;
                if let Err(e) = write_gff3(&transcripts, gff3_path) {
                    warn!("Failed to write GFF3 file: {}", e);
                } else {
                    info!("GFF3 output complete: {} transcripts written", transcripts.len());
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
                use crate::quant::tpm::{count_reads, compute_tpm, write_tpm_table, filter_by_tpm};

                let counts = count_reads(&transcripts, &records);
                let tpms = compute_tpm(&counts, &transcripts);

                // Filter transcripts by TPM if requested
                if min_tpm > 0.0 {
                    info!("Filtering transcripts with TPM < {}", min_tpm);
                    let (filtered_transcripts, filtered_tpms) = filter_by_tpm(&transcripts, &tpms, min_tpm);
                    let filtered_count = transcripts.len() - filtered_transcripts.len();
                    info!("Filtered out {} transcripts with low expression", filtered_count);

                    // Replace transcripts with filtered set
                    transcripts = filtered_transcripts;

                    // Update TPMs to match filtered transcripts
                    let tpm_path = format!("{}_isoform.tpm.tsv", output_path);
                    write_tpm_table(&transcripts, &filtered_tpms, &tpm_path);
                    info!("TPM values written to {}", tpm_path);

                    // Write transcript metrics to JSON if requested
                    if let Some(json_path) = &json_metadata {
                        info!("Writing transcript metrics to JSON: {}", json_path);
                        use crate::io::metadata::write_transcript_metrics;
                        if let Err(e) = write_transcript_metrics(&transcripts, &filtered_tpms, json_path) {
                            warn!("Failed to write transcript metrics to JSON: {}", e);
                        } else {
                            info!("Transcript metrics written to {}", json_path);
                        }
                    }
                } else {
                    let tpm_path = format!("{}_isoform.tpm.tsv", output_path);
                    write_tpm_table(&transcripts, &tpms, &tpm_path);
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
                    use std::collections::HashMap;
                    use crate::io::sam::parse_sam_transcript_hits;
                    use crate::quant::matrix::write_counts_matrix;

                    info!("Processing multi-sample data from {}", sample_file);
                    let sample_content = std::fs::read_to_string(sample_file)
                        .expect("Failed to read sample file");

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

                        info!("Processing sample: {} from alignment {}", sample_name, sam_path);

                        // Parse SAM file for transcript hits
                        let hits = match parse_sam_transcript_hits(sam_path) {
                            Ok(h) => h,
                            Err(e) => {
                                warn!("Failed to parse SAM file {}: {}", sam_path, e);
                                continue;
                            }
                        };

                        // Count hits per transcript
                        let sam_counts: Vec<usize> = transcripts.iter().map(|t| {
                            *hits.get(&format!("transcript_{}", t.id)).unwrap_or(&0)
                        }).collect();

                        // Compute TPM values
                        let sam_tpms = compute_tpm(&sam_counts, &transcripts);
                        sample_tpms.insert(sample_name, sam_tpms);
                    }

                    if !sample_tpms.is_empty() {
                        // Write the matrix output
                        let matrix_path = format!("{}_isoform.counts.matrix", output_path);
                        write_counts_matrix(&sample_tpms, &transcripts, &matrix_path);
                        info!("Multi-sample counts matrix written to {}", matrix_path);
                    } else {
                        warn!("No valid samples found in {}", sample_file);
                    }
                }
            }

            // Apply transcript polishing with long reads if requested
            if let Some(polish_sam_path) = &long_reads {
                info!("Polishing transcript sequences with long read alignments from {}", polish_sam_path);
                use crate::polish::longread::polish_transcripts;

                match polish_transcripts(&mut transcripts, polish_sam_path) {
                    Ok(count) => {
                        info!("Successfully polished {} transcripts", count);
                    },
                    Err(e) => {
                        warn!("Error during transcript polishing: {}", e);
                    }
                }
            }

            // Display transcript evaluation metrics
            let transcript_lengths: Vec<usize> = transcripts.iter()
                .map(|t| t.sequence.len())
                .collect();
            let stats = evaluate_lengths(&transcript_lengths);
            info!(
                "Transcript statistics: {} transcripts, Avg length: {:.1} bp, N50: {} bp",
                stats.total, stats.avg_length, stats.n50
            );
        }
        // Output GFA if requested
        if _output_gfa {
            let gfa_path = get_output_filename(output_path, "gfa");
            info!("Writing GFA to: {}", gfa_path);

            // Write GFA output
            let mut gfa_writer = GfaWriter::new(&gfa_path);

            if _use_rle {
                gfa_writer.write_rle_segments(&contigs).expect("Failed to write RLE GFA segments");
            } else {
                gfa_writer.write_segments(&contigs).expect("Failed to write GFA segments");
            }

            gfa_writer.write_links(&links).expect("Failed to write GFA links");
            gfa_writer.write_assembly_paths(&paths).expect("Failed to write GFA paths");

            info!("GFA output complete: {} segments, {} links, {} paths written",
                  contigs.len(), links.len(), paths.len());
        }

        // Output GFA2 if requested
        if _output_gfa2 {
            let gfa2_path = get_output_filename(output_path, "gfa2");
            info!("Writing GFA2 to: {}", gfa2_path);

            // Write GFA2 output
            let mut gfa2_writer = Gfa2Writer::new(&gfa2_path);
            gfa2_writer.write_segments(&contigs).expect("Failed to write GFA2 segments");
            gfa2_writer.write_links(&links).expect("Failed to write GFA2 links");
            gfa2_writer.write_paths(&paths).expect("Failed to write GFA2 paths");

            info!("GFA2 output complete: {} segments, {} links, {} paths written",
                   contigs.len(), links.len(), paths.len());
        }
    }
}

// Helper function to generate output filenames
fn get_output_filename(output_path: &str, extension: &str) -> String {
    let path = if output_path.ends_with(".gz") {
        output_path.strip_suffix(".gz").unwrap_or(output_path)
    } else {
        output_path
    };
    
    if path.ends_with(".fasta") || path.ends_with(".fa") {
        format!("{}.{}", &path[..path.rfind('.').unwrap_or(path.len())], extension)
    } else {
        format!("{}.{}", path, extension)
    }
} 