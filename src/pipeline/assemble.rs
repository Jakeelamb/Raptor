use crate::io::fastq::{open_fastq, read_fastq_records, stream_fastq_records};
use crate::io::fasta::FastaWriter;
use crate::io::gfa::GfaWriter;
use crate::io::gfa2::Gfa2Writer;
use crate::kmer::kmer::canonical_kmer;
use crate::kmer::variable_k::{optimal_k, kmer_coverage_histogram, select_best_k};
use crate::graph::assembler::greedy_assembly;
use crate::graph::overlap::find_overlaps;
use crate::graph::stitch::{OverlapGraphBuilder, Path};
use crate::dist::scheduler::{DistributedConfig, run_distributed_assembly};
use std::collections::HashMap;
use std::path::Path as FilePath;
use tracing::{info, warn};
use std::fs;
use std::io::Write;
pub fn assemble_reads(
    input_path: &str, 
    output_path: &str, 
    min_len: usize, 
    output_gfa: bool, 
    output_gfa2: bool,
    adaptive_k: bool, 
    use_rle: bool,
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
    samples: Option<String>,
    min_tpm: f64,
    polish_reads: Option<String>,
) { 
    assemble_reads_old(
        input_path, output_path, min_len, output_gfa, output_gfa2, 
        adaptive_k, use_rle, collapse_repeats, min_repeat_len, 
        polish, polish_window, streaming, export_metadata, 
        json_metadata, tsv_metadata, isoforms, gtf_path, gff3_path, 
        max_path_depth, min_confidence, compute_tpm, polish_isoforms,
        samples, min_tpm, polish_reads
    ); 
}

pub fn assemble_reads_old(
    input_path: &str, 
    output_path: &str, 
    min_len: usize, 
    output_gfa: bool, 
    output_gfa2: bool,
    adaptive_k: bool, 
    use_rle: bool,
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
    samples: Option<String>,
    min_tpm: f64,
    polish_reads: Option<String>,
) {
    info!("Starting assembly from: {}", input_path);
    
    // Determine k-mer size - either adaptive or fixed optimal
    let max_k = 41;
    let mut k = 31; // Default k-mer size
    
    let mut kmer_counts = HashMap::new();
    let mut records = Vec::new();
    
    if streaming {
        info!("Using streaming mode for assembly");
        
        // First pass: count k-mers to determine optimal k
        info!("First pass: counting k-mers for optimal k selection");
        let reader = open_fastq(input_path);
        let mut temp_sequences = Vec::new();
        
        // Sample some sequences for k-mer optimization
        for (i, record) in stream_fastq_records(reader).enumerate() {
            if i < 10000 { // Sample only the first 10,000 sequences
                temp_sequences.push(record.sequence.clone());
            }
            if i >= 10000 && !adaptive_k {
                break; // Only need samples for k determination if not using full dataset
            }
        }
        
        // Determine k-mer size
        if adaptive_k {
            let hist = kmer_coverage_histogram(&temp_sequences, max_k);
            k = select_best_k(&hist);
            info!("Using adaptive k-mer: {}", k);
        } else {
            k = optimal_k(&temp_sequences, max_k);
            info!("Using optimal k-mer size: {}", k);
        }
        
        // Second pass: count k-mers for assembly
        info!("Second pass: counting k-mers with k={}", k);
        let reader = open_fastq(input_path);
        
        // In streaming mode, accumulate minimal required data
        for record in stream_fastq_records(reader) {
            if record.sequence.len() < k {
                continue;
            }
            
            // Count k-mers for assembly
            for i in 0..=record.sequence.len() - k {
                if let Some(kmer) = canonical_kmer(&record.sequence[i..i + k]) {
                    *kmer_counts.entry(kmer).or_insert(0) += 1;
                }
            }
            
            // Only store records if polishing is needed or for isoform processing
            if polish || (isoforms && (polish_isoforms || compute_tpm)) {
                records.push(record);
            }
        }
    } else {
        info!("Using in-memory mode for assembly");
        
        // Read all records into memory
        let reader = open_fastq(input_path);
        records = read_fastq_records(reader).collect();
        
        // Extract sequences for k-mer optimization
        let sequences: Vec<String> = records.iter()
            .map(|record| record.sequence.clone())
            .collect();
        
        // Determine k-mer size
        if adaptive_k {
            let hist = kmer_coverage_histogram(&sequences, max_k);
            k = select_best_k(&hist);
            info!("Using adaptive k-mer: {}", k);
        } else {
            k = optimal_k(&sequences, max_k);
            info!("Using optimal k-mer size: {}", k);
        }
        
        // Count k-mers
        for record in &records {
            if record.sequence.len() < k {
                continue;
            }
            
            for i in 0..=record.sequence.len() - k {
                if let Some(kmer) = canonical_kmer(&record.sequence[i..i + k]) {
                    *kmer_counts.entry(kmer).or_insert(0) += 1;
                }
            }
        }
    }
    
    // Perform greedy assembly
    info!("Assembling contigs with minimum length: {}", min_len);
    let mut contigs = greedy_assembly(k, &kmer_counts, min_len);
    
    // Collapse repeats if requested
    if collapse_repeats {
        use crate::graph::simplify::collapse_repeats;
        let before_count = contigs.len();
        contigs = collapse_repeats(contigs, min_repeat_len);
        info!("Collapsed repeats: {} -> {} contigs (removed {})", 
              before_count, contigs.len(), before_count - contigs.len());
    }
    
    // Polish contigs if requested
    if polish {
        use crate::graph::polish::polish_contig;
        info!("Polishing contigs using aligned reads (window size: {})", polish_window);
        
        for contig in &mut contigs {
            contig.sequence = polish_contig(&contig.sequence, &records, polish_window);
        }
        
        info!("Completed contig polishing");
    }
    
    // Write FASTA output
    let mut writer = FastaWriter::new(output_path);
    for (i, contig) in contigs.iter().enumerate() {
        writer.write_contig(contig, i + 1).expect("Failed to write contig");
        
        // Write RLE version if requested
        if use_rle {
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
    if (output_gfa || output_gfa2 || isoforms) && !contigs.is_empty() {
        // Extract raw sequences for overlap finding
        let contig_seqs: Vec<String> = contigs.iter().map(|c| c.sequence.clone()).collect();
        
        // Build overlap graph
        let min_overlap = (k / 2).max(15); // Use at least half of k but minimum 15bp
        let max_mismatches = 3; // Allow up to 3 mismatches in the overlap
        
        let builder = OverlapGraphBuilder::new(min_overlap, max_mismatches, max_mismatches);
        let graph = builder.build_overlap_graph(&contig_seqs);
        let (_, paths) = builder.stitch_contigs(&graph);
        
        // Find links (simple edges) for GFA format
        let links = find_overlaps(&contigs, min_overlap, max_mismatches);
        

        // Process isoforms if requested
        if isoforms {
            use crate::pipeline::isoform_processor::process_isoforms;
            let mut transcripts = process_isoforms(&contigs, &links, &kmer_counts, output_path, gtf_path.as_deref(), max_path_depth, min_confidence, get_output_filename);
            
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
                if let Some(sample_file) = &samples {
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
            if let Some(polish_sam_path) = &polish_reads {
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
            use crate::eval::metrics::evaluate;
            let stats = evaluate(&transcripts);
            info!(
                "Transcript statistics: {} transcripts, Avg length: {:.1} bp, N50: {} bp",
                stats.total, stats.avg_length, stats.n50
            );
        }
        // Output GFA if requested
        if output_gfa {
            let gfa_path = get_output_filename(output_path, "gfa");
            info!("Writing GFA to: {}", gfa_path);
            
            // Write GFA output
            let mut gfa_writer = GfaWriter::new(&gfa_path);
            
            if use_rle {
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
        if output_gfa2 {
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

/// Perform distributed assembly by partitioning input sequences
pub fn distributed_assembly(
    input_path: &str, 
    output_path: &str, 
    min_len: usize, 
    output_gfa: bool, 
    adaptive_k: bool, 
    use_rle: bool,
    buckets: usize,
    threads: usize
) {
    info!("Starting distributed assembly from: {}", input_path);
    
    // Read input file
    let reader = open_fastq(input_path);
    let records: Vec<_> = read_fastq_records(reader).collect();
    
    // Extract sequences
    let sequences: Vec<String> = records.iter()
        .map(|record| record.sequence.clone())
        .collect();
    
    info!("Partitioning {} sequences into {} buckets", sequences.len(), buckets);
    
    // Create output directory structure
    let output_dir = FilePath::new(output_path).parent().unwrap_or(FilePath::new(".")).to_str().unwrap();
    let output_name = FilePath::new(output_path).file_stem().unwrap_or_default().to_str().unwrap();
    let distributed_dir = format!("{}/{}_distributed", output_dir, output_name);
    
    // Configure distributed assembly
    let config = DistributedConfig {
        buckets,
        threads_per_job: threads / buckets.max(1),
        use_slurm: false, // Use local execution by default
        ..Default::default()
    };
    
    // Run distributed assembly
    let output_files = run_distributed_assembly(&sequences, &distributed_dir, min_len, &config);
    
    info!("Distributed assembly complete. Partitioned results written to: {}", distributed_dir);
    info!("Output files: {:?}", output_files);
    
    // Optionally merge results
    // TODO: Implement merging of distributed results
} 