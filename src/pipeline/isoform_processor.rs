use crate::graph::assembler::Contig;
use crate::graph::isoform_graph::{build_isoform_graph, find_start_nodes, IsoformGraph};
use crate::graph::isoform_parallel::{parallel_path_discovery, parallel_transcript_assembly, parallel_path_filtering};
use crate::graph::transcript::{calculate_transcript_stats, Transcript};
use crate::io::transcript_io::{TranscriptFastaWriter, TranscriptGtfWriter, write_transcript_json, add_transcripts_to_gfa, write_transcript_stats};
use std::collections::HashMap;
use tracing::info;

/// Processes transcripts/isoforms from assembled contigs
pub fn process_isoforms(
    contigs: &[Contig],
    links: &[(usize, usize, usize)],
    kmer_counts: &HashMap<String, u32>,
    output_path: &str,
    gtf_path: Option<&str>,
    max_path_depth: usize,
    min_confidence: f64,
    get_output_filename: fn(&str, &str) -> String,
) -> Vec<Transcript> {
    info!("Performing isoform inference and transcript path export with parallel traversal");
    
    // Generate coverage data from k-mer counts
    let mut coverage_map = HashMap::new();
    for (i, contig) in contigs.iter().enumerate() {
        // Estimate coverage based on k-mer abundances of this contig's k-mers
        let mut coverage_sum = 0.0;
        let mut count = 0;
        
        for kmer in &contig.kmer_path {
            if let Some(&abundance) = kmer_counts.get(kmer) {
                coverage_sum += abundance as f64;
                count += 1;
            }
        }
        
        let coverage = if count > 0 {
            coverage_sum / count as f64
        } else {
            1.0 // Default coverage if no k-mers
        };
        
        coverage_map.insert(i, coverage);
    }
    
    // Build isoform graph using coverage information
    info!("Building isoform graph from {} contigs and {} links", contigs.len(), links.len());
    let isoform_graph = build_isoform_graph(contigs, links, &coverage_map);
    let start_nodes = find_start_nodes(&isoform_graph);
    
    info!("Found {} potential transcript start points", start_nodes.len());
    
    // Define traversal parameters - use user-provided values
    info!("Using max path depth: {}, min confidence: {:.2}", max_path_depth, min_confidence);
    
    // Enumerate transcript paths using parallel traversal
    info!("Enumerating transcript paths in parallel (max depth: {})", max_path_depth);
    let transcript_paths = parallel_path_discovery(&isoform_graph, &start_nodes, max_path_depth);
    info!("Found {} potential transcript paths", transcript_paths.len());
    
    // Filter paths by confidence in parallel
    let filtered_paths = parallel_path_filtering(&transcript_paths, min_confidence as f32);
    info!("Filtered to {} paths with confidence >= {}", filtered_paths.len(), min_confidence);
    
    // Process paths into transcripts in parallel
    info!("Assembling transcripts from paths in parallel");
    let batch_size = 50; // Process in batches for better performance
    let transcripts = parallel_transcript_assembly(&filtered_paths, contigs, links, batch_size);
    
    // Filter similar transcripts to reduce redundancy
    let similarity_threshold = 0.85;
    info!("Filtering similar transcripts (threshold: {})", similarity_threshold);
    let filtered_transcripts = crate::graph::isoform_filter::filter_similar_transcripts(
        &transcripts, similarity_threshold);
    
    info!("Filtered {} transcripts to {} unique isoforms", 
          transcripts.len(), filtered_transcripts.len());
    
    // Generate stats
    let stats = calculate_transcript_stats(&filtered_transcripts);
    
    // Output files - use default paths
    let fasta_path = get_output_filename(output_path, "isoforms.fasta");
    let default_gtf_path = get_output_filename(output_path, "isoforms.gtf");
    let gfa_path = get_output_filename(output_path, "isoforms.gfa");
    let json_path = get_output_filename(output_path, "isoforms.json");
    let stats_path = get_output_filename(output_path, "isoforms.stats.json");
    
    // Use custom GTF path if provided
    let final_gtf_path = gtf_path.unwrap_or(&default_gtf_path);
    
    // Write outputs
    info!("Writing transcript outputs");
    let mut fasta_writer = TranscriptFastaWriter::new(&fasta_path);
    fasta_writer.write_transcripts(&filtered_transcripts)
        .expect("Failed to write isoform FASTA");
    
    let mut gtf_writer = TranscriptGtfWriter::new(final_gtf_path);
    gtf_writer.write_transcripts(&filtered_transcripts)
        .expect("Failed to write isoform GTF");
    
    add_transcripts_to_gfa(&gfa_path, &filtered_transcripts)
        .expect("Failed to write isoform GFA");
        
    write_transcript_json(&json_path, &filtered_transcripts)
        .expect("Failed to write isoform JSON");
        
    write_transcript_stats(&stats_path, &stats)
        .expect("Failed to write isoform stats");
    
    info!("Isoform inference complete. Written to:");
    info!("  FASTA: {}", fasta_path);
    info!("  GTF: {}", final_gtf_path);
    info!("  GFA: {}", gfa_path);
    info!("  JSON: {}", json_path);
    info!("  Stats: {}", stats_path);
    
    filtered_transcripts
}
