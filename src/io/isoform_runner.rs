use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::time::Instant;
use log::{info, warn};

use crate::graph::isoform_graph::{build_isoform_graph, find_start_nodes, find_end_nodes};
use crate::graph::isoform_traverse::{find_directed_paths, filter_paths_by_confidence};
use crate::graph::transcript::{assemble_transcripts, calculate_transcript_stats};
use crate::graph::isoform_filter::{filter_similar_transcripts, merge_transcripts};
use crate::io::transcript_io::{TranscriptFastaWriter, TranscriptGtfWriter, add_transcripts_to_gfa, write_transcript_stats};
use crate::io::gfa::{read_gfa_contigs, read_gfa_links};

/// Load expression data from a TSV file
/// Format: contig_id, coverage
pub fn load_expression_data(expression_path: &str) -> HashMap<usize, f64> {
    let file = File::open(expression_path).unwrap_or_else(|_| panic!("Failed to open expression file: {}", expression_path));
    let reader = BufReader::new(file);
    
    let mut expression_data = HashMap::new();
    
    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            warn!("Invalid line in expression data: {}", line);
            continue;
        }
        
        let contig_id = parts[0].parse::<usize>().unwrap_or_else(|_| {
            panic!("Invalid contig id: {}", parts[0]);
        });
        
        let coverage = parts[1].parse::<f64>().unwrap_or_else(|_| {
            panic!("Invalid coverage value: {}", parts[1]);
        });
        
        expression_data.insert(contig_id, coverage);
    }
    
    info!("Loaded expression data for {} contigs", expression_data.len());
    expression_data
}

/// Run the isoform reconstruction pipeline
pub fn run_isoform_reconstruction(
    input_gfa: &str,
    expression_path: &str,
    output_prefix: &str,
    min_confidence: f64,
    max_depth: usize,
    formats: &str,
    output_stats: bool,
    similarity_threshold: Option<f64>,
    merge_similar: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let start_time = Instant::now();
    info!("Starting isoform reconstruction");
    
    // Load contigs and overlaps from GFA
    info!("Loading contigs from GFA");
    let contigs = read_gfa_contigs(input_gfa)?;
    info!("Loaded {} contigs", contigs.len());
    
    info!("Loading overlaps from GFA");
    let overlaps = read_gfa_links(input_gfa)?;
    info!("Loaded {} overlaps", overlaps.len());
    
    // Load expression data
    info!("Loading expression data");
    let expression_data = load_expression_data(expression_path);
    
    // Build isoform graph
    info!("Building isoform graph");
    let graph = build_isoform_graph(&contigs, &overlaps, &expression_data);
    let node_count = graph.node_count();
    let edge_count = graph.edge_count();
    info!("Built graph with {} nodes and {} edges", node_count, edge_count);
    
    // Find start and end nodes
    let start_nodes = find_start_nodes(&graph);
    let end_nodes = find_end_nodes(&graph);
    info!("Identified {} potential start nodes and {} potential end nodes", 
           start_nodes.len(), end_nodes.len());
    
    // Find transcript paths
    info!("Finding directed paths through the graph (max depth: {})", max_depth);
    let paths = find_directed_paths(&graph, &start_nodes, &end_nodes, max_depth);
    info!("Found {} potential transcript paths", paths.len());
    
    // Filter paths by confidence
    info!("Filtering paths by confidence (min: {})", min_confidence);
    let min_confidence_f32 = min_confidence as f32;
    let min_path_len = 50; // Minimum path length for standard filtering
    let high_conf_threshold = Some(0.9); // Higher threshold for shorter paths
    info!("Using minimum path length {} with high confidence threshold {:.2}", 
          min_path_len, high_conf_threshold.unwrap_or(1.0));
    let filtered_paths = filter_paths_by_confidence(&paths, min_confidence_f32, min_path_len, high_conf_threshold);
    info!("Retained {} transcript paths after filtering", filtered_paths.len());
    
    // Assemble transcripts
    info!("Assembling transcripts");
    let mut transcripts = assemble_transcripts(&filtered_paths, &contigs, &overlaps, None);
    info!("Assembled {} transcripts", transcripts.len());
    
    // Apply similarity filtering if requested
    if let Some(threshold) = similarity_threshold {
        info!("Filtering similar transcripts (similarity threshold: {})", threshold);
        if merge_similar {
            info!("Merging similar transcripts");
            transcripts = merge_transcripts(&transcripts, threshold);
            info!("After merging, {} transcripts remain", transcripts.len());
        } else {
            transcripts = filter_similar_transcripts(&transcripts, threshold);
            info!("After filtering, {} transcripts remain", transcripts.len());
        }
    }
    
    // Output transcripts in requested formats
    let format_list: Vec<&str> = formats.split(',').collect();
    
    if format_list.contains(&"fasta") {
        let fasta_path = format!("{}.transcripts.fa", output_prefix);
        info!("Writing transcripts to FASTA: {}", fasta_path);
        let mut fasta_writer = TranscriptFastaWriter::new(&fasta_path);
        fasta_writer.write_transcripts(&transcripts)?;
    }
    
    if format_list.contains(&"gtf") {
        let gtf_path = format!("{}.transcripts.gtf", output_prefix);
        info!("Writing transcripts to GTF: {}", gtf_path);
        let mut gtf_writer = TranscriptGtfWriter::new(&gtf_path);
        gtf_writer.write_transcripts(&transcripts)?;
    }
    
    if format_list.contains(&"gfa") {
        let gfa_out_path = format!("{}.transcripts.gfa", output_prefix);
        
        // If input GFA exists, copy it to output GFA first
        if Path::new(input_gfa).exists() && input_gfa != gfa_out_path {
            info!("Copying input GFA to output GFA");
            std::fs::copy(input_gfa, &gfa_out_path)?;
        }
        
        info!("Adding transcript paths to GFA: {}", gfa_out_path);
        add_transcripts_to_gfa(&gfa_out_path, &transcripts)?;
    }
    
    // Calculate and output statistics
    if output_stats {
        let stats_path = format!("{}.transcript_stats.json", output_prefix);
        info!("Calculating transcript statistics");
        let stats = calculate_transcript_stats(&transcripts);
        
        info!("Writing transcript statistics to: {}", stats_path);
        write_transcript_stats(&stats_path, &stats)?;
        
        // Print some key statistics
        println!("Transcript count: {}", stats.get("count").unwrap_or(&0.0));
        println!("Average length: {:.1}", stats.get("mean_length").unwrap_or(&0.0));
        println!("N50: {}", stats.get("n50").unwrap_or(&0.0));
        println!("Average confidence: {:.4}", stats.get("mean_confidence").unwrap_or(&0.0));
    }
    
    let elapsed = start_time.elapsed();
    info!("Isoform reconstruction completed in {:.2}s", elapsed.as_secs_f64());
    
    Ok(())
} 