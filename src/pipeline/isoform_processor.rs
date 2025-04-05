use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use tracing::info;
use crate::graph::assembler::Contig;
use crate::graph::transcript::Transcript;
use crate::graph::isoform_graph::build_isoform_graph;
use crate::graph::isoform_traverse::{find_directed_paths, filter_paths_by_confidence};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

/// Load contig coverages from an expression data file
pub fn load_contig_coverages(expression_path: &str) -> Result<HashMap<usize, f64>> {
    info!("Loading contig coverage data from {}", expression_path);
    let file = File::open(expression_path)?;
    let reader = BufReader::new(file);
    
    let mut coverage_map = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        
        // Skip comments and empty lines
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        
        // Parse contig ID and coverage
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            if let (Ok(id), Ok(coverage)) = (parts[0].parse::<usize>(), parts[1].parse::<f64>()) {
                coverage_map.insert(id, coverage);
            }
        }
    }
    
    info!("Loaded coverage data for {} contigs", coverage_map.len());
    Ok(coverage_map)
}

/// Process assembled contigs into transcripts
pub fn process_isoforms(
    contigs: &[Contig],
    links: &[(usize, usize, usize)],
    kmer_counts: &HashMap<usize, usize>,
    output_path: &str,
    gtf_path: Option<&str>,
    max_path_depth: usize,
    min_confidence: f64,
    min_tpm: Option<f64>,
    strand_aware: bool,
    bam_path: Option<&str>,
    long_reads: Option<&str>,
    get_output_filename: fn(&str, Option<&str>) -> String,
) -> Result<Vec<Transcript>> {
    info!("Processing isoforms from {} contigs with {} links", contigs.len(), links.len());
    
    // Build expression map from k-mer counts if available
    let mut expression_map = HashMap::new();
    for (contig_id, count) in kmer_counts {
        expression_map.insert(*contig_id, *count as f64);
    }
    
    // Convert contigs to a HashMap for easier lookup
    let mut contig_map = HashMap::new();
    for contig in contigs {
        contig_map.insert(contig.id, contig.sequence.clone());
    }
    
    // Build the isoform graph
    let overlaps_vec = links.to_vec();
    let graph = build_isoform_graph(&contig_map, &overlaps_vec, &expression_map);
    
    // Prepare start and end nodes (for simplicity, using all nodes as candidates)
    let all_nodes: Vec<usize> = graph.nodes().collect();
    let start_nodes = all_nodes.clone();
    let end_nodes = all_nodes;
    
    // Find paths through the graph
    let paths = find_directed_paths(&graph, &start_nodes, &end_nodes, max_path_depth);
    info!("Found {} raw transcript paths", paths.len());
    
    // Filter paths by confidence (using min_confidence as f32, 100 as min path length, and None for high threshold)
    let filtered_paths = filter_paths_by_confidence(&paths, min_confidence as f32, 100, None);
    info!("Filtered to {} high-confidence paths", filtered_paths.len());
    
    // Convert paths to transcripts
    let mut transcripts = Vec::new();
    for (i, path) in filtered_paths.iter().enumerate() {
        // Create transcript with proper ID, sequence, and path info
        let mut sequence = String::new();
        
        // Simple sequence stitching (can replace with actual transcript::stitch_isoform later)
        for node_id in &path.nodes {
            if let Some(node_seq) = contig_map.get(node_id) {
                sequence.push_str(node_seq);
            }
        }
        
        let transcript = Transcript::new(
            i,
            sequence,
            path.nodes.clone(),
            path.confidence as f64
        );
        
        transcripts.push(transcript);
    }
    
    info!("Generated {} transcript isoforms", transcripts.len());
    Ok(transcripts)
}
