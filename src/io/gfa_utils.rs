use crate::graph::navigation;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Write};
use tracing::info;

/// Traverses paths in a GFA file and exports them in various formats
pub fn traverse_paths(
    input: &str,
    segments: &str,
    output: &str,
    formats: &str,
    include_edges: bool,
    visualize: bool,
    metadata: bool,
) -> io::Result<usize> {
    info!("Traversing paths in GFA file: {}", input);
    
    // Load GFA file
    let gfa_content = std::fs::read_to_string(input)?;
    
    // Load segment sequences
    let segment_map = navigation::load_segment_sequences(segments)?;
    
    // Parse GFA paths
    let gfa_lines: Vec<String> = gfa_content.lines().map(String::from).collect();
    let paths = navigation::parse_gfa_paths(&gfa_lines);
    
    // Reconstruct path sequences
    let reconstructed_paths = navigation::reconstruct_paths(&paths, &segment_map);
    
    // Process and export based on requested formats
    let format_list: Vec<&str> = formats.split(',').collect();
    
    for format in format_list {
        match format.trim() {
            "fasta" => {
                let fasta_path = format!("{}.fasta", output);
                export_paths_to_fasta(&reconstructed_paths, &fasta_path)?;
                info!("Exported path sequences to FASTA: {}", fasta_path);
            },
            "json" => {
                if metadata {
                    let json_path = format!("{}.json", output);
                    export_paths_to_json(&reconstructed_paths, &paths, &json_path)?;
                    info!("Exported path metadata to JSON: {}", json_path);
                }
            },
            "dot" => {
                if visualize {
                    let dot_path = format!("{}.dot", output);
                    export_paths_to_dot(&paths, &segment_map, include_edges, &dot_path)?;
                    info!("Created DOT graph at {}", dot_path);
                }
            },
            _ => {
                info!("Unsupported export format: {}", format);
            }
        }
    }
    
    info!("Path traversal complete. Found {} paths, exported to requested formats.", paths.len());
    Ok(paths.len())
}

/// Export path sequences to FASTA format
fn export_paths_to_fasta(paths: &HashMap<String, String>, output_path: &str) -> io::Result<()> {
    let mut file = File::create(output_path)?;
    
    for (path_id, sequence) in paths {
        writeln!(file, ">{}", path_id)?;
        
        // Write sequence in chunks of 60 characters
        for i in (0..sequence.len()).step_by(60) {
            let end = (i + 60).min(sequence.len());
            writeln!(file, "{}", &sequence[i..end])?;
        }
    }
    
    Ok(())
}

/// Export path metadata to JSON
fn export_paths_to_json(
    reconstructed_paths: &HashMap<String, String>,
    paths: &HashMap<String, Vec<(String, char)>>,
    output_path: &str
) -> io::Result<()> {
    #[derive(serde::Serialize)]
    struct PathMetadata {
        path_id: String,
        segment_count: usize,
        unique_segments: usize,
        total_length: usize,
        has_inversions: bool,
    }
    
    let mut path_metadatas = Vec::new();
    
    for (id, sequence) in reconstructed_paths {
        // Retrieve the original path segments
        let empty_vec = Vec::new();
        let segments = paths.get(id).unwrap_or(&empty_vec);
        let segment_count = segments.len();
        
        // Count unique segments
        let mut unique_segment_ids = std::collections::HashSet::new();
        for (seg_id, _) in segments {
            unique_segment_ids.insert(seg_id);
        }
        
        // Check if path has inversions
        let has_inversions = segments.iter().any(|(_, dir)| *dir == '-');
        
        // Create metadata
        let meta = PathMetadata {
            path_id: id.clone(),
            segment_count,
            unique_segments: unique_segment_ids.len(),
            total_length: sequence.len(),
            has_inversions,
        };
        
        path_metadatas.push(meta);
    }
    
    // Serialize and write to file
    let json_data = serde_json::to_string_pretty(&path_metadatas)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
    
    std::fs::write(output_path, json_data)?;
    
    Ok(())
}

/// Export paths to DOT format for visualization
fn export_paths_to_dot(
    paths: &HashMap<String, Vec<(String, char)>>,
    segment_map: &HashMap<String, String>,
    include_edges: bool,
    output_path: &str
) -> io::Result<()> {
    let mut file = File::create(output_path)?;
    
    // Start DOT file
    writeln!(file, "digraph G {{")?;
    writeln!(file, "  rankdir=LR;")?;
    writeln!(file, "  node [shape=box style=filled];")?;
    
    // Create node set to avoid duplicates
    let mut nodes = std::collections::HashSet::new();
    let mut edges = std::collections::HashSet::new();
    
    // Process each path
    for (path_id, segments) in paths {
        // Add path as subgraph
        writeln!(file, "  subgraph cluster_{} {{", path_id.replace('-', "_"))?;
        writeln!(file, "    label=\"Path {}\"", path_id)?;
        writeln!(file, "    style=filled;")?;
        writeln!(file, "    color=lightgrey;")?;
        
        // Add nodes for each segment in this path
        for (i, (seg_id, orientation)) in segments.iter().enumerate() {
            let node_id = format!("{}_{}", seg_id, orientation);
            let label = if let Some(seq) = segment_map.get(seg_id) {
                format!("{} ({} bp)", seg_id, seq.len())
            } else {
                seg_id.clone()
            };
            
            let color = match orientation {
                '+' => "lightblue",
                '-' => "salmon",
                _ => "white",
            };
            
            if nodes.insert(node_id.clone()) {
                writeln!(file, "    \"{}\" [label=\"{}\" fillcolor=\"{}\"];", node_id, label, color)?;
            }
            
            // Add edges between consecutive segments
            if i > 0 && include_edges {
                let prev_node_id = format!("{}_{}", segments[i-1].0, segments[i-1].1);
                let edge_id = format!("{}->{}", prev_node_id, node_id);
                
                if edges.insert(edge_id) {
                    writeln!(file, "    \"{}\" -> \"{}\";", prev_node_id, node_id)?;
                }
            }
        }
        
        writeln!(file, "  }}")?;
    }
    
    // End DOT file
    writeln!(file, "}}")?;
    
    Ok(())
} 