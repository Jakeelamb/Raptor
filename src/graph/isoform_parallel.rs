use rayon::prelude::*;
use crate::graph::isoform_graph::IsoformGraph;
use crate::graph::isoform_traverse::TranscriptPath;
use crate::graph::transcript::{Transcript, stitch_isoform};
use crate::graph::assembler::Contig;
use std::collections::HashSet;
use petgraph::visit::EdgeRef;

/// Process start nodes in parallel to find all possible transcript paths
pub fn parallel_path_discovery(
    graph: &IsoformGraph,
    start_nodes: &[usize],
    max_depth: usize
) -> Vec<TranscriptPath> {
    start_nodes.par_iter()
        .flat_map(|&start| {
            // For each start node, find all possible paths
            let mut paths = Vec::new();
            let mut stack = Vec::new();
            let initial_path = vec![start];
            stack.push((initial_path, HashSet::new()));
            
            while let Some((current_path, mut visited)) = stack.pop() {
                let current_node = *current_path.last().unwrap();
                visited.insert(current_node);
                
                // We've reached max depth, add the path
                if current_path.len() >= max_depth {
                    // Calculate confidence as average of edge weights
                    let confidence = calculate_path_confidence(graph, &current_path);
                    paths.push(TranscriptPath {
                        nodes: current_path,
                        confidence,
                    });
                    continue;
                }
                
                // Try all neighbors
                let mut has_neighbors = false;
                for edge in graph.edges(current_node) {
                    let neighbor = edge.target();
                    
                    // Skip if we've already visited this node (avoid cycles)
                    if visited.contains(&neighbor) {
                        continue;
                    }
                    
                    has_neighbors = true;
                    
                    // Create new path with this neighbor
                    let mut new_path = current_path.clone();
                    new_path.push(neighbor);
                    
                    // Push to stack to continue DFS
                    stack.push((new_path, visited.clone()));
                }
                
                // If this is a leaf node (no outgoing edges), add the path
                if !has_neighbors && !current_path.is_empty() {
                    let confidence = calculate_path_confidence(graph, &current_path);
                    paths.push(TranscriptPath {
                        nodes: current_path,
                        confidence,
                    });
                }
            }
            
            paths
        })
        .collect()
}

/// Process transcript assembly in parallel
pub fn parallel_transcript_assembly(
    paths: &[TranscriptPath],
    contigs: &[Contig],
    links: &[(usize, usize, usize)],
    batch_size: usize
) -> Vec<Transcript> {
    // Process paths in parallel using batches for better performance
    let transcripts: Vec<_> = paths.par_chunks(batch_size.max(1))
        .flat_map(|chunk| {
            chunk.iter().enumerate().map(|(i, path)| {
                let sequence = stitch_isoform(contigs, &path.nodes, links);
                Transcript {
                    id: i, // Temporary ID, will be fixed later
                    sequence: sequence.clone(),
                    path: path.nodes.clone(),
                    confidence: path.confidence as f64,
                    length: sequence.len(),
                }
            }).collect::<Vec<_>>()
        })
        .collect();
    
    // Fix transcript IDs
    transcripts.into_iter()
        .enumerate()
        .map(|(i, mut t)| {
            t.id = i;
            t
        })
        .collect()
}

/// Filter transcript paths by confidence score in parallel
pub fn parallel_path_filtering(
    paths: &[TranscriptPath],
    min_confidence: f32
) -> Vec<TranscriptPath> {
    paths.par_iter()
        .filter(|path| path.confidence >= min_confidence)
        .cloned()
        .collect()
}

/// Calculate confidence for a path as the average of edge weights
fn calculate_path_confidence(graph: &IsoformGraph, path: &[usize]) -> f32 {
    if path.len() <= 1 {
        return 1.0; // Single node has perfect confidence
    }
    
    let mut sum_weights = 0.0;
    let mut count = 0;
    
    for i in 0..path.len() - 1 {
        if let Some(weight) = graph.edge_weight(path[i], path[i + 1]) {
            sum_weights += *weight;
            count += 1;
        }
    }
    
    if count == 0 {
        0.0
    } else {
        sum_weights / count as f32
    }
}
