use rayon::prelude::*;
use crate::graph::isoform_graph::IsoformGraph;
use crate::graph::isoform_traverse::{find_directed_paths, TranscriptPath};
use crate::graph::transcript::Transcript;
use std::collections::HashMap;
use std::sync::Arc;

/// Process paths in parallel using multiple threads
pub fn parallel_path_discovery(
    graph: &IsoformGraph,
    start_nodes: &[usize],
    end_nodes: &[usize],
    max_depth: usize
) -> Vec<TranscriptPath> {
    // Define chunk size based on number of start nodes
    let chunk_size = (start_nodes.len() / rayon::current_num_threads()).max(1);
    
    // Split start nodes into chunks for parallel processing
    let chunks: Vec<Vec<usize>> = start_nodes.chunks(chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect();
    
    // Process each chunk in parallel
    let results: Vec<Vec<TranscriptPath>> = chunks.par_iter()
        .map(|chunk| {
            find_directed_paths(graph, chunk, end_nodes, max_depth)
        })
        .collect();
    
    // Flatten results
    results.into_iter().flatten().collect()
}

/// Process transcript assembly in parallel
pub fn parallel_transcript_assembly(
    paths: &[TranscriptPath],
    contigs: &HashMap<usize, String>,
    min_confidence: f32
) -> Vec<Transcript> {
    // Filter paths by confidence
    let filtered_paths: Vec<&TranscriptPath> = paths.iter()
        .filter(|path| path.confidence >= min_confidence)
        .collect();
    
    // Convert to Arc for thread-safe sharing
    let contigs_arc = Arc::new(contigs.clone());
    
    // Process in parallel
    filtered_paths.par_iter()
        .enumerate()
        .map(|(idx, path)| {
            let contigs = Arc::clone(&contigs_arc);
            assemble_single_transcript(idx, path, &contigs)
        })
        .collect()
}

/// Assemble a single transcript from a path
fn assemble_single_transcript(
    id: usize,
    path: &TranscriptPath,
    contigs: &HashMap<usize, String>
) -> Transcript {
    let mut sequence = String::new();
    
    // Stitch together contigs along the path
    for (i, &node_id) in path.nodes.iter().enumerate() {
        if let Some(contig_seq) = contigs.get(&node_id) {
            if i == 0 {
                // First contig is added in full
                sequence.push_str(contig_seq);
            } else {
                // Subsequent contigs might overlap with previous
                // For simplicity, we're using a fixed overlap of 20 bp
                // In a real implementation, you'd use the actual overlaps
                let overlap = 20.min(contig_seq.len());
                
                if contig_seq.len() > overlap {
                    sequence.push_str(&contig_seq[overlap..]);
                }
            }
        }
    }
    
    // Create transcript with appropriate metadata
    Transcript::new(
        id,
        sequence,
        path.nodes.clone(),
        path.confidence as f64
    )
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

/// Calculate confidence for a path as the average of edge weights in a DiGraphMap
fn calculate_path_confidence(graph: &IsoformGraph, path: &[usize]) -> f32 {
    if path.len() <= 1 {
        return 1.0; // Single node has perfect confidence
    }
    
    let mut sum_weights = 0.0;
    let mut count = 0;
    
    for i in 0..path.len() - 1 {
        // For DiGraphMap, we use edge_weight directly with node IDs
        if let Some(&weight) = graph.edge_weight(path[i], path[i + 1]) {
            sum_weights += weight;
            count += 1;
        }
    }
    
    if count == 0 {
        0.0
    } else {
        sum_weights / count as f32
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::Graph;
    
    #[test]
    fn test_parallel_assembly() {
        // Create test graph
        let mut graph = Graph::<usize, f32>::new();
        let n1 = graph.add_node(1);
        let n2 = graph.add_node(2);
        let n3 = graph.add_node(3);
        
        graph.add_edge(n1, n2, 0.9);
        graph.add_edge(n2, n3, 0.8);
        
        // Create test contigs
        let mut contigs = HashMap::new();
        contigs.insert(1, "ATCGATCG".to_string());
        contigs.insert(2, "ATCGGCTA".to_string());
        contigs.insert(3, "GCTATAGC".to_string());
        
        // Create test paths
        let path1 = TranscriptPath {
            nodes: vec![1, 2],
            confidence: 0.9,
            length: 16,
        };
        
        let path2 = TranscriptPath {
            nodes: vec![2, 3],
            confidence: 0.8,
            length: 16,
        };
        
        let paths = vec![path1, path2];
        
        // Test parallel assembly
        let transcripts = parallel_transcript_assembly(&paths, &contigs, 0.7);
        
        assert_eq!(transcripts.len(), 2);
        assert_eq!(transcripts[0].path, vec![1, 2]);
        assert_eq!(transcripts[1].path, vec![2, 3]);
    }
}
