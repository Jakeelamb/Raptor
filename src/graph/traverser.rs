use std::collections::HashSet;
use rayon::prelude::*;
use petgraph::visit::EdgeRef;
use crate::graph::isoform_graph::IsoformGraph;
use crate::graph::isoform_traverse::TranscriptPath;

/// Parallel version of enumerate_paths that scales better for complex graphs
pub fn enumerate_paths_parallel(
    graph: &IsoformGraph,
    start_nodes: &[usize],
    max_depth: usize
) -> Vec<TranscriptPath> {
    // Process start nodes in parallel
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

/// Batch process paths in parallel for high performance
pub fn process_paths_parallel<F, T>(
    paths: &[TranscriptPath],
    batch_size: usize,
    processor: F
) -> Vec<T>
where
    F: Fn(&TranscriptPath) -> T + Sync,
    T: Send,
{
    paths.par_chunks(batch_size.max(1))
        .flat_map(|chunk| {
            chunk.iter()
                .map(|path| processor(path))
                .collect::<Vec<_>>()
        })
        .collect()
}

/// Filter a large number of paths in parallel by confidence threshold
pub fn filter_paths_parallel(
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

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::graphmap::DiGraphMap;
    
    #[test]
    fn test_enumerate_paths_parallel() {
        // Create a simple test graph
        let mut graph = DiGraphMap::new();
        
        // Add nodes and edges
        graph.add_node(0);
        graph.add_node(1);
        graph.add_node(2);
        graph.add_node(3);
        
        graph.add_edge(0, 1, 0.9);
        graph.add_edge(1, 2, 0.8);
        graph.add_edge(1, 3, 0.6);
        
        // Enumerate paths
        let paths = enumerate_paths_parallel(&graph, &[0], 3);
        
        // Should find two paths: 0->1->2 and 0->1->3
        assert_eq!(paths.len(), 2);
        
        // Check if paths are as expected
        let path_strings: Vec<String> = paths.iter()
            .map(|p| format!("{:?}", p.nodes))
            .collect();
            
        assert!(path_strings.contains(&"[0, 1, 2]".to_string()));
        assert!(path_strings.contains(&"[0, 1, 3]".to_string()));
        
        // Check confidence scores
        let path1 = paths.iter().find(|p| p.nodes == vec![0, 1, 2]).unwrap();
        let path2 = paths.iter().find(|p| p.nodes == vec![0, 1, 3]).unwrap();
        
        assert!((path1.confidence - 0.85).abs() < 0.01); // Average of 0.9 and 0.8
        assert!((path2.confidence - 0.75).abs() < 0.01); // Average of 0.9 and 0.6
    }
    
    #[test]
    fn test_filter_paths_parallel() {
        // Create test paths
        let paths = vec![
            TranscriptPath { nodes: vec![0, 1, 2], confidence: 0.9 },
            TranscriptPath { nodes: vec![0, 1, 3], confidence: 0.7 },
            TranscriptPath { nodes: vec![0, 2, 3], confidence: 0.3 },
        ];
        
        // Filter with threshold 0.6
        let filtered = filter_paths_parallel(&paths, 0.6);
        
        // Should keep first two paths only
        assert_eq!(filtered.len(), 2);
        assert!(filtered.iter().any(|p| p.nodes == vec![0, 1, 2]));
        assert!(filtered.iter().any(|p| p.nodes == vec![0, 1, 3]));
    }
}
