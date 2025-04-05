use crate::graph::isoform_graph::IsoformGraph;
use std::collections::{HashSet};
use tracing::debug;

/// Represents a path through the transcript graph
#[derive(Debug, Clone)]
pub struct TranscriptPath {
    /// Sequence of node IDs making up the path
    pub nodes: Vec<usize>,
    
    /// Confidence score for this path
    pub confidence: f32,
    
    /// Total path length in nucleotides
    pub length: usize,
}

/// Find all directed paths from start nodes to end nodes
pub fn find_directed_paths(
    graph: &IsoformGraph,
    start_nodes: &[usize],
    end_nodes: &[usize],
    max_depth: usize,
) -> Vec<TranscriptPath> {
    let mut paths = Vec::new();
    let end_node_set: HashSet<usize> = end_nodes.iter().cloned().collect();
    
    for &start in start_nodes {
        let initial_path = vec![start];
        let visited = HashSet::new();
        
        find_paths_dfs(
            graph, 
            start, 
            &end_node_set, 
            initial_path, 
            visited, 
            0, 
            max_depth, 
            &mut paths
        );
    }
    
    debug!("Found {} raw paths in graph traversal", paths.len());
    
    // Convert paths to TranscriptPath objects with confidence scores
    let mut transcript_paths = Vec::new();
    for path in paths {
        // Calculate confidence as average of edge weights
        let confidence = calculate_path_confidence(graph, &path);
        
        // Estimate path length (simplified - would depend on your data model)
        let length = path.len() * 25; // Assuming average node length of 25 bp
        
        transcript_paths.push(TranscriptPath {
            nodes: path,
            confidence,
            length,
        });
    }
    
    transcript_paths
}

/// Recursive DFS to find all paths
fn find_paths_dfs(
    graph: &IsoformGraph,
    current: usize,
    end_nodes: &HashSet<usize>,
    path: Vec<usize>,
    mut visited: HashSet<usize>,
    depth: usize,
    max_depth: usize,
    results: &mut Vec<Vec<usize>>,
) {
    // Check if we've reached an end node
    if end_nodes.contains(&current) && depth > 0 {
        results.push(path.clone());
        return;
    }
    
    // Check if we've reached max depth
    if depth >= max_depth {
        return;
    }
    
    // Mark current node as visited
    visited.insert(current);
    
    // Explore neighbors
    for neighbor in graph.neighbors(current) {
        if !visited.contains(&neighbor) {
            let mut new_path = path.clone();
            new_path.push(neighbor);
            
            find_paths_dfs(
                graph,
                neighbor,
                end_nodes,
                new_path,
                visited.clone(),
                depth + 1,
                max_depth,
                results,
            );
        }
    }
}

/// Calculate confidence score for a path
fn calculate_path_confidence(graph: &IsoformGraph, path: &[usize]) -> f32 {
    if path.len() <= 1 {
        return 1.0; // Single node path has perfect confidence
    }
    
    let mut sum_weight = 0.0;
    let mut count = 0;
    
    for i in 0..path.len() - 1 {
        if let Some(&weight) = graph.edge_weight(path[i], path[i + 1]) {
            sum_weight += weight;
            count += 1;
        }
    }
    
    if count == 0 {
        return 0.0;
    }
    
    sum_weight / count as f32
}

/// Filter transcript paths by confidence and optional length criteria
pub fn filter_paths_by_confidence(
    paths: &[TranscriptPath],
    min_confidence: f32,
    min_path_len: usize,
    high_confidence_threshold: Option<f32>
) -> Vec<TranscriptPath> {
    paths.iter()
        .filter(|path| {
            // Standard confidence check
            if path.confidence >= min_confidence {
                // Standard length check
                if path.length >= min_path_len {
                    return true;
                }
                
                // Short paths can pass if they have high confidence
                if let Some(high_threshold) = high_confidence_threshold {
                    if path.confidence >= high_threshold {
                        return true;
                    }
                }
            }
            
            false
        })
        .cloned()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::graphmap::DiGraphMap;
    
    #[test]
    fn test_find_directed_paths() {
        // Create a test graph with multiple paths
        let mut graph = DiGraphMap::new();
        
        // Add nodes and edges
        for i in 0..5 {
            graph.add_node(i);
        }
        
        // Path 1: 0 -> 1 -> 3 -> 4
        graph.add_edge(0, 1, 0.9);
        graph.add_edge(1, 3, 0.8);
        graph.add_edge(3, 4, 0.7);
        
        // Path 2: 0 -> 2 -> 4
        graph.add_edge(0, 2, 0.6);
        graph.add_edge(2, 4, 0.5);
        
        // Find paths from 0 to 4
        let paths = find_directed_paths(&graph, &[0], &[4], 5);
        
        // Should find two paths
        assert_eq!(paths.len(), 2);
        
        // Check path contents
        let path_strings: Vec<String> = paths.iter()
            .map(|p| format!("{:?}", p.nodes))
            .collect();
            
        assert!(path_strings.contains(&"[0, 1, 3, 4]".to_string()));
        assert!(path_strings.contains(&"[0, 2, 4]".to_string()));
        
        // Check confidence scores
        let path1 = paths.iter().find(|p| p.nodes == vec![0, 1, 3, 4]).unwrap();
        let path2 = paths.iter().find(|p| p.nodes == vec![0, 2, 4]).unwrap();
        
        assert!((path1.confidence - 0.8).abs() < 0.01); // Average of 0.9, 0.8, 0.7
        assert!((path2.confidence - 0.55).abs() < 0.01); // Average of 0.6 and 0.5
    }
} 