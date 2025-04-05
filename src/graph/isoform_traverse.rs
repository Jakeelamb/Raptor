use std::collections::{HashSet, VecDeque};
use petgraph::visit::EdgeRef;
use crate::graph::isoform_graph::IsoformGraph;

/// Represents a path through the isoform graph
#[derive(Debug, Clone)]
pub struct TranscriptPath {
    /// Vector of contig IDs in the path
    pub nodes: Vec<usize>,
    /// Confidence score for this path (0.0-1.0)
    pub confidence: f32,
}

/// Enumerate all possible paths from start nodes up to max_depth
pub fn enumerate_paths(
    graph: &IsoformGraph, 
    start_nodes: &[usize], 
    max_depth: usize
) -> Vec<TranscriptPath> {
    let mut paths = Vec::new();
    
    // Process each start node
    for &start in start_nodes {
        // Start a DFS from this node
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
    }
    
    paths
}

/// Find directed paths from start to end nodes
pub fn find_directed_paths(
    graph: &IsoformGraph, 
    start_nodes: &[usize], 
    end_nodes: &[usize],
    max_depth: usize
) -> Vec<TranscriptPath> {
    let mut paths = Vec::new();
    let end_set: HashSet<usize> = end_nodes.iter().cloned().collect();
    
    // Process each start node
    for &start in start_nodes {
        // Use BFS to enumerate paths
        let mut queue = VecDeque::new();
        queue.push_back((vec![start], HashSet::new()));
        
        while let Some((current_path, mut visited)) = queue.pop_front() {
            if current_path.len() > max_depth {
                continue; // Exceeded max depth
            }
            
            let current_node = *current_path.last().unwrap();
            visited.insert(current_node);
            
            // If we've reached an end node, add the path
            if end_set.contains(&current_node) && current_path.len() > 1 {
                let confidence = calculate_path_confidence(graph, &current_path);
                paths.push(TranscriptPath {
                    nodes: current_path,
                    confidence,
                });
                continue;
            }
            
            // Try all neighbors
            for edge in graph.edges(current_node) {
                let neighbor = edge.target();
                
                // Skip if we've already visited this node (avoid cycles)
                if visited.contains(&neighbor) {
                    continue;
                }
                
                // Create new path with this neighbor
                let mut new_path = current_path.clone();
                new_path.push(neighbor);
                
                // Add to queue to continue BFS
                queue.push_back((new_path, visited.clone()));
            }
        }
    }
    
    paths
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

/// Filter paths by minimum confidence threshold
pub fn filter_paths_by_confidence(
    paths: &[TranscriptPath],
    min_confidence: f32
) -> Vec<TranscriptPath> {
    paths.iter()
        .filter(|path| path.confidence >= min_confidence)
        .cloned()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::graphmap::DiGraphMap;
    
    #[test]
    fn test_enumerate_paths() {
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
        let paths = enumerate_paths(&graph, &[0], 3);
        
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