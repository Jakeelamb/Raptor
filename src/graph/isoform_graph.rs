use petgraph::graphmap::DiGraphMap;
use std::collections::HashMap;
use crate::graph::assembler::Contig;

/// Type definition for IsoformGraph - a directed graph where:
/// - Nodes are contig IDs
/// - Edge weights are confidence scores based on coverage
pub type IsoformGraph = DiGraphMap<usize, f32>;

/// Build a directed isoform graph from contigs, overlaps, and expression data
pub fn build_isoform_graph(
    contigs: &[Contig],
    overlaps: &[(usize, usize, usize)],
    expression_data: &HashMap<usize, f64>,
) -> IsoformGraph {
    let mut graph = DiGraphMap::new();
    
    // Add all contigs as nodes
    for contig in contigs {
        graph.add_node(contig.id);
    }
    
    // Add edges based on overlaps
    for &(from, to, _overlap) in overlaps {
        // For each overlap, calculate edge confidence based on expression
        // The confidence will be higher if both contigs have similar expression
        let from_expr = expression_data.get(&from).cloned().unwrap_or(1.0);
        let to_expr = expression_data.get(&to).cloned().unwrap_or(1.0);
        
        // Calculate confidence as normalized similarity in expression
        let max_expr = from_expr.max(to_expr);
        let min_expr = from_expr.min(to_expr);
        
        // Prevent division by zero
        let confidence = if max_expr > 0.0 {
            (min_expr / max_expr) as f32
        } else {
            0.0
        };
        
        // Only add edges with minimum confidence
        if confidence > 0.1 {
            graph.add_edge(from, to, confidence);
        }
    }
    
    // Remove low-weight edges that might be spurious
    prune_low_confidence_edges(&mut graph);
    
    graph
}

/// Remove low confidence edges, especially when better alternatives exist
fn prune_low_confidence_edges(graph: &mut IsoformGraph) {
    let nodes: Vec<usize> = graph.nodes().collect();
    
    for &node in &nodes {
        let mut outgoing_edges = Vec::new();
        
        // Collect all outgoing edges with their weights
        for neighbor in graph.neighbors(node) {
            if let Some(weight) = graph.edge_weight(node, neighbor) {
                outgoing_edges.push((neighbor, *weight));
            }
        }
        
        // If there are multiple outgoing edges, keep only the high confidence ones
        if outgoing_edges.len() > 1 {
            // Sort by descending confidence
            outgoing_edges.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
            
            // Get the highest confidence
            let max_confidence = outgoing_edges[0].1;
            
            // Keep only edges with confidence at least 60% of the best one
            let threshold = max_confidence * 0.6;
            
            for (neighbor, weight) in outgoing_edges {
                if weight < threshold {
                    graph.remove_edge(node, neighbor);
                }
            }
        }
    }
}

/// Find nodes with no incoming edges (potential transcript start points)
pub fn find_start_nodes(graph: &IsoformGraph) -> Vec<usize> {
    let mut start_nodes = Vec::new();
    
    for node in graph.nodes() {
        let has_incoming = graph.neighbors_directed(node, petgraph::Direction::Incoming).count() > 0;
        if !has_incoming {
            start_nodes.push(node);
        }
    }
    
    start_nodes
}

/// Find nodes with no outgoing edges (potential transcript end points)
pub fn find_end_nodes(graph: &IsoformGraph) -> Vec<usize> {
    let mut end_nodes = Vec::new();
    
    for node in graph.nodes() {
        let has_outgoing = graph.neighbors_directed(node, petgraph::Direction::Outgoing).count() > 0;
        if !has_outgoing {
            end_nodes.push(node);
        }
    }
    
    end_nodes
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_build_isoform_graph() {
        // Create sample contigs
        let contigs = vec![
            Contig { id: 0, sequence: "ATCG".to_string(), kmer_path: vec![] },
            Contig { id: 1, sequence: "TCGA".to_string(), kmer_path: vec![] },
            Contig { id: 2, sequence: "CGAT".to_string(), kmer_path: vec![] },
        ];
        
        // Create sample overlaps
        let overlaps = vec![
            (0, 1, 3), // ATCG -> TCGA (overlap of 3)
            (1, 2, 3), // TCGA -> CGAT (overlap of 3)
        ];
        
        // Create sample expression data
        let mut expression_data = HashMap::new();
        expression_data.insert(0, 10.0);
        expression_data.insert(1, 12.0);
        expression_data.insert(2, 8.0);
        
        // Build graph
        let graph = build_isoform_graph(&contigs, &overlaps, &expression_data);
        
        // Check nodes
        assert_eq!(graph.node_count(), 3);
        
        // Check edges
        assert_eq!(graph.edge_count(), 2);
        
        // Check edge weights
        let weight_0_1 = graph.edge_weight(0, 1).unwrap();
        let weight_1_2 = graph.edge_weight(1, 2).unwrap();
        
        assert!(*weight_0_1 > 0.8); // High confidence due to similar expression
        assert!(*weight_1_2 > 0.6); // Medium confidence
        
        // Check start/end nodes
        let start_nodes = find_start_nodes(&graph);
        let end_nodes = find_end_nodes(&graph);
        
        assert_eq!(start_nodes.len(), 1);
        assert_eq!(start_nodes[0], 0);
        
        assert_eq!(end_nodes.len(), 1);
        assert_eq!(end_nodes[0], 2);
    }
} 