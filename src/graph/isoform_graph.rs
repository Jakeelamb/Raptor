use std::collections::HashMap;
use petgraph::graphmap::DiGraphMap;
use petgraph::Direction;
use tracing::debug;

/// Type definition for IsoformGraph - a directed graph where:
/// - Nodes are contig IDs (usize)
/// - Edge weights are confidence scores based on coverage (f32)
pub type IsoformGraph = DiGraphMap<usize, f32>;

/// Build a graph of potential isoforms from contigs and their overlaps
pub fn build_isoform_graph(
    contigs: &HashMap<usize, String>,
    overlaps: &Vec<(usize, usize, usize)>,
    expression_data: &HashMap<usize, f64>
) -> IsoformGraph {
    // Create a new directed graph
    let mut graph = DiGraphMap::new();
    
    // Add nodes for each contig
    for &contig_id in contigs.keys() {
        graph.add_node(contig_id);
    }
    
    debug!("Added {} nodes to isoform graph", graph.node_count());
    
    // Add edges for overlaps
    for &(from_id, to_id, overlap_len) in overlaps {
        // Calculate edge weight based on overlap length and expression
        let weight = calculate_edge_weight(from_id, to_id, overlap_len, expression_data);
        
        // Add edge
        graph.add_edge(from_id, to_id, weight);
    }
    
    debug!("Added {} edges to isoform graph", graph.edge_count());
    
    graph
}

/// Calculate edge weight based on overlap quality and expression
fn calculate_edge_weight(
    from_id: usize,
    to_id: usize,
    overlap_len: usize,
    expression_data: &HashMap<usize, f64>
) -> f32 {
    // Base weight from overlap length (normalized)
    let overlap_factor = (overlap_len as f32).min(100.0) / 100.0;
    
    // Expression similarity factor
    let expr_similarity = match (expression_data.get(&from_id), expression_data.get(&to_id)) {
        (Some(&from_expr), Some(&to_expr)) => {
            let min_expr = from_expr.min(to_expr);
            let max_expr = from_expr.max(to_expr);
            
            if max_expr == 0.0 {
                0.5 // Default if no expression
            } else {
                (min_expr / max_expr) as f32
            }
        },
        _ => 0.5, // Default if expression data is missing
    };
    
    // Combine factors - higher is better
    let weight = 0.7 * overlap_factor + 0.3 * expr_similarity;
    
    // Ensure weight is between 0.0 and 1.0
    weight.max(0.0).min(1.0)
}

/// Find potential transcript start nodes in the graph
pub fn find_start_nodes(graph: &IsoformGraph) -> Vec<usize> {
    let mut start_nodes = Vec::new();
    
    for node in graph.nodes() {
        // Start nodes have no incoming edges or many more outgoing than incoming
        let in_degree = graph.neighbors_directed(node, Direction::Incoming).count();
        let out_degree = graph.neighbors_directed(node, Direction::Outgoing).count();
        
        if in_degree == 0 || (out_degree > 0 && in_degree > 0 && out_degree >= 2 * in_degree) {
            start_nodes.push(node);
        }
    }
    
    debug!("Found {} potential start nodes", start_nodes.len());
    start_nodes
}

/// Find potential transcript end nodes in the graph
pub fn find_end_nodes(graph: &IsoformGraph) -> Vec<usize> {
    let mut end_nodes = Vec::new();
    
    for node in graph.nodes() {
        // End nodes have no outgoing edges or many more incoming than outgoing
        let in_degree = graph.neighbors_directed(node, Direction::Incoming).count();
        let out_degree = graph.neighbors_directed(node, Direction::Outgoing).count();
        
        if out_degree == 0 || (in_degree > 0 && out_degree > 0 && in_degree >= 2 * out_degree) {
            end_nodes.push(node);
        }
    }
    
    debug!("Found {} potential end nodes", end_nodes.len());
    end_nodes
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_build_simple_graph() {
        // Create test contigs
        let mut contigs = HashMap::new();
        contigs.insert(1, "ATCGATCG".to_string());
        contigs.insert(2, "CGATCGAT".to_string());
        contigs.insert(3, "TCGATCGA".to_string());
        
        // Create test overlaps
        let overlaps = vec![
            (1, 2, 5),
            (2, 3, 6),
        ];
        
        // Create test expression data
        let mut expression = HashMap::new();
        expression.insert(1, 10.0);
        expression.insert(2, 12.0);
        expression.insert(3, 8.0);
        
        // Build graph
        let graph = build_isoform_graph(&contigs, &overlaps, &expression);
        
        // Verify graph structure
        assert_eq!(graph.node_count(), 3);
        assert_eq!(graph.edge_count(), 2);
        
        // Verify start and end nodes
        let start_nodes = find_start_nodes(&graph);
        let end_nodes = find_end_nodes(&graph);
        
        assert_eq!(start_nodes.len(), 1);  // Node 1 should be the only start
        assert_eq!(end_nodes.len(), 1);    // Node 3 should be the only end
    }
    
    #[test]
    fn test_edge_weight_calculation() {
        // Test basic weight calculation
        let mut expression = HashMap::new();
        expression.insert(1, 100.0);
        expression.insert(2, 80.0);
        
        let weight = calculate_edge_weight(1, 2, 50, &expression);
        
        // Expected: 0.7 * (50/100) + 0.3 * (80/100) = 0.35 + 0.24 = 0.59
        assert!((weight - 0.59).abs() < 0.01);
        
        // Test with missing expression data
        let weight_missing = calculate_edge_weight(1, 3, 50, &expression);
        assert!((weight_missing - 0.5).abs() < 0.01);
    }
} 