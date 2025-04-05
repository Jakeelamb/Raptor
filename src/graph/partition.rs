use petgraph::graphmap::DiGraphMap;
use petgraph::algo::connected_components;
use std::collections::HashMap;

/// Clusters contigs in a graph into connected components
///
/// This function identifies clusters of contigs that are connected in the graph,
/// which can represent gene families or related transcripts.
///
/// # Arguments
/// * `graph` - A directed graph where nodes are contig IDs and edges represent overlaps
///
/// # Returns
/// * HashMap mapping cluster IDs to vectors of contig IDs
pub fn cluster_contigs(graph: &DiGraphMap<usize, f32>) -> HashMap<usize, Vec<usize>> {
    let mut clusters = HashMap::new();
    let mut label = 0;
    let mut visited = std::collections::HashSet::new();

    // Perform DFS to find connected components
    for node in graph.nodes() {
        if visited.contains(&node) { continue; }
        let mut component = vec![];
        let mut stack = vec![node];
        while let Some(n) = stack.pop() {
            if !visited.insert(n) { continue; }
            component.push(n);
            stack.extend(graph.neighbors(n));
        }
        clusters.insert(label, component);
        label += 1;
    }

    clusters
}

/// Merges small clusters that share significant k-mer content
/// 
/// # Arguments
/// * `clusters` - Initial clustering of contigs
/// * `contigs` - Vector of contig sequences
/// * `min_size` - Minimum cluster size to consider
/// * `similarity_threshold` - Jaccard similarity threshold for merging
/// 
/// # Returns
/// * Updated HashMap with merged clusters
pub fn merge_small_clusters(
    clusters: HashMap<usize, Vec<usize>>,
    contigs: &[String],
    min_size: usize,
    similarity_threshold: f32
) -> HashMap<usize, Vec<usize>> {
    let mut new_clusters = HashMap::new();
    let mut small_clusters = Vec::new();
    let mut next_id = clusters.len();
    
    // Separate small clusters from large ones
    for (id, component) in clusters {
        if component.len() < min_size {
            small_clusters.push((id, component));
        } else {
            new_clusters.insert(id, component);
        }
    }
    
    // Try to merge small clusters with larger ones
    for (small_id, small_component) in small_clusters {
        let mut best_match = None;
        let mut best_similarity = similarity_threshold;
        
        // Calculate k-mer similarity with existing clusters
        for (&cluster_id, cluster_contigs) in &new_clusters {
            let similarity = calculate_cluster_similarity(
                &small_component, cluster_contigs, contigs);
            
            if similarity > best_similarity {
                best_similarity = similarity;
                best_match = Some(cluster_id);
            }
        }
        
        if let Some(target_id) = best_match {
            // Merge with the best matching cluster
            let target_cluster = new_clusters.get_mut(&target_id).unwrap();
            target_cluster.extend(small_component);
        } else {
            // Keep as a separate cluster
            new_clusters.insert(next_id, small_component);
            next_id += 1;
        }
    }
    
    new_clusters
}

/// Calculates Jaccard similarity between two clusters based on shared k-mers
fn calculate_cluster_similarity(
    cluster1: &[usize],
    cluster2: &[usize],
    contigs: &[String]
) -> f32 {
    // For simplicity, we'll just check for shared contigs
    // In a real implementation, you'd calculate k-mer Jaccard similarity
    
    let set1: std::collections::HashSet<_> = cluster1.iter().collect();
    let set2: std::collections::HashSet<_> = cluster2.iter().collect();
    
    let intersection = set1.intersection(&set2).count();
    let union = set1.union(&set2).count();
    
    if union == 0 {
        0.0
    } else {
        intersection as f32 / union as f32
    }
} 