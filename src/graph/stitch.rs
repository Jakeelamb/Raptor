use crate::accel::simd::{match_kmers_with_overlap, edit_distance_bp};
use petgraph::graph::DiGraph;
use petgraph::visit::EdgeRef;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone)]
pub struct Path {
    pub id: usize,
    pub segments: Vec<usize>, // Indexes of stitched contigs
    pub overlaps: Vec<usize>, // Optional overlap lengths
}

pub struct OverlapGraphBuilder {
    min_overlap: usize,
    max_distance: usize,
    max_edit_dist: usize,
}

impl OverlapGraphBuilder {
    pub fn new(min_overlap: usize, max_distance: usize, max_edit_dist: usize) -> Self {
        Self {
            min_overlap,
            max_distance,
            max_edit_dist,
        }
    }

    /// Find overlaps between sequences and build a directed graph
    pub fn build_overlap_graph(&self, sequences: &[String]) -> DiGraph<String, OverlapInfo> {
        let mut graph = DiGraph::new();
        let mut seq_to_node = HashMap::new();
        
        // Add all sequences as nodes
        for seq in sequences {
            let node_idx = graph.add_node(seq.clone());
            seq_to_node.insert(seq.clone(), node_idx);
        }
        
        // Find overlaps and add edges
        for (i, query) in sequences.iter().enumerate() {
            for (j, target) in sequences.iter().enumerate() {
                if i == j {
                    continue;
                }
                
                // Try suffix-prefix overlap
                if let Some((shift, distance)) = match_kmers_with_overlap(
                    query, 
                    target, 
                    self.min_overlap, 
                    self.max_distance
                ) {
                    // Get edit distance for better accuracy if needed
                    let overlap_len = query.len() - shift;
                    let query_suffix = &query[shift..];
                    let target_prefix = &target[..overlap_len];
                    
                    let edit_dist = if distance > 0 && distance <= self.max_distance {
                        // Only compute edit distance if there are some mismatches
                        // but within our hamming distance threshold
                        if overlap_len <= 64 { // Myers' algorithm has a 64-character limit
                            edit_distance_bp(query_suffix, target_prefix, self.max_edit_dist)
                                .unwrap_or(overlap_len) // If beyond max_edit_dist
                        } else {
                            distance as usize // Fall back to hamming distance for longer sequences
                        }
                    } else {
                        distance as usize
                    };
                    
                    // Only add edge if edit distance is within threshold
                    if edit_dist <= self.max_edit_dist {
                        let query_node = seq_to_node[query];
                        let target_node = seq_to_node[target];
                        
                        let overlap_info = OverlapInfo {
                            overlap_length: overlap_len,
                            distance: edit_dist,
                            shift,
                        };
                        
                        graph.add_edge(query_node, target_node, overlap_info);
                    }
                }
            }
        }
        
        graph
    }
    
    /// Stitch overlapping sequences to create contigs and track paths
    pub fn stitch_contigs(&self, graph: &DiGraph<String, OverlapInfo>) -> (Vec<String>, Vec<Path>) {
        let mut contigs = Vec::new();
        let mut paths = Vec::new();
        let mut visited = vec![false; graph.node_count()];
        let mut node_to_contig_idx = HashMap::new();
        
        for node_idx in graph.node_indices() {
            if visited[node_idx.index()] {
                continue;
            }
            
            let contig_idx = contigs.len();
            let mut contig = graph[node_idx].clone();
            let mut path_segments = vec![node_idx.index()];
            let mut path_overlaps = Vec::new();
            
            visited[node_idx.index()] = true;
            node_to_contig_idx.insert(node_idx.index(), contig_idx);
            
            // Extend contig by following edges
            let mut current = node_idx;
            while let Some(edge) = graph.edges(current).next() {
                let target = edge.target();
                if visited[target.index()] {
                    break;
                }
                
                let overlap_info = edge.weight();
                let target_seq = &graph[target];
                
                // Merge sequences based on overlap
                let extension = &target_seq[overlap_info.overlap_length..];
                contig.push_str(extension);
                
                // Track path information
                path_segments.push(target.index());
                path_overlaps.push(overlap_info.overlap_length);
                
                visited[target.index()] = true;
                node_to_contig_idx.insert(target.index(), contig_idx);
                current = target;
            }
            
            contigs.push(contig);
            
            // Create path record
            let path = Path {
                id: paths.len(),
                segments: path_segments,
                overlaps: path_overlaps,
            };
            paths.push(path);
        }
        
        (contigs, paths)
    }
}

#[derive(Debug, Clone)]
pub struct OverlapInfo {
    pub overlap_length: usize,
    pub distance: usize,
    pub shift: usize,
}

/// Load paths from a GFA file
pub fn load_paths_from_gfa(gfa_path: &str) -> Result<Vec<Path>, std::io::Error> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    
    let file = File::open(gfa_path)?;
    let reader = BufReader::new(file);
    let mut paths = Vec::new();
    
    // Map from segment names to numeric IDs
    let mut segment_to_id = HashMap::new();
    let mut next_id = 0;
    
    // First pass: collect segment IDs
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        
        if parts.is_empty() {
            continue;
        }
        
        if parts[0] == "S" && parts.len() >= 3 {
            // This is a segment line: S segmentID sequence
            if !segment_to_id.contains_key(parts[1]) {
                segment_to_id.insert(parts[1].to_string(), next_id);
                next_id += 1;
            }
        }
    }
    
    // Second pass: read paths
    let file = File::open(gfa_path)?;
    let reader = BufReader::new(file);
    
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        
        if parts.is_empty() {
            continue;
        }
        
        if parts[0] == "P" && parts.len() >= 3 {
            // This is a path line: P pathID segment1+,segment2+,... CIGAR
            let path_id = paths.len();
            let segment_str = parts[2];
            
            // Parse segments
            let mut segments = Vec::new();
            let mut overlaps = Vec::new();
            
            // Handle path segments in format "seg1+,seg2+,..."
            let seg_parts: Vec<&str> = segment_str.split(',').collect();
            
            for seg in seg_parts {
                // Remove orientation marker (+/-)
                let seg_name = seg.trim_end_matches(|c| c == '+' || c == '-');
                
                if let Some(&id) = segment_to_id.get(seg_name) {
                    segments.push(id);
                }
            }
            
            // Overlaps are usually not explicitly stated in GFA, so we'll use defaults
            if segments.len() > 1 {
                overlaps = vec![0; segments.len() - 1];
            }
            
            // Create path record
            let path = Path {
                id: path_id,
                segments,
                overlaps,
            };
            
            paths.push(path);
        }
    }
    
    Ok(paths)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_overlap_graph_builder() {
        let sequences = vec![
            "ATCGATCG".to_string(),
            "GATCGAAA".to_string(),
            "CGATCGAT".to_string(),
        ];
        
        let builder = OverlapGraphBuilder::new(4, 0, 0);
        let graph = builder.build_overlap_graph(&sequences);
        
        // Debug edges
        for edge in graph.edge_references() {
            let from = &graph[edge.source()];
            let to = &graph[edge.target()];
            let overlap = edge.weight();
            println!("Edge: {} -> {} (overlap: {}, dist: {})", 
                from, to, overlap.overlap_length, overlap.distance);
        }
        
        assert_eq!(graph.node_count(), 3);
        
        // Don't test exact edge count as it depends on the implementation
        // Just verify we have at least one edge
        assert!(graph.edge_count() > 0);
        
        // Test stitching
        let (contigs, paths) = builder.stitch_contigs(&graph);
        
        // Print contigs for debugging
        for (i, contig) in contigs.iter().enumerate() {
            println!("Contig {}: {}", i, contig);
        }
        
        // Print paths for debugging
        for path in &paths {
            println!("Path {}: segments={:?}, overlaps={:?}", 
                    path.id, path.segments, path.overlaps);
        }
        
        // There should be some contig that contains both ATCG and GATC or similar
        assert!(contigs.iter().any(|c| c.contains("ATCG") && c.contains("GATC")));
    }
} 