use std::collections::{HashMap, HashSet};
use crate::graph::transcript::Transcript;

/// Calculate the Levenshtein edit distance between two strings
fn edit_distance(a: &str, b: &str) -> usize {
    let a_len = a.len();
    let b_len = b.len();
    
    if a_len == 0 { return b_len; }
    if b_len == 0 { return a_len; }
    
    let mut matrix = vec![vec![0; b_len + 1]; a_len + 1];
    
    for i in 0..=a_len {
        matrix[i][0] = i;
    }
    
    for j in 0..=b_len {
        matrix[0][j] = j;
    }
    
    for i in 1..=a_len {
        for j in 1..=b_len {
            let cost = if a.chars().nth(i-1) == b.chars().nth(j-1) { 0 } else { 1 };
            
            matrix[i][j] = std::cmp::min(
                std::cmp::min(
                    matrix[i-1][j] + 1,     // deletion
                    matrix[i][j-1] + 1      // insertion
                ),
                matrix[i-1][j-1] + cost     // substitution
            );
        }
    }
    
    matrix[a_len][b_len]
}

/// Calculate the sequence similarity between two transcripts as a ratio
fn calculate_sequence_similarity(t1: &Transcript, t2: &Transcript) -> f64 {
    let dist = edit_distance(&t1.sequence, &t2.sequence);
    let max_len = std::cmp::max(t1.sequence.len(), t2.sequence.len());
    
    if max_len == 0 {
        return 1.0;
    }
    
    1.0 - (dist as f64 / max_len as f64)
}

/// Calculate the path similarity between two transcripts as a ratio
fn calculate_path_similarity(t1: &Transcript, t2: &Transcript) -> f64 {
    let t1_nodes: HashSet<usize> = t1.path.iter().cloned().collect();
    let t2_nodes: HashSet<usize> = t2.path.iter().cloned().collect();
    
    let intersection_size = t1_nodes.intersection(&t2_nodes).count();
    let union_size = t1_nodes.union(&t2_nodes).count();
    
    if union_size == 0 {
        return 1.0;
    }
    
    intersection_size as f64 / union_size as f64
}

/// Calculate the overall similarity between two transcripts
/// Uses a weighted combination of sequence similarity and path similarity
pub fn calculate_transcript_similarity(t1: &Transcript, t2: &Transcript) -> f64 {
    let seq_sim = calculate_sequence_similarity(t1, t2);
    let path_sim = calculate_path_similarity(t1, t2);
    
    // Weight sequence similarity higher than path similarity
    0.7 * seq_sim + 0.3 * path_sim
}

/// Cluster similar transcripts based on similarity threshold
pub fn cluster_similar_transcripts(
    transcripts: &[Transcript], 
    similarity_threshold: f64
) -> Vec<usize> {
    let n = transcripts.len();
    let mut clusters = vec![0; n];
    
    for i in 0..n {
        clusters[i] = i; // Initially each transcript is in its own cluster
    }
    
    // Union-find approach to cluster similar transcripts
    for i in 0..n {
        for j in (i+1)..n {
            let similarity = calculate_transcript_similarity(&transcripts[i], &transcripts[j]);
            
            if similarity >= similarity_threshold {
                // Find the root clusters
                let mut root_i = i;
                while clusters[root_i] != root_i {
                    root_i = clusters[root_i];
                }
                
                let mut root_j = j;
                while clusters[root_j] != root_j {
                    root_j = clusters[root_j];
                }
                
                // Merge clusters by making the smaller index the parent
                if root_i < root_j {
                    clusters[root_j] = root_i;
                } else {
                    clusters[root_i] = root_j;
                }
            }
        }
    }
    
    // Flatten clusters to their roots
    for i in 0..n {
        let mut root = i;
        while clusters[root] != root {
            root = clusters[root];
        }
        clusters[i] = root;
    }
    
    clusters
}

/// Filter similar transcripts to reduce redundancy
pub fn filter_similar_transcripts(
    transcripts: &[Transcript], 
    similarity_threshold: f64
) -> Vec<Transcript> {
    let clusters = cluster_similar_transcripts(transcripts, similarity_threshold);
    
    // Select one representative transcript from each cluster (the one with highest confidence)
    let mut cluster_reps: HashMap<usize, usize> = HashMap::new();
    
    for (i, &cluster) in clusters.iter().enumerate() {
        cluster_reps
            .entry(cluster)
            .and_modify(|best_idx| {
                if transcripts[i].confidence > transcripts[*best_idx].confidence {
                    *best_idx = i;
                }
            })
            .or_insert(i);
    }
    
    // Collect the representative transcripts
    cluster_reps.values()
        .map(|&idx| transcripts[idx].clone())
        .collect()
}

/// Merge similar transcripts into combined transcripts
pub fn merge_transcripts(
    transcripts: &[Transcript], 
    similarity_threshold: f64
) -> Vec<Transcript> {
    let clusters = cluster_similar_transcripts(transcripts, similarity_threshold);
    
    // Group transcripts by cluster
    let mut cluster_map: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, &cluster) in clusters.iter().enumerate() {
        cluster_map.entry(cluster).or_default().push(i);
    }
    
    // Merge transcripts in each cluster
    let mut merged_transcripts = Vec::new();
    for (cluster_id, indices) in cluster_map {
        if indices.len() == 1 {
            // Single transcript in cluster, no merging needed
            merged_transcripts.push(transcripts[indices[0]].clone());
            continue;
        }
        
        // Calculate representative sequence and path
        // For simplicity, choose the transcript with highest confidence
        let mut best_idx = indices[0];
        let mut best_confidence = transcripts[best_idx].confidence;
        
        for &idx in &indices[1..] {
            if transcripts[idx].confidence > best_confidence {
                best_idx = idx;
                best_confidence = transcripts[idx].confidence;
            }
        }
        
        // Create a merged transcript using the best transcript as base
        let best_transcript = &transcripts[best_idx];
        
        // Calculate merged confidence as average of all transcripts in the cluster
        let avg_confidence = indices.iter()
            .map(|&idx| transcripts[idx].confidence)
            .sum::<f64>() / indices.len() as f64;
        
        // Create merged transcript
        let mut merged = best_transcript.clone();
        merged.confidence = avg_confidence;
        merged.id = cluster_id; // Use cluster ID as merged transcript ID
        
        merged_transcripts.push(merged);
    }
    
    merged_transcripts
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_edit_distance() {
        assert_eq!(edit_distance("kitten", "sitting"), 3);
        assert_eq!(edit_distance("GATTACA", "GCATGCU"), 4);
        assert_eq!(edit_distance("", "abc"), 3);
        assert_eq!(edit_distance("abc", ""), 3);
    }
    
    #[test]
    fn test_calculate_sequence_similarity() {
        let t1 = Transcript {
            id: 1,
            sequence: "GATTACA".to_string(),
            path: vec![1, 2, 3],
            confidence: 0.9,
            length: 7,
        };
        
        let t2 = Transcript {
            id: 2,
            sequence: "GCTTACA".to_string(),
            path: vec![1, 2, 3],
            confidence: 0.8,
            length: 7,
        };
        
        // One substitution in a 7-char string = 1-1/7 = ~0.857 similarity
        assert!((calculate_sequence_similarity(&t1, &t2) - 0.857).abs() < 0.001);
        
        let t3 = Transcript {
            id: 3,
            sequence: "GATTACA".to_string(),
            path: vec![1, 2, 3],
            confidence: 0.7,
            length: 7,
        };
        
        // Identical sequences = 1.0 similarity
        assert_eq!(calculate_sequence_similarity(&t1, &t3), 1.0);
    }
    
    #[test]
    fn test_calculate_path_similarity() {
        let t1 = Transcript {
            id: 1,
            sequence: "GATTACA".to_string(),
            path: vec![1, 2, 3],
            confidence: 0.9,
            length: 7,
        };
        
        let t2 = Transcript {
            id: 2,
            sequence: "GCTTACA".to_string(),
            path: vec![1, 2, 4],
            confidence: 0.8,
            length: 7,
        };
        
        // 2 common nodes out of 4 unique nodes = 0.5 similarity
        assert_eq!(calculate_path_similarity(&t1, &t2), 0.5);
        
        let t3 = Transcript {
            id: 3,
            sequence: "GATTACA".to_string(),
            path: vec![1, 2, 3],
            confidence: 0.7,
            length: 7,
        };
        
        // Identical paths = 1.0 similarity
        assert_eq!(calculate_path_similarity(&t1, &t3), 1.0);
    }
    
    #[test]
    fn test_filter_similar_transcripts() {
        let t1 = Transcript {
            id: 1,
            sequence: "GATTACA".to_string(),
            path: vec![1, 2, 3],
            confidence: 0.9,
            length: 7,
        };
        
        let t2 = Transcript {
            id: 2,
            sequence: "GCTTACA".to_string(), // 1 diff
            path: vec![1, 2, 4],
            confidence: 0.8,
            length: 7,
        };
        
        let t3 = Transcript {
            id: 3,
            sequence: "GATTACA".to_string(), // identical to t1
            path: vec![1, 2, 3],
            confidence: 0.7,
            length: 7,
        };
        
        let transcripts = vec![t1.clone(), t2.clone(), t3.clone()];
        
        // High threshold - all transcripts should be kept
        let filtered_high = filter_similar_transcripts(&transcripts, 0.95);
        assert_eq!(filtered_high.len(), 2);
        
        // Low threshold - only one transcript should remain
        let filtered_low = filter_similar_transcripts(&transcripts, 0.5);
        assert_eq!(filtered_low.len(), 1);
        
        // The remaining transcript should be t1 (highest confidence among similar ones)
        assert_eq!(filtered_low[0].id, 1);
    }
    
    #[test]
    fn test_merge_transcripts() {
        let t1 = Transcript {
            id: 1,
            sequence: "GATTACA".to_string(),
            path: vec![1, 2, 3],
            confidence: 0.9,
            length: 7,
        };
        
        let t2 = Transcript {
            id: 2,
            sequence: "GCTTACA".to_string(), // 1 diff
            path: vec![1, 2, 4],
            confidence: 0.8,
            length: 7,
        };
        
        let t3 = Transcript {
            id: 3,
            sequence: "GATTACA".to_string(), // identical to t1
            path: vec![1, 2, 3],
            confidence: 0.7,
            length: 7,
        };
        
        let transcripts = vec![t1.clone(), t2.clone(), t3.clone()];
        
        // High threshold - less merging
        let merged_high = merge_transcripts(&transcripts, 0.95);
        assert_eq!(merged_high.len(), 2);
        
        // Low threshold - more merging
        let merged_low = merge_transcripts(&transcripts, 0.5);
        assert_eq!(merged_low.len(), 1);
        
        // The merged transcript should have averaged confidence
        assert!((merged_low[0].confidence - 0.8).abs() < 0.001);
    }
} 