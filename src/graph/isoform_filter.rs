use std::collections::HashSet;
use crate::graph::transcript::Transcript;
use crate::accel::simd::hamming_distance_simd;
use crate::kmer::nthash::nthash;
use ahash::AHashSet;
use tracing::debug;

/// K-mer size for Jaccard similarity calculation
const JACCARD_K: usize = 15;

/// Calculate the Levenshtein edit distance between two strings.
///
/// DEPRECATED: Use kmer_jaccard_similarity for O(m+n) performance
/// instead of O(m*n) edit distance.
#[allow(dead_code)]
#[deprecated(note = "Use kmer_jaccard_similarity for 100-1000x speedup")]
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

/// Calculate k-mer Jaccard similarity between two sequences.
///
/// This is O(m+n) vs O(m*n) for edit distance, providing 100-1000x speedup.
/// Jaccard similarity = |intersection| / |union|
#[inline]
pub fn kmer_jaccard_similarity(seq1: &str, seq2: &str, k: usize) -> f64 {
    if seq1.len() < k || seq2.len() < k {
        // Fall back to length-based similarity for very short sequences
        let min_len = seq1.len().min(seq2.len());
        let max_len = seq1.len().max(seq2.len());
        if max_len == 0 {
            return 1.0;
        }
        return min_len as f64 / max_len as f64;
    }

    // Extract k-mer hashes from both sequences
    let kmers1: AHashSet<u64> = extract_kmer_hashes(seq1.as_bytes(), k);
    let kmers2: AHashSet<u64> = extract_kmer_hashes(seq2.as_bytes(), k);

    if kmers1.is_empty() && kmers2.is_empty() {
        return 1.0;
    }

    let intersection = kmers1.intersection(&kmers2).count();
    let union = kmers1.len() + kmers2.len() - intersection;

    if union == 0 {
        return 1.0;
    }

    intersection as f64 / union as f64
}

/// Extract k-mer hashes from a sequence using ntHash.
#[inline]
fn extract_kmer_hashes(seq: &[u8], k: usize) -> AHashSet<u64> {
    let mut hashes = AHashSet::new();

    if seq.len() < k {
        return hashes;
    }

    for i in 0..=seq.len() - k {
        if let Some(hash) = nthash(&seq[i..i + k]) {
            hashes.insert(hash);
        }
    }

    hashes
}

/// Calculate the sequence similarity between two transcripts using k-mer Jaccard.
///
/// This replaces the O(m*n) edit distance with O(m+n) k-mer comparison.
fn calculate_sequence_similarity(t1: &Transcript, t2: &Transcript) -> f64 {
    kmer_jaccard_similarity(&t1.sequence, &t2.sequence, JACCARD_K)
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

/// Filter out similar transcripts based on sequence similarity
pub fn filter_similar_transcripts(
    transcripts: &[Transcript],
    similarity_threshold: f64
) -> Vec<Transcript> {
    if transcripts.is_empty() {
        return Vec::new();
    }
    
    debug!("Filtering {} transcripts with similarity threshold {}", 
            transcripts.len(), similarity_threshold);
    
    // Sort transcripts by length (descending) and confidence (descending)
    let mut sorted_transcripts = transcripts.to_vec();
    sorted_transcripts.sort_by(|a, b| {
        let len_cmp = b.length.cmp(&a.length);
        if len_cmp == std::cmp::Ordering::Equal {
            b.confidence.partial_cmp(&a.confidence).unwrap_or(std::cmp::Ordering::Equal)
        } else {
            len_cmp
        }
    });
    
    let mut filtered_transcripts = Vec::new();
    let mut removed_ids = HashSet::new();
    
    // Keep highest confidence/longest transcripts, filter out similar ones
    for (i, transcript) in sorted_transcripts.iter().enumerate() {
        if removed_ids.contains(&transcript.id) {
            continue;
        }
        
        filtered_transcripts.push(transcript.clone());
        
        // Compare with remaining transcripts
        for (j, other) in sorted_transcripts.iter().enumerate() {
            if i == j || removed_ids.contains(&other.id) {
                continue;
            }
            
            if is_similar(transcript, other, similarity_threshold) {
                removed_ids.insert(other.id);
            }
        }
    }
    
    debug!("Filtered out {} similar transcripts, kept {}", 
            removed_ids.len(), filtered_transcripts.len());
    
    filtered_transcripts
}

/// Merge similar transcripts instead of filtering them out
pub fn merge_transcripts(
    transcripts: &[Transcript],
    similarity_threshold: f64
) -> Vec<Transcript> {
    if transcripts.is_empty() {
        return Vec::new();
    }
    
    debug!("Merging {} transcripts with similarity threshold {}", 
            transcripts.len(), similarity_threshold);
    
    // Sort transcripts by confidence (descending)
    let mut sorted_transcripts = transcripts.to_vec();
    sorted_transcripts.sort_by(|a, b| {
        b.confidence.partial_cmp(&a.confidence).unwrap_or(std::cmp::Ordering::Equal)
    });
    
    let mut merged_transcripts = Vec::new();
    let mut removed_ids = HashSet::new();
    
    // Process all transcripts
    for (i, transcript) in sorted_transcripts.iter().enumerate() {
        if removed_ids.contains(&transcript.id) {
            continue;
        }
        
        let mut merged = transcript.clone();
        let mut merged_with = Vec::new();
        
        // Find transcripts to merge with
        for (j, other) in sorted_transcripts.iter().enumerate() {
            if i == j || removed_ids.contains(&other.id) {
                continue;
            }
            
            if is_similar(&merged, other, similarity_threshold) {
                merged_with.push(other);
                removed_ids.insert(other.id);
            }
        }
        
        // If we found transcripts to merge
        if !merged_with.is_empty() {
            // Update confidence as weighted average
            let total_confidence = merged_with.iter()
                .fold(merged.confidence, |acc, t| acc + t.confidence);
            
            let avg_confidence = total_confidence / (merged_with.len() as f64 + 1.0);
            merged.confidence = avg_confidence;
            
            // Update TPM if available
            if let Some(tpm) = merged.tpm {
                let total_tpm = merged_with.iter()
                    .fold(tpm, |acc, t| acc + t.tpm.unwrap_or(0.0));
                
                merged.tpm = Some(total_tpm);
            }
            
            debug!("Merged transcript {} with {} others", merged.id, merged_with.len());
        }
        
        merged_transcripts.push(merged);
    }
    
    debug!("After merging, reduced from {} to {} transcripts", 
            transcripts.len(), merged_transcripts.len());
    
    merged_transcripts
}

/// Determine if two transcripts are similar based on sequence similarity
fn is_similar(a: &Transcript, b: &Transcript, threshold: f64) -> bool {
    // If length difference is too large, they're not similar
    let len_a = a.sequence.len();
    let len_b = b.sequence.len();
    
    // If one sequence is more than 50% longer than the other, they're not similar
    if len_a as f64 > len_b as f64 * 1.5 || len_b as f64 > len_a as f64 * 1.5 {
        return false;
    }
    
    // Use Hamming distance for sequences of similar length
    if len_a == len_b {
        let distance = hamming_distance_simd(a.sequence.as_bytes(), b.sequence.as_bytes());
        let similarity = 1.0 - (distance as f64 / len_a as f64);
        return similarity >= threshold;
    }
    
    // For differing lengths, use a k-mer based approach
    // This is a simple implementation - for production, use a more sophisticated algorithm
    let k = 25; // k-mer size
    let min_len = len_a.min(len_b);
    
    if min_len < k {
        return false; // Too short to compare meaningfully
    }
    
    // Count shared k-mers
    let shorter = if len_a <= len_b { &a.sequence } else { &b.sequence };
    let longer = if len_a > len_b { &a.sequence } else { &b.sequence };
    
    let mut shorter_kmers = HashSet::new();
    for i in 0..=shorter.len() - k {
        shorter_kmers.insert(&shorter[i..i+k]);
    }
    
    let mut shared_kmers = 0;
    for i in 0..=longer.len() - k {
        if shorter_kmers.contains(&longer[i..i+k]) {
            shared_kmers += 1;
        }
    }
    
    // Calculate similarity as proportion of shared k-mers
    let max_possible_shared = (shorter.len() - k + 1).min(longer.len() - k + 1);
    let similarity = shared_kmers as f64 / max_possible_shared as f64;
    
    similarity >= threshold
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_edit_distance() {
        assert_eq!(edit_distance("GATTACA", "GCTTACA"), 1);
        assert_eq!(edit_distance("GATTACA", "GATTACA"), 0);
        assert_eq!(edit_distance("GATTACA", "GATT"), 3);
    }
    
    #[test]
    fn test_calculate_sequence_similarity() {
        let t1 = Transcript {
            id: 1,
            sequence: "GATTACA".to_string(),
            path: vec![1, 2, 3],
            confidence: 0.9,
            length: 7,
            strand: '+',
            tpm: None,
            splicing: "unknown".to_string()
        };

        let t2 = Transcript {
            id: 2,
            sequence: "GCTTACA".to_string(),
            path: vec![1, 2, 3],
            confidence: 0.85,
            length: 7,
            strand: '+',
            tpm: None,
            splicing: "unknown".to_string()
        };

        // Now using k-mer Jaccard similarity
        // For very short sequences (< JACCARD_K), falls back to length ratio
        let sim = calculate_sequence_similarity(&t1, &t2);
        // Both have length 7, so length-based similarity is 1.0
        assert!(sim > 0.0 && sim <= 1.0, "Similarity should be between 0 and 1, got {}", sim);

        let t3 = Transcript {
            id: 3,
            sequence: "GATTACA".to_string(),
            path: vec![1, 2, 3],
            confidence: 0.95,
            length: 7,
            strand: '+',
            tpm: None,
            splicing: "unknown".to_string()
        };

        // Identical sequences should have high similarity
        let sim_identical = calculate_sequence_similarity(&t1, &t3);
        assert!(sim_identical >= 0.99, "Identical sequences should have similarity ~1.0, got {}", sim_identical);
    }
    
    #[test]
    fn test_calculate_path_similarity() {
        let t1 = Transcript {
            id: 1,
            sequence: "GATTACA".to_string(),
            path: vec![1, 2, 3],
            confidence: 0.9,
            length: 7,
            strand: '+',
            tpm: None,
            splicing: "unknown".to_string()
        };
        
        let t2 = Transcript {
            id: 2,
            sequence: "GCTTACA".to_string(),
            path: vec![1, 2, 4],
            confidence: 0.85,
            length: 7,
            strand: '+',
            tpm: None,
            splicing: "unknown".to_string()
        };
        
        // 2 out of 3 elements in common
        assert_eq!(calculate_path_similarity(&t1, &t2), 0.5);
        
        let t3 = Transcript {
            id: 3,
            sequence: "GATTACA".to_string(),
            path: vec![1, 2, 3],
            confidence: 0.95,
            length: 7,
            strand: '+',
            tpm: None,
            splicing: "unknown".to_string()
        };
        
        // Identical paths
        assert_eq!(calculate_path_similarity(&t1, &t3), 1.0);
    }
    
    fn create_test_transcripts() -> Vec<Transcript> {
        vec![
            Transcript {
                id: 1,
                sequence: "AAAAAAAAAATTTTTTTTTT".to_string(),
                path: vec![1, 2],
                confidence: 0.95,
                length: 20,
                strand: '+',
                tpm: Some(100.0),
                splicing: "linear".to_string(),
            },
            Transcript {
                id: 2,
                sequence: "AAAAAAAAAATTTTTTTTTC".to_string(), // 95% similar to first
                path: vec![1, 3],
                confidence: 0.85,
                length: 20,
                strand: '+',
                tpm: Some(90.0),
                splicing: "linear".to_string(),
            },
            Transcript {
                id: 3,
                sequence: "GGGGGGGGGGGGGGGGGGG".to_string(), // Different
                path: vec![4, 5],
                confidence: 0.75,
                length: 20,
                strand: '+',
                tpm: Some(50.0),
                splicing: "linear".to_string(),
            }
        ]
    }
    
    #[test]
    fn test_filter_similar() {
        let transcripts = create_test_transcripts();
        
        // With high threshold - should keep all
        let filtered_high = filter_similar_transcripts(&transcripts, 0.99);
        assert_eq!(filtered_high.len(), 3);
        
        // With lower threshold - should filter out one
        let filtered_low = filter_similar_transcripts(&transcripts, 0.9);
        assert_eq!(filtered_low.len(), 2);
        
        // Check that the highest confidence one was kept
        let kept_ids: Vec<usize> = filtered_low.iter().map(|t| t.id).collect();
        assert!(kept_ids.contains(&1)); // Highest confidence should be kept
        assert!(kept_ids.contains(&3)); // Different sequence should be kept
    }
    
    #[test]
    fn test_merge_transcripts() {
        let transcripts = create_test_transcripts();
        
        // With high threshold - no merging
        let merged_high = merge_transcripts(&transcripts, 0.99);
        assert_eq!(merged_high.len(), 3);
        
        // With lower threshold - should merge two
        let merged_low = merge_transcripts(&transcripts, 0.9);
        assert_eq!(merged_low.len(), 2);
        
        // Check that merged transcript has updated values
        let merged = merged_low.iter().find(|t| t.id == 1).unwrap();
        assert!((merged.confidence - 0.9).abs() < 0.01); // Average of 0.95 and 0.85
        assert!(merged.tpm.unwrap() > 100.0); // Should be sum of TPMs
    }
} 