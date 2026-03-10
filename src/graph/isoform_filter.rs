use crate::accel::simd::hamming_distance_simd;
use crate::graph::transcript::Transcript;
use crate::kmer::nthash::nthash;
use ahash::AHashSet;
use std::cmp::Ordering;
use std::collections::HashSet;
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

    if a_len == 0 {
        return b_len;
    }
    if b_len == 0 {
        return a_len;
    }

    let mut matrix = vec![vec![0; b_len + 1]; a_len + 1];

    for i in 0..=a_len {
        matrix[i][0] = i;
    }

    for j in 0..=b_len {
        matrix[0][j] = j;
    }

    for i in 1..=a_len {
        for j in 1..=b_len {
            let cost = if a.chars().nth(i - 1) == b.chars().nth(j - 1) {
                0
            } else {
                1
            };

            matrix[i][j] = std::cmp::min(
                std::cmp::min(
                    matrix[i - 1][j] + 1, // deletion
                    matrix[i][j - 1] + 1, // insertion
                ),
                matrix[i - 1][j - 1] + cost, // substitution
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
    if seq.len() < k {
        return AHashSet::new();
    }

    let mut hashes = AHashSet::with_capacity(seq.len() - k + 1);
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
    similarity_threshold: f64,
) -> Vec<usize> {
    let n = transcripts.len();
    let mut clusters = vec![0; n];

    for i in 0..n {
        clusters[i] = i; // Initially each transcript is in its own cluster
    }

    // Union-find approach to cluster similar transcripts
    for i in 0..n {
        for j in (i + 1)..n {
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

#[inline]
fn normalized_confidence(confidence: f64) -> f64 {
    if confidence.is_finite() {
        confidence
    } else {
        f64::NEG_INFINITY
    }
}

#[inline]
fn cmp_confidence_desc(a: f64, b: f64) -> Ordering {
    normalized_confidence(b).total_cmp(&normalized_confidence(a))
}

/// Filter out similar transcripts based on sequence similarity
pub fn filter_similar_transcripts(
    transcripts: &[Transcript],
    similarity_threshold: f64,
) -> Vec<Transcript> {
    if transcripts.is_empty() {
        return Vec::new();
    }

    debug!(
        "Filtering {} transcripts with similarity threshold {}",
        transcripts.len(),
        similarity_threshold
    );

    // Sort transcripts by length (descending) and confidence (descending)
    let mut sorted_transcripts = transcripts.to_vec();
    sorted_transcripts.sort_by(|a, b| {
        let len_cmp = b.length.cmp(&a.length);
        if len_cmp == std::cmp::Ordering::Equal {
            let confidence_cmp = cmp_confidence_desc(a.confidence, b.confidence);
            if confidence_cmp == Ordering::Equal {
                a.id.cmp(&b.id)
            } else {
                confidence_cmp
            }
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

    debug!(
        "Filtered out {} similar transcripts, kept {}",
        removed_ids.len(),
        filtered_transcripts.len()
    );

    filtered_transcripts
}

/// Merge similar transcripts instead of filtering them out
pub fn merge_transcripts(transcripts: &[Transcript], similarity_threshold: f64) -> Vec<Transcript> {
    if transcripts.is_empty() {
        return Vec::new();
    }

    debug!(
        "Merging {} transcripts with similarity threshold {}",
        transcripts.len(),
        similarity_threshold
    );

    // Sort transcripts by confidence (descending)
    let mut sorted_transcripts = transcripts.to_vec();
    sorted_transcripts.sort_by(|a, b| {
        let confidence_cmp = cmp_confidence_desc(a.confidence, b.confidence);
        if confidence_cmp != Ordering::Equal {
            return confidence_cmp;
        }

        let len_cmp = b.length.cmp(&a.length);
        if len_cmp != Ordering::Equal {
            return len_cmp;
        }

        a.id.cmp(&b.id)
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
            let mut confidence_sum = 0.0;
            let mut confidence_count = 0usize;
            for transcript in std::iter::once(&merged).chain(merged_with.iter().copied()) {
                if transcript.confidence.is_finite() {
                    confidence_sum += transcript.confidence;
                    confidence_count += 1;
                }
            }
            if confidence_count > 0 {
                merged.confidence = confidence_sum / confidence_count as f64;
            }

            // Update TPM if available, ignoring non-finite values.
            let mut total_tpm = 0.0;
            let mut has_finite_tpm = false;
            if let Some(tpm) = merged.tpm.filter(|v| v.is_finite()) {
                total_tpm += tpm;
                has_finite_tpm = true;
            }
            for transcript in &merged_with {
                if let Some(tpm) = transcript.tpm.filter(|v| v.is_finite()) {
                    total_tpm += tpm;
                    has_finite_tpm = true;
                }
            }
            merged.tpm = has_finite_tpm.then_some(total_tpm);

            debug!(
                "Merged transcript {} with {} others",
                merged.id,
                merged_with.len()
            );
        }

        merged_transcripts.push(merged);
    }

    debug!(
        "After merging, reduced from {} to {} transcripts",
        transcripts.len(),
        merged_transcripts.len()
    );

    merged_transcripts
}

/// Determine if two transcripts are similar based on sequence similarity
fn is_similar(a: &Transcript, b: &Transcript, threshold: f64) -> bool {
    // If length difference is too large, they're not similar
    let len_a = a.sequence.len();
    let len_b = b.sequence.len();
    if len_a == 0 || len_b == 0 {
        return len_a == len_b;
    }

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

    // Count shared unique k-mers (Jaccard-style).
    // Using unique sets avoids over-counting repeated k-mers in one sequence,
    // which could otherwise inflate similarity above 1.0.
    let shorter = if len_a <= len_b {
        &a.sequence
    } else {
        &b.sequence
    };
    let longer = if len_a > len_b {
        &a.sequence
    } else {
        &b.sequence
    };

    let shorter_bytes = shorter.as_bytes();
    let longer_bytes = longer.as_bytes();

    let mut shorter_kmers: AHashSet<&[u8]> = AHashSet::with_capacity(shorter_bytes.len() - k + 1);
    for i in 0..=shorter_bytes.len() - k {
        shorter_kmers.insert(&shorter_bytes[i..i + k]);
    }

    let mut longer_unique: AHashSet<&[u8]> = AHashSet::with_capacity(longer_bytes.len() - k + 1);
    let mut shared_unique_kmers = 0usize;
    for i in 0..=longer_bytes.len() - k {
        let kmer = &longer_bytes[i..i + k];
        if longer_unique.insert(kmer) && shorter_kmers.contains(kmer) {
            shared_unique_kmers += 1;
        }
    }

    let union_size = shorter_kmers.len() + longer_unique.len() - shared_unique_kmers;
    if union_size == 0 {
        return true;
    }
    let similarity = shared_unique_kmers as f64 / union_size as f64;

    similarity >= threshold
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[allow(deprecated)]
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
            splicing: "unknown".to_string(),
        };

        let t2 = Transcript {
            id: 2,
            sequence: "GCTTACA".to_string(),
            path: vec![1, 2, 3],
            confidence: 0.85,
            length: 7,
            strand: '+',
            tpm: None,
            splicing: "unknown".to_string(),
        };

        // Now using k-mer Jaccard similarity
        // For very short sequences (< JACCARD_K), falls back to length ratio
        let sim = calculate_sequence_similarity(&t1, &t2);
        // Both have length 7, so length-based similarity is 1.0
        assert!(
            sim > 0.0 && sim <= 1.0,
            "Similarity should be between 0 and 1, got {}",
            sim
        );

        let t3 = Transcript {
            id: 3,
            sequence: "GATTACA".to_string(),
            path: vec![1, 2, 3],
            confidence: 0.95,
            length: 7,
            strand: '+',
            tpm: None,
            splicing: "unknown".to_string(),
        };

        // Identical sequences should have high similarity
        let sim_identical = calculate_sequence_similarity(&t1, &t3);
        assert!(
            sim_identical >= 0.99,
            "Identical sequences should have similarity ~1.0, got {}",
            sim_identical
        );
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
            splicing: "unknown".to_string(),
        };

        let t2 = Transcript {
            id: 2,
            sequence: "GCTTACA".to_string(),
            path: vec![1, 2, 4],
            confidence: 0.85,
            length: 7,
            strand: '+',
            tpm: None,
            splicing: "unknown".to_string(),
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
            splicing: "unknown".to_string(),
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
            },
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

    #[test]
    fn test_filter_similar_is_deterministic_with_nan_and_tied_confidence() {
        let transcripts = vec![
            Transcript {
                id: 11,
                sequence: "AAAACCCCTTTT".to_string(),
                path: vec![1, 2, 3],
                confidence: f64::NAN,
                length: 12,
                strand: '+',
                tpm: Some(5.0),
                splicing: "linear".to_string(),
            },
            Transcript {
                id: 2,
                sequence: "AAAACCCCTTTT".to_string(),
                path: vec![1, 2, 3],
                confidence: 0.9,
                length: 12,
                strand: '+',
                tpm: Some(10.0),
                splicing: "linear".to_string(),
            },
            Transcript {
                id: 5,
                sequence: "AAAACCCCTTTT".to_string(),
                path: vec![1, 2, 3],
                confidence: 0.9,
                length: 12,
                strand: '+',
                tpm: Some(20.0),
                splicing: "linear".to_string(),
            },
        ];

        let mut reversed = transcripts.clone();
        reversed.reverse();

        let baseline = filter_similar_transcripts(&transcripts, 0.99);
        let reversed_result = filter_similar_transcripts(&reversed, 0.99);

        let baseline_ids: Vec<usize> = baseline.iter().map(|t| t.id).collect();
        let reversed_ids: Vec<usize> = reversed_result.iter().map(|t| t.id).collect();
        assert_eq!(baseline_ids, reversed_ids);
        assert_eq!(baseline_ids, vec![2]);
    }

    #[test]
    fn test_merge_transcripts_ignores_non_finite_values() {
        let transcripts = vec![
            Transcript {
                id: 1,
                sequence: "AAAACCCCTTTTGGGG".to_string(),
                path: vec![1, 2],
                confidence: 1.0,
                length: 16,
                strand: '+',
                tpm: Some(10.0),
                splicing: "linear".to_string(),
            },
            Transcript {
                id: 2,
                sequence: "AAAACCCCTTTTGGGG".to_string(),
                path: vec![1, 3],
                confidence: f64::INFINITY,
                length: 16,
                strand: '+',
                tpm: Some(f64::NAN),
                splicing: "linear".to_string(),
            },
            Transcript {
                id: 3,
                sequence: "AAAACCCCTTTTGGGG".to_string(),
                path: vec![1, 4],
                confidence: 0.5,
                length: 16,
                strand: '+',
                tpm: Some(2.5),
                splicing: "linear".to_string(),
            },
        ];

        let merged = merge_transcripts(&transcripts, 0.99);
        assert_eq!(merged.len(), 1);
        let only = &merged[0];
        assert!((only.confidence - 0.75).abs() < 1e-12);
        assert_eq!(only.tpm, Some(12.5));
    }

    #[test]
    fn test_is_similar_treats_empty_sequences_as_equal_only_when_both_empty() {
        let empty_a = Transcript {
            id: 1,
            sequence: String::new(),
            path: vec![],
            confidence: 1.0,
            length: 0,
            strand: '+',
            tpm: None,
            splicing: "linear".to_string(),
        };
        let empty_b = Transcript {
            id: 2,
            sequence: String::new(),
            path: vec![],
            confidence: 1.0,
            length: 0,
            strand: '+',
            tpm: None,
            splicing: "linear".to_string(),
        };
        let non_empty = Transcript {
            id: 3,
            sequence: "A".to_string(),
            path: vec![1],
            confidence: 1.0,
            length: 1,
            strand: '+',
            tpm: None,
            splicing: "linear".to_string(),
        };

        assert!(is_similar(&empty_a, &empty_b, 0.99));
        assert!(!is_similar(&empty_a, &non_empty, 0.99));
    }

    #[test]
    fn test_is_similar_length_mismatch_uses_unique_kmers_for_similarity() {
        let repetitive_short = Transcript {
            id: 1,
            sequence: format!("{}C", "A".repeat(49)),
            path: vec![1, 2],
            confidence: 1.0,
            length: 50,
            strand: '+',
            tpm: None,
            splicing: "linear".to_string(),
        };
        let repetitive_long = Transcript {
            id: 2,
            sequence: "A".repeat(75),
            path: vec![1, 3],
            confidence: 1.0,
            length: 75,
            strand: '+',
            tpm: None,
            splicing: "linear".to_string(),
        };

        // With k=25 this pair shares only one unique k-mer ("A"*25), so
        // similarity should be low (1 / 26), not inflated by repeated windows.
        assert!(!is_similar(&repetitive_short, &repetitive_long, 0.9));
        assert!(is_similar(&repetitive_short, &repetitive_long, 0.03));
    }
}
