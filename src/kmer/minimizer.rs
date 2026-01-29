// src/kmer/minimizer.rs
//! Minimizer extraction for efficient sequence indexing.
//!
//! Minimizers reduce the number of indexed positions by 1/w while maintaining
//! sensitivity for overlap detection. This provides 10x less memory usage
//! compared to indexing all suffixes.
//!
//! Reference: Roberts, M., Hayes, W., Hunt, B. R., Mount, S. M., & Yorke, J. A. (2004).
//! Reducing storage requirements for biological sequence comparison.
//! Bioinformatics, 20(18), 3363-3369.

use crate::kmer::nthash::NtHashIterator;

/// A minimizer: a k-mer with minimum hash value in a window.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Minimizer {
    /// Hash value of the minimizer
    pub hash: u64,
    /// Position in the original sequence
    pub position: usize,
}

/// Extract minimizers from a sequence.
///
/// A minimizer is the k-mer with the minimum hash value in a sliding window
/// of w consecutive k-mers.
///
/// # Arguments
/// * `seq` - The DNA sequence as bytes
/// * `k` - K-mer size
/// * `w` - Window size (number of k-mers per window)
///
/// # Returns
/// Vector of (hash, position) minimizers
///
/// # Example
/// ```ignore
/// let minimizers = get_minimizers(b"ACGTACGTACGT", 4, 3);
/// // Returns minimizers for each window of 3 consecutive 4-mers
/// ```
pub fn get_minimizers(seq: &[u8], k: usize, w: usize) -> Vec<Minimizer> {
    if seq.len() < k || w == 0 {
        return Vec::new();
    }

    let mut minimizers = Vec::new();
    let mut window: Vec<(u64, usize)> = Vec::with_capacity(w);
    let mut last_minimizer: Option<(u64, usize)> = None;

    // Get all k-mer hashes
    for (pos, hash) in NtHashIterator::new(seq, k) {
        // Add to window
        window.push((hash, pos));

        // Remove old entries outside window
        if window.len() > w {
            window.remove(0);
        }

        // Only process once we have a full window
        if window.len() == w {
            // Find minimum in window
            let min_entry = window.iter().min_by_key(|(h, _)| h).copied();

            if let Some((min_hash, min_pos)) = min_entry {
                // Only add if this is a new minimizer
                if last_minimizer.map(|(_, p)| p) != Some(min_pos) {
                    minimizers.push(Minimizer {
                        hash: min_hash,
                        position: min_pos,
                    });
                    last_minimizer = Some((min_hash, min_pos));
                }
            }
        }
    }

    minimizers
}

/// Extract minimizers with their sequence context.
///
/// Returns minimizers along with the actual k-mer sequence for verification.
pub fn get_minimizers_with_seq<'a>(seq: &'a [u8], k: usize, w: usize) -> Vec<(Minimizer, &'a [u8])> {
    get_minimizers(seq, k, w)
        .into_iter()
        .filter(|m| m.position + k <= seq.len())
        .map(|m| (m, &seq[m.position..m.position + k]))
        .collect()
}

/// Build a minimizer index for a set of sequences.
///
/// Maps each minimizer hash to the list of (sequence_id, position) where it occurs.
pub struct MinimizerIndex {
    /// Maps hash -> Vec<(seq_id, position)>
    index: ahash::AHashMap<u64, Vec<(usize, usize)>>,
    /// K-mer size
    k: usize,
    /// Window size
    w: usize,
}

impl MinimizerIndex {
    /// Create a new minimizer index.
    pub fn new(k: usize, w: usize) -> Self {
        Self {
            index: ahash::AHashMap::new(),
            k,
            w,
        }
    }

    /// Create with pre-allocated capacity.
    pub fn with_capacity(k: usize, w: usize, capacity: usize) -> Self {
        Self {
            index: ahash::AHashMap::with_capacity(capacity),
            k,
            w,
        }
    }

    /// Add a sequence to the index.
    pub fn add_sequence(&mut self, seq_id: usize, seq: &[u8]) {
        for minimizer in get_minimizers(seq, self.k, self.w) {
            self.index
                .entry(minimizer.hash)
                .or_default()
                .push((seq_id, minimizer.position));
        }
    }

    /// Build index from multiple sequences.
    pub fn build(sequences: &[impl AsRef<[u8]>], k: usize, w: usize) -> Self {
        // Estimate capacity: roughly n/w minimizers per sequence
        let estimated_minimizers: usize = sequences
            .iter()
            .map(|s| s.as_ref().len().saturating_sub(k) / w.max(1) + 1)
            .sum();

        let mut index = Self::with_capacity(k, w, estimated_minimizers);

        for (seq_id, seq) in sequences.iter().enumerate() {
            index.add_sequence(seq_id, seq.as_ref());
        }

        index
    }

    /// Query the index for sequences sharing minimizers with the query.
    ///
    /// Returns a map from sequence_id to list of matching minimizer positions.
    pub fn query(&self, query: &[u8]) -> ahash::AHashMap<usize, Vec<(usize, usize)>> {
        let mut hits: ahash::AHashMap<usize, Vec<(usize, usize)>> = ahash::AHashMap::new();

        for minimizer in get_minimizers(query, self.k, self.w) {
            if let Some(matches) = self.index.get(&minimizer.hash) {
                for &(seq_id, target_pos) in matches {
                    hits.entry(seq_id)
                        .or_default()
                        .push((minimizer.position, target_pos));
                }
            }
        }

        hits
    }

    /// Find candidate overlaps between the query and indexed sequences.
    ///
    /// Returns (seq_id, estimated_overlap_quality) pairs sorted by quality.
    pub fn find_overlap_candidates(
        &self,
        query: &[u8],
        min_shared_minimizers: usize,
    ) -> Vec<(usize, usize)> {
        let hits = self.query(query);

        let mut candidates: Vec<(usize, usize)> = hits
            .into_iter()
            .filter(|(_, matches)| matches.len() >= min_shared_minimizers)
            .map(|(seq_id, matches)| (seq_id, matches.len()))
            .collect();

        // Sort by number of shared minimizers (descending)
        candidates.sort_by(|a, b| b.1.cmp(&a.1));

        candidates
    }

    /// Get the number of unique minimizers in the index.
    pub fn num_minimizers(&self) -> usize {
        self.index.len()
    }

    /// Get the total number of (minimizer, position) entries.
    pub fn num_entries(&self) -> usize {
        self.index.values().map(|v| v.len()).sum()
    }

    /// Get k-mer size.
    pub fn k(&self) -> usize {
        self.k
    }

    /// Get window size.
    pub fn w(&self) -> usize {
        self.w
    }
}

/// Compute the density of minimizers (ratio of minimizers to k-mers).
///
/// Theoretical density is approximately 2/(w+1) for random sequences.
pub fn minimizer_density(seq: &[u8], k: usize, w: usize) -> f64 {
    let num_kmers = seq.len().saturating_sub(k - 1);
    if num_kmers == 0 {
        return 0.0;
    }

    let num_minimizers = get_minimizers(seq, k, w).len();
    num_minimizers as f64 / num_kmers as f64
}

/// Syncmer extraction (an alternative to minimizers with more even spacing).
///
/// A syncmer is a k-mer where the minimum s-mer occurs at a specific position
/// (e.g., the first or last s-mer position).
///
/// Reference: Edgar, R. (2021). Syncmers are more sensitive than minimizers
/// for selecting conserved k-mers in biological sequences. PeerJ, 9, e10805.
pub fn get_syncmers(seq: &[u8], k: usize, s: usize, t: usize) -> Vec<Minimizer> {
    if seq.len() < k || s > k || t >= k - s + 1 {
        return Vec::new();
    }

    let mut syncmers = Vec::new();

    // Get all k-mer hashes
    for (kmer_pos, kmer_hash) in NtHashIterator::new(seq, k) {
        // Extract s-mers within this k-mer and find minimum position
        let kmer_seq = &seq[kmer_pos..kmer_pos + k];

        let mut min_hash = u64::MAX;
        let mut min_pos = 0;

        for (smer_pos, smer_hash) in NtHashIterator::new(kmer_seq, s) {
            if smer_hash < min_hash {
                min_hash = smer_hash;
                min_pos = smer_pos;
            }
        }

        // Check if minimum s-mer is at target position
        if min_pos == t {
            syncmers.push(Minimizer {
                hash: kmer_hash,
                position: kmer_pos,
            });
        }
    }

    syncmers
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_minimizers_basic() {
        let seq = b"ACGTACGTACGTACGT";
        let minimizers = get_minimizers(seq, 4, 3);

        // Should get some minimizers
        assert!(!minimizers.is_empty());

        // Each minimizer position should be valid
        for m in &minimizers {
            assert!(m.position + 4 <= seq.len());
        }
    }

    #[test]
    fn test_minimizer_density() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let density = minimizer_density(seq, 8, 5);

        // Density should be roughly 2/(w+1) = 2/6 ≈ 0.33
        assert!(density > 0.1 && density < 0.8, "Density: {}", density);
    }

    #[test]
    fn test_minimizer_index() {
        let sequences = vec![
            b"ACGTACGTACGT".to_vec(),
            b"ACGTACGTTTTTT".to_vec(),
            b"GGGGACGTACGT".to_vec(),
        ];

        let index = MinimizerIndex::build(&sequences, 4, 2);

        // Query with first sequence
        let candidates = index.find_overlap_candidates(b"ACGTACGT", 1);

        // Should find at least the first sequence
        assert!(!candidates.is_empty());
    }

    #[test]
    fn test_syncmers() {
        let seq = b"ACGTACGTACGTACGT";
        let syncmers = get_syncmers(seq, 8, 4, 0); // Open syncmers (t=0)

        // Should get some syncmers
        // Not all k-mers will be syncmers
        assert!(syncmers.len() <= seq.len() - 8 + 1);

        // Each position should be valid
        for s in &syncmers {
            assert!(s.position + 8 <= seq.len());
        }
    }

    #[test]
    fn test_empty_sequence() {
        assert!(get_minimizers(b"", 4, 3).is_empty());
        assert!(get_minimizers(b"ACG", 4, 3).is_empty()); // Too short
    }
}
