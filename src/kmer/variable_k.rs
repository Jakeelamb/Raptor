use std::collections::HashMap;
use ahash::AHashSet;
use crate::kmer::nthash::NtHashIterator;

/// Generate a histogram of unique k-mer counts for different k values.
/// Uses ntHash for fast O(1) rolling hash computation.
pub fn kmer_coverage_histogram(sequences: &[String], max_k: usize) -> HashMap<usize, usize> {
    let mut counts: HashMap<usize, usize> = HashMap::new();

    // Choose k values to test up to max_k
    let k_values: Vec<usize> = if max_k <= 21 {
        vec![15, 17, 19, 21]
    } else if max_k <= 31 {
        vec![21, 25, 27, 31]
    } else {
        vec![21, 25, 31, 35, 41].into_iter().filter(|&k| k <= max_k).collect()
    };

    for &k in &k_values {
        // Use u64 hash set for memory efficiency
        let mut seen: AHashSet<u64> = AHashSet::new();

        for seq in sequences {
            let bytes = seq.as_bytes();
            // Use ntHash for O(1) rolling hash per k-mer
            for (_, hash) in NtHashIterator::new(bytes, k) {
                seen.insert(hash);
            }
        }

        counts.insert(k, seen.len());
    }

    counts
}

/// Selects the best k value based on the k-mer histogram
/// This uses a heuristic that balances uniqueness with k-mer length
pub fn select_best_k(hist: &HashMap<usize, usize>) -> usize {
    if hist.is_empty() {
        return 31; // Default if no data
    }
    
    // Apply a simple heuristic: maximize uniqueness but penalize very large k values
    let best_k = hist.iter()
        .max_by_key(|(&k, &count)| {
            // Balance between uniqueness and reasonable k-size
            // For very large k, uniqueness might go up but utility goes down
            if k > 35 {
                count / (k - 30) // Penalize k > 35
            } else {
                count
            }
        })
        .map(|(&k, _)| k)
        .unwrap_or(31);
    
    best_k
}

/// Calculate the optimal k value for a set of sequences
pub fn optimal_k(sequences: &[String], max_k: usize) -> usize {
    let hist = kmer_coverage_histogram(sequences, max_k);
    select_best_k(&hist)
}
