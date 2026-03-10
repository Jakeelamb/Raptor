use crate::kmer::nthash::NtHashIterator;
use ahash::AHashSet;
use std::collections::HashMap;

const NTHASH_MAX_K: usize = 32;

#[inline]
fn kmer_score(k: usize, count: usize) -> usize {
    if k > 35 {
        count / (k - 30)
    } else {
        count
    }
}

#[inline]
fn candidate_k_values(max_k: usize) -> Vec<usize> {
    let capped_max = max_k.min(NTHASH_MAX_K);
    if capped_max == 0 {
        return Vec::new();
    }

    let template: &[usize] = if capped_max <= 21 {
        &[15, 17, 19, 21]
    } else if capped_max <= 31 {
        &[21, 25, 27, 31]
    } else {
        &[21, 25, 31, 32]
    };

    let mut k_values: Vec<usize> = template
        .iter()
        .copied()
        .filter(|&k| k <= capped_max)
        .collect();

    // Always include the requested upper bound (after ntHash capping).
    if k_values.last().copied() != Some(capped_max) {
        k_values.push(capped_max);
    }

    k_values.sort_unstable();
    k_values.dedup();
    k_values
}

#[inline]
fn kmer_coverage_histogram_for_candidates(
    sequences: &[String],
    k_values: &[usize],
) -> HashMap<usize, usize> {
    let mut counts: HashMap<usize, usize> = HashMap::new();

    for &k in k_values {
        // Use u64 hash set for memory efficiency.
        let mut seen: AHashSet<u64> = AHashSet::new();

        for seq in sequences {
            let bytes = seq.as_bytes();
            // Use ntHash for O(1) rolling hash per k-mer.
            for (_, hash) in NtHashIterator::new(bytes, k) {
                seen.insert(hash);
            }
        }

        counts.insert(k, seen.len());
    }

    counts
}

/// Generate a histogram of unique k-mer counts for different k values.
/// Uses ntHash for fast O(1) rolling hash computation.
pub fn kmer_coverage_histogram(sequences: &[String], max_k: usize) -> HashMap<usize, usize> {
    let k_values = candidate_k_values(max_k);
    kmer_coverage_histogram_for_candidates(sequences, &k_values)
}

#[inline]
fn select_best_k_with_default(hist: &HashMap<usize, usize>, default_k: usize) -> usize {
    hist.iter()
        .max_by(|(k_a, count_a), (k_b, count_b)| {
            let score_a = kmer_score(**k_a, **count_a);
            let score_b = kmer_score(**k_b, **count_b);

            // Deterministic tie-breaker: prefer smaller k in score ties.
            score_a.cmp(&score_b).then_with(|| k_b.cmp(k_a))
        })
        .map(|(&k, _)| k)
        .unwrap_or(default_k)
}

/// Selects the best k value based on the k-mer histogram
/// This uses a heuristic that balances uniqueness with k-mer length
pub fn select_best_k(hist: &HashMap<usize, usize>) -> usize {
    select_best_k_with_default(hist, 31)
}

/// Calculate the optimal k value for a set of sequences
pub fn optimal_k(sequences: &[String], max_k: usize) -> usize {
    let hist = kmer_coverage_histogram(sequences, max_k);
    let fallback_k = candidate_k_values(max_k).into_iter().next().unwrap_or(1);
    select_best_k_with_default(&hist, fallback_k)
}

#[cfg(test)]
mod tests {
    use super::{candidate_k_values, kmer_coverage_histogram, optimal_k, select_best_k};
    use std::collections::HashMap;

    #[test]
    fn candidate_values_are_bounded_and_include_requested_upper_bound() {
        assert_eq!(candidate_k_values(0), Vec::<usize>::new());
        assert_eq!(candidate_k_values(10), vec![10]);
        assert_eq!(candidate_k_values(22), vec![21, 22]);
        assert_eq!(candidate_k_values(32), vec![21, 25, 31, 32]);
        assert_eq!(candidate_k_values(41), vec![21, 25, 31, 32]);
    }

    #[test]
    fn select_best_k_prefers_smaller_k_when_scores_tie() {
        let mut hist = HashMap::new();
        hist.insert(31usize, 42usize);
        hist.insert(25usize, 42usize);
        hist.insert(21usize, 42usize);

        assert_eq!(select_best_k(&hist), 21);
    }

    #[test]
    fn optimal_k_respects_upper_bound_for_short_inputs() {
        let sequences = vec!["ACGTACGT".to_string(), "TGCATGCA".to_string()];

        assert_eq!(optimal_k(&sequences, 10), 10);
    }

    #[test]
    fn histogram_uses_only_admissible_k_values() {
        let sequences = vec!["A".repeat(64), "C".repeat(64)];

        let hist = kmer_coverage_histogram(&sequences, 41);
        let mut keys: Vec<usize> = hist.keys().copied().collect();
        keys.sort_unstable();

        assert_eq!(keys, vec![21, 25, 31, 32]);
    }
}
