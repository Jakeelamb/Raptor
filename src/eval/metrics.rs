pub struct TranscriptStats {
    pub total: usize,
    pub total_bases: usize,
    pub avg_length: f64,
    pub median_length: f64,
    pub n50: usize,
    pub n75: usize,
    pub n90: usize,
    pub n95: usize,
    pub n99: usize,
    pub l50: usize,
    pub l75: usize,
    pub l90: usize,
    pub l95: usize,
    pub l99: usize,
    pub au_n: f64,
    pub longest: usize,
    pub contigs_ge_1kb: usize,
    pub contigs_ge_10kb: usize,
    pub contigs_ge_50kb: usize,
    pub contigs_ge_100kb: usize,
    pub bases_ge_1kb: usize,
    pub bases_ge_10kb: usize,
    pub bases_ge_50kb: usize,
    pub bases_ge_100kb: usize,
}

#[inline]
fn empty_transcript_stats() -> TranscriptStats {
    TranscriptStats {
        total: 0,
        total_bases: 0,
        avg_length: 0.0,
        median_length: 0.0,
        n50: 0,
        n75: 0,
        n90: 0,
        n95: 0,
        n99: 0,
        l50: 0,
        l75: 0,
        l90: 0,
        l95: 0,
        l99: 0,
        au_n: 0.0,
        longest: 0,
        contigs_ge_1kb: 0,
        contigs_ge_10kb: 0,
        contigs_ge_50kb: 0,
        contigs_ge_100kb: 0,
        bases_ge_1kb: 0,
        bases_ge_10kb: 0,
        bases_ge_50kb: 0,
        bases_ge_100kb: 0,
    }
}

#[derive(Debug, Clone, Copy, Default, Eq, PartialEq)]
pub struct BaseComposition {
    pub gc_bases: usize,
    pub acgt_bases: usize,
    pub n_bases: usize,
    pub ambiguous_bases: usize,
}

impl BaseComposition {
    #[inline]
    pub fn add_sequence(&mut self, seq: &[u8]) {
        for &base in seq {
            match base {
                b'A' | b'a' | b'T' | b't' | b'U' | b'u' => self.acgt_bases += 1,
                b'G' | b'g' | b'C' | b'c' => {
                    self.gc_bases += 1;
                    self.acgt_bases += 1;
                }
                b'N' | b'n' => self.n_bases += 1,
                _ => self.ambiguous_bases += 1,
            }
        }
    }

    #[inline]
    pub fn gc_content(self) -> f64 {
        if self.acgt_bases > 0 {
            self.gc_bases as f64 / self.acgt_bases as f64
        } else {
            0.0
        }
    }

    #[inline]
    pub fn n_content(self, total_bases: usize) -> f64 {
        if total_bases > 0 {
            self.n_bases as f64 / total_bases as f64
        } else {
            0.0
        }
    }

    #[inline]
    pub fn ambiguous_content(self, total_bases: usize) -> f64 {
        if total_bases > 0 {
            self.ambiguous_bases as f64 / total_bases as f64
        } else {
            0.0
        }
    }
}

#[inline]
fn nx_lx(
    sorted_desc_lengths: &[usize],
    total_len: usize,
    numerator: usize,
    denominator: usize,
) -> (usize, usize) {
    if sorted_desc_lengths.is_empty() || total_len == 0 {
        return (0, 0);
    }

    let threshold =
        ((total_len as u128 * numerator as u128) + denominator as u128 - 1) / denominator as u128;
    let mut acc = 0u128;
    for (idx, &len) in sorted_desc_lengths.iter().enumerate() {
        acc += len as u128;
        if acc >= threshold {
            return (len, idx + 1);
        }
    }

    (0, sorted_desc_lengths.len())
}

#[inline]
fn compute_au_n(lengths: &[usize], total_len: usize) -> f64 {
    if lengths.is_empty() || total_len == 0 {
        return 0.0;
    }

    let sum_squares: u128 = lengths
        .iter()
        .map(|&len| {
            let len128 = len as u128;
            len128 * len128
        })
        .sum();
    sum_squares as f64 / total_len as f64
}

#[inline]
fn compute_median_sorted_desc(sorted_lengths: &[usize]) -> f64 {
    if sorted_lengths.is_empty() {
        return 0.0;
    }

    let mid = sorted_lengths.len() / 2;
    if sorted_lengths.len() % 2 == 1 {
        sorted_lengths[mid] as f64
    } else {
        (sorted_lengths[mid - 1] as f64 + sorted_lengths[mid] as f64) / 2.0
    }
}

#[derive(Debug, Clone, Copy, Default, Eq, PartialEq)]
struct LengthBuckets {
    contigs_ge_1kb: usize,
    contigs_ge_10kb: usize,
    contigs_ge_50kb: usize,
    contigs_ge_100kb: usize,
    bases_ge_1kb: usize,
    bases_ge_10kb: usize,
    bases_ge_50kb: usize,
    bases_ge_100kb: usize,
}

#[inline]
fn compute_length_buckets(lengths: &[usize]) -> LengthBuckets {
    let mut buckets = LengthBuckets::default();
    for &len in lengths {
        if len >= 1_000 {
            buckets.contigs_ge_1kb += 1;
            buckets.bases_ge_1kb = buckets.bases_ge_1kb.saturating_add(len);

            if len >= 10_000 {
                buckets.contigs_ge_10kb += 1;
                buckets.bases_ge_10kb = buckets.bases_ge_10kb.saturating_add(len);

                if len >= 50_000 {
                    buckets.contigs_ge_50kb += 1;
                    buckets.bases_ge_50kb = buckets.bases_ge_50kb.saturating_add(len);

                    if len >= 100_000 {
                        buckets.contigs_ge_100kb += 1;
                        buckets.bases_ge_100kb = buckets.bases_ge_100kb.saturating_add(len);
                    }
                }
            }
        }
    }
    buckets
}

/// Evaluate contig/transcript lengths that are already sorted descending.
///
/// This avoids an internal sort and is preferred for hot paths that already
/// maintain descending length order.
pub fn evaluate_lengths_sorted_desc(sorted_lengths: &[usize]) -> TranscriptStats {
    if sorted_lengths.is_empty() {
        return empty_transcript_stats();
    }

    let total_len = sorted_lengths
        .iter()
        .fold(0usize, |acc, &len| acc.saturating_add(len));
    let avg = total_len as f64 / sorted_lengths.len() as f64;
    let median_length = compute_median_sorted_desc(sorted_lengths);
    let (n50, l50) = nx_lx(sorted_lengths, total_len, 1, 2);
    let (n75, l75) = nx_lx(sorted_lengths, total_len, 3, 4);
    let (n90, l90) = nx_lx(sorted_lengths, total_len, 9, 10);
    let (n95, l95) = nx_lx(sorted_lengths, total_len, 19, 20);
    let (n99, l99) = nx_lx(sorted_lengths, total_len, 99, 100);
    let au_n = compute_au_n(sorted_lengths, total_len);
    let longest = sorted_lengths.first().copied().unwrap_or(0);
    let length_buckets = compute_length_buckets(sorted_lengths);

    TranscriptStats {
        total: sorted_lengths.len(),
        total_bases: total_len,
        avg_length: avg,
        median_length,
        n50,
        n75,
        n90,
        n95,
        n99,
        l50,
        l75,
        l90,
        l95,
        l99,
        au_n,
        longest,
        contigs_ge_1kb: length_buckets.contigs_ge_1kb,
        contigs_ge_10kb: length_buckets.contigs_ge_10kb,
        contigs_ge_50kb: length_buckets.contigs_ge_50kb,
        contigs_ge_100kb: length_buckets.contigs_ge_100kb,
        bases_ge_1kb: length_buckets.bases_ge_1kb,
        bases_ge_10kb: length_buckets.bases_ge_10kb,
        bases_ge_50kb: length_buckets.bases_ge_50kb,
        bases_ge_100kb: length_buckets.bases_ge_100kb,
    }
}

pub fn evaluate_lengths(lengths: &[usize]) -> TranscriptStats {
    if lengths.is_empty() {
        return empty_transcript_stats();
    }

    let mut sorted_lengths = lengths.to_vec();
    sorted_lengths.sort_unstable_by(|a, b| b.cmp(a));
    evaluate_lengths_sorted_desc(&sorted_lengths)
}

#[cfg(test)]
mod tests {
    use super::{evaluate_lengths, evaluate_lengths_sorted_desc, nx_lx, BaseComposition};

    #[test]
    fn evaluate_lengths_reports_nx_metrics() {
        let stats = evaluate_lengths(&[20, 24, 4]);
        assert_eq!(stats.total, 3);
        assert_eq!(stats.total_bases, 48);
        assert_eq!(stats.avg_length, 16.0);
        assert_eq!(stats.median_length, 20.0);
        assert_eq!(stats.n50, 24);
        assert_eq!(stats.n75, 20);
        assert_eq!(stats.n90, 20);
        assert_eq!(stats.n95, 4);
        assert_eq!(stats.n99, 4);
        assert_eq!(stats.l50, 1);
        assert_eq!(stats.l75, 2);
        assert_eq!(stats.l90, 2);
        assert_eq!(stats.l95, 3);
        assert_eq!(stats.l99, 3);
        assert!((stats.au_n - 20.6666666667).abs() < 1e-6);
        assert_eq!(stats.longest, 24);
        assert_eq!(stats.contigs_ge_1kb, 0);
        assert_eq!(stats.contigs_ge_10kb, 0);
        assert_eq!(stats.contigs_ge_50kb, 0);
        assert_eq!(stats.contigs_ge_100kb, 0);
        assert_eq!(stats.bases_ge_1kb, 0);
        assert_eq!(stats.bases_ge_10kb, 0);
        assert_eq!(stats.bases_ge_50kb, 0);
        assert_eq!(stats.bases_ge_100kb, 0);
    }

    #[test]
    fn nx_lx_handles_large_totals_without_overflow() {
        let very_large = usize::MAX - 3;
        let lengths = vec![very_large];
        let (n90, l90) = nx_lx(&lengths, very_large, 9, 10);
        assert_eq!(n90, very_large);
        assert_eq!(l90, 1);
    }

    #[test]
    fn evaluate_lengths_reports_length_bucket_counts_and_spans() {
        let stats = evaluate_lengths(&[50_000, 10_000, 1_000, 999]);
        assert_eq!(stats.contigs_ge_1kb, 3);
        assert_eq!(stats.contigs_ge_10kb, 2);
        assert_eq!(stats.contigs_ge_50kb, 1);
        assert_eq!(stats.contigs_ge_100kb, 0);
        assert_eq!(stats.bases_ge_1kb, 61_000);
        assert_eq!(stats.bases_ge_10kb, 60_000);
        assert_eq!(stats.bases_ge_50kb, 50_000);
        assert_eq!(stats.bases_ge_100kb, 0);
    }

    #[test]
    fn evaluate_lengths_reports_100kb_bucket_counts_and_spans() {
        let stats = evaluate_lengths(&[100_000, 99_999, 1_000]);
        assert_eq!(stats.contigs_ge_100kb, 1);
        assert_eq!(stats.bases_ge_100kb, 100_000);
    }

    #[test]
    fn evaluate_lengths_reports_even_count_median_as_midpoint() {
        let stats = evaluate_lengths(&[10, 8, 6, 4]);
        assert_eq!(stats.median_length, 7.0);
    }

    #[test]
    fn base_composition_classifies_bases_consistently() {
        let mut composition = BaseComposition::default();
        composition.add_sequence(b"GgCcAaTtUuNnRY");

        assert_eq!(composition.gc_bases, 4);
        assert_eq!(composition.acgt_bases, 10);
        assert_eq!(composition.n_bases, 2);
        assert_eq!(composition.ambiguous_bases, 2);
        assert!((composition.gc_content() - 0.4).abs() < 1e-12);
        assert!((composition.n_content(14) - (2.0 / 14.0)).abs() < 1e-12);
        assert!((composition.ambiguous_content(14) - (2.0 / 14.0)).abs() < 1e-12);
    }

    #[test]
    fn evaluate_lengths_sorted_desc_matches_unsorted_entrypoint() {
        let input = vec![999, 50_000, 1_000, 10_000];
        let expected = evaluate_lengths(&input);

        let mut sorted = input.clone();
        sorted.sort_unstable_by(|a, b| b.cmp(a));
        let observed = evaluate_lengths_sorted_desc(&sorted);

        assert_eq!(observed.total, expected.total);
        assert_eq!(observed.total_bases, expected.total_bases);
        assert!((observed.median_length - expected.median_length).abs() < 1e-12);
        assert_eq!(observed.n50, expected.n50);
        assert_eq!(observed.n75, expected.n75);
        assert_eq!(observed.n90, expected.n90);
        assert_eq!(observed.n95, expected.n95);
        assert_eq!(observed.n99, expected.n99);
        assert_eq!(observed.l50, expected.l50);
        assert_eq!(observed.l75, expected.l75);
        assert_eq!(observed.l90, expected.l90);
        assert_eq!(observed.l95, expected.l95);
        assert_eq!(observed.l99, expected.l99);
        assert!((observed.au_n - expected.au_n).abs() < 1e-12);
        assert_eq!(observed.longest, expected.longest);
        assert_eq!(observed.contigs_ge_1kb, expected.contigs_ge_1kb);
        assert_eq!(observed.contigs_ge_10kb, expected.contigs_ge_10kb);
        assert_eq!(observed.contigs_ge_50kb, expected.contigs_ge_50kb);
        assert_eq!(observed.contigs_ge_100kb, expected.contigs_ge_100kb);
        assert_eq!(observed.bases_ge_1kb, expected.bases_ge_1kb);
        assert_eq!(observed.bases_ge_10kb, expected.bases_ge_10kb);
        assert_eq!(observed.bases_ge_50kb, expected.bases_ge_50kb);
        assert_eq!(observed.bases_ge_100kb, expected.bases_ge_100kb);
    }

    #[test]
    fn evaluate_lengths_handles_all_zero_lengths() {
        let stats = evaluate_lengths(&[0, 0, 0, 0]);
        assert_eq!(stats.total, 4);
        assert_eq!(stats.total_bases, 0);
        assert_eq!(stats.avg_length, 0.0);
        assert_eq!(stats.median_length, 0.0);
        assert_eq!(stats.n50, 0);
        assert_eq!(stats.n75, 0);
        assert_eq!(stats.n90, 0);
        assert_eq!(stats.n95, 0);
        assert_eq!(stats.n99, 0);
        assert_eq!(stats.l50, 0);
        assert_eq!(stats.l75, 0);
        assert_eq!(stats.l90, 0);
        assert_eq!(stats.l95, 0);
        assert_eq!(stats.l99, 0);
        assert_eq!(stats.au_n, 0.0);
        assert_eq!(stats.longest, 0);
        assert_eq!(stats.contigs_ge_1kb, 0);
        assert_eq!(stats.contigs_ge_10kb, 0);
        assert_eq!(stats.contigs_ge_50kb, 0);
        assert_eq!(stats.contigs_ge_100kb, 0);
        assert_eq!(stats.bases_ge_1kb, 0);
        assert_eq!(stats.bases_ge_10kb, 0);
        assert_eq!(stats.bases_ge_50kb, 0);
        assert_eq!(stats.bases_ge_100kb, 0);
    }
}
