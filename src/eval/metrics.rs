pub struct TranscriptStats {
    pub total: usize,
    pub total_bases: usize,
    pub avg_length: f64,
    pub n50: usize,
    pub n75: usize,
    pub n90: usize,
    pub l50: usize,
    pub l90: usize,
    pub au_n: f64,
    pub longest: usize,
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

pub fn evaluate_lengths(lengths: &[usize]) -> TranscriptStats {
    if lengths.is_empty() {
        return TranscriptStats {
            total: 0,
            total_bases: 0,
            avg_length: 0.0,
            n50: 0,
            n75: 0,
            n90: 0,
            l50: 0,
            l90: 0,
            au_n: 0.0,
            longest: 0,
        };
    }

    let mut sorted_lengths = lengths.to_vec();
    sorted_lengths.sort_unstable_by(|a, b| b.cmp(a));

    let total_len = sorted_lengths
        .iter()
        .fold(0usize, |acc, &len| acc.saturating_add(len));
    let avg = total_len as f64 / sorted_lengths.len() as f64;
    let (n50, l50) = nx_lx(&sorted_lengths, total_len, 1, 2);
    let (n75, _) = nx_lx(&sorted_lengths, total_len, 3, 4);
    let (n90, l90) = nx_lx(&sorted_lengths, total_len, 9, 10);
    let au_n = compute_au_n(&sorted_lengths, total_len);
    let longest = sorted_lengths.first().copied().unwrap_or(0);

    TranscriptStats {
        total: sorted_lengths.len(),
        total_bases: total_len,
        avg_length: avg,
        n50,
        n75,
        n90,
        l50,
        l90,
        au_n,
        longest,
    }
}

#[cfg(test)]
mod tests {
    use super::{evaluate_lengths, nx_lx};

    #[test]
    fn evaluate_lengths_reports_nx_metrics() {
        let stats = evaluate_lengths(&[20, 24, 4]);
        assert_eq!(stats.total, 3);
        assert_eq!(stats.total_bases, 48);
        assert_eq!(stats.avg_length, 16.0);
        assert_eq!(stats.n50, 24);
        assert_eq!(stats.n75, 20);
        assert_eq!(stats.n90, 20);
        assert_eq!(stats.l50, 1);
        assert_eq!(stats.l90, 2);
        assert!((stats.au_n - 20.6666666667).abs() < 1e-6);
        assert_eq!(stats.longest, 24);
    }

    #[test]
    fn nx_lx_handles_large_totals_without_overflow() {
        let very_large = usize::MAX - 3;
        let lengths = vec![very_large];
        let (n90, l90) = nx_lx(&lengths, very_large, 9, 10);
        assert_eq!(n90, very_large);
        assert_eq!(l90, 1);
    }
}
