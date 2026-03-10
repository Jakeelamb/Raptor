pub struct TranscriptStats {
    pub total: usize,
    pub total_bases: usize,
    pub avg_length: f64,
    pub n50: usize,
    pub n90: usize,
    pub l50: usize,
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

    let threshold = (total_len * numerator).div_ceil(denominator);
    let mut acc = 0usize;
    for (idx, &len) in sorted_desc_lengths.iter().enumerate() {
        acc += len;
        if acc >= threshold {
            return (len, idx + 1);
        }
    }

    (0, sorted_desc_lengths.len())
}

pub fn evaluate_lengths(lengths: &[usize]) -> TranscriptStats {
    if lengths.is_empty() {
        return TranscriptStats {
            total: 0,
            total_bases: 0,
            avg_length: 0.0,
            n50: 0,
            n90: 0,
            l50: 0,
            longest: 0,
        };
    }

    let mut sorted_lengths = lengths.to_vec();
    sorted_lengths.sort_unstable_by(|a, b| b.cmp(a));

    let total_len: usize = sorted_lengths.iter().sum();
    let avg = total_len as f64 / sorted_lengths.len() as f64;
    let (n50, l50) = nx_lx(&sorted_lengths, total_len, 1, 2);
    let (n90, _) = nx_lx(&sorted_lengths, total_len, 9, 10);
    let longest = sorted_lengths.first().copied().unwrap_or(0);

    TranscriptStats {
        total: sorted_lengths.len(),
        total_bases: total_len,
        avg_length: avg,
        n50,
        n90,
        l50,
        longest,
    }
}

#[cfg(test)]
mod tests {
    use super::evaluate_lengths;

    #[test]
    fn evaluate_lengths_reports_nx_metrics() {
        let stats = evaluate_lengths(&[20, 24, 4]);
        assert_eq!(stats.total, 3);
        assert_eq!(stats.total_bases, 48);
        assert_eq!(stats.avg_length, 16.0);
        assert_eq!(stats.n50, 24);
        assert_eq!(stats.n90, 20);
        assert_eq!(stats.l50, 1);
        assert_eq!(stats.longest, 24);
    }
}
