pub struct TranscriptStats {
    pub total: usize,
    pub avg_length: f64,
    pub n50: usize,
}

pub fn evaluate_lengths(lengths: &[usize]) -> TranscriptStats {
    if lengths.is_empty() {
        return TranscriptStats {
            total: 0,
            avg_length: 0.0,
            n50: 0,
        };
    }
    
    let mut sorted_lengths = lengths.to_vec();
    sorted_lengths.sort_unstable();
    
    let total_len: usize = sorted_lengths.iter().sum();
    let avg = total_len as f64 / sorted_lengths.len() as f64;

    let mut acc = 0;
    let n50 = sorted_lengths.iter().rev().find(|&&l| {
        acc += l;
        acc >= total_len / 2
    }).copied().unwrap_or(0);

    TranscriptStats {
        total: sorted_lengths.len(),
        avg_length: avg,
        n50,
    }
} 