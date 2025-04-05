use crate::graph::transcript::Transcript;

pub struct TranscriptStats {
    pub total: usize,
    pub avg_length: f64,
    pub n50: usize,
}

pub fn evaluate(transcripts: &[Transcript]) -> TranscriptStats {
    let mut lengths: Vec<usize> = transcripts.iter().map(|t| t.sequence.len()).collect();
    lengths.sort_unstable();
    let total_len: usize = lengths.iter().sum();
    let avg = total_len as f64 / lengths.len() as f64;

    let mut acc = 0;
    let n50 = lengths.iter().rev().find(|&&l| {
        acc += l;
        acc >= total_len / 2
    }).copied().unwrap_or(0);

    TranscriptStats {
        total: lengths.len(),
        avg_length: avg,
        n50,
    }
} 