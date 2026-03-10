use crate::graph::transcript::Transcript;
use crate::io::fastq::FastqRecord;

/// Naïve count of exact matches to transcript sequence
pub fn count_reads(transcripts: &[Transcript], reads: &[FastqRecord]) -> Vec<usize> {
    transcripts
        .iter()
        .map(|t| {
            reads
                .iter()
                .filter(|r| t.sequence.contains(&r.sequence))
                .count()
        })
        .collect()
}

/// Compute TPM values from raw counts and transcript lengths
pub fn compute_tpm(counts: &[usize], transcripts: &[Transcript]) -> Vec<f64> {
    let mut norm_counts = Vec::with_capacity(counts.len().min(transcripts.len()));
    for (&count, transcript) in counts.iter().zip(transcripts.iter()) {
        let len_bp = transcript.sequence.len();
        if len_bp == 0 {
            norm_counts.push(0.0);
            continue;
        }

        let len_kb = len_bp as f64 / 1000.0;
        norm_counts.push(count as f64 / len_kb);
    }

    let sum: f64 = norm_counts.iter().sum();
    if sum <= f64::EPSILON || !sum.is_finite() {
        return vec![0.0; norm_counts.len()];
    }

    norm_counts
        .into_iter()
        .map(|x| {
            if x.is_finite() {
                (x / sum) * 1_000_000.0
            } else {
                0.0
            }
        })
        .collect()
}

/// Write TPM output table
pub fn write_tpm_table(transcripts: &[Transcript], tpms: &[f64], output: &str) {
    let mut f = std::fs::File::create(output).expect("Failed to write TPM table");
    use std::io::Write;
    writeln!(f, "transcript_id\tlength\ttpm").unwrap();
    for (t, &v) in transcripts.iter().zip(tpms.iter()) {
        writeln!(f, "transcript_{}\t{}\t{:.2}", t.id, t.sequence.len(), v).unwrap();
    }
}

/// Filter transcripts based on a minimum TPM threshold
pub fn filter_by_tpm(
    transcripts: &[Transcript],
    tpms: &[f64],
    min_tpm: f64,
) -> (Vec<Transcript>, Vec<f64>) {
    let mut kept_tx = vec![];
    let mut kept_tpms = vec![];

    for (t, &tpm) in transcripts.iter().zip(tpms.iter()) {
        if tpm >= min_tpm {
            kept_tx.push(t.clone());
            kept_tpms.push(tpm);
        }
    }

    (kept_tx, kept_tpms)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_transcript(id: usize, sequence: &str) -> Transcript {
        Transcript {
            id,
            sequence: sequence.to_string(),
            path: vec![id],
            confidence: 1.0,
            length: sequence.len(),
            strand: '+',
            tpm: None,
            splicing: "linear".to_string(),
        }
    }

    #[test]
    fn compute_tpm_handles_zero_length_transcripts_without_nan_or_inf() {
        let transcripts = vec![make_transcript(1, ""), make_transcript(2, "ACGT")];
        let tpms = compute_tpm(&[10, 5], &transcripts);

        assert_eq!(tpms.len(), 2);
        assert_eq!(tpms[0], 0.0);
        assert!(tpms[1].is_finite());
        assert!((tpms.iter().sum::<f64>() - 1_000_000.0).abs() < 1e-6);
    }

    #[test]
    fn compute_tpm_returns_zeroes_for_all_zero_counts() {
        let transcripts = vec![make_transcript(1, "AAAA"), make_transcript(2, "CCCC")];
        let tpms = compute_tpm(&[0, 0], &transcripts);
        assert_eq!(tpms, vec![0.0, 0.0]);
    }
}
