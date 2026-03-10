use crate::graph::transcript::Transcript;
use crate::io::fastq::FastqRecord;
use std::io::{BufWriter, Write};

#[inline]
fn ordered_transcript_rows<'a>(
    transcripts: &'a [Transcript],
    tpms: &[f64],
) -> Vec<(&'a Transcript, f64)> {
    let mut rows: Vec<(usize, &Transcript, f64)> = transcripts
        .iter()
        .enumerate()
        .map(|(idx, tx)| (idx, tx, tpms.get(idx).copied().unwrap_or(0.0)))
        .collect();
    rows.sort_unstable_by(|a, b| a.1.id.cmp(&b.1.id).then_with(|| a.0.cmp(&b.0)));
    rows.into_iter().map(|(_, tx, tpm)| (tx, tpm)).collect()
}

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
    let mut norm_counts = Vec::with_capacity(transcripts.len());
    for (idx, transcript) in transcripts.iter().enumerate() {
        let count = counts.get(idx).copied().unwrap_or(0);
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
        return vec![0.0; transcripts.len()];
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
pub fn write_tpm_table(
    transcripts: &[Transcript],
    tpms: &[f64],
    output: &str,
) -> std::io::Result<()> {
    let file = std::fs::File::create(output)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "transcript_id\tlength\ttpm")?;
    for (tx, tpm) in ordered_transcript_rows(transcripts, tpms) {
        writeln!(
            writer,
            "transcript_{}\t{}\t{:.2}",
            tx.id,
            tx.sequence.len(),
            tpm
        )?;
    }
    writer.flush()
}

/// Filter transcripts based on a minimum TPM threshold
pub fn filter_by_tpm(
    transcripts: &[Transcript],
    tpms: &[f64],
    min_tpm: f64,
) -> (Vec<Transcript>, Vec<f64>) {
    let mut kept_tx = vec![];
    let mut kept_tpms = vec![];

    for (idx, t) in transcripts.iter().enumerate() {
        let tpm = tpms.get(idx).copied().unwrap_or(0.0);
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
    use tempfile::NamedTempFile;

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

    #[test]
    fn compute_tpm_pads_missing_counts_to_transcript_count() {
        let transcripts = vec![
            make_transcript(5, "AAAAAA"),
            make_transcript(2, "CCCCCC"),
            make_transcript(9, "GGGGGG"),
        ];
        let tpms = compute_tpm(&[10], &transcripts);

        assert_eq!(tpms.len(), transcripts.len());
        assert!(tpms[0] > 0.0);
        assert_eq!(tpms[1], 0.0);
        assert_eq!(tpms[2], 0.0);
        assert!((tpms.iter().sum::<f64>() - 1_000_000.0).abs() < 1e-6);
    }

    #[test]
    fn filter_by_tpm_treats_missing_values_as_zero() {
        let transcripts = vec![
            make_transcript(1, "AAAA"),
            make_transcript(2, "CCCC"),
            make_transcript(3, "GGGG"),
        ];

        let (filtered_tx, filtered_tpms) = filter_by_tpm(&transcripts, &[5.0], 1.0);
        let filtered_ids: Vec<usize> = filtered_tx.iter().map(|tx| tx.id).collect();

        assert_eq!(filtered_ids, vec![1]);
        assert_eq!(filtered_tpms, vec![5.0]);
    }

    #[test]
    fn write_tpm_table_is_sorted_by_transcript_id_and_pads_missing_values() {
        let transcripts = vec![
            make_transcript(9, "AAA"),
            make_transcript(2, "CCCC"),
            make_transcript(5, "GG"),
        ];
        let output = NamedTempFile::new().unwrap();

        write_tpm_table(&transcripts, &[7.0], output.path().to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(output.path()).unwrap();
        let lines: Vec<&str> = content.lines().collect();

        assert_eq!(lines[0], "transcript_id\tlength\ttpm");
        assert_eq!(lines[1], "transcript_2\t4\t0.00");
        assert_eq!(lines[2], "transcript_5\t2\t0.00");
        assert_eq!(lines[3], "transcript_9\t3\t7.00");
    }
}
