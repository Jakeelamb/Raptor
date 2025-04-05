use crate::graph::transcript::Transcript;
use crate::io::fastq::FastqRecord;

/// Naïve count of exact matches to transcript sequence
pub fn count_reads(transcripts: &[Transcript], reads: &[FastqRecord]) -> Vec<usize> {
    transcripts.iter().map(|t| {
        reads.iter().filter(|r| t.sequence.contains(&r.sequence)).count()
    }).collect()
}

/// Compute TPM values from raw counts and transcript lengths
pub fn compute_tpm(counts: &[usize], transcripts: &[Transcript]) -> Vec<f64> {
    let mut norm_counts = vec![];
    for (i, &count) in counts.iter().enumerate() {
        let len_kb = transcripts[i].sequence.len() as f64 / 1000.0;
        norm_counts.push(count as f64 / len_kb);
    }

    let sum: f64 = norm_counts.iter().sum();
    norm_counts.into_iter().map(|x| (x / sum) * 1_000_000.0).collect()
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