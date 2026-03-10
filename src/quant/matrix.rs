use crate::graph::transcript::Transcript;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/// Generate matrix from isoform -> TPM mapping per sample
pub fn write_counts_matrix(
    samples: &HashMap<String, Vec<f64>>,
    transcripts: &[Transcript],
    output: &str,
) -> std::io::Result<()> {
    let file = std::fs::File::create(output)?;
    let mut writer = BufWriter::new(file);
    let mut sample_names: Vec<&str> = samples.keys().map(|name| name.as_str()).collect();
    sample_names.sort_unstable();

    // Header
    write!(writer, "transcript_id")?;
    for sample in &sample_names {
        write!(writer, "\t{}", sample)?;
    }
    writeln!(writer)?;

    for (tx_idx, tx) in transcripts.iter().enumerate() {
        write!(writer, "transcript_{}", tx.id)?;
        for sample in &sample_names {
            let tpm = samples
                .get(*sample)
                .and_then(|tpms| tpms.get(tx_idx))
                .copied()
                .unwrap_or(0.0);
            write!(writer, "\t{:.2}", tpm)?;
        }
        writeln!(writer)?;
    }

    writer.flush()
}

/// Read a TPM matrix file into a map of sample name -> transcript TPM values
pub fn read_tpm_matrix(path: &str) -> std::io::Result<HashMap<String, Vec<f64>>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Read header to get sample names
    let header = match lines.next() {
        Some(Ok(line)) => line,
        _ => {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Empty or invalid matrix file",
            ))
        }
    };

    let headers: Vec<&str> = header.split('\t').collect();
    if headers.len() < 2 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Invalid matrix format - no sample columns found",
        ));
    }

    let sample_headers = &headers[1..];
    let mut seen_samples = std::collections::HashSet::with_capacity(sample_headers.len());
    for &sample in sample_headers {
        if !seen_samples.insert(sample) {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Duplicate sample column: {}", sample),
            ));
        }
    }

    // Initialize result map
    let mut result: HashMap<String, Vec<f64>> = HashMap::new();
    for &h in sample_headers {
        // Skip transcript_id column
        result.insert(h.to_string(), Vec::new());
    }

    // Process each data row
    for line in lines {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let mut fields = line.split('\t');
        let transcript_id = fields.next().unwrap_or("").trim();
        if transcript_id.is_empty() {
            continue;
        }

        // Parse each sample's TPM value.
        // Missing/invalid values are padded with 0.0 to keep row alignment deterministic.
        for &sample in sample_headers {
            let parsed = fields
                .next()
                .and_then(|raw| raw.parse::<f64>().ok())
                .filter(|value| value.is_finite())
                .unwrap_or(0.0);
            result
                .get_mut(sample)
                .expect("sample header exists")
                .push(parsed);
        }
    }

    Ok(result)
}

/// Writes transcript information to a counts matrix format
///
/// # Arguments
/// * `transcripts` - Vector of transcripts to write
/// * `out_path` - Output file path
///
/// # Returns
/// * IO Result
pub fn write_isoform_counts_matrix(
    transcripts: &[Transcript],
    out_path: &str,
) -> std::io::Result<()> {
    let file = File::create(out_path)?;
    let mut w = BufWriter::new(file);
    let mut ordered: Vec<&Transcript> = transcripts.iter().collect();
    ordered.sort_unstable_by_key(|tx| tx.id);

    writeln!(w, "transcript_id\tlength\ttpm\tconfidence")?;

    for tx in ordered {
        writeln!(
            w,
            "transcript_{}\t{}\t{:.3}\t{:.2}",
            tx.id,
            tx.sequence.len(),
            tx.tpm.unwrap_or(0.0),
            tx.confidence
        )?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use tempfile::NamedTempFile;

    fn make_transcript(id: usize) -> Transcript {
        Transcript {
            id,
            sequence: "ACGT".to_string(),
            path: vec![id],
            confidence: 0.9,
            length: 4,
            strand: '+',
            tpm: None,
            splicing: "linear".to_string(),
        }
    }

    #[test]
    fn write_counts_matrix_sorts_headers_and_is_insertion_order_independent() {
        let transcripts = vec![make_transcript(1), make_transcript(2)];

        let mut samples_a = HashMap::new();
        samples_a.insert("zeta".to_string(), vec![2.0, 4.0]);
        samples_a.insert("alpha".to_string(), vec![1.0, 3.0]);

        let mut samples_b = HashMap::new();
        samples_b.insert("alpha".to_string(), vec![1.0, 3.0]);
        samples_b.insert("zeta".to_string(), vec![2.0, 4.0]);

        let out_a = NamedTempFile::new().unwrap();
        let out_b = NamedTempFile::new().unwrap();

        write_counts_matrix(&samples_a, &transcripts, out_a.path().to_str().unwrap()).unwrap();
        write_counts_matrix(&samples_b, &transcripts, out_b.path().to_str().unwrap()).unwrap();

        let content_a = std::fs::read_to_string(out_a.path()).unwrap();
        let content_b = std::fs::read_to_string(out_b.path()).unwrap();
        assert_eq!(content_a, content_b);

        let mut lines = content_a.lines();
        assert_eq!(lines.next().unwrap(), "transcript_id\talpha\tzeta");
        assert_eq!(lines.next().unwrap(), "transcript_1\t1.00\t2.00");
        assert_eq!(lines.next().unwrap(), "transcript_2\t3.00\t4.00");
    }

    #[test]
    fn write_counts_matrix_pads_missing_values_with_zero() {
        let transcripts = vec![make_transcript(1), make_transcript(2), make_transcript(3)];

        let mut samples = HashMap::new();
        samples.insert("sample".to_string(), vec![5.0]);

        let out = NamedTempFile::new().unwrap();
        write_counts_matrix(&samples, &transcripts, out.path().to_str().unwrap()).unwrap();

        let content = std::fs::read_to_string(out.path()).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines[0], "transcript_id\tsample");
        assert_eq!(lines[1], "transcript_1\t5.00");
        assert_eq!(lines[2], "transcript_2\t0.00");
        assert_eq!(lines[3], "transcript_3\t0.00");
    }

    #[test]
    fn read_tpm_matrix_pads_missing_and_invalid_values_per_row() {
        let file = NamedTempFile::new().unwrap();
        std::fs::write(
            file.path(),
            "transcript_id\ta\tb\n\
             transcript_1\t1.25\t2.50\n\
             transcript_2\t3.00\n\
             transcript_3\tnan\t4.75\n\
             transcript_4\tbad\t\n",
        )
        .unwrap();

        let parsed = read_tpm_matrix(file.path().to_str().unwrap()).unwrap();
        assert_eq!(parsed.get("a"), Some(&vec![1.25, 3.0, 0.0, 0.0]));
        assert_eq!(parsed.get("b"), Some(&vec![2.5, 0.0, 4.75, 0.0]));
    }

    #[test]
    fn read_tpm_matrix_rejects_duplicate_sample_columns() {
        let file = NamedTempFile::new().unwrap();
        std::fs::write(
            file.path(),
            "transcript_id\ts1\ts1\ntranscript_1\t1.0\t2.0\n",
        )
        .unwrap();

        let err = read_tpm_matrix(file.path().to_str().unwrap()).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("Duplicate sample column"));
    }

    #[test]
    fn read_tpm_matrix_ignores_blank_and_comment_rows() {
        let file = NamedTempFile::new().unwrap();
        std::fs::write(
            file.path(),
            "transcript_id\ta\tb\n\
             \n\
             # this is a comment\n\
             transcript_1\t1.0\t2.0\n\
             \t\n\
             transcript_2\t3.5\t4.5\n",
        )
        .unwrap();

        let parsed = read_tpm_matrix(file.path().to_str().unwrap()).unwrap();
        assert_eq!(parsed.get("a"), Some(&vec![1.0, 3.5]));
        assert_eq!(parsed.get("b"), Some(&vec![2.0, 4.5]));
    }

    fn round_two(v: f64) -> f64 {
        (v * 100.0).round() / 100.0
    }

    #[test]
    fn write_isoform_counts_matrix_is_sorted_by_transcript_id() {
        let out = NamedTempFile::new().unwrap();
        let transcripts = vec![make_transcript(9), make_transcript(2), make_transcript(5)];

        write_isoform_counts_matrix(&transcripts, out.path().to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(out.path()).unwrap();
        let lines: Vec<&str> = content.lines().collect();

        assert_eq!(lines[0], "transcript_id\tlength\ttpm\tconfidence");
        assert!(lines[1].starts_with("transcript_2\t"));
        assert!(lines[2].starts_with("transcript_5\t"));
        assert!(lines[3].starts_with("transcript_9\t"));
    }

    proptest! {
        #[test]
        fn write_then_read_counts_matrix_round_trips_to_two_decimal_precision(
            tx_count in 1usize..20,
            sample_data in prop::collection::btree_map("[a-z]{1,6}", prop::collection::vec(-1_000.0f64..1_000.0, 1..20), 1..6)
        ) {
            let transcripts: Vec<Transcript> = (0..tx_count).map(make_transcript).collect();
            let mut samples = HashMap::new();

            for (name, mut values) in sample_data {
                values.truncate(tx_count);
                while values.len() < tx_count {
                    values.push(0.0);
                }
                samples.insert(name, values);
            }

            let out = NamedTempFile::new().unwrap();
            write_counts_matrix(&samples, &transcripts, out.path().to_str().unwrap()).unwrap();
            let parsed = read_tpm_matrix(out.path().to_str().unwrap()).unwrap();

            for (name, expected_values) in samples {
                let observed = parsed.get(&name).expect("sample present after round-trip");
                prop_assert_eq!(observed.len(), tx_count);
                for (observed_v, expected_v) in observed.iter().zip(expected_values.iter()) {
                    prop_assert!((observed_v - round_two(*expected_v)).abs() < 1e-9);
                }
            }
        }
    }
}
