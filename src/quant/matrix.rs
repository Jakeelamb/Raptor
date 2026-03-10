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

    // Initialize result map
    let mut result: HashMap<String, Vec<f64>> = HashMap::new();
    for &h in headers.iter().skip(1) {
        // Skip transcript_id column
        result.insert(h.to_string(), Vec::new());
    }

    // Process each data row
    for line in lines {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < headers.len() {
            continue; // Skip invalid lines
        }

        // Parse each sample's TPM value
        for (idx, &sample) in headers.iter().skip(1).enumerate() {
            let col_idx = idx + 1; // +1 because we skipped the first column
            if col_idx < fields.len() {
                if let Ok(tpm) = fields[col_idx].parse::<f64>() {
                    result.get_mut(sample).unwrap().push(tpm);
                } else {
                    // Add 0.0 for invalid/empty values
                    result.get_mut(sample).unwrap().push(0.0);
                }
            }
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

    writeln!(w, "transcript_id\tlength\ttpm\tconfidence")?;

    for tx in transcripts {
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
}
