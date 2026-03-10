use crate::graph::transcript::Transcript;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Error, ErrorKind, Result};

/// Estimates transcript expression in TPM from a BAM alignment file
///
/// # Arguments
/// * `bam_path` - Path to the BAM file containing alignments
/// * `transcripts` - Vector of transcript objects
///
/// # Returns
/// * Vector of TPM values corresponding to each transcript
pub fn estimate_tpm_from_bam(_bam_path: &str, transcripts: &[Transcript]) -> Result<Vec<f64>> {
    // In a real implementation, we would use rust-htslib to parse the BAM file
    // For this example, we'll simulate the BAM reading process

    // Mock counts for demonstration purposes
    let mut counts = vec![0usize; transcripts.len()];

    // In an actual implementation, this would read the BAM file:
    /*
    let mut reader = Reader::from_path(bam_path)?;
    for record in reader.records().filter_map(Result::ok) {
        if let Some(tid) = record.tid().checked_into::<usize>().ok() {
            if tid < transcripts.len() {
                counts[tid] += 1;
            }
        }
    }
    */

    // For demonstration, populate with random counts
    for i in 0..counts.len() {
        // Random count proportional to transcript length (longer transcripts get more reads)
        counts[i] = transcripts[i].length / 100 + i % 10;
    }

    // Calculate TPM values
    let lengths: Vec<f64> = transcripts.iter().map(|t| t.length as f64).collect();
    let rpk: Vec<f64> = counts
        .iter()
        .zip(lengths.iter())
        .map(|(&c, &l)| if l > 0.0 { c as f64 * 1000.0 / l } else { 0.0 })
        .collect();

    let scaling_factor: f64 = rpk.iter().sum::<f64>() / 1_000_000.0;

    let tpm = if scaling_factor > 0.0 {
        rpk.iter().map(|&r| r / scaling_factor).collect()
    } else {
        vec![0.0; transcripts.len()]
    };

    Ok(tpm)
}

/// Updates transcripts with TPM values from alignment file
///
/// # Arguments
/// * `transcripts` - Mutable vector of transcripts to update
/// * `bam_path` - Path to BAM file
///
/// # Returns
/// * Result with number of transcripts updated
pub fn update_transcripts_with_tpm(
    transcripts: &mut [Transcript],
    bam_path: &str,
) -> Result<usize> {
    let tpm_values = estimate_tpm_from_bam(bam_path, transcripts)?;

    for (transcript, tpm) in transcripts.iter_mut().zip(tpm_values.iter()) {
        transcript.tpm = Some(*tpm);
    }

    Ok(transcripts.len())
}

/// Filters transcripts by TPM threshold
///
/// # Arguments
/// * `transcripts` - Vector of transcripts to filter
/// * `min_tpm` - Minimum TPM value to keep a transcript
///
/// # Returns
/// * Filtered vector of transcripts
pub fn filter_by_tpm(transcripts: Vec<Transcript>, min_tpm: f64) -> Vec<Transcript> {
    transcripts
        .into_iter()
        .filter(|t| t.tpm.unwrap_or(0.0) >= min_tpm)
        .collect()
}

/// Creates an expression matrix for multiple samples
///
/// # Arguments
/// * `transcripts` - Vector of transcripts
/// * `sample_bams` - HashMap of sample names to BAM file paths
///
/// # Returns
/// * HashMap where keys are transcript IDs and values are HashMaps of sample to TPM
pub fn create_expression_matrix(
    transcripts: &[Transcript],
    sample_bams: &HashMap<String, String>,
) -> Result<HashMap<usize, HashMap<String, f64>>> {
    let mut matrix = HashMap::new();

    // Initialize matrix entries
    for transcript in transcripts {
        let mut sample_values = HashMap::new();
        for (sample, _) in sample_bams {
            sample_values.insert(sample.clone(), 0.0);
        }
        matrix.insert(transcript.id, sample_values);
    }

    // Fill in TPM values for each sample
    for (sample, bam_path) in sample_bams {
        let mut transcript_copies = transcripts.to_vec();
        update_transcripts_with_tpm(&mut transcript_copies, bam_path)?;

        for transcript in transcript_copies {
            if let Some(tpm) = transcript.tpm {
                if let Some(sample_map) = matrix.get_mut(&transcript.id) {
                    sample_map.insert(sample.clone(), tpm);
                }
            }
        }
    }

    Ok(matrix)
}

/// Loads a CSV file with alignment counts for multiple samples
///
/// # Arguments
/// * `csv_path` - Path to CSV file with counts
/// * `transcripts` - Vector of transcript objects
///
/// # Returns
/// * HashMap mapping sample names to vectors of counts
pub fn load_counts_matrix(
    csv_path: &str,
    transcripts: &[Transcript],
) -> Result<HashMap<String, Vec<f64>>> {
    let file = File::open(csv_path)?;
    let reader = BufReader::new(file);
    let transcript_index: HashMap<usize, usize> = transcripts
        .iter()
        .enumerate()
        .map(|(idx, tx)| (tx.id, idx))
        .collect();

    let mut counts_matrix = HashMap::new();
    let mut delimiter = ',';
    let mut sample_names: Vec<String> = Vec::new();
    let mut saw_header = false;

    for (line_idx, line_result) in reader.lines().enumerate() {
        let line_number = line_idx + 1;
        let line = line_result?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if !saw_header {
            delimiter = if trimmed.contains('\t') { '\t' } else { ',' };
            let fields: Vec<&str> = trimmed.split(delimiter).map(str::trim).collect();
            if fields.len() < 2 {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    format!(
                        "counts matrix header must include transcript id and at least one sample at line {}",
                        line_number
                    ),
                ));
            }

            for sample_name in fields.into_iter().skip(1) {
                if sample_name.is_empty() {
                    return Err(Error::new(
                        ErrorKind::InvalidData,
                        format!(
                            "counts matrix header contains an empty sample name at line {}",
                            line_number
                        ),
                    ));
                }
                sample_names.push(sample_name.to_string());
                counts_matrix.insert(sample_name.to_string(), vec![0.0; transcripts.len()]);
            }

            saw_header = true;
            continue;
        }

        let fields: Vec<&str> = trimmed.split(delimiter).map(str::trim).collect();
        if fields.is_empty() || fields[0].is_empty() {
            return Err(Error::new(
                ErrorKind::InvalidData,
                format!("missing transcript id at line {}", line_number),
            ));
        }

        let transcript_id = fields[0].parse::<usize>().map_err(|_| {
            Error::new(
                ErrorKind::InvalidData,
                format!(
                    "invalid transcript id '{}' at line {}",
                    fields[0], line_number
                ),
            )
        })?;

        let Some(&tx_idx) = transcript_index.get(&transcript_id) else {
            continue;
        };

        for (sample_idx, sample_name) in sample_names.iter().enumerate() {
            let raw = fields.get(sample_idx + 1).copied().unwrap_or("");
            if raw.is_empty() {
                continue;
            }

            let value = raw.parse::<f64>().map_err(|_| {
                Error::new(
                    ErrorKind::InvalidData,
                    format!(
                        "invalid count '{}' for sample '{}' at line {}",
                        raw, sample_name, line_number
                    ),
                )
            })?;
            if !value.is_finite() || value < 0.0 {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    format!(
                        "count must be finite and non-negative for sample '{}' at line {}",
                        sample_name, line_number
                    ),
                ));
            }

            if let Some(sample_counts) = counts_matrix.get_mut(sample_name) {
                sample_counts[tx_idx] += value;
            }
        }
    }

    if !saw_header {
        return Err(Error::new(
            ErrorKind::InvalidData,
            "counts matrix is missing a header row",
        ));
    }

    Ok(counts_matrix)
}

/// Generate counts matrix from TPM values
///
/// # Arguments
/// * `transcripts` - Vector of transcript objects
/// * `tpm_values` - Vector of TPM values
///
/// # Returns
/// * Vector of count values
pub fn tpm_to_counts(transcripts: &[Transcript], tpm_values: &[f64]) -> Vec<u32> {
    assert_eq!(transcripts.len(), tpm_values.len(), "Mismatched lengths");

    let total_reads = 1_000_000; // Simulated total read count
    let effective_lengths: Vec<f64> = transcripts.iter().map(|t| t.length.max(1) as f64).collect();

    // Convert TPM to expected counts
    let counts: Vec<u32> = tpm_values
        .iter()
        .zip(effective_lengths.iter())
        .map(|(&tpm, &len)| {
            let fraction = tpm / 1_000_000.0;
            let count = fraction * total_reads as f64 * len / 1000.0;
            count.round() as u32
        })
        .collect();

    counts
}

/// Counts the number of transcripts with TPM above a threshold
///
/// # Arguments
/// * `transcripts` - Vector of transcript objects
/// * `tpm_values` - Vector of TPM values
/// * `threshold` - TPM threshold for counting
///
/// # Returns
/// * Number of transcripts above threshold
pub fn count_expressed_transcripts(
    transcripts: &[Transcript],
    tpm_values: &[f64],
    threshold: f64,
) -> Result<usize> {
    assert_eq!(transcripts.len(), tpm_values.len(), "Mismatched lengths");

    let count = tpm_values.iter().filter(|&&tpm| tpm >= threshold).count();

    Ok(count)
}

/// Estimates transcript abundances from an alignment file
///
/// # Arguments
/// * `transcripts` - Vector of transcript objects
/// * `bam_path` - Path to the BAM file with alignments
///
/// # Returns
/// * Vector of abundance values
pub fn estimate_transcript_abundances_from_alignments(
    transcripts: &[Transcript],
    _bam_path: &str,
) -> Result<Vec<f64>> {
    // This would normally parse the BAM file, but for demonstration we'll use the other function
    estimate_tpm_from_bam(_bam_path, transcripts)
}

#[cfg(test)]
mod tests {
    use super::load_counts_matrix;
    use crate::graph::transcript::Transcript;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_transcript(id: usize) -> Transcript {
        Transcript {
            id,
            sequence: "ACGT".to_string(),
            path: vec![id],
            confidence: 1.0,
            length: 4,
            strand: '+',
            tpm: None,
            splicing: "linear".to_string(),
        }
    }

    #[test]
    fn load_counts_matrix_parses_csv_and_aggregates_duplicate_ids() {
        let transcripts = vec![make_transcript(42), make_transcript(7), make_transcript(1)];
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "transcript_id,sample_a,sample_b").unwrap();
        writeln!(file, "7,1.5,2.0").unwrap();
        writeln!(file, "1,0.0,3.0").unwrap();
        writeln!(file, "42,4.0,5.0").unwrap();
        writeln!(file, "7,0.5,1.0").unwrap();
        writeln!(file, "999,100.0,100.0").unwrap();

        let observed = load_counts_matrix(file.path().to_str().unwrap(), &transcripts).unwrap();
        assert_eq!(observed["sample_a"], vec![4.0, 2.0, 0.0]);
        assert_eq!(observed["sample_b"], vec![5.0, 3.0, 3.0]);
    }

    #[test]
    fn load_counts_matrix_parses_tsv_and_treats_missing_values_as_zero() {
        let transcripts = vec![make_transcript(42), make_transcript(7)];
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "# comment").unwrap();
        writeln!(file, "transcript_id\ts1\ts2").unwrap();
        writeln!(file, "7\t1.0\t").unwrap();
        writeln!(file, "42\t2.0").unwrap();

        let observed = load_counts_matrix(file.path().to_str().unwrap(), &transcripts).unwrap();
        assert_eq!(observed["s1"], vec![2.0, 1.0]);
        assert_eq!(observed["s2"], vec![0.0, 0.0]);
    }

    #[test]
    fn load_counts_matrix_rejects_invalid_counts() {
        let transcripts = vec![make_transcript(0)];
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "transcript_id,sample").unwrap();
        writeln!(file, "0,NaN").unwrap();

        let err = load_counts_matrix(file.path().to_str().unwrap(), &transcripts).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    }
}
