use crate::graph::transcript::Transcript;
use std::collections::HashMap;
use std::io::Result;

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
    let rpk: Vec<f64> = counts.iter().zip(lengths.iter())
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
pub fn update_transcripts_with_tpm(transcripts: &mut [Transcript], bam_path: &str) -> Result<usize> {
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
    transcripts.into_iter()
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
    sample_bams: &HashMap<String, String>
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
pub fn load_counts_matrix(csv_path: &str, transcripts: &[Transcript]) -> Result<HashMap<String, Vec<f64>>> {
    let mut counts_matrix = HashMap::new();
    
    // In a real implementation, this would parse the CSV file
    // For demonstration, we'll create random counts for samples
    let sample_names = vec!["sample1", "sample2", "sample3"];
    
    for sample in &sample_names {
        let counts: Vec<f64> = transcripts.iter()
            .map(|t| {
                // Random count based on transcript length
                (t.length as f64 / 100.0) * ((t.id as f64 % 5.0) + 0.5)
            })
            .collect();
        
        counts_matrix.insert(sample.to_string(), counts);
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
    let effective_lengths: Vec<f64> = transcripts.iter()
        .map(|t| t.length.max(1) as f64)
        .collect();
    
    // Convert TPM to expected counts
    let counts: Vec<u32> = tpm_values.iter()
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
pub fn count_expressed_transcripts(transcripts: &[Transcript], tpm_values: &[f64], threshold: f64) -> Result<usize> {
    assert_eq!(transcripts.len(), tpm_values.len(), "Mismatched lengths");
    
    let count = tpm_values.iter()
        .filter(|&&tpm| tpm >= threshold)
        .count();
    
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
pub fn estimate_transcript_abundances_from_alignments(transcripts: &[Transcript], _bam_path: &str) -> Result<Vec<f64>> {
    // This would normally parse the BAM file, but for demonstration we'll use the other function
    estimate_tpm_from_bam(_bam_path, transcripts)
} 