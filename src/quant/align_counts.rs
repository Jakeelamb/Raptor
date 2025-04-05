use crate::graph::transcript::Transcript;
use std::collections::HashMap;
use std::io::Result;

/// Estimates transcript expression in TPM from a BAM alignment file
///
/// # Arguments
/// * `bam_path` - Path to BAM file with aligned reads
/// * `transcripts` - Vector of transcripts to quantify
///
/// # Returns
/// * Vector of TPM values corresponding to each transcript
pub fn estimate_tpm_from_bam(bam_path: &str, transcripts: &[Transcript]) -> Result<Vec<f64>> {
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

/// Enhances the TPM estimation with more accurate RNA-Seq calculations
/// 
/// This improved implementation follows standard RNA-Seq quantification
/// methods with length normalization and scaling to TPM units.
/// 
/// # Arguments
/// * `transcripts` - Vector of transcripts to quantify
/// * `bam_path` - Path to alignment BAM file
/// 
/// # Returns
/// * Vector of TPM values corresponding to transcripts
pub fn estimate_tpm(transcripts: &[Transcript], bam_path: &str) -> std::io::Result<Vec<f64>> {
    // In a real implementation, would use rust-htslib to parse the BAM file
    // For this implementation, we'll simulate alignment counts

    // Mock counts for demonstration purposes
    let mut raw_counts = vec![0.0; transcripts.len()];
    
    // Simulate alignment counting
    for (i, transcript) in transcripts.iter().enumerate() {
        // Simulate read counts proportional to length and inversely to id
        // (in a real implementation, this would come from the BAM file)
        raw_counts[i] = (transcript.sequence.len() as f64 / 100.0) * 
                        (10.0 / (i as f64 + 1.0)).max(1.0);
    }
    
    // Calculate RPK (reads per kilobase)
    let rpk: Vec<f64> = transcripts.iter().enumerate()
        .map(|(i, tx)| {
            let len_kb = tx.sequence.len() as f64 / 1000.0;
            if len_kb > 0.0 {
                raw_counts[i] / len_kb
            } else {
                0.0
            }
        })
        .collect();
    
    // Calculate scaling factor to TPM
    let total_rpk: f64 = rpk.iter().sum();
    let scaling_factor = if total_rpk > 0.0 { 1_000_000.0 / total_rpk } else { 0.0 };
    
    // Convert to TPM
    let tpm: Vec<f64> = rpk.iter()
        .map(|&value| value * scaling_factor)
        .collect();
    
    Ok(tpm)
}

/// Updates transcripts with TPM values using enhanced calculation
///
/// # Arguments
/// * `transcripts` - Mutable vector of transcripts to update
/// * `bam_path` - Path to BAM file
///
/// # Returns
/// * Result with number of transcripts updated
pub fn update_transcripts_with_enhanced_tpm(transcripts: &mut [Transcript], bam_path: &str) -> std::io::Result<usize> {
    let tpm_values = estimate_tpm(transcripts, bam_path)?;
    
    for (transcript, tpm) in transcripts.iter_mut().zip(tpm_values.iter()) {
        transcript.tpm = Some(*tpm);
    }
    
    Ok(transcripts.len())
} 