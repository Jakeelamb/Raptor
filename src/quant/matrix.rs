use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use crate::graph::transcript::Transcript;

/// Generate matrix from isoform -> TPM mapping per sample
pub fn write_counts_matrix(
    samples: &HashMap<String, Vec<f64>>,
    transcripts: &[Transcript],
    output: &str,
) {
    let mut file = std::fs::File::create(output).expect("Unable to write matrix");
    use std::io::Write;

    // Header
    write!(file, "transcript_id").unwrap();
    for sample in samples.keys() {
        write!(file, "\t{}", sample).unwrap();
    }
    writeln!(file).unwrap();

    for (i, tx) in transcripts.iter().enumerate() {
        write!(file, "transcript_{}", tx.id).unwrap();
        for tpms in samples.values() {
            write!(file, "\t{:.2}", tpms[i]).unwrap();
        }
        writeln!(file).unwrap();
    }
}

/// Read a TPM matrix file into a map of sample name -> transcript TPM values
pub fn read_tpm_matrix(path: &str) -> Result<HashMap<String, Vec<f64>>, std::io::Error> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    
    // Read header to get sample names
    let header = match lines.next() {
        Some(Ok(line)) => line,
        _ => return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "Empty or invalid matrix file")),
    };
    
    let headers: Vec<&str> = header.split('\t').collect();
    if headers.len() < 2 {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid matrix format - no sample columns found"));
    }
    
    // Initialize result map
    let mut result: HashMap<String, Vec<f64>> = HashMap::new();
    for &h in headers.iter().skip(1) { // Skip transcript_id column
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