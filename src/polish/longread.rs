use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::fs::File;
use crate::graph::transcript::Transcript;

/// Parse mapped reads from SAM (transcript_id → sequences)
pub fn group_alignments_by_transcript(sam_path: &str) -> Result<HashMap<String, Vec<String>>, std::io::Error> {
    let mut map: HashMap<String, Vec<String>> = HashMap::new();
    let file = File::open(sam_path)?;
    
    for line in BufReader::new(file).lines() {
        let line = line?;
        if line.starts_with('@') { continue; } // Skip header lines
        
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 || fields[2] == "*" { continue; } // Skip unmapped or invalid alignments
        
        let tx_id = fields[2].to_string();
        let seq = fields[9].to_string();
        
        if !seq.is_empty() {
            map.entry(tx_id).or_default().push(seq);
        }
    }
    
    Ok(map)
}

/// Polish a transcript by consensus majority vote per base
pub fn polish_transcript(seq: &str, reads: &[String]) -> String {
    // Early return if no reads are available
    if reads.is_empty() {
        return seq.to_string();
    }

    // Create count matrix for each base position (A, C, G, T)
    let mut counts = vec![[0u32; 4]; seq.len()];
    
    // Count nucleotide occurrences at each position
    for read in reads {
        let len = read.len().min(seq.len());
        for (i, c) in read[..len].bytes().enumerate() {
            match c {
                b'A' | b'a' => counts[i][0] += 1,
                b'C' | b'c' => counts[i][1] += 1,
                b'G' | b'g' => counts[i][2] += 1,
                b'T' | b't' => counts[i][3] += 1,
                _ => {}
            }
        }
    }

    // Build polished sequence by choosing the most frequent base at each position
    let polished: String = counts.iter().enumerate().map(|(i, arr)| {
        // If no reads cover this position, keep original base
        if arr.iter().sum::<u32>() == 0 {
            return seq.chars().nth(i).unwrap_or('N');
        }
        
        // Otherwise choose the most frequent base
        match arr.iter().enumerate().max_by_key(|&(_, &v)| v) {
            Some((0, _)) => 'A',
            Some((1, _)) => 'C',
            Some((2, _)) => 'G',
            Some((3, _)) => 'T',
            _ => 'N', // Fallback (should not happen)
        }
    }).collect();

    polished
}

/// Polish a set of transcripts using long read alignments
pub fn polish_transcripts(
    transcripts: &mut [Transcript],
    sam_path: &str
) -> Result<usize, std::io::Error> {
    let alignments = group_alignments_by_transcript(sam_path)?;
    let mut polished_count = 0;
    
    for t in transcripts {
        let tx_id = format!("transcript_{}", t.id);
        if let Some(reads) = alignments.get(&tx_id) {
            if !reads.is_empty() {
                t.sequence = polish_transcript(&t.sequence, reads);
                polished_count += 1;
            }
        }
    }
    
    Ok(polished_count)
} 