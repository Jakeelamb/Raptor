use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::fs::File;
use crate::graph::transcript::Transcript;
use crate::io::fastq::FastqRecord;

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

/// Polishes transcripts using long read data
///
/// This function aligns long reads to transcripts and uses them to correct
/// potential errors in the short-read assembly.
///
/// # Arguments
/// * `transcripts` - Vector of transcripts to polish
/// * `long_reads` - Vector of long read sequences
///
/// # Returns
/// * Vector of polished transcripts
pub fn polish_transcripts_with_long_reads(
    transcripts: &[Transcript], 
    long_reads: &[FastqRecord]
) -> Vec<Transcript> {
    // Clone the transcripts for modification
    let mut polished = transcripts.to_vec();
    
    // For each transcript, try to find long reads that align to it
    for transcript in &mut polished {
        // Find long reads that match this transcript
        let aligned_reads = find_aligned_reads(&transcript.sequence, long_reads);
        
        if !aligned_reads.is_empty() {
            // Polish the sequence using aligned reads
            transcript.sequence = polish_sequence(&transcript.sequence, &aligned_reads);
        }
    }
    
    polished
}

/// Finds long reads that align to a transcript
///
/// # Arguments
/// * `transcript_seq` - The transcript sequence
/// * `long_reads` - Vector of long read sequences
///
/// # Returns
/// * Vector of FastqRecord that align to the transcript
fn find_aligned_reads(transcript_seq: &str, long_reads: &[FastqRecord]) -> Vec<FastqRecord> {
    // In a real implementation, we would use a proper alignment algorithm
    // For this example, we'll use a simple substring match as a placeholder
    
    long_reads.iter()
        .filter(|read| {
            // Use a sliding window to check for partial matches
            let min_overlap = std::cmp::min(read.sequence.len(), transcript_seq.len()) / 2;
            for window_size in (min_overlap..=read.sequence.len()).rev() {
                let read_window = &read.sequence[0..window_size];
                if transcript_seq.contains(read_window) {
                    return true;
                }
            }
            false
        })
        .cloned()
        .collect()
}

/// Polishes a sequence using aligned reads
///
/// # Arguments
/// * `seq` - The sequence to polish
/// * `reads` - Aligned reads
///
/// # Returns
/// * Polished sequence
fn polish_sequence(seq: &str, reads: &[FastqRecord]) -> String {
    // Convert sequence to nucleotide array for consensus building
    let seq_len = seq.len();
    let mut consensus_counts = vec![[0u32; 5]; seq_len]; // A, C, G, T, N
    
    // For each aligned read, update the consensus counts
    for read in reads {
        // In a real implementation, we would use the actual alignment positions
        // For this example, find approximate position using substring match
        if let Some(pos) = find_approximate_position(seq, &read.sequence) {
            update_consensus_counts(&mut consensus_counts, pos, &read.sequence);
        }
    }
    
    // Build consensus sequence
    build_consensus_sequence(&consensus_counts)
}

/// Finds approximate starting position of a read in a reference sequence
///
/// # Arguments
/// * `reference` - Reference sequence
/// * `read` - Read sequence
///
/// # Returns
/// * Optional starting position
fn find_approximate_position(reference: &str, read: &str) -> Option<usize> {
    // Simple sliding window search for a match
    let window_size = std::cmp::min(30, read.len());
    let read_prefix = &read[0..window_size];
    
    reference.find(read_prefix)
}

/// Updates consensus counts based on an aligned read
///
/// # Arguments
/// * `counts` - Consensus count matrix
/// * `start_pos` - Starting position in the reference
/// * `read_seq` - Read sequence
fn update_consensus_counts(counts: &mut Vec<[u32; 5]>, start_pos: usize, read_seq: &str) {
    for (i, c) in read_seq.bytes().enumerate() {
        let ref_pos = start_pos + i;
        if ref_pos >= counts.len() {
            break;
        }
        
        let idx = match c {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 4, // N or any other character
        };
        
        counts[ref_pos][idx] += 1;
    }
}

/// Builds consensus sequence from nucleotide counts
///
/// # Arguments
/// * `counts` - Consensus count matrix
///
/// # Returns
/// * Consensus sequence
fn build_consensus_sequence(counts: &Vec<[u32; 5]>) -> String {
    counts.iter().map(|pos_counts| {
        let max_idx = pos_counts.iter()
            .enumerate()
            .max_by_key(|(_, &count)| count)
            .map(|(idx, _)| idx)
            .unwrap_or(4);
        
        match max_idx {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => 'N',
        }
    }).collect()
} 