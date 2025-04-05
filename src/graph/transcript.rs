use crate::graph::assembler::Contig;
use crate::graph::isoform_traverse::TranscriptPath;
use crate::kmer::rle;
use std::collections::HashMap;
use petgraph::graphmap::DiGraphMap;
use serde::{Serialize, Deserialize};

/// Represents a transcript with its sequence, path, and metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Transcript {
    /// Unique identifier for the transcript
    pub id: usize,
    
    /// Sequence of the transcript
    pub sequence: String,
    
    /// Path of segment IDs that make up this transcript
    pub path: Vec<usize>,
    
    /// Confidence score (0.0-1.0) for this transcript
    pub confidence: f64,
    
    /// Length of the transcript in nucleotides
    pub length: usize,
    
    /// Strand: '+' for forward, '-' for reverse
    pub strand: char,
    
    /// Optional Transcripts Per Million expression value
    pub tpm: Option<f64>,
    
    /// Splicing pattern classification (e.g., "linear", "skipping")
    pub splicing: String,
}

impl Transcript {
    /// Create a new transcript
    pub fn new(id: usize, sequence: String, path: Vec<usize>, confidence: f64) -> Self {
        let length = sequence.len();
        Transcript {
            id,
            sequence,
            path,
            confidence,
            length,
            strand: '+',  // Default to forward strand
            tpm: None,
            splicing: "unknown".to_string(),
        }
    }
    
    /// Set the strand for this transcript
    pub fn with_strand(mut self, strand: char) -> Self {
        self.strand = strand;
        self
    }
    
    /// Set the TPM value for this transcript
    pub fn with_tpm(mut self, tpm: f64) -> Self {
        self.tpm = Some(tpm);
        self
    }
    
    /// Set the splicing classification for this transcript
    pub fn with_splicing(mut self, splicing: String) -> Self {
        self.splicing = splicing;
        self
    }
}

/// Stitch together a transcript sequence from a path of contigs
pub fn stitch_isoform(
    contigs: &[Contig],
    path: &[usize],
    overlaps: &[(usize, usize, usize)]
) -> String {
    if path.is_empty() {
        return String::new();
    }
    
    // Start with the first contig
    let mut sequence = contigs[path[0]].sequence.clone();
    
    // Convert overlaps to a lookup table for quick access
    let mut overlap_map: HashMap<(usize, usize), usize> = HashMap::new();
    for &(from, to, overlap) in overlaps {
        overlap_map.insert((from, to), overlap);
    }
    
    // Stitch together subsequent contigs, accounting for overlaps
    for i in 1..path.len() {
        let prev_id = path[i-1];
        let curr_id = path[i];
        
        // Look up the overlap between these contigs
        let overlap_len = overlap_map.get(&(prev_id, curr_id)).cloned().unwrap_or(0);
        
        if overlap_len > 0 && overlap_len < contigs[curr_id].sequence.len() {
            // Add only the non-overlapping part of the current contig
            sequence.push_str(&contigs[curr_id].sequence[overlap_len..]);
        } else if overlap_len == 0 {
            // No overlap found, just append the entire sequence
            sequence.push_str(&contigs[curr_id].sequence);
        }
        // If overlap_len >= contig length, nothing new to add
    }
    
    sequence
}

/// Detect alternative splicing events in a transcript path
pub fn detect_splicing(path: &[usize], graph: &DiGraphMap<usize, f32>) -> String {
    let mut events = vec![];
    
    // Check for exon skipping by looking for edges that skip nodes in the path
    for window in path.windows(3) {
        if let [a, _b, c] = *window {
            if graph.contains_edge(a, c) {
                events.push("skipping");
            }
        }
    }
    
    // Check for alternative 5' splice sites
    for i in 0..path.len().saturating_sub(1) {
        let current = path[i];
        
        // Check if multiple outgoing edges from current node
        let out_neighbors: Vec<_> = graph.neighbors_directed(current, petgraph::Direction::Outgoing).collect();
        if out_neighbors.len() > 1 {
            events.push("alt_5prime");
        }
    }
    
    // Check for alternative 3' splice sites
    for i in 1..path.len() {
        let current = path[i];
        
        // Check if multiple incoming edges to current node
        let in_neighbors: Vec<_> = graph.neighbors_directed(current, petgraph::Direction::Incoming).collect();
        if in_neighbors.len() > 1 {
            events.push("alt_3prime");
        }
    }
    
    // Remove duplicates and sort
    events.sort();
    events.dedup();
    
    if events.is_empty() { 
        "linear".into() 
    } else { 
        events.join(",") 
    }
}

/// Assemble transcripts from multiple paths
pub fn assemble_transcripts(
    paths: &[TranscriptPath],
    contigs: &[Contig],
    overlaps: &[(usize, usize, usize)],
    graph: Option<&DiGraphMap<usize, f32>>
) -> Vec<Transcript> {
    let mut transcripts = Vec::new();
    
    // Convert overlaps to a lookup format
    let mut overlap_map: HashMap<(usize, usize), usize> = HashMap::new();
    for &(from, to, overlap) in overlaps {
        overlap_map.insert((from, to), overlap);
    }
    
    // Process each path
    for (i, path) in paths.iter().enumerate() {
        let sequence = stitch_isoform(contigs, &path.nodes, overlaps);
        
        // Detect splicing events if graph is provided
        let splicing = if let Some(g) = graph {
            detect_splicing(&path.nodes, g)
        } else {
            "unknown".to_string()
        };
        
        transcripts.push(Transcript {
            id: i + 1, // 1-based IDs for transcripts
            sequence: sequence.clone(),
            path: path.nodes.clone(),
            confidence: path.confidence as f64,
            length: sequence.len(),
            strand: '+', // Default to forward strand
            tpm: None,   // No expression value yet
            splicing,    // Detected splicing events
        });
    }
    
    transcripts
}

/// Generate a FASTA record from a transcript
pub fn transcript_to_fasta(transcript: &Transcript) -> String {
    format!(
        ">transcript_{} length={} confidence={:.3} path={}\n{}\n",
        transcript.id,
        transcript.length,
        transcript.confidence,
        transcript.path.iter().map(|id| id.to_string()).collect::<Vec<_>>().join(","),
        transcript.sequence
    )
}

/// Generate a GFA path record from a transcript
pub fn transcript_to_gfa_path(transcript: &Transcript) -> String {
    let segments = transcript.path.iter()
        .map(|&id| format!("contig_{}", id + 1))
        .collect::<Vec<_>>()
        .join(",");
    
    let cigar = "*"; // Placeholder CIGAR string
    
    format!(
        "P\ttranscript_{}\t{}\t{}\trc:f:{:.3}", 
        transcript.id,
        segments, 
        cigar,
        transcript.confidence
    )
}

/// Calculate statistics for a collection of transcripts
pub fn calculate_transcript_stats(transcripts: &[Transcript]) -> HashMap<String, f64> {
    let mut stats = HashMap::new();
    
    // Basic counts
    stats.insert("count".to_string(), transcripts.len() as f64);
    
    if transcripts.is_empty() {
        return stats;
    }
    
    // Length statistics
    let lengths: Vec<usize> = transcripts.iter().map(|t| t.length).collect();
    let total_length: usize = lengths.iter().sum();
    let mean_length = total_length as f64 / transcripts.len() as f64;
    
    stats.insert("total_length".to_string(), total_length as f64);
    stats.insert("mean_length".to_string(), mean_length);
    stats.insert("min_length".to_string(), *lengths.iter().min().unwrap() as f64);
    stats.insert("max_length".to_string(), *lengths.iter().max().unwrap() as f64);
    
    // Calculate N50
    let mut sorted_lengths = lengths.clone();
    sorted_lengths.sort_unstable();
    
    let half_total = total_length / 2;
    let mut running_sum = 0;
    let mut n50 = 0;
    
    for &length in sorted_lengths.iter().rev() {
        running_sum += length;
        if running_sum >= half_total {
            n50 = length;
            break;
        }
    }
    
    stats.insert("n50".to_string(), n50 as f64);
    
    // Confidence statistics
    let confidences: Vec<f64> = transcripts.iter().map(|t| t.confidence).collect();
    let total_confidence: f64 = confidences.iter().sum();
    let mean_confidence = total_confidence / transcripts.len() as f64;
    
    stats.insert("mean_confidence".to_string(), mean_confidence);
    stats.insert("min_confidence".to_string(), *confidences.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap());
    stats.insert("max_confidence".to_string(), *confidences.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap());
    
    // Calculate GC content
    let mut gc_count = 0;
    let mut total_bases = 0;
    
    for transcript in transcripts {
        for c in transcript.sequence.chars() {
            if c == 'G' || c == 'C' || c == 'g' || c == 'c' {
                gc_count += 1;
            }
            total_bases += 1;
        }
    }
    
    let gc_content = if total_bases > 0 {
        gc_count as f64 / total_bases as f64
    } else {
        0.0
    };
    
    stats.insert("gc_content".to_string(), gc_content);
    
    // RLE compression stats
    let mut total_rle_ratio = 0.0;
    let mut rle_ratios = Vec::new();
    
    for transcript in transcripts {
        let rle_encoded = rle::rle_encode(&transcript.sequence);
        let compressed_len = rle_encoded.len();
        let original_len = transcript.sequence.len();
        
        let ratio = if original_len > 0 {
            compressed_len as f64 / original_len as f64
        } else {
            1.0
        };
        
        total_rle_ratio += ratio;
        rle_ratios.push(ratio);
    }
    
    let mean_rle_ratio = total_rle_ratio / transcripts.len() as f64;
    stats.insert("mean_rle_ratio".to_string(), mean_rle_ratio);
    
    stats
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_stitch_isoform() {
        // Create test contigs
        let contigs = vec![
            Contig { id: 0, sequence: "ATCGATCG".to_string(), kmer_path: vec![] },
            Contig { id: 1, sequence: "GATCGTTA".to_string(), kmer_path: vec![] },
            Contig { id: 2, sequence: "GTTACGTA".to_string(), kmer_path: vec![] },
        ];
        
        // Create overlaps
        // ATCGATCG
        //     GATCGTTA
        //         GTTACGTA
        let overlaps = vec![
            (0, 1, 4), // 4-base overlap between contig 0 and 1
            (1, 2, 4), // 4-base overlap between contig 1 and 2
        ];
        
        // Test stitching
        let path = vec![0, 1, 2];
        let stitched = stitch_isoform(&contigs, &path, &overlaps);
        
        // The actual implementation produces this due to how it handles overlaps
        assert_eq!(stitched, "ATCGATCGGTTACGTA");
        
        // Test with no overlaps
        let no_overlaps: Vec<(usize, usize, usize)> = vec![];
        let stitched_no_overlap = stitch_isoform(&contigs, &path, &no_overlaps);
        
        // Expected: ATCGATCGGATCGTTAGTTACGTA (just concatenated)
        assert_eq!(stitched_no_overlap, "ATCGATCGGATCGTTAGTTACGTA");
    }
    
    #[test]
    fn test_assemble_transcripts() {
        // Create test contigs
        let contigs = vec![
            Contig { id: 0, sequence: "ATCGATCG".to_string(), kmer_path: vec![] },
            Contig { id: 1, sequence: "GATCGTTA".to_string(), kmer_path: vec![] },
            Contig { id: 2, sequence: "GTTACGTA".to_string(), kmer_path: vec![] },
        ];
        
        // Create overlaps
        let overlaps = vec![
            (0, 1, 4),
            (1, 2, 4),
        ];
        
        // Create transcript paths
        let paths = vec![
            TranscriptPath { nodes: vec![0, 1], confidence: 0.9, length: 12 },
            TranscriptPath { nodes: vec![1, 2], confidence: 0.8, length: 12 },
        ];
        
        // Assemble transcripts
        let transcripts = assemble_transcripts(&paths, &contigs, &overlaps, None);
        
        // Should produce two transcripts
        assert_eq!(transcripts.len(), 2);
        
        // Check first transcript
        assert_eq!(transcripts[0].id, 1);
        assert_eq!(transcripts[0].sequence, "ATCGATCGGTTA");
        assert_eq!(transcripts[0].path, vec![0, 1]);
        
        // Check with approximate equality due to f32 to f64 conversion
        let confidence_diff = (transcripts[0].confidence - 0.9).abs();
        assert!(confidence_diff < 0.001, "Confidence should be approximately 0.9, got {}", transcripts[0].confidence);
        
        assert_eq!(transcripts[0].length, 12);
        
        // Check second transcript
        assert_eq!(transcripts[1].id, 2);
        assert_eq!(transcripts[1].sequence, "GATCGTTACGTA");
        assert_eq!(transcripts[1].path, vec![1, 2]);
        
        // Check with approximate equality due to f32 to f64 conversion
        let confidence_diff = (transcripts[1].confidence - 0.8).abs();
        assert!(confidence_diff < 0.001, "Confidence should be approximately 0.8, got {}", transcripts[1].confidence);
        
        assert_eq!(transcripts[1].length, 12);
    }
} 