use crate::graph::assembler::Contig;
use crate::graph::transcript::Transcript;
use crate::kmer::rle::rle_encode;
use serde::Serialize;
use std::fs::File;
use std::io;

/// Represents metadata for a contig
#[derive(Serialize, Debug)]
pub struct ContigMetadata {
    pub id: usize,
    pub length: usize,
    pub rle_compression: f64,
    pub gc_content: f64,
}

/// Represents metadata and metrics for a transcript
#[derive(Serialize, Debug)]
pub struct TranscriptMetrics {
    pub id: String,
    pub length: usize,
    pub confidence: f64,
    pub tpm: f64,
}

/// Generate metadata for a list of contigs
pub fn generate_metadata(contigs: &[Contig]) -> Vec<ContigMetadata> {
    contigs.iter().enumerate().map(|(i, c)| {
        let rle = rle_encode(&c.sequence);
        let rle_len = rle.len();
        let orig_len = c.sequence.len();
        
        // Calculate GC content
        let gc_count = c.sequence.bytes()
            .filter(|&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
            .count();
        
        ContigMetadata {
            id: i + 1,
            length: orig_len,
            rle_compression: if orig_len > 0 { 1.0 - (rle_len as f64 / orig_len as f64) } else { 0.0 },
            gc_content: if orig_len > 0 { gc_count as f64 / orig_len as f64 } else { 0.0 },
        }
    }).collect()
}

/// Write transcript metrics to a JSON file
pub fn write_transcript_metrics(
    transcripts: &[Transcript],
    tpms: &[f64],
    output: &str,
) -> io::Result<()> {
    let mut metrics = Vec::with_capacity(transcripts.len());
    for (tx, &tpm) in transcripts.iter().zip(tpms) {
        metrics.push(TranscriptMetrics {
            id: format!("transcript_{}", tx.id),
            length: tx.sequence.len(),
            confidence: tx.confidence,
            tpm,
        });
    }
    
    let file = File::create(output)?;
    serde_json::to_writer_pretty(file, &metrics)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_generate_metadata() {
        let contigs = vec![
            Contig {
                id: 0,
                sequence: "AAAACCCCCGGGGTTTT".to_string(),
                kmer_path: vec!["AAAAC".to_string()],
            },
            Contig {
                id: 1,
                sequence: "ATATATAT".to_string(),
                kmer_path: vec!["ATATA".to_string()],
            },
        ];
        
        let metadata = generate_metadata(&contigs);
        
        assert_eq!(metadata.len(), 2);
        
        // First contig
        assert_eq!(metadata[0].id, 1);
        assert_eq!(metadata[0].length, 17);
        
        // Check RLE compression: Original is "AAAACCCCCGGGGTTTT", 
        // RLE is [(A,4), (C,5), (G,4), (T,4)], so 4 elements vs 17 original
        // AAAACCCCCGGGGTTTT has 4 elements in RLE: A4, C5, G4, T4
        let expected_compression = 1.0 - (4.0 / 17.0);
        println!("Expected: {}, Actual: {}", expected_compression, metadata[0].rle_compression);
        assert!((metadata[0].rle_compression - expected_compression).abs() < 0.001);
        
        // Check GC content: "AAAACCCCCGGGGTTTT" has 9 G/C out of 17 total
        let expected_gc = 9.0 / 17.0;
        assert!((metadata[0].gc_content - expected_gc).abs() < 0.001);
        
        // Second contig
        assert_eq!(metadata[1].id, 2);
        assert_eq!(metadata[1].length, 8);
    }

    // Create and return mock contigs for testing
    fn create_test_contigs() -> Vec<Contig> {
        vec![
            Contig {
                id: 0,
                sequence: "ATCGATCGATCG".to_string(),
                kmer_path: vec!["ATC".to_string(), "TCG".to_string()],
            },
            Contig {
                id: 1,
                sequence: "GCTAGCTAGCT".to_string(),
                kmer_path: vec!["GCT".to_string(), "CTA".to_string()],
            },
        ]
    }
} 