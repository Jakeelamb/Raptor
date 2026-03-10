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
    contigs
        .iter()
        .enumerate()
        .map(|(i, c)| {
            let rle = rle_encode(&c.sequence);
            let rle_len = rle.len();
            let orig_len = c.sequence.len();

            // Calculate GC content
            let gc_count = c
                .sequence
                .bytes()
                .filter(|&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
                .count();

            ContigMetadata {
                id: i + 1,
                length: orig_len,
                rle_compression: if orig_len > 0 {
                    1.0 - (rle_len as f64 / orig_len as f64)
                } else {
                    0.0
                },
                gc_content: if orig_len > 0 {
                    gc_count as f64 / orig_len as f64
                } else {
                    0.0
                },
            }
        })
        .collect()
}

/// Write transcript metrics to a JSON file
pub fn write_transcript_metrics(
    transcripts: &[Transcript],
    tpms: &[f64],
    output: &str,
) -> io::Result<()> {
    let mut ordered: Vec<(usize, &Transcript)> = transcripts.iter().enumerate().collect();
    ordered.sort_unstable_by(|a, b| a.1.id.cmp(&b.1.id).then_with(|| a.0.cmp(&b.0)));

    let mut metrics = Vec::with_capacity(transcripts.len());
    for (idx, tx) in ordered {
        let tpm = tpms.get(idx).copied().unwrap_or(0.0);
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
    use crate::graph::transcript::Transcript;
    use crate::kmer::kmer::encode_kmer;
    use tempfile::NamedTempFile;

    #[test]
    fn test_generate_metadata() {
        let contigs = vec![
            Contig {
                id: 0,
                sequence: "AAAACCCCCGGGGTTTT".to_string(),
                kmer_path: vec![encode_kmer("AAAAC").unwrap()],
            },
            Contig {
                id: 1,
                sequence: "ATATATAT".to_string(),
                kmer_path: vec![encode_kmer("ATATA").unwrap()],
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
        println!(
            "Expected: {}, Actual: {}",
            expected_compression, metadata[0].rle_compression
        );
        assert!((metadata[0].rle_compression - expected_compression).abs() < 0.001);

        // Check GC content: "AAAACCCCCGGGGTTTT" has 9 G/C out of 17 total
        let expected_gc = 9.0 / 17.0;
        assert!((metadata[0].gc_content - expected_gc).abs() < 0.001);

        // Second contig
        assert_eq!(metadata[1].id, 2);
        assert_eq!(metadata[1].length, 8);
    }

    // Create and return mock contigs for testing
    #[allow(dead_code)]
    fn create_test_contigs() -> Vec<Contig> {
        vec![
            Contig {
                id: 0,
                sequence: "ATCGATCGATCG".to_string(),
                kmer_path: vec![encode_kmer("ATC").unwrap(), encode_kmer("TCG").unwrap()],
            },
            Contig {
                id: 1,
                sequence: "GCTAGCTAGCT".to_string(),
                kmer_path: vec![encode_kmer("GCT").unwrap(), encode_kmer("CTA").unwrap()],
            },
        ]
    }

    fn make_transcript(id: usize, sequence: &str, confidence: f64) -> Transcript {
        Transcript {
            id,
            sequence: sequence.to_string(),
            path: vec![id],
            confidence,
            length: sequence.len(),
            strand: '+',
            tpm: None,
            splicing: "linear".to_string(),
        }
    }

    #[test]
    fn write_transcript_metrics_is_sorted_by_id_and_pads_missing_tpms() {
        let transcripts = vec![
            make_transcript(8, "ACGT", 0.9),
            make_transcript(3, "AAAAAA", 0.8),
            make_transcript(5, "CC", 0.7),
        ];
        let output = NamedTempFile::new().unwrap();

        write_transcript_metrics(&transcripts, &[42.5], output.path().to_str().unwrap()).unwrap();
        let json = std::fs::read_to_string(output.path()).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&json).unwrap();
        let rows = parsed.as_array().unwrap();

        assert_eq!(rows.len(), 3);
        assert_eq!(rows[0]["id"], "transcript_3");
        assert_eq!(rows[0]["tpm"], 0.0);
        assert_eq!(rows[1]["id"], "transcript_5");
        assert_eq!(rows[1]["tpm"], 0.0);
        assert_eq!(rows[2]["id"], "transcript_8");
        assert_eq!(rows[2]["tpm"], 42.5);
    }
}
