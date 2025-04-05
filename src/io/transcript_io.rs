use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use serde_json::{json, Value as JsonValue};

use crate::graph::transcript::Transcript;

/// Writer for transcript FASTA files
pub struct TranscriptFastaWriter {
    writer: BufWriter<File>,
}

impl TranscriptFastaWriter {
    /// Create a new FASTA writer for transcripts
    pub fn new(path: &str) -> Self {
        let file = File::create(path).expect("Failed to create FASTA file");
        TranscriptFastaWriter {
            writer: BufWriter::new(file),
        }
    }
    
    /// Write transcripts to FASTA format
    pub fn write_transcripts(&mut self, transcripts: &[Transcript]) -> io::Result<()> {
        for transcript in transcripts {
            writeln!(
                self.writer,
                ">transcript_{} length={} confidence={:.4}",
                transcript.id,
                transcript.length,
                transcript.confidence
            )?;
            
            // Write sequence with 60 characters per line
            let seq = &transcript.sequence;
            for i in (0..seq.len()).step_by(60) {
                let end = std::cmp::min(i + 60, seq.len());
                writeln!(self.writer, "{}", &seq[i..end])?;
            }
        }
        
        self.writer.flush()?;
        Ok(())
    }
}

/// Writer for transcript GTF files
pub struct TranscriptGtfWriter {
    writer: BufWriter<File>,
}

impl TranscriptGtfWriter {
    /// Create a new GTF writer for transcripts
    pub fn new(path: &str) -> Self {
        let file = File::create(path).expect("Failed to create GTF file");
        TranscriptGtfWriter {
            writer: BufWriter::new(file),
        }
    }
    
    /// Write transcripts to GTF format
    pub fn write_transcripts(&mut self, transcripts: &[Transcript]) -> io::Result<()> {
        for transcript in transcripts {
            let transcript_id = format!("transcript_{}", transcript.id);
            let gene_id = format!("gene_{}", transcript.id);
            
            // Write transcript entry
            writeln!(
                self.writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                transcript_id, // seqname
                "RNAtools",  // source
                "transcript", // feature
                1,            // start
                transcript.length, // end
                transcript.confidence, // score
                transcript.strand,  // strand (now using the strand field)
                ".",          // frame
                format!("gene_id \"{}\"; transcript_id \"{}\"; confidence \"{:.4}\"; tpm \"{:.3}\"; splicing \"{}\";",
                    gene_id,
                    transcript_id,
                    transcript.confidence,
                    transcript.tpm.unwrap_or(0.0),
                    transcript.splicing
                )
            )?;
            
            // Write exon entry (simplified - just one exon for the whole transcript)
            writeln!(
                self.writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                transcript_id, // seqname
                "RNAtools",  // source
                "exon",       // feature
                1,            // start
                transcript.length, // end
                transcript.confidence, // score
                transcript.strand,   // strand - use transcript strand
                ".",          // frame
                format!("gene_id \"{}\"; transcript_id \"{}\"; exon_number \"1\";", 
                    gene_id, 
                    transcript_id
                )
            )?;
        }
        
        self.writer.flush()?;
        Ok(())
    }
}

/// Add transcript paths to an existing GFA file
pub fn add_transcripts_to_gfa(gfa_path: &str, transcripts: &[Transcript]) -> io::Result<()> {
    // Check if the GFA file exists - if not, create it
    let file = if Path::new(gfa_path).exists() {
        std::fs::OpenOptions::new()
            .append(true)
            .open(gfa_path)?
    } else {
        // Create a new file with header
        let f = File::create(gfa_path)?;
        {
            let mut writer = BufWriter::new(&f);
            writeln!(writer, "H\tVN:Z:1.0")?;
            writer.flush()?;
        }
        f
    };
    
    let mut writer = BufWriter::new(file);
    
    // Add path lines for each transcript
    for transcript in transcripts {
        let path_id = format!("transcript_{}", transcript.id);
        let path_nodes = transcript.path.iter()
            .map(|&id| format!("{}+", id))
            .collect::<Vec<_>>()
            .join(",");
        
        // GFA path line format: P <path_id> <segment_ids> <cigar>
        // Using a placeholder CIGAR string for simplicity
        writeln!(writer, "P\t{}\t{}\t*", path_id, path_nodes)?;
    }
    
    writer.flush()?;
    Ok(())
}

/// Write transcript statistics to JSON
pub fn write_transcript_stats(path: &str, stats: &std::collections::HashMap<String, f64>) -> io::Result<()> {
    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, stats)?;
    Ok(())
}

/// Write transcript metadata to JSON
pub fn write_transcript_json(path: &str, transcripts: &[Transcript]) -> io::Result<()> {
    let file = File::create(path)?;
    let writer = BufWriter::new(file);
    
    // Create JSON array of transcript objects
    let json_array: Vec<JsonValue> = transcripts.iter().map(|t| {
        json!({
            "id": format!("transcript_{}", t.id),
            "length": t.length,
            "confidence": t.confidence,
            "path": t.path
        })
    }).collect();
    
    serde_json::to_writer_pretty(writer, &json_array)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::NamedTempFile;
    
    fn create_test_transcripts() -> Vec<Transcript> {
        vec![
            Transcript {
                id: 1,
                sequence: "GATTACA".to_string(),
                path: vec![1, 2, 3],
                confidence: 0.95,
                length: 7,
                strand: '+',
                tpm: Some(100.0),
                splicing: "linear".to_string(),
            },
            Transcript {
                id: 2,
                sequence: "ACGTACGT".to_string(),
                path: vec![4, 5],
                confidence: 0.85,
                length: 8,
                strand: '-',
                tpm: Some(75.5),
                splicing: "skipping".to_string(),
            },
        ]
    }
    
    #[test]
    fn test_write_fasta() {
        let transcripts = create_test_transcripts();
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path().to_str().unwrap();
        
        let mut writer = TranscriptFastaWriter::new(path);
        writer.write_transcripts(&transcripts).unwrap();
        
        // Check file contents
        let mut file = File::open(path).unwrap();
        let mut contents = String::new();
        file.read_to_string(&mut contents).unwrap();
        
        assert!(contents.contains(">transcript_1 length=7 confidence=0.95"));
        assert!(contents.contains("GATTACA"));
        assert!(contents.contains(">transcript_2 length=8 confidence=0.85"));
        assert!(contents.contains("ACGTACGT"));
    }
    
    #[test]
    fn test_write_gtf() {
        let transcripts = create_test_transcripts();
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path().to_str().unwrap();
        
        let mut writer = TranscriptGtfWriter::new(path);
        writer.write_transcripts(&transcripts).unwrap();
        
        // Check file contents
        let mut file = File::open(path).unwrap();
        let mut contents = String::new();
        file.read_to_string(&mut contents).unwrap();
        
        assert!(contents.contains("transcript_1\traptor\ttranscript"));
        assert!(contents.contains("transcript_2\traptor\ttranscript"));
        assert!(contents.contains("exon"));
        assert!(contents.contains("transcript_id \"transcript_1\""));
    }
    
    #[test]
    fn test_write_transcript_json() {
        let transcripts = create_test_transcripts();
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path().to_str().unwrap();
        
        write_transcript_json(path, &transcripts).unwrap();
        
        // Check file contents
        let mut file = File::open(path).unwrap();
        let mut contents = String::new();
        file.read_to_string(&mut contents).unwrap();
        
        assert!(contents.contains("\"id\": \"transcript_1\""));
        assert!(contents.contains("\"confidence\": 0.95"));
        assert!(contents.contains("\"path\": [1, 2, 3]"));
    }
} 