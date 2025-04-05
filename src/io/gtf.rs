use crate::graph::transcript::Transcript;
use std::fs::File;
use std::io::{BufWriter, Write, Result};

/// Write transcripts to GTF format
///
/// # Arguments
/// * `transcripts` - Vector of transcripts to write
/// * `out_path` - Output file path
///
/// # Returns
/// * IO Result
pub fn write_gtf(transcripts: &[Transcript], out_path: &str) -> Result<()> {
    let file = File::create(out_path)?;
    let mut writer = BufWriter::new(file);
    
    for tx in transcripts {
        let tx_id = format!("transcript_{}", tx.id);
        let gene_id = format!("gene_{}", tx.id);
        
        // Write transcript feature
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            tx_id,         // seqname
            "RNAtools",    // source
            "transcript",  // feature
            1,             // start
            tx.sequence.len(), // end
            tx.confidence, // score
            tx.strand,     // strand
            ".",          // frame
            format!(
                "gene_id \"{}\"; transcript_id \"{}\"; tpm \"{:.3}\"; confidence \"{:.3}\"; splicing \"{}\";",
                gene_id, tx_id, tx.tpm.unwrap_or(0.0), tx.confidence, tx.splicing
            )
        )?;
        
        // Write exon feature (single exon for simple transcripts)
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            tx_id,         // seqname
            "RNAtools",    // source
            "exon",        // feature
            1,             // start
            tx.sequence.len(), // end
            tx.confidence, // score
            tx.strand,     // strand
            ".",          // frame
            format!(
                "gene_id \"{}\"; transcript_id \"{}\"; exon_number \"1\";",
                gene_id, tx_id
            )
        )?;
    }
    
    Ok(())
}

/// Write transcripts to GTF format with detailed exon structure
///
/// # Arguments
/// * `transcripts` - Vector of transcripts to write
/// * `out_path` - Output file path
/// * `exon_coords` - Optional vector of exon coordinates for each transcript
///
/// # Returns
/// * IO Result
pub fn write_detailed_gtf(
    transcripts: &[Transcript], 
    out_path: &str,
    exon_coords: Option<&Vec<Vec<(usize, usize)>>>
) -> Result<()> {
    let file = File::create(out_path)?;
    let mut writer = BufWriter::new(file);
    
    for (i, tx) in transcripts.iter().enumerate() {
        let tx_id = format!("transcript_{}", tx.id);
        let gene_id = format!("gene_{}", tx.id);
        
        // Write transcript feature
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            tx_id,         // seqname
            "RNAtools",    // source
            "transcript",  // feature
            1,             // start
            tx.sequence.len(), // end
            tx.confidence, // score
            tx.strand,     // strand
            ".",          // frame
            format!(
                "gene_id \"{}\"; transcript_id \"{}\"; tpm \"{:.3}\"; confidence \"{:.3}\"; splicing \"{}\";",
                gene_id, tx_id, tx.tpm.unwrap_or(0.0), tx.confidence, tx.splicing
            )
        )?;
        
        // Check if we have detailed exon coordinates
        if let Some(coords) = exon_coords {
            if i < coords.len() && !coords[i].is_empty() {
                // Write each exon with its coordinates
                for (exon_idx, (start, end)) in coords[i].iter().enumerate() {
                    writeln!(
                        writer,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        tx_id,         // seqname
                        "RNAtools",    // source
                        "exon",        // feature
                        start,         // start
                        end,           // end
                        tx.confidence, // score
                        tx.strand,     // strand
                        ".",          // frame
                        format!(
                            "gene_id \"{}\"; transcript_id \"{}\"; exon_number \"{}\";",
                            gene_id, tx_id, exon_idx + 1
                        )
                    )?;
                }
                continue;
            }
        }
        
        // Write a single exon if no detailed coordinates
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            tx_id,         // seqname
            "RNAtools",    // source
            "exon",        // feature
            1,             // start
            tx.sequence.len(), // end
            tx.confidence, // score
            tx.strand,     // strand
            ".",          // frame
            format!(
                "gene_id \"{}\"; transcript_id \"{}\"; exon_number \"1\";",
                gene_id, tx_id
            )
        )?;
    }
    
    Ok(())
} 