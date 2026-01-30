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
            "{}\tRNAtools\ttranscript\t1\t{}\t{}\t{}\t.\tgene_id \"{}\"; transcript_id \"{}\"; tpm \"{:.3}\"; confidence \"{:.3}\"; splicing \"{}\";",
            tx_id,
            tx.sequence.len(),
            tx.confidence,
            tx.strand,
            gene_id, tx_id, tx.tpm.unwrap_or(0.0), tx.confidence, tx.splicing
        )?;

        // Write exon feature (single exon for simple transcripts)
        writeln!(
            writer,
            "{}\tRNAtools\texon\t1\t{}\t{}\t{}\t.\tgene_id \"{}\"; transcript_id \"{}\"; exon_number \"1\";",
            tx_id,
            tx.sequence.len(),
            tx.confidence,
            tx.strand,
            gene_id, tx_id
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
            "{}\tRNAtools\ttranscript\t1\t{}\t{}\t{}\t.\tgene_id \"{}\"; transcript_id \"{}\"; tpm \"{:.3}\"; confidence \"{:.3}\"; splicing \"{}\";",
            tx_id,
            tx.sequence.len(),
            tx.confidence,
            tx.strand,
            gene_id, tx_id, tx.tpm.unwrap_or(0.0), tx.confidence, tx.splicing
        )?;

        // Check if we have detailed exon coordinates
        if let Some(coords) = exon_coords {
            if i < coords.len() && !coords[i].is_empty() {
                // Write each exon with its coordinates
                for (exon_idx, (start, end)) in coords[i].iter().enumerate() {
                    writeln!(
                        writer,
                        "{}\tRNAtools\texon\t{}\t{}\t{}\t{}\t.\tgene_id \"{}\"; transcript_id \"{}\"; exon_number \"{}\";",
                        tx_id,
                        start,
                        end,
                        tx.confidence,
                        tx.strand,
                        gene_id, tx_id, exon_idx + 1
                    )?;
                }
                continue;
            }
        }

        // Write a single exon if no detailed coordinates
        writeln!(
            writer,
            "{}\tRNAtools\texon\t1\t{}\t{}\t{}\t.\tgene_id \"{}\"; transcript_id \"{}\"; exon_number \"1\";",
            tx_id,
            tx.sequence.len(),
            tx.confidence,
            tx.strand,
            gene_id, tx_id
        )?;
    }

    Ok(())
}
