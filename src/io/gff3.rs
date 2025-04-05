use crate::graph::transcript::Transcript;
use std::fs::File;
use std::io::{BufWriter, Write, Result};

/// Write transcripts to GFF3 format
pub fn write_gff3(transcripts: &[Transcript], output: &str) -> Result<()> {
    let f = File::create(output)?;
    let mut writer = BufWriter::new(f);
    
    // Write GFF3 header
    writeln!(writer, "##gff-version 3")?;
    
    for tx in transcripts {
        let start = 1;
        let end = tx.sequence.len();
        let gene_id = format!("gene_{}", tx.id);
        let transcript_id = format!("transcript_{}", tx.id);
        
        // Write gene feature
        writeln!(
            writer,
            "{}\tAssembler\tgene\t{}\t{}\t.\t+\t.\tID={};Name={}",
            tx.id, start, end, gene_id, gene_id
        )?;
        
        // Write mRNA feature
        writeln!(
            writer,
            "{}\tAssembler\tmRNA\t{}\t{}\t.\t+\t.\tID={};Parent={};confidence={:.3}",
            tx.id, start, end, transcript_id, gene_id, tx.confidence
        )?;
        
        // Write exon feature
        writeln!(
            writer,
            "{}\tAssembler\texon\t{}\t{}\t.\t+\t.\tParent={}",
            tx.id, start, end, transcript_id
        )?;
    }
    
    Ok(())
}

/// Write transcripts to GFF3 format with multiple exons based on contig paths
pub fn write_detailed_gff3(
    transcripts: &[Transcript], 
    output: &str,
    contigs: &[usize],  // Positions where contigs start in each transcript
) -> Result<()> {
    let f = File::create(output)?;
    let mut writer = BufWriter::new(f);
    
    // Write GFF3 header
    writeln!(writer, "##gff-version 3")?;
    
    for tx in transcripts {
        let gene_id = format!("gene_{}", tx.id);
        let transcript_id = format!("transcript_{}", tx.id);
        
        // Write gene feature
        writeln!(
            writer,
            "{}\tAssembler\tgene\t1\t{}\t.\t+\t.\tID={};Name={}",
            tx.id, tx.length, gene_id, gene_id
        )?;
        
        // Write mRNA feature
        writeln!(
            writer,
            "{}\tAssembler\tmRNA\t1\t{}\t.\t+\t.\tID={};Parent={};confidence={:.3}",
            tx.id, tx.length, transcript_id, gene_id, tx.confidence
        )?;
        
        // If we have contig positions, write multiple exons
        if !contigs.is_empty() {
            // Create exon coordinates from contig positions
            let mut exon_starts = vec![1];
            exon_starts.extend(contigs);
            
            let mut exon_ends = contigs.to_vec();
            exon_ends.push(tx.length);
            
            // Write each exon
            for (i, (start, end)) in exon_starts.iter().zip(exon_ends.iter()).enumerate() {
                writeln!(
                    writer,
                    "{}\tAssembler\texon\t{}\t{}\t.\t+\t.\tID=exon_{}_{};Parent={}",
                    tx.id, start, end, tx.id, i+1, transcript_id
                )?;
            }
        } else {
            // Write a single exon for the entire transcript
            writeln!(
                writer,
                "{}\tAssembler\texon\t1\t{}\t.\t+\t.\tID=exon_{}_1;Parent={}",
                tx.id, tx.length, tx.id, transcript_id
            )?;
        }
    }
    
    Ok(())
} 