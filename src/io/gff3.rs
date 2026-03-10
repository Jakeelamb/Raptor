use crate::graph::transcript::Transcript;
use std::fs::File;
use std::io::{BufWriter, Result, Write};

#[inline]
fn ordered_transcripts(transcripts: &[Transcript]) -> Vec<&Transcript> {
    let mut ordered: Vec<&Transcript> = transcripts.iter().collect();
    ordered.sort_unstable_by_key(|tx| tx.id);
    ordered
}

#[inline]
fn transcript_span(tx: &Transcript) -> usize {
    let seq_len = tx.sequence.len();
    if seq_len > 0 {
        seq_len
    } else {
        tx.length
    }
}

#[inline]
fn normalized_strand(strand: char) -> char {
    match strand {
        '+' | '-' | '.' | '?' => strand,
        _ => '.',
    }
}

/// Write transcripts to GFF3 format
pub fn write_gff3(transcripts: &[Transcript], output: &str) -> Result<()> {
    let f = File::create(output)?;
    let mut writer = BufWriter::new(f);

    // Write GFF3 header
    writeln!(writer, "##gff-version 3")?;

    for tx in ordered_transcripts(transcripts) {
        let start = 1;
        let end = transcript_span(tx);
        let strand = normalized_strand(tx.strand);
        let gene_id = format!("gene_{}", tx.id);
        let transcript_id = format!("transcript_{}", tx.id);

        // Write gene feature
        writeln!(
            writer,
            "{}\tRaptor\tgene\t{}\t{}\t.\t{}\t.\tID={};Name={}",
            tx.id, start, end, strand, gene_id, gene_id
        )?;

        // Write mRNA feature
        writeln!(
            writer,
            "{}\tRaptor\tmRNA\t{}\t{}\t.\t{}\t.\tID={};Parent={};confidence={:.3}",
            tx.id, start, end, strand, transcript_id, gene_id, tx.confidence
        )?;

        // Write exon feature
        writeln!(
            writer,
            "{}\tRaptor\texon\t{}\t{}\t.\t{}\t.\tParent={}",
            tx.id, start, end, strand, transcript_id
        )?;
    }

    Ok(())
}

/// Write transcripts to GFF3 format with multiple exons based on contig paths
pub fn write_detailed_gff3(
    transcripts: &[Transcript],
    output: &str,
    contigs: &[usize], // Positions where contigs start in each transcript
) -> Result<()> {
    let f = File::create(output)?;
    let mut writer = BufWriter::new(f);

    // Write GFF3 header
    writeln!(writer, "##gff-version 3")?;

    for tx in ordered_transcripts(transcripts) {
        let span = transcript_span(tx);
        let strand = normalized_strand(tx.strand);
        let gene_id = format!("gene_{}", tx.id);
        let transcript_id = format!("transcript_{}", tx.id);

        // Write gene feature
        writeln!(
            writer,
            "{}\tRaptor\tgene\t1\t{}\t.\t{}\t.\tID={};Name={}",
            tx.id, span, strand, gene_id, gene_id
        )?;

        // Write mRNA feature
        writeln!(
            writer,
            "{}\tRaptor\tmRNA\t1\t{}\t.\t{}\t.\tID={};Parent={};confidence={:.3}",
            tx.id, span, strand, transcript_id, gene_id, tx.confidence
        )?;

        // If we have contig positions, write multiple exons
        if !contigs.is_empty() {
            // Create exon coordinates from contig positions
            let mut exon_starts = vec![1];
            exon_starts.extend(contigs);

            let mut exon_ends = contigs.to_vec();
            exon_ends.push(span);

            // Write each exon
            for (i, (start, end)) in exon_starts.iter().zip(exon_ends.iter()).enumerate() {
                writeln!(
                    writer,
                    "{}\tRaptor\texon\t{}\t{}\t.\t{}\t.\tID=exon_{}_{};Parent={}",
                    tx.id,
                    start,
                    end,
                    strand,
                    tx.id,
                    i + 1,
                    transcript_id
                )?;
            }
        } else {
            // Write a single exon for the entire transcript
            writeln!(
                writer,
                "{}\tRaptor\texon\t1\t{}\t.\t{}\t.\tID=exon_{}_1;Parent={}",
                tx.id, span, strand, tx.id, transcript_id
            )?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{write_detailed_gff3, write_gff3};
    use crate::graph::transcript::Transcript;
    use std::fs;
    use tempfile::NamedTempFile;

    fn tx(id: usize, seq: &str, strand: char, confidence: f64) -> Transcript {
        Transcript {
            id,
            sequence: seq.to_string(),
            path: vec![id],
            confidence,
            length: seq.len(),
            strand,
            tpm: None,
            splicing: "linear".to_string(),
        }
    }

    #[test]
    fn write_gff3_is_sorted_by_transcript_id_and_uses_transcript_strand() {
        let transcripts = vec![
            tx(10, "AAAA", '-', 0.7),
            tx(4, "CC", 'x', 0.9),
            tx(7, "GGG", '+', 0.8),
        ];
        let output = NamedTempFile::new().unwrap();

        write_gff3(&transcripts, output.path().to_str().unwrap()).unwrap();
        let content = fs::read_to_string(output.path()).unwrap();
        let lines: Vec<&str> = content.lines().collect();

        assert_eq!(lines[0], "##gff-version 3");
        assert!(lines[1].starts_with("4\tRaptor\tgene\t"));
        assert!(lines[4].starts_with("7\tRaptor\tgene\t"));
        assert!(lines[7].starts_with("10\tRaptor\tgene\t"));

        // Invalid strand gets normalized to '.'
        assert_eq!(lines[1].split('\t').nth(6), Some("."));
        assert_eq!(lines[4].split('\t').nth(6), Some("+"));
        assert_eq!(lines[7].split('\t').nth(6), Some("-"));
    }

    #[test]
    fn write_detailed_gff3_keeps_sorted_order_and_strand() {
        let transcripts = vec![tx(9, "AAAAAA", '-', 0.6), tx(3, "CCCC", '+', 0.9)];
        let output = NamedTempFile::new().unwrap();

        write_detailed_gff3(&transcripts, output.path().to_str().unwrap(), &[2]).unwrap();
        let content = fs::read_to_string(output.path()).unwrap();
        let lines: Vec<&str> = content.lines().collect();

        assert_eq!(lines[0], "##gff-version 3");
        // transcript 3 comes first (sorted)
        assert!(lines[1].starts_with("3\tRaptor\tgene\t"));
        assert_eq!(lines[1].split('\t').nth(6), Some("+"));

        // transcript 9 section preserves '-' strand
        assert!(lines
            .iter()
            .any(|line| line.starts_with("9\tRaptor\tgene\t1\t6\t.\t-\t.")));
    }
}
