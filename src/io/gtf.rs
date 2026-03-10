use crate::graph::transcript::Transcript;
use std::fs::File;
use std::io::{BufWriter, Result, Write};

#[inline]
fn ordered_transcripts(transcripts: &[Transcript]) -> Vec<(usize, &Transcript)> {
    let mut ordered: Vec<(usize, &Transcript)> = transcripts.iter().enumerate().collect();
    ordered.sort_unstable_by(|a, b| a.1.id.cmp(&b.1.id).then_with(|| a.0.cmp(&b.0)));
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

    for (_, tx) in ordered_transcripts(transcripts) {
        let tx_id = format!("transcript_{}", tx.id);
        let gene_id = format!("gene_{}", tx.id);
        let span = transcript_span(tx);
        let strand = normalized_strand(tx.strand);

        // Write transcript feature
        writeln!(
            writer,
            "{}\tRNAtools\ttranscript\t1\t{}\t{}\t{}\t.\tgene_id \"{}\"; transcript_id \"{}\"; tpm \"{:.3}\"; confidence \"{:.3}\"; splicing \"{}\";",
            tx_id,
            span,
            tx.confidence,
            strand,
            gene_id, tx_id, tx.tpm.unwrap_or(0.0), tx.confidence, tx.splicing
        )?;

        // Write exon feature (single exon for simple transcripts)
        writeln!(
            writer,
            "{}\tRNAtools\texon\t1\t{}\t{}\t{}\t.\tgene_id \"{}\"; transcript_id \"{}\"; exon_number \"1\";",
            tx_id,
            span,
            tx.confidence,
            strand,
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
    exon_coords: Option<&[Vec<(usize, usize)>]>,
) -> Result<()> {
    let file = File::create(out_path)?;
    let mut writer = BufWriter::new(file);

    for (original_idx, tx) in ordered_transcripts(transcripts) {
        let tx_id = format!("transcript_{}", tx.id);
        let gene_id = format!("gene_{}", tx.id);
        let span = transcript_span(tx);
        let strand = normalized_strand(tx.strand);

        // Write transcript feature
        writeln!(
            writer,
            "{}\tRNAtools\ttranscript\t1\t{}\t{}\t{}\t.\tgene_id \"{}\"; transcript_id \"{}\"; tpm \"{:.3}\"; confidence \"{:.3}\"; splicing \"{}\";",
            tx_id,
            span,
            tx.confidence,
            strand,
            gene_id, tx_id, tx.tpm.unwrap_or(0.0), tx.confidence, tx.splicing
        )?;

        // Check if we have detailed exon coordinates
        if let Some(coords) = exon_coords {
            if original_idx < coords.len() && !coords[original_idx].is_empty() {
                // Write each exon with its coordinates
                for (exon_idx, (start, end)) in coords[original_idx].iter().enumerate() {
                    writeln!(
                        writer,
                        "{}\tRNAtools\texon\t{}\t{}\t{}\t{}\t.\tgene_id \"{}\"; transcript_id \"{}\"; exon_number \"{}\";",
                        tx_id,
                        start,
                        end,
                        tx.confidence,
                        strand,
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
            span,
            tx.confidence,
            strand,
            gene_id, tx_id
        )?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{write_detailed_gtf, write_gtf};
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
            tpm: Some(1.0),
            splicing: "linear".to_string(),
        }
    }

    #[test]
    fn write_gtf_is_sorted_by_transcript_id_and_normalizes_strand() {
        let transcripts = vec![
            tx(9, "AAAA", '-', 0.7),
            tx(2, "CCC", 'x', 0.9),
            tx(5, "GG", '+', 0.8),
        ];
        let output = NamedTempFile::new().unwrap();

        write_gtf(&transcripts, output.path().to_str().unwrap()).unwrap();
        let content = fs::read_to_string(output.path()).unwrap();
        let lines: Vec<&str> = content.lines().collect();

        // Transcript rows should be sorted by transcript ID: 2,5,9.
        assert!(lines[0].starts_with("transcript_2\t"));
        assert!(lines[2].starts_with("transcript_5\t"));
        assert!(lines[4].starts_with("transcript_9\t"));

        // Invalid strand is normalized to '.'.
        assert_eq!(lines[0].split('\t').nth(6), Some("."));
        assert_eq!(lines[2].split('\t').nth(6), Some("+"));
        assert_eq!(lines[4].split('\t').nth(6), Some("-"));
    }

    #[test]
    fn write_detailed_gtf_uses_exon_coords_from_original_transcript_index() {
        let transcripts = vec![tx(8, "AAAA", '+', 0.9), tx(3, "CCCCCC", '-', 0.8)];
        let coords = vec![vec![(10, 20)], vec![(1, 3), (5, 6)]];
        let output = NamedTempFile::new().unwrap();

        write_detailed_gtf(&transcripts, output.path().to_str().unwrap(), Some(&coords)).unwrap();
        let content = fs::read_to_string(output.path()).unwrap();
        let lines: Vec<&str> = content.lines().collect();

        // Sorted output starts with transcript_3, but should still use coords[1].
        assert!(lines[0].starts_with("transcript_3\t"));
        assert!(lines[1].contains("\texon\t1\t3\t"));
        assert!(lines[2].contains("\texon\t5\t6\t"));
        assert!(lines[3].starts_with("transcript_8\t"));
        assert!(lines[4].contains("\texon\t10\t20\t"));
    }
}
