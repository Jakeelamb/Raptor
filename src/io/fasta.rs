// src/io/fasta.rs
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Result, Write};
use std::path::Path;

use crate::graph::assembler::Contig;
use crate::kmer::rle::rle_encode;

pub enum FastaWriter {
    Plain(BufWriter<File>),
    Compressed(BufWriter<GzEncoder<File>>),
}

/// Open a FASTA file for reading, handles gzipped files automatically
pub fn open_fasta(path: &str) -> Box<dyn BufRead> {
    let file = File::open(path).expect("Unable to open FASTA file");
    if path.ends_with(".gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    }
}

impl FastaWriter {
    pub fn new(path: &str) -> Self {
        let path_ref = Path::new(path);
        if let Some(parent) = path_ref.parent().filter(|p| !p.as_os_str().is_empty()) {
            fs::create_dir_all(parent).expect("Unable to create FASTA parent directory");
        }
        let file = File::create(path_ref).expect("Unable to create FASTA file");
        if path.ends_with(".gz") {
            let encoder = GzEncoder::new(file, Compression::default());
            FastaWriter::Compressed(BufWriter::new(encoder))
        } else {
            FastaWriter::Plain(BufWriter::new(file))
        }
    }

    pub fn write_record(&mut self, header: &str, sequence: &str) -> Result<()> {
        match self {
            FastaWriter::Plain(writer) => {
                writeln!(writer, ">{}", header.trim_start_matches('@'))?;
                writeln!(writer, "{}", sequence)?;
            }
            FastaWriter::Compressed(writer) => {
                writeln!(writer, ">{}", header.trim_start_matches('@'))?;
                writeln!(writer, "{}", sequence)?;
            }
        };
        Ok(())
    }

    pub fn write_contig(&mut self, contig: &Contig, id: usize) -> Result<()> {
        match self {
            FastaWriter::Plain(writer) => {
                writeln!(writer, ">contig_{}", id)?;
                writeln!(writer, "{}", contig.sequence)?;
            }
            FastaWriter::Compressed(writer) => {
                writeln!(writer, ">contig_{}", id)?;
                writeln!(writer, "{}", contig.sequence)?;
            }
        };
        Ok(())
    }

    /// Write contig with RLE-encoded sequence in comment or custom format
    pub fn write_rle_contig(&mut self, contig: &Contig, id: usize) -> Result<()> {
        let rle = rle_encode(&contig.sequence);
        let encoded = rle
            .iter()
            .map(|(b, c)| format!("{}{}", *b as char, c))
            .collect::<Vec<_>>()
            .join("");

        match self {
            FastaWriter::Plain(writer) => {
                writeln!(writer, ">contig_{}_RLE", id)?;
                writeln!(writer, "{}", encoded)?;
            }
            FastaWriter::Compressed(writer) => {
                writeln!(writer, ">contig_{}_RLE", id)?;
                writeln!(writer, "{}", encoded)?;
            }
        };
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::FastaWriter;
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn test_fasta_writer_creates_parent_directories() {
        let temp_dir = TempDir::new().unwrap();
        let output = temp_dir.path().join("nested/output/test.fa");

        let mut writer = FastaWriter::new(output.to_str().unwrap());
        writer.write_record("seq1", "ACGT").unwrap();
        drop(writer);

        assert!(output.exists());
        assert!(fs::read_to_string(output).unwrap().contains(">seq1"));
    }
}
