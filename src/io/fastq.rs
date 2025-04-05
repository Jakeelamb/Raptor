// src/io/fastq.rs
#[derive(Debug, Clone)]
pub struct FastqRecord {
    pub header: String,
    pub sequence: String,
    pub plus: String,
    pub quality: String,
}

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

pub fn open_fastq(path: &str) -> Box<dyn BufRead> {
    let file = File::open(path).expect("Unable to open FASTQ file");
    if path.ends_with(".gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    }
}

pub fn read_fastq_records<R: BufRead>(reader: R) -> impl Iterator<Item = FastqRecord> {
    reader
        .lines()
        .collect::<std::io::Result<Vec<_>>>()
        .expect("Failed to read lines")
        .chunks(4)
        .map(|chunk| FastqRecord {
            header: chunk[0].clone(),
            sequence: chunk[1].clone(),
            plus: chunk[2].clone(),
            quality: chunk[3].clone(),
        })
        .collect::<Vec<_>>()
        .into_iter()
}

/// Stream FASTQ records for memory-efficient processing
/// 
/// This function processes FASTQ records in a streaming fashion rather than
/// loading the entire file into memory first. This allows handling of very large
/// FASTQ files with bounded memory usage.
pub fn stream_fastq_records<R: BufRead>(reader: R) -> impl Iterator<Item = FastqRecord> {
    let lines = reader.lines();
    FastqStreamParser { lines }
}

/// Iterator adaptor to handle streaming FASTQ parsing
pub struct FastqStreamParser<I>
where
    I: Iterator<Item = io::Result<String>>,
{
    lines: I,
}

impl<I> Iterator for FastqStreamParser<I>
where
    I: Iterator<Item = io::Result<String>>,
{
    type Item = FastqRecord;

    fn next(&mut self) -> Option<Self::Item> {
        let header = match self.lines.next() {
            Some(Ok(line)) => line,
            _ => return None,
        };
        
        let sequence = match self.lines.next() {
            Some(Ok(line)) => line,
            _ => return None,
        };
        
        let plus = match self.lines.next() {
            Some(Ok(line)) => line,
            _ => return None,
        };
        
        let quality = match self.lines.next() {
            Some(Ok(line)) => line,
            _ => return None,
        };
        
        Some(FastqRecord {
            header,
            sequence,
            plus,
            quality,
        })
    }
}

pub fn read_paired_fastq_records<R1: BufRead, R2: BufRead>(
    reader1: R1,
    reader2: R2,
) -> impl Iterator<Item = (FastqRecord, FastqRecord)> {
    let r1_lines = reader1.lines().collect::<Result<Vec<_>, _>>().unwrap();
    let r2_lines = reader2.lines().collect::<Result<Vec<_>, _>>().unwrap();

    r1_lines
        .chunks(4)
        .zip(r2_lines.chunks(4))
        .map(|(r1, r2)| {
            (
                FastqRecord {
                    header: r1[0].clone(),
                    sequence: r1[1].clone(),
                    plus: r1[2].clone(),
                    quality: r1[3].clone(),
                },
                FastqRecord {
                    header: r2[0].clone(),
                    sequence: r2[1].clone(),
                    plus: r2[2].clone(),
                    quality: r2[3].clone(),
                },
            )
        })
        .collect::<Vec<_>>()
        .into_iter()
}

/// Stream paired FASTQ records for memory-efficient processing
pub fn stream_paired_fastq_records<R1: BufRead, R2: BufRead>(
    reader1: R1,
    reader2: R2,
) -> impl Iterator<Item = (FastqRecord, FastqRecord)> {
    let r1_parser = FastqStreamParser { lines: reader1.lines() };
    let r2_parser = FastqStreamParser { lines: reader2.lines() };
    r1_parser.zip(r2_parser)
}

pub enum FastqWriter {
    Plain(BufWriter<File>),
    Compressed(BufWriter<GzEncoder<File>>),
}

impl FastqWriter {
    pub fn new(path: &str) -> Self {
        let file = File::create(path).expect("Unable to create output FASTQ file");
        if path.ends_with(".gz") {
            let encoder = GzEncoder::new(file, Compression::default());
            FastqWriter::Compressed(BufWriter::new(encoder))
        } else {
            FastqWriter::Plain(BufWriter::new(file))
        }
    }

    pub fn write_record(&mut self, record: &FastqRecord) -> io::Result<()> {
        match self {
            FastqWriter::Plain(writer) => {
                writeln!(writer, "{}", record.header)?;
                writeln!(writer, "{}", record.sequence)?;
                writeln!(writer, "{}", record.plus)?;
                writeln!(writer, "{}", record.quality)?;
            },
            FastqWriter::Compressed(writer) => {
                writeln!(writer, "{}", record.header)?;
                writeln!(writer, "{}", record.sequence)?;
                writeln!(writer, "{}", record.plus)?;
                writeln!(writer, "{}", record.quality)?;
            },
        };
        Ok(())
    }
}
