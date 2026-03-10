// src/io/fastq.rs
#[derive(Debug, Clone)]
pub struct FastqRecord {
    pub header: String,
    pub sequence: String,
    pub plus: String,
    pub quality: String,
}

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

pub fn open_fastq(path: &str) -> Box<dyn BufRead> {
    let file = File::open(path).expect("Unable to open FASTQ file");
    if path.ends_with(".gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    }
}

/// DEPRECATED: Use stream_fastq_records() instead for memory efficiency.
/// This function loads the entire file twice (once for lines, once for records).
#[deprecated(note = "Use stream_fastq_records() for 50-66% memory reduction")]
pub fn read_fastq_records<R: BufRead>(reader: R) -> impl Iterator<Item = FastqRecord> {
    // Delegate to streaming implementation for backwards compatibility
    stream_fastq_records(reader)
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

/// Stream FASTQ records with strict structural validation.
///
/// Unlike `stream_fastq_records`, this parser returns `io::Result` items and
/// fails fast on:
/// - truncated/incomplete records
/// - headers not starting with '@'
/// - plus lines not starting with '+'
/// - sequence/quality length mismatch
pub fn stream_fastq_records_checked<R: BufRead>(
    reader: R,
) -> impl Iterator<Item = io::Result<FastqRecord>> {
    FastqCheckedStreamParser {
        lines: reader.lines(),
        record_index: 0,
        finished: false,
    }
}

/// Iterator adaptor to handle streaming FASTQ parsing
pub struct FastqStreamParser<I>
where
    I: Iterator<Item = io::Result<String>>,
{
    lines: I,
}

pub struct FastqCheckedStreamParser<I>
where
    I: Iterator<Item = io::Result<String>>,
{
    lines: I,
    record_index: usize,
    finished: bool,
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

impl<I> FastqCheckedStreamParser<I>
where
    I: Iterator<Item = io::Result<String>>,
{
    fn invalid_data(record_index: usize, message: &str) -> io::Error {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("FASTQ record {}: {}", record_index, message),
        )
    }

    fn unexpected_eof(record_index: usize, missing_field: &str) -> io::Error {
        io::Error::new(
            io::ErrorKind::UnexpectedEof,
            format!(
                "FASTQ record {}: truncated input, missing {} line",
                record_index, missing_field
            ),
        )
    }

    fn read_required_line(
        lines: &mut I,
        record_index: usize,
        field_name: &str,
    ) -> io::Result<String> {
        match lines.next() {
            Some(Ok(line)) => Ok(line),
            Some(Err(err)) => Err(err),
            None => Err(Self::unexpected_eof(record_index, field_name)),
        }
    }
}

impl<I> Iterator for FastqCheckedStreamParser<I>
where
    I: Iterator<Item = io::Result<String>>,
{
    type Item = io::Result<FastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        let header = match self.lines.next() {
            Some(Ok(line)) => line,
            Some(Err(err)) => {
                self.finished = true;
                return Some(Err(err));
            }
            None => return None,
        };

        self.record_index += 1;
        let record_index = self.record_index;

        let sequence = match Self::read_required_line(&mut self.lines, record_index, "sequence") {
            Ok(line) => line,
            Err(err) => {
                self.finished = true;
                return Some(Err(err));
            }
        };
        let plus = match Self::read_required_line(&mut self.lines, record_index, "plus") {
            Ok(line) => line,
            Err(err) => {
                self.finished = true;
                return Some(Err(err));
            }
        };
        let quality = match Self::read_required_line(&mut self.lines, record_index, "quality") {
            Ok(line) => line,
            Err(err) => {
                self.finished = true;
                return Some(Err(err));
            }
        };

        if !header.starts_with('@') {
            self.finished = true;
            return Some(Err(Self::invalid_data(
                record_index,
                "header must start with '@'",
            )));
        }
        if !plus.starts_with('+') {
            self.finished = true;
            return Some(Err(Self::invalid_data(
                record_index,
                "plus line must start with '+'",
            )));
        }
        if sequence.len() != quality.len() {
            self.finished = true;
            return Some(Err(Self::invalid_data(
                record_index,
                "sequence and quality lengths differ",
            )));
        }

        Some(Ok(FastqRecord {
            header,
            sequence,
            plus,
            quality,
        }))
    }
}

/// DEPRECATED: Use stream_paired_fastq_records() instead for memory efficiency.
/// This function loads both files entirely into memory twice.
#[deprecated(note = "Use stream_paired_fastq_records() for 50-66% memory reduction")]
pub fn read_paired_fastq_records<R1: BufRead, R2: BufRead>(
    reader1: R1,
    reader2: R2,
) -> impl Iterator<Item = (FastqRecord, FastqRecord)> {
    // Delegate to streaming implementation for backwards compatibility
    stream_paired_fastq_records(reader1, reader2)
}

/// Stream paired FASTQ records for memory-efficient processing
pub fn stream_paired_fastq_records<R1: BufRead, R2: BufRead>(
    reader1: R1,
    reader2: R2,
) -> impl Iterator<Item = (FastqRecord, FastqRecord)> {
    let r1_parser = FastqStreamParser {
        lines: reader1.lines(),
    };
    let r2_parser = FastqStreamParser {
        lines: reader2.lines(),
    };
    r1_parser.zip(r2_parser)
}

/// Stream paired FASTQ records with strict validation.
///
/// This parser validates both input streams with `stream_fastq_records_checked`
/// and additionally errors when the two files have different record counts.
pub fn stream_paired_fastq_records_checked<R1: BufRead, R2: BufRead>(
    reader1: R1,
    reader2: R2,
) -> impl Iterator<Item = io::Result<(FastqRecord, FastqRecord)>> {
    PairedFastqCheckedParser {
        r1: stream_fastq_records_checked(reader1),
        r2: stream_fastq_records_checked(reader2),
        pair_index: 0,
        finished: false,
    }
}

pub struct PairedFastqCheckedParser<I1, I2>
where
    I1: Iterator<Item = io::Result<FastqRecord>>,
    I2: Iterator<Item = io::Result<FastqRecord>>,
{
    r1: I1,
    r2: I2,
    pair_index: usize,
    finished: bool,
}

impl<I1, I2> Iterator for PairedFastqCheckedParser<I1, I2>
where
    I1: Iterator<Item = io::Result<FastqRecord>>,
    I2: Iterator<Item = io::Result<FastqRecord>>,
{
    type Item = io::Result<(FastqRecord, FastqRecord)>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        let left = self.r1.next();
        let right = self.r2.next();

        match (left, right) {
            (None, None) => None,
            (Some(Err(err)), _) => {
                self.finished = true;
                Some(Err(err))
            }
            (_, Some(Err(err))) => {
                self.finished = true;
                Some(Err(err))
            }
            (Some(Ok(r1)), Some(Ok(r2))) => {
                self.pair_index += 1;
                Some(Ok((r1, r2)))
            }
            (Some(Ok(_)), None) => {
                self.finished = true;
                Some(Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "Paired FASTQ inputs have different record counts: missing mate in R2 at pair {}",
                        self.pair_index + 1
                    ),
                )))
            }
            (None, Some(Ok(_))) => {
                self.finished = true;
                Some(Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "Paired FASTQ inputs have different record counts: missing mate in R1 at pair {}",
                        self.pair_index + 1
                    ),
                )))
            }
        }
    }
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
            }
            FastqWriter::Compressed(writer) => {
                writeln!(writer, "{}", record.header)?;
                writeln!(writer, "{}", record.sequence)?;
                writeln!(writer, "{}", record.plus)?;
                writeln!(writer, "{}", record.quality)?;
            }
        };
        Ok(())
    }
}

/// Read long reads from a FASTQ file
///
/// This function reads all records from a FASTQ file and returns them as a vector.
/// It's specifically intended for loading long reads for transcript polishing.
///
/// # Arguments
/// * `path` - Path to the FASTQ file (can be gzipped)
///
/// # Returns
/// * Result containing a vector of FastqRecord on success, or an io::Error on failure
pub fn read_long_reads(path: &str) -> io::Result<Vec<FastqRecord>> {
    let reader = open_fastq(path);
    let mut records = Vec::new();
    for record in stream_fastq_records_checked(reader) {
        records.push(record?);
    }

    if records.is_empty() {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "No reads found in the file",
        ))
    } else {
        Ok(records)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn checked_fastq_parser_reports_truncated_record() {
        let input = b"@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\n";
        let mut iter = stream_fastq_records_checked(Cursor::new(&input[..]));

        assert!(iter.next().expect("first record missing").is_ok());
        let err = iter
            .next()
            .expect("expected parse error for truncated record")
            .expect_err("expected truncated FASTQ error");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
        assert!(err.to_string().contains("record 2"));
    }

    #[test]
    fn checked_fastq_parser_rejects_length_mismatch() {
        let input = b"@r1\nACGT\n+\nIII\n";
        let mut iter = stream_fastq_records_checked(Cursor::new(&input[..]));

        let err = iter
            .next()
            .expect("expected parse error")
            .expect_err("expected length mismatch error");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("lengths differ"));
    }

    #[test]
    fn checked_paired_fastq_parser_rejects_mismatched_record_counts() {
        let r1 = b"@r1/1\nAAAA\n+\nIIII\n@r2/1\nCCCC\n+\nIIII\n";
        let r2 = b"@r1/2\nTTTT\n+\nIIII\n";
        let mut iter =
            stream_paired_fastq_records_checked(Cursor::new(&r1[..]), Cursor::new(&r2[..]));

        assert!(iter.next().expect("first pair missing").is_ok());
        let err = iter
            .next()
            .expect("expected mismatch error")
            .expect_err("expected paired mismatch error");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("different record counts"));
    }

    #[test]
    fn read_long_reads_rejects_malformed_fastq() {
        let mut file = NamedTempFile::new().expect("temp fastq");
        writeln!(file, "@r1").expect("header");
        writeln!(file, "ACGT").expect("sequence");
        writeln!(file, "+").expect("plus");
        file.flush().expect("flush");

        let err = read_long_reads(file.path().to_str().expect("path"))
            .expect_err("expected strict parser error");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn read_long_reads_accepts_well_formed_fastq() {
        let mut file = NamedTempFile::new().expect("temp fastq");
        writeln!(file, "@r1").expect("header");
        writeln!(file, "ACGT").expect("sequence");
        writeln!(file, "+").expect("plus");
        writeln!(file, "IIII").expect("quality");
        file.flush().expect("flush");

        let records = read_long_reads(file.path().to_str().expect("path")).expect("valid fastq");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].sequence, "ACGT");
    }
}
