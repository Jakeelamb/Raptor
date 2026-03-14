#![doc = include_str!("../README.md")]

use flate2::read::MultiGzDecoder;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MaskMode {
    N,
    X,
    Lowercase,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RepeatMaskerAnnotation {
    pub score: usize,
    pub percent_divergence_tenths: u16,
    pub percent_deletions_tenths: u16,
    pub percent_insertions_tenths: u16,
    pub query_name: String,
    pub query_begin: usize,
    pub query_end: usize,
    pub query_left: usize,
    pub orientation: char,
    pub subject_name: String,
    pub subject_class: String,
    pub subject_begin: usize,
    pub subject_end: usize,
    pub subject_left: usize,
    pub annotation_id: String,
    pub lineage_id: String,
    pub overlap_marker: Option<String>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct RepeatMaskerAlignmentMetrics {
    pub aligned_query_bases: usize,
    pub aligned_subject_bases: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub matrix_name: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub kimura_divergence: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub raw_kimura_divergence: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cpg_sites: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub transition_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub transversion_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gap_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub average_gap_size: Option<f64>,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct RepeatMaskerAlignmentRecord {
    pub score: usize,
    pub percent_divergence_tenths: u16,
    pub percent_deletions_tenths: u16,
    pub percent_insertions_tenths: u16,
    pub query_name: String,
    pub query_begin: usize,
    pub query_end: usize,
    pub query_left: usize,
    pub orientation: char,
    pub subject_name: String,
    pub subject_begin: usize,
    pub subject_end: usize,
    pub subject_left: usize,
    pub annotation_id: String,
    pub lineage_id: String,
    pub query_sequence: String,
    pub subject_sequence: String,
    pub metrics: RepeatMaskerAlignmentMetrics,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct RepeatMaskerAlignmentBlock {
    pub query_begin: usize,
    pub query_end: usize,
    pub subject_begin: usize,
    pub subject_end: usize,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct RepeatMaskerAlignmentEncoding {
    pub cigar: String,
    pub caf: String,
    pub match_bases: usize,
    pub mismatch_bases: usize,
    pub insertion_bases: usize,
    pub deletion_bases: usize,
    pub blocks: Vec<RepeatMaskerAlignmentBlock>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MaskInterval {
    pub begin: usize,
    pub end: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FastaRecord {
    pub id: String,
    pub description: String,
    pub sequence: Vec<u8>,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct MaskStats {
    pub sequence_count: usize,
    pub total_bases: usize,
    pub non_ambiguous_bases_excluding_long_runs: usize,
    pub gc_fraction: f64,
    pub masked_bases: usize,
    pub annotation_count: usize,
    pub interval_count: usize,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct AnnotationStats {
    pub annotation_count: usize,
    pub total_aligned_bases: usize,
    pub occupied_bases: usize,
    pub sequence_count: usize,
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize)]
pub struct RepeatMaskerMetadata {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub total_sequences: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub total_length: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub total_non_mask_bases: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub total_non_sub_bases: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub version: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub engine: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub library: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reported_masked_bases: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reported_gc_percent_hundredths: Option<u16>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reported_masked_gc_percent_hundredths: Option<u16>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub batch_overlap_boundaries: BTreeMap<String, Vec<usize>>,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct RepeatMaskerCatalog {
    pub annotations: Vec<RepeatMaskerAnnotation>,
    pub metadata: RepeatMaskerMetadata,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct RepeatSummaryRow {
    pub key: String,
    pub annotation_count: usize,
    pub aligned_bases: usize,
    pub occupied_bases: usize,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct RepeatMaskerAnnotationChain {
    pub chain_id: usize,
    pub query_name: String,
    pub orientation: char,
    pub subject_name: String,
    pub subject_class: String,
    pub query_begin: usize,
    pub query_end: usize,
    pub subject_begin: usize,
    pub subject_end: usize,
    pub annotation_count: usize,
    pub annotation_ids: Vec<String>,
    pub lineage_ids: Vec<String>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct RepeatMaskerCatalogStats {
    #[serde(flatten)]
    pub annotation_stats: AnnotationStats,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub raw_annotation_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub suppressed_annotation_count: Option<usize>,
    pub element_chain_count: usize,
    pub repeat_class_summary: Vec<RepeatSummaryRow>,
    pub repeat_family_summary: Vec<RepeatSummaryRow>,
    #[serde(default)]
    pub metadata: RepeatMaskerMetadata,
}

type AlignmentLinkKey = (String, String, String, String);

#[derive(Debug)]
pub enum RepeatMaskerError {
    Io(io::Error),
    InvalidAnnotationLine {
        line_number: usize,
        line: String,
        message: String,
    },
    InvalidUtf8Path(String),
}

impl fmt::Display for RepeatMaskerError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(err) => write!(f, "{err}"),
            Self::InvalidAnnotationLine {
                line_number,
                message,
                ..
            } => write!(
                f,
                "invalid RepeatMasker .out annotation at line {line_number}: {message}"
            ),
            Self::InvalidUtf8Path(path) => write!(f, "path is not valid UTF-8: {path}"),
        }
    }
}

impl std::error::Error for RepeatMaskerError {}

impl From<io::Error> for RepeatMaskerError {
    fn from(value: io::Error) -> Self {
        Self::Io(value)
    }
}

fn path_to_str(path: &Path) -> Result<&str, RepeatMaskerError> {
    path.to_str()
        .ok_or_else(|| RepeatMaskerError::InvalidUtf8Path(path.display().to_string()))
}

fn open_maybe_gzip(path: &Path) -> Result<Box<dyn BufRead>, RepeatMaskerError> {
    let file = File::open(path)?;
    if path_to_str(path)?.ends_with(".gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn parse_parenthesized_usize(token: &str) -> Result<usize, String> {
    token
        .trim_matches(|c| c == '(' || c == ')')
        .parse::<usize>()
        .map_err(|err| format!("invalid parenthesized integer '{token}': {err}"))
}

fn is_parenthesized_integer(token: &str) -> bool {
    token.starts_with('(')
        && token.ends_with(')')
        && token
            .trim_matches(|c| c == '(' || c == ')')
            .parse::<isize>()
            .is_ok()
}

fn is_percent_token(token: &str) -> bool {
    token.parse::<f64>().is_ok()
}

fn is_positive_integer_token(token: &str) -> bool {
    token.parse::<usize>().is_ok()
}

fn looks_like_annotation_prefix(fields: &[&str]) -> bool {
    fields.len() >= 8
        && is_positive_integer_token(fields[0])
        && is_percent_token(fields[1])
        && is_percent_token(fields[2])
        && is_percent_token(fields[3])
        && is_positive_integer_token(fields[5])
        && is_positive_integer_token(fields[6])
        && is_parenthesized_integer(fields[7])
}

fn looks_like_out_annotation_fields(fields: &[&str]) -> bool {
    looks_like_annotation_prefix(fields)
        && fields.len() >= 15
        && match fields[8] {
            "+" => is_parenthesized_integer(fields[13]),
            "C" | "c" => is_parenthesized_integer(fields[11]),
            _ => false,
        }
}

fn looks_like_cat_annotation_fields(fields: &[&str]) -> bool {
    looks_like_annotation_prefix(fields)
        && fields.len() >= 12
        && match fields[8] {
            "+" => false,
            "C" | "c" => is_parenthesized_integer(fields[10]),
            _ => is_parenthesized_integer(fields[11]),
        }
}

fn looks_like_align_annotation_fields(fields: &[&str]) -> bool {
    looks_like_annotation_prefix(fields)
        && fields.len() >= 13
        && !looks_like_out_annotation_fields(fields)
        && match fields[8] {
            "+" => false,
            "C" | "c" => is_parenthesized_integer(fields[10]),
            _ => is_parenthesized_integer(fields[11]),
        }
}

fn split_subject_token(token: &str) -> (String, String) {
    match token.rsplit_once('#') {
        Some((name, class_name)) => (name.to_string(), class_name.to_string()),
        None => (token.to_string(), String::new()),
    }
}

fn parse_batch_overlap_boundaries(line: &str) -> Option<(String, Vec<usize>)> {
    let trimmed = line.strip_prefix("##").map(str::trim_start).unwrap_or(line);
    let (name, values) = trimmed.split_once(char::is_whitespace)?;
    let values = values.trim();
    if name.is_empty() || !values.contains(',') {
        return None;
    }
    let boundaries: Option<Vec<usize>> = values
        .split(',')
        .map(|value| value.trim().parse::<usize>().ok())
        .collect();
    let boundaries = boundaries?;
    if boundaries.is_empty() {
        return None;
    }
    Some((name.to_string(), boundaries))
}

fn update_metadata_from_line(line: &str, metadata: &mut RepeatMaskerMetadata) {
    let trimmed = line.trim();
    if trimmed.is_empty() {
        return;
    }

    if let Some((name, boundaries)) = parse_batch_overlap_boundaries(trimmed) {
        metadata.batch_overlap_boundaries.insert(name, boundaries);
        return;
    }

    let normalized = trimmed
        .strip_prefix("##")
        .map(str::trim_start)
        .unwrap_or(trimmed);

    if let Some(value) = normalized.strip_prefix("Total Sequences:") {
        metadata.total_sequences = value.trim().parse::<usize>().ok();
        return;
    }
    if let Some(value) = normalized.strip_prefix("Total Length:") {
        metadata.total_length = value.trim().parse::<usize>().ok();
        return;
    }
    if let Some(value) = normalized.strip_prefix("Total NonMask") {
        if let Some((_, count)) = value.split_once(':') {
            metadata.total_non_mask_bases = count.trim().parse::<usize>().ok();
        }
        return;
    }
    if let Some(value) = normalized.strip_prefix("Total NonSub") {
        if let Some((_, count)) = value.split_once(':') {
            metadata.total_non_sub_bases = count.trim().parse::<usize>().ok();
        }
        return;
    }
    if normalized.starts_with("RepeatMasker version ") {
        metadata.version = Some(normalized.to_string());
        return;
    }
    if normalized.starts_with("run with ") {
        metadata.engine = Some(normalized.to_string());
        return;
    }
    if let Some(value) = normalized.strip_prefix("RM Library:") {
        metadata.library = Some(value.trim().to_string());
        return;
    }
    if normalized.starts_with("RepBase ") {
        metadata.library = Some(normalized.to_string());
        return;
    }

    let fields: Vec<&str> = normalized.split_whitespace().collect();
    if fields.len() >= 5
        && fields[1].starts_with("sequence")
        && fields[2] == "with"
        && fields[4] == "bp"
    {
        metadata.total_sequences = fields[0].parse::<usize>().ok();
        metadata.total_length = fields[3].parse::<usize>().ok();
        if let Some(start) = normalized.find('(') {
            if let Some(end) = normalized[start + 1..].find(" bp") {
                metadata.total_non_mask_bases = normalized[start + 1..start + 1 + end]
                    .trim()
                    .parse::<usize>()
                    .ok();
            }
        }
        return;
    }
    if fields.len() == 3 && fields[1] == "bp" && fields[2] == "masked" {
        metadata.reported_masked_bases = fields[0].parse::<usize>().ok();
        return;
    }
    if normalized.starts_with("GC level") {
        let values: Vec<u16> = normalized
            .split_whitespace()
            .map(|token| token.trim_matches(|c| c == '(' || c == ')' || c == '%'))
            .filter_map(|token| parse_percent_hundredths(token).ok())
            .collect();
        metadata.reported_gc_percent_hundredths = values.first().copied();
        metadata.reported_masked_gc_percent_hundredths = values.get(1).copied();
    }
}

fn format_percent_tenths(value: u16) -> String {
    format!("{}.{}", value / 10, value % 10)
}

fn format_percent_hundredths(value: u16) -> String {
    format!("{}.{:02}", value / 100, value % 100)
}

fn gff_escape(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for byte in value.bytes() {
        let safe = byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b':' | b'_' | b'-');
        if safe {
            escaped.push(byte as char);
        } else {
            use std::fmt::Write as _;
            let _ = write!(escaped, "%{byte:02X}");
        }
    }
    escaped
}

fn parse_percent_tenths(token: &str) -> Result<u16, String> {
    let trimmed = token.trim();
    if trimmed.is_empty() {
        return Err("empty percentage field".to_string());
    }
    let negative = trimmed.starts_with('-');
    if negative {
        return Err(format!("negative percentage is not supported: '{token}'"));
    }

    let (whole, frac) = match trimmed.split_once('.') {
        Some((whole, frac)) => (whole, frac),
        None => (trimmed, "0"),
    };

    let whole_value = whole
        .parse::<u16>()
        .map_err(|err| format!("invalid percentage '{token}': {err}"))?;
    let frac_digit = frac
        .chars()
        .next()
        .unwrap_or('0')
        .to_digit(10)
        .ok_or_else(|| format!("invalid percentage '{token}'"))? as u16;
    Ok(whole_value.saturating_mul(10).saturating_add(frac_digit))
}

fn parse_percent_hundredths(token: &str) -> Result<u16, String> {
    let trimmed = token.trim().trim_end_matches('%');
    if trimmed.is_empty() {
        return Err("empty percentage field".to_string());
    }
    let negative = trimmed.starts_with('-');
    if negative {
        return Err(format!("negative percentage is not supported: '{token}'"));
    }

    let (whole, frac) = match trimmed.split_once('.') {
        Some((whole, frac)) => (whole, frac),
        None => (trimmed, "0"),
    };

    let whole_value = whole
        .parse::<u16>()
        .map_err(|err| format!("invalid percentage '{token}': {err}"))?;
    let mut digits = frac.chars();
    let first = digits.next().unwrap_or('0');
    let second = digits.next().unwrap_or('0');
    let first = first
        .to_digit(10)
        .ok_or_else(|| format!("invalid percentage '{token}'"))? as u16;
    let second = second
        .to_digit(10)
        .ok_or_else(|| format!("invalid percentage '{token}'"))? as u16;
    Ok(whole_value
        .saturating_mul(100)
        .saturating_add(first.saturating_mul(10))
        .saturating_add(second))
}

fn alignment_link_key(
    query_name: &str,
    subject_name: &str,
    lineage_id: &str,
    annotation_id: &str,
) -> AlignmentLinkKey {
    (
        query_name.to_string(),
        subject_name.to_string(),
        lineage_id.to_string(),
        annotation_id.to_string(),
    )
}

fn annotation_alignment_key(annotation: &RepeatMaskerAnnotation) -> AlignmentLinkKey {
    alignment_link_key(
        annotation.query_name.as_str(),
        annotation.subject_name.as_str(),
        annotation.lineage_id.as_str(),
        annotation.annotation_id.as_str(),
    )
}

fn alignment_record_key(record: &RepeatMaskerAlignmentRecord) -> AlignmentLinkKey {
    alignment_link_key(
        record.query_name.as_str(),
        record.subject_name.as_str(),
        record.lineage_id.as_str(),
        record.annotation_id.as_str(),
    )
}

fn reverse_complement_iupac(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'Y' => b'R',
            b'R' => b'Y',
            b'M' => b'K',
            b'K' => b'M',
            b'H' => b'D',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            b'y' => b'r',
            b'r' => b'y',
            b'm' => b'k',
            b'k' => b'm',
            b'h' => b'd',
            b'b' => b'v',
            b'v' => b'b',
            b'd' => b'h',
            other => other,
        })
        .map(char::from)
        .collect()
}

fn parse_alignment_metadata_f64(
    line: &str,
    prefix: &str,
) -> Result<Option<f64>, RepeatMaskerError> {
    let Some(value) = line.trim().strip_prefix(prefix) else {
        return Ok(None);
    };
    let value = value.trim();
    if value.eq_ignore_ascii_case("Unknown") {
        return Ok(None);
    }
    value
        .parse::<f64>()
        .map(Some)
        .map_err(|err| RepeatMaskerError::InvalidAnnotationLine {
            line_number: 0,
            line: line.to_string(),
            message: format!("invalid floating-point alignment metadata '{value}': {err}"),
        })
}

fn count_non_gap_bases(sequence: &str) -> usize {
    sequence.bytes().filter(|base| *base != b'-').count()
}

fn parse_align_annotation_line(
    line: &str,
    line_number: usize,
) -> Result<Option<RepeatMaskerAlignmentRecord>, RepeatMaskerError> {
    let trimmed = line.trim();
    if trimmed.is_empty() {
        return Ok(None);
    }
    let first = trimmed.as_bytes()[0];
    if !first.is_ascii_digit() {
        return Ok(None);
    }

    let fields: Vec<&str> = trimmed.split_whitespace().collect();
    if !looks_like_align_annotation_fields(&fields) {
        return Ok(None);
    }

    let score =
        fields[0]
            .parse::<usize>()
            .map_err(|err| RepeatMaskerError::InvalidAnnotationLine {
                line_number,
                line: line.to_string(),
                message: format!("invalid score '{}': {err}", fields[0]),
            })?;
    let percent_divergence_tenths = parse_percent_tenths(fields[1]).map_err(|message| {
        RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message,
        }
    })?;
    let percent_deletions_tenths = parse_percent_tenths(fields[2]).map_err(|message| {
        RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message,
        }
    })?;
    let percent_insertions_tenths = parse_percent_tenths(fields[3]).map_err(|message| {
        RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message,
        }
    })?;
    let query_begin =
        fields[5]
            .parse::<usize>()
            .map_err(|err| RepeatMaskerError::InvalidAnnotationLine {
                line_number,
                line: line.to_string(),
                message: format!("invalid query begin '{}': {err}", fields[5]),
            })?;
    let query_end =
        fields[6]
            .parse::<usize>()
            .map_err(|err| RepeatMaskerError::InvalidAnnotationLine {
                line_number,
                line: line.to_string(),
                message: format!("invalid query end '{}': {err}", fields[6]),
            })?;
    if query_begin == 0 || query_end < query_begin {
        return Err(RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message: format!("invalid query interval {query_begin}-{query_end}"),
        });
    }
    let query_left = parse_parenthesized_usize(fields[7]).map_err(|message| {
        RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message,
        }
    })?;

    let (orientation, subject_name, subject_coord1, subject_coord2, subject_coord3, tail_index) =
        if fields[8].eq_ignore_ascii_case("C") {
            ('C', fields[9], fields[10], fields[11], fields[12], 13usize)
        } else {
            ('+', fields[8], fields[9], fields[10], fields[11], 12usize)
        };

    let (subject_begin, subject_end, subject_left) = if orientation == 'C' {
        (
            subject_coord3.parse::<usize>().map_err(|err| {
                RepeatMaskerError::InvalidAnnotationLine {
                    line_number,
                    line: line.to_string(),
                    message: format!("invalid reverse subject begin '{subject_coord3}': {err}"),
                }
            })?,
            subject_coord2.parse::<usize>().map_err(|err| {
                RepeatMaskerError::InvalidAnnotationLine {
                    line_number,
                    line: line.to_string(),
                    message: format!("invalid reverse subject end '{subject_coord2}': {err}"),
                }
            })?,
            parse_parenthesized_usize(subject_coord1).map_err(|message| {
                RepeatMaskerError::InvalidAnnotationLine {
                    line_number,
                    line: line.to_string(),
                    message,
                }
            })?,
        )
    } else {
        (
            subject_coord1.parse::<usize>().map_err(|err| {
                RepeatMaskerError::InvalidAnnotationLine {
                    line_number,
                    line: line.to_string(),
                    message: format!("invalid subject begin '{subject_coord1}': {err}"),
                }
            })?,
            subject_coord2.parse::<usize>().map_err(|err| {
                RepeatMaskerError::InvalidAnnotationLine {
                    line_number,
                    line: line.to_string(),
                    message: format!("invalid subject end '{subject_coord2}': {err}"),
                }
            })?,
            parse_parenthesized_usize(subject_coord3).map_err(|message| {
                RepeatMaskerError::InvalidAnnotationLine {
                    line_number,
                    line: line.to_string(),
                    message,
                }
            })?,
        )
    };
    if subject_begin == 0 || subject_end < subject_begin {
        return Err(RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message: format!("invalid subject interval {subject_begin}-{subject_end}"),
        });
    }

    let tail = &fields[tail_index..];
    let (lineage_id, annotation_id) = match tail {
        [lineage, annotation, ..] => ((*lineage).to_string(), (*annotation).to_string()),
        [lineage] => ((*lineage).to_string(), (*lineage).to_string()),
        [] => (String::new(), format!("align-{line_number}")),
    };

    Ok(Some(RepeatMaskerAlignmentRecord {
        score,
        percent_divergence_tenths,
        percent_deletions_tenths,
        percent_insertions_tenths,
        query_name: fields[4].to_string(),
        query_begin,
        query_end,
        query_left,
        orientation,
        subject_name: subject_name.to_string(),
        subject_begin,
        subject_end,
        subject_left,
        annotation_id,
        lineage_id,
        query_sequence: String::new(),
        subject_sequence: String::new(),
        metrics: RepeatMaskerAlignmentMetrics::default(),
    }))
}

fn finalize_alignment_record(
    current: &mut Option<RepeatMaskerAlignmentRecord>,
    query_complemented: &mut bool,
    records: &mut Vec<RepeatMaskerAlignmentRecord>,
) {
    let Some(mut record) = current.take() else {
        return;
    };
    if *query_complemented {
        record.query_sequence = reverse_complement_iupac(&record.query_sequence);
        record.subject_sequence = reverse_complement_iupac(&record.subject_sequence);
    }
    record.metrics.aligned_query_bases = count_non_gap_bases(&record.query_sequence);
    record.metrics.aligned_subject_bases = count_non_gap_bases(&record.subject_sequence);
    records.push(record);
    *query_complemented = false;
}

pub fn load_repeatmasker_alignments(
    reader: impl BufRead,
) -> Result<Vec<RepeatMaskerAlignmentRecord>, RepeatMaskerError> {
    let mut records = Vec::new();
    let mut current: Option<RepeatMaskerAlignmentRecord> = None;
    let mut query_complemented = false;
    let mut expect_query_line = true;

    for (index, line) in reader.lines().enumerate() {
        let line = line?;
        let line_number = index + 1;

        if let Some(record) = parse_align_annotation_line(&line, line_number)? {
            finalize_alignment_record(&mut current, &mut query_complemented, &mut records);
            current = Some(record);
            expect_query_line = true;
            continue;
        }

        let Some(active) = current.as_mut() else {
            continue;
        };

        if let Some(matrix_name) = line.trim().strip_prefix("Matrix = ") {
            if !matrix_name.eq_ignore_ascii_case("Unknown") {
                active.metrics.matrix_name = Some(matrix_name.trim().to_string());
            }
            continue;
        }
        if let Some(kimura) = parse_alignment_metadata_f64(&line, "Kimura (with divCpGMod) = ")? {
            active.metrics.kimura_divergence = Some(kimura);
            continue;
        }
        if let Some(values) = line.trim().strip_prefix("CpG sites = ") {
            if let Some((cpg_sites, raw_kimura)) = values.split_once(", Kimura (unadjusted) = ") {
                active.metrics.cpg_sites = cpg_sites.trim().parse::<usize>().ok();
                active.metrics.raw_kimura_divergence = raw_kimura.trim().parse::<f64>().ok();
            }
            continue;
        }
        if let Some(values) = line.trim().strip_prefix("Transitions / transversions = ") {
            if let Some((_, counts)) = values.split_once('(') {
                if let Some(counts) = counts.strip_suffix(')') {
                    let mut fields = counts.split('/');
                    active.metrics.transition_count = fields
                        .next()
                        .and_then(|value| value.trim().parse::<usize>().ok());
                    active.metrics.transversion_count = fields
                        .next()
                        .and_then(|value| value.trim().parse::<usize>().ok());
                }
            }
            continue;
        }
        if let Some(values) = line.trim().strip_prefix("Gap_init rate = ") {
            let mut fields = values.split(", avg. gap size = ");
            if let Some(left) = fields.next() {
                if let Some((_, gap_counts)) = left.split_once('(') {
                    if let Some(gap_counts) = gap_counts.strip_suffix(')') {
                        active.metrics.gap_count = gap_counts
                            .split('/')
                            .next()
                            .and_then(|value| value.trim().parse::<usize>().ok());
                    }
                }
            }
            if let Some(right) = fields.next() {
                active.metrics.average_gap_size = right
                    .split_whitespace()
                    .next()
                    .and_then(|value| value.parse::<f64>().ok());
            }
            finalize_alignment_record(&mut current, &mut query_complemented, &mut records);
            expect_query_line = true;
            continue;
        }

        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let mut sequence_fields = trimmed.split_whitespace();
        let first = sequence_fields.next();
        let second = sequence_fields.next();
        let third = sequence_fields.next();
        let fourth = sequence_fields.next();
        if let (Some(name_or_c), Some(name_or_start), Some(seq_or_start), Some(seq_or_end)) =
            (first, second, third, fourth)
        {
            let (is_complemented, _name, start, sequence, end) = if name_or_c == "C" {
                let Some(end) = sequence_fields.next() else {
                    continue;
                };
                (true, name_or_start, seq_or_start, seq_or_end, end)
            } else {
                (false, name_or_c, name_or_start, seq_or_start, seq_or_end)
            };
            if start.parse::<usize>().is_ok()
                && end.parse::<usize>().is_ok()
                && sequence.bytes().all(|byte| !byte.is_ascii_whitespace())
            {
                if expect_query_line {
                    if is_complemented {
                        query_complemented = true;
                    }
                    active.query_sequence.push_str(sequence);
                } else {
                    active.subject_sequence.push_str(sequence);
                }
                expect_query_line = !expect_query_line;
            }
        }
    }

    finalize_alignment_record(&mut current, &mut query_complemented, &mut records);
    Ok(records)
}

pub fn load_repeatmasker_alignments_path(
    path: &Path,
) -> Result<Vec<RepeatMaskerAlignmentRecord>, RepeatMaskerError> {
    load_repeatmasker_alignments(open_maybe_gzip(path)?)
}

pub fn encode_repeatmasker_alignment(
    record: &RepeatMaskerAlignmentRecord,
) -> RepeatMaskerAlignmentEncoding {
    let mut cigar = String::new();
    let mut prev_op = '\0';
    let mut op_len = 0usize;

    let mut caf = String::new();
    let mut in_deletion = false;
    let mut in_insertion = false;

    let mut match_bases = 0usize;
    let mut mismatch_bases = 0usize;
    let mut insertion_bases = 0usize;
    let mut deletion_bases = 0usize;

    let mut blocks = Vec::new();
    let mut current_block_query_begin = 0usize;
    let mut current_block_subject_begin = 0usize;
    let mut current_block_len = 0usize;

    let mut query_pos = record.query_begin;
    let mut subject_pos = if record.orientation == 'C' {
        record.subject_end
    } else {
        record.subject_begin
    };

    let flush_cigar = |cigar: &mut String, op: char, len: usize| {
        if len > 0 {
            use std::fmt::Write as _;
            let _ = write!(cigar, "{len}{op}");
        }
    };

    let flush_block = |blocks: &mut Vec<RepeatMaskerAlignmentBlock>,
                       query_begin: usize,
                       subject_begin: usize,
                       len: usize,
                       orientation: char| {
        if len == 0 {
            return;
        }
        let query_end = query_begin + len - 1;
        let (subject_begin, subject_end) = if orientation == 'C' {
            (subject_begin - len + 1, subject_begin)
        } else {
            (subject_begin, subject_begin + len - 1)
        };
        blocks.push(RepeatMaskerAlignmentBlock {
            query_begin,
            query_end,
            subject_begin,
            subject_end,
        });
    };

    for (query_char, subject_char) in record
        .query_sequence
        .bytes()
        .zip(record.subject_sequence.bytes())
    {
        let op = if query_char != b'-' && subject_char != b'-' {
            'M'
        } else if query_char == b'-' {
            'I'
        } else {
            'D'
        };
        if prev_op == op {
            op_len += 1;
        } else {
            flush_cigar(&mut cigar, prev_op, op_len);
            prev_op = op;
            op_len = 1;
        }

        match op {
            'M' => {
                if in_deletion {
                    caf.push('-');
                    in_deletion = false;
                } else if in_insertion {
                    caf.push('+');
                    in_insertion = false;
                }
                if query_char.eq_ignore_ascii_case(&subject_char) {
                    caf.push(query_char as char);
                    match_bases += 1;
                } else {
                    caf.push(query_char as char);
                    caf.push('/');
                    caf.push(subject_char as char);
                    mismatch_bases += 1;
                }

                if current_block_len == 0 {
                    current_block_query_begin = query_pos;
                    current_block_subject_begin = subject_pos;
                }
                current_block_len += 1;
            }
            'I' => {
                if current_block_len > 0 {
                    flush_block(
                        &mut blocks,
                        current_block_query_begin,
                        current_block_subject_begin,
                        current_block_len,
                        record.orientation,
                    );
                    current_block_len = 0;
                }
                if in_deletion {
                    caf.push('-');
                    in_deletion = false;
                }
                if !in_insertion {
                    caf.push('+');
                    in_insertion = true;
                }
                caf.push(subject_char as char);
                insertion_bases += 1;
            }
            'D' => {
                if current_block_len > 0 {
                    flush_block(
                        &mut blocks,
                        current_block_query_begin,
                        current_block_subject_begin,
                        current_block_len,
                        record.orientation,
                    );
                    current_block_len = 0;
                }
                if in_insertion {
                    caf.push('+');
                    in_insertion = false;
                }
                if !in_deletion {
                    caf.push('-');
                    in_deletion = true;
                }
                caf.push(query_char as char);
                deletion_bases += 1;
            }
            _ => unreachable!(),
        }

        if query_char != b'-' {
            query_pos += 1;
        }
        if subject_char != b'-' {
            if record.orientation == 'C' {
                subject_pos = subject_pos.saturating_sub(1);
            } else {
                subject_pos += 1;
            }
        }
    }

    if current_block_len > 0 {
        flush_block(
            &mut blocks,
            current_block_query_begin,
            current_block_subject_begin,
            current_block_len,
            record.orientation,
        );
    }
    if in_deletion {
        caf.push('-');
    } else if in_insertion {
        caf.push('+');
    }
    flush_cigar(&mut cigar, prev_op, op_len);

    RepeatMaskerAlignmentEncoding {
        cigar,
        caf,
        match_bases,
        mismatch_bases,
        insertion_bases,
        deletion_bases,
        blocks,
    }
}

fn format_alignment_blocks(blocks: &[RepeatMaskerAlignmentBlock], query_axis: bool) -> String {
    blocks
        .iter()
        .map(|block| {
            if query_axis {
                format!("{}-{}", block.query_begin, block.query_end)
            } else {
                format!("{}-{}", block.subject_begin, block.subject_end)
            }
        })
        .collect::<Vec<_>>()
        .join(";")
}

pub fn write_alignment_tsv_path(
    path: &Path,
    alignments: &[RepeatMaskerAlignmentRecord],
) -> Result<(), RepeatMaskerError> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "annotation_id\tlineage_id\tquery_name\tquery_begin\tquery_end\tsubject_name\tsubject_begin\tsubject_end\tstrand\taligned_query_bases\taligned_subject_bases\tmatch_bases\tmismatch_bases\tinsertion_bases\tdeletion_bases\tcigar\tcaf\tblock_count\tquery_blocks\tsubject_blocks\tmatrix_name\tkimura_divergence\traw_kimura_divergence\tcpg_sites\ttransition_count\ttransversion_count\tgap_count\taverage_gap_size"
    )?;
    for record in alignments {
        let encoding = encode_repeatmasker_alignment(record);
        writeln!(
            writer,
            "{annotation_id}\t{lineage_id}\t{query_name}\t{query_begin}\t{query_end}\t{subject_name}\t{subject_begin}\t{subject_end}\t{strand}\t{aligned_query_bases}\t{aligned_subject_bases}\t{match_bases}\t{mismatch_bases}\t{insertion_bases}\t{deletion_bases}\t{cigar}\t{caf}\t{block_count}\t{query_blocks}\t{subject_blocks}\t{matrix_name}\t{kimura_divergence}\t{raw_kimura_divergence}\t{cpg_sites}\t{transition_count}\t{transversion_count}\t{gap_count}\t{average_gap_size}",
            annotation_id = record.annotation_id,
            lineage_id = record.lineage_id,
            query_name = record.query_name,
            query_begin = record.query_begin,
            query_end = record.query_end,
            subject_name = record.subject_name,
            subject_begin = record.subject_begin,
            subject_end = record.subject_end,
            strand = if record.orientation == 'C' { "-" } else { "+" },
            aligned_query_bases = record.metrics.aligned_query_bases,
            aligned_subject_bases = record.metrics.aligned_subject_bases,
            match_bases = encoding.match_bases,
            mismatch_bases = encoding.mismatch_bases,
            insertion_bases = encoding.insertion_bases,
            deletion_bases = encoding.deletion_bases,
            cigar = encoding.cigar,
            caf = encoding.caf,
            block_count = encoding.blocks.len(),
            query_blocks = format_alignment_blocks(&encoding.blocks, true),
            subject_blocks = format_alignment_blocks(&encoding.blocks, false),
            matrix_name = record.metrics.matrix_name.as_deref().unwrap_or(""),
            kimura_divergence = record.metrics.kimura_divergence.map(|value| value.to_string()).unwrap_or_default(),
            raw_kimura_divergence = record.metrics.raw_kimura_divergence.map(|value| value.to_string()).unwrap_or_default(),
            cpg_sites = record.metrics.cpg_sites.map(|value| value.to_string()).unwrap_or_default(),
            transition_count = record.metrics.transition_count.map(|value| value.to_string()).unwrap_or_default(),
            transversion_count = record.metrics.transversion_count.map(|value| value.to_string()).unwrap_or_default(),
            gap_count = record.metrics.gap_count.map(|value| value.to_string()).unwrap_or_default(),
            average_gap_size = record.metrics.average_gap_size.map(|value| value.to_string()).unwrap_or_default(),
        )?;
    }
    writer.flush()?;
    Ok(())
}

pub fn parse_out_annotation_line(
    line: &str,
    line_number: usize,
) -> Result<Option<RepeatMaskerAnnotation>, RepeatMaskerError> {
    let trimmed = line.trim();
    if trimmed.is_empty() {
        return Ok(None);
    }
    let first = trimmed.as_bytes()[0];
    if !first.is_ascii_digit() {
        return Ok(None);
    }

    let fields: Vec<&str> = trimmed.split_whitespace().collect();
    if !looks_like_annotation_prefix(&fields) {
        return Ok(None);
    }
    if fields.len() < 15 {
        return Err(RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message: format!(
                "expected at least 15 whitespace-separated fields, got {}",
                fields.len()
            ),
        });
    }
    if !matches!(fields[8], "+" | "C" | "c") {
        return Err(RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message: format!("invalid orientation field '{}'", fields[8]),
        });
    }
    if !looks_like_out_annotation_fields(&fields) {
        return Err(RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message: "line does not match RepeatMasker .out field layout".to_string(),
        });
    }

    let score =
        fields[0]
            .parse::<usize>()
            .map_err(|err| RepeatMaskerError::InvalidAnnotationLine {
                line_number,
                line: line.to_string(),
                message: format!("invalid score '{}': {err}", fields[0]),
            })?;
    let percent_divergence_tenths = parse_percent_tenths(fields[1]).map_err(|message| {
        RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message,
        }
    })?;
    let percent_deletions_tenths = parse_percent_tenths(fields[2]).map_err(|message| {
        RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message,
        }
    })?;
    let percent_insertions_tenths = parse_percent_tenths(fields[3]).map_err(|message| {
        RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message,
        }
    })?;
    let query_begin =
        fields[5]
            .parse::<usize>()
            .map_err(|err| RepeatMaskerError::InvalidAnnotationLine {
                line_number,
                line: line.to_string(),
                message: format!("invalid query begin '{}': {err}", fields[5]),
            })?;
    let query_end =
        fields[6]
            .parse::<usize>()
            .map_err(|err| RepeatMaskerError::InvalidAnnotationLine {
                line_number,
                line: line.to_string(),
                message: format!("invalid query end '{}': {err}", fields[6]),
            })?;
    if query_begin == 0 || query_end < query_begin {
        return Err(RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message: format!("invalid query interval {query_begin}-{query_end}"),
        });
    }
    let query_left = parse_parenthesized_usize(fields[7]).map_err(|message| {
        RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message,
        }
    })?;

    let orientation = match fields[8] {
        "+" => '+',
        "C" | "c" => 'C',
        other => {
            return Err(RepeatMaskerError::InvalidAnnotationLine {
                line_number,
                line: line.to_string(),
                message: format!("unsupported orientation '{other}'"),
            })
        }
    };

    let (subject_begin, subject_end, subject_left) = if orientation == 'C' {
        let subject_left = parse_parenthesized_usize(fields[11]).map_err(|message| {
            RepeatMaskerError::InvalidAnnotationLine {
                line_number,
                line: line.to_string(),
                message,
            }
        })?;
        let subject_end = fields[12].parse::<usize>().map_err(|err| {
            RepeatMaskerError::InvalidAnnotationLine {
                line_number,
                line: line.to_string(),
                message: format!("invalid reverse subject end '{}': {err}", fields[12]),
            }
        })?;
        let subject_begin = fields[13].parse::<usize>().map_err(|err| {
            RepeatMaskerError::InvalidAnnotationLine {
                line_number,
                line: line.to_string(),
                message: format!("invalid reverse subject begin '{}': {err}", fields[13]),
            }
        })?;
        (subject_begin, subject_end, subject_left)
    } else {
        let subject_begin = fields[11].parse::<usize>().map_err(|err| {
            RepeatMaskerError::InvalidAnnotationLine {
                line_number,
                line: line.to_string(),
                message: format!("invalid subject begin '{}': {err}", fields[11]),
            }
        })?;
        let subject_end = fields[12].parse::<usize>().map_err(|err| {
            RepeatMaskerError::InvalidAnnotationLine {
                line_number,
                line: line.to_string(),
                message: format!("invalid subject end '{}': {err}", fields[12]),
            }
        })?;
        let subject_left = parse_parenthesized_usize(fields[13]).map_err(|message| {
            RepeatMaskerError::InvalidAnnotationLine {
                line_number,
                line: line.to_string(),
                message,
            }
        })?;
        (subject_begin, subject_end, subject_left)
    };

    let tail = &fields[15..];
    let (lineage_id, overlap_marker) = match tail {
        [] => (String::new(), None),
        [marker] if *marker == "*" => (String::new(), Some((*marker).to_string())),
        [lineage] => ((*lineage).to_string(), None),
        [lineage, overlap, ..] => ((*lineage).to_string(), Some((*overlap).to_string())),
    };

    Ok(Some(RepeatMaskerAnnotation {
        score,
        percent_divergence_tenths,
        percent_deletions_tenths,
        percent_insertions_tenths,
        query_name: fields[4].to_string(),
        query_begin,
        query_end,
        query_left,
        orientation,
        subject_name: fields[9].to_string(),
        subject_class: fields[10].to_string(),
        subject_begin,
        subject_end,
        subject_left,
        annotation_id: fields[14].to_string(),
        lineage_id,
        overlap_marker,
    }))
}

pub fn parse_cat_annotation_line(
    line: &str,
    line_number: usize,
) -> Result<Option<RepeatMaskerAnnotation>, RepeatMaskerError> {
    let trimmed = line.trim();
    if trimmed.is_empty() {
        return Ok(None);
    }
    let first = trimmed.as_bytes()[0];
    if !first.is_ascii_digit() {
        return Ok(None);
    }

    let fields: Vec<&str> = trimmed.split_whitespace().collect();
    if !looks_like_annotation_prefix(&fields) || looks_like_out_annotation_fields(&fields) {
        return Ok(None);
    }
    if !looks_like_cat_annotation_fields(&fields) {
        return Err(RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message: "line does not match RepeatMasker .cat field layout".to_string(),
        });
    }

    let score =
        fields[0]
            .parse::<usize>()
            .map_err(|err| RepeatMaskerError::InvalidAnnotationLine {
                line_number,
                line: line.to_string(),
                message: format!("invalid score '{}': {err}", fields[0]),
            })?;
    let percent_divergence_tenths = parse_percent_tenths(fields[1]).map_err(|message| {
        RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message,
        }
    })?;
    let percent_deletions_tenths = parse_percent_tenths(fields[2]).map_err(|message| {
        RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message,
        }
    })?;
    let percent_insertions_tenths = parse_percent_tenths(fields[3]).map_err(|message| {
        RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message,
        }
    })?;
    let query_begin =
        fields[5]
            .parse::<usize>()
            .map_err(|err| RepeatMaskerError::InvalidAnnotationLine {
                line_number,
                line: line.to_string(),
                message: format!("invalid query begin '{}': {err}", fields[5]),
            })?;
    let query_end =
        fields[6]
            .parse::<usize>()
            .map_err(|err| RepeatMaskerError::InvalidAnnotationLine {
                line_number,
                line: line.to_string(),
                message: format!("invalid query end '{}': {err}", fields[6]),
            })?;
    if query_begin == 0 || query_end < query_begin {
        return Err(RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message: format!("invalid query interval {query_begin}-{query_end}"),
        });
    }
    let query_left = parse_parenthesized_usize(fields[7]).map_err(|message| {
        RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message,
        }
    })?;

    let (orientation, subject_token, subject_coord1, subject_coord2, subject_coord3, tail_index) =
        if fields[8].eq_ignore_ascii_case("C") {
            ('C', fields[9], fields[10], fields[11], fields[12], 13usize)
        } else {
            ('+', fields[8], fields[9], fields[10], fields[11], 12usize)
        };

    let (subject_name, subject_class) = split_subject_token(subject_token);
    let (subject_begin, subject_end, subject_left) = if orientation == 'C' {
        (
            subject_coord3.parse::<usize>().map_err(|err| {
                RepeatMaskerError::InvalidAnnotationLine {
                    line_number,
                    line: line.to_string(),
                    message: format!("invalid reverse subject begin '{subject_coord3}': {err}"),
                }
            })?,
            subject_coord2.parse::<usize>().map_err(|err| {
                RepeatMaskerError::InvalidAnnotationLine {
                    line_number,
                    line: line.to_string(),
                    message: format!("invalid reverse subject end '{subject_coord2}': {err}"),
                }
            })?,
            parse_parenthesized_usize(subject_coord1).map_err(|message| {
                RepeatMaskerError::InvalidAnnotationLine {
                    line_number,
                    line: line.to_string(),
                    message,
                }
            })?,
        )
    } else {
        (
            subject_coord1.parse::<usize>().map_err(|err| {
                RepeatMaskerError::InvalidAnnotationLine {
                    line_number,
                    line: line.to_string(),
                    message: format!("invalid subject begin '{subject_coord1}': {err}"),
                }
            })?,
            subject_coord2.parse::<usize>().map_err(|err| {
                RepeatMaskerError::InvalidAnnotationLine {
                    line_number,
                    line: line.to_string(),
                    message: format!("invalid subject end '{subject_coord2}': {err}"),
                }
            })?,
            parse_parenthesized_usize(subject_coord3).map_err(|message| {
                RepeatMaskerError::InvalidAnnotationLine {
                    line_number,
                    line: line.to_string(),
                    message,
                }
            })?,
        )
    };
    if subject_begin == 0 || subject_end < subject_begin {
        return Err(RepeatMaskerError::InvalidAnnotationLine {
            line_number,
            line: line.to_string(),
            message: format!("invalid subject interval {subject_begin}-{subject_end}"),
        });
    }

    let tail = &fields[tail_index..];
    let (lineage_id, annotation_id, overlap_marker) = match tail {
        [lineage, annotation, overlap, ..] => (
            (*lineage).to_string(),
            (*annotation).to_string(),
            Some((*overlap).to_string()),
        ),
        [lineage, annotation] => ((*lineage).to_string(), (*annotation).to_string(), None),
        [lineage] => ((*lineage).to_string(), (*lineage).to_string(), None),
        [] => (String::new(), format!("cat-{line_number}"), None),
    };

    Ok(Some(RepeatMaskerAnnotation {
        score,
        percent_divergence_tenths,
        percent_deletions_tenths,
        percent_insertions_tenths,
        query_name: fields[4].to_string(),
        query_begin,
        query_end,
        query_left,
        orientation,
        subject_name,
        subject_class,
        subject_begin,
        subject_end,
        subject_left,
        annotation_id,
        lineage_id,
        overlap_marker,
    }))
}

pub fn load_out_annotations(
    reader: impl BufRead,
) -> Result<Vec<RepeatMaskerAnnotation>, RepeatMaskerError> {
    let mut annotations = Vec::new();
    for (index, line) in reader.lines().enumerate() {
        let line = line?;
        if let Some(annotation) = parse_out_annotation_line(&line, index + 1)? {
            annotations.push(annotation);
        }
    }
    Ok(annotations)
}

pub fn load_out_annotations_path(
    path: &Path,
) -> Result<Vec<RepeatMaskerAnnotation>, RepeatMaskerError> {
    load_out_annotations(open_maybe_gzip(path)?)
}

pub fn load_repeatmasker_catalog(
    reader: impl BufRead,
) -> Result<RepeatMaskerCatalog, RepeatMaskerError> {
    let mut annotations = Vec::new();
    let mut metadata = RepeatMaskerMetadata::default();

    for (index, line) in reader.lines().enumerate() {
        let line = line?;
        let line_number = index + 1;
        let fields: Vec<&str> = line.split_whitespace().collect();

        if looks_like_out_annotation_fields(&fields) {
            if let Some(annotation) = parse_out_annotation_line(&line, line_number)? {
                annotations.push(annotation);
            }
            continue;
        }
        if looks_like_cat_annotation_fields(&fields) {
            if let Some(annotation) = parse_cat_annotation_line(&line, line_number)? {
                annotations.push(annotation);
            }
            continue;
        }
        if let Some(annotation) = parse_out_annotation_line(&line, line_number)? {
            annotations.push(annotation);
            continue;
        }
        if let Some(annotation) = parse_cat_annotation_line(&line, line_number)? {
            annotations.push(annotation);
            continue;
        }
        update_metadata_from_line(&line, &mut metadata);
    }

    Ok(RepeatMaskerCatalog {
        annotations,
        metadata,
    })
}

pub fn load_repeatmasker_catalog_path(
    path: &Path,
) -> Result<RepeatMaskerCatalog, RepeatMaskerError> {
    load_repeatmasker_catalog(open_maybe_gzip(path)?)
}

fn annotation_sort_key(
    left: &RepeatMaskerAnnotation,
    right: &RepeatMaskerAnnotation,
) -> std::cmp::Ordering {
    left.query_name
        .cmp(&right.query_name)
        .then(left.query_begin.cmp(&right.query_begin))
        .then(right.query_end.cmp(&left.query_end))
        .then(right.score.cmp(&left.score))
        .then(left.subject_class.cmp(&right.subject_class))
        .then(left.lineage_id.cmp(&right.lineage_id))
        .then(left.subject_name.cmp(&right.subject_name))
        .then(left.annotation_id.cmp(&right.annotation_id))
}

fn alignment_metrics_index(
    alignments: &[RepeatMaskerAlignmentRecord],
) -> BTreeMap<AlignmentLinkKey, RepeatMaskerAlignmentMetrics> {
    alignments
        .iter()
        .map(|record| (alignment_record_key(record), record.metrics.clone()))
        .collect()
}

fn cmp_optional_f64(left: Option<f64>, right: Option<f64>) -> std::cmp::Ordering {
    match (left, right) {
        (Some(left), Some(right)) => left
            .partial_cmp(&right)
            .unwrap_or(std::cmp::Ordering::Equal),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => std::cmp::Ordering::Equal,
    }
}

fn alignment_metric_is_no_worse(
    previous: &RepeatMaskerAlignmentMetrics,
    current: &RepeatMaskerAlignmentMetrics,
) -> bool {
    previous.aligned_query_bases >= current.aligned_query_bases
        && previous.aligned_subject_bases >= current.aligned_subject_bases
        && !matches!(
            cmp_optional_f64(previous.kimura_divergence, current.kimura_divergence),
            std::cmp::Ordering::Greater
        )
        && !matches!(
            cmp_optional_f64(previous.average_gap_size, current.average_gap_size),
            std::cmp::Ordering::Greater
        )
}

fn annotation_sort_key_with_alignments(
    left: &RepeatMaskerAnnotation,
    right: &RepeatMaskerAnnotation,
    alignment_metrics: &BTreeMap<AlignmentLinkKey, RepeatMaskerAlignmentMetrics>,
) -> std::cmp::Ordering {
    left.query_name
        .cmp(&right.query_name)
        .then(left.query_begin.cmp(&right.query_begin))
        .then(right.query_end.cmp(&left.query_end))
        .then(right.score.cmp(&left.score))
        .then_with(|| {
            let left_metrics = alignment_metrics.get(&annotation_alignment_key(left));
            let right_metrics = alignment_metrics.get(&annotation_alignment_key(right));
            cmp_optional_f64(
                left_metrics.and_then(|metrics| metrics.kimura_divergence),
                right_metrics.and_then(|metrics| metrics.kimura_divergence),
            )
        })
        .then_with(|| {
            let left_metrics = alignment_metrics.get(&annotation_alignment_key(left));
            let right_metrics = alignment_metrics.get(&annotation_alignment_key(right));
            cmp_optional_f64(
                left_metrics.and_then(|metrics| metrics.average_gap_size),
                right_metrics.and_then(|metrics| metrics.average_gap_size),
            )
        })
        .then_with(|| {
            let left_bases = alignment_metrics
                .get(&annotation_alignment_key(left))
                .map(|metrics| metrics.aligned_subject_bases)
                .unwrap_or(0);
            let right_bases = alignment_metrics
                .get(&annotation_alignment_key(right))
                .map(|metrics| metrics.aligned_subject_bases)
                .unwrap_or(0);
            right_bases.cmp(&left_bases)
        })
        .then(left.subject_class.cmp(&right.subject_class))
        .then(left.lineage_id.cmp(&right.lineage_id))
        .then(left.subject_name.cmp(&right.subject_name))
        .then(left.annotation_id.cmp(&right.annotation_id))
}

fn is_short_fragment(annotation: &RepeatMaskerAnnotation) -> bool {
    let length = annotation.query_end - annotation.query_begin + 1;
    length < 5 || (!annotation.subject_class.contains("Simple") && length < 10)
}

fn signed_interval_overlap(
    left_begin: usize,
    left_end: usize,
    right_begin: usize,
    right_end: usize,
) -> isize {
    if right_begin <= left_end && left_begin <= right_end {
        left_end.min(right_end) as isize - left_begin.max(right_begin) as isize + 1
    } else if left_begin <= right_begin {
        left_end as isize - right_begin as isize + 1
    } else {
        right_end as isize - left_begin as isize + 1
    }
}

fn signed_interval_overlap_i64(
    left_begin: i64,
    left_end: i64,
    right_begin: i64,
    right_end: i64,
) -> i64 {
    if right_begin <= left_end && left_begin <= right_end {
        left_end.min(right_end) - left_begin.max(right_begin) + 1
    } else if left_begin <= right_begin {
        left_end - right_begin + 1
    } else {
        right_end - left_begin + 1
    }
}

fn annotations_can_chain(
    previous: &RepeatMaskerAnnotation,
    current: &RepeatMaskerAnnotation,
) -> bool {
    if annotation_ids_imply_same_element(previous, current) {
        return true;
    }

    if previous.query_name != current.query_name
        || previous.orientation != current.orientation
        || previous.subject_name != current.subject_name
        || previous.subject_class != current.subject_class
    {
        return false;
    }

    let query_overlap = signed_interval_overlap(
        previous.query_begin,
        previous.query_end,
        current.query_begin,
        current.query_end,
    );
    let consensus_overlap = signed_interval_overlap(
        previous.subject_begin,
        previous.subject_end,
        current.subject_begin,
        current.subject_end,
    );
    let strict_join = query_overlap >= -10
        && consensus_overlap >= -100
        && (consensus_overlap <= 21
            || (query_overlap < 0 && consensus_overlap - query_overlap <= 20));
    if strict_join {
        return true;
    }

    let divergence_delta = previous
        .percent_divergence_tenths
        .abs_diff(current.percent_divergence_tenths);
    if divergence_delta >= 100 {
        return false;
    }

    let gap_scale = (100 - divergence_delta) as isize;
    let gap_max = 3750isize * gap_scale / 100;
    query_overlap > -gap_max && consensus_overlap < 33 && consensus_overlap > -gap_max
}

fn annotation_ids_imply_same_element(
    previous: &RepeatMaskerAnnotation,
    current: &RepeatMaskerAnnotation,
) -> bool {
    previous.query_name == current.query_name
        && previous.orientation == current.orientation
        && previous.subject_class == current.subject_class
        && !previous.annotation_id.is_empty()
        && previous.annotation_id == current.annotation_id
        // Old `.cat` fixtures use small numeric stage values here rather than
        // stable element IDs, so only trust IDs on `.out`-style annotations.
        && previous.lineage_id.is_empty()
        && current.lineage_id.is_empty()
}

fn alignment_encoding_index(
    alignments: &[RepeatMaskerAlignmentRecord],
) -> BTreeMap<AlignmentLinkKey, RepeatMaskerAlignmentEncoding> {
    alignments
        .iter()
        .map(|record| {
            (
                alignment_record_key(record),
                encode_repeatmasker_alignment(record),
            )
        })
        .collect()
}

fn annotation_query_span(
    annotation: &RepeatMaskerAnnotation,
    encoding: Option<&RepeatMaskerAlignmentEncoding>,
) -> (usize, usize) {
    if let Some(encoding) = encoding {
        if let (Some(first), Some(last)) = (encoding.blocks.first(), encoding.blocks.last()) {
            return (first.query_begin, last.query_end);
        }
    }
    (annotation.query_begin, annotation.query_end)
}

fn annotation_oriented_subject_span(
    annotation: &RepeatMaskerAnnotation,
    encoding: Option<&RepeatMaskerAlignmentEncoding>,
) -> (i64, i64) {
    if let Some(encoding) = encoding {
        if let (Some(first), Some(last)) = (encoding.blocks.first(), encoding.blocks.last()) {
            return if annotation.orientation == 'C' {
                (-(first.subject_end as i64), -(last.subject_begin as i64))
            } else {
                (first.subject_begin as i64, last.subject_end as i64)
            };
        }
    }

    if annotation.orientation == 'C' {
        (
            -(annotation.subject_end as i64),
            -(annotation.subject_begin as i64),
        )
    } else {
        (
            annotation.subject_begin as i64,
            annotation.subject_end as i64,
        )
    }
}

fn annotations_can_chain_with_encoding_index(
    previous: &RepeatMaskerAnnotation,
    current: &RepeatMaskerAnnotation,
    encodings: &BTreeMap<AlignmentLinkKey, RepeatMaskerAlignmentEncoding>,
) -> bool {
    if annotations_can_chain(previous, current) {
        return true;
    }
    if previous.query_name != current.query_name
        || previous.orientation != current.orientation
        || previous.subject_name != current.subject_name
        || previous.subject_class != current.subject_class
    {
        return false;
    }

    let previous_encoding = encodings.get(&annotation_alignment_key(previous));
    let current_encoding = encodings.get(&annotation_alignment_key(current));
    let (previous_query_begin, previous_query_end) =
        annotation_query_span(previous, previous_encoding);
    let (current_query_begin, current_query_end) = annotation_query_span(current, current_encoding);
    let (previous_subject_begin, previous_subject_end) =
        annotation_oriented_subject_span(previous, previous_encoding);
    let (current_subject_begin, current_subject_end) =
        annotation_oriented_subject_span(current, current_encoding);

    let query_overlap = signed_interval_overlap_i64(
        previous_query_begin as i64,
        previous_query_end as i64,
        current_query_begin as i64,
        current_query_end as i64,
    );
    let consensus_overlap = signed_interval_overlap_i64(
        previous_subject_begin,
        previous_subject_end,
        current_subject_begin,
        current_subject_end,
    );
    let strict_join = query_overlap >= -10
        && consensus_overlap >= -100
        && (consensus_overlap <= 21
            || (query_overlap < 0 && consensus_overlap - query_overlap <= 20));
    if strict_join {
        return true;
    }

    let divergence_delta = previous
        .percent_divergence_tenths
        .abs_diff(current.percent_divergence_tenths);
    if divergence_delta >= 100 {
        return false;
    }

    let gap_scale = (100 - divergence_delta) as i64;
    let gap_max = 3750i64 * gap_scale / 100;
    query_overlap > -gap_max && consensus_overlap < 33 && consensus_overlap > -gap_max
}

pub fn adjudicate_annotations(
    annotations: &[RepeatMaskerAnnotation],
) -> Vec<RepeatMaskerAnnotation> {
    let mut sorted = annotations.to_vec();
    sorted.sort_unstable_by(annotation_sort_key);

    let mut adjudicated = Vec::with_capacity(sorted.len());
    let mut previous_kept: Option<&RepeatMaskerAnnotation> = None;

    for annotation in sorted {
        if is_short_fragment(&annotation) {
            continue;
        }

        let suppress_as_duplicate = previous_kept.is_some_and(|previous| {
            annotation.query_name == previous.query_name
                && annotation.query_end <= previous.query_end
                && annotation.score <= previous.score
                && annotation.subject_class == previous.subject_class
                && annotation.lineage_id == previous.lineage_id
        });
        if suppress_as_duplicate {
            continue;
        }

        adjudicated.push(annotation);
        previous_kept = adjudicated.last();
    }

    adjudicated
}

pub fn adjudicate_annotations_with_alignments(
    annotations: &[RepeatMaskerAnnotation],
    alignments: &[RepeatMaskerAlignmentRecord],
) -> Vec<RepeatMaskerAnnotation> {
    let alignment_metrics = alignment_metrics_index(alignments);
    let mut sorted = annotations.to_vec();
    sorted.sort_unstable_by(|left, right| {
        annotation_sort_key_with_alignments(left, right, &alignment_metrics)
    });

    let mut adjudicated = Vec::with_capacity(sorted.len());
    let mut previous_kept: Option<&RepeatMaskerAnnotation> = None;

    for annotation in sorted {
        if is_short_fragment(&annotation) {
            continue;
        }

        let suppress_as_duplicate = previous_kept.is_some_and(|previous| {
            if annotation.query_name != previous.query_name
                || annotation.query_begin < previous.query_begin
                || annotation.query_end > previous.query_end
                || annotation.subject_name != previous.subject_name
                || annotation.subject_class != previous.subject_class
                || annotation.lineage_id != previous.lineage_id
            {
                return false;
            }

            let current_metrics = alignment_metrics.get(&annotation_alignment_key(&annotation));
            let previous_metrics = alignment_metrics.get(&annotation_alignment_key(previous));

            previous.score >= annotation.score
                && match (previous_metrics, current_metrics) {
                    (Some(previous_metrics), Some(current_metrics)) => {
                        alignment_metric_is_no_worse(previous_metrics, current_metrics)
                    }
                    _ => true,
                }
        });
        if suppress_as_duplicate {
            continue;
        }

        adjudicated.push(annotation);
        previous_kept = adjudicated.last();
    }

    adjudicated
}

pub fn build_annotation_chains(
    annotations: &[RepeatMaskerAnnotation],
) -> Vec<RepeatMaskerAnnotationChain> {
    struct ChainAccumulator {
        chain: RepeatMaskerAnnotationChain,
        last_annotation: RepeatMaskerAnnotation,
    }

    let mut sorted = annotations.to_vec();
    sorted.sort_unstable_by(annotation_sort_key);

    let mut chains: Vec<ChainAccumulator> = Vec::new();
    for annotation in sorted {
        if let Some(previous_chain) = chains.last_mut() {
            if annotations_can_chain(&previous_chain.last_annotation, &annotation) {
                previous_chain.chain.query_end =
                    previous_chain.chain.query_end.max(annotation.query_end);
                previous_chain.chain.subject_end =
                    previous_chain.chain.subject_end.max(annotation.subject_end);
                previous_chain.chain.annotation_count += 1;
                previous_chain
                    .chain
                    .annotation_ids
                    .push(annotation.annotation_id.clone());
                if !annotation.lineage_id.is_empty()
                    && previous_chain.chain.lineage_ids.last() != Some(&annotation.lineage_id)
                {
                    previous_chain
                        .chain
                        .lineage_ids
                        .push(annotation.lineage_id.clone());
                }
                previous_chain.last_annotation = annotation;
                continue;
            }
        }

        let mut lineage_ids = Vec::new();
        if !annotation.lineage_id.is_empty() {
            lineage_ids.push(annotation.lineage_id.clone());
        }
        chains.push(ChainAccumulator {
            chain: RepeatMaskerAnnotationChain {
                chain_id: chains.len() + 1,
                query_name: annotation.query_name.clone(),
                orientation: annotation.orientation,
                subject_name: annotation.subject_name.clone(),
                subject_class: annotation.subject_class.clone(),
                query_begin: annotation.query_begin,
                query_end: annotation.query_end,
                subject_begin: annotation.subject_begin,
                subject_end: annotation.subject_end,
                annotation_count: 1,
                annotation_ids: vec![annotation.annotation_id.clone()],
                lineage_ids,
            },
            last_annotation: annotation,
        });
    }
    chains
        .into_iter()
        .map(|accumulator| accumulator.chain)
        .collect()
}

pub fn build_annotation_chains_with_alignments(
    annotations: &[RepeatMaskerAnnotation],
    alignments: &[RepeatMaskerAlignmentRecord],
) -> Vec<RepeatMaskerAnnotationChain> {
    struct ChainAccumulator {
        chain: RepeatMaskerAnnotationChain,
        last_annotation: RepeatMaskerAnnotation,
    }

    let alignment_metrics = alignment_metrics_index(alignments);
    let alignment_encodings = alignment_encoding_index(alignments);

    let mut sorted = annotations.to_vec();
    sorted.sort_unstable_by(|left, right| {
        annotation_sort_key_with_alignments(left, right, &alignment_metrics)
    });

    let mut chains: Vec<ChainAccumulator> = Vec::new();
    for annotation in sorted {
        if let Some(previous_chain) = chains.last_mut() {
            if annotations_can_chain_with_encoding_index(
                &previous_chain.last_annotation,
                &annotation,
                &alignment_encodings,
            ) {
                previous_chain.chain.query_end =
                    previous_chain.chain.query_end.max(annotation.query_end);
                previous_chain.chain.subject_end =
                    previous_chain.chain.subject_end.max(annotation.subject_end);
                previous_chain.chain.annotation_count += 1;
                previous_chain
                    .chain
                    .annotation_ids
                    .push(annotation.annotation_id.clone());
                if !annotation.lineage_id.is_empty()
                    && previous_chain.chain.lineage_ids.last() != Some(&annotation.lineage_id)
                {
                    previous_chain
                        .chain
                        .lineage_ids
                        .push(annotation.lineage_id.clone());
                }
                previous_chain.last_annotation = annotation;
                continue;
            }
        }

        let mut lineage_ids = Vec::new();
        if !annotation.lineage_id.is_empty() {
            lineage_ids.push(annotation.lineage_id.clone());
        }
        chains.push(ChainAccumulator {
            chain: RepeatMaskerAnnotationChain {
                chain_id: chains.len() + 1,
                query_name: annotation.query_name.clone(),
                orientation: annotation.orientation,
                subject_name: annotation.subject_name.clone(),
                subject_class: annotation.subject_class.clone(),
                query_begin: annotation.query_begin,
                query_end: annotation.query_end,
                subject_begin: annotation.subject_begin,
                subject_end: annotation.subject_end,
                annotation_count: 1,
                annotation_ids: vec![annotation.annotation_id.clone()],
                lineage_ids,
            },
            last_annotation: annotation,
        });
    }

    chains
        .into_iter()
        .map(|accumulator| accumulator.chain)
        .collect()
}

pub fn annotation_stats(annotations: &[RepeatMaskerAnnotation]) -> AnnotationStats {
    AnnotationStats {
        annotation_count: annotations.len(),
        total_aligned_bases: annotations
            .iter()
            .map(|annotation| annotation.query_end - annotation.query_begin + 1)
            .sum(),
        occupied_bases: occupied_bases(annotations),
        sequence_count: annotations
            .iter()
            .map(|annotation| annotation.query_name.as_str())
            .collect::<std::collections::BTreeSet<_>>()
            .len(),
    }
}

fn occupied_bases(annotations: &[RepeatMaskerAnnotation]) -> usize {
    consolidate_mask_intervals(annotations)
        .values()
        .flatten()
        .map(|interval| interval.end - interval.begin + 1)
        .sum()
}

fn summarize_annotations_by<F>(
    annotations: &[RepeatMaskerAnnotation],
    mut key_fn: F,
) -> Vec<RepeatSummaryRow>
where
    F: FnMut(&RepeatMaskerAnnotation) -> &str,
{
    let mut grouped: BTreeMap<String, Vec<RepeatMaskerAnnotation>> = BTreeMap::new();
    for annotation in annotations {
        grouped
            .entry(key_fn(annotation).to_string())
            .or_default()
            .push(annotation.clone());
    }

    let mut rows: Vec<RepeatSummaryRow> = grouped
        .into_iter()
        .map(|(key, grouped_annotations)| RepeatSummaryRow {
            key,
            annotation_count: grouped_annotations.len(),
            aligned_bases: grouped_annotations
                .iter()
                .map(|annotation| annotation.query_end - annotation.query_begin + 1)
                .sum(),
            occupied_bases: occupied_bases(&grouped_annotations),
        })
        .collect();
    rows.sort_unstable_by(|left, right| {
        right
            .occupied_bases
            .cmp(&left.occupied_bases)
            .then(right.aligned_bases.cmp(&left.aligned_bases))
            .then(right.annotation_count.cmp(&left.annotation_count))
            .then(left.key.cmp(&right.key))
    });
    rows
}

fn out_style_annotation_id_count(annotations: &[RepeatMaskerAnnotation]) -> Option<usize> {
    if annotations
        .iter()
        .all(|annotation| !annotation.annotation_id.is_empty() && annotation.lineage_id.is_empty())
    {
        Some(
            annotations
                .iter()
                .map(|annotation| annotation.annotation_id.as_str())
                .collect::<BTreeSet<_>>()
                .len(),
        )
    } else {
        None
    }
}

fn default_element_count(annotations: &[RepeatMaskerAnnotation]) -> usize {
    out_style_annotation_id_count(annotations)
        .unwrap_or_else(|| build_annotation_chains(annotations).len())
}

fn default_element_count_with_alignments(
    annotations: &[RepeatMaskerAnnotation],
    alignments: &[RepeatMaskerAlignmentRecord],
) -> usize {
    out_style_annotation_id_count(annotations)
        .unwrap_or_else(|| build_annotation_chains_with_alignments(annotations, alignments).len())
}

fn catalog_stats_for_annotations_with_chain_count(
    annotations: &[RepeatMaskerAnnotation],
    raw_annotation_count: Option<usize>,
    metadata: RepeatMaskerMetadata,
    element_chain_count: usize,
) -> RepeatMaskerCatalogStats {
    RepeatMaskerCatalogStats {
        annotation_stats: annotation_stats(annotations),
        raw_annotation_count,
        suppressed_annotation_count: raw_annotation_count.map(|count| count - annotations.len()),
        element_chain_count,
        repeat_class_summary: summarize_annotations_by(annotations, |annotation| {
            if annotation.subject_class.is_empty() {
                "unclassified"
            } else {
                annotation.subject_class.as_str()
            }
        }),
        repeat_family_summary: summarize_annotations_by(annotations, |annotation| {
            annotation.subject_name.as_str()
        }),
        metadata,
    }
}

pub fn catalog_stats_for_annotations(
    annotations: &[RepeatMaskerAnnotation],
    raw_annotation_count: Option<usize>,
    metadata: RepeatMaskerMetadata,
) -> RepeatMaskerCatalogStats {
    catalog_stats_for_annotations_with_chain_count(
        annotations,
        raw_annotation_count,
        metadata,
        default_element_count(annotations),
    )
}

pub fn catalog_stats_for_annotations_with_alignments(
    annotations: &[RepeatMaskerAnnotation],
    raw_annotation_count: Option<usize>,
    metadata: RepeatMaskerMetadata,
    alignments: &[RepeatMaskerAlignmentRecord],
) -> RepeatMaskerCatalogStats {
    catalog_stats_for_annotations_with_chain_count(
        annotations,
        raw_annotation_count,
        metadata,
        default_element_count_with_alignments(annotations, alignments),
    )
}

pub fn catalog_stats(catalog: &RepeatMaskerCatalog) -> RepeatMaskerCatalogStats {
    catalog_stats_for_annotations(&catalog.annotations, None, catalog.metadata.clone())
}

pub fn adjudicated_catalog_stats(catalog: &RepeatMaskerCatalog) -> RepeatMaskerCatalogStats {
    let adjudicated = adjudicate_annotations(&catalog.annotations);
    catalog_stats_for_annotations(
        &adjudicated,
        Some(catalog.annotations.len()),
        catalog.metadata.clone(),
    )
}

pub fn adjudicated_catalog_stats_with_alignments(
    catalog: &RepeatMaskerCatalog,
    alignments: &[RepeatMaskerAlignmentRecord],
) -> RepeatMaskerCatalogStats {
    let adjudicated = adjudicate_annotations_with_alignments(&catalog.annotations, alignments);
    catalog_stats_for_annotations_with_alignments(
        &adjudicated,
        Some(catalog.annotations.len()),
        catalog.metadata.clone(),
        alignments,
    )
}

pub fn consolidate_mask_intervals(
    annotations: &[RepeatMaskerAnnotation],
) -> BTreeMap<String, Vec<MaskInterval>> {
    let mut grouped: BTreeMap<String, Vec<(usize, usize)>> = BTreeMap::new();
    for annotation in annotations {
        grouped
            .entry(annotation.query_name.clone())
            .or_default()
            .push((annotation.query_begin, annotation.query_end));
    }

    let mut intervals = BTreeMap::new();
    for (query_name, mut ranges) in grouped {
        ranges.sort_unstable_by(|left, right| left.0.cmp(&right.0).then(left.1.cmp(&right.1)));

        let mut consolidated: Vec<MaskInterval> = Vec::with_capacity(ranges.len());
        for (begin, end) in ranges {
            if let Some(last) = consolidated.last_mut() {
                if begin <= last.end.saturating_add(1) {
                    if end > last.end {
                        last.end = end;
                    }
                    continue;
                }
            }
            consolidated.push(MaskInterval { begin, end });
        }
        intervals.insert(query_name, consolidated);
    }
    intervals
}

#[inline(always)]
fn count_non_ambiguous_bases_excluding_long_runs(sequence: &[u8]) -> usize {
    let mut remaining = sequence.len();
    let mut index = 0usize;
    while index < sequence.len() {
        let byte = sequence[index];
        if byte != b'N' && byte != b'n' && byte != b'X' && byte != b'x' {
            index += 1;
            continue;
        }
        let run_start = index;
        while index < sequence.len() {
            let current = sequence[index];
            if current != b'N' && current != b'n' && current != b'X' && current != b'x' {
                break;
            }
            index += 1;
        }
        let run_len = index - run_start;
        if run_len >= 20 {
            remaining = remaining.saturating_sub(run_len);
        }
    }
    remaining
}

#[inline(always)]
fn gc_fraction(sequence: &[u8]) -> usize {
    sequence
        .iter()
        .filter(|base| matches!(**base, b'G' | b'g' | b'C' | b'c'))
        .count()
}

pub fn read_fasta_records(reader: impl BufRead) -> Result<Vec<FastaRecord>, RepeatMaskerError> {
    let mut records = Vec::new();
    let mut current_id = String::new();
    let mut current_description = String::new();
    let mut current_sequence = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if !current_id.is_empty() {
                records.push(FastaRecord {
                    id: std::mem::take(&mut current_id),
                    description: std::mem::take(&mut current_description),
                    sequence: std::mem::take(&mut current_sequence),
                });
            }
            let trimmed = header.trim();
            let mut parts = trimmed.splitn(2, char::is_whitespace);
            current_id = parts.next().unwrap_or_default().to_string();
            current_description = parts.next().unwrap_or_default().trim().to_string();
        } else {
            current_sequence.extend_from_slice(line.trim().as_bytes());
        }
    }

    if !current_id.is_empty() {
        records.push(FastaRecord {
            id: current_id,
            description: current_description,
            sequence: current_sequence,
        });
    }

    Ok(records)
}

pub fn read_fasta_records_path(path: &Path) -> Result<Vec<FastaRecord>, RepeatMaskerError> {
    read_fasta_records(open_maybe_gzip(path)?)
}

fn write_wrapped_record(
    writer: &mut dyn Write,
    record: &FastaRecord,
    wrap_width: usize,
) -> Result<(), RepeatMaskerError> {
    if record.description.is_empty() {
        writeln!(writer, ">{}", record.id)?;
    } else {
        writeln!(writer, ">{} {}", record.id, record.description)?;
    }

    let mut start = 0usize;
    while start < record.sequence.len() {
        let end = (start + wrap_width).min(record.sequence.len());
        writer.write_all(&record.sequence[start..end])?;
        writer.write_all(b"\n")?;
        start = end;
    }
    if record.sequence.is_empty() {
        writer.write_all(b"\n")?;
    }
    Ok(())
}

pub fn write_fasta_records_path(
    path: &Path,
    records: &[FastaRecord],
    wrap_width: usize,
) -> Result<(), RepeatMaskerError> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    for record in records {
        write_wrapped_record(&mut writer, record, wrap_width)?;
    }
    writer.flush()?;
    Ok(())
}

#[inline]
pub fn mask_fasta_records(
    records: &mut [FastaRecord],
    intervals: &BTreeMap<String, Vec<MaskInterval>>,
    mask_mode: MaskMode,
) -> MaskStats {
    let mut total_bases = 0usize;
    let mut gc_bases = 0usize;
    let mut masked_bases = 0usize;
    let mut non_ambiguous_bases_excluding_long_runs = 0usize;
    let interval_count = intervals.values().map(Vec::len).sum();

    for record in records.iter_mut() {
        total_bases += record.sequence.len();
        gc_bases += gc_fraction(&record.sequence);
        non_ambiguous_bases_excluding_long_runs +=
            count_non_ambiguous_bases_excluding_long_runs(&record.sequence);

        if let Some(mask_ranges) = intervals.get(&record.id) {
            for range in mask_ranges {
                if range.begin == 0 || range.begin > record.sequence.len() {
                    continue;
                }
                let start = range.begin - 1;
                let end = range.end.min(record.sequence.len());
                if end <= start {
                    continue;
                }
                masked_bases += end - start;
                match mask_mode {
                    MaskMode::N => record.sequence[start..end].fill(b'N'),
                    MaskMode::X => record.sequence[start..end].fill(b'X'),
                    MaskMode::Lowercase => {
                        for base in &mut record.sequence[start..end] {
                            base.make_ascii_lowercase();
                        }
                    }
                }
            }
        }
    }

    MaskStats {
        sequence_count: records.len(),
        total_bases,
        non_ambiguous_bases_excluding_long_runs,
        gc_fraction: if total_bases == 0 {
            0.0
        } else {
            gc_bases as f64 / total_bases as f64
        },
        masked_bases,
        annotation_count: interval_count,
        interval_count,
    }
}

pub fn mask_repeatmasker_out_to_fasta(
    annotation_path: &Path,
    fasta_path: &Path,
    output_path: &Path,
    mask_mode: MaskMode,
    wrap_width: usize,
) -> Result<MaskStats, RepeatMaskerError> {
    let annotations = load_out_annotations_path(annotation_path)?;
    mask_repeatmasker_annotations_to_fasta(
        &annotations,
        fasta_path,
        output_path,
        mask_mode,
        wrap_width,
    )
}

pub fn mask_repeatmasker_annotations_to_fasta(
    annotations: &[RepeatMaskerAnnotation],
    fasta_path: &Path,
    output_path: &Path,
    mask_mode: MaskMode,
    wrap_width: usize,
) -> Result<MaskStats, RepeatMaskerError> {
    compute_repeatmasker_mask_stats_from_fasta(
        annotations,
        fasta_path,
        mask_mode,
        Some((output_path, wrap_width)),
    )
}

pub fn compute_repeatmasker_mask_stats_from_fasta(
    annotations: &[RepeatMaskerAnnotation],
    fasta_path: &Path,
    mask_mode: MaskMode,
    masked_output: Option<(&Path, usize)>,
) -> Result<MaskStats, RepeatMaskerError> {
    let intervals = consolidate_mask_intervals(annotations);
    let mut records = read_fasta_records_path(fasta_path)?;
    let mut stats = mask_fasta_records(&mut records, &intervals, mask_mode);
    stats.annotation_count = annotations.len();
    if let Some((output_path, wrap_width)) = masked_output {
        write_fasta_records_path(output_path, &records, wrap_width)?;
    }
    Ok(stats)
}

fn query_sequence_lengths(annotations: &[RepeatMaskerAnnotation]) -> BTreeMap<&str, usize> {
    let mut lengths = BTreeMap::new();
    for annotation in annotations {
        let inferred = annotation.query_end + annotation.query_left;
        lengths
            .entry(annotation.query_name.as_str())
            .and_modify(|length| {
                if inferred > *length {
                    *length = inferred;
                }
            })
            .or_insert(inferred);
    }
    lengths
}

pub fn write_out_annotations_path(
    path: &Path,
    annotations: &[RepeatMaskerAnnotation],
) -> Result<(), RepeatMaskerError> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    for annotation in annotations {
        let orientation = if annotation.orientation == 'C' {
            "C"
        } else {
            "+"
        };
        let (subject_coord1, subject_coord2, subject_coord3) = if annotation.orientation == 'C' {
            (
                format!("({})", annotation.subject_left),
                annotation.subject_end.to_string(),
                annotation.subject_begin.to_string(),
            )
        } else {
            (
                annotation.subject_begin.to_string(),
                annotation.subject_end.to_string(),
                format!("({})", annotation.subject_left),
            )
        };

        writeln!(
            writer,
            "{score:>6} {div:>4} {del:>4} {ins:>4} {query_name:<17} {query_begin:>8} {query_end:>8} {query_left:>8} {orientation:>1} {subject_name:<15} {subject_class:<15} {subject_coord1:>7} {subject_coord2:>7} {subject_coord3:>7} {annotation_id:<5} {lineage_id:>3}{overlap}",
            score = annotation.score,
            div = format_percent_tenths(annotation.percent_divergence_tenths),
            del = format_percent_tenths(annotation.percent_deletions_tenths),
            ins = format_percent_tenths(annotation.percent_insertions_tenths),
            query_name = annotation.query_name,
            query_begin = annotation.query_begin,
            query_end = annotation.query_end,
            query_left = format!("({})", annotation.query_left),
            orientation = orientation,
            subject_name = annotation.subject_name,
            subject_class = annotation.subject_class,
            subject_coord1 = subject_coord1,
            subject_coord2 = subject_coord2,
            subject_coord3 = subject_coord3,
            annotation_id = annotation.annotation_id,
            lineage_id = annotation.lineage_id,
            overlap = annotation
                .overlap_marker
                .as_ref()
                .map(|value| format!(" {value}"))
                .unwrap_or_default(),
        )?;
    }
    writer.flush()?;
    Ok(())
}

pub fn write_annotations_tsv_path(
    path: &Path,
    annotations: &[RepeatMaskerAnnotation],
) -> Result<(), RepeatMaskerError> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "query_name\tquery_begin\tquery_end\tquery_left\tscore\tpercent_divergence\tpercent_deletions\tpercent_insertions\tstrand\tsubject_name\tsubject_class\tsubject_begin\tsubject_end\tsubject_left\tannotation_id\tlineage_id\toverlap_marker"
    )?;
    for annotation in annotations {
        writeln!(
            writer,
            "{query_name}\t{query_begin}\t{query_end}\t{query_left}\t{score}\t{div}\t{del}\t{ins}\t{strand}\t{subject_name}\t{subject_class}\t{subject_begin}\t{subject_end}\t{subject_left}\t{annotation_id}\t{lineage_id}\t{overlap}",
            query_name = annotation.query_name,
            query_begin = annotation.query_begin,
            query_end = annotation.query_end,
            query_left = annotation.query_left,
            score = annotation.score,
            div = format_percent_tenths(annotation.percent_divergence_tenths),
            del = format_percent_tenths(annotation.percent_deletions_tenths),
            ins = format_percent_tenths(annotation.percent_insertions_tenths),
            strand = if annotation.orientation == 'C' { "-" } else { "+" },
            subject_name = annotation.subject_name,
            subject_class = annotation.subject_class,
            subject_begin = annotation.subject_begin,
            subject_end = annotation.subject_end,
            subject_left = annotation.subject_left,
            annotation_id = annotation.annotation_id,
            lineage_id = annotation.lineage_id,
            overlap = annotation.overlap_marker.as_deref().unwrap_or(""),
        )?;
    }
    writer.flush()?;
    Ok(())
}

pub fn write_gff3_path(
    path: &Path,
    annotations: &[RepeatMaskerAnnotation],
) -> Result<(), RepeatMaskerError> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "##gff-version 3")?;

    for (query_name, length) in query_sequence_lengths(annotations) {
        writeln!(writer, "##sequence-region {query_name} 1 {length}")?;
    }

    for annotation in annotations {
        let strand = if annotation.orientation == 'C' {
            '-'
        } else {
            '+'
        };
        writeln!(
            writer,
            "{query_name}\tRepeatMasker-Rust\tdispersed_repeat\t{query_begin}\t{query_end}\t{score}\t{strand}\t.\tID={id};Name={name};Class={class_name};Target=Motif:{target_name} {subject_begin} {subject_end}",
            query_name = annotation.query_name,
            query_begin = annotation.query_begin,
            query_end = annotation.query_end,
            score = format_percent_tenths(annotation.percent_divergence_tenths),
            strand = strand,
            id = gff_escape(&annotation.annotation_id),
            name = gff_escape(&annotation.subject_name),
            class_name = gff_escape(&annotation.subject_class),
            target_name = gff_escape(&annotation.subject_name),
            subject_begin = annotation.subject_begin,
            subject_end = annotation.subject_end,
        )?;
    }

    writer.flush()?;
    Ok(())
}

pub fn write_repeat_summary_tsv_path(
    path: &Path,
    catalog: &RepeatMaskerCatalog,
) -> Result<(), RepeatMaskerError> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    let annotations = adjudicate_annotations(&catalog.annotations);
    writeln!(
        writer,
        "scope\tkey\tannotation_count\taligned_bases\toccupied_bases"
    )?;
    for row in summarize_annotations_by(&annotations, |annotation| {
        if annotation.subject_class.is_empty() {
            "unclassified"
        } else {
            annotation.subject_class.as_str()
        }
    }) {
        writeln!(
            writer,
            "class\t{}\t{}\t{}\t{}",
            row.key, row.annotation_count, row.aligned_bases, row.occupied_bases
        )?;
    }
    for row in summarize_annotations_by(&annotations, |annotation| annotation.subject_name.as_str())
    {
        writeln!(
            writer,
            "family\t{}\t{}\t{}\t{}",
            row.key, row.annotation_count, row.aligned_bases, row.occupied_bases
        )?;
    }
    writer.flush()?;
    Ok(())
}

pub fn write_chain_tsv_path(
    path: &Path,
    annotations: &[RepeatMaskerAnnotation],
) -> Result<(), RepeatMaskerError> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    let chains = build_annotation_chains(annotations);

    writeln!(
        writer,
        "chain_id\tquery_name\torientation\tsubject_name\tsubject_class\tannotation_count\tquery_begin\tquery_end\tsubject_begin\tsubject_end\tannotation_ids\tlineage_ids"
    )?;
    for chain in chains {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            chain.chain_id,
            chain.query_name,
            chain.orientation,
            chain.subject_name,
            chain.subject_class,
            chain.annotation_count,
            chain.query_begin,
            chain.query_end,
            chain.subject_begin,
            chain.subject_end,
            chain.annotation_ids.join(";"),
            chain.lineage_ids.join(";"),
        )?;
    }

    writer.flush()?;
    Ok(())
}

pub fn write_chain_tsv_with_alignments_path(
    path: &Path,
    annotations: &[RepeatMaskerAnnotation],
    alignments: &[RepeatMaskerAlignmentRecord],
) -> Result<(), RepeatMaskerError> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    let chains = build_annotation_chains_with_alignments(annotations, alignments);

    writeln!(
        writer,
        "chain_id\tquery_name\torientation\tsubject_name\tsubject_class\tannotation_count\tquery_begin\tquery_end\tsubject_begin\tsubject_end\tannotation_ids\tlineage_ids"
    )?;
    for chain in chains {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            chain.chain_id,
            chain.query_name,
            chain.orientation,
            chain.subject_name,
            chain.subject_class,
            chain.annotation_count,
            chain.query_begin,
            chain.query_end,
            chain.subject_begin,
            chain.subject_end,
            chain.annotation_ids.join(";"),
            chain.lineage_ids.join(";"),
        )?;
    }

    writer.flush()?;
    Ok(())
}

fn summarize_class_filter_with_counter<F, C>(
    annotations: &[RepeatMaskerAnnotation],
    mut predicate: F,
    chain_counter: &C,
) -> (usize, usize)
where
    F: FnMut(&str) -> bool,
    C: Fn(&[RepeatMaskerAnnotation]) -> usize,
{
    let filtered: Vec<RepeatMaskerAnnotation> = annotations
        .iter()
        .filter(|annotation| predicate(annotation.subject_class.as_str()))
        .cloned()
        .collect();
    let element_count = chain_counter(&filtered);
    let intervals = consolidate_mask_intervals(&filtered);
    (
        element_count,
        intervals
            .values()
            .flatten()
            .map(|interval| interval.end - interval.begin + 1)
            .sum(),
    )
}

fn percent_of(occupied_bases: usize, total_bases: usize) -> f64 {
    if total_bases == 0 {
        0.0
    } else {
        occupied_bases as f64 * 100.0 / total_bases as f64
    }
}

fn write_tbl_row(
    writer: &mut dyn Write,
    label: &str,
    annotation_count: usize,
    occupied_bases: usize,
    total_bases: usize,
) -> Result<(), RepeatMaskerError> {
    writeln!(
        writer,
        "{label:<20}{annotation_count:>6}{occupied_bases:>13} bp   {percent:>5.2} %",
        percent = percent_of(occupied_bases, total_bases),
    )?;
    Ok(())
}

fn class_is_small_rna(class_name: &str) -> bool {
    matches!(
        class_name,
        "RNA" | "rRNA" | "scRNA" | "snRNA" | "srpRNA" | "tRNA"
    ) || class_name.ends_with("RNA")
        || class_name.contains("/RNA")
}

fn class_is_interspersed_other(class_name: &str) -> bool {
    !class_name.starts_with("SINE")
        && !class_name.starts_with("PLE")
        && !class_name.starts_with("Penelope")
        && !class_name.starts_with("LINE")
        && !class_name.starts_with("LTR")
        && !class_name.starts_with("DNA")
        && !class_name.starts_with("Satellite")
        && !class_name.contains("Simple_repeat")
        && !class_name.contains("Low_complexity")
        && !class_is_small_rna(class_name)
}

fn class_is_interspersed(class_name: &str) -> bool {
    class_name.starts_with("SINE")
        || class_name.starts_with("PLE")
        || class_name.starts_with("Penelope")
        || class_name.starts_with("LINE")
        || class_name.starts_with("LTR")
        || class_name.starts_with("DNA")
        || class_is_interspersed_other(class_name)
}

fn catalog_requires_report_adjudication(catalog: &RepeatMaskerCatalog) -> bool {
    catalog.metadata != RepeatMaskerMetadata::default()
}

fn write_repeatmasker_tbl_with_chain_counter<C>(
    path: &Path,
    source_name: &str,
    annotations: &[RepeatMaskerAnnotation],
    metadata: &RepeatMaskerMetadata,
    mask_stats: Option<&MaskStats>,
    chain_counter: &C,
) -> Result<(), RepeatMaskerError>
where
    C: Fn(&[RepeatMaskerAnnotation]) -> usize,
{
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    let sequence_count = mask_stats
        .map(|stats| stats.sequence_count)
        .or(metadata.total_sequences)
        .unwrap_or_else(|| annotation_stats(annotations).sequence_count);
    let total_length_display = mask_stats
        .map(|stats| stats.total_bases)
        .or(metadata.total_length)
        .or(metadata.total_non_mask_bases)
        .unwrap_or_else(|| annotation_stats(annotations).occupied_bases);
    let percentage_denominator_bases = total_length_display.max(1);
    let non_mask_bases = mask_stats
        .map(|stats| stats.non_ambiguous_bases_excluding_long_runs)
        .or(metadata.total_non_mask_bases)
        .or(Some(total_length_display));
    let non_mask_note = non_mask_bases
        .map(|bases| format!("({bases} bp excl N/X-runs)"))
        .unwrap_or_default();
    let gc_level = mask_stats
        .map(|stats| format!("{:.2}", stats.gc_fraction * 100.0))
        .or_else(|| {
            metadata
                .reported_gc_percent_hundredths
                .map(format_percent_hundredths)
        })
        .unwrap_or_else(|| "Unknown".to_string());
    let masked_bases = mask_stats
        .map(|stats| stats.masked_bases)
        .or(metadata.reported_masked_bases)
        .unwrap_or_else(|| annotation_stats(annotations).occupied_bases);

    let (sines_count, sines_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.starts_with("SINE"),
        chain_counter,
    );
    let (alus_count, alus_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.starts_with("SINE/Alu"),
        chain_counter,
    );
    let (mirs_count, mirs_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.starts_with("SINE/MIR"),
        chain_counter,
    );
    let (penelope_count, penelope_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.starts_with("PLE") || class_name.starts_with("Penelope"),
        chain_counter,
    );
    let (lines_count, lines_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.starts_with("LINE"),
        chain_counter,
    );
    let (line1_count, line1_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.starts_with("LINE/L1"),
        chain_counter,
    );
    let (line2_count, line2_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.starts_with("LINE/L2"),
        chain_counter,
    );
    let (linecr1_count, linecr1_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.starts_with("LINE/CR1"),
        chain_counter,
    );
    let (ltr_count, ltr_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.starts_with("LTR"),
        chain_counter,
    );
    let (malr_count, malr_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.contains("MaLR"),
        chain_counter,
    );
    let (ervl_count, ervl_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.starts_with("LTR/ERVL") && !class_name.contains("MaLR"),
        chain_counter,
    );
    let (erv1_count, erv1_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.starts_with("LTR/ERV1"),
        chain_counter,
    );
    let (ervk_count, ervk_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.starts_with("LTR/ERVK"),
        chain_counter,
    );
    let (dna_count, dna_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.starts_with("DNA"),
        chain_counter,
    );
    let (mer1_count, mer1_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.contains("MER1_type"),
        chain_counter,
    );
    let (mer2_count, mer2_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.contains("MER2_type"),
        chain_counter,
    );
    let (other_count, other_bases) = summarize_class_filter_with_counter(
        annotations,
        class_is_interspersed_other,
        chain_counter,
    );
    let (rna_count, rna_bases) =
        summarize_class_filter_with_counter(annotations, class_is_small_rna, chain_counter);
    let (satellite_count, satellite_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.starts_with("Satellite"),
        chain_counter,
    );
    let (simple_count, simple_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.contains("Simple_repeat"),
        chain_counter,
    );
    let (low_complexity_count, low_complexity_bases) = summarize_class_filter_with_counter(
        annotations,
        |class_name| class_name.contains("Low_complexity"),
        chain_counter,
    );
    let (_, total_interspersed_bases) =
        summarize_class_filter_with_counter(annotations, class_is_interspersed, chain_counter);

    writeln!(writer, "==================================================")?;
    writeln!(writer, "file name: {source_name}")?;
    writeln!(writer, "sequences: {:>13}", sequence_count)?;
    writeln!(
        writer,
        "total length: {:>10} bp  {}",
        total_length_display, non_mask_note
    )?;
    writeln!(writer, "GC level: {:>13} %", gc_level)?;
    writeln!(
        writer,
        "bases masked: {:>10} bp ( {:>5.2} %)",
        masked_bases,
        percent_of(masked_bases, percentage_denominator_bases)
    )?;
    writeln!(writer, "==================================================")?;
    writeln!(writer, "               number of      length   percentage")?;
    writeln!(writer, "               elements*    occupied  of sequence")?;
    writeln!(writer, "--------------------------------------------------")?;
    write_tbl_row(
        &mut writer,
        "SINEs:",
        sines_count,
        sines_bases,
        percentage_denominator_bases,
    )?;
    write_tbl_row(
        &mut writer,
        "      ALUs",
        alus_count,
        alus_bases,
        percentage_denominator_bases,
    )?;
    write_tbl_row(
        &mut writer,
        "      MIRs",
        mirs_count,
        mirs_bases,
        percentage_denominator_bases,
    )?;
    writeln!(writer)?;
    if penelope_count > 0 || penelope_bases > 0 {
        write_tbl_row(
            &mut writer,
            "Penelope",
            penelope_count,
            penelope_bases,
            percentage_denominator_bases,
        )?;
        writeln!(writer)?;
    }
    write_tbl_row(
        &mut writer,
        "LINEs:",
        lines_count,
        lines_bases,
        percentage_denominator_bases,
    )?;
    write_tbl_row(
        &mut writer,
        "      LINE1",
        line1_count,
        line1_bases,
        percentage_denominator_bases,
    )?;
    write_tbl_row(
        &mut writer,
        "      LINE2",
        line2_count,
        line2_bases,
        percentage_denominator_bases,
    )?;
    write_tbl_row(
        &mut writer,
        "      L3/CR1",
        linecr1_count,
        linecr1_bases,
        percentage_denominator_bases,
    )?;
    writeln!(writer)?;
    write_tbl_row(
        &mut writer,
        "LTR elements:",
        ltr_count,
        ltr_bases,
        percentage_denominator_bases,
    )?;
    write_tbl_row(
        &mut writer,
        "      MaLRs",
        malr_count,
        malr_bases,
        percentage_denominator_bases,
    )?;
    write_tbl_row(
        &mut writer,
        "      ERVL",
        ervl_count,
        ervl_bases,
        percentage_denominator_bases,
    )?;
    write_tbl_row(
        &mut writer,
        "      ERV_classI",
        erv1_count,
        erv1_bases,
        percentage_denominator_bases,
    )?;
    write_tbl_row(
        &mut writer,
        "      ERV_classII",
        ervk_count,
        ervk_bases,
        percentage_denominator_bases,
    )?;
    writeln!(writer)?;
    write_tbl_row(
        &mut writer,
        "DNA elements:",
        dna_count,
        dna_bases,
        percentage_denominator_bases,
    )?;
    write_tbl_row(
        &mut writer,
        "      MER1_type",
        mer1_count,
        mer1_bases,
        percentage_denominator_bases,
    )?;
    write_tbl_row(
        &mut writer,
        "      MER2_type",
        mer2_count,
        mer2_bases,
        percentage_denominator_bases,
    )?;
    writeln!(writer)?;
    write_tbl_row(
        &mut writer,
        "Unclassified:",
        other_count,
        other_bases,
        percentage_denominator_bases,
    )?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Total interspersed repeats:{:>10} bp   {:>5.2} %",
        total_interspersed_bases,
        percent_of(total_interspersed_bases, percentage_denominator_bases)
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    write_tbl_row(
        &mut writer,
        "Small RNA:",
        rna_count,
        rna_bases,
        percentage_denominator_bases,
    )?;
    writeln!(writer)?;
    write_tbl_row(
        &mut writer,
        "Satellites:",
        satellite_count,
        satellite_bases,
        percentage_denominator_bases,
    )?;
    write_tbl_row(
        &mut writer,
        "Simple repeats:",
        simple_count,
        simple_bases,
        percentage_denominator_bases,
    )?;
    write_tbl_row(
        &mut writer,
        "Low complexity:",
        low_complexity_count,
        low_complexity_bases,
        percentage_denominator_bases,
    )?;
    writeln!(writer, "==================================================")?;
    writeln!(writer)?;
    writeln!(
        writer,
        "* most repeats fragmented by insertions or deletions"
    )?;
    writeln!(writer, "  have been counted as one element")?;
    writeln!(writer)?;
    if let Some(version) = &metadata.version {
        writeln!(writer, "{version}")?;
    }
    if let Some(engine) = &metadata.engine {
        writeln!(writer, "{engine}")?;
    }
    if let Some(library) = &metadata.library {
        writeln!(writer, "{library}")?;
    }
    writer.flush()?;
    Ok(())
}

pub fn write_repeatmasker_tbl_path(
    path: &Path,
    source_name: &str,
    catalog: &RepeatMaskerCatalog,
    mask_stats: Option<&MaskStats>,
) -> Result<(), RepeatMaskerError> {
    let annotations = if catalog_requires_report_adjudication(catalog) {
        adjudicate_annotations(&catalog.annotations)
    } else {
        catalog.annotations.clone()
    };
    write_repeatmasker_tbl_with_chain_counter(
        path,
        source_name,
        &annotations,
        &catalog.metadata,
        mask_stats,
        &default_element_count,
    )
}

pub fn write_repeatmasker_tbl_with_alignments_path(
    path: &Path,
    source_name: &str,
    catalog: &RepeatMaskerCatalog,
    mask_stats: Option<&MaskStats>,
    alignments: &[RepeatMaskerAlignmentRecord],
) -> Result<(), RepeatMaskerError> {
    let annotations = if catalog_requires_report_adjudication(catalog) {
        adjudicate_annotations_with_alignments(&catalog.annotations, alignments)
    } else {
        catalog.annotations.clone()
    };
    write_repeatmasker_tbl_with_chain_counter(
        path,
        source_name,
        &annotations,
        &catalog.metadata,
        mask_stats,
        &|filtered| default_element_count_with_alignments(filtered, alignments),
    )
}

#[cfg(test)]
mod tests {
    use super::{
        adjudicate_annotations, adjudicate_annotations_with_alignments, adjudicated_catalog_stats,
        build_annotation_chains, build_annotation_chains_with_alignments, catalog_stats,
        compute_repeatmasker_mask_stats_from_fasta, consolidate_mask_intervals,
        encode_repeatmasker_alignment, load_out_annotations, load_repeatmasker_alignments,
        load_repeatmasker_catalog, mask_fasta_records, mask_repeatmasker_annotations_to_fasta,
        mask_repeatmasker_out_to_fasta, parse_cat_annotation_line, parse_out_annotation_line,
        write_alignment_tsv_path, write_annotations_tsv_path, write_chain_tsv_path,
        write_chain_tsv_with_alignments_path, write_gff3_path, write_out_annotations_path,
        write_repeat_summary_tsv_path, write_repeatmasker_tbl_path,
        write_repeatmasker_tbl_with_alignments_path, FastaRecord, MaskInterval, MaskMode,
    };
    use std::collections::BTreeMap;
    use std::fs;
    use std::io::Cursor;
    use tempfile::TempDir;

    #[test]
    fn parse_out_annotation_line_skips_headers() {
        assert!(parse_out_annotation_line("   SW  perc perc", 1)
            .unwrap()
            .is_none());
        assert!(parse_out_annotation_line("", 2).unwrap().is_none());
    }

    #[test]
    fn parse_out_annotation_line_parses_forward_and_reverse_annotations() {
        let forward = parse_out_annotation_line(
            "  307  6.2  0.0  0.0  chr1  10  20 (80) + AluJo SINE/Alu 1 11 (0) 6 m_b1s1i0 *",
            1,
        )
        .unwrap()
        .unwrap();
        assert_eq!(forward.query_name, "chr1");
        assert_eq!(forward.query_begin, 10);
        assert_eq!(forward.query_end, 20);
        assert_eq!(forward.orientation, '+');
        assert_eq!(forward.subject_class, "SINE/Alu");
        assert_eq!(forward.subject_begin, 1);
        assert_eq!(forward.subject_end, 11);
        assert_eq!(forward.subject_left, 0);

        let reverse = parse_out_annotation_line(
            "  512 12.3  1.4  0.5  chr2 100 140 (860) C L1 LINE/L1 (10) 220 180 12 m_b2s4i9",
            2,
        )
        .unwrap()
        .unwrap();
        assert_eq!(reverse.query_name, "chr2");
        assert_eq!(reverse.orientation, 'C');
        assert_eq!(reverse.query_left, 860);
        assert_eq!(reverse.subject_begin, 180);
        assert_eq!(reverse.subject_end, 220);
        assert_eq!(reverse.subject_left, 10);
        assert_eq!(reverse.annotation_id, "12");
        assert_eq!(reverse.lineage_id, "m_b2s4i9");
    }

    #[test]
    fn parse_out_annotation_line_treats_trailing_star_as_overlap_marker() {
        let annotation = parse_out_annotation_line(
            "  307  6.2  0.0  0.0  chr1  10  20 (80) + AluJo SINE/Alu 1 11 (0) 6 *",
            1,
        )
        .unwrap()
        .unwrap();
        assert_eq!(annotation.annotation_id, "6");
        assert!(annotation.lineage_id.is_empty());
        assert_eq!(annotation.overlap_marker.as_deref(), Some("*"));
    }

    #[test]
    fn parse_cat_annotation_line_parses_forward_and_reverse_annotations() {
        let forward = parse_cat_annotation_line(
            " 2134 10.70 3.30 1.32 seq1 1809 2111 (21) AluSx#SINE/Alu 3 311 (1) 1",
            1,
        )
        .unwrap()
        .unwrap();
        assert_eq!(forward.query_name, "seq1");
        assert_eq!(forward.orientation, '+');
        assert_eq!(forward.subject_name, "AluSx");
        assert_eq!(forward.subject_class, "SINE/Alu");
        assert_eq!(forward.subject_begin, 3);
        assert_eq!(forward.subject_end, 311);
        assert_eq!(forward.subject_left, 1);
        assert_eq!(forward.annotation_id, "1");

        let reverse = parse_cat_annotation_line(
            " 333 30.39 10.80 5.56 seq2 709 1032 (24) C L1MC4a_3end#LINE/L1 (2103) 452 112 m_b1s2 5",
            2,
        )
        .unwrap()
        .unwrap();
        assert_eq!(reverse.orientation, 'C');
        assert_eq!(reverse.subject_name, "L1MC4a_3end");
        assert_eq!(reverse.subject_class, "LINE/L1");
        assert_eq!(reverse.subject_begin, 112);
        assert_eq!(reverse.subject_end, 452);
        assert_eq!(reverse.subject_left, 2103);
        assert_eq!(reverse.lineage_id, "m_b1s2");
        assert_eq!(reverse.annotation_id, "5");
    }

    #[test]
    fn consolidate_mask_intervals_merges_overlaps_and_adjacency() {
        let annotations = load_out_annotations(Cursor::new(
            "\
  300  1.0  0.0  0.0 seq1  5 10 (90) + Alu SINE/Alu 1 6 (0) 1 a\n\
  290  1.0  0.0  0.0 seq1  8 12 (88) + Alu SINE/Alu 7 11 (0) 2 a\n\
  280  1.0  0.0  0.0 seq1 13 15 (85) + Alu SINE/Alu 12 14 (0) 3 a\n\
  270  1.0  0.0  0.0 seq2  1  2 (10) + Alu SINE/Alu 1 2 (0) 4 a\n",
        ))
        .unwrap();
        let intervals = consolidate_mask_intervals(&annotations);
        assert_eq!(
            intervals.get("seq1").unwrap(),
            &vec![MaskInterval { begin: 5, end: 15 }]
        );
        assert_eq!(
            intervals.get("seq2").unwrap(),
            &vec![MaskInterval { begin: 1, end: 2 }]
        );
    }

    #[test]
    fn mask_fasta_records_supports_n_x_and_lowercase_modes() {
        let mut records = vec![FastaRecord {
            id: "seq1".to_string(),
            description: String::new(),
            sequence: b"ACGTACGT".to_vec(),
        }];
        let intervals =
            BTreeMap::from([("seq1".to_string(), vec![MaskInterval { begin: 3, end: 6 }])]);

        let stats = mask_fasta_records(&mut records, &intervals, MaskMode::N);
        assert_eq!(records[0].sequence, b"ACNNNNGT");
        assert_eq!(stats.masked_bases, 4);

        let mut lower_records = vec![FastaRecord {
            id: "seq1".to_string(),
            description: String::new(),
            sequence: b"ACGTACGT".to_vec(),
        }];
        mask_fasta_records(&mut lower_records, &intervals, MaskMode::Lowercase);
        assert_eq!(lower_records[0].sequence, b"ACgtacGT");

        let mut x_records = vec![FastaRecord {
            id: "seq1".to_string(),
            description: String::new(),
            sequence: b"ACGTACGT".to_vec(),
        }];
        mask_fasta_records(&mut x_records, &intervals, MaskMode::X);
        assert_eq!(x_records[0].sequence, b"ACXXXXGT");
    }

    #[test]
    fn mask_repeatmasker_out_to_fasta_masks_real_files_and_reports_stats() {
        let temp_dir = TempDir::new().unwrap();
        let out_path = temp_dir.path().join("sample.out");
        let fasta_path = temp_dir.path().join("sample.fa");
        let masked_path = temp_dir.path().join("sample.fa.masked");

        fs::write(
            &out_path,
            "\
header line\n\
  450  2.0  0.0  0.0 seq1  2  4 (4) + Alu SINE/Alu 1 3 (0) 1 m_b1s1i0\n\
  430  2.0  0.0  0.0 seq1  4  6 (2) + Alu SINE/Alu 3 5 (0) 2 m_b1s1i1\n\
  410  2.0  0.0  0.0 seq2  1  3 (3) C L1 LINE/L1 (1) 10 8 3 m_b2s1i0\n",
        )
        .unwrap();
        fs::write(&fasta_path, ">seq1 first\nACGTACGT\n>seq2\nGGGGTT\n").unwrap();

        let stats = mask_repeatmasker_out_to_fasta(
            &out_path,
            &fasta_path,
            &masked_path,
            MaskMode::Lowercase,
            50,
        )
        .unwrap();

        let masked = fs::read_to_string(masked_path).unwrap();
        assert!(masked.contains(">seq1 first\nAcgtacGT\n"));
        assert!(masked.contains(">seq2\ngggGTT\n"));
        assert_eq!(stats.sequence_count, 2);
        assert_eq!(stats.annotation_count, 3);
        assert_eq!(stats.masked_bases, 8);
    }

    #[test]
    fn load_repeatmasker_catalog_parses_cat_metadata_and_summary_rows() {
        let catalog = load_repeatmasker_catalog(Cursor::new(
            "\
 2134 10.70 3.30 1.32 seq1 1809 2111 (21) AluSx#SINE/Alu 3 311 (1) 1\n\
 333 30.39 10.80 5.56 seq2 709 1032 (24) C L1MC4a_3end#LINE/L1 (2103) 452 112 5\n\
## Batch Overlap Boundaries\n\
##   seq1  10, 20, 30, 40\n\
## Total Sequences: 2\n\
## Total Length: 4266\n\
## Total NonMask ( excluding >20bp runs of N/X bases ): 3490\n\
## Total NonSub ( excluding all non ACGT bases ):4266\n\
RepeatMasker version 4.1.9, default mode\n\
run with cross_match version 0.990329\n\
RM Library: CONS-Dfam_withRBRM_3.9\n",
        ))
        .unwrap();

        assert_eq!(catalog.annotations.len(), 2);
        assert_eq!(catalog.metadata.total_sequences, Some(2));
        assert_eq!(catalog.metadata.total_length, Some(4266));
        assert_eq!(catalog.metadata.total_non_mask_bases, Some(3490));
        assert_eq!(catalog.metadata.total_non_sub_bases, Some(4266));
        assert_eq!(
            catalog
                .metadata
                .batch_overlap_boundaries
                .get("seq1")
                .unwrap(),
            &vec![10, 20, 30, 40]
        );

        let summary = catalog_stats(&catalog);
        assert_eq!(summary.annotation_stats.annotation_count, 2);
        assert_eq!(summary.annotation_stats.occupied_bases, 627);
        assert_eq!(summary.repeat_class_summary[0].key, "LINE/L1");
        assert_eq!(summary.repeat_class_summary[0].aligned_bases, 324);
        assert_eq!(summary.repeat_class_summary[0].occupied_bases, 324);
        assert_eq!(summary.repeat_family_summary[1].key, "AluSx");
    }

    #[test]
    fn load_repeatmasker_catalog_parses_legacy_cat_footer_lines() {
        let catalog = load_repeatmasker_catalog(Cursor::new(
            "\
 2134 10.70 3.30 1.32 seq1 1809 2111 (21) AluSx#SINE/Alu 3 311 (1) 1\n\
2 sequences with 4266 bp (3490 bp not contained in strings of more than 20 Ns)\n\
1254 bp masked\n\
GC level 47.87 % (46.64 % after masking)\n\
RepeatMasker version development, default mode\n\
run with cross_match version 0.990329\n\
RepBase Update 8.12, RM database version 20040306\n",
        ))
        .unwrap();

        assert_eq!(catalog.metadata.total_sequences, Some(2));
        assert_eq!(catalog.metadata.total_length, Some(4266));
        assert_eq!(catalog.metadata.total_non_mask_bases, Some(3490));
        assert_eq!(catalog.metadata.reported_masked_bases, Some(1254));
        assert_eq!(catalog.metadata.reported_gc_percent_hundredths, Some(4787));
        assert_eq!(
            catalog.metadata.reported_masked_gc_percent_hundredths,
            Some(4664)
        );
        assert_eq!(
            catalog.metadata.library.as_deref(),
            Some("RepBase Update 8.12, RM database version 20040306")
        );
    }

    #[test]
    fn load_repeatmasker_alignments_parses_synthetic_record_and_metrics() {
        let alignments = load_repeatmasker_alignments(Cursor::new(
            "\
  500  1.0  0.0  0.0 seq1 10 21 (79) AluY 1 12 (0) m_b1s1i0 2\n\
\n\
  seq1          10 ACGTACGTACGT 21\n\
                           i   v      \n\
  AluY           1 ACGTTCGTAC-T 11\n\
\n\
Matrix = simple.matrix\n\
Kimura (with divCpGMod) = 2.50\n\
CpG sites = 3, Kimura (unadjusted) = 2.10\n\
Transitions / transversions = 2.00 (4/2)\n\
Gap_init rate = 0.10 (1 / 10), avg. gap size = 1.00 (1 / 1)\n\
",
        ))
        .unwrap();

        assert_eq!(alignments.len(), 1);
        let alignment = &alignments[0];
        assert_eq!(alignment.query_name, "seq1");
        assert_eq!(alignment.subject_name, "AluY");
        assert_eq!(alignment.annotation_id, "2");
        assert_eq!(alignment.lineage_id, "m_b1s1i0");
        assert_eq!(alignment.query_sequence, "ACGTACGTACGT");
        assert_eq!(alignment.subject_sequence, "ACGTTCGTAC-T");
        assert_eq!(alignment.metrics.aligned_query_bases, 12);
        assert_eq!(alignment.metrics.aligned_subject_bases, 11);
        assert_eq!(
            alignment.metrics.matrix_name.as_deref(),
            Some("simple.matrix")
        );
        assert_eq!(alignment.metrics.cpg_sites, Some(3));
        assert_eq!(alignment.metrics.transition_count, Some(4));
        assert_eq!(alignment.metrics.transversion_count, Some(2));
        assert_eq!(alignment.metrics.gap_count, Some(1));
        assert_eq!(alignment.metrics.average_gap_size, Some(1.0));
        assert_eq!(alignment.metrics.kimura_divergence, Some(2.5));
        assert_eq!(alignment.metrics.raw_kimura_divergence, Some(2.1));
    }

    #[test]
    fn load_repeatmasker_alignments_reverse_complements_query_oriented_records() {
        let alignments = load_repeatmasker_alignments(Cursor::new(
            "\
  400  1.0  0.0  0.0 seq2 30 33 (10) C L1MC4 (20) 8 5 m_b1s2 7\n\
\n\
C seq2          33 ACGA 30\n\
                           \n\
  L1MC4          8 TCGT 5\n\
\n\
Gap_init rate = 0.00 (0 / 3), avg. gap size = 0.00 (0 / 0)\n\
",
        ))
        .unwrap();

        assert_eq!(alignments.len(), 1);
        let alignment = &alignments[0];
        assert_eq!(alignment.orientation, 'C');
        assert_eq!(alignment.query_sequence, "TCGT");
        assert_eq!(alignment.subject_sequence, "ACGA");
    }

    #[test]
    fn adjudicate_annotations_suppresses_contained_lower_score_duplicates_and_short_fragments() {
        let annotations = load_out_annotations(Cursor::new(
            "\
  450  2.0  0.0  0.0 seq1 10 40 (60) + Alu SINE/Alu 1 31 (0) 1 a\n\
  440  2.0  0.0  0.0 seq1 15 30 (70) + Alu SINE/Alu 6 21 (0) 2 a\n\
  430  2.0  0.0  0.0 seq1 50 57 (43) + L1 LINE/L1 1 8 (0) 3 a\n\
  420  2.0  0.0  0.0 seq1 70 75 (25) + (CA)n Simple_repeat 1 6 (0) 4 a\n",
        ))
        .unwrap();

        let adjudicated = adjudicate_annotations(&annotations);
        assert_eq!(adjudicated.len(), 2);
        assert_eq!(adjudicated[0].annotation_id, "1");
        assert_eq!(adjudicated[1].annotation_id, "4");
    }

    #[test]
    fn adjudicated_catalog_stats_record_suppressed_counts() {
        let catalog = load_repeatmasker_catalog(Cursor::new(
            "\
  450  2.0  0.0  0.0 seq1 10 40 (60) + Alu SINE/Alu 1 31 (0) 1 a\n\
  440  2.0  0.0  0.0 seq1 15 30 (70) + Alu SINE/Alu 6 21 (0) 2 a\n",
        ))
        .unwrap();

        let summary = adjudicated_catalog_stats(&catalog);
        assert_eq!(summary.annotation_stats.annotation_count, 1);
        assert_eq!(summary.raw_annotation_count, Some(2));
        assert_eq!(summary.suppressed_annotation_count, Some(1));
    }

    #[test]
    fn build_annotation_chains_merges_exact_family_fragments() {
        let annotations = load_out_annotations(Cursor::new(
            "\
  450  2.0  0.0  0.0 seq1 10 40 (60) + AluY SINE/Alu 1 31 (0) 1 a\n\
  440  2.5  0.0  0.0 seq1 44 72 (28) + AluY SINE/Alu 32 60 (0) 2 a\n\
  430  2.0  0.0  0.0 seq1 90 120 (0) + AluSx SINE/Alu 1 31 (0) 3 a\n",
        ))
        .unwrap();

        let chains = build_annotation_chains(&annotations);
        assert_eq!(chains.len(), 2);
        assert_eq!(chains[0].subject_name, "AluY");
        assert_eq!(chains[0].annotation_count, 2);
        assert_eq!(chains[0].query_begin, 10);
        assert_eq!(chains[0].query_end, 72);
        assert_eq!(chains[0].annotation_ids, vec!["1", "2"]);
        assert_eq!(chains[1].annotation_ids, vec!["3"]);
    }

    #[test]
    fn alignments_refine_fragment_chain_boundaries() {
        let annotations = load_out_annotations(Cursor::new(
            "\
  500  1.0  0.0  0.0 seq1 10 19 (81) + AluY SINE/Alu 1 50 (0) 1 m_b1s1i0\n\
  490  1.0  0.0  0.0 seq1 20 59 (41) + AluY SINE/Alu 18 60 (0) 2 m_b1s1i0\n\
",
        ))
        .unwrap();
        let alignments = load_repeatmasker_alignments(Cursor::new(
            "\
  500  1.0  0.0  0.0 seq1 10 19 (81) AluY 1 50 (0) m_b1s1i0 1\n\
\n\
  seq1          10 ACGTACGTAA---------------------------------------- 19\n\
                                                                   \n\
  AluY           1 ACGTACGTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 50\n\
\n\
Gap_init rate = 0.50 (1 / 2), avg. gap size = 41.00 (41 / 1)\n\
  490  1.0  0.0  0.0 seq1 20 59 (41) AluY 21 60 (0) m_b1s1i0 2\n\
\n\
  seq1          20 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT 59\n\
                                                                   \n\
  AluY          21 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT 60\n\
\n\
Gap_init rate = 0.00 (0 / 39), avg. gap size = 0.00 (0 / 0)\n\
",
        ))
        .unwrap();

        let plain_chains = build_annotation_chains(&annotations);
        let refined_chains = build_annotation_chains_with_alignments(&annotations, &alignments);
        assert_eq!(plain_chains.len(), 2);
        assert_eq!(refined_chains.len(), 1);
        assert_eq!(refined_chains[0].annotation_ids, vec!["1", "2"]);
    }

    #[test]
    fn chain_counts_flow_into_catalog_stats_and_tbl() {
        let temp_dir = TempDir::new().unwrap();
        let tbl_path = temp_dir.path().join("chained.tbl");
        let catalog = load_repeatmasker_catalog(Cursor::new(
            "\
  450  2.0  0.0  0.0 sample.fa 10 40 (60) + AluY SINE/Alu 1 31 (0) 1 a\n\
  440  2.5  0.0  0.0 sample.fa 44 72 (28) + AluY SINE/Alu 32 60 (0) 2 a\n\
1 sequences with 200 bp (200 bp not contained in strings of more than 20 Ns)\n\
63 bp masked\n\
GC level 50.00 % (50.00 % after masking)\n\
",
        ))
        .unwrap();

        let summary = catalog_stats(&catalog);
        assert_eq!(summary.annotation_stats.annotation_count, 2);
        assert_eq!(summary.element_chain_count, 1);

        write_repeatmasker_tbl_path(&tbl_path, "sample.fa", &catalog, None).unwrap();
        let tbl_text = fs::read_to_string(tbl_path).unwrap();
        assert!(tbl_text.contains("SINEs:"));
        assert!(tbl_text.contains("     1           60 bp"));
        assert!(!tbl_text.contains("     2           60 bp"));
    }

    #[test]
    fn out_tables_count_global_annotation_ids_across_split_query_sequences() {
        let temp_dir = TempDir::new().unwrap();
        let tbl_path = temp_dir.path().join("split.tbl");
        let catalog = load_repeatmasker_catalog(Cursor::new(
            "\
  450  2.0  0.0  0.0 part1 10 20 (80) + L2a LINE/L2 1 11 (0) 1 *\n\
  440  2.0  0.0  0.0 part2 30 40 (60) + L2a LINE/L2 12 22 (0) 1\n\
",
        ))
        .unwrap();

        write_repeatmasker_tbl_path(&tbl_path, "split.fa", &catalog, None).unwrap();
        let tbl_text = fs::read_to_string(tbl_path).unwrap();
        assert!(tbl_text.contains("LINEs:"));
        assert!(tbl_text.contains("     1           22 bp"));
        assert!(!tbl_text.contains("     2           22 bp"));
    }

    #[test]
    fn alignments_flow_into_chain_counts_and_tbl() {
        let temp_dir = TempDir::new().unwrap();
        let plain_tbl_path = temp_dir.path().join("plain.tbl");
        let refined_tbl_path = temp_dir.path().join("refined.tbl");
        let catalog = load_repeatmasker_catalog(Cursor::new(
            "\
  500  1.0  0.0  0.0 sample.fa 10 19 (81) + AluY SINE/Alu 1 50 (0) 1 m_b1s1i0\n\
  490  1.0  0.0  0.0 sample.fa 20 59 (41) + AluY SINE/Alu 18 60 (0) 2 m_b1s1i0\n\
1 sequences with 200 bp (200 bp not contained in strings of more than 20 Ns)\n\
50 bp masked\n\
GC level 50.00 % (50.00 % after masking)\n\
",
        ))
        .unwrap();
        let alignments = load_repeatmasker_alignments(Cursor::new(
            "\
  500  1.0  0.0  0.0 sample.fa 10 19 (81) AluY 1 50 (0) m_b1s1i0 1\n\
\n\
  sample.fa     10 ACGTACGTAA---------------------------------------- 19\n\
                                                                   \n\
  AluY           1 ACGTACGTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 50\n\
\n\
Gap_init rate = 0.50 (1 / 2), avg. gap size = 41.00 (41 / 1)\n\
  490  1.0  0.0  0.0 sample.fa 20 59 (41) AluY 21 60 (0) m_b1s1i0 2\n\
\n\
  sample.fa     20 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT 59\n\
                                                                   \n\
  AluY          21 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT 60\n\
\n\
Gap_init rate = 0.00 (0 / 39), avg. gap size = 0.00 (0 / 0)\n\
",
        ))
        .unwrap();

        write_repeatmasker_tbl_path(&plain_tbl_path, "sample.fa", &catalog, None).unwrap();
        write_repeatmasker_tbl_with_alignments_path(
            &refined_tbl_path,
            "sample.fa",
            &catalog,
            None,
            &alignments,
        )
        .unwrap();

        let plain_tbl_text = fs::read_to_string(plain_tbl_path).unwrap();
        let refined_tbl_text = fs::read_to_string(refined_tbl_path).unwrap();
        assert!(plain_tbl_text.contains("SINEs:"));
        assert!(plain_tbl_text.contains("     2           50 bp"));
        assert!(refined_tbl_text.contains("     1           50 bp"));
    }

    #[test]
    fn alignments_break_ties_for_identical_family_hits() {
        let catalog = load_repeatmasker_catalog(Cursor::new(
            "\
  500  1.0  0.0  0.0 seq1 10 21 (79) AluY#SINE/Alu 1 12 (0) m_b1s1i0 1\n\
  500  1.0  0.0  0.0 seq1 10 21 (79) AluY#SINE/Alu 1 12 (0) m_b1s1i0 2\n\
",
        ))
        .unwrap();
        let alignments = load_repeatmasker_alignments(Cursor::new(
            "\
  500  1.0  0.0  0.0 seq1 10 21 (79) AluY 1 12 (0) m_b1s1i0 1\n\
\n\
  seq1          10 ACGTACGTACGT 21\n\
                           i        \n\
  AluY           1 ACGTTCGTACGT 12\n\
\n\
Kimura (with divCpGMod) = 8.00\n\
Gap_init rate = 0.00 (0 / 11), avg. gap size = 0.00 (0 / 0)\n\
  500  1.0  0.0  0.0 seq1 10 21 (79) AluY 1 12 (0) m_b1s1i0 2\n\
\n\
  seq1          10 ACGTACGTACGT 21\n\
                                   \n\
  AluY           1 ACGTACGTACGT 12\n\
\n\
Kimura (with divCpGMod) = 2.00\n\
Gap_init rate = 0.00 (0 / 11), avg. gap size = 0.00 (0 / 0)\n\
",
        ))
        .unwrap();

        let adjudicated = adjudicate_annotations_with_alignments(&catalog.annotations, &alignments);
        assert_eq!(adjudicated.len(), 1);
        assert_eq!(adjudicated[0].annotation_id, "2");
    }

    #[test]
    fn write_chain_tsv_exports_fragment_chains() {
        let temp_dir = TempDir::new().unwrap();
        let chain_path = temp_dir.path().join("sample.chains.tsv");
        let annotations = load_out_annotations(Cursor::new(
            "\
  450  2.0  0.0  0.0 seq1 10 40 (60) + AluY SINE/Alu 1 31 (0) 1 a\n\
  440  2.5  0.0  0.0 seq1 44 72 (28) + AluY SINE/Alu 32 60 (0) 2 a\n",
        ))
        .unwrap();

        write_chain_tsv_path(&chain_path, &annotations).unwrap();
        let text = fs::read_to_string(chain_path).unwrap();
        assert!(text.contains("chain_id\tquery_name\torientation"));
        assert!(text.contains("\tseq1\t+\tAluY\tSINE/Alu\t2\t10\t72\t1\t60\t1;2\t"));
    }

    #[test]
    fn write_chain_tsv_with_alignments_exports_refined_fragment_chains() {
        let temp_dir = TempDir::new().unwrap();
        let chain_path = temp_dir.path().join("sample.refined.chains.tsv");
        let annotations = load_out_annotations(Cursor::new(
            "\
  500  1.0  0.0  0.0 seq1 10 19 (81) + AluY SINE/Alu 1 50 (0) 1 m_b1s1i0\n\
  490  1.0  0.0  0.0 seq1 20 59 (41) + AluY SINE/Alu 18 60 (0) 2 m_b1s1i0\n\
",
        ))
        .unwrap();
        let alignments = load_repeatmasker_alignments(Cursor::new(
            "\
  500  1.0  0.0  0.0 seq1 10 19 (81) AluY 1 50 (0) m_b1s1i0 1\n\
\n\
  seq1          10 ACGTACGTAA---------------------------------------- 19\n\
                                                                   \n\
  AluY           1 ACGTACGTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 50\n\
\n\
Gap_init rate = 0.50 (1 / 2), avg. gap size = 41.00 (41 / 1)\n\
  490  1.0  0.0  0.0 seq1 20 59 (41) AluY 21 60 (0) m_b1s1i0 2\n\
\n\
  seq1          20 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT 59\n\
                                                                   \n\
  AluY          21 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT 60\n\
\n\
Gap_init rate = 0.00 (0 / 39), avg. gap size = 0.00 (0 / 0)\n\
",
        ))
        .unwrap();

        write_chain_tsv_with_alignments_path(&chain_path, &annotations, &alignments).unwrap();
        let text = fs::read_to_string(chain_path).unwrap();
        assert!(text.contains("chain_id\tquery_name\torientation"));
        assert!(text.contains("\tseq1\t+\tAluY\tSINE/Alu\t2\t10\t59\t1\t60\t1;2\t"));
    }

    #[test]
    fn encode_repeatmasker_alignment_emits_cigar_caf_and_blocks() {
        let alignment = load_repeatmasker_alignments(Cursor::new(
            "\
  500  1.0  0.0  0.0 seq1 10 16 (79) AluY 1 7 (0) m_b1s1i0 2\n\
\n\
  seq1          10 AAGCTA--A 16\n\
                           \n\
  AluY           1 AA--TAGGA 7\n\
\n\
Gap_init rate = 0.20 (2 / 10), avg. gap size = 1.50 (3 / 2)\n\
",
        ))
        .unwrap()
        .pop()
        .unwrap();

        let encoding = encode_repeatmasker_alignment(&alignment);
        assert_eq!(encoding.cigar, "2M2D2M2I1M");
        assert_eq!(encoding.caf, "AA-GC-TA+GG+A");
        assert_eq!(encoding.match_bases, 5);
        assert_eq!(encoding.mismatch_bases, 0);
        assert_eq!(encoding.insertion_bases, 2);
        assert_eq!(encoding.deletion_bases, 2);
        assert_eq!(encoding.blocks.len(), 3);
        assert_eq!(
            encoding.blocks[0],
            super::RepeatMaskerAlignmentBlock {
                query_begin: 10,
                query_end: 11,
                subject_begin: 1,
                subject_end: 2,
            }
        );
        assert_eq!(
            encoding.blocks[1],
            super::RepeatMaskerAlignmentBlock {
                query_begin: 14,
                query_end: 15,
                subject_begin: 3,
                subject_end: 4,
            }
        );
        assert_eq!(
            encoding.blocks[2],
            super::RepeatMaskerAlignmentBlock {
                query_begin: 16,
                query_end: 16,
                subject_begin: 7,
                subject_end: 7,
            }
        );
    }

    #[test]
    fn write_alignment_tsv_exports_alignment_encodings() {
        let temp_dir = TempDir::new().unwrap();
        let tsv_path = temp_dir.path().join("sample.align.tsv");
        let alignments = load_repeatmasker_alignments(Cursor::new(
            "\
  500  1.0  0.0  0.0 seq1 10 16 (79) AluY 1 7 (0) m_b1s1i0 2\n\
\n\
  seq1          10 AAGCTA--A 16\n\
                           \n\
  AluY           1 AA--TAGGA 7\n\
\n\
Kimura (with divCpGMod) = 2.00\n\
Gap_init rate = 0.20 (2 / 10), avg. gap size = 1.50 (3 / 2)\n\
",
        ))
        .unwrap();

        write_alignment_tsv_path(&tsv_path, &alignments).unwrap();
        let text = fs::read_to_string(tsv_path).unwrap();
        assert!(text.contains("annotation_id\tlineage_id\tquery_name"));
        assert!(text.contains("\t2M2D2M2I1M\t"));
        assert!(text.contains("\tAA-GC-TA+GG+A\t"));
        assert!(text.contains("\t10-11;14-15;16-16\t"));
        assert!(text.contains("\t1-2;3-4;7-7\t"));
    }

    #[test]
    fn write_repeatmasker_tbl_reports_core_repeatmasker_sections() {
        let temp_dir = TempDir::new().unwrap();
        let tbl_path = temp_dir.path().join("sample.tbl");
        let catalog = load_repeatmasker_catalog(Cursor::new(
            "\
 2134 10.70 3.30 1.32 sample.fa 10 312 (688) AluSx#SINE/Alu 3 311 (1) 1\n\
 586 29.17 9.75 4.69 sample.fa 400 676 (324) L1ME4a_3end#LINE/L1 164 454 (414) 5\n\
 297 0.00 0.00 0.00 sample.fa 700 732 (268) (CAAAA)n#Simple_repeat 3 35 (145) 7\n\
2 sequences with 1000 bp (900 bp not contained in strings of more than 20 Ns)\n\
420 bp masked\n\
GC level 47.87 % (46.64 % after masking)\n\
RepeatMasker version development, default mode\n\
run with cross_match version 0.990329\n\
RepBase Update 8.12, RM database version 20040306\n",
        ))
        .unwrap();

        write_repeatmasker_tbl_path(&tbl_path, "sample.fa", &catalog, None).unwrap();
        let tbl_text = fs::read_to_string(tbl_path).unwrap();
        assert!(tbl_text.contains("file name: sample.fa"));
        assert!(tbl_text.contains("GC level:         47.87 %"));
        assert!(tbl_text.contains("bases masked:        420 bp"));
        assert!(tbl_text.contains("SINEs:"));
        assert!(tbl_text.contains("LINEs:"));
        assert!(tbl_text.contains("Simple repeats:"));
        assert!(tbl_text.contains("Total interspersed repeats:"));
        assert!(tbl_text.contains("RepBase Update 8.12, RM database version 20040306"));
    }

    #[test]
    fn summary_and_tbl_outputs_use_adjudicated_annotations() {
        let temp_dir = TempDir::new().unwrap();
        let summary_path = temp_dir.path().join("sample.summary.tsv");
        let tbl_path = temp_dir.path().join("sample.tbl");
        let catalog = load_repeatmasker_catalog(Cursor::new(
            "\
  450  2.0  0.0  0.0 sample.fa 10 40 (60) + Alu SINE/Alu 1 31 (0) 1 a\n\
  440  2.0  0.0  0.0 sample.fa 15 30 (70) + Alu SINE/Alu 6 21 (0) 2 a\n\
1 sequences with 100 bp (100 bp not contained in strings of more than 20 Ns)\n\
31 bp masked\n\
GC level 50.00 % (50.00 % after masking)\n\
",
        ))
        .unwrap();

        write_repeat_summary_tsv_path(&summary_path, &catalog).unwrap();
        write_repeatmasker_tbl_path(&tbl_path, "sample.fa", &catalog, None).unwrap();

        let summary_text = fs::read_to_string(summary_path).unwrap();
        assert!(summary_text.contains("class\tSINE/Alu\t1\t31\t31"));
        assert!(summary_text.contains("family\tAlu\t1\t31\t31"));

        let tbl_text = fs::read_to_string(tbl_path).unwrap();
        assert!(tbl_text.contains("SINEs:"));
        assert!(tbl_text.contains("     1           31 bp"));
        assert!(!tbl_text.contains("     2           31 bp"));
    }

    #[test]
    fn fasta_driven_tbl_stats_use_total_length_and_exact_gc_without_mask_output() {
        let temp_dir = TempDir::new().unwrap();
        let fasta_path = temp_dir.path().join("sample.fa");
        let tbl_path = temp_dir.path().join("sample.tbl");
        fs::write(&fasta_path, ">sample.fa\nGGGGGNNNNNNNNNNNNNNNNNNNNAAAAA\n").unwrap();

        let catalog = load_repeatmasker_catalog(Cursor::new(
            "\
 2134 10.70 3.30 1.32 sample.fa 1 5 (25) AluSx#SINE/Alu 3 7 (1) 1\n\
2 sequences with 1000 bp (900 bp not contained in strings of more than 20 Ns)\n\
420 bp masked\n\
GC level 47.87 % (46.64 % after masking)\n\
RepeatMasker version development, default mode\n\
run with cross_match version 0.990329\n",
        ))
        .unwrap();

        let stats = compute_repeatmasker_mask_stats_from_fasta(
            &catalog.annotations,
            &fasta_path,
            MaskMode::N,
            None,
        )
        .unwrap();
        assert_eq!(stats.total_bases, 30);
        assert_eq!(stats.non_ambiguous_bases_excluding_long_runs, 10);
        assert_eq!(stats.masked_bases, 5);
        assert!((stats.gc_fraction - (5.0 / 30.0)).abs() < 1e-9);

        write_repeatmasker_tbl_path(&tbl_path, "sample.fa", &catalog, Some(&stats)).unwrap();
        let tbl_text = fs::read_to_string(tbl_path).unwrap();
        assert!(tbl_text.contains("total length:         30 bp  (10 bp excl N/X-runs)"));
        assert!(tbl_text.contains("GC level:         16.67 %"));
        assert!(tbl_text.contains("bases masked:          5 bp ( 16.67 %)"));
    }

    #[test]
    fn write_gff_and_out_exports_match_expected_core_fields() {
        let temp_dir = TempDir::new().unwrap();
        let out_path = temp_dir.path().join("sample.out");
        let gff_path = temp_dir.path().join("sample.gff3");
        let tsv_path = temp_dir.path().join("sample.tsv");
        let summary_path = temp_dir.path().join("sample.summary.tsv");

        let annotations = load_out_annotations(Cursor::new(
            "\
  450  2.0  0.0  0.0 seq1  2 14 (4) + Alu SINE/Alu 1 13 (0) 1 m_b1s1i0\n\
  410  2.0  0.0  0.0 seq2  1 13 (3) C L1 LINE/L1 (1) 20 8 3 m_b2s1i0\n",
        ))
        .unwrap();

        write_out_annotations_path(&out_path, &annotations).unwrap();
        write_gff3_path(&gff_path, &annotations).unwrap();
        write_annotations_tsv_path(&tsv_path, &annotations).unwrap();
        write_repeat_summary_tsv_path(
            &summary_path,
            &super::RepeatMaskerCatalog {
                annotations: annotations.clone(),
                metadata: Default::default(),
            },
        )
        .unwrap();

        let out_text = fs::read_to_string(out_path).unwrap();
        assert!(out_text.contains("seq1"));
        assert!(out_text.contains("SINE/Alu"));
        assert!(out_text.contains("C L1"));

        let gff_text = fs::read_to_string(gff_path).unwrap();
        assert!(gff_text.contains("##gff-version 3"));
        assert!(gff_text.contains("##sequence-region seq1 1 18"));
        assert!(gff_text.contains("RepeatMasker-Rust\tdispersed_repeat\t2\t14\t2.0\t+"));
        assert!(gff_text.contains("Target=Motif:Alu 1 13"));
        assert!(gff_text.contains("Target=Motif:L1 8 20"));

        let tsv_text = fs::read_to_string(tsv_path).unwrap();
        assert!(tsv_text.contains("query_name\tquery_begin"));
        assert!(tsv_text.contains("seq2\t1\t13\t3\t410\t2.0\t0.0\t0.0\t-"));

        let summary_text = fs::read_to_string(summary_path).unwrap();
        assert!(
            summary_text.contains("scope\tkey\tannotation_count\taligned_bases\toccupied_bases")
        );
        assert!(summary_text.contains("class\tSINE/Alu\t1\t13\t13"));
        assert!(summary_text.contains("family\tAlu\t1\t13\t13"));
    }

    #[test]
    fn mask_repeatmasker_annotations_to_fasta_masks_cat_annotations() {
        let temp_dir = TempDir::new().unwrap();
        let fasta_path = temp_dir.path().join("sample.fa");
        let masked_path = temp_dir.path().join("sample.fa.masked");

        fs::write(&fasta_path, ">seq1\nACGTACGT\n").unwrap();

        let annotations = vec![parse_cat_annotation_line(
            " 2134 10.70 3.30 1.32 seq1 2 5 (3) AluSx#SINE/Alu 3 6 (1) 1",
            1,
        )
        .unwrap()
        .unwrap()];

        let stats = mask_repeatmasker_annotations_to_fasta(
            &annotations,
            &fasta_path,
            &masked_path,
            MaskMode::Lowercase,
            50,
        )
        .unwrap();

        assert_eq!(
            fs::read_to_string(masked_path).unwrap(),
            ">seq1\nAcgtaCGT\n"
        );
        assert_eq!(stats.masked_bases, 4);
    }
}
