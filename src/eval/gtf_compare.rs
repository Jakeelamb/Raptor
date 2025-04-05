use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Result as IoResult};

/// A transcript location with chromosome, start, and end positions
#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct TranscriptLocation {
    pub chromosome: String,
    pub start: usize,
    pub end: usize,
    pub transcript_id: String,
}

/// Load transcript start-end ranges from GTF
pub fn load_gtf_ranges(path: &str) -> IoResult<HashSet<TranscriptLocation>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut ranges = HashSet::new();

    for line_result in reader.lines() {
        let line = line_result?;
        if line.starts_with('#') {
            continue; // Skip comment lines
        }

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue; // Skip malformed lines
        }

        // We only care about exon features for this comparison
        if cols[2] != "exon" {
            continue;
        }

        // Parse the attributes field to extract transcript_id
        let attrs = cols[8];
        let transcript_id = if let Some(pos) = attrs.find("transcript_id") {
            let remaining = &attrs[pos + 14..]; // Skip "transcript_id \""
            if let Some(end_pos) = remaining.find('"') {
                remaining[..end_pos].to_string()
            } else {
                continue; // Malformed attribute
            }
        } else {
            continue; // No transcript ID found
        };

        // Parse start and end positions
        let start = match cols[3].parse::<usize>() {
            Ok(s) => s,
            Err(_) => continue,
        };

        let end = match cols[4].parse::<usize>() {
            Ok(e) => e,
            Err(_) => continue,
        };

        ranges.insert(TranscriptLocation {
            chromosome: cols[0].to_string(),
            start,
            end,
            transcript_id,
        });
    }

    Ok(ranges)
}

/// Compare truth GTF with predicted GTF
/// Returns a tuple of (true_positives, false_positives, false_negatives)
pub fn compare_gtf(truth: &str, predicted: &str) -> IoResult<(usize, usize, usize)> {
    let true_ranges = load_gtf_ranges(truth)?;
    let pred_ranges = load_gtf_ranges(predicted)?;

    // Count true positives, false positives, and false negatives
    let tp = pred_ranges.intersection(&true_ranges).count();
    let fp = pred_ranges.difference(&true_ranges).count();
    let fn_ = true_ranges.difference(&pred_ranges).count();

    Ok((tp, fp, fn_))
}

/// Calculate precision, recall, and F1 score
pub fn calculate_metrics(tp: usize, fp: usize, fn_: usize) -> (f64, f64, f64) {
    let precision = if tp + fp > 0 {
        tp as f64 / (tp + fp) as f64
    } else {
        0.0
    };

    let recall = if tp + fn_ > 0 {
        tp as f64 / (tp + fn_) as f64
    } else {
        0.0
    };

    let f1 = if precision + recall > 0.0 {
        2.0 * precision * recall / (precision + recall)
    } else {
        0.0
    };

    (precision, recall, f1)
} 