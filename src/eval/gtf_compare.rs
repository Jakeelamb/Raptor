use std::collections::{HashSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, Result as IoResult};
use bio::data_structures::interval_tree::IntervalTree;

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

/// Load exon information from a GTF file into interval trees for more complex overlap analysis
pub fn load_exon_intervals(
    path: &str,
) -> IoResult<HashMap<String, IntervalTree<usize, String>>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut interval_map: HashMap<String, IntervalTree<usize, String>> = HashMap::new();

    for line_result in reader.lines() {
        let line = line_result?;
        if line.starts_with('#') {
            continue; // Skip comment lines
        }

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue; // Skip malformed lines
        }

        // We only care about exon features
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

        let chrom = cols[0].to_string();
        interval_map
            .entry(chrom)
            .or_insert_with(IntervalTree::new)
            .insert(start..end, transcript_id);
    }

    Ok(interval_map)
}

/// Advanced GTF comparison using interval trees to find overlapping exons
pub fn compare_gtf_advanced(
    truth: &str,
    predicted: &str,
) -> IoResult<(usize, usize, usize, f64)> {
    // Load interval trees for both files
    let true_intervals = load_exon_intervals(truth)?;
    let pred_intervals = load_exon_intervals(predicted)?;

    let mut tp = 0;
    let mut fp = 0;
    let mut fn_ = 0;
    let mut overlap_sum = 0.0;
    let mut overlap_count = 0;

    // Compute false positives and true positives
    for (chrom, intervals) in &pred_intervals {
        if let Some(true_tree) = true_intervals.get(chrom) {
            // For each predicted exon
            for interval in intervals.intervals() {
                let pred_start = interval.start;
                let pred_end = interval.end;
                let pred_len = pred_end - pred_start;

                // Find overlapping true exons
                let overlaps = true_tree.find(pred_start..pred_end);
                if overlaps.is_empty() {
                    fp += 1; // No overlapping exon found
                } else {
                    tp += 1; // Found at least one overlapping exon
                    
                    // Calculate overlap percentage for evaluation
                    for overlap in overlaps {
                        if let Some(range) = overlap.range() {
                            let overlap_start = range.start.max(pred_start);
                            let overlap_end = range.end.min(pred_end);
                            let overlap_len = overlap_end - overlap_start;
                            let overlap_percent = overlap_len as f64 / pred_len as f64;
                            
                            overlap_sum += overlap_percent;
                            overlap_count += 1;
                        }
                    }
                }
            }
        } else {
            // All exons in a chromosome not in reference are FP
            fp += intervals.intervals().count();
        }
    }

    // Compute false negatives (true exons not overlapped by any predicted exon)
    for (chrom, intervals) in &true_intervals {
        if let Some(pred_tree) = pred_intervals.get(chrom) {
            for interval in intervals.intervals() {
                let true_start = interval.start;
                let true_end = interval.end;
                
                let overlaps = pred_tree.find(true_start..true_end);
                if overlaps.is_empty() {
                    fn_ += 1; // No predicted exon overlapping this true exon
                }
            }
        } else {
            // All exons in a chromosome not in predictions are FN
            fn_ += intervals.intervals().count();
        }
    }

    let avg_overlap = if overlap_count > 0 {
        overlap_sum / overlap_count as f64
    } else {
        0.0
    };

    Ok((tp, fp, fn_, avg_overlap))
} 