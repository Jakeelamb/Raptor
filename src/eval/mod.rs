//! Evaluation module - optional benchmarking and comparison features
#![allow(dead_code)]

pub mod metrics;
pub mod gtf_compare;

use std::collections::HashSet;
use std::fs;
use std::io::{self, BufRead, BufReader};
use tracing::info;

/// Evaluates assembled transcripts against ground truth
pub fn evaluate_assembly(
    truth_path: &str,
    pred_path: &str,
    output_path: Option<&str>,
) -> io::Result<()> {
    info!("Evaluating assembled transcripts against ground truth");
    info!("Truth file: {}", truth_path);
    info!("Predicted file: {}", pred_path);
    
    // Compare GTF files
    let (tp, fp, fn_) = compare_gtf(truth_path, pred_path)?;
    let (precision, recall, f1) = calculate_metrics(tp, fp, fn_);
    
    // Display results
    println!("Evaluation Results:");
    println!("  True Positives: {}", tp);
    println!("  False Positives: {}", fp);
    println!("  False Negatives: {}", fn_);
    println!("  Precision: {:.4}", precision);
    println!("  Recall: {:.4}", recall);
    println!("  F1 Score: {:.4}", f1);
    
    // Write output if requested
    if let Some(output_path) = output_path {
        let content = format!(
            "metric\tvalue\ntp\t{}\nfp\t{}\nfn\t{}\nprecision\t{:.4}\nrecall\t{:.4}\nf1\t{:.4}\n",
            tp, fp, fn_, precision, recall, f1
        );
        
        fs::write(output_path, content)?;
        info!("Evaluation results written to {}", output_path);
    }
    
    Ok(())
}

/// Compare truth and predicted GTF files and count TP, FP, and FN
pub fn compare_gtf(truth_path: &str, pred_path: &str) -> io::Result<(usize, usize, usize)> {
    let truth_ranges = load_gtf_ranges(truth_path)?;
    let pred_ranges = load_gtf_ranges(pred_path)?;
    
    let tp = pred_ranges.intersection(&truth_ranges).count();
    let fp = pred_ranges.difference(&truth_ranges).count();
    let fn_ = truth_ranges.difference(&pred_ranges).count();
    
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

/// Load GTF ranges for comparison
fn load_gtf_ranges(path: &str) -> io::Result<HashSet<(String, String, usize, usize)>> {
    let mut ranges = HashSet::new();
    let file = fs::File::open(path)?;
    let reader = BufReader::new(file);
    
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
        
        // Add to the ranges set
        ranges.insert((cols[0].to_string(), transcript_id, start, end));
    }
    
    Ok(ranges)
}

/// GFA comparison module for assembly evaluation
pub mod gfa_compare {
    use super::*;
    use std::collections::HashSet;
    
    /// Compare two GFA files and count shared paths
    pub fn compare_gfa(truth_path: &str, pred_path: &str) -> io::Result<(usize, usize, usize)> {
        let truth_paths = load_gfa_paths(truth_path)?;
        let pred_paths = load_gfa_paths(pred_path)?;
        
        let truth_set: HashSet<String> = truth_paths.into_iter().collect();
        let pred_set: HashSet<String> = pred_paths.into_iter().collect();
        
        let tp = truth_set.intersection(&pred_set).count();
        let fp = pred_set.difference(&truth_set).count();
        let fn_ = truth_set.difference(&pred_set).count();
        
        Ok((tp, fp, fn_))
    }
    
    /// Load paths from GFA file
    fn load_gfa_paths(path: &str) -> io::Result<Vec<String>> {
        let mut paths = Vec::new();
        let file = fs::File::open(path)?;
        let reader = BufReader::new(file);
        
        for line_result in reader.lines() {
            let line = line_result?;
            if line.starts_with('P') {
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() >= 3 {
                    // Extract path name and segment list
                    let path_name = parts[1].to_string();
                    let segments = parts[2].to_string();
                    paths.push(format!("{}:{}", path_name, segments));
                }
            }
        }
        
        Ok(paths)
    }
} 