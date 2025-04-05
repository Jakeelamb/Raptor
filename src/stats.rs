use crate::io::fasta::open_fasta;
use std::fs::File;
use std::io::{BufReader, BufRead};
use serde::{Serialize};

#[derive(Serialize)]
pub struct Stats {
    pub total_contigs: usize,
    pub total_length: usize,
    pub average_length: f64,
    pub n50: usize,
    // Graph-related stats
    pub path_count: Option<usize>,
    pub avg_path_length: Option<f64>,
    pub branch_count: Option<usize>,
}

pub fn calculate_stats(path: &str) -> Stats {
    let reader = open_fasta(path);
    let mut lengths = vec![];
    let mut total = 0;
    let mut in_sequence = false;
    let mut current_sequence = String::new();

    for line in reader.lines().flatten() {
        if line.starts_with('>') {
            // If we were in a sequence, add the complete sequence to lengths
            if in_sequence && !current_sequence.is_empty() {
                let len = current_sequence.len();
                total += len;
                lengths.push(len);
                current_sequence.clear();
            }
            in_sequence = true;
        } else if in_sequence {
            // Add this line to the current sequence
            current_sequence.push_str(line.trim());
        }
    }

    // Add the last sequence if there is one
    if in_sequence && !current_sequence.is_empty() {
        let len = current_sequence.len();
        total += len;
        lengths.push(len);
    }

    lengths.sort_unstable();
    let total_contigs = lengths.len();
    let avg = if total_contigs > 0 { total as f64 / total_contigs as f64 } else { 0.0 };

    // Calculate N50
    let mut acc = 0;
    let half_total = total / 2;
    let n50 = lengths.iter().rev().find(|&&len| {
        acc += len;
        acc >= half_total
    }).copied().unwrap_or(0);

    Stats {
        total_contigs,
        total_length: total,
        average_length: avg,
        n50,
        path_count: None,
        avg_path_length: None,
        branch_count: None,
    }
}

/// Generate and display graph complexity stats from a GFA file
pub fn calculate_graph_stats(path: &str) -> Option<crate::graph::complexity::PathStats> {
    use crate::graph::stitch::load_paths_from_gfa;
    use crate::graph::complexity::compute_path_stats;
    
    match load_paths_from_gfa(path) {
        Ok(paths) => {
            let stats = compute_path_stats(&paths);
            Some(stats)
        },
        Err(e) => {
            eprintln!("Error loading paths from GFA: {}", e);
            None
        }
    }
}

/// Update Stats object with graph complexity information
pub fn update_with_graph_stats(stats: &mut Stats, graph_stats: &crate::graph::complexity::PathStats) {
    stats.path_count = Some(graph_stats.total_paths);
    stats.avg_path_length = Some(graph_stats.average_length);
    stats.branch_count = Some(graph_stats.branch_count);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_calculate_stats() {
        // Create a temporary FASTA file
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">contig_1").unwrap();
        writeln!(file, "ATCGATCGATCGATCGATCG").unwrap(); // 20 bp
        writeln!(file, ">contig_2").unwrap();
        writeln!(file, "GCTAGCTAGCTAGCTAGCTAGCTA").unwrap(); // 24 bp
        writeln!(file, ">contig_3").unwrap();
        writeln!(file, "ATCG").unwrap(); // 4 bp

        let stats = calculate_stats(file.path().to_str().unwrap());
        
        assert_eq!(stats.total_contigs, 3);
        assert_eq!(stats.total_length, 48);
        assert_eq!(stats.average_length, 16.0);
        assert_eq!(stats.n50, 24); // N50 should be 24
    }
} 