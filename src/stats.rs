use crate::eval::metrics::evaluate_lengths;
use crate::io::fasta::open_fasta;
use serde::Serialize;
use std::io::BufRead;

#[derive(Serialize)]
pub struct Stats {
    pub total_contigs: usize,
    pub total_length: usize,
    pub average_length: f64,
    pub n50: usize,
    pub n75: usize,
    pub n90: usize,
    pub n95: usize,
    pub l50: usize,
    pub l90: usize,
    pub l95: usize,
    pub au_n: f64,
    pub longest_contig: usize,
    // Graph-related stats
    pub path_count: Option<usize>,
    pub avg_path_length: Option<f64>,
    pub branch_count: Option<usize>,
}

pub fn calculate_stats(path: &str) -> Stats {
    let reader = open_fasta(path);
    let mut lengths = vec![];
    let mut in_sequence = false;
    let mut current_len = 0usize;

    for line in reader.lines().map_while(Result::ok) {
        if line.starts_with('>') {
            // If we were in a sequence, add its final length.
            if in_sequence && current_len > 0 {
                lengths.push(current_len);
            }
            in_sequence = true;
            current_len = 0;
        } else if in_sequence {
            current_len += line.trim().len();
        }
    }

    // Add the last sequence if there is one
    if in_sequence && current_len > 0 {
        lengths.push(current_len);
    }

    let length_stats = evaluate_lengths(&lengths);

    Stats {
        total_contigs: length_stats.total,
        total_length: length_stats.total_bases,
        average_length: length_stats.avg_length,
        n50: length_stats.n50,
        n75: length_stats.n75,
        n90: length_stats.n90,
        n95: length_stats.n95,
        l50: length_stats.l50,
        l90: length_stats.l90,
        l95: length_stats.l95,
        au_n: length_stats.au_n,
        longest_contig: length_stats.longest,
        path_count: None,
        avg_path_length: None,
        branch_count: None,
    }
}

/// Generate and display graph complexity stats from a GFA file
pub fn calculate_graph_stats(path: &str) -> Option<crate::graph::complexity::PathStats> {
    // Check if the file exists and has a GFA extension
    if !path.to_lowercase().ends_with(".gfa") {
        return None;
    }

    // Use the new compute_path_stats function that takes a GFA file path
    match crate::graph::complexity::compute_path_stats(path) {
        Ok(stats) => Some(stats),
        Err(e) => {
            eprintln!("Error computing graph stats: {}", e);
            None
        }
    }
}

/// Update Stats object with graph complexity information
pub fn update_with_graph_stats(
    stats: &mut Stats,
    graph_stats: &crate::graph::complexity::PathStats,
) {
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
        assert_eq!(stats.n50, 24);
        assert_eq!(stats.n75, 20);
        assert_eq!(stats.n90, 20);
        assert_eq!(stats.n95, 4);
        assert_eq!(stats.l50, 1);
        assert_eq!(stats.l90, 2);
        assert_eq!(stats.l95, 3);
        assert!((stats.au_n - 20.6666666667).abs() < 1e-6);
        assert_eq!(stats.longest_contig, 24);
    }

    #[test]
    fn test_calculate_stats_multiline_sequences() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">contig_1").unwrap();
        writeln!(file, "ATCGATCG").unwrap();
        writeln!(file, "ATCG").unwrap(); // 12 bp total
        writeln!(file, ">contig_2").unwrap();
        writeln!(file, "GCTA").unwrap(); // 4 bp

        let stats = calculate_stats(file.path().to_str().unwrap());
        assert_eq!(stats.total_contigs, 2);
        assert_eq!(stats.total_length, 16);
        assert_eq!(stats.n50, 12);
        assert_eq!(stats.n75, 12);
        assert_eq!(stats.n90, 4);
        assert_eq!(stats.n95, 4);
        assert_eq!(stats.l50, 1);
        assert_eq!(stats.l90, 2);
        assert_eq!(stats.l95, 2);
        assert!((stats.au_n - 10.0).abs() < 1e-6);
        assert_eq!(stats.longest_contig, 12);
    }
}
