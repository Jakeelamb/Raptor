use log::{info, warn};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::time::Instant;

use crate::graph::isoform_filter::{filter_similar_transcripts, merge_transcripts};
use crate::graph::isoform_graph::{build_isoform_graph, find_end_nodes, find_start_nodes};
use crate::graph::isoform_traverse::{filter_paths_by_confidence, find_directed_paths};
use crate::graph::transcript::{assemble_transcripts, calculate_transcript_stats};
use crate::io::gfa::{read_gfa_contigs, read_gfa_links};
use crate::io::transcript_io::{
    add_transcripts_to_gfa, write_transcript_stats, TranscriptFastaWriter, TranscriptGtfWriter,
};

/// Load expression data from a TSV file
/// Format: contig_id, coverage
pub fn load_expression_data(expression_path: &str) -> io::Result<HashMap<usize, f64>> {
    let file = File::open(expression_path)?;
    let reader = BufReader::new(file);

    let mut expression_data = HashMap::new();

    for (line_number, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        let trimmed = line.trim();
        if trimmed.starts_with('#') || trimmed.is_empty() {
            continue;
        }

        let mut parts = trimmed.split_whitespace();
        let id_token = parts.next().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid expression row at line {}: missing contig id",
                    line_number + 1
                ),
            )
        })?;
        let coverage_token = parts.next().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid expression row at line {}: missing coverage value",
                    line_number + 1
                ),
            )
        })?;

        let contig_id = id_token.parse::<usize>().map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid expression row at line {}: contig id '{}' is not a usize",
                    line_number + 1,
                    id_token
                ),
            )
        })?;

        let coverage = coverage_token.parse::<f64>().map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid expression row at line {}: coverage '{}' is not a finite float",
                    line_number + 1,
                    coverage_token
                ),
            )
        })?;
        if !coverage.is_finite() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid expression row at line {}: coverage must be finite",
                    line_number + 1
                ),
            ));
        }
        if coverage < 0.0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid expression row at line {}: coverage must be non-negative",
                    line_number + 1
                ),
            ));
        }
        if parts.next().is_some() {
            warn!(
                "Ignoring extra columns in expression row at line {}",
                line_number + 1
            );
        }

        expression_data.insert(contig_id, coverage);
    }

    info!(
        "Loaded expression data for {} contigs",
        expression_data.len()
    );
    Ok(expression_data)
}

#[derive(Debug, Clone, Copy, Default)]
struct OutputFormats {
    fasta: bool,
    gtf: bool,
    gfa: bool,
}

fn parse_output_formats(formats: &str) -> io::Result<OutputFormats> {
    let mut parsed = OutputFormats::default();

    for format in formats
        .split(',')
        .map(str::trim)
        .filter(|value| !value.is_empty())
    {
        match format.to_ascii_lowercase().as_str() {
            "fasta" => parsed.fasta = true,
            "gtf" => parsed.gtf = true,
            "gfa" => parsed.gfa = true,
            unsupported => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "unsupported isoform output format '{unsupported}' (supported: fasta, gtf, gfa)"
                    ),
                ))
            }
        }
    }

    if !parsed.fasta && !parsed.gtf && !parsed.gfa {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "at least one isoform output format is required (fasta, gtf, gfa)",
        ));
    }

    Ok(parsed)
}

/// Run the isoform reconstruction pipeline
pub fn run_isoform_reconstruction(
    input_gfa: &str,
    expression_path: &str,
    output_prefix: &str,
    min_confidence: f64,
    max_depth: usize,
    formats: &str,
    output_stats: bool,
    similarity_threshold: Option<f64>,
    merge_similar: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let start_time = Instant::now();
    info!("Starting isoform reconstruction");

    // Load contigs and overlaps from GFA
    info!("Loading contigs from GFA");
    let contigs = read_gfa_contigs(input_gfa)?;
    info!("Loaded {} contigs", contigs.len());

    info!("Loading overlaps from GFA");
    let overlaps = read_gfa_links(input_gfa)?;
    info!("Loaded {} overlaps", overlaps.len());

    // Load expression data
    info!("Loading expression data");
    let expression_data = load_expression_data(expression_path)?;

    // Build isoform graph
    info!("Building isoform graph");
    let mut contig_map = HashMap::new();
    for contig in &contigs {
        contig_map.insert(contig.id, contig.sequence.clone());
    }
    let overlaps_vec = overlaps.to_vec();
    let graph = build_isoform_graph(&contig_map, &overlaps_vec, &expression_data);
    let node_count = graph.node_count();
    let edge_count = graph.edge_count();
    info!(
        "Built graph with {} nodes and {} edges",
        node_count, edge_count
    );

    // Find start and end nodes
    let start_nodes = find_start_nodes(&graph);
    let end_nodes = find_end_nodes(&graph);
    info!(
        "Identified {} potential start nodes and {} potential end nodes",
        start_nodes.len(),
        end_nodes.len()
    );

    // Find transcript paths
    info!(
        "Finding directed paths through the graph (max depth: {})",
        max_depth
    );
    let paths = find_directed_paths(&graph, &start_nodes, &end_nodes, max_depth);
    info!("Found {} potential transcript paths", paths.len());

    // Filter paths by confidence
    info!("Filtering paths by confidence (min: {})", min_confidence);
    let min_confidence_f32 = min_confidence as f32;
    let min_path_len = 50; // Minimum path length for standard filtering
    let high_conf_threshold = Some(0.9); // Higher threshold for shorter paths
    info!(
        "Using minimum path length {} with high confidence threshold {:.2}",
        min_path_len,
        high_conf_threshold.unwrap_or(1.0)
    );
    let filtered_paths = filter_paths_by_confidence(
        &paths,
        min_confidence_f32,
        min_path_len,
        high_conf_threshold,
    );
    info!(
        "Retained {} transcript paths after filtering",
        filtered_paths.len()
    );

    // Assemble transcripts
    info!("Assembling transcripts");
    let mut transcripts = assemble_transcripts(&filtered_paths, &contigs, &overlaps, None);
    info!("Assembled {} transcripts", transcripts.len());

    // Apply similarity filtering if requested
    if let Some(threshold) = similarity_threshold {
        info!(
            "Filtering similar transcripts (similarity threshold: {})",
            threshold
        );
        if merge_similar {
            info!("Merging similar transcripts");
            transcripts = merge_transcripts(&transcripts, threshold);
            info!("After merging, {} transcripts remain", transcripts.len());
        } else {
            transcripts = filter_similar_transcripts(&transcripts, threshold);
            info!("After filtering, {} transcripts remain", transcripts.len());
        }
    }

    // Output transcripts in requested formats
    let output_formats = parse_output_formats(formats)?;

    if output_formats.fasta {
        let fasta_path = format!("{}.transcripts.fa", output_prefix);
        info!("Writing transcripts to FASTA: {}", fasta_path);
        let mut fasta_writer = TranscriptFastaWriter::new(&fasta_path);
        fasta_writer.write_transcripts(&transcripts)?;
    }

    if output_formats.gtf {
        let gtf_path = format!("{}.transcripts.gtf", output_prefix);
        info!("Writing transcripts to GTF: {}", gtf_path);
        let mut gtf_writer = TranscriptGtfWriter::new(&gtf_path);
        gtf_writer.write_transcripts(&transcripts)?;
    }

    if output_formats.gfa {
        let gfa_out_path = format!("{}.transcripts.gfa", output_prefix);

        // If input GFA exists, copy it to output GFA first
        if Path::new(input_gfa).exists() && input_gfa != gfa_out_path {
            info!("Copying input GFA to output GFA");
            std::fs::copy(input_gfa, &gfa_out_path)?;
        }

        info!("Adding transcript paths to GFA: {}", gfa_out_path);
        add_transcripts_to_gfa(&gfa_out_path, &transcripts)?;
    }

    // Calculate and output statistics
    if output_stats {
        let stats_path = format!("{}.transcript_stats.json", output_prefix);
        info!("Calculating transcript statistics");
        let stats = calculate_transcript_stats(&transcripts);

        info!("Writing transcript statistics to: {}", stats_path);
        write_transcript_stats(&stats_path, &stats)?;

        // Print some key statistics
        println!("Transcript count: {}", stats.get("count").unwrap_or(&0.0));
        println!(
            "Average length: {:.1}",
            stats.get("mean_length").unwrap_or(&0.0)
        );
        println!("N50: {}", stats.get("n50").unwrap_or(&0.0));
        println!(
            "Average confidence: {:.4}",
            stats.get("mean_confidence").unwrap_or(&0.0)
        );
    }

    let elapsed = start_time.elapsed();
    info!(
        "Isoform reconstruction completed in {:.2}s",
        elapsed.as_secs_f64()
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{load_expression_data, parse_output_formats};
    use std::io::{ErrorKind, Write};
    use tempfile::NamedTempFile;

    #[test]
    fn load_expression_data_parses_valid_rows_and_ignores_comments() {
        let mut file = NamedTempFile::new().expect("create temp expression file");
        writeln!(file, "# comment").expect("write comment");
        writeln!(file, "0\t10.5").expect("write tsv row");
        writeln!(file, "  ").expect("write empty row");
        writeln!(file, "1 8.25 extra").expect("write row with extra columns");
        writeln!(file, "1\t9.75").expect("write duplicate id overwrite row");
        file.flush().expect("flush expression file");

        let expression = load_expression_data(file.path().to_str().expect("utf8 path"))
            .expect("valid expression file should parse");
        assert_eq!(expression.len(), 2);
        assert_eq!(expression.get(&0), Some(&10.5));
        assert_eq!(expression.get(&1), Some(&9.75));
    }

    #[test]
    fn load_expression_data_rejects_invalid_rows_with_context() {
        let mut file = NamedTempFile::new().expect("create temp expression file");
        writeln!(file, "abc\t10.0").expect("write invalid id row");
        file.flush().expect("flush expression file");

        let err = load_expression_data(file.path().to_str().expect("utf8 path"))
            .expect_err("invalid row must fail");
        assert_eq!(err.kind(), ErrorKind::InvalidData);
        assert!(err.to_string().contains("line 1"));
    }

    #[test]
    fn load_expression_data_rejects_non_finite_and_negative_coverage() {
        let mut non_finite = NamedTempFile::new().expect("create temp expression file");
        writeln!(non_finite, "0\tNaN").expect("write non-finite row");
        non_finite.flush().expect("flush expression file");
        let err = load_expression_data(non_finite.path().to_str().expect("utf8 path"))
            .expect_err("NaN coverage must fail");
        assert_eq!(err.kind(), ErrorKind::InvalidData);
        assert!(err.to_string().contains("finite"));

        let mut negative = NamedTempFile::new().expect("create temp expression file");
        writeln!(negative, "0\t-1.0").expect("write negative row");
        negative.flush().expect("flush expression file");
        let err = load_expression_data(negative.path().to_str().expect("utf8 path"))
            .expect_err("negative coverage must fail");
        assert_eq!(err.kind(), ErrorKind::InvalidData);
        assert!(err.to_string().contains("non-negative"));
    }

    #[test]
    fn parse_output_formats_normalizes_case_and_whitespace() {
        let parsed = parse_output_formats("  FASTA, gTf ,gfa ")
            .expect("mixed-case formats should parse successfully");
        assert!(parsed.fasta);
        assert!(parsed.gtf);
        assert!(parsed.gfa);
    }

    #[test]
    fn parse_output_formats_rejects_unknown_and_empty_input() {
        let unknown = parse_output_formats("fasta,bed")
            .expect_err("unsupported format should return an error");
        assert_eq!(unknown.kind(), ErrorKind::InvalidInput);

        let empty =
            parse_output_formats("  , ").expect_err("empty format string should return an error");
        assert_eq!(empty.kind(), ErrorKind::InvalidInput);
    }
}
