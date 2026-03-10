use crate::graph::assembler::Contig;
use crate::graph::isoform_graph::{build_isoform_graph, find_end_nodes, find_start_nodes};
use crate::graph::isoform_traverse::{filter_paths_by_confidence, find_directed_paths};
use crate::graph::transcript::{assemble_transcripts, Transcript};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use tracing::info;

type Result<T> = std::result::Result<T, Box<dyn Error>>;

/// Load contig coverages from an expression data file
pub fn load_contig_coverages(expression_path: &str) -> Result<HashMap<usize, f64>> {
    info!("Loading contig coverage data from {}", expression_path);
    let file = File::open(expression_path)?;
    let reader = BufReader::new(file);

    let mut coverage_map = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();

        // Skip comments and empty lines
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        // Parse contig ID and coverage
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            if let (Ok(id), Ok(coverage)) = (parts[0].parse::<usize>(), parts[1].parse::<f64>()) {
                coverage_map.insert(id, coverage);
            }
        }
    }

    info!("Loaded coverage data for {} contigs", coverage_map.len());
    Ok(coverage_map)
}

/// Process assembled contigs into transcripts
pub fn process_isoforms(
    contigs: &[Contig],
    links: &[(usize, usize, usize)],
    kmer_counts: &HashMap<usize, usize>,
    _output_path: &str,
    _gtf_path: Option<&str>,
    max_path_depth: usize,
    min_confidence: f64,
    _min_tpm: Option<f64>,
    _strand_aware: bool,
    _bam_path: Option<&str>,
    _long_reads: Option<&str>,
    _get_output_filename: fn(&str, Option<&str>) -> String,
) -> Result<Vec<Transcript>> {
    info!(
        "Processing isoforms from {} contigs with {} links",
        contigs.len(),
        links.len()
    );

    // Build expression map from k-mer counts if available
    let expression_map: HashMap<usize, f64> = kmer_counts
        .iter()
        .map(|(&contig_id, &count)| (contig_id, count as f64))
        .collect();

    // Convert contigs to a HashMap for easier lookup
    let mut contig_map = HashMap::new();
    for contig in contigs {
        contig_map.insert(contig.id, contig.sequence.clone());
    }

    // Build the isoform graph
    let graph = build_isoform_graph(&contig_map, links, &expression_map);

    // Determine start/end nodes using graph structure; fall back to all nodes if degenerate.
    let mut all_nodes: Vec<usize> = graph.nodes().collect();
    all_nodes.sort_unstable();

    let mut start_nodes = find_start_nodes(&graph);
    let mut end_nodes = find_end_nodes(&graph);

    if start_nodes.is_empty() {
        start_nodes = all_nodes.clone();
    }
    if end_nodes.is_empty() {
        end_nodes = all_nodes.clone();
    }

    start_nodes.sort_unstable();
    end_nodes.sort_unstable();

    // Find paths through the graph
    let paths = find_directed_paths(&graph, &start_nodes, &end_nodes, max_path_depth);
    info!("Found {} raw transcript paths", paths.len());

    // Filter paths by confidence (using min_confidence as f32, 100 as min path length, and None for high threshold)
    let mut filtered_paths = filter_paths_by_confidence(&paths, min_confidence as f32, 100, None);
    filtered_paths.sort_unstable_by(|a, b| {
        b.confidence
            .total_cmp(&a.confidence)
            .then_with(|| b.length.cmp(&a.length))
            .then_with(|| a.nodes.cmp(&b.nodes))
    });
    info!("Filtered to {} high-confidence paths", filtered_paths.len());

    // Build overlap-aware transcript sequences using the shared deterministic assembler.
    let transcripts = assemble_transcripts(&filtered_paths, contigs, links, Some(&graph));

    info!("Generated {} transcript isoforms", transcripts.len());
    Ok(transcripts)
}

#[cfg(test)]
mod tests {
    use super::process_isoforms;
    use crate::graph::assembler::Contig;
    use crate::graph::transcript::Transcript;
    use std::collections::HashMap;

    fn long_contig(id: usize, base: char) -> Contig {
        Contig {
            id,
            sequence: base.to_string().repeat(30),
            kmer_path: vec![],
        }
    }

    fn transcript_signature(transcripts: &[Transcript]) -> Vec<(usize, Vec<usize>, String, u64)> {
        transcripts
            .iter()
            .map(|tx| {
                (
                    tx.id,
                    tx.path.clone(),
                    tx.sequence.clone(),
                    tx.confidence.to_bits(),
                )
            })
            .collect()
    }

    fn passthrough_filename(base: &str, _ext: Option<&str>) -> String {
        base.to_string()
    }

    #[test]
    fn process_isoforms_uses_overlap_aware_stitching() {
        let contigs = vec![
            long_contig(0, 'A'),
            long_contig(1, 'C'),
            long_contig(2, 'G'),
            long_contig(3, 'T'),
        ];
        let links = vec![(0, 1, 10), (1, 2, 10), (2, 3, 10)];
        let expression = HashMap::from([(0usize, 12usize), (1, 11), (2, 10), (3, 9)]);

        let transcripts = process_isoforms(
            &contigs,
            &links,
            &expression,
            "unused",
            None,
            8,
            0.0,
            None,
            false,
            None,
            None,
            passthrough_filename,
        )
        .expect("isoform processing should succeed");

        assert_eq!(transcripts.len(), 1);
        assert_eq!(transcripts[0].id, 1);
        assert_eq!(transcripts[0].path, vec![0, 1, 2, 3]);

        let expected = [
            "A".repeat(30),
            "C".repeat(20),
            "G".repeat(20),
            "T".repeat(20),
        ]
        .concat();
        assert_eq!(transcripts[0].sequence, expected);
    }

    #[test]
    fn process_isoforms_is_stable_under_contig_and_link_permutation() {
        let contigs_a = vec![
            long_contig(0, 'A'),
            long_contig(1, 'C'),
            long_contig(2, 'G'),
            long_contig(3, 'T'),
            long_contig(4, 'N'),
        ];
        let contigs_b = vec![
            long_contig(3, 'T'),
            long_contig(1, 'C'),
            long_contig(4, 'N'),
            long_contig(0, 'A'),
            long_contig(2, 'G'),
        ];

        let links_a = vec![(0, 1, 10), (0, 2, 10), (1, 3, 10), (2, 3, 10), (3, 4, 10)];
        let links_b = vec![(3, 4, 10), (2, 3, 10), (0, 2, 10), (1, 3, 10), (0, 1, 10)];

        let expression = HashMap::from([(0usize, 10usize), (1, 10), (2, 10), (3, 10), (4, 10)]);

        let transcripts_a = process_isoforms(
            &contigs_a,
            &links_a,
            &expression,
            "unused",
            None,
            8,
            0.0,
            None,
            false,
            None,
            None,
            passthrough_filename,
        )
        .expect("isoform processing should succeed for ordering A");

        let transcripts_b = process_isoforms(
            &contigs_b,
            &links_b,
            &expression,
            "unused",
            None,
            8,
            0.0,
            None,
            false,
            None,
            None,
            passthrough_filename,
        )
        .expect("isoform processing should succeed for ordering B");

        assert_eq!(
            transcript_signature(&transcripts_a),
            transcript_signature(&transcripts_b)
        );
    }
}
