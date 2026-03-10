use crate::graph::navigation;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, Write};
use tracing::info;

#[inline]
fn sorted_entries<V>(map: &HashMap<String, V>) -> Vec<(&String, &V)> {
    let mut entries: Vec<(&String, &V)> = map.iter().collect();
    entries.sort_unstable_by(|(a, _), (b, _)| a.cmp(b));
    entries
}

/// Traverses paths in a GFA file and exports them in various formats
pub fn traverse_paths(
    input: &str,
    segments: &str,
    output: &str,
    formats: &str,
    include_edges: bool,
    visualize: bool,
    metadata: bool,
) -> io::Result<usize> {
    info!("Traversing paths in GFA file: {}", input);

    // Load segment sequences
    let segment_map = navigation::load_segment_sequences(segments)?;

    // Parse GFA paths without loading the entire file into memory.
    let gfa_file = File::open(input)?;
    let paths = navigation::parse_gfa_paths_reader(BufReader::new(gfa_file))?;

    // Reconstruct path sequences
    let reconstructed_paths = navigation::reconstruct_paths(&paths, &segment_map);

    // Process and export based on requested formats
    let format_list: Vec<&str> = formats.split(',').collect();

    for format in format_list {
        match format.trim() {
            "fasta" => {
                let fasta_path = format!("{}.fasta", output);
                export_paths_to_fasta(&reconstructed_paths, &fasta_path)?;
                info!("Exported path sequences to FASTA: {}", fasta_path);
            }
            "json" => {
                if metadata {
                    let json_path = format!("{}.json", output);
                    export_paths_to_json(&reconstructed_paths, &paths, &json_path)?;
                    info!("Exported path metadata to JSON: {}", json_path);
                }
            }
            "dot" => {
                if visualize {
                    let dot_path = format!("{}.dot", output);
                    export_paths_to_dot(&paths, &segment_map, include_edges, &dot_path)?;
                    info!("Created DOT graph at {}", dot_path);
                }
            }
            _ => {
                info!("Unsupported export format: {}", format);
            }
        }
    }

    info!(
        "Path traversal complete. Found {} paths, exported to requested formats.",
        paths.len()
    );
    Ok(paths.len())
}

/// Export path sequences to FASTA format
fn export_paths_to_fasta(paths: &HashMap<String, String>, output_path: &str) -> io::Result<()> {
    let mut file = File::create(output_path)?;

    for (path_id, sequence) in sorted_entries(paths) {
        writeln!(file, ">{}", path_id)?;

        // Write sequence in chunks of 60 characters
        for i in (0..sequence.len()).step_by(60) {
            let end = (i + 60).min(sequence.len());
            writeln!(file, "{}", &sequence[i..end])?;
        }
    }

    Ok(())
}

/// Export path metadata to JSON
fn export_paths_to_json(
    reconstructed_paths: &HashMap<String, String>,
    paths: &HashMap<String, Vec<(String, char)>>,
    output_path: &str,
) -> io::Result<()> {
    #[derive(serde::Serialize)]
    struct PathMetadata {
        path_id: String,
        segment_count: usize,
        unique_segments: usize,
        total_length: usize,
        has_inversions: bool,
    }

    let mut path_metadatas = Vec::new();

    for (id, sequence) in sorted_entries(reconstructed_paths) {
        // Retrieve the original path segments
        let Some(segments) = paths.get(id) else {
            continue;
        };
        let segment_count = segments.len();

        // Count unique segments
        let mut unique_segment_ids = std::collections::HashSet::new();
        for (seg_id, _) in segments {
            unique_segment_ids.insert(seg_id);
        }

        // Check if path has inversions
        let has_inversions = segments.iter().any(|(_, dir)| *dir == '-');

        // Create metadata
        let meta = PathMetadata {
            path_id: (*id).clone(),
            segment_count,
            unique_segments: unique_segment_ids.len(),
            total_length: sequence.len(),
            has_inversions,
        };

        path_metadatas.push(meta);
    }

    // Serialize and write to file
    let json_data = serde_json::to_string_pretty(&path_metadatas).map_err(io::Error::other)?;

    std::fs::write(output_path, json_data)?;

    Ok(())
}

/// Export paths to DOT format for visualization
fn export_paths_to_dot(
    paths: &HashMap<String, Vec<(String, char)>>,
    segment_map: &HashMap<String, String>,
    include_edges: bool,
    output_path: &str,
) -> io::Result<()> {
    let mut file = File::create(output_path)?;

    // Start DOT file
    writeln!(file, "digraph G {{")?;
    writeln!(file, "  rankdir=LR;")?;
    writeln!(file, "  node [shape=box style=filled];")?;

    // Create node set to avoid duplicates
    let mut nodes = std::collections::HashSet::new();
    let mut edges = std::collections::HashSet::new();

    // Process each path
    for (path_id, segments) in sorted_entries(paths) {
        // Add path as subgraph
        writeln!(file, "  subgraph cluster_{} {{", path_id.replace('-', "_"))?;
        writeln!(file, "    label=\"Path {}\"", path_id)?;
        writeln!(file, "    style=filled;")?;
        writeln!(file, "    color=lightgrey;")?;

        // Add nodes for each segment in this path
        for (i, (seg_id, orientation)) in segments.iter().enumerate() {
            let node_id = format!("{}_{}", seg_id, orientation);
            let label = if let Some(seq) = segment_map.get(seg_id) {
                format!("{} ({} bp)", seg_id, seq.len())
            } else {
                seg_id.clone()
            };

            let color = match orientation {
                '+' => "lightblue",
                '-' => "salmon",
                _ => "white",
            };

            if nodes.insert(node_id.clone()) {
                writeln!(
                    file,
                    "    \"{}\" [label=\"{}\" fillcolor=\"{}\"];",
                    node_id, label, color
                )?;
            }

            // Add edges between consecutive segments
            if i > 0 && include_edges {
                let prev_node_id = format!("{}_{}", segments[i - 1].0, segments[i - 1].1);
                let edge_id = format!("{}->{}", prev_node_id, node_id);

                if edges.insert(edge_id) {
                    writeln!(file, "    \"{}\" -> \"{}\";", prev_node_id, node_id)?;
                }
            }
        }

        writeln!(file, "  }}")?;
    }

    // End DOT file
    writeln!(file, "}}")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::NamedTempFile;

    #[test]
    fn export_paths_to_fasta_is_sorted_by_path_id() {
        let mut paths = HashMap::new();
        paths.insert("path2".to_string(), "CCCC".to_string());
        paths.insert("path1".to_string(), "AAAA".to_string());

        let file = NamedTempFile::new().unwrap();
        export_paths_to_fasta(&paths, file.path().to_str().unwrap()).unwrap();

        let mut contents = String::new();
        File::open(file.path())
            .unwrap()
            .read_to_string(&mut contents)
            .unwrap();
        let lines: Vec<&str> = contents.lines().collect();
        assert_eq!(lines[0], ">path1");
        assert_eq!(lines[2], ">path2");
    }

    #[test]
    fn export_paths_to_json_is_deterministic_under_map_insertion_order() {
        let mut reconstructed_a = HashMap::new();
        reconstructed_a.insert("path2".to_string(), "CCCC".to_string());
        reconstructed_a.insert("path1".to_string(), "AAAA".to_string());
        let mut paths_a = HashMap::new();
        paths_a.insert("path2".to_string(), vec![("2".to_string(), '+')]);
        paths_a.insert(
            "path1".to_string(),
            vec![("1".to_string(), '+'), ("2".to_string(), '-')],
        );

        let mut reconstructed_b = HashMap::new();
        reconstructed_b.insert("path1".to_string(), "AAAA".to_string());
        reconstructed_b.insert("path2".to_string(), "CCCC".to_string());
        let mut paths_b = HashMap::new();
        paths_b.insert(
            "path1".to_string(),
            vec![("1".to_string(), '+'), ("2".to_string(), '-')],
        );
        paths_b.insert("path2".to_string(), vec![("2".to_string(), '+')]);

        let out_a = NamedTempFile::new().unwrap();
        let out_b = NamedTempFile::new().unwrap();
        export_paths_to_json(&reconstructed_a, &paths_a, out_a.path().to_str().unwrap()).unwrap();
        export_paths_to_json(&reconstructed_b, &paths_b, out_b.path().to_str().unwrap()).unwrap();

        let mut json_a = String::new();
        let mut json_b = String::new();
        File::open(out_a.path())
            .unwrap()
            .read_to_string(&mut json_a)
            .unwrap();
        File::open(out_b.path())
            .unwrap()
            .read_to_string(&mut json_b)
            .unwrap();

        assert_eq!(json_a, json_b);
        assert!(
            json_a.find("\"path_id\": \"path1\"").unwrap()
                < json_a.find("\"path_id\": \"path2\"").unwrap()
        );
    }
}
