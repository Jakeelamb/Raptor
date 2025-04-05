use crate::graph::navigation::{PathMetadata, reconstruct_paths};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

/// Export reconstructed isoform paths to a FASTA file
pub fn export_paths_to_fasta(
    isoforms: &HashMap<String, String>,
    output_path: &str,
) -> io::Result<()> {
    let path = Path::new(output_path);
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for (id, sequence) in isoforms {
        writeln!(writer, ">{}", id)?;
        
        // Write sequence in lines of 80 characters
        for chunk in sequence.as_bytes().chunks(80) {
            if let Ok(line) = std::str::from_utf8(chunk) {
                writeln!(writer, "{}", line)?;
            }
        }
    }

    Ok(())
}

/// Export path metadata to a TSV file for analysis
pub fn export_path_metadata(
    metadata: &[PathMetadata],
    output_path: &str,
) -> io::Result<()> {
    let path = Path::new(output_path);
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    // Write header
    writeln!(
        writer,
        "path_id\tsegment_count\tunique_segment_count\thas_inversions"
    )?;

    // Write data rows
    for path in metadata {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}",
            path.id,
            path.segment_count,
            path.unique_segment_count,
            path.has_inversions
        )?;
    }

    Ok(())
}

/// Create ODGI-compatible path definitions for visualization
pub fn export_odgi_paths(
    path_map: &HashMap<String, Vec<(String, char)>>,
    output_path: &str,
) -> io::Result<()> {
    let path = Path::new(output_path);
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for (id, segments) in path_map {
        let segment_str: Vec<String> = segments
            .iter()
            .map(|(seg, dir)| format!("{}{}", seg, dir))
            .collect();
        
        writeln!(writer, "P\t{}\t{}\t*", id, segment_str.join(","))?;
    }

    Ok(())
}

/// Export path visualization in DOT format for tools like Graphviz
pub fn export_dot_graph(
    path_map: &HashMap<String, Vec<(String, char)>>,
    output_path: &str,
) -> io::Result<()> {
    let path = Path::new(output_path);
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "digraph G {{")?;
    writeln!(writer, "  node [shape=box];")?;
    
    // Create subgraphs for each path
    for (id, segments) in path_map {
        writeln!(writer, "  subgraph cluster_{} {{", id.replace("-", "_"))?;
        writeln!(writer, "    label=\"{}\";", id)?;
        
        // Add path nodes
        for (i, (seg, dir)) in segments.iter().enumerate() {
            let node_id = format!("{}_{}", id.replace("-", "_"), i);
            writeln!(
                writer,
                "    {} [label=\"{}{}\"];",
                node_id, seg, dir
            )?;
            
            // Add edges between consecutive nodes
            if i > 0 {
                let prev_node = format!("{}_{}", id.replace("-", "_"), i-1);
                writeln!(
                    writer,
                    "    {} -> {};",
                    prev_node, node_id
                )?;
            }
        }
        
        writeln!(writer, "  }}")?;
    }
    
    writeln!(writer, "}}")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::NamedTempFile;

    #[test]
    fn test_export_paths_to_fasta() {
        let mut isoforms = HashMap::new();
        isoforms.insert("iso1".to_string(), "ACGTACGTACGT".to_string());
        isoforms.insert("iso2".to_string(), "TTTTGGGGCCCCAAAA".to_string());
        
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path().to_str().unwrap();
        
        export_paths_to_fasta(&isoforms, path).unwrap();
        
        let mut file = File::open(path).unwrap();
        let mut contents = String::new();
        file.read_to_string(&mut contents).unwrap();
        
        assert!(contents.contains(">iso1"));
        assert!(contents.contains("ACGTACGTACGT"));
        assert!(contents.contains(">iso2"));
        assert!(contents.contains("TTTTGGGGCCCCAAAA"));
    }
    
    #[test]
    fn test_export_path_metadata() {
        let metadata = vec![
            PathMetadata {
                id: "path1".to_string(),
                segment_count: 3,
                unique_segment_count: 3,
                has_inversions: false,
            },
            PathMetadata {
                id: "path2".to_string(),
                segment_count: 2,
                unique_segment_count: 2,
                has_inversions: true,
            },
        ];
        
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path().to_str().unwrap();
        
        export_path_metadata(&metadata, path).unwrap();
        
        let mut file = File::open(path).unwrap();
        let mut contents = String::new();
        file.read_to_string(&mut contents).unwrap();
        
        assert!(contents.contains("path_id\tsegment_count\tunique_segment_count\thas_inversions"));
        assert!(contents.contains("path1\t3\t3\tfalse"));
        assert!(contents.contains("path2\t2\t2\ttrue"));
    }
} 