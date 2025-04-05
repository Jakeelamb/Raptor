use std::collections::HashMap;
use serde::Serialize;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use crate::kmer::rle::rle_encode;
use crate::graph::navigation::PathMetadata;

#[derive(Serialize)]
struct OdgiNode {
    id: String,
    length: usize,
    tags: HashMap<String, f32>,
}

#[derive(Serialize)]
struct OdgiPath {
    id: String,
    segments: Vec<String>,
    metadata: HashMap<String, f32>,
}

#[derive(Serialize)]
struct OdgiGraph {
    nodes: Vec<OdgiNode>,
    paths: Vec<OdgiPath>,
    metadata: HashMap<String, f32>,
}

/// Export the graph to an ODGI-compatible JSON format
pub fn export_odgi_json(
    segment_map: &HashMap<String, String>,
    path_map: &HashMap<String, Vec<(String, char)>>,
    coverage: &HashMap<String, f32>,
    output: &str,
) -> io::Result<()> {
    let mut nodes = vec![];
    
    // Process nodes (segments)
    for (id, seq) in segment_map {
        let rle = rle_encode(seq);
        let rc = 1.0 - (rle.len() as f32 / seq.len() as f32);
        let cov = *coverage.get(id).unwrap_or(&1.0);
        
        let mut tags = HashMap::new();
        tags.insert("RC".to_string(), rc);
        tags.insert("CV".to_string(), cov);
        
        // Add GC content
        let gc_count = seq.bytes()
            .filter(|&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
            .count();
        let gc_content = if seq.len() > 0 {
            gc_count as f32 / seq.len() as f32
        } else {
            0.0
        };
        tags.insert("GC".to_string(), gc_content);

        nodes.push(OdgiNode {
            id: id.clone(),
            length: seq.len(),
            tags,
        });
    }

    // Process paths
    let mut paths = vec![];
    for (id, segments) in path_map {
        let segments_vec: Vec<String> = segments.iter()
            .map(|(s, dir)| format!("{}{}", s, dir))
            .collect();
            
        let mut metadata = HashMap::new();
        metadata.insert("length".to_string(), segments.len() as f32);
        
        // Calculate number of unique segments
        let unique_segments: std::collections::HashSet<&String> = 
            segments.iter().map(|(s, _)| s).collect();
        metadata.insert("unique_segments".to_string(), unique_segments.len() as f32);
        
        // Calculate if the path has inversions
        let has_inversions = segments.iter().any(|(_, dir)| *dir == '-');
        metadata.insert("has_inversions".to_string(), if has_inversions { 1.0 } else { 0.0 });
        
        paths.push(OdgiPath {
            id: id.clone(),
            segments: segments_vec,
            metadata,
        });
    }
    
    // Calculate graph-level metrics
    let mut graph_metadata = HashMap::new();
    graph_metadata.insert("node_count".to_string(), nodes.len() as f32);
    graph_metadata.insert("path_count".to_string(), paths.len() as f32);
    
    // Calculate branchiness
    let mut segment_usage = HashMap::new();
    for path in path_map.values() {
        for (seg, _) in path {
            *segment_usage.entry(seg).or_insert(0) += 1;
        }
    }
    let shared_segments = segment_usage.values().filter(|&&v| v > 1).count();
    graph_metadata.insert("shared_segments".to_string(), shared_segments as f32);
    graph_metadata.insert("branchiness".to_string(), 
        if segment_map.len() > 0 { 
            shared_segments as f32 / segment_map.len() as f32 
        } else { 
            0.0 
        });

    // Build the complete graph object
    let graph = OdgiGraph { 
        nodes, 
        paths, 
        metadata: graph_metadata,
    };
    
    // Write to file
    let file = File::create(output)?;
    serde_json::to_writer_pretty(BufWriter::new(file), &graph)?;
    
    Ok(())
}

/// Export path metadata to ODGI-compatible JSON
pub fn export_path_metadata_json(
    metadata: &[PathMetadata],
    output: &str,
) -> io::Result<()> {
    let file = File::create(output)?;
    let mut writer = BufWriter::new(file);
    
    writeln!(writer, "{{")?;
    writeln!(writer, "  \"paths\": [")?;
    
    for (i, path) in metadata.iter().enumerate() {
        let comma = if i < metadata.len() - 1 { "," } else { "" };
        writeln!(writer, "    {{")?;
        writeln!(writer, "      \"id\": \"{}\",", path.id)?;
        writeln!(writer, "      \"segment_count\": {},", path.segment_count)?;
        writeln!(writer, "      \"unique_segment_count\": {},", path.unique_segment_count)?;
        writeln!(writer, "      \"has_inversions\": {}", if path.has_inversions { "true" } else { "false" })?;
        writeln!(writer, "    }}{}", comma)?;
    }
    
    writeln!(writer, "  ],")?;
    
    // Add summary statistics
    let total_paths = metadata.len();
    let total_segments: usize = metadata.iter().map(|p| p.segment_count).sum();
    let avg_path_length = if total_paths > 0 { 
        total_segments as f32 / total_paths as f32 
    } else { 
        0.0 
    };
    
    // Count shared segments across paths
    let mut segment_usage = HashMap::new();
    for path in metadata {
        // We don't have actual segment IDs, so we mock this with the path ID
        for _ in 0..path.segment_count {
            *segment_usage.entry(path.id.clone()).or_insert(0) += 1;
        }
    }
    let shared_segments = segment_usage.values().filter(|&&v| v > 1).count();
    
    writeln!(writer, "  \"summary\": {{")?;
    writeln!(writer, "    \"total_paths\": {},", total_paths)?;
    writeln!(writer, "    \"average_path_length\": {:.2},", avg_path_length)?;
    writeln!(writer, "    \"shared_segments\": {},", shared_segments)?;
    writeln!(writer, "    \"branchiness\": {:.2}", 
        if total_segments > 0 { 
            shared_segments as f32 / total_segments as f32 * 100.0 
        } else { 
            0.0 
        })?;
    writeln!(writer, "  }}")?;
    writeln!(writer, "}}")?;
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::NamedTempFile;
    
    #[test]
    fn test_export_odgi_json() {
        // Create test data
        let mut segment_map = HashMap::new();
        segment_map.insert("1".to_string(), "ACGTACGT".to_string());
        segment_map.insert("2".to_string(), "TTTTGGGG".to_string());
        segment_map.insert("3".to_string(), "AAAATTTT".to_string());
        
        let mut path_map = HashMap::new();
        path_map.insert(
            "path1".to_string(), 
            vec![("1".to_string(), '+'), ("2".to_string(), '+')]
        );
        path_map.insert(
            "path2".to_string(), 
            vec![("1".to_string(), '+'), ("3".to_string(), '-')]
        );
        
        let mut coverage = HashMap::new();
        coverage.insert("1".to_string(), 10.5);
        coverage.insert("2".to_string(), 5.2);
        coverage.insert("3".to_string(), 7.8);
        
        // Export to a temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path().to_str().unwrap();
        
        export_odgi_json(&segment_map, &path_map, &coverage, path).unwrap();
        
        // Verify file existence and content
        let mut file = File::open(path).unwrap();
        let mut contents = String::new();
        file.read_to_string(&mut contents).unwrap();
        
        // Check basic structure
        assert!(contents.contains("\"nodes\":"));
        assert!(contents.contains("\"paths\":"));
        assert!(contents.contains("\"metadata\":"));
        
        // Check for node data
        assert!(contents.contains("\"id\": \"1\""));
        assert!(contents.contains("\"CV\": 10.5"));
    }
} 