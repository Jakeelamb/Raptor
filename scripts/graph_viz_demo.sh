#!/bin/bash
set -e

# Graph Visualization Demo Script
echo "Graph Visualization Demo Script"
echo "=============================="

# Create directory for test data if it doesn't exist
mkdir -p vis_demo
cd vis_demo

# Create a sample GFA file with paths and node annotations
cat > sample_graph.gfa << 'EOF'
H	VN:Z:1.0
S	1	ACGTACGTACGT	RC:f:0.25	CV:f:10.5
S	2	TTTTGGGGCCCC	RC:f:0.75	CV:f:5.2
S	3	AAAATTTTGGGG	RC:f:0.33	CV:f:7.8
S	4	CCCCAAAATTTT	RC:f:0.42	CV:f:12.3
P	transcript_1	1+,2+,3+	*
P	transcript_2	1+,3-,4+	*
P	transcript_3	2+,3+,4-	*
P	transcript_4	1+,2+,4+	*
EOF

echo "Created sample GFA file with node annotations and paths."

# Create segment sequences file
cat > segment_sequences.tsv << 'EOF'
1	ACGTACGTACGT
2	TTTTGGGGCCCC
3	AAAATTTTGGGG
4	CCCCAAAATTTT
EOF

echo "Created segment sequences file."

# Create coverage data file
cat > coverage.tsv << 'EOF'
1	10.5
2	5.2
3	7.8
4	12.3
EOF

echo "Created coverage data file."

# Prepare path map and segment map for ODGI export
cat > prepare_odgi_data.py << 'EOF'
import json
import sys

def main():
    # Read GFA file
    with open('sample_graph.gfa', 'r') as f:
        lines = f.readlines()
    
    # Parse segment sequences
    segments = {}
    with open('segment_sequences.tsv', 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                segments[parts[0]] = parts[1]
    
    # Parse coverage
    coverage = {}
    with open('coverage.tsv', 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                coverage[parts[0]] = float(parts[1])
    
    # Parse paths
    paths = {}
    for line in lines:
        if line.startswith('P'):
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                path_id = parts[1]
                path_segments = []
                for seg in parts[2].split(','):
                    if seg.endswith('+') or seg.endswith('-'):
                        seg_id = seg[:-1]
                        dir_char = seg[-1]
                        path_segments.append([seg_id, dir_char])
                paths[path_id] = path_segments
    
    # Convert to expected format
    path_map = {}
    for path_id, segments_list in paths.items():
        path_map[path_id] = [[s, d] for s, d in segments_list]
    
    # Create output data
    data = {
        "segment_map": segments,
        "path_map": path_map,
        "coverage": coverage
    }
    
    # Write to JSON file
    with open('graph_data.json', 'w') as f:
        json.dump(data, f, indent=2)
    
    print("Created graph_data.json with segment_map, path_map and coverage.")

if __name__ == "__main__":
    main()
EOF

# Run the python script to prepare data
python prepare_odgi_data.py

# Create Rust program to generate ODGI JSON
cat > generate_odgi_json.rs << 'EOF'
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::Read;
use serde::{Deserialize, Serialize};
use serde_json::Value;

#[derive(Deserialize)]
struct GraphData {
    segment_map: HashMap<String, String>,
    path_map: HashMap<String, Vec<Vec<String>>>,
    coverage: HashMap<String, f32>,
}

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

fn main() -> Result<(), Box<dyn Error>> {
    // Read the JSON data
    let mut file = File::open("graph_data.json")?;
    let mut content = String::new();
    file.read_to_string(&mut content)?;
    
    let data: GraphData = serde_json::from_str(&content)?;
    
    // Convert path_map format to expected structure
    let mut converted_path_map = HashMap::new();
    for (id, segments) in &data.path_map {
        let converted = segments.iter()
            .map(|seg| (seg[0].clone(), seg[1].chars().next().unwrap()))
            .collect::<Vec<_>>();
        converted_path_map.insert(id.clone(), converted);
    }
    
    // Generate ODGI JSON
    let mut nodes = vec![];
    
    // Process nodes (segments)
    for (id, seq) in &data.segment_map {
        let rle_len = rle_encode(seq).len();
        let rc = 1.0 - (rle_len as f32 / seq.len() as f32);
        let cov = *data.coverage.get(id).unwrap_or(&1.0);
        
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
    for (id, segments) in &converted_path_map {
        let segments_vec: Vec<String> = segments.iter()
            .map(|(s, dir)| format!("{}{}", s, dir))
            .collect();
            
        let mut metadata = HashMap::new();
        metadata.insert("length".to_string(), segments.len() as f32);
        
        // Calculate number of unique segments
        let mut unique_segments = std::collections::HashSet::new();
        for (s, _) in segments {
            unique_segments.insert(s);
        }
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
    for segments in converted_path_map.values() {
        for (seg, _) in segments {
            *segment_usage.entry(seg).or_insert(0) += 1;
        }
    }
    let shared_segments = segment_usage.values().filter(|&&v| v > 1).count();
    graph_metadata.insert("shared_segments".to_string(), shared_segments as f32);
    graph_metadata.insert("branchiness".to_string(), 
        if data.segment_map.len() > 0 { 
            shared_segments as f32 / data.segment_map.len() as f32 
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
    let file = File::create("odgi_graph.json")?;
    serde_json::to_writer_pretty(file, &graph)?;
    
    println!("Generated ODGI-compatible JSON file: odgi_graph.json");
    
    // Calculate and print the graph complexity
    println!("\nGraph Complexity Analysis:");
    println!("Total segments: {}", data.segment_map.len());
    println!("Total paths: {}", converted_path_map.len());
    println!("Shared segments: {} ({}%)", 
        shared_segments,
        (shared_segments as f32 / data.segment_map.len() as f32 * 100.0) as usize
    );
    println!("Graph branchiness: {:.1}%", 
        shared_segments as f32 / data.segment_map.len() as f32 * 100.0
    );
    
    Ok(())
}

// Simple RLE encoder for demo
fn rle_encode(seq: &str) -> Vec<(u8, u8)> {
    let mut result = Vec::new();
    let mut chars = seq.bytes().peekable();

    while let Some(b) = chars.next() {
        let mut count = 1;
        while let Some(&next) = chars.peek() {
            if next == b {
                count += 1;
                chars.next();
            } else {
                break;
            }
        }
        result.push((b, count));
    }

    result
}
EOF

echo "Compile and run the Rust program to generate ODGI JSON"
rustc generate_odgi_json.rs && ./generate_odgi_json

# Create a GFA2 file with segment annotations
cat > gfa2_annotated.rs << 'EOF'
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use serde::Deserialize;

#[derive(Deserialize)]
struct GraphData {
    segment_map: HashMap<String, String>,
    path_map: HashMap<String, Vec<Vec<String>>>,
    coverage: HashMap<String, f32>,
}

fn main() -> Result<(), Box<dyn Error>> {
    // Read the JSON data
    let mut file = File::open("graph_data.json")?;
    let mut content = String::new();
    file.read_to_string(&mut content)?;
    
    let data: GraphData = serde_json::from_str(&content)?;
    
    // Create GFA2 file
    let gfa2_file = File::create("annotated_graph.gfa2")?;
    let mut writer = BufWriter::new(gfa2_file);
    
    // Write header
    writeln!(writer, "H\tVN:Z:2.0")?;
    
    // Write segments with annotations
    for (id, sequence) in &data.segment_map {
        // Calculate RLE compression
        let rle = rle_encode(sequence);
        let rle_ratio = if sequence.len() > 0 {
            1.0 - (rle.len() as f32 / sequence.len() as f32)
        } else {
            0.0
        };
        
        // Get coverage
        let coverage = data.coverage.get(id).copied().unwrap_or(1.0);
        
        // Determine color based on RLE ratio for BandageNG
        let color = if rle_ratio > 0.7 {
            "red"
        } else if rle_ratio > 0.4 {
            "orange"
        } else {
            "green"
        };
        
        // Write segment with additional tags
        writeln!(
            writer,
            "S\t{}\t{}\t*\tRC:f:{:.3}\tCV:f:{:.2}\tCL:Z:{}",
            id, sequence.len(), rle_ratio, coverage, color
        )?;
        writeln!(writer, "a\t{}\tseq\t{}", id, sequence)?;
    }
    
    // Write paths
    for (id, segments) in &data.path_map {
        let path_str: Vec<String> = segments.iter()
            .map(|s| format!("{}{}", s[0], s[1]))
            .collect();
        
        writeln!(writer, "O\t{}\t{}\t*", id, path_str.join(" "))?;
    }
    
    println!("Generated annotated GFA2 file: annotated_graph.gfa2");
    Ok(())
}

// Simple RLE encoder for demo
fn rle_encode(seq: &str) -> Vec<(u8, u8)> {
    let mut result = Vec::new();
    let mut chars = seq.bytes().peekable();

    while let Some(b) = chars.next() {
        let mut count = 1;
        while let Some(&next) = chars.peek() {
            if next == b {
                count += 1;
                chars.next();
            } else {
                break;
            }
        }
        result.push((b, count));
    }

    result
}
EOF

echo "Compile and run the Rust program to generate annotated GFA2"
rustc gfa2_annotated.rs && ./gfa2_annotated

echo ""
echo "Generated files:"
ls -l

echo ""
echo "Graph complexity metrics were printed above."
echo ""
echo "To visualize the GFA2 file, you can use BandageNG or ODGI:"
echo "- BandageNG: Load annotated_graph.gfa2"
echo "- ODGI: Convert GFA to ODGI format with 'odgi build' and visualize"
echo ""
echo "For web-based visualization, you can use:"
echo "- OGDI Viz: http://ogdi-viz.example.org (load odgi_graph.json)"
echo "- Graphviz for DOT format visualization" 