use std::collections::HashMap;
use std::error::Error;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::Path;

fn main() -> Result<(), Box<dyn Error>> {
    println!("GFA Path Traversal Example");
    println!("=========================");
    
    // Create sample GFA file with path definitions
    let gfa_content = r#"H	VN:Z:1.0
S	1	ACGTACGTACGT
S	2	TTTTGGGGCCCC
S	3	AAAATTTTGGGG
S	4	CCCCAAAATTTT
P	iso1	1+,2+,3+	*
P	iso2	1+,3-,4+	*
P	iso3	2+,3+,4-	*
"#;
    
    let gfa_path = Path::new("example.gfa");
    fs::write(gfa_path, gfa_content)?;
    println!("Created sample GFA file: {}", gfa_path.display());
    
    // Create segment sequences file
    let segments_content = r#"1	ACGTACGTACGT
2	TTTTGGGGCCCC
3	AAAATTTTGGGG
4	CCCCAAAATTTT
"#;
    
    let segments_path = Path::new("segments.tsv");
    fs::write(segments_path, segments_content)?;
    println!("Created segment sequences file: {}", segments_path.display());
    
    // Load segments
    let mut segments = HashMap::new();
    for line in segments_content.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            segments.insert(parts[0].to_string(), parts[1].to_string());
        }
    }
    
    // Parse GFA paths
    let mut paths = HashMap::new();
    for line in gfa_content.lines() {
        if line.starts_with('P') {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 3 {
                let id = fields[1].to_string();
                let segments_str = fields[2];
                
                let path_segments = segments_str
                    .split(',')
                    .map(|s| {
                        let (seg, dir) = s.split_at(s.len() - 1);
                        (seg.to_string(), dir.chars().next().unwrap_or('+'))
                    })
                    .collect::<Vec<_>>();
                
                paths.insert(id, path_segments);
            }
        }
    }
    
    println!("\nParsed {} paths from GFA file", paths.len());
    
    // Generate reverse complement function
    fn revcomp(seq: &str) -> String {
        seq.chars()
            .rev()
            .map(|c| match c {
                'A' | 'a' => 'T',
                'T' | 't' => 'A',
                'C' | 'c' => 'G',
                'G' | 'g' => 'C',
                _ => 'N',
            })
            .collect()
    }
    
    // Traverse paths and reconstruct sequences
    println!("\nReconstructed path sequences:");
    
    for (id, path) in &paths {
        let mut sequence = String::new();
        
        println!("Path {}: ", id);
        for (i, (seg_id, orientation)) in path.iter().enumerate() {
            print!("{}{}→", seg_id, orientation);
            
            if let Some(seg_seq) = segments.get(seg_id) {
                let oriented_seq = match orientation {
                    '+' => seg_seq.clone(),
                    '-' => revcomp(seg_seq),
                    _ => seg_seq.clone(),
                };
                
                sequence.push_str(&oriented_seq);
            }
        }
        println!();
        
        println!("Sequence: {}", sequence);
        println!();
    }
    
    // Export to FASTA
    let fasta_path = Path::new("paths.fasta");
    let file = File::create(fasta_path)?;
    let mut writer = BufWriter::new(file);
    
    for (id, path) in &paths {
        writeln!(writer, ">{}", id)?;
        
        let mut sequence = String::new();
        for (seg_id, orientation) in path {
            if let Some(seg_seq) = segments.get(seg_id) {
                let oriented_seq = match orientation {
                    '+' => seg_seq.clone(),
                    '-' => revcomp(seg_seq),
                    _ => seg_seq.clone(),
                };
                
                sequence.push_str(&oriented_seq);
            }
        }
        
        // Write sequence in lines of 80 characters
        for chunk in sequence.as_bytes().chunks(80) {
            writeln!(writer, "{}", std::str::from_utf8(chunk).unwrap_or(""))?;
        }
    }
    
    println!("Exported sequences to FASTA file: {}", fasta_path.display());
    println!("\nTo use the CLI tool for path traversal:");
    println!("cargo run --release -- traverse -i example.gfa -s segments.tsv -o output --formats fasta,dot --metadata");
    
    Ok(())
} 