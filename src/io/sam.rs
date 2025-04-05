use std::collections::HashMap;
use std::io::{BufRead, BufReader, Result as IoResult};
use std::fs::File;

/// Parse SAM file to count transcript hits
/// Returns a HashMap of transcript IDs to hit counts
pub fn parse_sam_transcript_hits(path: &str) -> IoResult<HashMap<String, usize>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut hits = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('@') { continue; } // Skip header lines
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() > 2 {
            let tid = fields[2].to_string();
            if tid != "*" { // Skip unmapped reads
                *hits.entry(tid).or_insert(0) += 1;
            }
        }
    }

    Ok(hits)
} 