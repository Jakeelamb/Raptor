use crate::graph::stitch::Path;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Traverse GFA path nodes for export or visualization
pub fn traverse_path(path: &Path, include_edges: bool) -> Vec<String> {
    let mut nav = vec![];

    for (i, &id) in path.segments.iter().enumerate() {
        nav.push(format!("contig_{}+", id + 1));
        if include_edges && i < path.segments.len() - 1 {
            nav.push("->".to_string());
        }
    }

    nav
}

/// Traverse GFA `P` lines into segment ID chains
pub fn parse_gfa_paths(lines: &[String]) -> HashMap<String, Vec<(String, char)>> {
    let mut paths = HashMap::new();

    for line in lines {
        if line.starts_with("P") {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 3 {
                let id = fields[1].to_string();
                let segments = fields[2]
                    .split(',')
                    .filter_map(|s| {
                        let (seg, dir) = s.split_at(s.len() - 1);
                        Some((seg.to_string(), dir.chars().next()?))
                    })
                    .collect::<Vec<_>>();
                paths.insert(id, segments);
            }
        }
    }

    paths
}

/// Load all `P`-lines from GFA file and parse them into path navigation
pub fn load_and_parse_gfa_paths(path: &str) -> Result<HashMap<String, Vec<(String, char)>>, std::io::Error> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines().filter_map(Result::ok).collect();
    Ok(parse_gfa_paths(&lines))
}

/// Load segment sequences from a TSV file (segment_id\tsequence)
pub fn load_segment_sequences(path: &str) -> Result<HashMap<String, String>, std::io::Error> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut segments = HashMap::new();
    
    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 2 {
            segments.insert(fields[0].to_string(), fields[1].to_string());
        }
    }
    
    Ok(segments)
}

/// Return node sequence with correct orientation
pub fn get_oriented_sequence(segment: &str, dir: char, map: &HashMap<String, String>) -> String {
    let seq = map.get(segment).unwrap_or(&"".to_string()).clone();
    match dir {
        '+' => seq,
        '-' => revcomp(&seq),
        _ => seq,
    }
}

/// Generate reverse complement of a DNA sequence
pub fn revcomp(seq: &str) -> String {
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

/// Traverse and reconstruct full paths using segment sequences
pub fn reconstruct_paths(
    path_map: &HashMap<String, Vec<(String, char)>>,
    segment_map: &HashMap<String, String>,
) -> HashMap<String, String> {
    let mut isoforms = HashMap::new();
    for (id, path) in path_map {
        let mut seq = String::new();
        for (seg, dir) in path {
            let part = get_oriented_sequence(seg, *dir, segment_map);
            seq.push_str(&part);
        }
        isoforms.insert(id.clone(), seq);
    }
    isoforms
}

/// Extract path information from GFA for visualization tools like Bandage/ODGI
pub fn extract_path_metadata(
    path_map: &HashMap<String, Vec<(String, char)>>,
) -> Vec<PathMetadata> {
    path_map
        .iter()
        .map(|(id, segments)| {
            let length = segments.len();
            let unique_segments: HashSet<&String> = segments.iter().map(|(seg, _)| seg).collect();
            
            PathMetadata {
                id: id.clone(),
                segment_count: length,
                unique_segment_count: unique_segments.len(),
                has_inversions: segments.iter().any(|(_, dir)| *dir == '-'),
            }
        })
        .collect()
}

/// Metadata about a GFA path for analysis and visualization
#[derive(Debug, Clone)]
pub struct PathMetadata {
    pub id: String,
    pub segment_count: usize,
    pub unique_segment_count: usize,
    pub has_inversions: bool,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::stitch::Path;
    
    #[test]
    fn test_path_traversal() {
        let path = Path {
            id: 0,
            segments: vec![0, 2, 4],
            overlaps: vec![5, 10],
        };
        
        // Test with edges
        let result = traverse_path(&path, true);
        assert_eq!(result, vec!["contig_1+", "->", "contig_3+", "->", "contig_5+"]);
        
        // Test without edges
        let result = traverse_path(&path, false);
        assert_eq!(result, vec!["contig_1+", "contig_3+", "contig_5+"]);
    }
    
    #[test]
    fn test_parse_gfa_paths() {
        let lines = vec![
            "H\tVN:Z:1.0".to_string(),
            "S\t1\tACGT".to_string(), 
            "S\t2\tGGGG".to_string(),
            "S\t3\tTTTT".to_string(),
            "P\tpath1\t1+,2+,3+\t*".to_string(),
            "P\tpath2\t1+,3-\t*".to_string(),
        ];
        
        let paths = parse_gfa_paths(&lines);
        assert_eq!(paths.len(), 2);
        assert_eq!(paths["path1"], vec![("1".to_string(), '+'), ("2".to_string(), '+'), ("3".to_string(), '+')]);
        assert_eq!(paths["path2"], vec![("1".to_string(), '+'), ("3".to_string(), '-')]);
    }
    
    #[test]
    fn test_revcomp() {
        assert_eq!(revcomp("ACGT"), "ACGT");
        assert_eq!(revcomp("AAAAAA"), "TTTTTT");
        assert_eq!(revcomp("GATTACA"), "TGTAATC");
    }
} 