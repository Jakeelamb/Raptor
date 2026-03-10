use crate::graph::stitch::Path;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader};

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
    let mut paths = HashMap::with_capacity(lines.len());

    for line in lines {
        if let Some((path_name, segments)) = parse_path_record(line) {
            paths.insert(path_name, segments);
        }
    }

    paths
}

/// Parse GFA `P` lines from a buffered reader.
///
/// Unlike `parse_gfa_paths`, this API does not allocate an intermediate
/// collection of all lines and preserves I/O errors from the source reader.
pub fn parse_gfa_paths_reader<R: BufRead>(
    reader: R,
) -> io::Result<HashMap<String, Vec<(String, char)>>> {
    let mut paths = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if let Some((path_name, segments)) = parse_path_record(&line) {
            paths.insert(path_name, segments);
        }
    }

    Ok(paths)
}

#[inline]
fn parse_path_record(line: &str) -> Option<(String, Vec<(String, char)>)> {
    if !line.starts_with("P\t") {
        return None;
    }

    let mut parts = line.split('\t');
    let _record_type = parts.next()?;
    let path_name = parts.next()?.to_string();
    let segment_str = parts.next()?;

    // Parse segment string (e.g., "1+,2-,3+")
    let mut segments = Vec::new();
    for seg in segment_str.split(',') {
        if seg.len() < 2 {
            continue;
        }

        let orientation = seg.chars().last().unwrap_or('+');
        if orientation != '+' && orientation != '-' {
            continue;
        }

        let segment_id = seg[..seg.len() - 1].to_string();
        if segment_id.is_empty() {
            continue;
        }
        segments.push((segment_id, orientation));
    }

    Some((path_name, segments))
}

/// Load all `P`-lines from GFA file and parse them into path navigation
pub fn load_and_parse_gfa_paths(
    path: &str,
) -> Result<HashMap<String, Vec<(String, char)>>, std::io::Error> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    parse_gfa_paths_reader(reader)
}

/// Load segment sequences from a TSV file (segment_id\tsequence)
pub fn load_segment_sequences(path: &str) -> io::Result<HashMap<String, String>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut segments = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();

        if parts.len() >= 2 {
            let segment_id = parts[0].to_string();
            let sequence = parts[1].to_string();
            segments.insert(segment_id, sequence);
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
    paths: &HashMap<String, Vec<(String, char)>>,
    segment_map: &HashMap<String, String>,
) -> HashMap<String, String> {
    let mut reconstructed = HashMap::new();

    for (path_id, segments) in paths {
        let mut sequence = String::new();

        for segment in segments {
            let (segment_id, orientation) = segment;
            if let Some(seg_seq) = segment_map.get(segment_id) {
                if *orientation == '+' {
                    sequence.push_str(seg_seq);
                } else {
                    // Reverse complement for negative orientation
                    sequence.push_str(&reverse_complement(seg_seq));
                }
            }
        }

        reconstructed.insert(path_id.clone(), sequence);
    }

    reconstructed
}

/// Generate reverse complement of a DNA sequence
fn reverse_complement(sequence: &str) -> String {
    let mut result = String::with_capacity(sequence.len());

    for c in sequence.chars().rev() {
        let complement = match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'G' | 'g' => 'C',
            'C' | 'c' => 'G',
            'N' | 'n' => 'N',
            _ => c, // Preserve non-standard characters
        };
        result.push(complement);
    }

    result
}

/// Extract path information from GFA for visualization tools like Bandage/ODGI
pub fn extract_path_metadata(path_map: &HashMap<String, Vec<(String, char)>>) -> Vec<PathMetadata> {
    let mut metadata: Vec<PathMetadata> = path_map
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
        .collect();
    metadata.sort_unstable_by(|a, b| a.id.cmp(&b.id));
    metadata
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
    use std::io::Cursor;

    #[test]
    fn test_path_traversal() {
        let path = Path {
            id: 0,
            segments: vec![0, 2, 4],
            overlaps: vec![5, 10],
        };

        // Test with edges
        let result = traverse_path(&path, true);
        assert_eq!(
            result,
            vec!["contig_1+", "->", "contig_3+", "->", "contig_5+"]
        );

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
        assert_eq!(
            paths["path1"],
            vec![
                ("1".to_string(), '+'),
                ("2".to_string(), '+'),
                ("3".to_string(), '+')
            ]
        );
        assert_eq!(
            paths["path2"],
            vec![("1".to_string(), '+'), ("3".to_string(), '-')]
        );
    }

    #[test]
    fn test_parse_gfa_paths_reader_matches_slice_parser() {
        let lines = vec![
            "H\tVN:Z:1.0".to_string(),
            "P\tpath2\t2+,3-\t*".to_string(),
            "P\tpath1\t1+,2+,3+\t*".to_string(),
        ];
        let expected = parse_gfa_paths(&lines);
        let raw = lines.join("\n");
        let observed = parse_gfa_paths_reader(Cursor::new(raw)).unwrap();
        assert_eq!(observed, expected);
    }

    #[test]
    fn test_revcomp() {
        assert_eq!(revcomp("ACGT"), "ACGT");
        assert_eq!(revcomp("AAAAAA"), "TTTTTT");
        assert_eq!(revcomp("GATTACA"), "TGTAATC");
    }

    #[test]
    fn extract_path_metadata_is_sorted_by_id() {
        let mut paths = HashMap::new();
        paths.insert(
            "path2".to_string(),
            vec![("1".to_string(), '+'), ("2".to_string(), '-')],
        );
        paths.insert("path1".to_string(), vec![("3".to_string(), '+')]);

        let metadata = extract_path_metadata(&paths);
        let ids: Vec<&str> = metadata.iter().map(|m| m.id.as_str()).collect();
        assert_eq!(ids, vec!["path1", "path2"]);
    }
}
