use crate::graph::assembler::Contig;
use crate::graph::navigation::traverse_path;
use crate::graph::stitch::Path;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Result, Write};

pub struct GfaWriter {
    writer: BufWriter<File>,
}

#[derive(Debug, Clone, Copy)]
struct SegmentRecord<'a> {
    id: &'a str,
    sequence: &'a str,
    rle_tag: Option<&'a str>,
}

#[inline]
fn parse_segment_line(line: &str) -> Option<SegmentRecord<'_>> {
    if !line.starts_with("S\t") {
        return None;
    }

    let mut fields = line.split('\t');
    let record_type = fields.next()?;
    debug_assert_eq!(record_type, "S");
    let id = fields.next()?;
    let sequence = fields.next()?;
    let rle_tag = fields.find_map(|field| field.strip_prefix("RN:Z:"));

    Some(SegmentRecord {
        id,
        sequence,
        rle_tag,
    })
}

#[inline]
fn parse_link_line(line: &str) -> Option<(&str, &str, &str)> {
    if !line.starts_with("L\t") {
        return None;
    }

    let mut fields = line.split('\t');
    let record_type = fields.next()?;
    debug_assert_eq!(record_type, "L");
    let from_id = fields.next()?;
    fields.next()?;
    let to_id = fields.next()?;
    fields.next()?;
    let cigar = fields.next()?;
    Some((from_id, to_id, cigar))
}

#[inline]
fn parse_overlap_size(cigar: &str) -> Option<usize> {
    if cigar == "*" {
        return Some(0);
    }

    let bytes = cigar.as_bytes();
    let mut overlap = 0usize;
    let mut idx = 0usize;

    while idx < bytes.len() {
        let byte = bytes[idx];
        if byte.is_ascii_digit() {
            overlap = overlap
                .checked_mul(10)?
                .checked_add((byte - b'0') as usize)?;
            idx += 1;
        } else {
            break;
        }
    }

    if idx == 0 {
        return None;
    }

    match bytes.get(idx).copied() {
        Some(b'M') | Some(b'=') => Some(overlap),
        _ => None,
    }
}

#[inline]
fn decode_rle_tag(encoded: &str) -> Option<String> {
    if encoded.is_empty() {
        return Some(String::new());
    }

    let bytes = encoded.as_bytes();
    let mut idx = 0usize;
    let mut runs = Vec::new();
    let mut total_len = 0usize;

    while idx < bytes.len() {
        let base = *bytes.get(idx)?;
        if !base.is_ascii_alphabetic() {
            return None;
        }
        idx += 1;

        let digit_start = idx;
        let mut count = 0usize;
        while idx < bytes.len() && bytes[idx].is_ascii_digit() {
            count = count
                .checked_mul(10)?
                .checked_add((bytes[idx] - b'0') as usize)?;
            idx += 1;
        }

        if idx == digit_start || count == 0 {
            return None;
        }

        total_len = total_len.checked_add(count)?;
        runs.push((base, count));
    }

    let mut decoded = String::with_capacity(total_len);
    for (base, count) in runs {
        decoded.extend(std::iter::repeat(base as char).take(count));
    }
    Some(decoded)
}

fn collect_segment_id_map(gfa_path: &str) -> Result<HashMap<String, usize>> {
    let file = File::open(gfa_path)?;
    let reader = BufReader::new(file);
    let mut id_map = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let Some(segment) = parse_segment_line(&line) else {
            continue;
        };

        let fallback_idx = id_map.len();
        if let Entry::Vacant(slot) = id_map.entry(segment.id.to_string()) {
            slot.insert(fallback_idx);
        }
    }

    Ok(id_map)
}

impl GfaWriter {
    pub fn new(output_path: &str) -> Self {
        let file = File::create(output_path).expect("Could not create GFA file");
        Self {
            writer: BufWriter::new(file),
        }
    }

    /// Write segments (contigs)
    pub fn write_segments(&mut self, contigs: &[Contig]) -> Result<()> {
        // Write header
        writeln!(self.writer, "H\tVN:Z:1.0")?;

        for (i, contig) in contigs.iter().enumerate() {
            writeln!(self.writer, "S\tcontig_{}\t{}", i + 1, contig.sequence)?;
        }
        Ok(())
    }

    /// Write overlaps/links between contigs
    pub fn write_links(&mut self, links: &[(usize, usize, usize)]) -> Result<()> {
        for (from, to, overlap) in links {
            writeln!(
                self.writer,
                "L\tcontig_{}\t+\tcontig_{}\t+\t{}M",
                from + 1,
                to + 1,
                overlap
            )?;
        }
        Ok(())
    }

    /// Write paths for traversal visualization (using contig k-mer paths)
    pub fn write_paths(&mut self, contigs: &[Contig]) -> Result<()> {
        for (i, contig) in contigs.iter().enumerate() {
            let segments = contig
                .kmer_path
                .iter()
                .map(|_| format!("contig_{}", i + 1))
                .collect::<Vec<_>>()
                .join(",");
            writeln!(self.writer, "P\tpath_{}\t{}\t*", i + 1, segments)?;
        }
        Ok(())
    }

    /// Write assembly paths from Path objects
    pub fn write_assembly_paths(&mut self, paths: &[Path]) -> Result<()> {
        for path in paths {
            // Use our navigation module to get ODGI-style path representation
            let nav = traverse_path(path, false); // Don't include edges in GFA format
            let segments = nav.join(",");

            writeln!(self.writer, "P\tpath_{}\t{}\t*", path.id + 1, segments)?;
        }
        Ok(())
    }

    /// Write segments with RLE in tag field
    pub fn write_rle_segments(&mut self, contigs: &[Contig]) -> Result<()> {
        // Write header
        writeln!(self.writer, "H\tVN:Z:1.0")?;

        for (i, contig) in contigs.iter().enumerate() {
            let rle = crate::kmer::rle::rle_encode(&contig.sequence);
            let encoded = rle
                .iter()
                .map(|(b, c)| format!("{}{}", *b as char, c))
                .collect::<Vec<_>>()
                .join("");
            writeln!(self.writer, "S\tcontig_{}\t*\tRN:Z:{}", i + 1, encoded)?;
        }
        Ok(())
    }
}

/// Read contigs from a GFA file
pub fn read_gfa_contigs(gfa_path: &str) -> Result<Vec<Contig>> {
    let file = File::open(gfa_path)?;
    let reader = BufReader::new(file);
    let mut id_map = HashMap::new();
    let mut contigs = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let Some(segment) = parse_segment_line(&line) else {
            continue;
        };

        let sequence = if segment.sequence == "*" {
            segment.rle_tag.and_then(decode_rle_tag).unwrap_or_default()
        } else {
            segment.sequence.to_string()
        };

        let next_idx = id_map.len();
        if let Entry::Vacant(slot) = id_map.entry(segment.id.to_string()) {
            slot.insert(next_idx);
            contigs.push(Contig {
                id: next_idx,
                sequence,
                kmer_path: Vec::new(),
            });
        }
    }

    Ok(contigs)
}

/// Read links from a GFA file
pub fn read_gfa_links(gfa_path: &str) -> Result<Vec<(usize, usize, usize)>> {
    let id_map = collect_segment_id_map(gfa_path)?;
    if id_map.is_empty() {
        return Ok(Vec::new());
    }

    let mut links = Vec::new();

    // Read links and keep only those whose segment IDs exist.
    let file = File::open(gfa_path)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line?;
        let Some((from_id, to_id, cigar)) = parse_link_line(&line) else {
            continue;
        };

        let Some(from_idx) = id_map.get(from_id).copied() else {
            continue;
        };
        let Some(to_idx) = id_map.get(to_id).copied() else {
            continue;
        };
        let Some(overlap_size) = parse_overlap_size(cigar) else {
            continue;
        };

        links.push((from_idx, to_idx, overlap_size));
    }

    // Keep link output deterministic regardless of record order in input GFA.
    links.sort_unstable_by(|a, b| {
        a.0.cmp(&b.0)
            .then_with(|| a.1.cmp(&b.1))
            .then_with(|| a.2.cmp(&b.2))
    });

    Ok(links)
}

#[cfg(test)]
mod tests {
    use super::{decode_rle_tag, parse_overlap_size, read_gfa_contigs, read_gfa_links, GfaWriter};
    use crate::graph::assembler::Contig;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn read_gfa_contigs_assigns_dense_ids_for_malformed_and_sparse_ids() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "H\tVN:Z:1.0").unwrap();
        writeln!(file, "S\tcontig_x\tAAAA").unwrap();
        writeln!(file, "S\tcontig_0\tCCCC").unwrap();
        writeln!(file, "S\tcontig_42\tGGGG").unwrap();

        let contigs = read_gfa_contigs(file.path().to_str().unwrap()).unwrap();
        assert_eq!(contigs.len(), 3);
        assert_eq!(contigs[0].id, 0);
        assert_eq!(contigs[1].id, 1);
        assert_eq!(contigs[2].id, 2);
    }

    #[test]
    fn read_gfa_links_handles_malformed_contig_ids_deterministically() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "H\tVN:Z:1.0").unwrap();
        writeln!(file, "S\tcontig_x\tAAAA").unwrap();
        writeln!(file, "S\tcustom\tCCCC").unwrap();
        writeln!(file, "S\tcontig_2\tGGGG").unwrap();
        writeln!(file, "L\tcontig_x\t+\tcustom\t+\t3M").unwrap();
        writeln!(file, "L\tcontig_2\t+\tcontig_x\t+\t2M").unwrap();

        let links = read_gfa_links(file.path().to_str().unwrap()).unwrap();
        assert_eq!(links, vec![(0, 1, 3), (2, 0, 2)]);
    }

    #[test]
    fn read_gfa_links_ignores_unknown_segments_and_invalid_cigar() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "H\tVN:Z:1.0").unwrap();
        writeln!(file, "S\tseg_a\tAAAA").unwrap();
        writeln!(file, "S\tseg_b\tCCCC").unwrap();
        writeln!(file, "L\tseg_a\t+\tseg_b\t+\t12M1I").unwrap();
        writeln!(file, "L\tseg_a\t+\tmissing\t+\t5M").unwrap();
        writeln!(file, "L\tmissing\t+\tseg_b\t+\t5M").unwrap();
        writeln!(file, "L\tseg_b\t+\tseg_a\t+\tZZM").unwrap();

        let links = read_gfa_links(file.path().to_str().unwrap()).unwrap();
        assert_eq!(links, vec![(0, 1, 12)]);
    }

    #[test]
    fn read_gfa_links_is_stable_under_link_record_reordering() {
        let mut file_a = NamedTempFile::new().unwrap();
        writeln!(file_a, "H\tVN:Z:1.0").unwrap();
        writeln!(file_a, "S\tseg_a\tAAAA").unwrap();
        writeln!(file_a, "S\tseg_b\tCCCC").unwrap();
        writeln!(file_a, "S\tseg_c\tGGGG").unwrap();
        writeln!(file_a, "L\tseg_c\t+\tseg_a\t+\t4M").unwrap();
        writeln!(file_a, "L\tseg_a\t+\tseg_b\t+\t2M").unwrap();
        writeln!(file_a, "L\tseg_a\t+\tseg_c\t+\t3M").unwrap();

        let mut file_b = NamedTempFile::new().unwrap();
        writeln!(file_b, "H\tVN:Z:1.0").unwrap();
        writeln!(file_b, "S\tseg_a\tAAAA").unwrap();
        writeln!(file_b, "S\tseg_b\tCCCC").unwrap();
        writeln!(file_b, "S\tseg_c\tGGGG").unwrap();
        writeln!(file_b, "L\tseg_a\t+\tseg_c\t+\t3M").unwrap();
        writeln!(file_b, "L\tseg_c\t+\tseg_a\t+\t4M").unwrap();
        writeln!(file_b, "L\tseg_a\t+\tseg_b\t+\t2M").unwrap();

        let links_a = read_gfa_links(file_a.path().to_str().unwrap()).unwrap();
        let links_b = read_gfa_links(file_b.path().to_str().unwrap()).unwrap();
        assert_eq!(links_a, links_b);
        assert_eq!(links_a, vec![(0, 1, 2), (0, 2, 3), (2, 0, 4)]);
    }

    #[test]
    fn read_gfa_contigs_deduplicates_segments_by_first_occurrence() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "H\tVN:Z:1.0").unwrap();
        writeln!(file, "S\tdup\tAAAA").unwrap();
        writeln!(file, "S\tdup\tCCCC").unwrap();
        writeln!(file, "S\tother\tGGGG").unwrap();

        let contigs = read_gfa_contigs(file.path().to_str().unwrap()).unwrap();
        assert_eq!(contigs.len(), 2);
        assert_eq!(contigs[0].id, 0);
        assert_eq!(contigs[0].sequence, "AAAA");
        assert_eq!(contigs[1].id, 1);
        assert_eq!(contigs[1].sequence, "GGGG");
    }

    #[test]
    fn parse_overlap_size_requires_match_operator_after_length_prefix() {
        assert_eq!(parse_overlap_size("*"), Some(0));
        assert_eq!(parse_overlap_size("12M"), Some(12));
        assert_eq!(parse_overlap_size("12M1I"), Some(12));
        assert_eq!(parse_overlap_size("7=2X"), Some(7));

        assert_eq!(parse_overlap_size("12"), None);
        assert_eq!(parse_overlap_size("12I"), None);
        assert_eq!(parse_overlap_size("12S10M"), None);
        assert_eq!(parse_overlap_size("M12"), None);
    }

    #[test]
    fn decode_rle_tag_parses_multi_digit_runs() {
        assert_eq!(
            decode_rle_tag("A4C3T12"),
            Some("AAAACCCTTTTTTTTTTTT".to_string())
        );
        assert_eq!(decode_rle_tag("N1"), Some("N".to_string()));
        assert_eq!(decode_rle_tag(""), Some(String::new()));
    }

    #[test]
    fn decode_rle_tag_rejects_malformed_inputs() {
        assert_eq!(decode_rle_tag("A"), None);
        assert_eq!(decode_rle_tag("A0"), None);
        assert_eq!(decode_rle_tag("3A"), None);
        assert_eq!(decode_rle_tag("A2-"), None);
    }

    #[test]
    fn read_gfa_contigs_decodes_rle_sequence_tag_when_sequence_is_unknown() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "H\tVN:Z:1.0").unwrap();
        writeln!(file, "S\tseg_1\t*\tRN:Z:A4C2T1").unwrap();
        writeln!(file, "S\tseg_2\tACGT").unwrap();

        let contigs = read_gfa_contigs(file.path().to_str().unwrap()).unwrap();
        assert_eq!(contigs.len(), 2);
        assert_eq!(contigs[0].sequence, "AAAACCT");
        assert_eq!(contigs[1].sequence, "ACGT");
    }

    #[test]
    fn read_gfa_contigs_rle_writer_round_trip_preserves_sequences() {
        let contigs = vec![
            Contig {
                id: 10,
                sequence: "AAAACCCCG".to_string(),
                kmer_path: Vec::new(),
            },
            Contig {
                id: 42,
                sequence: "TTTTGG".to_string(),
                kmer_path: Vec::new(),
            },
            Contig {
                id: 100,
                sequence: String::new(),
                kmer_path: Vec::new(),
            },
        ];

        let file = NamedTempFile::new().unwrap();
        let path = file.path().to_str().unwrap();
        let mut writer = GfaWriter::new(path);
        writer.write_rle_segments(&contigs).unwrap();
        drop(writer);

        let observed = read_gfa_contigs(path).unwrap();
        let observed_sequences: Vec<String> = observed.into_iter().map(|c| c.sequence).collect();
        let expected_sequences: Vec<String> = contigs.into_iter().map(|c| c.sequence).collect();
        assert_eq!(observed_sequences, expected_sequences);
    }
}
