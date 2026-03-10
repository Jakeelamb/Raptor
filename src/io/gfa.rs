use crate::graph::assembler::Contig;
use crate::graph::navigation::traverse_path;
use crate::graph::stitch::Path;
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Result, Write};

pub struct GfaWriter {
    writer: BufWriter<File>,
}

#[inline]
fn parse_contig_index(id: &str) -> Option<usize> {
    let raw = id.strip_prefix("contig_")?;
    let one_based = raw.parse::<usize>().ok()?;
    one_based.checked_sub(1)
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
    let mut contigs = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('S') {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 3 {
                continue;
            }

            let id_str = parts[1];
            let id = parse_contig_index(id_str).unwrap_or(contigs.len());

            let sequence = parts[2].to_string();

            // Create a contig with empty kmer_path for now
            let contig = Contig {
                id,
                sequence,
                kmer_path: Vec::new(),
            };

            contigs.push(contig);
        }
    }

    Ok(contigs)
}

/// Read links from a GFA file
pub fn read_gfa_links(gfa_path: &str) -> Result<Vec<(usize, usize, usize)>> {
    let file = File::open(gfa_path)?;
    let reader = BufReader::new(file);
    let mut links = Vec::new();
    let mut id_map = HashMap::new();

    // First pass: build id mapping if needed
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('S') {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 3 {
                continue;
            }

            let id_str = parts[1];
            if parse_contig_index(id_str).is_none() {
                let fallback_idx = id_map.len();
                if let Entry::Vacant(slot) = id_map.entry(id_str.to_string()) {
                    slot.insert(fallback_idx);
                }
            }
        }
    }

    // Second pass: read links
    let file = File::open(gfa_path)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('L') {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 6 {
                continue;
            }

            let from_id = parts[1];
            let to_id = parts[3];

            // Extract overlap size from CIGAR (format like "10M")
            let cigar = parts[5];
            let overlap_size = cigar
                .trim_end_matches(['M', 'm'])
                .parse::<usize>()
                .unwrap_or(0);

            // Convert IDs to numeric indices
            let from_idx = parse_contig_index(from_id)
                .or_else(|| id_map.get(from_id).copied())
                .unwrap_or(0);
            let to_idx = parse_contig_index(to_id)
                .or_else(|| id_map.get(to_id).copied())
                .unwrap_or(0);

            links.push((from_idx, to_idx, overlap_size));
        }
    }

    Ok(links)
}

#[cfg(test)]
mod tests {
    use super::{read_gfa_contigs, read_gfa_links};
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn read_gfa_contigs_handles_malformed_contig_ids_without_underflow() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "H\tVN:Z:1.0").unwrap();
        writeln!(file, "S\tcontig_x\tAAAA").unwrap();
        writeln!(file, "S\tcontig_0\tCCCC").unwrap();
        writeln!(file, "S\tcontig_2\tGGGG").unwrap();

        let contigs = read_gfa_contigs(file.path().to_str().unwrap()).unwrap();
        assert_eq!(contigs.len(), 3);
        assert_eq!(contigs[0].id, 0);
        assert_eq!(contigs[1].id, 1);
        assert_eq!(contigs[2].id, 1);
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
        assert_eq!(links, vec![(0, 1, 3), (1, 0, 2)]);
    }
}
