use std::fs::File;
use std::io::{BufReader, BufRead, BufWriter, Write, Result};
use std::collections::HashMap;
use crate::graph::assembler::Contig;
use crate::graph::stitch::Path;
use crate::graph::navigation::traverse_path;

pub struct GfaWriter {
    writer: BufWriter<File>,
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
            writeln!(self.writer, "L\tcontig_{}\t+\tcontig_{}\t+\t{}M", from + 1, to + 1, overlap)?;
        }
        Ok(())
    }
    
    /// Write paths for traversal visualization (using contig k-mer paths)
    pub fn write_paths(&mut self, contigs: &[Contig]) -> Result<()> {
        for (i, contig) in contigs.iter().enumerate() {
            let segments = contig.kmer_path.iter()
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
            let encoded = rle.iter()
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
            let id = if id_str.starts_with("contig_") {
                id_str[7..].parse::<usize>().unwrap_or(contigs.len()) - 1
            } else {
                contigs.len()
            };
            
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
            if !id_str.starts_with("contig_") {
                id_map.insert(id_str.to_string(), id_map.len());
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
                .trim_end_matches(|c| c == 'M' || c == 'm')
                .parse::<usize>()
                .unwrap_or(0);
            
            // Convert IDs to numeric indices
            let from_idx = if from_id.starts_with("contig_") {
                from_id[7..].parse::<usize>().unwrap_or(0) - 1
            } else {
                *id_map.get(from_id).unwrap_or(&0)
            };
            
            let to_idx = if to_id.starts_with("contig_") {
                to_id[7..].parse::<usize>().unwrap_or(0) - 1
            } else {
                *id_map.get(to_id).unwrap_or(&0)
            };
            
            links.push((from_idx, to_idx, overlap_size));
        }
    }
    
    Ok(links)
}
