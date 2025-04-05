use std::fs::File;
use std::io::{BufWriter, Write, Result};
use crate::graph::assembler::Contig;
use crate::graph::stitch::Path;
use crate::kmer::rle::rle_encode;
use crate::graph::navigation::traverse_path;
use std::collections::HashMap;

pub struct Gfa2Writer {
    writer: BufWriter<File>,
    coverage: Option<HashMap<usize, f32>>,
}

impl Gfa2Writer {
    pub fn new(path: &str) -> Self {
        let file = File::create(path).expect("Failed to create GFA2 file");
        let mut writer = BufWriter::new(file);
        writeln!(writer, "H\tVN:Z:2.0").unwrap();
        Self { 
            writer,
            coverage: None,
        }
    }

    pub fn with_coverage(path: &str, coverage: HashMap<usize, f32>) -> Self {
        let file = File::create(path).expect("Failed to create GFA2 file");
        let mut writer = BufWriter::new(file);
        writeln!(writer, "H\tVN:Z:2.0").unwrap();
        Self { 
            writer,
            coverage: Some(coverage),
        }
    }

    /// Annotate segment with coverage + compression
    pub fn annotate_segment(
        &mut self,
        name: &str,
        sequence: &str,
        coverage: Option<f32>,
    ) -> Result<()> {
        let rle = rle_encode(sequence);
        let rle_ratio = if sequence.len() > 0 {
            1.0 - (rle.len() as f32 / sequence.len() as f32)
        } else {
            0.0
        };
        let cov = coverage.unwrap_or(1.0);

        // Determine color based on RLE ratio for BandageNG
        let color = if rle_ratio > 0.7 {
            "red"
        } else if rle_ratio > 0.4 {
            "orange"
        } else {
            "green"
        };

        // Write segment with additional tags for ODGI/BandageNG visualization
        writeln!(
            self.writer,
            "S\t{}\t{}\t*\tRC:f:{:.3}\tCV:f:{:.2}\tCL:Z:{}",
            name, sequence.len(), rle_ratio, cov, color
        )?;
        writeln!(self.writer, "a\t{}\tseq\t{}", name, sequence)?;

        Ok(())
    }

    pub fn write_segments(&mut self, contigs: &[Contig]) -> Result<()> {
        for (i, contig) in contigs.iter().enumerate() {
            let contig_id = i + 1;
            let name = format!("contig_{}", contig_id);
            
            // Get coverage for this contig if available
            let coverage = self.coverage.as_ref()
                .and_then(|cov_map| cov_map.get(&contig.id).or_else(|| cov_map.get(&contig_id)).copied());
            
            self.annotate_segment(&name, &contig.sequence, coverage)?;
        }
        Ok(())
    }

    pub fn write_links(&mut self, links: &[(usize, usize, usize)]) -> Result<()> {
        for (from, to, overlap) in links {
            writeln!(self.writer, "E\tedge_{}_{}\tcontig_{}\t+\tcontig_{}\t+\t0\t{}\t0\t{}\t{}M",
                     from + 1, to + 1, from + 1, to + 1, overlap, overlap, overlap)?;
        }
        Ok(())
    }

    pub fn write_paths(&mut self, paths: &[Path]) -> Result<()> {
        for path in paths {
            // Generate path using the navigation module, with -> between nodes
            let nav = traverse_path(path, true);
            let segs = nav.join(" ");
            
            // Calculate the total path length including overlaps
            let _total_length: usize = path.segments.iter()
                .zip(path.overlaps.iter().chain(std::iter::once(&0)))
                .map(|(_, overlap)| *overlap)
                .sum();
                
            writeln!(self.writer, "O\tpath_{}\t{}\t*", path.id + 1, segs)?;
        }
        Ok(())
    }
} 