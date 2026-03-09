//! Contig Polisher - Correct errors in assembled contigs using raw reads
//!
//! Algorithm:
//! 1. Index contigs using minimizers
//! 2. Map reads to contigs
//! 3. Build pileup at each position
//! 4. Call consensus based on base frequencies
//! 5. Output polished contigs

use crate::io::fasta::{open_fasta, FastaWriter};
use crate::io::fastq::{open_fastq, stream_fastq_records};
use ahash::AHashMap;
use std::io::{BufRead, Result};
use tracing::info;

/// Statistics from polishing
#[derive(Debug, Clone, Default)]
pub struct PolishStats {
    pub contigs_input: usize,
    pub reads_processed: u64,
    pub reads_mapped: u64,
    pub corrections: usize,
    pub insertions: usize,
    pub deletions: usize,
    pub substitutions: usize,
}

/// Pileup column for consensus calling
#[derive(Debug, Default, Clone)]
struct PileupColumn {
    bases: [u32; 5], // A, C, G, T, gap
    depth: u32,
}

impl PileupColumn {
    fn add_base(&mut self, base: u8) {
        let idx = match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 4, // gap or N
        };
        self.bases[idx] += 1;
        self.depth += 1;
    }

    fn consensus(&self, min_freq: f64) -> Option<u8> {
        if self.depth == 0 {
            return None;
        }

        let max_idx = self
            .bases
            .iter()
            .enumerate()
            .max_by_key(|(_, &count)| count)
            .map(|(idx, _)| idx)?;

        let freq = self.bases[max_idx] as f64 / self.depth as f64;
        if freq >= min_freq {
            match max_idx {
                0 => Some(b'A'),
                1 => Some(b'C'),
                2 => Some(b'G'),
                3 => Some(b'T'),
                _ => None, // Gap - deletion
            }
        } else {
            None
        }
    }
}

/// Minimizer index for read mapping
struct MinimizerIndex {
    index: AHashMap<u64, Vec<(usize, usize)>>, // minimizer -> [(contig_id, position)]
    k: usize,
    w: usize,
}

impl MinimizerIndex {
    fn new(k: usize, w: usize) -> Self {
        Self {
            index: AHashMap::new(),
            k,
            w,
        }
    }

    fn add_contig(&mut self, contig_id: usize, sequence: &[u8]) {
        for (pos, minimizer) in self.extract_minimizers(sequence) {
            self.index
                .entry(minimizer)
                .or_default()
                .push((contig_id, pos));
        }
    }

    fn extract_minimizers(&self, seq: &[u8]) -> Vec<(usize, u64)> {
        let mut minimizers = Vec::new();
        if seq.len() < self.k + self.w - 1 {
            return minimizers;
        }

        for window_start in 0..=(seq.len() - self.k - self.w + 1) {
            let mut min_hash = u64::MAX;
            let mut min_pos = 0;

            for i in 0..self.w {
                let pos = window_start + i;
                if pos + self.k <= seq.len() {
                    let hash = hash_kmer(&seq[pos..pos + self.k]);
                    if hash < min_hash {
                        min_hash = hash;
                        min_pos = pos;
                    }
                }
            }

            if minimizers.is_empty() || minimizers.last().map(|(p, _)| *p) != Some(min_pos) {
                minimizers.push((min_pos, min_hash));
            }
        }

        minimizers
    }

    fn map_read(&self, sequence: &[u8]) -> Option<(usize, usize, bool)> {
        let mut hits: AHashMap<(usize, usize), usize> = AHashMap::new();

        // Forward mapping
        for (read_pos, minimizer) in self.extract_minimizers(sequence) {
            if let Some(entries) = self.index.get(&minimizer) {
                for &(contig_id, contig_pos) in entries {
                    let start = if contig_pos >= read_pos {
                        contig_pos - read_pos
                    } else {
                        0
                    };
                    *hits.entry((contig_id, start)).or_default() += 1;
                }
            }
        }

        // Reverse complement mapping
        let rc = reverse_complement(sequence);
        for (read_pos, minimizer) in self.extract_minimizers(&rc) {
            if let Some(entries) = self.index.get(&minimizer) {
                for &(contig_id, contig_pos) in entries {
                    let start = if contig_pos >= read_pos {
                        contig_pos - read_pos
                    } else {
                        0
                    };
                    // Use negative positions to mark reverse mappings (we'll handle this later)
                    *hits.entry((contig_id, start | 0x80000000)).or_default() += 1;
                }
            }
        }

        // Find best hit with at least 3 minimizer matches
        hits.into_iter()
            .filter(|(_, count)| *count >= 3)
            .max_by_key(|(_, count)| *count)
            .map(|((contig_id, pos), _)| {
                let is_reverse = (pos & 0x80000000) != 0;
                let actual_pos = pos & 0x7FFFFFFF;
                (contig_id, actual_pos, is_reverse)
            })
    }
}

fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut hash = 0u64;
    for &base in kmer {
        let val = match base {
            b'A' | b'a' => 0u64,
            b'C' | b'c' => 1u64,
            b'G' | b'g' => 2u64,
            b'T' | b't' => 3u64,
            _ => continue,
        };
        hash = hash.wrapping_mul(4).wrapping_add(val);
    }
    hash
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        })
        .collect()
}

/// Read contigs from FASTA file
fn read_contigs(path: &str) -> Result<Vec<(String, Vec<u8>)>> {
    let reader = open_fasta(path);
    let mut contigs = Vec::new();
    let mut current_header = String::new();
    let mut current_seq = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_header.is_empty() {
                contigs.push((current_header.clone(), current_seq.clone()));
            }
            current_header = line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            current_seq.clear();
        } else {
            current_seq.extend(line.trim().bytes());
        }
    }

    if !current_header.is_empty() {
        contigs.push((current_header, current_seq));
    }

    Ok(contigs)
}

/// Main polishing function
pub fn polish_contigs(
    contigs_path: &str,
    reads1_path: &str,
    reads2_path: Option<&str>,
    output_path: &str,
    iterations: usize,
) -> Result<PolishStats> {
    let mut stats = PolishStats::default();

    // Read initial contigs
    let mut contigs = read_contigs(contigs_path)?;
    stats.contigs_input = contigs.len();
    info!("Loaded {} contigs for polishing", contigs.len());

    for iter in 0..iterations {
        info!("Polishing iteration {}/{}", iter + 1, iterations);

        // Build minimizer index
        let mut index = MinimizerIndex::new(15, 10);
        for (i, (_header, seq)) in contigs.iter().enumerate() {
            index.add_contig(i, seq);
        }

        // Build pileups for each contig
        let mut pileups: Vec<Vec<PileupColumn>> = contigs
            .iter()
            .map(|(_, seq)| vec![PileupColumn::default(); seq.len()])
            .collect();

        // Process reads from first file
        info!("  Mapping reads from {}...", reads1_path);
        let reader1 = open_fastq(reads1_path);
        for record in stream_fastq_records(reader1) {
            stats.reads_processed += 1;
            if let Some((contig_id, start, is_reverse)) = index.map_read(record.sequence.as_bytes())
            {
                stats.reads_mapped += 1;
                let read_seq = if is_reverse {
                    reverse_complement(record.sequence.as_bytes())
                } else {
                    record.sequence.into_bytes()
                };
                add_to_pileup(&mut pileups[contig_id], start, &read_seq);
            }
        }

        // Process reads from second file if provided
        if let Some(reads2) = reads2_path {
            info!("  Mapping reads from {}...", reads2);
            let reader2 = open_fastq(reads2);
            for record in stream_fastq_records(reader2) {
                stats.reads_processed += 1;
                if let Some((contig_id, start, is_reverse)) =
                    index.map_read(record.sequence.as_bytes())
                {
                    stats.reads_mapped += 1;
                    let read_seq = if is_reverse {
                        reverse_complement(record.sequence.as_bytes())
                    } else {
                        record.sequence.into_bytes()
                    };
                    add_to_pileup(&mut pileups[contig_id], start, &read_seq);
                }
            }
        }

        // Call consensus and apply corrections
        info!("  Calling consensus...");
        let iter_corrections = apply_corrections(&mut contigs, &pileups, &mut stats);
        info!("  Made {} corrections", iter_corrections);

        if iter_corrections == 0 {
            info!("  No more corrections needed, stopping early");
            break;
        }
    }

    // Write polished contigs
    info!("Writing polished contigs to {}", output_path);
    let mut writer = FastaWriter::new(output_path);
    for (header, seq) in &contigs {
        let seq_str = String::from_utf8_lossy(seq);
        writer.write_record(header, &seq_str)?;
    }

    info!(
        "Polishing complete: {} corrections ({} substitutions, {} insertions, {} deletions)",
        stats.corrections, stats.substitutions, stats.insertions, stats.deletions
    );

    Ok(stats)
}

/// Add a read alignment to the pileup
fn add_to_pileup(pileup: &mut [PileupColumn], start: usize, read: &[u8]) {
    for (i, &base) in read.iter().enumerate() {
        let pos = start + i;
        if pos < pileup.len() {
            pileup[pos].add_base(base);
        }
    }
}

/// Apply corrections based on pileup consensus
fn apply_corrections(
    contigs: &mut [(String, Vec<u8>)],
    pileups: &[Vec<PileupColumn>],
    stats: &mut PolishStats,
) -> usize {
    let mut total_corrections = 0;
    const MIN_DEPTH: u32 = 5;
    const MIN_FREQ: f64 = 0.6;

    for (contig_idx, (_, seq)) in contigs.iter_mut().enumerate() {
        let pileup = &pileups[contig_idx];
        let mut corrections = Vec::new();

        for (pos, column) in pileup.iter().enumerate() {
            if column.depth < MIN_DEPTH {
                continue;
            }

            if let Some(consensus_base) = column.consensus(MIN_FREQ) {
                if pos < seq.len() && seq[pos] != consensus_base {
                    corrections.push((pos, consensus_base));
                }
            }
        }

        // Apply corrections (in reverse order to preserve positions)
        for (pos, base) in corrections.into_iter().rev() {
            if pos < seq.len() {
                let old_base = seq[pos];
                seq[pos] = base;
                total_corrections += 1;
                stats.corrections += 1;

                // Categorize the correction
                if old_base == b'-' || old_base == b'N' {
                    stats.insertions += 1;
                } else if base == b'-' {
                    stats.deletions += 1;
                } else {
                    stats.substitutions += 1;
                }
            }
        }
    }

    total_corrections
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pileup_column() {
        let mut col = PileupColumn::default();
        col.add_base(b'A');
        col.add_base(b'A');
        col.add_base(b'A');
        col.add_base(b'C');

        assert_eq!(col.depth, 4);
        assert_eq!(col.consensus(0.5), Some(b'A'));
        assert_eq!(col.consensus(0.9), None);
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
    }

    #[test]
    fn test_minimizer_index() {
        let mut index = MinimizerIndex::new(11, 5);
        index.add_contig(0, b"ACGTACGTACGTACGTACGTACGTACGTACGT");

        // Should find the contig
        let result = index.map_read(b"ACGTACGTACGTACGTACGT");
        assert!(result.is_some());
    }
}
