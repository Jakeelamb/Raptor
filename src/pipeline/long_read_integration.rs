//! Long Read Integration - Hybrid assembly with PacBio/Nanopore reads
//!
//! Strategy:
//! 1. Use short reads to build accurate de Bruijn graph
//! 2. Map long reads to contigs using minimizers
//! 3. Use long reads to span repeat regions
//! 4. Anchor long reads to extend/connect contigs

use crate::io::fasta::{open_fasta, FastaWriter};
use crate::io::fastq::{open_fastq, stream_fastq_records};
use ahash::{AHashMap, AHashSet};
use std::io::{BufRead, Result};
use tracing::info;

/// Statistics from long read integration
#[derive(Debug, Clone, Default)]
pub struct LongReadStats {
    pub contigs_input: usize,
    pub long_reads_processed: u64,
    pub long_reads_mapped: u64,
    pub contigs_extended: usize,
    pub contigs_joined: usize,
    pub total_extension_bp: usize,
    pub output_contigs: usize,
}

/// Configuration for long read integration
#[derive(Clone, Debug)]
pub struct LongReadConfig {
    /// Minimum minimizer matches to anchor a long read
    pub min_anchor_matches: usize,
    /// Minimum overlap between long read and contig
    pub min_overlap: usize,
    /// Minimum long read length to use
    pub min_read_length: usize,
    /// K-mer size for minimizers
    pub minimizer_k: usize,
    /// Window size for minimizers
    pub minimizer_w: usize,
}

impl Default for LongReadConfig {
    fn default() -> Self {
        Self {
            min_anchor_matches: 5,
            min_overlap: 500,
            min_read_length: 1000,
            minimizer_k: 15,
            minimizer_w: 10,
        }
    }
}

/// A mapped region of a long read to a contig
#[derive(Debug, Clone)]
struct ReadMapping {
    contig_id: usize,
    contig_start: usize,
    contig_end: usize,
    read_start: usize,
    read_end: usize,
    is_reverse: bool,
    score: usize, // Number of matching minimizers
}

/// Minimizer index for fast contig mapping
struct MinimizerIndex {
    index: AHashMap<u64, Vec<(usize, usize, bool)>>, // minimizer -> [(contig_id, pos, is_rc)]
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
        // Forward strand
        for (pos, minimizer) in self.extract_minimizers(sequence) {
            self.index
                .entry(minimizer)
                .or_default()
                .push((contig_id, pos, false));
        }

        // Reverse complement
        let rc = reverse_complement(sequence);
        for (pos, minimizer) in self.extract_minimizers(&rc) {
            let orig_pos = sequence.len().saturating_sub(pos + self.k);
            self.index
                .entry(minimizer)
                .or_default()
                .push((contig_id, orig_pos, true));
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

    /// Map a long read to contigs
    fn map_read(&self, sequence: &[u8]) -> Vec<ReadMapping> {
        let mut hits: AHashMap<(usize, bool), Vec<(usize, usize)>> = AHashMap::new();

        // Collect minimizer hits
        for (read_pos, minimizer) in self.extract_minimizers(sequence) {
            if let Some(entries) = self.index.get(&minimizer) {
                for &(contig_id, contig_pos, is_rc) in entries {
                    hits.entry((contig_id, is_rc))
                        .or_default()
                        .push((read_pos, contig_pos));
                }
            }
        }

        // Convert hits to mappings
        let mut mappings = Vec::new();
        for ((contig_id, is_reverse), positions) in hits {
            if positions.len() < 3 {
                continue;
            }

            // Find the range of the mapping
            let read_start = positions.iter().map(|(r, _)| *r).min().unwrap_or(0);
            let read_end = positions.iter().map(|(r, _)| *r).max().unwrap_or(0) + self.k;
            let contig_start = positions.iter().map(|(_, c)| *c).min().unwrap_or(0);
            let contig_end = positions.iter().map(|(_, c)| *c).max().unwrap_or(0) + self.k;

            mappings.push(ReadMapping {
                contig_id,
                contig_start,
                contig_end,
                read_start,
                read_end,
                is_reverse,
                score: positions.len(),
            });
        }

        // Sort by score
        mappings.sort_by(|a, b| b.score.cmp(&a.score));
        mappings
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

/// Read contigs from FASTA
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

/// A link between two contigs discovered from long reads
#[derive(Debug, Clone)]
struct ContigLink {
    contig_a: usize,
    contig_b: usize,
    orientation_a: bool, // false = forward, true = reverse
    orientation_b: bool,
    gap_estimate: i32,
    supporting_reads: usize,
    bridging_sequences: Vec<Vec<u8>>,
}

/// Main function for long read integration
pub fn integrate_long_reads(
    contigs_path: &str,
    long_reads_path: &str,
    output_path: &str,
    config: LongReadConfig,
) -> Result<LongReadStats> {
    let mut stats = LongReadStats::default();

    info!("Loading contigs from {}", contigs_path);
    let mut contigs = read_contigs(contigs_path)?;
    stats.contigs_input = contigs.len();
    info!("Loaded {} contigs", contigs.len());

    if contigs.is_empty() {
        return Ok(stats);
    }

    // Build minimizer index
    info!("Building minimizer index...");
    let mut index = MinimizerIndex::new(config.minimizer_k, config.minimizer_w);
    for (i, (_header, seq)) in contigs.iter().enumerate() {
        index.add_contig(i, seq);
    }

    // Process long reads and find links
    info!("Processing long reads from {}...", long_reads_path);
    let mut links: AHashMap<(usize, usize), ContigLink> = AHashMap::new();
    let mut extensions: AHashMap<usize, Vec<(bool, Vec<u8>)>> = AHashMap::new(); // contig_id -> [(is_right, extension_seq)]

    let reader = open_fastq(long_reads_path);
    for record in stream_fastq_records(reader) {
        stats.long_reads_processed += 1;

        if record.sequence.len() < config.min_read_length {
            continue;
        }

        let mappings = index.map_read(record.sequence.as_bytes());

        if mappings.is_empty() {
            continue;
        }

        stats.long_reads_mapped += 1;

        // Check for read spanning multiple contigs (potential links)
        if mappings.len() >= 2 {
            let best = &mappings[0];
            let second = &mappings[1];

            if best.score >= config.min_anchor_matches && second.score >= config.min_anchor_matches
            {
                let key = if best.contig_id < second.contig_id {
                    (best.contig_id, second.contig_id)
                } else {
                    (second.contig_id, best.contig_id)
                };

                let link = links.entry(key).or_insert_with(|| ContigLink {
                    contig_a: key.0,
                    contig_b: key.1,
                    orientation_a: best.is_reverse,
                    orientation_b: second.is_reverse,
                    gap_estimate: 0,
                    supporting_reads: 0,
                    bridging_sequences: Vec::new(),
                });

                link.supporting_reads += 1;

                // Extract bridging sequence if available
                if best.read_end < second.read_start && second.read_start - best.read_end < 10000 {
                    let bridge = record.sequence[best.read_end..second.read_start]
                        .as_bytes()
                        .to_vec();
                    if bridge.len() < 5000 {
                        link.bridging_sequences.push(bridge);
                    }
                }
            }
        }

        // Check for contig extensions (read extends past contig end)
        if let Some(best) = mappings.first() {
            if best.score >= config.min_anchor_matches {
                let contig_len = contigs[best.contig_id].1.len();

                // Check for left extension
                if best.contig_start < 50 && best.read_start > config.min_overlap {
                    let extension = record.sequence[..best.read_start].as_bytes().to_vec();
                    extensions
                        .entry(best.contig_id)
                        .or_default()
                        .push((false, extension));
                }

                // Check for right extension
                if best.contig_end > contig_len - 50
                    && record.sequence.len() - best.read_end > config.min_overlap
                {
                    let extension = record.sequence[best.read_end..].as_bytes().to_vec();
                    extensions
                        .entry(best.contig_id)
                        .or_default()
                        .push((true, extension));
                }
            }
        }

        if stats.long_reads_processed % 100_000 == 0 {
            info!(
                "  Processed {} long reads, {} mapped...",
                stats.long_reads_processed, stats.long_reads_mapped
            );
        }
    }

    // Apply extensions to contigs
    info!("Applying {} contig extensions...", extensions.len());
    for (contig_id, exts) in &extensions {
        let contig_seq = &mut contigs[*contig_id].1;

        // Find consensus extension for each end
        for is_right in [false, true] {
            let relevant: Vec<&Vec<u8>> = exts
                .iter()
                .filter(|(r, _)| *r == is_right)
                .map(|(_, seq)| seq)
                .collect();

            if relevant.len() >= 2 {
                // Simple consensus: use the most common sequence or the longest
                if let Some(best_ext) = relevant.iter().max_by_key(|s| s.len()) {
                    if is_right {
                        contig_seq.extend_from_slice(best_ext);
                    } else {
                        let mut new_seq = (*best_ext).clone();
                        new_seq.extend_from_slice(contig_seq);
                        *contig_seq = new_seq;
                    }
                    stats.contigs_extended += 1;
                    stats.total_extension_bp += best_ext.len();
                }
            }
        }
    }

    // Join contigs using links with bridging sequences
    info!("Processing {} potential contig links...", links.len());
    let mut joined_contigs: AHashSet<usize> = AHashSet::new();
    let mut new_contigs: Vec<(String, Vec<u8>)> = Vec::new();

    let mut sorted_links: Vec<_> = links
        .values()
        .filter(|l| l.supporting_reads >= 3 && !l.bridging_sequences.is_empty())
        .collect();
    sorted_links.sort_by(|a, b| b.supporting_reads.cmp(&a.supporting_reads));

    for link in sorted_links {
        if joined_contigs.contains(&link.contig_a) || joined_contigs.contains(&link.contig_b) {
            continue;
        }

        // Use the most common bridging sequence
        if let Some(bridge) = link.bridging_sequences.first() {
            let mut joined = contigs[link.contig_a].1.clone();
            joined.extend_from_slice(bridge);
            joined.extend_from_slice(&contigs[link.contig_b].1);

            let header = format!(
                "joined_{}_{}",
                contigs[link.contig_a].0, contigs[link.contig_b].0
            );
            new_contigs.push((header, joined));

            joined_contigs.insert(link.contig_a);
            joined_contigs.insert(link.contig_b);
            stats.contigs_joined += 2;
        }
    }

    // Collect final contigs
    let mut final_contigs: Vec<(String, Vec<u8>)> = contigs
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !joined_contigs.contains(i))
        .map(|(_, c)| c)
        .collect();
    final_contigs.extend(new_contigs);

    stats.output_contigs = final_contigs.len();

    // Write output
    info!("Writing {} contigs to {}", final_contigs.len(), output_path);
    let mut writer = FastaWriter::new(output_path);
    for (header, seq) in &final_contigs {
        let seq_str = String::from_utf8_lossy(seq);
        writer.write_record(header, &seq_str)?;
    }

    info!("Long read integration complete:");
    info!("  Contigs extended: {}", stats.contigs_extended);
    info!("  Contigs joined: {}", stats.contigs_joined / 2);
    info!("  Total extension: {} bp", stats.total_extension_bp);

    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"GCGC"), b"GCGC");
    }

    #[test]
    fn test_minimizer_index() {
        let mut index = MinimizerIndex::new(11, 5);
        index.add_contig(0, b"ACGTACGTACGTACGTACGTACGTACGTACGT");

        let mappings = index.map_read(b"ACGTACGTACGTACGTACGT");
        assert!(!mappings.is_empty());
        assert_eq!(mappings[0].contig_id, 0);
    }

    #[test]
    fn test_hash_kmer() {
        let h1 = hash_kmer(b"ACGT");
        let h2 = hash_kmer(b"ACGT");
        assert_eq!(h1, h2);

        let h3 = hash_kmer(b"TGCA");
        assert_ne!(h1, h3);
    }

    #[test]
    fn test_long_read_config_default() {
        let config = LongReadConfig::default();
        assert_eq!(config.min_anchor_matches, 5);
        assert_eq!(config.min_read_length, 1000);
    }
}
