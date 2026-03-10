//! Long Read Integration - Hybrid assembly with PacBio/Nanopore reads
//!
//! Strategy:
//! 1. Use short reads to build accurate de Bruijn graph
//! 2. Map long reads to contigs using minimizers
//! 3. Use long reads to span repeat regions
//! 4. Anchor long reads to extend/connect contigs

use crate::io::fasta::{open_fasta, FastaWriter};
use crate::io::fastq::{stream_fastq_records_checked, try_open_fastq};
use ahash::{AHashMap, AHashSet};
use std::cmp::Ordering;
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

#[derive(Debug, Clone, Copy)]
struct MappingHitStats {
    count: usize,
    read_min: usize,
    read_max: usize,
    contig_min: usize,
    contig_max: usize,
}

impl MappingHitStats {
    #[inline]
    fn new(read_pos: usize, contig_pos: usize) -> Self {
        Self {
            count: 1,
            read_min: read_pos,
            read_max: read_pos,
            contig_min: contig_pos,
            contig_max: contig_pos,
        }
    }

    #[inline]
    fn observe(&mut self, read_pos: usize, contig_pos: usize) {
        self.count += 1;
        self.read_min = self.read_min.min(read_pos);
        self.read_max = self.read_max.max(read_pos);
        self.contig_min = self.contig_min.min(contig_pos);
        self.contig_max = self.contig_max.max(contig_pos);
    }
}

#[inline]
fn compare_read_mappings(a: &ReadMapping, b: &ReadMapping) -> Ordering {
    b.score
        .cmp(&a.score)
        .then_with(|| a.contig_id.cmp(&b.contig_id))
        .then_with(|| a.is_reverse.cmp(&b.is_reverse))
        .then_with(|| a.contig_start.cmp(&b.contig_start))
        .then_with(|| a.contig_end.cmp(&b.contig_end))
        .then_with(|| a.read_start.cmp(&b.read_start))
        .then_with(|| a.read_end.cmp(&b.read_end))
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
        let mut hits: AHashMap<(usize, bool), MappingHitStats> = AHashMap::new();

        // Collect minimizer hits
        for (read_pos, minimizer) in self.extract_minimizers(sequence) {
            if let Some(entries) = self.index.get(&minimizer) {
                for &(contig_id, contig_pos, is_rc) in entries {
                    hits.entry((contig_id, is_rc))
                        .and_modify(|stats| stats.observe(read_pos, contig_pos))
                        .or_insert_with(|| MappingHitStats::new(read_pos, contig_pos));
                }
            }
        }

        // Convert hits to mappings
        let mut mappings = Vec::new();
        for ((contig_id, is_reverse), stats) in hits {
            if stats.count < 3 {
                continue;
            }

            mappings.push(ReadMapping {
                contig_id,
                contig_start: stats.contig_min,
                contig_end: stats.contig_max + self.k,
                read_start: stats.read_min,
                read_end: stats.read_max + self.k,
                is_reverse,
                score: stats.count,
            });
        }

        mappings.sort_unstable_by(compare_read_mappings);
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

#[inline]
fn compare_extension_sequences(left: &[u8], right: &[u8]) -> Ordering {
    left.len().cmp(&right.len()).then_with(|| right.cmp(left))
}

fn choose_consensus_extension<'a>(extensions: &[&'a Vec<u8>]) -> Option<&'a Vec<u8>> {
    if extensions.is_empty() {
        return None;
    }

    let mut sorted = extensions.to_vec();
    sorted.sort_unstable_by(|left, right| left.as_slice().cmp(right.as_slice()));

    let mut best_seq = sorted[0];
    let mut best_count = 1usize;

    let mut run_seq = sorted[0];
    let mut run_count = 1usize;
    for seq in sorted.iter().copied().skip(1) {
        if seq.as_slice() == run_seq.as_slice() {
            run_count += 1;
            continue;
        }

        if run_count > best_count
            || (run_count == best_count
                && compare_extension_sequences(run_seq.as_slice(), best_seq.as_slice())
                    == Ordering::Greater)
        {
            best_seq = run_seq;
            best_count = run_count;
        }
        run_seq = seq;
        run_count = 1;
    }

    if run_count > best_count
        || (run_count == best_count
            && compare_extension_sequences(run_seq.as_slice(), best_seq.as_slice())
                == Ordering::Greater)
    {
        best_seq = run_seq;
    }

    Some(best_seq)
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

#[inline]
fn compare_contig_links(a: &ContigLink, b: &ContigLink) -> Ordering {
    b.supporting_reads
        .cmp(&a.supporting_reads)
        .then_with(|| a.contig_a.cmp(&b.contig_a))
        .then_with(|| a.contig_b.cmp(&b.contig_b))
        .then_with(|| a.orientation_a.cmp(&b.orientation_a))
        .then_with(|| a.orientation_b.cmp(&b.orientation_b))
        .then_with(|| a.gap_estimate.cmp(&b.gap_estimate))
        .then_with(|| b.bridging_sequences.len().cmp(&a.bridging_sequences.len()))
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

    let reader = try_open_fastq(long_reads_path)?;
    for record in stream_fastq_records_checked(reader) {
        let record = record?;
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
                if best.contig_end > contig_len.saturating_sub(50)
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
                // Consensus uses support first, with deterministic tie-breakers.
                if let Some(best_ext) = choose_consensus_extension(&relevant) {
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
    sorted_links.sort_unstable_by(|a, b| compare_contig_links(a, b));

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
    let mut writer = FastaWriter::try_new(output_path)?;
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
    use std::io::{self, Write};
    use tempfile::{NamedTempFile, TempDir};

    fn write_long_reads_fastq(path: &std::path::Path, reads: &[&str]) {
        let mut file = std::fs::File::create(path).expect("create FASTQ");
        for (idx, read) in reads.iter().enumerate() {
            writeln!(file, "@read_{idx}").expect("write read header");
            writeln!(file, "{read}").expect("write read sequence");
            writeln!(file, "+").expect("write plus line");
            writeln!(file, "{}", "I".repeat(read.len())).expect("write quality");
        }
    }

    fn read_single_output_sequence(path: &std::path::Path) -> String {
        let content = std::fs::read_to_string(path).expect("read output FASTA");
        let mut seq = String::new();
        for line in content.lines() {
            if !line.starts_with('>') {
                seq.push_str(line);
            }
        }
        seq
    }

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
    fn map_read_reports_full_span_and_score_for_exact_match() {
        let mut index = MinimizerIndex::new(3, 2);
        let seq = b"ACGTACGTACGTAC";
        index.add_contig(0, seq);

        let mappings = index.map_read(seq);
        let best = mappings.first().expect("expected at least one mapping");

        assert_eq!(best.contig_id, 0);
        assert!(!best.is_reverse);
        assert_eq!(best.read_start, 0);
        assert_eq!(best.contig_start, 0);
        assert_eq!(best.read_end, best.contig_end);
        assert!(best.read_end >= seq.len().saturating_sub(1));
        assert!(best.read_end <= seq.len());
        assert!(best.score >= 3);
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

    #[test]
    fn map_read_tie_breaks_are_deterministic() {
        let mut mappings = vec![
            ReadMapping {
                contig_id: 2,
                contig_start: 10,
                contig_end: 20,
                read_start: 1,
                read_end: 11,
                is_reverse: false,
                score: 7,
            },
            ReadMapping {
                contig_id: 1,
                contig_start: 11,
                contig_end: 21,
                read_start: 2,
                read_end: 12,
                is_reverse: true,
                score: 7,
            },
            ReadMapping {
                contig_id: 1,
                contig_start: 8,
                contig_end: 18,
                read_start: 0,
                read_end: 10,
                is_reverse: false,
                score: 7,
            },
            ReadMapping {
                contig_id: 0,
                contig_start: 0,
                contig_end: 10,
                read_start: 0,
                read_end: 10,
                is_reverse: false,
                score: 9,
            },
        ];

        mappings.sort_unstable_by(compare_read_mappings);
        let ordering: Vec<(usize, bool, usize)> = mappings
            .iter()
            .map(|m| (m.contig_id, m.is_reverse, m.score))
            .collect();

        assert_eq!(
            ordering,
            vec![(0, false, 9), (1, false, 7), (1, true, 7), (2, false, 7)]
        );
    }

    #[test]
    fn contig_link_tie_breaks_are_deterministic() {
        let mut links = vec![
            ContigLink {
                contig_a: 5,
                contig_b: 8,
                orientation_a: false,
                orientation_b: false,
                gap_estimate: 0,
                supporting_reads: 4,
                bridging_sequences: vec![b"AAA".to_vec()],
            },
            ContigLink {
                contig_a: 2,
                contig_b: 3,
                orientation_a: false,
                orientation_b: false,
                gap_estimate: 0,
                supporting_reads: 4,
                bridging_sequences: vec![b"CCC".to_vec()],
            },
            ContigLink {
                contig_a: 1,
                contig_b: 2,
                orientation_a: false,
                orientation_b: false,
                gap_estimate: 0,
                supporting_reads: 5,
                bridging_sequences: vec![b"GGG".to_vec()],
            },
        ];

        links.sort_unstable_by(compare_contig_links);
        let ordering: Vec<(usize, usize, usize)> = links
            .iter()
            .map(|l| (l.contig_a, l.contig_b, l.supporting_reads))
            .collect();

        assert_eq!(ordering, vec![(1, 2, 5), (2, 3, 4), (5, 8, 4)]);
    }

    #[test]
    fn integrate_long_reads_handles_short_contigs_without_underflow() {
        let mut contigs = NamedTempFile::new().expect("temp contig fasta");
        writeln!(contigs, ">contig_1").expect("write header");
        writeln!(contigs, "ACGTACGTACGT").expect("write sequence");

        let mut reads = NamedTempFile::new().expect("temp long-read fastq");
        let read_seq = "ACGTACGTACGTAAAA";
        writeln!(reads, "@read_1").expect("write read id");
        writeln!(reads, "{}", read_seq).expect("write read seq");
        writeln!(reads, "+").expect("write plus line");
        writeln!(reads, "{}", "I".repeat(read_seq.len())).expect("write quality");

        let output = NamedTempFile::new().expect("temp output fasta");
        let config = LongReadConfig {
            min_anchor_matches: 3,
            min_overlap: 1,
            min_read_length: 8,
            minimizer_k: 3,
            minimizer_w: 2,
        };

        let stats = integrate_long_reads(
            contigs.path().to_str().expect("utf8 contigs path"),
            reads.path().to_str().expect("utf8 reads path"),
            output.path().to_str().expect("utf8 output path"),
            config,
        )
        .expect("integration should succeed");

        assert_eq!(stats.contigs_input, 1);
        assert_eq!(stats.long_reads_processed, 1);
        assert_eq!(stats.long_reads_mapped, 1);
        assert_eq!(stats.output_contigs, 1);
    }

    #[test]
    fn integrate_long_reads_returns_error_for_invalid_output_path() {
        let mut contigs = NamedTempFile::new().expect("temp contig fasta");
        writeln!(contigs, ">contig_1").expect("write header");
        writeln!(contigs, "ACGTACGTACGT").expect("write sequence");

        let mut reads = NamedTempFile::new().expect("temp long-read fastq");
        writeln!(reads, "@read_1").expect("write read id");
        writeln!(reads, "ACGTACGT").expect("write read seq");
        writeln!(reads, "+").expect("write plus line");
        writeln!(reads, "IIIIIIII").expect("write quality");

        let output_dir = TempDir::new().expect("temp output directory");
        let result = integrate_long_reads(
            contigs.path().to_str().expect("utf8 contigs path"),
            reads.path().to_str().expect("utf8 reads path"),
            output_dir.path().to_str().expect("utf8 output dir"),
            LongReadConfig::default(),
        );
        assert!(result.is_err());
    }

    #[test]
    fn integrate_long_reads_rejects_truncated_fastq() {
        let mut contigs = NamedTempFile::new().expect("temp contig fasta");
        writeln!(contigs, ">contig_1").expect("write header");
        writeln!(contigs, "ACGTACGTACGT").expect("write sequence");

        let mut reads = NamedTempFile::new().expect("temp long-read fastq");
        writeln!(reads, "@read_1").expect("write read id");
        writeln!(reads, "ACGTACGT").expect("write read seq");
        writeln!(reads, "+").expect("write plus line");
        reads.flush().expect("flush reads");

        let output = NamedTempFile::new().expect("temp output fasta");
        let err = integrate_long_reads(
            contigs.path().to_str().expect("utf8 contigs path"),
            reads.path().to_str().expect("utf8 reads path"),
            output.path().to_str().expect("utf8 output path"),
            LongReadConfig::default(),
        )
        .expect_err("truncated FASTQ must return an error");

        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn choose_consensus_extension_prefers_support_then_length_then_lexicographic_order() {
        let a = b"AAA".to_vec();
        let b = b"AAAA".to_vec();
        let c = b"TTTT".to_vec();
        let d = b"CCCC".to_vec();

        let best_supported = choose_consensus_extension(&[&a, &b, &a, &c, &d, &c, &c]);
        assert_eq!(best_supported.expect("consensus exists"), &c);

        let best_tie_by_length = choose_consensus_extension(&[&a, &b]);
        assert_eq!(best_tie_by_length.expect("consensus exists"), &b);

        let best_tie_by_lex = choose_consensus_extension(&[&c, &d]);
        assert_eq!(best_tie_by_lex.expect("consensus exists"), &d);
    }

    #[test]
    fn integrate_long_reads_consensus_extension_is_order_invariant_and_majority_driven() {
        let mut contigs = NamedTempFile::new().expect("temp contig fasta");
        writeln!(contigs, ">contig_1").expect("write header");
        writeln!(contigs, "ACGTACGTACGT").expect("write sequence");

        let temp_dir = TempDir::new().expect("temp directory");
        let reads_a = temp_dir.path().join("reads_a.fastq");
        let reads_b = temp_dir.path().join("reads_b.fastq");
        let out_a = temp_dir.path().join("out_a.fasta");
        let out_b = temp_dir.path().join("out_b.fasta");

        let majority_short = "ACGTACGTACGTAAA";
        let minority_long = "ACGTACGTACGTTTTTTT";
        write_long_reads_fastq(
            &reads_a,
            &[
                majority_short,
                minority_long,
                majority_short,
                majority_short,
                majority_short,
            ],
        );
        write_long_reads_fastq(
            &reads_b,
            &[
                minority_long,
                majority_short,
                majority_short,
                majority_short,
                majority_short,
            ],
        );

        let config = LongReadConfig {
            min_anchor_matches: 3,
            min_overlap: 1,
            min_read_length: 8,
            minimizer_k: 3,
            minimizer_w: 2,
        };

        integrate_long_reads(
            contigs.path().to_str().expect("utf8 contig path"),
            reads_a.to_str().expect("utf8 reads_a path"),
            out_a.to_str().expect("utf8 out_a path"),
            config.clone(),
        )
        .expect("first integration should succeed");

        integrate_long_reads(
            contigs.path().to_str().expect("utf8 contig path"),
            reads_b.to_str().expect("utf8 reads_b path"),
            out_b.to_str().expect("utf8 out_b path"),
            config,
        )
        .expect("second integration should succeed");

        let assembled_a = read_single_output_sequence(&out_a);
        let assembled_b = read_single_output_sequence(&out_b);
        assert_eq!(assembled_a, assembled_b);
        assert!(
            assembled_a.ends_with("AA"),
            "expected majority extension suffix to preserve A-supported tail, observed {assembled_a}"
        );
        assert!(
            !assembled_a.ends_with("TTTT"),
            "longest minority extension should not win, observed {assembled_a}"
        );
    }
}
