//! Contig Polisher - Correct errors in assembled contigs using raw reads
//!
//! Algorithm:
//! 1. Index contigs using minimizers
//! 2. Map reads to contigs
//! 3. Build pileup at each position
//! 4. Call consensus based on base frequencies
//! 5. Output polished contigs

use crate::io::fasta::{open_fasta, FastaWriter};
use crate::io::fastq::{for_each_fastq_sequence_checked, try_open_fastq};
use ahash::AHashMap;
use std::cmp::Ordering;
use std::io::{BufRead, Result};
use tracing::info;

const REVERSE_STRAND_FLAG: usize = 1usize << (usize::BITS - 1);
const FAST_MINIMIZER_WINDOW_LIMIT: usize = 255;
pub const DEFAULT_READ_MAPPING_K: usize = 15;
pub const DEFAULT_READ_MAPPING_W: usize = 10;
pub const DEFAULT_MIN_PRIMARY_MATCHES: usize = 3;
pub const DEFAULT_MIN_SCAFFOLD_MATCHES: usize = 2;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadMappingConfig {
    pub minimizer_k: usize,
    pub minimizer_w: usize,
    pub min_primary_matches: usize,
    pub min_scaffold_matches: usize,
}

impl Default for ReadMappingConfig {
    fn default() -> Self {
        Self {
            minimizer_k: DEFAULT_READ_MAPPING_K,
            minimizer_w: DEFAULT_READ_MAPPING_W,
            min_primary_matches: DEFAULT_MIN_PRIMARY_MATCHES,
            min_scaffold_matches: DEFAULT_MIN_SCAFFOLD_MATCHES,
        }
    }
}

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
pub(crate) struct PileupColumn {
    bases: [u32; 5], // A, C, G, T, gap
    depth: u32,
}

impl PileupColumn {
    #[inline]
    fn base_to_index(base: u8) -> usize {
        match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 4, // gap or ambiguous base
        }
    }

    #[inline]
    fn index_to_base(idx: usize) -> Option<u8> {
        match idx {
            0 => Some(b'A'),
            1 => Some(b'C'),
            2 => Some(b'G'),
            3 => Some(b'T'),
            _ => None, // Gap - deletion or unknown
        }
    }

    fn add_base(&mut self, base: u8) {
        let idx = Self::base_to_index(base);
        self.bases[idx] += 1;
        self.depth += 1;
    }

    fn consensus(&self, min_freq: f64, preferred_base: Option<u8>) -> Option<u8> {
        if self.depth == 0 {
            return None;
        }

        let mut max_count = 0u32;
        let mut tied = [false; 5];
        for (idx, &count) in self.bases.iter().enumerate() {
            if count > max_count {
                max_count = count;
                tied = [false; 5];
                tied[idx] = true;
            } else if count == max_count {
                tied[idx] = true;
            }
        }

        if max_count == 0 {
            return None;
        }

        let max_idx = preferred_base
            .map(Self::base_to_index)
            .filter(|&idx| tied[idx])
            .unwrap_or_else(|| tied.iter().position(|is_tied| *is_tied).unwrap_or(4));

        let freq = self.bases[max_idx] as f64 / self.depth as f64;
        if freq >= min_freq {
            Self::index_to_base(max_idx)
        } else {
            None
        }
    }
}

/// Minimizer index for read mapping
pub(crate) struct MinimizerIndex {
    index: AHashMap<u64, Vec<(usize, usize)>>, // minimizer -> [(contig_id, position)]
    k: usize,
    w: usize,
    contig_count: usize,
}

impl MinimizerIndex {
    pub(crate) fn new(k: usize, w: usize) -> Self {
        Self {
            index: AHashMap::new(),
            k,
            w,
            contig_count: 0,
        }
    }

    pub(crate) fn add_contig(&mut self, contig_id: usize, sequence: &[u8]) {
        self.contig_count = self.contig_count.max(contig_id.saturating_add(1));
        Self::for_each_minimizer(self.k, self.w, sequence, |pos, minimizer| {
            self.index
                .entry(minimizer)
                .or_default()
                .push((contig_id, pos));
        });
    }

    fn for_each_minimizer<F>(k: usize, w: usize, seq: &[u8], mut emit: F)
    where
        F: FnMut(usize, u64),
    {
        if k == 0 || w == 0 {
            return;
        }
        let Some(min_seq_len) = k.checked_add(w).and_then(|sum| sum.checked_sub(1)) else {
            return;
        };
        if seq.len() < min_seq_len {
            return;
        }
        if k > 32 || w > FAST_MINIMIZER_WINDOW_LIMIT {
            Self::for_each_minimizer_reference(k, w, seq, emit);
            return;
        }

        let rolling_mask = if k == 32 {
            u64::MAX
        } else {
            (1u64 << (k * 2)) - 1
        };
        let mut rolling = 0u64;
        let mut valid_run = 0usize;
        let mut last_emitted_pos = usize::MAX;
        let mut deque_positions = [0usize; FAST_MINIMIZER_WINDOW_LIMIT];
        let mut deque_hashes = [0u64; FAST_MINIMIZER_WINDOW_LIMIT];
        let mut deque_head = 0usize;
        let mut deque_len = 0usize;

        for (idx, &base) in seq.iter().enumerate() {
            let current_entry = if let Some(base_bits) = encode_dna_base_2bit(base) {
                rolling = ((rolling << 2) | base_bits) & rolling_mask;
                valid_run += 1;
                if valid_run >= k {
                    Some((idx + 1 - k, rolling))
                } else {
                    None
                }
            } else {
                rolling = 0;
                valid_run = 0;
                None
            };

            if idx + 1 >= min_seq_len {
                let window_start = idx + 1 - min_seq_len;
                while deque_len > 0 && deque_positions[deque_head] < window_start {
                    deque_head = (deque_head + 1) % FAST_MINIMIZER_WINDOW_LIMIT;
                    deque_len -= 1;
                }
            }

            if let Some((kmer_pos, hash)) = current_entry {
                while deque_len > 0 {
                    let back_idx = (deque_head + deque_len - 1) % FAST_MINIMIZER_WINDOW_LIMIT;
                    if deque_hashes[back_idx] > hash {
                        deque_len -= 1;
                    } else {
                        break;
                    }
                }

                let insert_idx = (deque_head + deque_len) % FAST_MINIMIZER_WINDOW_LIMIT;
                deque_positions[insert_idx] = kmer_pos;
                deque_hashes[insert_idx] = hash;
                deque_len += 1;
            }

            if idx + 1 < min_seq_len || deque_len == 0 {
                continue;
            }

            let min_pos = deque_positions[deque_head];
            if min_pos != last_emitted_pos {
                emit(min_pos, deque_hashes[deque_head]);
                last_emitted_pos = min_pos;
            }
        }
    }

    fn for_each_minimizer_reference<F>(k: usize, w: usize, seq: &[u8], mut emit: F)
    where
        F: FnMut(usize, u64),
    {
        let Some(min_seq_len) = k.checked_add(w).and_then(|sum| sum.checked_sub(1)) else {
            return;
        };
        if k == 0 || w == 0 || seq.len() < min_seq_len {
            return;
        }

        let mut last_emitted_pos = usize::MAX;
        for window_start in 0..=(seq.len() - min_seq_len) {
            let mut min_hash = 0u64;
            let mut min_pos = 0usize;
            let mut found_valid_kmer = false;

            for offset in 0..w {
                let pos = window_start + offset;
                let end = pos + k;
                if end > seq.len() {
                    break;
                }
                if let Some(hash) = hash_kmer(&seq[pos..end]) {
                    if !found_valid_kmer || hash < min_hash {
                        min_hash = hash;
                        min_pos = pos;
                        found_valid_kmer = true;
                    }
                }
            }

            if found_valid_kmer && min_pos != last_emitted_pos {
                emit(min_pos, min_hash);
                last_emitted_pos = min_pos;
            }
        }
    }

    fn map_read(&self, sequence: &[u8]) -> Option<(usize, usize, bool)> {
        let mut hits = AHashMap::new();
        let mut contig_hits = ContigHitScratch::default();
        let mut reverse_scratch = Vec::new();
        self.map_read_with_scaffold_support(
            sequence,
            &mut hits,
            &mut contig_hits,
            &mut reverse_scratch,
        )
        .0
    }

    pub(crate) fn map_read_with_scratch(
        &self,
        sequence: &[u8],
        hits: &mut AHashMap<(usize, usize), usize>,
        reverse_scratch: &mut Vec<u8>,
    ) -> Option<(usize, usize, bool)> {
        let mut contig_hits = ContigHitScratch::default();
        self.map_read_with_scaffold_support(sequence, hits, &mut contig_hits, reverse_scratch)
            .0
    }

    pub(crate) fn map_read_with_scaffold_support(
        &self,
        sequence: &[u8],
        hits: &mut AHashMap<(usize, usize), usize>,
        contig_hits: &mut ContigHitScratch,
        reverse_scratch: &mut Vec<u8>,
    ) -> (Option<(usize, usize, bool)>, Option<(usize, usize, bool)>) {
        self.map_read_with_scaffold_support_thresholds(
            sequence,
            hits,
            contig_hits,
            reverse_scratch,
            DEFAULT_MIN_PRIMARY_MATCHES,
            DEFAULT_MIN_SCAFFOLD_MATCHES,
        )
    }

    pub(crate) fn map_read_with_scaffold_support_thresholds(
        &self,
        sequence: &[u8],
        hits: &mut AHashMap<(usize, usize), usize>,
        contig_hits: &mut ContigHitScratch,
        reverse_scratch: &mut Vec<u8>,
        min_primary_matches: usize,
        min_scaffold_matches: usize,
    ) -> (Option<(usize, usize, bool)>, Option<(usize, usize, bool)>) {
        hits.clear();
        contig_hits.prepare(self.contig_count);

        // Forward mapping
        Self::for_each_minimizer(self.k, self.w, sequence, |read_pos, minimizer| {
            if let Some(entries) = self.index.get(&minimizer) {
                for &(contig_id, contig_pos) in entries {
                    let start = contig_pos.saturating_sub(read_pos);
                    *hits.entry((contig_id, start)).or_default() += 1;
                }
            }
        });

        // Reverse complement mapping
        reverse_complement_into(sequence, reverse_scratch);
        Self::for_each_minimizer(self.k, self.w, reverse_scratch, |read_pos, minimizer| {
            if let Some(entries) = self.index.get(&minimizer) {
                for &(contig_id, contig_pos) in entries {
                    let start = contig_pos.saturating_sub(read_pos);
                    *hits
                        .entry((contig_id, start | REVERSE_STRAND_FLAG))
                        .or_default() += 1;
                }
            }
        });

        let mut best_primary: Option<((usize, usize), usize)> = None;
        for (&(contig_id, pos), &count) in hits.iter() {
            let is_reverse = (pos & REVERSE_STRAND_FLAG) != 0;
            contig_hits.add(contig_id, is_reverse, count);

            if count < min_primary_matches {
                continue;
            }

            let candidate = ((contig_id, pos), count);
            match best_primary {
                None => best_primary = Some(candidate),
                Some(current) => {
                    if compare_mapping_hits(&candidate, &current).is_gt() {
                        best_primary = Some(candidate);
                    }
                }
            }
        }

        (
            best_primary.map(|((contig_id, pos), _)| {
                let (actual_pos, is_reverse) = decode_mapping_position(pos);
                (contig_id, actual_pos, is_reverse)
            }),
            contig_hits.best_hit(min_scaffold_matches),
        )
    }
}

#[derive(Debug, Default)]
pub(crate) struct ContigHitScratch {
    counts: Vec<u32>,
    touched: Vec<usize>,
}

impl ContigHitScratch {
    fn prepare(&mut self, contig_count: usize) {
        for &idx in &self.touched {
            self.counts[idx] = 0;
        }
        self.touched.clear();

        let needed = contig_count.saturating_mul(2);
        if self.counts.len() < needed {
            self.counts.resize(needed, 0);
        }
    }

    #[inline]
    fn add(&mut self, contig_id: usize, is_reverse: bool, count: usize) {
        let idx = contig_id
            .saturating_mul(2)
            .saturating_add(usize::from(is_reverse));
        if self.counts[idx] == 0 {
            self.touched.push(idx);
        }
        self.counts[idx] = self.counts[idx].saturating_add(count as u32);
    }

    fn best_hit(&self, min_matches: usize) -> Option<(usize, usize, bool)> {
        self.touched
            .iter()
            .copied()
            .filter_map(|idx| {
                let count = self.counts[idx] as usize;
                (count >= min_matches).then(|| {
                    let contig_id = idx / 2;
                    let is_reverse = (idx & 1) != 0;
                    ((contig_id, is_reverse), count)
                })
            })
            .max_by(compare_contig_hits)
            .map(|((contig_id, is_reverse), count)| (contig_id, count, is_reverse))
    }
}

#[inline]
fn decode_mapping_position(pos: usize) -> (usize, bool) {
    let is_reverse = (pos & REVERSE_STRAND_FLAG) != 0;
    let actual_pos = pos & !REVERSE_STRAND_FLAG;
    (actual_pos, is_reverse)
}

#[inline]
fn compare_mapping_hits(
    left: &((usize, usize), usize),
    right: &((usize, usize), usize),
) -> Ordering {
    let ((left_contig, left_pos), left_count) = left;
    let ((right_contig, right_pos), right_count) = right;
    let (left_actual_pos, left_is_reverse) = decode_mapping_position(*left_pos);
    let (right_actual_pos, right_is_reverse) = decode_mapping_position(*right_pos);

    left_count
        .cmp(right_count)
        // Prefer lower contig IDs when support ties.
        .then_with(|| right_contig.cmp(left_contig))
        // Prefer earlier start positions when contig IDs tie.
        .then_with(|| right_actual_pos.cmp(&left_actual_pos))
        // Prefer forward mappings over reverse-complement mappings.
        .then_with(|| right_is_reverse.cmp(&left_is_reverse))
}

#[inline]
fn compare_contig_hits(left: &((usize, bool), usize), right: &((usize, bool), usize)) -> Ordering {
    let ((left_contig, left_is_reverse), left_count) = left;
    let ((right_contig, right_is_reverse), right_count) = right;

    left_count
        .cmp(right_count)
        .then_with(|| right_contig.cmp(left_contig))
        .then_with(|| right_is_reverse.cmp(left_is_reverse))
}

#[inline]
fn select_best_mapping_hit<I>(hits: I, min_matches: usize) -> Option<(usize, usize, bool)>
where
    I: IntoIterator<Item = ((usize, usize), usize)>,
{
    hits.into_iter()
        .filter(|(_, count)| *count >= min_matches)
        .max_by(compare_mapping_hits)
        .map(|((contig_id, pos), _)| {
            let (actual_pos, is_reverse) = decode_mapping_position(pos);
            (contig_id, actual_pos, is_reverse)
        })
}

#[inline]
fn select_best_contig_hit<I>(hits: I, min_matches: usize) -> Option<(usize, usize, bool)>
where
    I: IntoIterator<Item = ((usize, bool), usize)>,
{
    hits.into_iter()
        .filter(|(_, count)| *count >= min_matches)
        .max_by(compare_contig_hits)
        .map(|((contig_id, is_reverse), count)| (contig_id, count, is_reverse))
}

#[inline]
fn populate_contig_hits_from_positional_hits(
    hits: &AHashMap<(usize, usize), usize>,
    contig_hits: &mut AHashMap<(usize, bool), usize>,
) {
    contig_hits.clear();
    for (&(contig_id, pos), &count) in hits {
        let is_reverse = (pos & REVERSE_STRAND_FLAG) != 0;
        *contig_hits.entry((contig_id, is_reverse)).or_default() += count;
    }
}

fn hash_kmer(kmer: &[u8]) -> Option<u64> {
    let mut hash = 0u64;
    for &base in kmer {
        let val = match base {
            b'A' | b'a' => 0u64,
            b'C' | b'c' => 1u64,
            b'G' | b'g' => 2u64,
            b'T' | b't' => 3u64,
            _ => return None,
        };
        hash = hash.wrapping_mul(4).wrapping_add(val);
    }
    Some(hash)
}

#[inline]
fn encode_dna_base_2bit(base: u8) -> Option<u64> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn reverse_complement_into(seq: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.reserve(seq.len());
    out.extend(seq.iter().rev().map(|&b| match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        _ => b'N',
    }));
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut rc = Vec::with_capacity(seq.len());
    reverse_complement_into(seq, &mut rc);
    rc
}

/// Read contigs from FASTA file
pub(crate) fn read_contigs(path: &str) -> Result<Vec<(String, Vec<u8>)>> {
    let reader = open_fasta(path);
    let mut contigs = Vec::new();
    let mut current_header = String::new();
    let mut current_seq = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_header.is_empty() {
                contigs.push((
                    std::mem::take(&mut current_header),
                    std::mem::take(&mut current_seq),
                ));
            }
            current_header = line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
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
    mapping_config: ReadMappingConfig,
) -> Result<PolishStats> {
    let mut stats = PolishStats::default();

    // Read initial contigs
    let mut contigs = read_contigs(contigs_path)?;
    stats.contigs_input = contigs.len();
    info!("Loaded {} contigs for polishing", contigs.len());

    for iter in 0..iterations {
        info!("Polishing iteration {}/{}", iter + 1, iterations);

        // Build minimizer index
        let mut index = MinimizerIndex::new(mapping_config.minimizer_k, mapping_config.minimizer_w);
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
        map_reads_into_pileups(
            reads1_path,
            &index,
            &mut pileups,
            &mut stats,
            mapping_config,
        )?;

        // Process reads from second file if provided
        if let Some(reads2) = reads2_path {
            info!("  Mapping reads from {}...", reads2);
            map_reads_into_pileups(reads2, &index, &mut pileups, &mut stats, mapping_config)?;
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

    info!("Writing polished contigs to {}", output_path);
    write_contigs(&contigs, output_path)?;

    info!(
        "Polishing complete: {} corrections ({} substitutions, {} insertions, {} deletions)",
        stats.corrections, stats.substitutions, stats.insertions, stats.deletions
    );

    Ok(stats)
}

/// Add a read alignment to the pileup
pub(crate) fn add_to_pileup(pileup: &mut [PileupColumn], start: usize, read: &[u8]) {
    for (i, &base) in read.iter().enumerate() {
        let pos = start + i;
        if pos < pileup.len() {
            pileup[pos].add_base(base);
        }
    }
}

fn map_reads_into_pileups(
    reads_path: &str,
    index: &MinimizerIndex,
    pileups: &mut [Vec<PileupColumn>],
    stats: &mut PolishStats,
    mapping_config: ReadMappingConfig,
) -> Result<()> {
    let reader = try_open_fastq(reads_path)?;
    let mut hits = AHashMap::new();
    let mut contig_hits = ContigHitScratch::default();
    let mut reverse_scratch = Vec::new();

    for_each_fastq_sequence_checked(reader, |sequence| {
        stats.reads_processed += 1;

        if let Some((contig_id, start, is_reverse)) = index
            .map_read_with_scaffold_support_thresholds(
                sequence,
                &mut hits,
                &mut contig_hits,
                &mut reverse_scratch,
                mapping_config.min_primary_matches,
                mapping_config.min_scaffold_matches,
            )
            .0
        {
            stats.reads_mapped += 1;
            let read_seq = if is_reverse {
                reverse_scratch.as_slice()
            } else {
                sequence
            };
            add_to_pileup(&mut pileups[contig_id], start, read_seq);
        }
        Ok(())
    })?;

    Ok(())
}

/// Apply corrections based on pileup consensus
pub(crate) fn apply_corrections(
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

            if let Some(consensus_base) = column.consensus(MIN_FREQ, seq.get(pos).copied()) {
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
                if old_base == b'-' {
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

pub(crate) fn write_contigs(contigs: &[(String, Vec<u8>)], output_path: &str) -> Result<()> {
    let mut writer = FastaWriter::try_new(output_path)?;
    for (header, seq) in contigs {
        let seq_str = String::from_utf8_lossy(seq);
        writer.write_record(header, &seq_str)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::bool::ANY;
    use proptest::prelude::*;
    use std::io::Write;
    use tempfile::{NamedTempFile, TempDir};

    #[test]
    fn test_pileup_column() {
        let mut col = PileupColumn::default();
        col.add_base(b'A');
        col.add_base(b'A');
        col.add_base(b'A');
        col.add_base(b'C');

        assert_eq!(col.depth, 4);
        assert_eq!(col.consensus(0.5, None), Some(b'A'));
        assert_eq!(col.consensus(0.9, None), None);
    }

    #[test]
    fn test_pileup_column_tie_prefers_requested_reference_base() {
        let mut col = PileupColumn::default();
        col.add_base(b'A');
        col.add_base(b'A');
        col.add_base(b'C');
        col.add_base(b'C');

        assert_eq!(col.consensus(0.5, None), Some(b'A'));
        assert_eq!(col.consensus(0.5, Some(b'C')), Some(b'C'));
        assert_eq!(col.consensus(0.75, Some(b'C')), None);
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

    fn collect_minimizers<F>(collector: F, k: usize, w: usize, sequence: &[u8]) -> Vec<(usize, u64)>
    where
        F: Fn(usize, usize, &[u8], &mut dyn FnMut(usize, u64)),
    {
        let mut minimizers = Vec::new();
        collector(k, w, sequence, &mut |pos, hash| {
            minimizers.push((pos, hash))
        });
        minimizers
    }

    #[test]
    fn fast_minimizer_path_matches_reference_with_ambiguous_bases_and_ties() {
        let sequence = b"ACGTACGTNNACGTAAAACGT";
        let observed = collect_minimizers(
            |k, w, seq, emit| MinimizerIndex::for_each_minimizer(k, w, seq, emit),
            5,
            4,
            sequence,
        );
        let expected = collect_minimizers(
            |k, w, seq, emit| MinimizerIndex::for_each_minimizer_reference(k, w, seq, emit),
            5,
            4,
            sequence,
        );

        assert_eq!(observed, expected);
    }

    #[test]
    fn hash_kmer_rejects_ambiguous_bases() {
        assert_eq!(hash_kmer(b"ACGT"), Some(27));
        assert_eq!(hash_kmer(b"ACGN"), None);
    }

    #[test]
    fn minimizer_index_skips_windows_with_only_ambiguous_kmers() {
        let mut index = MinimizerIndex::new(3, 2);
        index.add_contig(0, b"AAAAAAAAAAAA");

        // Ambiguous read should not map via collapsed hash collisions.
        assert_eq!(index.map_read(b"NNNNNNNNNNNN"), None);
    }

    #[test]
    fn select_best_mapping_hit_prefers_contig_position_and_forward_strand_on_ties() {
        let hits = vec![
            ((1, 10), 7),
            ((0, 20), 7),
            ((0, 5 | REVERSE_STRAND_FLAG), 7),
            ((0, 5), 7),
        ];

        let best = select_best_mapping_hit(hits.iter().copied(), 3);
        assert_eq!(best, Some((0, 5, false)));

        let mut reversed = hits.clone();
        reversed.reverse();
        let best_reversed = select_best_mapping_hit(reversed.iter().copied(), 3);
        assert_eq!(best_reversed, Some((0, 5, false)));
    }

    #[test]
    fn populate_contig_hits_from_positional_hits_aggregates_all_positions_per_strand() {
        let mut hits = AHashMap::new();
        hits.insert((0, 10), 2);
        hits.insert((0, 20), 1);
        hits.insert((0, 30 | REVERSE_STRAND_FLAG), 3);
        hits.insert((1, 15), 3);
        hits.insert((1, 25 | REVERSE_STRAND_FLAG), 2);

        let mut contig_hits = AHashMap::new();
        populate_contig_hits_from_positional_hits(&hits, &mut contig_hits);

        assert_eq!(contig_hits.get(&(0, false)), Some(&3));
        assert_eq!(contig_hits.get(&(0, true)), Some(&3));
        assert_eq!(contig_hits.get(&(1, false)), Some(&3));
        assert_eq!(contig_hits.get(&(1, true)), Some(&2));
        assert_eq!(
            select_best_contig_hit(contig_hits.iter().map(|(k, v)| (*k, *v)), 3),
            Some((0, 3, false))
        );
    }

    #[test]
    fn contig_hit_scratch_aggregates_and_reuses_storage() {
        let mut scratch = ContigHitScratch::default();
        scratch.prepare(2);
        scratch.add(0, false, 2);
        scratch.add(0, false, 1);
        scratch.add(0, true, 3);
        scratch.add(1, false, 3);
        scratch.add(1, true, 2);
        assert_eq!(scratch.best_hit(3), Some((0, 3, false)));

        scratch.prepare(2);
        assert_eq!(scratch.best_hit(1), None);
        scratch.add(1, true, 4);
        assert_eq!(scratch.best_hit(3), Some((1, 4, true)));
    }

    #[test]
    fn apply_corrections_counts_n_to_acgt_as_substitution() {
        let mut contigs = vec![("contig_1".to_string(), b"NAAAA".to_vec())];
        let mut pileups = vec![vec![PileupColumn::default(); 5]];
        for _ in 0..6 {
            pileups[0][0].add_base(b'A');
        }

        let mut stats = PolishStats::default();
        let corrections = apply_corrections(&mut contigs, &pileups, &mut stats);

        assert_eq!(corrections, 1);
        assert_eq!(contigs[0].1[0], b'A');
        assert_eq!(stats.substitutions, 1);
        assert_eq!(stats.insertions, 0);
        assert_eq!(stats.deletions, 0);
    }

    fn canonical_best_mapping_hit(
        hits: &[((usize, usize), usize)],
        min_matches: usize,
    ) -> Option<(usize, usize, bool)> {
        let mut filtered: Vec<((usize, usize), usize)> = hits
            .iter()
            .copied()
            .filter(|(_, count)| *count >= min_matches)
            .collect();
        filtered.sort_unstable_by(|a, b| compare_mapping_hits(b, a));
        filtered.first().map(|((contig_id, pos), _)| {
            let (actual_pos, is_reverse) = decode_mapping_position(*pos);
            (*contig_id, actual_pos, is_reverse)
        })
    }

    proptest! {
        #[test]
        fn fast_minimizer_path_matches_reference_under_random_input(
            sequence in prop::collection::vec(
                prop_oneof![
                    Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T'),
                    Just(b'N'), Just(b'R'), Just(b'a'), Just(b't')
                ],
                0..96
            ),
            k in 1usize..20,
            w in 1usize..16
        ) {
            let observed = collect_minimizers(
                |k, w, seq, emit| MinimizerIndex::for_each_minimizer(k, w, seq, emit),
                k,
                w,
                &sequence,
            );
            let expected = collect_minimizers(
                |k, w, seq, emit| MinimizerIndex::for_each_minimizer_reference(k, w, seq, emit),
                k,
                w,
                &sequence,
            );
            prop_assert_eq!(observed, expected);
        }

        #[test]
        fn select_best_mapping_hit_is_invariant_to_candidate_order(
            raw_hits in prop::collection::vec(((0usize..8, 0usize..256usize, ANY), 0usize..10usize), 1..128)
        ) {
            let hits: Vec<((usize, usize), usize)> = raw_hits
                .iter()
                .map(|((contig, pos, is_reverse), count)| {
                    let encoded_pos = if *is_reverse { *pos | REVERSE_STRAND_FLAG } else { *pos };
                    ((*contig, encoded_pos), *count)
                })
                .collect();

            let expected = canonical_best_mapping_hit(&hits, 3);
            let observed = select_best_mapping_hit(hits.iter().copied(), 3);
            prop_assert_eq!(observed, expected);

            let mut reversed = hits.clone();
            reversed.reverse();
            let observed_reversed = select_best_mapping_hit(reversed.iter().copied(), 3);
            prop_assert_eq!(observed_reversed, expected);
        }
    }

    #[test]
    fn polish_contigs_returns_error_for_invalid_output_path() {
        let mut contigs = NamedTempFile::new().unwrap();
        writeln!(contigs, ">contig_1").unwrap();
        writeln!(contigs, "ACGTACGTACGT").unwrap();

        let mut reads = NamedTempFile::new().unwrap();
        writeln!(reads, "@read_1").unwrap();
        writeln!(reads, "ACGTACGT").unwrap();
        writeln!(reads, "+").unwrap();
        writeln!(reads, "IIIIIIII").unwrap();

        let output_dir = TempDir::new().unwrap();
        let result = polish_contigs(
            contigs.path().to_str().unwrap(),
            reads.path().to_str().unwrap(),
            None,
            output_dir.path().to_str().unwrap(),
            1,
            ReadMappingConfig::default(),
        );
        assert!(result.is_err());
    }
}
