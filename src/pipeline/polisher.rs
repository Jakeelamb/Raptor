//! Contig Polisher - Correct errors in assembled contigs using raw reads
//!
//! Algorithm:
//! 1. Index contigs using minimizers
//! 2. Map reads to contigs
//! 3. Build pileup at each position
//! 4. Call consensus based on base frequencies
//! 5. Output polished contigs

use crate::io::fasta::{open_fasta, FastaWriter};
use crate::io::fastq::{stream_fastq_records_checked, try_open_fastq};
use ahash::AHashMap;
use std::cmp::Ordering;
use std::io::{BufRead, Result};
use tracing::info;

const REVERSE_STRAND_FLAG: usize = 1usize << (usize::BITS - 1);

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
        if seq.len() < k + w - 1 {
            return;
        }

        let mut last_emitted_pos = usize::MAX;
        for window_start in 0..=(seq.len() - k - w + 1) {
            let mut min_hash = u64::MAX;
            let mut min_pos = 0;

            for i in 0..w {
                let pos = window_start + i;
                if pos + k <= seq.len() {
                    let hash = hash_kmer(&seq[pos..pos + k]);
                    if hash < min_hash {
                        min_hash = hash;
                        min_pos = pos;
                    }
                }
            }

            if min_pos != last_emitted_pos {
                emit(min_pos, min_hash);
                last_emitted_pos = min_pos;
            }
        }
    }

    fn map_read(&self, sequence: &[u8]) -> Option<(usize, usize, bool)> {
        let mut hits: AHashMap<(usize, usize), usize> = AHashMap::new();

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
        let rc = reverse_complement(sequence);
        Self::for_each_minimizer(self.k, self.w, &rc, |read_pos, minimizer| {
            if let Some(entries) = self.index.get(&minimizer) {
                for &(contig_id, contig_pos) in entries {
                    let start = contig_pos.saturating_sub(read_pos);
                    *hits
                        .entry((contig_id, start | REVERSE_STRAND_FLAG))
                        .or_default() += 1;
                }
            }
        });

        // Find best hit with at least 3 minimizer matches using explicit,
        // insertion-order-independent tie-breaking.
        select_best_mapping_hit(hits, 3)
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
        let reader1 = try_open_fastq(reads1_path)?;
        for record in stream_fastq_records_checked(reader1) {
            let record = record?;
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
            let reader2 = try_open_fastq(reads2)?;
            for record in stream_fastq_records_checked(reader2) {
                let record = record?;
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
    let mut writer = FastaWriter::try_new(output_path)?;
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
        );
        assert!(result.is_err());
    }
}
