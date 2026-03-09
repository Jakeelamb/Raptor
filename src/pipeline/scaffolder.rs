//! Scaffolder - Link contigs into scaffolds using paired-end reads
//!
//! Algorithm:
//! 1. Index contigs using minimizers
//! 2. Map paired reads to contigs
//! 3. Build scaffold graph from read pair links
//! 4. Filter edges by minimum support
//! 5. Extract scaffold paths
//! 6. Output scaffolds with gap estimates

use crate::io::fasta::{open_fasta, FastaWriter};
use crate::io::fastq::{open_fastq, stream_paired_fastq_records};
use ahash::{AHashMap, AHashSet};
use std::io::{BufRead, Result};
use tracing::info;

/// Statistics from scaffolding
#[derive(Debug, Clone, Default)]
pub struct ScaffoldStats {
    pub contigs_input: usize,
    pub reads_processed: u64,
    pub reads_mapped: u64,
    pub links_found: usize,
    pub links_after_filter: usize,
    pub num_scaffolds: usize,
    pub scaffold_n50: usize,
    pub total_length: usize,
    pub total_gaps: usize,
}

/// Orientation of a contig in a scaffold
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Orientation {
    Forward,
    Reverse,
}

/// A link between two contigs supported by read pairs
#[derive(Debug, Clone)]
struct ContigLink {
    contig_a: usize,
    contig_b: usize,
    orientation_a: Orientation,
    orientation_b: Orientation,
    gap_estimates: Vec<i32>,
    support_count: usize,
}

impl ContigLink {
    fn median_gap(&self) -> i32 {
        if self.gap_estimates.is_empty() {
            return 0;
        }
        let mut sorted = self.gap_estimates.clone();
        sorted.sort_unstable();
        sorted[sorted.len() / 2]
    }
}

/// A scaffold is a list of oriented contigs with gaps
#[derive(Debug, Clone)]
pub struct Scaffold {
    pub contigs: Vec<(usize, Orientation)>,
    pub gaps: Vec<i32>,
}

/// Minimizer index for fast contig mapping
struct MinimizerIndex {
    /// Map from minimizer hash to (contig_id, position, is_reverse)
    index: AHashMap<u64, Vec<(usize, usize, bool)>>,
    k: usize,
    w: usize, // window size
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
        if sequence.len() < self.k {
            return;
        }

        // Extract minimizers from forward strand
        for (pos, minimizer) in self.extract_minimizers(sequence) {
            self.index
                .entry(minimizer)
                .or_default()
                .push((contig_id, pos, false));
        }

        // Extract minimizers from reverse complement
        let rc = reverse_complement(sequence);
        for (pos, minimizer) in self.extract_minimizers(&rc) {
            let orig_pos = sequence.len() - pos - self.k;
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

    fn query(&self, sequence: &[u8]) -> Vec<(usize, usize, bool)> {
        let mut hits: AHashMap<(usize, bool), usize> = AHashMap::new();

        for (_, minimizer) in self.extract_minimizers(sequence) {
            if let Some(entries) = self.index.get(&minimizer) {
                for &(contig_id, _pos, is_rev) in entries {
                    *hits.entry((contig_id, is_rev)).or_default() += 1;
                }
            }
        }

        // Return hits with sufficient support (at least 2 minimizer matches)
        hits.into_iter()
            .filter(|(_, count)| *count >= 2)
            .map(|((contig_id, is_rev), count)| (contig_id, count, is_rev))
            .collect()
    }
}

/// Hash a k-mer using simple polynomial rolling hash
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

/// Compute reverse complement
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

/// Read contigs from a FASTA file
fn read_contigs(path: &str) -> Result<Vec<(String, String)>> {
    let reader = open_fasta(path);
    let mut contigs = Vec::new();
    let mut current_header = String::new();
    let mut current_seq = String::new();

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
            current_seq.push_str(line.trim());
        }
    }

    if !current_header.is_empty() {
        contigs.push((current_header, current_seq));
    }

    Ok(contigs)
}

/// Main scaffolding function
pub fn scaffold_contigs(
    contigs_path: &str,
    reads1_path: &str,
    reads2_path: &str,
    output_path: &str,
    min_links: usize,
) -> Result<ScaffoldStats> {
    let mut stats = ScaffoldStats::default();

    info!("Reading contigs from {}", contigs_path);
    let contigs = read_contigs(contigs_path)?;
    stats.contigs_input = contigs.len();
    info!("Loaded {} contigs", contigs.len());

    if contigs.is_empty() {
        return Ok(stats);
    }

    // Build minimizer index
    info!("Building minimizer index...");
    let mut index = MinimizerIndex::new(15, 10);
    for (i, (_header, seq)) in contigs.iter().enumerate() {
        index.add_contig(i, seq.as_bytes());
    }

    // Map read pairs and collect links
    info!("Mapping paired-end reads...");
    let mut links: AHashMap<(usize, usize), ContigLink> = AHashMap::new();

    let reader1 = open_fastq(reads1_path);
    let reader2 = open_fastq(reads2_path);

    // Estimate insert size from first batch of reads
    let mut insert_sizes: Vec<i32> = Vec::new();
    const INSERT_SAMPLE_SIZE: usize = 10000;

    for (r1, r2) in stream_paired_fastq_records(reader1, reader2) {
        stats.reads_processed += 2;

        let hits1 = index.query(r1.sequence.as_bytes());
        let hits2 = index.query(r2.sequence.as_bytes());

        if hits1.is_empty() || hits2.is_empty() {
            continue;
        }

        stats.reads_mapped += 2;

        // Get best hit for each read
        let best1 = hits1.iter().max_by_key(|(_, count, _)| count);
        let best2 = hits2.iter().max_by_key(|(_, count, _)| count);

        if let (Some(&(contig1, _, rev1)), Some(&(contig2, _, rev2))) = (best1, best2) {
            // Skip if both reads map to the same contig (used for insert size estimation)
            if contig1 == contig2 && insert_sizes.len() < INSERT_SAMPLE_SIZE {
                // Rough insert size estimate based on read lengths
                insert_sizes.push((r1.sequence.len() + r2.sequence.len()) as i32);
                continue;
            }

            // Record link between different contigs
            if contig1 != contig2 {
                let key = if contig1 < contig2 {
                    (contig1, contig2)
                } else {
                    (contig2, contig1)
                };

                let link = links.entry(key).or_insert_with(|| ContigLink {
                    contig_a: key.0,
                    contig_b: key.1,
                    orientation_a: if rev1 {
                        Orientation::Reverse
                    } else {
                        Orientation::Forward
                    },
                    orientation_b: if rev2 {
                        Orientation::Reverse
                    } else {
                        Orientation::Forward
                    },
                    gap_estimates: Vec::new(),
                    support_count: 0,
                });

                // Estimate gap (negative means overlap)
                let estimated_insert = if !insert_sizes.is_empty() {
                    let mut sorted = insert_sizes.clone();
                    sorted.sort_unstable();
                    sorted[sorted.len() / 2]
                } else {
                    500 // Default insert size
                };

                let gap = estimated_insert - (r1.sequence.len() + r2.sequence.len()) as i32;
                link.gap_estimates.push(gap);
                link.support_count += 1;
            }
        }

        if stats.reads_processed % 1_000_000 == 0 {
            info!(
                "  Processed {} million read pairs...",
                stats.reads_processed / 2_000_000
            );
        }
    }

    stats.links_found = links.len();
    info!("Found {} potential links", links.len());

    // Filter links by minimum support
    let filtered_links: Vec<ContigLink> = links
        .into_values()
        .filter(|link| link.support_count >= min_links)
        .collect();
    stats.links_after_filter = filtered_links.len();
    info!(
        "Retained {} links with >= {} supporting pairs",
        filtered_links.len(),
        min_links
    );

    // Build scaffold graph and extract paths
    let scaffolds = build_scaffolds(&contigs, &filtered_links);
    stats.num_scaffolds = scaffolds.len();

    // Write output
    info!("Writing {} scaffolds to {}", scaffolds.len(), output_path);
    let mut writer = FastaWriter::new(output_path);
    let mut scaffold_lengths: Vec<usize> = Vec::new();

    for (i, scaffold) in scaffolds.iter().enumerate() {
        let (sequence, gaps) = assemble_scaffold(scaffold, &contigs);
        stats.total_gaps += gaps;
        stats.total_length += sequence.len();
        scaffold_lengths.push(sequence.len());

        let header = format!(
            "scaffold_{} contigs={} gaps={}",
            i + 1,
            scaffold.contigs.len(),
            gaps
        );
        writer.write_record(&header, &sequence)?;
    }

    // Calculate N50
    scaffold_lengths.sort_unstable_by(|a, b| b.cmp(a));
    let half = stats.total_length / 2;
    let mut cumsum = 0;
    for len in &scaffold_lengths {
        cumsum += len;
        if cumsum >= half {
            stats.scaffold_n50 = *len;
            break;
        }
    }

    info!(
        "Scaffolding complete: {} scaffolds, N50 = {} bp",
        stats.num_scaffolds, stats.scaffold_n50
    );
    Ok(stats)
}

/// Build scaffolds from filtered links using greedy path extension
fn build_scaffolds(contigs: &[(String, String)], links: &[ContigLink]) -> Vec<Scaffold> {
    let n = contigs.len();
    let mut used: AHashSet<usize> = AHashSet::new();
    let mut scaffolds = Vec::new();

    // Build adjacency list
    let mut adj: AHashMap<usize, Vec<&ContigLink>> = AHashMap::new();
    for link in links {
        adj.entry(link.contig_a).or_default().push(link);
        adj.entry(link.contig_b).or_default().push(link);
    }

    // Sort contigs by length (start with longest)
    let mut contig_order: Vec<usize> = (0..n).collect();
    contig_order.sort_unstable_by_key(|&i| std::cmp::Reverse(contigs[i].1.len()));

    for start in contig_order {
        if used.contains(&start) {
            continue;
        }

        let mut scaffold = Scaffold {
            contigs: vec![(start, Orientation::Forward)],
            gaps: Vec::new(),
        };
        used.insert(start);

        // Extend right
        let mut current = start;
        loop {
            let next_link = adj.get(&current).and_then(|links| {
                links
                    .iter()
                    .filter(|link| {
                        let other = if link.contig_a == current {
                            link.contig_b
                        } else {
                            link.contig_a
                        };
                        !used.contains(&other)
                    })
                    .max_by_key(|link| link.support_count)
            });

            match next_link {
                Some(link) => {
                    let other = if link.contig_a == current {
                        link.contig_b
                    } else {
                        link.contig_a
                    };
                    let orientation = if link.contig_a == current {
                        link.orientation_b
                    } else {
                        link.orientation_a
                    };
                    scaffold.contigs.push((other, orientation));
                    scaffold.gaps.push(link.median_gap().max(1)); // At least 1 N for gap
                    used.insert(other);
                    current = other;
                }
                None => break,
            }
        }

        scaffolds.push(scaffold);
    }

    scaffolds
}

/// Assemble a scaffold sequence from its contigs
fn assemble_scaffold(scaffold: &Scaffold, contigs: &[(String, String)]) -> (String, usize) {
    let mut sequence = String::new();
    let mut gap_count = 0;

    for (i, &(contig_id, orientation)) in scaffold.contigs.iter().enumerate() {
        let contig_seq = &contigs[contig_id].1;

        let seq = match orientation {
            Orientation::Forward => contig_seq.clone(),
            Orientation::Reverse => {
                String::from_utf8(reverse_complement(contig_seq.as_bytes())).unwrap()
            }
        };

        sequence.push_str(&seq);

        // Add gap if not the last contig
        if i < scaffold.gaps.len() {
            let gap_size = scaffold.gaps[i].max(1) as usize;
            sequence.push_str(&"N".repeat(gap_size));
            gap_count += 1;
        }
    }

    (sequence, gap_count)
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

        let hits = index.query(b"ACGTACGTACGTACGTACGT");
        assert!(!hits.is_empty());
    }

    #[test]
    fn test_hash_kmer() {
        let h1 = hash_kmer(b"ACGT");
        let h2 = hash_kmer(b"ACGT");
        assert_eq!(h1, h2);

        let h3 = hash_kmer(b"TGCA");
        assert_ne!(h1, h3);
    }
}
