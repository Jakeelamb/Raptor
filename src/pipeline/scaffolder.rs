//! Scaffolder - Link contigs into scaffolds using paired-end reads
//!
//! Algorithm:
//! 1. Index contigs using minimizers
//! 2. Map paired reads to contigs
//! 3. Build scaffold graph from read pair links
//! 4. Filter edges by minimum support
//! 5. Extract scaffold paths
//! 6. Output scaffolds with gap estimates

use crate::eval::metrics::evaluate_lengths_sorted_desc;
use crate::io::fasta::{open_fasta, FastaWriter};
use crate::io::fastq::{stream_paired_fastq_records_checked, try_open_fastq};
use crate::pipeline::polisher::{
    add_to_pileup, apply_corrections, read_contigs as read_contigs_bytes, write_contigs,
    MinimizerIndex as PolishMinimizerIndex, PileupColumn, PolishStats, ReadMappingConfig,
};
use ahash::{AHashMap, AHashSet};
use std::cmp::Ordering;
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
    pub scaffold_n25: usize,
    pub scaffold_n50: usize,
    pub scaffold_n75: usize,
    pub scaffold_n90: usize,
    pub scaffold_n95: usize,
    pub scaffold_n99: usize,
    pub scaffold_l25: usize,
    pub scaffold_l50: usize,
    pub scaffold_l75: usize,
    pub scaffold_l90: usize,
    pub scaffold_l95: usize,
    pub scaffold_l99: usize,
    pub scaffold_avg_len: f64,
    pub scaffold_median_len: f64,
    pub scaffold_au_n: f64,
    pub longest_scaffold: usize,
    pub scaffolds_ge_1kb: usize,
    pub scaffolds_ge_10kb: usize,
    pub scaffolds_ge_50kb: usize,
    pub scaffolds_ge_100kb: usize,
    pub scaffolds_ge_1mb: usize,
    pub bases_ge_1kb: usize,
    pub bases_ge_10kb: usize,
    pub bases_ge_50kb: usize,
    pub bases_ge_100kb: usize,
    pub bases_ge_1mb: usize,
    pub scaffolds_ge_1kb_frac: f64,
    pub scaffolds_ge_10kb_frac: f64,
    pub scaffolds_ge_50kb_frac: f64,
    pub scaffolds_ge_100kb_frac: f64,
    pub scaffolds_ge_1mb_frac: f64,
    pub bases_ge_1kb_frac: f64,
    pub bases_ge_10kb_frac: f64,
    pub bases_ge_50kb_frac: f64,
    pub bases_ge_100kb_frac: f64,
    pub bases_ge_1mb_frac: f64,
    pub total_length: usize,
    pub total_gaps: usize,
}

/// Orientation of a contig in a scaffold
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
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
        median_i32(&self.gap_estimates)
    }
}

#[inline]
fn midpoint_i32(a: i32, b: i32) -> i32 {
    ((a as i64 + b as i64) / 2) as i32
}

#[inline]
fn median_i32(values: &[i32]) -> i32 {
    if values.is_empty() {
        return 0;
    }
    let mut scratch = values.to_vec();
    median_i32_in_place(&mut scratch)
}

#[inline]
fn median_i32_in_place(values: &mut [i32]) -> i32 {
    debug_assert!(!values.is_empty());
    let mid = values.len() / 2;
    let upper = {
        let (_, upper, _) = values.select_nth_unstable(mid);
        *upper
    };
    if values.len() % 2 == 1 {
        upper
    } else {
        let lower = values[..mid].iter().copied().max().unwrap_or(upper);
        midpoint_i32(lower, upper)
    }
}

#[inline]
fn other_contig_id(link: &ContigLink, current: usize) -> usize {
    if link.contig_a == current {
        link.contig_b
    } else {
        link.contig_a
    }
}

#[inline]
fn compare_extension_candidates(current: usize, a: &ContigLink, b: &ContigLink) -> Ordering {
    a.support_count
        .cmp(&b.support_count)
        .then_with(|| other_contig_id(b, current).cmp(&other_contig_id(a, current)))
        .then_with(|| a.orientation_a.cmp(&b.orientation_a))
        .then_with(|| a.orientation_b.cmp(&b.orientation_b))
}

#[inline]
fn compare_links(a: &ContigLink, b: &ContigLink) -> Ordering {
    a.contig_a
        .cmp(&b.contig_a)
        .then_with(|| a.contig_b.cmp(&b.contig_b))
        .then_with(|| b.support_count.cmp(&a.support_count))
        .then_with(|| a.orientation_a.cmp(&b.orientation_a))
        .then_with(|| a.orientation_b.cmp(&b.orientation_b))
}

#[inline]
fn select_next_link<'a>(
    current: usize,
    links: &'a [&'a ContigLink],
    used: &AHashSet<usize>,
) -> Option<&'a ContigLink> {
    links
        .iter()
        .filter(|link| {
            let other = other_contig_id(link, current);
            !used.contains(&other)
        })
        .copied()
        .max_by(|a, b| compare_extension_candidates(current, a, b))
}

fn update_scaffold_continuity(stats: &mut ScaffoldStats, scaffold_lengths: &mut [usize]) {
    scaffold_lengths.sort_unstable_by(|a, b| b.cmp(a));
    if scaffold_lengths.is_empty() {
        return;
    }

    let continuity = evaluate_lengths_sorted_desc(scaffold_lengths);
    stats.scaffold_n25 = continuity.n25;
    stats.scaffold_n50 = continuity.n50;
    stats.scaffold_n75 = continuity.n75;
    stats.scaffold_n90 = continuity.n90;
    stats.scaffold_n95 = continuity.n95;
    stats.scaffold_n99 = continuity.n99;
    stats.scaffold_l25 = continuity.l25;
    stats.scaffold_l50 = continuity.l50;
    stats.scaffold_l75 = continuity.l75;
    stats.scaffold_l90 = continuity.l90;
    stats.scaffold_l95 = continuity.l95;
    stats.scaffold_l99 = continuity.l99;
    stats.scaffold_avg_len = continuity.avg_length;
    stats.scaffold_median_len = continuity.median_length;
    stats.scaffold_au_n = continuity.au_n;
    stats.longest_scaffold = continuity.longest;
    stats.scaffolds_ge_1kb = continuity.contigs_ge_1kb;
    stats.scaffolds_ge_10kb = continuity.contigs_ge_10kb;
    stats.scaffolds_ge_50kb = continuity.contigs_ge_50kb;
    stats.scaffolds_ge_100kb = continuity.contigs_ge_100kb;
    stats.scaffolds_ge_1mb = continuity.contigs_ge_1mb;
    stats.bases_ge_1kb = continuity.bases_ge_1kb;
    stats.bases_ge_10kb = continuity.bases_ge_10kb;
    stats.bases_ge_50kb = continuity.bases_ge_50kb;
    stats.bases_ge_100kb = continuity.bases_ge_100kb;
    stats.bases_ge_1mb = continuity.bases_ge_1mb;
    stats.scaffolds_ge_1kb_frac = continuity.contigs_ge_1kb_frac;
    stats.scaffolds_ge_10kb_frac = continuity.contigs_ge_10kb_frac;
    stats.scaffolds_ge_50kb_frac = continuity.contigs_ge_50kb_frac;
    stats.scaffolds_ge_100kb_frac = continuity.contigs_ge_100kb_frac;
    stats.scaffolds_ge_1mb_frac = continuity.contigs_ge_1mb_frac;
    stats.bases_ge_1kb_frac = continuity.bases_ge_1kb_frac;
    stats.bases_ge_10kb_frac = continuity.bases_ge_10kb_frac;
    stats.bases_ge_50kb_frac = continuity.bases_ge_50kb_frac;
    stats.bases_ge_100kb_frac = continuity.bases_ge_100kb_frac;
    stats.bases_ge_1mb_frac = continuity.bases_ge_1mb_frac;
}

/// A scaffold is a list of oriented contigs with gaps
#[derive(Debug, Clone)]
pub struct Scaffold {
    pub contigs: Vec<(usize, Orientation)>,
    pub gaps: Vec<i32>,
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
    mapping_config: ReadMappingConfig,
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
    let mut index =
        PolishMinimizerIndex::new(mapping_config.minimizer_k, mapping_config.minimizer_w);
    for (i, (_header, seq)) in contigs.iter().enumerate() {
        index.add_contig(i, seq.as_bytes());
    }

    // Map read pairs and collect links
    info!("Mapping paired-end reads...");
    let mut links: AHashMap<(usize, usize), ContigLink> = AHashMap::new();

    let reader1 = try_open_fastq(reads1_path)?;
    let reader2 = try_open_fastq(reads2_path)?;

    // Estimate insert size from first batch of reads
    let mut insert_sizes: Vec<i32> = Vec::new();
    const INSERT_SAMPLE_SIZE: usize = 10000;
    const INSERT_MEDIAN_REFRESH: usize = 256;
    let mut estimated_insert = 500i32;
    let mut hits = AHashMap::new();
    let mut contig_hits = AHashMap::new();
    let mut reverse_scratch = Vec::new();

    for pair in stream_paired_fastq_records_checked(reader1, reader2) {
        let (r1, r2) = pair?;
        stats.reads_processed += 2;

        let sequence1 = r1.sequence.into_bytes();
        let sequence2 = r2.sequence.into_bytes();
        let scaffold_hit1 = index
            .map_read_with_scaffold_support_thresholds(
                &sequence1,
                &mut hits,
                &mut contig_hits,
                &mut reverse_scratch,
                mapping_config.min_primary_matches,
                mapping_config.min_scaffold_matches,
            )
            .1;
        let scaffold_hit2 = index
            .map_read_with_scaffold_support_thresholds(
                &sequence2,
                &mut hits,
                &mut contig_hits,
                &mut reverse_scratch,
                mapping_config.min_primary_matches,
                mapping_config.min_scaffold_matches,
            )
            .1;

        if let (Some((contig1, _support1, rev1)), Some((contig2, _support2, rev2))) =
            (scaffold_hit1, scaffold_hit2)
        {
            stats.reads_mapped += 2;

            // Skip if both reads map to the same contig (used for insert size estimation)
            if contig1 == contig2 && insert_sizes.len() < INSERT_SAMPLE_SIZE {
                // Rough insert size estimate based on read lengths
                insert_sizes.push((sequence1.len() + sequence2.len()) as i32);
                if insert_sizes.len() == 1
                    || insert_sizes.len() % INSERT_MEDIAN_REFRESH == 0
                    || insert_sizes.len() == INSERT_SAMPLE_SIZE
                {
                    estimated_insert = median_i32(&insert_sizes);
                }
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
                let gap = estimated_insert - (sequence1.len() + sequence2.len()) as i32;
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
    let mut filtered_links: Vec<ContigLink> = links
        .into_values()
        .filter(|link| link.support_count >= min_links)
        .collect();
    filtered_links.sort_unstable_by(compare_links);
    stats.links_after_filter = filtered_links.len();
    info!(
        "Retained {} links with >= {} supporting pairs",
        filtered_links.len(),
        min_links
    );

    write_scaffold_output(&contigs, &filtered_links, output_path, &mut stats)?;
    Ok(stats)
}

fn write_scaffold_output<T: AsRef<[u8]>>(
    contigs: &[(String, T)],
    filtered_links: &[ContigLink],
    output_path: &str,
    stats: &mut ScaffoldStats,
) -> Result<()> {
    let scaffolds = build_scaffolds(contigs, filtered_links);
    stats.num_scaffolds = scaffolds.len();

    info!("Writing {} scaffolds to {}", scaffolds.len(), output_path);
    let mut writer = FastaWriter::try_new(output_path)?;
    let mut scaffold_lengths = Vec::with_capacity(scaffolds.len());

    for (i, scaffold) in scaffolds.iter().enumerate() {
        let (sequence, gaps) = assemble_scaffold(scaffold, contigs);
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

    update_scaffold_continuity(stats, &mut scaffold_lengths);
    info!(
        "Scaffolding complete: {} scaffolds, N50 = {} bp, N90 = {} bp, N95 = {} bp",
        stats.num_scaffolds, stats.scaffold_n50, stats.scaffold_n90, stats.scaffold_n95
    );
    Ok(())
}

pub fn scaffold_and_polish_contigs(
    contigs_path: &str,
    reads1_path: &str,
    reads2_path: &str,
    scaffold_output_path: &str,
    polish_output_path: &str,
    min_links: usize,
    mapping_config: ReadMappingConfig,
) -> Result<(ScaffoldStats, PolishStats)> {
    let mut scaffold_stats = ScaffoldStats::default();
    let mut polish_stats = PolishStats::default();

    info!("Reading contigs from {}", contigs_path);
    let mut contigs = read_contigs_bytes(contigs_path)?;
    scaffold_stats.contigs_input = contigs.len();
    polish_stats.contigs_input = contigs.len();
    info!("Loaded {} contigs", contigs.len());

    if contigs.is_empty() {
        write_scaffold_output(&contigs, &[], scaffold_output_path, &mut scaffold_stats)?;
        info!("Writing polished contigs to {}", polish_output_path);
        write_contigs(&contigs, polish_output_path)?;
        return Ok((scaffold_stats, polish_stats));
    }

    info!("Building shared minimizer index...");
    let mut index =
        PolishMinimizerIndex::new(mapping_config.minimizer_k, mapping_config.minimizer_w);
    for (contig_id, (_, sequence)) in contigs.iter().enumerate() {
        index.add_contig(contig_id, sequence);
    }

    let mut pileups: Vec<Vec<PileupColumn>> = contigs
        .iter()
        .map(|(_, sequence)| vec![PileupColumn::default(); sequence.len()])
        .collect();
    let mut links: AHashMap<(usize, usize), ContigLink> = AHashMap::new();

    let reader1 = try_open_fastq(reads1_path)?;
    let reader2 = try_open_fastq(reads2_path)?;
    let mut hits = AHashMap::new();
    let mut contig_hits = AHashMap::new();
    let mut reverse_scratch = Vec::new();

    let mut insert_sizes = Vec::new();
    const INSERT_SAMPLE_SIZE: usize = 10000;
    const INSERT_MEDIAN_REFRESH: usize = 256;
    let mut estimated_insert = 500i32;

    info!("Mapping paired-end reads once for scaffolding and polishing...");
    for pair in stream_paired_fastq_records_checked(reader1, reader2) {
        let (r1, r2) = pair?;
        scaffold_stats.reads_processed += 2;
        polish_stats.reads_processed += 2;

        let sequence1 = r1.sequence.into_bytes();
        let (mapping1, scaffold_hit1) = index.map_read_with_scaffold_support_thresholds(
            &sequence1,
            &mut hits,
            &mut contig_hits,
            &mut reverse_scratch,
            mapping_config.min_primary_matches,
            mapping_config.min_scaffold_matches,
        );
        if let Some((contig_id, start, is_reverse)) = mapping1 {
            polish_stats.reads_mapped += 1;
            let read_sequence = if is_reverse {
                reverse_scratch.as_slice()
            } else {
                sequence1.as_slice()
            };
            add_to_pileup(&mut pileups[contig_id], start, read_sequence);
        }

        let sequence2 = r2.sequence.into_bytes();
        let (mapping2, scaffold_hit2) = index.map_read_with_scaffold_support_thresholds(
            &sequence2,
            &mut hits,
            &mut contig_hits,
            &mut reverse_scratch,
            mapping_config.min_primary_matches,
            mapping_config.min_scaffold_matches,
        );
        if let Some((contig_id, start, is_reverse)) = mapping2 {
            polish_stats.reads_mapped += 1;
            let read_sequence = if is_reverse {
                reverse_scratch.as_slice()
            } else {
                sequence2.as_slice()
            };
            add_to_pileup(&mut pileups[contig_id], start, read_sequence);
        }

        if let (Some((contig1, _support1, rev1)), Some((contig2, _support2, rev2))) =
            (scaffold_hit1, scaffold_hit2)
        {
            scaffold_stats.reads_mapped += 2;

            if contig1 == contig2 && insert_sizes.len() < INSERT_SAMPLE_SIZE {
                insert_sizes.push((sequence1.len() + sequence2.len()) as i32);
                if insert_sizes.len() == 1
                    || insert_sizes.len() % INSERT_MEDIAN_REFRESH == 0
                    || insert_sizes.len() == INSERT_SAMPLE_SIZE
                {
                    estimated_insert = median_i32(&insert_sizes);
                }
                continue;
            }

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

                let gap = estimated_insert - (sequence1.len() + sequence2.len()) as i32;
                link.gap_estimates.push(gap);
                link.support_count += 1;
            }
        }

        if scaffold_stats.reads_processed % 1_000_000 == 0 {
            info!(
                "  Processed {} million read pairs...",
                scaffold_stats.reads_processed / 2_000_000
            );
        }
    }

    scaffold_stats.links_found = links.len();
    info!("Found {} potential links", links.len());

    let mut filtered_links: Vec<ContigLink> = links
        .into_values()
        .filter(|link| link.support_count >= min_links)
        .collect();
    filtered_links.sort_unstable_by(compare_links);
    scaffold_stats.links_after_filter = filtered_links.len();
    info!(
        "Retained {} links with >= {} supporting pairs",
        filtered_links.len(),
        min_links
    );

    write_scaffold_output(
        &contigs,
        &filtered_links,
        scaffold_output_path,
        &mut scaffold_stats,
    )?;

    info!("Calling consensus...");
    let corrections = apply_corrections(&mut contigs, &pileups, &mut polish_stats);
    info!("Made {} corrections", corrections);

    info!("Writing polished contigs to {}", polish_output_path);
    write_contigs(&contigs, polish_output_path)?;
    info!(
        "Polishing complete: {} corrections ({} substitutions, {} insertions, {} deletions)",
        polish_stats.corrections,
        polish_stats.substitutions,
        polish_stats.insertions,
        polish_stats.deletions
    );

    Ok((scaffold_stats, polish_stats))
}

/// Build scaffolds from filtered links using greedy path extension
fn build_scaffolds<T: AsRef<[u8]>>(contigs: &[(String, T)], links: &[ContigLink]) -> Vec<Scaffold> {
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
    contig_order.sort_unstable_by(|&a, &b| {
        contigs[b]
            .1
            .as_ref()
            .len()
            .cmp(&contigs[a].1.as_ref().len())
            .then_with(|| a.cmp(&b))
    });

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
            let next_link = adj
                .get(&current)
                .and_then(|links| select_next_link(current, links, &used));

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

        // Extend left from the same seed to avoid splitting linear chains
        // when the chosen seed lies in the middle of the true scaffold.
        let mut left_contigs_rev: Vec<(usize, Orientation)> = Vec::new();
        let mut left_gaps_rev: Vec<i32> = Vec::new();
        current = start;
        loop {
            let next_link = adj
                .get(&current)
                .and_then(|links| select_next_link(current, links, &used));

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
                    left_contigs_rev.push((other, orientation));
                    left_gaps_rev.push(link.median_gap().max(1)); // At least 1 N for gap
                    used.insert(other);
                    current = other;
                }
                None => break,
            }
        }

        if !left_contigs_rev.is_empty() {
            left_contigs_rev.reverse();
            left_gaps_rev.reverse();
            left_contigs_rev.extend(scaffold.contigs);
            left_gaps_rev.extend(scaffold.gaps);
            scaffold.contigs = left_contigs_rev;
            scaffold.gaps = left_gaps_rev;
        }

        scaffolds.push(scaffold);
    }

    scaffolds
}

/// Assemble a scaffold sequence from its contigs
fn assemble_scaffold<T: AsRef<[u8]>>(
    scaffold: &Scaffold,
    contigs: &[(String, T)],
) -> (String, usize) {
    let mut sequence = String::new();
    let mut gap_count = 0;

    for (i, &(contig_id, orientation)) in scaffold.contigs.iter().enumerate() {
        let contig_seq = contigs[contig_id].1.as_ref();
        match orientation {
            Orientation::Forward => sequence.push_str(
                std::str::from_utf8(contig_seq).expect("scaffold contigs must remain ASCII DNA"),
            ),
            Orientation::Reverse => {
                let rc = reverse_complement(contig_seq);
                sequence.push_str(
                    std::str::from_utf8(&rc).expect("reverse complement must remain ASCII DNA"),
                );
            }
        }

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
    use rand::rngs::StdRng;
    use rand::seq::SliceRandom;
    use rand::SeedableRng;
    use std::io::Write;
    use tempfile::{NamedTempFile, TempDir};

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"GCGC"), b"GCGC");
    }

    #[test]
    fn test_median_i32_even_sample_uses_midpoint() {
        let values = [100, 200, 300, 400];
        assert_eq!(median_i32(&values), 250);
    }

    #[test]
    fn test_contig_link_median_gap_even_sample_uses_midpoint() {
        let link = ContigLink {
            contig_a: 0,
            contig_b: 1,
            orientation_a: Orientation::Forward,
            orientation_b: Orientation::Forward,
            gap_estimates: vec![100, 200, 300, 400],
            support_count: 4,
        };
        assert_eq!(link.median_gap(), 250);
    }

    #[test]
    fn test_update_scaffold_continuity_populates_extended_metrics() {
        let mut stats = ScaffoldStats::default();
        let mut lengths = vec![100, 50, 25];
        update_scaffold_continuity(&mut stats, &mut lengths);

        assert_eq!(stats.scaffold_n25, 100);
        assert_eq!(stats.scaffold_n50, 100);
        assert_eq!(stats.scaffold_n75, 50);
        assert_eq!(stats.scaffold_n90, 25);
        assert_eq!(stats.scaffold_n95, 25);
        assert_eq!(stats.scaffold_n99, 25);
        assert_eq!(stats.scaffold_l25, 1);
        assert_eq!(stats.scaffold_l50, 1);
        assert_eq!(stats.scaffold_l75, 2);
        assert_eq!(stats.scaffold_l90, 3);
        assert_eq!(stats.scaffold_l95, 3);
        assert_eq!(stats.scaffold_l99, 3);
        assert!((stats.scaffold_avg_len - 58.3333333333).abs() < 1e-9);
        assert_eq!(stats.scaffold_median_len, 50.0);
        assert_eq!(stats.longest_scaffold, 100);
        assert!((stats.scaffold_au_n - 75.0).abs() < 1e-12);
        assert_eq!(stats.scaffolds_ge_1kb, 0);
        assert_eq!(stats.bases_ge_1kb, 0);
        assert_eq!(stats.scaffolds_ge_1kb_frac, 0.0);
        assert_eq!(stats.bases_ge_1kb_frac, 0.0);
        assert_eq!(stats.scaffolds_ge_1mb, 0);
        assert_eq!(stats.bases_ge_1mb, 0);
        assert_eq!(stats.scaffolds_ge_1mb_frac, 0.0);
        assert_eq!(stats.bases_ge_1mb_frac, 0.0);
    }

    #[test]
    fn test_update_scaffold_continuity_reports_length_bucket_metrics() {
        let mut stats = ScaffoldStats::default();
        let mut lengths = vec![2_000, 1_200, 800, 400];
        update_scaffold_continuity(&mut stats, &mut lengths);

        assert_eq!(stats.scaffolds_ge_1kb, 2);
        assert_eq!(stats.scaffolds_ge_10kb, 0);
        assert_eq!(stats.scaffolds_ge_50kb, 0);
        assert_eq!(stats.scaffolds_ge_100kb, 0);
        assert_eq!(stats.scaffolds_ge_1mb, 0);
        assert_eq!(stats.bases_ge_1kb, 3_200);
        assert_eq!(stats.bases_ge_10kb, 0);
        assert_eq!(stats.bases_ge_50kb, 0);
        assert_eq!(stats.bases_ge_100kb, 0);
        assert_eq!(stats.bases_ge_1mb, 0);
        assert!((stats.scaffolds_ge_1kb_frac - 0.5).abs() < 1e-12);
        assert_eq!(stats.scaffolds_ge_10kb_frac, 0.0);
        assert_eq!(stats.scaffolds_ge_50kb_frac, 0.0);
        assert_eq!(stats.scaffolds_ge_100kb_frac, 0.0);
        assert_eq!(stats.scaffolds_ge_1mb_frac, 0.0);
        assert!((stats.bases_ge_1kb_frac - (3_200.0 / 4_400.0)).abs() < 1e-12);
        assert_eq!(stats.bases_ge_10kb_frac, 0.0);
        assert_eq!(stats.bases_ge_50kb_frac, 0.0);
        assert_eq!(stats.bases_ge_100kb_frac, 0.0);
        assert_eq!(stats.bases_ge_1mb_frac, 0.0);
    }

    fn make_link(a: usize, b: usize, support_count: usize) -> ContigLink {
        ContigLink {
            contig_a: a.min(b),
            contig_b: a.max(b),
            orientation_a: Orientation::Forward,
            orientation_b: Orientation::Forward,
            gap_estimates: vec![100],
            support_count,
        }
    }

    fn scaffold_fingerprint(scaffolds: &[Scaffold]) -> Vec<Vec<usize>> {
        scaffolds
            .iter()
            .map(|scaffold| scaffold.contigs.iter().map(|(id, _)| *id).collect())
            .collect()
    }

    #[test]
    fn build_scaffolds_is_stable_under_randomized_link_order() {
        let contigs = vec![
            ("c0".to_string(), "A".repeat(200)),
            ("c1".to_string(), "C".repeat(150)),
            ("c2".to_string(), "G".repeat(150)),
            ("c3".to_string(), "T".repeat(120)),
        ];

        let base_links = vec![
            make_link(0, 1, 10),
            make_link(0, 2, 10),
            make_link(1, 3, 8),
            make_link(2, 3, 8),
        ];

        let baseline = scaffold_fingerprint(&build_scaffolds(&contigs, &base_links));

        for seed in 0u64..64 {
            let mut rng = StdRng::seed_from_u64(seed);
            let mut shuffled = base_links.clone();
            shuffled.shuffle(&mut rng);
            let observed = scaffold_fingerprint(&build_scaffolds(&contigs, &shuffled));
            assert_eq!(observed, baseline);
        }
    }

    #[test]
    fn build_scaffolds_extends_both_directions_from_internal_seed() {
        let contigs = vec![
            ("c0".to_string(), "A".repeat(200)),
            ("c1".to_string(), "C".repeat(500)),
            ("c2".to_string(), "G".repeat(180)),
        ];
        let links = vec![make_link(0, 1, 10), make_link(1, 2, 8)];

        let scaffolds = build_scaffolds(&contigs, &links);
        assert_eq!(scaffolds.len(), 1);
        assert_eq!(
            scaffolds[0].contigs,
            vec![
                (2, Orientation::Forward),
                (1, Orientation::Forward),
                (0, Orientation::Forward)
            ]
        );
        assert_eq!(scaffolds[0].gaps, vec![100, 100]);
    }

    #[test]
    fn scaffold_contigs_returns_error_for_invalid_output_path() {
        let mut contigs = NamedTempFile::new().unwrap();
        writeln!(contigs, ">contig_1").unwrap();
        writeln!(contigs, "ACGTACGTACGT").unwrap();

        let mut reads1 = NamedTempFile::new().unwrap();
        writeln!(reads1, "@r1").unwrap();
        writeln!(reads1, "ACGTACGT").unwrap();
        writeln!(reads1, "+").unwrap();
        writeln!(reads1, "IIIIIIII").unwrap();

        let mut reads2 = NamedTempFile::new().unwrap();
        writeln!(reads2, "@r2").unwrap();
        writeln!(reads2, "ACGTACGT").unwrap();
        writeln!(reads2, "+").unwrap();
        writeln!(reads2, "IIIIIIII").unwrap();

        let output_dir = TempDir::new().unwrap();
        let result = scaffold_contigs(
            contigs.path().to_str().unwrap(),
            reads1.path().to_str().unwrap(),
            reads2.path().to_str().unwrap(),
            output_dir.path().to_str().unwrap(),
            1,
            ReadMappingConfig::default(),
        );
        assert!(result.is_err());
    }

    #[test]
    fn scaffold_contigs_respects_configured_scaffold_support_threshold() {
        let contigs = NamedTempFile::new().unwrap();
        std::fs::write(
            contigs.path(),
            b">contig_1\nACGTTGCATGCAAGTCGATCGTACCGTTAAGGCTAACGTA\n>contig_2\nTTAACCGGATCCGTTAGGCCAATTCGATGGCCTTAAGGCC\n",
        )
        .unwrap();

        let mut reads1 = NamedTempFile::new().unwrap();
        let mut reads2 = NamedTempFile::new().unwrap();
        let quality = "I".repeat(30);
        for idx in 0..2 {
            writeln!(reads1, "@r{}_1", idx).unwrap();
            writeln!(reads1, "ACGTTGCATGCAAGTCGATCGTACCGTTAA").unwrap();
            writeln!(reads1, "+").unwrap();
            writeln!(reads1, "{}", quality).unwrap();

            writeln!(reads2, "@r{}_2", idx).unwrap();
            writeln!(reads2, "TTAACCGGATCCGTTAGGCCAATTCGATGG").unwrap();
            writeln!(reads2, "+").unwrap();
            writeln!(reads2, "{}", quality).unwrap();
        }

        let low_threshold_output = NamedTempFile::new().unwrap();
        let high_threshold_output = NamedTempFile::new().unwrap();
        let low_threshold_stats = scaffold_contigs(
            contigs.path().to_str().unwrap(),
            reads1.path().to_str().unwrap(),
            reads2.path().to_str().unwrap(),
            low_threshold_output.path().to_str().unwrap(),
            1,
            ReadMappingConfig::default(),
        )
        .unwrap();
        let high_threshold_stats = scaffold_contigs(
            contigs.path().to_str().unwrap(),
            reads1.path().to_str().unwrap(),
            reads2.path().to_str().unwrap(),
            high_threshold_output.path().to_str().unwrap(),
            1,
            ReadMappingConfig {
                min_scaffold_matches: 64,
                ..ReadMappingConfig::default()
            },
        )
        .unwrap();

        assert!(low_threshold_stats.links_after_filter > 0);
        assert_eq!(high_threshold_stats.links_after_filter, 0);
    }

    #[test]
    fn scaffold_and_polish_contigs_matches_separate_passes_on_simple_dataset() {
        let contig1_correct = b"ACGTTGCATGCAAGTCGATCGTACCGTTAAGGCTAACGTA".to_vec();
        let mut contig1_mutated = contig1_correct.clone();
        contig1_mutated[10] = b'T';
        let contig2 = b"TTAACCGGATCCGTTAGGCCAATTCGATGGCCTTAAGGCC".to_vec();

        let mut contigs = NamedTempFile::new().unwrap();
        writeln!(contigs, ">contig_1").unwrap();
        writeln!(contigs, "{}", String::from_utf8_lossy(&contig1_mutated)).unwrap();
        writeln!(contigs, ">contig_2").unwrap();
        writeln!(contigs, "{}", String::from_utf8_lossy(&contig2)).unwrap();

        let mut reads1 = NamedTempFile::new().unwrap();
        let mut reads2 = NamedTempFile::new().unwrap();
        let read1 = &contig1_correct[..30];
        let read2 = &contig2[..30];
        let quality = "I".repeat(read1.len());
        for idx in 0..5 {
            writeln!(reads1, "@r{}_1", idx).unwrap();
            writeln!(reads1, "{}", String::from_utf8_lossy(read1)).unwrap();
            writeln!(reads1, "+").unwrap();
            writeln!(reads1, "{}", quality).unwrap();

            writeln!(reads2, "@r{}_2", idx).unwrap();
            writeln!(reads2, "{}", String::from_utf8_lossy(read2)).unwrap();
            writeln!(reads2, "+").unwrap();
            writeln!(reads2, "{}", quality).unwrap();
        }

        let separate_scaffold = NamedTempFile::new().unwrap();
        let separate_polish = NamedTempFile::new().unwrap();
        let combined_scaffold = NamedTempFile::new().unwrap();
        let combined_polish = NamedTempFile::new().unwrap();

        let separate_scaffold_stats = scaffold_contigs(
            contigs.path().to_str().unwrap(),
            reads1.path().to_str().unwrap(),
            reads2.path().to_str().unwrap(),
            separate_scaffold.path().to_str().unwrap(),
            3,
            ReadMappingConfig::default(),
        )
        .unwrap();
        let separate_polish_stats = crate::pipeline::polisher::polish_contigs(
            contigs.path().to_str().unwrap(),
            reads1.path().to_str().unwrap(),
            Some(reads2.path().to_str().unwrap()),
            separate_polish.path().to_str().unwrap(),
            1,
            ReadMappingConfig::default(),
        )
        .unwrap();

        let (combined_scaffold_stats, combined_polish_stats) = scaffold_and_polish_contigs(
            contigs.path().to_str().unwrap(),
            reads1.path().to_str().unwrap(),
            reads2.path().to_str().unwrap(),
            combined_scaffold.path().to_str().unwrap(),
            combined_polish.path().to_str().unwrap(),
            3,
            ReadMappingConfig::default(),
        )
        .unwrap();

        assert_eq!(
            combined_scaffold_stats.links_after_filter,
            separate_scaffold_stats.links_after_filter
        );
        assert_eq!(
            combined_scaffold_stats.num_scaffolds,
            separate_scaffold_stats.num_scaffolds
        );
        assert_eq!(
            combined_polish_stats.corrections,
            separate_polish_stats.corrections
        );
        assert_eq!(
            combined_polish_stats.reads_mapped,
            separate_polish_stats.reads_mapped
        );
        assert!(combined_polish_stats.reads_mapped > 0);

        let separate_scaffold_output = std::fs::read_to_string(separate_scaffold.path()).unwrap();
        let combined_scaffold_output = std::fs::read_to_string(combined_scaffold.path()).unwrap();
        assert_eq!(combined_scaffold_output, separate_scaffold_output);

        let separate_polish_output = std::fs::read_to_string(separate_polish.path()).unwrap();
        let combined_polish_output = std::fs::read_to_string(combined_polish.path()).unwrap();
        assert_eq!(combined_polish_output, separate_polish_output);
        assert!(combined_polish_output.contains(">contig_1"));
    }
}
