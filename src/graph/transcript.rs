use crate::eval::metrics::{evaluate_lengths, BaseComposition};
use crate::graph::assembler::Contig;
use crate::graph::isoform_traverse::TranscriptPath;
use crate::kmer::rle;
use petgraph::graphmap::DiGraphMap;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Represents a transcript with its sequence, path, and metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Transcript {
    /// Unique identifier for the transcript
    pub id: usize,

    /// Sequence of the transcript
    pub sequence: String,

    /// Path of segment IDs that make up this transcript
    pub path: Vec<usize>,

    /// Confidence score (0.0-1.0) for this transcript
    pub confidence: f64,

    /// Length of the transcript in nucleotides
    pub length: usize,

    /// Strand: '+' for forward, '-' for reverse
    pub strand: char,

    /// Optional Transcripts Per Million expression value
    pub tpm: Option<f64>,

    /// Splicing pattern classification (e.g., "linear", "skipping")
    pub splicing: String,
}

impl Transcript {
    /// Create a new transcript
    pub fn new(id: usize, sequence: String, path: Vec<usize>, confidence: f64) -> Self {
        let length = sequence.len();
        Transcript {
            id,
            sequence,
            path,
            confidence,
            length,
            strand: '+', // Default to forward strand
            tpm: None,
            splicing: "unknown".to_string(),
        }
    }

    /// Set the strand for this transcript
    pub fn with_strand(mut self, strand: char) -> Self {
        self.strand = strand;
        self
    }

    /// Set the TPM value for this transcript
    pub fn with_tpm(mut self, tpm: f64) -> Self {
        self.tpm = Some(tpm);
        self
    }

    /// Set the splicing classification for this transcript
    pub fn with_splicing(mut self, splicing: String) -> Self {
        self.splicing = splicing;
        self
    }
}

/// Stitch together a transcript sequence from a path of contigs
pub fn stitch_isoform(
    contigs: &[Contig],
    path: &[usize],
    overlaps: &[(usize, usize, usize)],
) -> String {
    let contig_sequences: HashMap<usize, &str> = contigs
        .iter()
        .map(|contig| (contig.id, contig.sequence.as_str()))
        .collect();
    let overlap_map = build_overlap_map(overlaps);

    stitch_isoform_with_maps(&contig_sequences, path, &overlap_map)
}

#[inline]
fn build_overlap_map(overlaps: &[(usize, usize, usize)]) -> HashMap<(usize, usize), usize> {
    let mut overlap_map: HashMap<(usize, usize), usize> = HashMap::with_capacity(overlaps.len());
    for &(from, to, overlap) in overlaps {
        overlap_map
            .entry((from, to))
            .and_modify(|existing| *existing = (*existing).max(overlap))
            .or_insert(overlap);
    }
    overlap_map
}

#[inline]
fn stitch_isoform_with_maps(
    contig_sequences: &HashMap<usize, &str>,
    path: &[usize],
    overlaps: &HashMap<(usize, usize), usize>,
) -> String {
    if path.is_empty() {
        return String::new();
    }

    // Start with the first path node that resolves to a contig sequence.
    let mut first_idx = None;
    let mut sequence = String::new();
    for (idx, &node_id) in path.iter().enumerate() {
        if let Some(contig_seq) = contig_sequences.get(&node_id) {
            sequence.push_str(contig_seq);
            first_idx = Some(idx);
            break;
        }
    }
    let Some(start_idx) = first_idx else {
        return String::new();
    };
    let mut prev_id = path[start_idx];

    // Stitch together subsequent contigs, accounting for overlaps
    for &curr_id in path.iter().skip(start_idx + 1) {
        let Some(curr_seq) = contig_sequences.get(&curr_id) else {
            continue;
        };

        // Look up the overlap between these contigs
        let overlap_len = overlaps.get(&(prev_id, curr_id)).copied().unwrap_or(0);

        if overlap_len > 0 && overlap_len < curr_seq.len() {
            // Add only the non-overlapping part of the current contig
            sequence.push_str(&curr_seq[overlap_len..]);
        } else if overlap_len == 0 {
            // No overlap found, just append the entire sequence
            sequence.push_str(curr_seq);
        }
        // If overlap_len >= contig length, nothing new to add
        prev_id = curr_id;
    }

    sequence
}

/// Detect alternative splicing events in a transcript path
pub fn detect_splicing(path: &[usize], graph: &DiGraphMap<usize, f32>) -> String {
    let mut events = vec![];

    // Check for exon skipping by looking for edges that skip nodes in the path
    for window in path.windows(3) {
        if let [a, _b, c] = *window {
            if graph.contains_edge(a, c) {
                events.push("skipping");
            }
        }
    }

    // Check for alternative 5' splice sites
    for i in 0..path.len().saturating_sub(1) {
        let current = path[i];

        // Check if multiple outgoing edges from current node
        let out_neighbors: Vec<_> = graph
            .neighbors_directed(current, petgraph::Direction::Outgoing)
            .collect();
        if out_neighbors.len() > 1 {
            events.push("alt_5prime");
        }
    }

    // Check for alternative 3' splice sites
    for i in 1..path.len() {
        let current = path[i];

        // Check if multiple incoming edges to current node
        let in_neighbors: Vec<_> = graph
            .neighbors_directed(current, petgraph::Direction::Incoming)
            .collect();
        if in_neighbors.len() > 1 {
            events.push("alt_3prime");
        }
    }

    // Remove duplicates and sort
    events.sort();
    events.dedup();

    if events.is_empty() {
        "linear".into()
    } else {
        events.join(",")
    }
}

/// Assemble transcripts from multiple paths
pub fn assemble_transcripts(
    paths: &[TranscriptPath],
    contigs: &[Contig],
    overlaps: &[(usize, usize, usize)],
    graph: Option<&DiGraphMap<usize, f32>>,
) -> Vec<Transcript> {
    let mut transcripts = Vec::new();

    let contig_sequences: HashMap<usize, &str> = contigs
        .iter()
        .map(|contig| (contig.id, contig.sequence.as_str()))
        .collect();
    let overlap_map = build_overlap_map(overlaps);

    // Process each path
    for (i, path) in paths.iter().enumerate() {
        let sequence = stitch_isoform_with_maps(&contig_sequences, &path.nodes, &overlap_map);

        // Detect splicing events if graph is provided
        let splicing = if let Some(g) = graph {
            detect_splicing(&path.nodes, g)
        } else {
            "unknown".to_string()
        };

        transcripts.push(Transcript {
            id: i + 1, // 1-based IDs for transcripts
            sequence: sequence.clone(),
            path: path.nodes.clone(),
            confidence: path.confidence as f64,
            length: sequence.len(),
            strand: '+', // Default to forward strand
            tpm: None,   // No expression value yet
            splicing,    // Detected splicing events
        });
    }

    transcripts
}

/// Generate a FASTA record from a transcript
pub fn transcript_to_fasta(transcript: &Transcript) -> String {
    format!(
        ">transcript_{} length={} confidence={:.3} path={}\n{}\n",
        transcript.id,
        transcript.length,
        transcript.confidence,
        transcript
            .path
            .iter()
            .map(|id| id.to_string())
            .collect::<Vec<_>>()
            .join(","),
        transcript.sequence
    )
}

/// Generate a GFA path record from a transcript
pub fn transcript_to_gfa_path(transcript: &Transcript) -> String {
    let segments = transcript
        .path
        .iter()
        .map(|&id| format!("contig_{}", id + 1))
        .collect::<Vec<_>>()
        .join(",");

    let cigar = "*"; // Placeholder CIGAR string

    format!(
        "P\ttranscript_{}\t{}\t{}\trc:f:{:.3}",
        transcript.id, segments, cigar, transcript.confidence
    )
}

/// Calculate statistics for a collection of transcripts
pub fn calculate_transcript_stats(transcripts: &[Transcript]) -> HashMap<String, f64> {
    let mut stats = HashMap::with_capacity(40);

    // Basic counts
    stats.insert("count".to_string(), transcripts.len() as f64);

    if transcripts.is_empty() {
        for key in [
            "total_length",
            "mean_length",
            "min_length",
            "max_length",
            "n50",
            "n75",
            "n90",
            "n95",
            "n99",
            "l50",
            "l75",
            "l90",
            "l95",
            "l99",
            "au_n",
            "contigs_ge_1kb",
            "contigs_ge_10kb",
            "contigs_ge_50kb",
            "contigs_ge_100kb",
            "bases_ge_1kb",
            "bases_ge_10kb",
            "bases_ge_50kb",
            "bases_ge_100kb",
            "non_finite_confidence_count",
            "mean_confidence",
            "min_confidence",
            "max_confidence",
            "gc_content",
            "gc_content_acgt",
            "n_content",
            "ambiguous_content",
            "acgt_bases",
            "n_bases",
            "ambiguous_bases",
            "mean_rle_ratio",
        ] {
            stats.insert(key.to_string(), 0.0);
        }
        return stats;
    }

    let mut lengths: Vec<usize> = Vec::with_capacity(transcripts.len());
    let mut min_length = usize::MAX;
    let mut max_length = 0usize;
    let mut finite_confidence_count = 0usize;
    let mut finite_confidence_sum = 0.0;
    let mut min_confidence = f64::INFINITY;
    let mut max_confidence = f64::NEG_INFINITY;
    let mut composition = BaseComposition::default();
    let mut total_sequence_bases = 0usize;
    let mut total_rle_ratio = 0.0;

    for transcript in transcripts {
        let length = transcript.length;
        lengths.push(length);
        min_length = min_length.min(length);
        max_length = max_length.max(length);

        let confidence = transcript.confidence;
        if confidence.is_finite() {
            finite_confidence_count += 1;
            finite_confidence_sum += confidence;
            min_confidence = min_confidence.min(confidence);
            max_confidence = max_confidence.max(confidence);
        }

        let sequence = transcript.sequence.as_bytes();
        composition.add_sequence(sequence);
        total_sequence_bases = total_sequence_bases.saturating_add(sequence.len());

        let compressed_len = rle::rle_encode(&transcript.sequence).len();
        total_rle_ratio += if sequence.is_empty() {
            1.0
        } else {
            compressed_len as f64 / sequence.len() as f64
        };
    }

    let length_metrics = evaluate_lengths(&lengths);
    stats.insert(
        "total_length".to_string(),
        length_metrics.total_bases as f64,
    );
    stats.insert("mean_length".to_string(), length_metrics.avg_length);
    stats.insert("min_length".to_string(), min_length as f64);
    stats.insert("max_length".to_string(), max_length as f64);
    stats.insert("n50".to_string(), length_metrics.n50 as f64);
    stats.insert("n75".to_string(), length_metrics.n75 as f64);
    stats.insert("n90".to_string(), length_metrics.n90 as f64);
    stats.insert("n95".to_string(), length_metrics.n95 as f64);
    stats.insert("n99".to_string(), length_metrics.n99 as f64);
    stats.insert("l50".to_string(), length_metrics.l50 as f64);
    stats.insert("l75".to_string(), length_metrics.l75 as f64);
    stats.insert("l90".to_string(), length_metrics.l90 as f64);
    stats.insert("l95".to_string(), length_metrics.l95 as f64);
    stats.insert("l99".to_string(), length_metrics.l99 as f64);
    stats.insert("au_n".to_string(), length_metrics.au_n);
    stats.insert(
        "contigs_ge_1kb".to_string(),
        length_metrics.contigs_ge_1kb as f64,
    );
    stats.insert(
        "contigs_ge_10kb".to_string(),
        length_metrics.contigs_ge_10kb as f64,
    );
    stats.insert(
        "contigs_ge_50kb".to_string(),
        length_metrics.contigs_ge_50kb as f64,
    );
    stats.insert(
        "contigs_ge_100kb".to_string(),
        length_metrics.contigs_ge_100kb as f64,
    );
    stats.insert(
        "bases_ge_1kb".to_string(),
        length_metrics.bases_ge_1kb as f64,
    );
    stats.insert(
        "bases_ge_10kb".to_string(),
        length_metrics.bases_ge_10kb as f64,
    );
    stats.insert(
        "bases_ge_50kb".to_string(),
        length_metrics.bases_ge_50kb as f64,
    );
    stats.insert(
        "bases_ge_100kb".to_string(),
        length_metrics.bases_ge_100kb as f64,
    );

    let non_finite_confidence_count = transcripts.len() - finite_confidence_count;
    stats.insert(
        "non_finite_confidence_count".to_string(),
        non_finite_confidence_count as f64,
    );

    if finite_confidence_count == 0 {
        stats.insert("mean_confidence".to_string(), 0.0);
        stats.insert("min_confidence".to_string(), 0.0);
        stats.insert("max_confidence".to_string(), 0.0);
    } else {
        let mean_confidence = finite_confidence_sum / finite_confidence_count as f64;

        stats.insert("mean_confidence".to_string(), mean_confidence);
        stats.insert("min_confidence".to_string(), min_confidence);
        stats.insert("max_confidence".to_string(), max_confidence);
    }

    let gc_content = if total_sequence_bases > 0 {
        composition.gc_bases as f64 / total_sequence_bases as f64
    } else {
        0.0
    };
    stats.insert("gc_content".to_string(), gc_content);
    stats.insert("gc_content_acgt".to_string(), composition.gc_content());
    stats.insert(
        "n_content".to_string(),
        composition.n_content(total_sequence_bases),
    );
    stats.insert(
        "ambiguous_content".to_string(),
        composition.ambiguous_content(total_sequence_bases),
    );
    stats.insert("acgt_bases".to_string(), composition.acgt_bases as f64);
    stats.insert("n_bases".to_string(), composition.n_bases as f64);
    stats.insert(
        "ambiguous_bases".to_string(),
        composition.ambiguous_bases as f64,
    );

    let mean_rle_ratio = total_rle_ratio / transcripts.len() as f64;
    stats.insert("mean_rle_ratio".to_string(), mean_rle_ratio);

    stats
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stitch_isoform() {
        // Create test contigs
        let contigs = vec![
            Contig {
                id: 0,
                sequence: "ATCGATCG".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 1,
                sequence: "GATCGTTA".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 2,
                sequence: "GTTACGTA".to_string(),
                kmer_path: vec![],
            },
        ];

        // Create overlaps
        // ATCGATCG
        //     GATCGTTA
        //         GTTACGTA
        let overlaps = vec![
            (0, 1, 4), // 4-base overlap between contig 0 and 1
            (1, 2, 4), // 4-base overlap between contig 1 and 2
        ];

        // Test stitching
        let path = vec![0, 1, 2];
        let stitched = stitch_isoform(&contigs, &path, &overlaps);

        // The actual implementation produces this due to how it handles overlaps
        assert_eq!(stitched, "ATCGATCGGTTACGTA");

        // Test with no overlaps
        let no_overlaps: Vec<(usize, usize, usize)> = vec![];
        let stitched_no_overlap = stitch_isoform(&contigs, &path, &no_overlaps);

        // Expected: ATCGATCGGATCGTTAGTTACGTA (just concatenated)
        assert_eq!(stitched_no_overlap, "ATCGATCGGATCGTTAGTTACGTA");
    }

    #[test]
    fn test_assemble_transcripts() {
        // Create test contigs
        let contigs = vec![
            Contig {
                id: 0,
                sequence: "ATCGATCG".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 1,
                sequence: "GATCGTTA".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 2,
                sequence: "GTTACGTA".to_string(),
                kmer_path: vec![],
            },
        ];

        // Create overlaps
        let overlaps = vec![(0, 1, 4), (1, 2, 4)];

        // Create transcript paths
        let paths = vec![
            TranscriptPath {
                nodes: vec![0, 1],
                confidence: 0.9,
                length: 12,
            },
            TranscriptPath {
                nodes: vec![1, 2],
                confidence: 0.8,
                length: 12,
            },
        ];

        // Assemble transcripts
        let transcripts = assemble_transcripts(&paths, &contigs, &overlaps, None);

        // Should produce two transcripts
        assert_eq!(transcripts.len(), 2);

        // Check first transcript
        assert_eq!(transcripts[0].id, 1);
        assert_eq!(transcripts[0].sequence, "ATCGATCGGTTA");
        assert_eq!(transcripts[0].path, vec![0, 1]);

        // Check with approximate equality due to f32 to f64 conversion
        let confidence_diff = (transcripts[0].confidence - 0.9).abs();
        assert!(
            confidence_diff < 0.001,
            "Confidence should be approximately 0.9, got {}",
            transcripts[0].confidence
        );

        assert_eq!(transcripts[0].length, 12);

        // Check second transcript
        assert_eq!(transcripts[1].id, 2);
        assert_eq!(transcripts[1].sequence, "GATCGTTACGTA");
        assert_eq!(transcripts[1].path, vec![1, 2]);

        // Check with approximate equality due to f32 to f64 conversion
        let confidence_diff = (transcripts[1].confidence - 0.8).abs();
        assert!(
            confidence_diff < 0.001,
            "Confidence should be approximately 0.8, got {}",
            transcripts[1].confidence
        );

        assert_eq!(transcripts[1].length, 12);
    }

    #[test]
    fn test_stitch_isoform_handles_sparse_contig_ids() {
        let contigs = vec![
            Contig {
                id: 0,
                sequence: "ATCG".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 2,
                sequence: "CGTT".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 5,
                sequence: "TTAA".to_string(),
                kmer_path: vec![],
            },
        ];
        let overlaps = vec![(0, 2, 2), (2, 5, 2)];
        let path = vec![0, 2, 5];

        let stitched = stitch_isoform(&contigs, &path, &overlaps);
        assert_eq!(stitched, "ATCGTTAA");
    }

    #[test]
    fn test_assemble_transcripts_handles_sparse_contig_ids() {
        let contigs = vec![
            Contig {
                id: 1,
                sequence: "AAAC".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 4,
                sequence: "ACGT".to_string(),
                kmer_path: vec![],
            },
            Contig {
                id: 9,
                sequence: "GTTT".to_string(),
                kmer_path: vec![],
            },
        ];
        let overlaps = vec![(1, 4, 2), (4, 9, 2)];
        let paths = vec![TranscriptPath {
            nodes: vec![1, 4, 9],
            confidence: 0.75,
            length: 10,
        }];

        let transcripts = assemble_transcripts(&paths, &contigs, &overlaps, None);
        assert_eq!(transcripts.len(), 1);
        assert_eq!(transcripts[0].sequence, "AAACGTTT");
        assert_eq!(transcripts[0].path, vec![1, 4, 9]);
    }

    #[test]
    fn test_calculate_transcript_stats_ignores_non_finite_confidence_values() {
        let transcripts = vec![
            Transcript::new(1, "AAAA".to_string(), vec![0], 0.8),
            Transcript::new(2, "CCCC".to_string(), vec![1], f64::NAN),
            Transcript::new(3, "GGGG".to_string(), vec![2], 0.2),
        ];

        let stats = calculate_transcript_stats(&transcripts);
        assert_eq!(stats.get("non_finite_confidence_count").copied(), Some(1.0));
        assert_eq!(stats.get("min_confidence").copied(), Some(0.2));
        assert_eq!(stats.get("max_confidence").copied(), Some(0.8));
        assert_eq!(stats.get("mean_confidence").copied(), Some(0.5));
    }

    #[test]
    fn test_calculate_transcript_stats_handles_all_non_finite_confidence_values() {
        let transcripts = vec![
            Transcript::new(1, "AAAA".to_string(), vec![0], f64::NAN),
            Transcript::new(2, "CCCC".to_string(), vec![1], f64::INFINITY),
            Transcript::new(3, "GGGG".to_string(), vec![2], f64::NEG_INFINITY),
        ];

        let stats = calculate_transcript_stats(&transcripts);
        assert_eq!(stats.get("non_finite_confidence_count").copied(), Some(3.0));
        assert_eq!(stats.get("min_confidence").copied(), Some(0.0));
        assert_eq!(stats.get("max_confidence").copied(), Some(0.0));
        assert_eq!(stats.get("mean_confidence").copied(), Some(0.0));
    }

    #[test]
    fn test_calculate_transcript_stats_empty_input_has_complete_zero_schema() {
        let stats = calculate_transcript_stats(&[]);

        for key in [
            "count",
            "total_length",
            "mean_length",
            "min_length",
            "max_length",
            "n50",
            "n75",
            "n90",
            "n95",
            "n99",
            "l50",
            "l75",
            "l90",
            "l95",
            "l99",
            "au_n",
            "contigs_ge_1kb",
            "contigs_ge_10kb",
            "contigs_ge_50kb",
            "contigs_ge_100kb",
            "bases_ge_1kb",
            "bases_ge_10kb",
            "bases_ge_50kb",
            "bases_ge_100kb",
            "non_finite_confidence_count",
            "mean_confidence",
            "min_confidence",
            "max_confidence",
            "gc_content",
            "gc_content_acgt",
            "n_content",
            "ambiguous_content",
            "acgt_bases",
            "n_bases",
            "ambiguous_bases",
            "mean_rle_ratio",
        ] {
            assert_eq!(
                stats.get(key).copied(),
                Some(0.0),
                "expected metric `{}` to be present with zero default",
                key
            );
        }
    }

    #[test]
    fn test_calculate_transcript_stats_reports_base_composition_breakdown() {
        let transcripts = vec![Transcript::new(1, "GCNRYat".to_string(), vec![0], 0.5)];
        let stats = calculate_transcript_stats(&transcripts);

        assert_eq!(stats.get("acgt_bases").copied(), Some(4.0));
        assert_eq!(stats.get("n_bases").copied(), Some(1.0));
        assert_eq!(stats.get("ambiguous_bases").copied(), Some(2.0));
        assert!((stats.get("gc_content").copied().unwrap_or(0.0) - (2.0 / 7.0)).abs() < 1e-12);
        assert!((stats.get("gc_content_acgt").copied().unwrap_or(0.0) - 0.5).abs() < 1e-12);
        assert!((stats.get("n_content").copied().unwrap_or(0.0) - (1.0 / 7.0)).abs() < 1e-12);
        assert!(
            (stats.get("ambiguous_content").copied().unwrap_or(0.0) - (2.0 / 7.0)).abs() < 1e-12
        );
    }
}
