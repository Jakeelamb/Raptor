use crate::graph::isoform_graph::IsoformGraph;
use crate::graph::isoform_traverse::{find_directed_paths, TranscriptPath};
use crate::graph::transcript::Transcript;
use rayon::prelude::*;
use std::collections::HashMap;

const START_NODE_CHUNK_SIZE: usize = 64;

/// Process paths in parallel using multiple threads
pub fn parallel_path_discovery(
    graph: &IsoformGraph,
    start_nodes: &[usize],
    end_nodes: &[usize],
    max_depth: usize,
) -> Vec<TranscriptPath> {
    if start_nodes.is_empty() {
        return Vec::new();
    }

    // Sort starts once so output is deterministic regardless of caller order.
    let mut ordered_starts = start_nodes.to_vec();
    ordered_starts.sort_unstable();
    ordered_starts.dedup();

    // Fixed-size chunking keeps output independent of runtime thread count.
    let mut chunked_results: Vec<Vec<TranscriptPath>> = ordered_starts
        .par_chunks(START_NODE_CHUNK_SIZE)
        .with_min_len(1)
        .map(|chunk| find_directed_paths(graph, chunk, end_nodes, max_depth))
        .collect();

    let total_paths: usize = chunked_results.iter().map(Vec::len).sum();
    let mut all_paths = Vec::with_capacity(total_paths);
    for paths in chunked_results.iter_mut() {
        all_paths.append(paths);
    }
    all_paths
}

/// Process transcript assembly in parallel
pub fn parallel_transcript_assembly(
    paths: &[TranscriptPath],
    contigs: &HashMap<usize, String>,
    min_confidence: f32,
) -> Vec<Transcript> {
    // Filter paths by confidence
    let filtered_paths: Vec<&TranscriptPath> = paths
        .iter()
        .filter(|path| path.confidence >= min_confidence)
        .collect();

    // Use reference directly - no need to clone into Arc
    // The contigs map is only read, not modified
    let contigs_ref = contigs;

    // Process in parallel - use reference directly
    filtered_paths
        .par_iter()
        .enumerate()
        .map(|(idx, path)| assemble_single_transcript(idx, path, contigs_ref))
        .collect()
}

/// Assemble a single transcript from a path
fn assemble_single_transcript(
    id: usize,
    path: &TranscriptPath,
    contigs: &HashMap<usize, String>,
) -> Transcript {
    let mut sequence = String::new();

    // Stitch together contigs along the path
    for (i, &node_id) in path.nodes.iter().enumerate() {
        if let Some(contig_seq) = contigs.get(&node_id) {
            if i == 0 {
                // First contig is added in full
                sequence.push_str(contig_seq);
            } else {
                // Subsequent contigs might overlap with previous
                // For simplicity, we're using a fixed overlap of 20 bp
                // In a real implementation, you'd use the actual overlaps
                let overlap = 20.min(contig_seq.len());

                if contig_seq.len() > overlap {
                    sequence.push_str(&contig_seq[overlap..]);
                }
            }
        }
    }

    // Create transcript with appropriate metadata
    Transcript::new(id, sequence, path.nodes.clone(), path.confidence as f64)
}

/// Filter transcript paths by confidence score in parallel
pub fn parallel_path_filtering(
    paths: &[TranscriptPath],
    min_confidence: f32,
) -> Vec<TranscriptPath> {
    paths
        .par_iter()
        .filter(|path| path.confidence >= min_confidence)
        .cloned()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::graphmap::DiGraphMap;
    use rand::rngs::StdRng;
    use rand::seq::SliceRandom;
    use rand::SeedableRng;
    use rayon::ThreadPoolBuilder;

    fn summarize_paths(paths: Vec<TranscriptPath>) -> Vec<(Vec<usize>, u32, usize)> {
        paths
            .into_iter()
            .map(|path| (path.nodes, path.confidence.to_bits(), path.length))
            .collect()
    }

    #[test]
    fn test_parallel_assembly() {
        // Create test contigs
        let mut contigs = HashMap::new();
        contigs.insert(1, "ATCGATCG".to_string());
        contigs.insert(2, "ATCGGCTA".to_string());
        contigs.insert(3, "GCTATAGC".to_string());

        // Create test paths
        let path1 = TranscriptPath {
            nodes: vec![1, 2],
            confidence: 0.9,
            length: 16,
        };

        let path2 = TranscriptPath {
            nodes: vec![2, 3],
            confidence: 0.8,
            length: 16,
        };

        let paths = vec![path1, path2];

        // Test parallel assembly
        let transcripts = parallel_transcript_assembly(&paths, &contigs, 0.7);

        assert_eq!(transcripts.len(), 2);
        assert_eq!(transcripts[0].path, vec![1, 2]);
        assert_eq!(transcripts[1].path, vec![2, 3]);
    }

    #[test]
    fn parallel_path_discovery_is_stable_under_start_permutation_and_thread_count() {
        let mut graph = DiGraphMap::new();
        let edges = [
            (0, 1, 0.90),
            (0, 2, 0.85),
            (1, 3, 0.80),
            (2, 3, 0.70),
            (2, 4, 0.75),
            (3, 5, 0.95),
            (4, 5, 0.65),
            (4, 6, 0.60),
        ];
        for &(u, v, w) in &edges {
            graph.add_edge(u, v, w);
        }

        let start_nodes = vec![4, 0, 2, 1];
        let end_nodes = vec![5, 6];

        let baseline = ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .expect("single-thread pool should build")
            .install(|| {
                summarize_paths(parallel_path_discovery(&graph, &start_nodes, &end_nodes, 8))
            });

        let multithread = ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("multi-thread pool should build")
            .install(|| {
                summarize_paths(parallel_path_discovery(&graph, &start_nodes, &end_nodes, 8))
            });

        assert_eq!(multithread, baseline);

        let mut rng = StdRng::seed_from_u64(0xDE7E_0123);
        for _ in 0..64 {
            let mut shuffled = start_nodes.clone();
            shuffled.shuffle(&mut rng);

            let observed = ThreadPoolBuilder::new()
                .num_threads(3)
                .build()
                .expect("thread pool should build")
                .install(|| {
                    summarize_paths(parallel_path_discovery(&graph, &shuffled, &end_nodes, 8))
                });
            assert_eq!(observed, baseline);
        }
    }

    #[test]
    fn parallel_path_discovery_deduplicates_repeated_start_nodes() {
        let mut graph = DiGraphMap::new();
        graph.add_edge(0, 1, 0.9);
        graph.add_edge(1, 2, 0.8);
        graph.add_edge(2, 3, 0.7);

        let unique = summarize_paths(parallel_path_discovery(&graph, &[0], &[3], 8));
        let repeated = summarize_paths(parallel_path_discovery(&graph, &[0, 0, 0, 0], &[3], 8));

        assert_eq!(repeated, unique);
    }
}
