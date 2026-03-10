use crate::graph::assembler::Contig;
use std::cmp::Ordering;

/// Collapse contigs with identical or highly similar RLE-encoded sequences
pub fn collapse_repeats(contigs: Vec<Contig>, min_repeat_len: usize) -> Vec<Contig> {
    if contigs.len() <= 1 {
        return contigs;
    }

    let signatures: Vec<Vec<(u8, u8)>> = contigs
        .iter()
        .map(|c| normalized_rle_signature(&c.sequence))
        .collect();
    let eligible: Vec<bool> = signatures
        .iter()
        .map(|sig| min_repeat_len == 0 || sig.len() >= min_repeat_len)
        .collect();

    let mut groups = DisjointSet::new(contigs.len());
    for i in 0..contigs.len() {
        if !eligible[i] {
            continue;
        }
        for j in (i + 1)..contigs.len() {
            if !eligible[j] {
                continue;
            }

            let same_or_similar = contigs[i]
                .sequence
                .eq_ignore_ascii_case(&contigs[j].sequence)
                || is_rle_similar_encoded(&signatures[i], &signatures[j], 0.9);

            if same_or_similar {
                groups.union(i, j);
            }
        }
    }

    let mut members_by_root: Vec<Vec<usize>> = vec![Vec::new(); contigs.len()];
    for idx in 0..contigs.len() {
        let root = groups.find(idx);
        members_by_root[root].push(idx);
    }

    let mut collapsed = Vec::new();
    for members in members_by_root.into_iter().filter(|m| !m.is_empty()) {
        let best_idx = members
            .into_iter()
            .min_by(|&a, &b| contig_preference(&contigs[a], &contigs[b]))
            .expect("non-empty component must have at least one member");
        collapsed.push(contigs[best_idx].clone());
    }

    // Keep output deterministic even if caller provides contigs in different orders.
    collapsed.sort_unstable_by(contig_preference);
    // Downstream overlap/isoform paths assume contig IDs are dense and index-like.
    for (new_id, contig) in collapsed.iter_mut().enumerate() {
        contig.id = new_id;
    }
    collapsed
}

/// Determine if two sequences are similar based on their run-length encoding
/// Compares the base types and their frequencies
fn is_rle_similar(seq1: &str, seq2: &str, threshold: f64) -> bool {
    let rle1 = normalized_rle_signature(seq1);
    let rle2 = normalized_rle_signature(seq2);
    is_rle_similar_encoded(&rle1, &rle2, threshold)
}

#[inline]
fn normalized_rle_signature(seq: &str) -> Vec<(u8, u8)> {
    let bytes = seq.as_bytes();
    if bytes.is_empty() {
        return Vec::new();
    }

    let mut encoded = Vec::new();
    let mut current = bytes[0].to_ascii_uppercase();
    let mut count: u8 = 1;

    for &next_raw in &bytes[1..] {
        let next = next_raw.to_ascii_uppercase();
        if next == current {
            if count == u8::MAX {
                encoded.push((current, u8::MAX));
                count = 1;
            } else {
                count += 1;
            }
        } else {
            encoded.push((current, count));
            current = next;
            count = 1;
        }
    }

    encoded.push((current, count));
    encoded
}

/// Determine if two pre-encoded RLE signatures are similar.
fn is_rle_similar_encoded(rle1: &[(u8, u8)], rle2: &[(u8, u8)], threshold: f64) -> bool {
    #[inline]
    fn rel_diff(lhs: usize, rhs: usize) -> f64 {
        let denom = lhs.max(rhs).max(1) as f64;
        lhs.abs_diff(rhs) as f64 / denom
    }

    // If lengths are too different, consider them dissimilar
    if rle1.is_empty() || rle2.is_empty() || rel_diff(rle1.len(), rle2.len()) > 0.3 {
        return false;
    }

    // Compare the base types and their frequencies
    let mut matches = 0;
    let total = rle1.len().max(rle2.len());

    for i in 0..rle1.len().min(rle2.len()) {
        if rle1[i].0 == rle2[i].0 && rel_diff(rle1[i].1 as usize, rle2[i].1 as usize) <= 0.2 {
            matches += 1;
        }
    }

    (matches as f64 / total as f64) >= threshold
}

fn contig_preference(a: &Contig, b: &Contig) -> Ordering {
    b.sequence
        .len()
        .cmp(&a.sequence.len())
        .then_with(|| a.sequence.cmp(&b.sequence))
        .then_with(|| a.kmer_path.cmp(&b.kmer_path))
        .then_with(|| a.id.cmp(&b.id))
}

#[derive(Debug, Clone)]
struct DisjointSet {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl DisjointSet {
    fn new(size: usize) -> Self {
        Self {
            parent: (0..size).collect(),
            rank: vec![0; size],
        }
    }

    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            let root = self.find(self.parent[x]);
            self.parent[x] = root;
        }
        self.parent[x]
    }

    fn union(&mut self, a: usize, b: usize) {
        let mut root_a = self.find(a);
        let mut root_b = self.find(b);
        if root_a == root_b {
            return;
        }

        if self.rank[root_a] < self.rank[root_b] {
            std::mem::swap(&mut root_a, &mut root_b);
        }
        self.parent[root_b] = root_a;
        if self.rank[root_a] == self.rank[root_b] {
            self.rank[root_a] = self.rank[root_a].saturating_add(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::kmer::encode_kmer;

    #[test]
    fn test_is_rle_similar() {
        // The RLE of "ATATATATATAT" might be [(A,1), (T,1)] * 6
        // The RLE of "ATATATATATAT" should be similar to itself
        assert!(is_rle_similar("ATATATATATAT", "ATATATATATAT", 0.9));

        // Should be similar with small changes
        assert!(is_rle_similar("ATATATATATAT", "ATATATATCTAT", 0.8));

        // Should not be similar with large changes
        assert!(!is_rle_similar("ATATATATATAT", "GCGCGCGCGCGC", 0.5));
    }

    #[test]
    fn test_is_rle_similar_is_case_insensitive() {
        assert!(is_rle_similar("aaaattttcccc", "AAAATTTTCCCC", 1.0));
    }

    #[test]
    fn test_collapse_repeats() {
        let contig1 = Contig {
            id: 1,
            sequence: "ATATATATATAT".to_string(),
            kmer_path: vec![encode_kmer("ATG").unwrap(), encode_kmer("TGC").unwrap()],
        };

        let contig2 = Contig {
            id: 2,
            sequence: "ATATATATATAT".to_string(),
            kmer_path: vec![encode_kmer("GCA").unwrap(), encode_kmer("CAT").unwrap()],
        };

        let contig3 = Contig {
            id: 3,
            sequence: "GCGCGCGCGCGC".to_string(),
            kmer_path: vec![encode_kmer("GCG").unwrap(), encode_kmer("CGC").unwrap()],
        };

        let contigs = vec![contig1, contig2, contig3];

        // Collapse with high similarity threshold (0.9)
        let collapsed = collapse_repeats(contigs.clone(), 0);

        // Should merge the first two contigs
        assert_eq!(collapsed.len(), 2);

        // Verify that one of the collapsed contigs is the similarity one
        let contains_similar = collapsed.iter().any(|c| c.sequence == "ATATATATATAT");
        assert!(contains_similar);

        // Verify that the non-similar contig is still there
        let contains_different = collapsed.iter().any(|c| c.sequence == "GCGCGCGCGCGC");
        assert!(contains_different);
    }

    #[test]
    fn test_collapse_repeats_respects_min_repeat_len_threshold() {
        let contig1 = Contig {
            id: 1,
            sequence: "ATATAT".to_string(),
            kmer_path: vec![encode_kmer("ATA").unwrap()],
        };
        let contig2 = Contig {
            id: 2,
            sequence: "ATATAT".to_string(),
            kmer_path: vec![encode_kmer("TAT").unwrap()],
        };

        // RLE tuple length is 6, so with threshold above 6 no collapse should happen.
        let no_collapse = collapse_repeats(vec![contig1.clone(), contig2.clone()], 7);
        assert_eq!(no_collapse.len(), 2);

        let collapsed = collapse_repeats(vec![contig1, contig2], 6);
        assert_eq!(collapsed.len(), 1);
    }

    #[test]
    fn test_collapse_repeats_is_invariant_to_input_order() {
        let contig_a = Contig {
            id: 10,
            sequence: "ATATATATATAT".to_string(),
            kmer_path: vec![encode_kmer("ATA").unwrap()],
        };
        let contig_b = Contig {
            id: 11,
            sequence: "ATATATATCTAT".to_string(),
            kmer_path: vec![encode_kmer("TAT").unwrap()],
        };
        let contig_c = Contig {
            id: 12,
            sequence: "GCGCGCGCGCGC".to_string(),
            kmer_path: vec![encode_kmer("GCG").unwrap()],
        };

        let first = collapse_repeats(
            vec![contig_a.clone(), contig_b.clone(), contig_c.clone()],
            0,
        );
        let second = collapse_repeats(vec![contig_c, contig_b, contig_a], 0);

        assert_eq!(first.len(), second.len());
        let first_sequences: Vec<&str> = first.iter().map(|c| c.sequence.as_str()).collect();
        let second_sequences: Vec<&str> = second.iter().map(|c| c.sequence.as_str()).collect();
        assert_eq!(first_sequences, second_sequences);
        let first_ids: Vec<usize> = first.iter().map(|c| c.id).collect();
        let second_ids: Vec<usize> = second.iter().map(|c| c.id).collect();
        assert_eq!(first_ids, second_ids);
    }

    #[test]
    fn is_rle_similarity_is_symmetric_for_near_threshold_run_lengths() {
        let shorter = "AAAACCCCGGGGTTTTAAAA";
        let longer = "AAAAACCCCCGGGGGTTTTTAAAAA";

        let left_to_right = is_rle_similar(shorter, longer, 0.9);
        let right_to_left = is_rle_similar(longer, shorter, 0.9);

        assert_eq!(left_to_right, right_to_left);
        assert!(left_to_right);
    }

    #[test]
    fn collapse_repeats_is_order_invariant_for_run_length_deltas() {
        let contig_short = Contig {
            id: 1,
            sequence: "AAAACCCCGGGGTTTTAAAA".to_string(),
            kmer_path: vec![encode_kmer("AAA").unwrap()],
        };
        let contig_long = Contig {
            id: 2,
            sequence: "AAAAACCCCCGGGGGTTTTTAAAAA".to_string(),
            kmer_path: vec![encode_kmer("AAC").unwrap()],
        };

        let first = collapse_repeats(vec![contig_short.clone(), contig_long.clone()], 0);
        let second = collapse_repeats(vec![contig_long, contig_short], 0);

        let first_fingerprint: Vec<(usize, &str, &[u64])> = first
            .iter()
            .map(|c| (c.id, c.sequence.as_str(), c.kmer_path.as_slice()))
            .collect();
        let second_fingerprint: Vec<(usize, &str, &[u64])> = second
            .iter()
            .map(|c| (c.id, c.sequence.as_str(), c.kmer_path.as_slice()))
            .collect();
        assert_eq!(first_fingerprint, second_fingerprint);
        assert_eq!(first.len(), 1);
    }

    #[test]
    fn collapse_repeats_reindexes_ids_to_dense_order() {
        let contigs = vec![
            Contig {
                id: 41,
                sequence: "ACGTACGT".to_string(),
                kmer_path: vec![encode_kmer("ACG").unwrap()],
            },
            Contig {
                id: 99,
                sequence: "ACGTACGT".to_string(),
                kmer_path: vec![encode_kmer("CGT").unwrap()],
            },
            Contig {
                id: 7,
                sequence: "TTTTCCCC".to_string(),
                kmer_path: vec![encode_kmer("TTT").unwrap()],
            },
        ];

        let collapsed = collapse_repeats(contigs, 0);
        let ids: Vec<usize> = collapsed.iter().map(|c| c.id).collect();
        assert_eq!(ids, vec![0, 1]);
    }

    #[test]
    fn collapse_repeats_deduplicates_case_only_sequence_variants() {
        let contigs = vec![
            Contig {
                id: 5,
                sequence: "aaaattttcccc".to_string(),
                kmer_path: vec![encode_kmer("AAA").unwrap()],
            },
            Contig {
                id: 6,
                sequence: "AAAATTTTCCCC".to_string(),
                kmer_path: vec![encode_kmer("AAC").unwrap()],
            },
        ];

        let collapsed = collapse_repeats(contigs, 0);
        assert_eq!(collapsed.len(), 1);
        assert_eq!(collapsed[0].id, 0);
    }
}
