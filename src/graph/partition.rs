use crate::graph::assembler::Contig;
use serde::Serialize;
use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct RaptorComponent {
    pub id: usize,
    pub contig_ids: Vec<usize>,
    pub contig_count: usize,
    pub total_bases: usize,
}

pub fn merge_clusters(
    clusters: Vec<Vec<usize>>,
    _distance_matrix: &[Vec<f32>],
    _contigs: &[String],
) -> Vec<Vec<usize>> {
    clusters
}

pub fn cluster_contigs_by_shared_sequence(
    contigs: &[Contig],
    min_shared_bases: usize,
) -> Vec<RaptorComponent> {
    if contigs.is_empty() {
        return Vec::new();
    }

    let threshold = min_shared_bases.max(1);
    let mut parent: Vec<usize> = (0..contigs.len()).collect();

    for left_idx in 0..contigs.len() {
        for right_idx in left_idx + 1..contigs.len() {
            if longest_common_substring_len(
                contigs[left_idx].sequence.as_bytes(),
                contigs[right_idx].sequence.as_bytes(),
            ) >= threshold
            {
                union(&mut parent, left_idx, right_idx);
            }
        }
    }

    let mut grouped: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for idx in 0..contigs.len() {
        let root = find(&mut parent, idx);
        grouped.entry(root).or_default().push(idx);
    }

    grouped
        .into_values()
        .enumerate()
        .map(|(component_id, member_indexes)| {
            let mut contig_ids: Vec<usize> = member_indexes
                .into_iter()
                .map(|idx| contigs[idx].id)
                .collect();
            contig_ids.sort_unstable();
            let total_bases = contig_ids
                .iter()
                .filter_map(|id| contigs.iter().find(|contig| contig.id == *id))
                .map(|contig| contig.sequence.len())
                .sum();
            RaptorComponent {
                id: component_id,
                contig_count: contig_ids.len(),
                contig_ids,
                total_bases,
            }
        })
        .collect()
}

fn find(parent: &mut [usize], idx: usize) -> usize {
    if parent[idx] != idx {
        parent[idx] = find(parent, parent[idx]);
    }
    parent[idx]
}

fn union(parent: &mut [usize], left: usize, right: usize) {
    let left_root = find(parent, left);
    let right_root = find(parent, right);
    if left_root == right_root {
        return;
    }
    if left_root < right_root {
        parent[right_root] = left_root;
    } else {
        parent[left_root] = right_root;
    }
}

fn longest_common_substring_len(left: &[u8], right: &[u8]) -> usize {
    if left.is_empty() || right.is_empty() {
        return 0;
    }

    let mut previous = vec![0usize; right.len() + 1];
    let mut best = 0usize;
    for &left_base in left {
        let mut current = vec![0usize; right.len() + 1];
        for (idx, &right_base) in right.iter().enumerate() {
            if left_base.eq_ignore_ascii_case(&right_base) {
                let value = previous[idx] + 1;
                current[idx + 1] = value;
                best = best.max(value);
            }
        }
        previous = current;
    }
    best
}

#[cfg(test)]
mod tests {
    use super::cluster_contigs_by_shared_sequence;
    use crate::graph::assembler::Contig;

    fn contig(id: usize, sequence: &str) -> Contig {
        Contig {
            id,
            sequence: sequence.to_string(),
            kmer_path: Vec::new(),
        }
    }

    #[test]
    fn shared_sequence_clusters_related_transcripts() {
        let contigs = vec![
            contig(10, "AAAACCCCGGGGTTTT"),
            contig(20, "CCCCGGGGAAAATTTT"),
            contig(30, "TATATATATATA"),
        ];

        let components = cluster_contigs_by_shared_sequence(&contigs, 8);
        assert_eq!(components.len(), 2);
        assert_eq!(components[0].contig_ids, vec![10, 20]);
        assert_eq!(components[0].contig_count, 2);
        assert_eq!(components[1].contig_ids, vec![30]);
    }

    #[test]
    fn component_output_is_stable_under_disconnected_inputs() {
        let contigs = vec![contig(2, "AAAA"), contig(0, "CCCC"), contig(1, "GGGG")];
        let components = cluster_contigs_by_shared_sequence(&contigs, 5);
        let ids: Vec<Vec<usize>> = components
            .iter()
            .map(|component| component.contig_ids.clone())
            .collect();
        assert_eq!(ids, vec![vec![2], vec![0], vec![1]]);
    }
}
