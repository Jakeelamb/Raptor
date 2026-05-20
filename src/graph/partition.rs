use crate::graph::assembler::Contig;
use crate::kmer::kmer::reverse_complement;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct RaptorComponent {
    pub id: usize,
    pub contig_ids: Vec<usize>,
    pub contig_count: usize,
    pub total_bases: usize,
    pub assigned_read_count: usize,
    pub assigned_pair_count: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct RaptorComponentGraph {
    pub component_id: usize,
    pub nodes: Vec<RaptorComponentNode>,
    pub edges: Vec<RaptorComponentEdge>,
    pub node_count: usize,
    pub edge_count: usize,
    pub assigned_read_count: usize,
    pub assigned_pair_count: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct RaptorComponentNode {
    pub contig_id: usize,
    pub length: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct RaptorComponentEdge {
    pub left_contig_id: usize,
    pub right_contig_id: usize,
    pub shared_bases: usize,
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
                assigned_read_count: 0,
                assigned_pair_count: 0,
            }
        })
        .collect()
}

pub fn assign_read_evidence_to_components(
    contigs: &[Contig],
    components: &mut [RaptorComponent],
    reads1: &[String],
    reads2: Option<&[String]>,
) {
    let contig_to_component = contig_component_index(components);
    for read in reads1 {
        for component_idx in matching_components(contigs, &contig_to_component, read) {
            components[component_idx].assigned_read_count += 1;
        }
    }

    let Some(reads2) = reads2 else {
        return;
    };

    for (read1, read2) in reads1.iter().zip(reads2.iter()) {
        let mut pair_components = matching_components(contigs, &contig_to_component, read1);
        for component_idx in matching_components(contigs, &contig_to_component, read2) {
            pair_components.insert(component_idx);
            components[component_idx].assigned_read_count += 1;
        }
        for component_idx in pair_components {
            components[component_idx].assigned_pair_count += 1;
        }
    }
}

pub fn build_component_graphs(
    contigs: &[Contig],
    components: &[RaptorComponent],
    min_shared_bases: usize,
) -> Vec<RaptorComponentGraph> {
    let contigs_by_id: BTreeMap<usize, &Contig> =
        contigs.iter().map(|contig| (contig.id, contig)).collect();
    let threshold = min_shared_bases.max(1);

    components
        .iter()
        .map(|component| {
            let nodes: Vec<RaptorComponentNode> = component
                .contig_ids
                .iter()
                .filter_map(|contig_id| {
                    contigs_by_id
                        .get(contig_id)
                        .map(|contig| RaptorComponentNode {
                            contig_id: *contig_id,
                            length: contig.sequence.len(),
                        })
                })
                .collect();
            let mut edges = Vec::new();
            for left_idx in 0..component.contig_ids.len() {
                for right_idx in left_idx + 1..component.contig_ids.len() {
                    let left_id = component.contig_ids[left_idx];
                    let right_id = component.contig_ids[right_idx];
                    let Some(left) = contigs_by_id.get(&left_id) else {
                        continue;
                    };
                    let Some(right) = contigs_by_id.get(&right_id) else {
                        continue;
                    };
                    let shared_bases = longest_common_substring_len(
                        left.sequence.as_bytes(),
                        right.sequence.as_bytes(),
                    );
                    if shared_bases >= threshold {
                        edges.push(RaptorComponentEdge {
                            left_contig_id: left_id,
                            right_contig_id: right_id,
                            shared_bases,
                        });
                    }
                }
            }
            RaptorComponentGraph {
                component_id: component.id,
                node_count: nodes.len(),
                edge_count: edges.len(),
                nodes,
                edges,
                assigned_read_count: component.assigned_read_count,
                assigned_pair_count: component.assigned_pair_count,
            }
        })
        .collect()
}

fn contig_component_index(components: &[RaptorComponent]) -> BTreeMap<usize, usize> {
    let mut index = BTreeMap::new();
    for (component_idx, component) in components.iter().enumerate() {
        for &contig_id in &component.contig_ids {
            index.insert(contig_id, component_idx);
        }
    }
    index
}

fn matching_components(
    contigs: &[Contig],
    contig_to_component: &BTreeMap<usize, usize>,
    read: &str,
) -> BTreeSet<usize> {
    let read_rc = reverse_complement(read);
    let mut matches = BTreeSet::new();
    for contig in contigs {
        if sequence_contains(&contig.sequence, read)
            || (!read_rc.eq_ignore_ascii_case(read)
                && sequence_contains(&contig.sequence, &read_rc))
        {
            if let Some(&component_idx) = contig_to_component.get(&contig.id) {
                matches.insert(component_idx);
            }
        }
    }
    matches
}

fn sequence_contains(haystack: &str, needle: &str) -> bool {
    let haystack = haystack.as_bytes();
    let needle = needle.as_bytes();
    if needle.is_empty() || needle.len() > haystack.len() {
        return false;
    }
    haystack
        .windows(needle.len())
        .any(|window| window.eq_ignore_ascii_case(needle))
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
    use super::{
        assign_read_evidence_to_components, build_component_graphs,
        cluster_contigs_by_shared_sequence,
    };
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

    #[test]
    fn read_evidence_counts_component_support_once_per_pair() {
        let contigs = vec![
            contig(0, "AAAACCCCGGGGTTTT"),
            contig(1, "CCCCGGGGAAAATTTT"),
            contig(2, "TATATATATATA"),
        ];
        let mut components = cluster_contigs_by_shared_sequence(&contigs, 8);

        let reads1 = vec!["AAAACCCC".to_string(), "TATATATA".to_string()];
        let reads2 = vec!["AAAACCCC".to_string(), "TATATATA".to_string()];
        assign_read_evidence_to_components(&contigs, &mut components, &reads1, Some(&reads2));

        assert_eq!(components[0].contig_ids, vec![0, 1]);
        assert_eq!(components[0].assigned_read_count, 2);
        assert_eq!(components[0].assigned_pair_count, 1);
        assert_eq!(components[1].contig_ids, vec![2]);
        assert_eq!(components[1].assigned_read_count, 2);
        assert_eq!(components[1].assigned_pair_count, 1);
    }

    #[test]
    fn component_graphs_preserve_shared_sequence_edges() {
        let contigs = vec![
            contig(10, "AAAACCCCGGGGTTTT"),
            contig(20, "CCCCGGGGAAAATTTT"),
            contig(30, "TATATATATATA"),
        ];
        let components = cluster_contigs_by_shared_sequence(&contigs, 8);
        let graphs = build_component_graphs(&contigs, &components, 8);

        assert_eq!(graphs.len(), 2);
        assert_eq!(graphs[0].component_id, 0);
        assert_eq!(graphs[0].node_count, 2);
        assert_eq!(graphs[0].edge_count, 1);
        assert_eq!(graphs[0].edges[0].left_contig_id, 10);
        assert_eq!(graphs[0].edges[0].right_contig_id, 20);
        assert_eq!(graphs[0].edges[0].shared_bases, 8);
        assert_eq!(graphs[1].node_count, 1);
        assert_eq!(graphs[1].edge_count, 0);
    }
}
