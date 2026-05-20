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
    pub read_kmer_k: usize,
    pub read_kmer_node_count: usize,
    pub read_kmer_edge_count: usize,
    pub read_kmer_nodes: Vec<String>,
    pub read_kmer_edges: Vec<RaptorReadKmerEdge>,
    pub read_kmer_edges_sample: Vec<RaptorReadKmerEdge>,
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
    pub shared_read_count: usize,
    pub shared_pair_count: usize,
    pub shared_observed_kmer_count: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct RaptorReadKmerEdge {
    pub from: String,
    pub to: String,
    pub support: usize,
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
    cluster_contigs(contigs, |left, right| {
        longest_common_substring_len(left.sequence.as_bytes(), right.sequence.as_bytes())
            >= min_shared_bases.max(1)
    })
}

pub fn cluster_contigs_by_sequence_or_read_kmers(
    contigs: &[Contig],
    min_shared_bases: usize,
    reads1: &[String],
    reads2: Option<&[String]>,
    read_k: usize,
    min_shared_observed_kmers: usize,
) -> Vec<RaptorComponent> {
    let threshold = min_shared_bases.max(1);
    let kmer_threshold = min_shared_observed_kmers.max(1);
    cluster_contigs(contigs, |left, right| {
        longest_common_substring_len(left.sequence.as_bytes(), right.sequence.as_bytes())
            >= threshold
            || shared_observed_kmer_count(left, right, reads1, reads2, read_k) >= kmer_threshold
    })
}

fn cluster_contigs<F>(contigs: &[Contig], mut should_link: F) -> Vec<RaptorComponent>
where
    F: FnMut(&Contig, &Contig) -> bool,
{
    if contigs.is_empty() {
        return Vec::new();
    }

    let mut parent: Vec<usize> = (0..contigs.len()).collect();

    for left_idx in 0..contigs.len() {
        for right_idx in left_idx + 1..contigs.len() {
            if should_link(&contigs[left_idx], &contigs[right_idx]) {
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
    reads1: &[String],
    reads2: Option<&[String]>,
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
                            shared_read_count: shared_read_support(left, right, reads1, reads2),
                            shared_pair_count: shared_pair_support(left, right, reads1, reads2),
                            shared_observed_kmer_count: shared_observed_kmer_count(
                                left,
                                right,
                                reads1,
                                reads2,
                                25.min(shared_bases),
                            ),
                        });
                    }
                }
            }
            let read_kmer_k = component_read_kmer_k(contigs, component, reads1, reads2, 25);
            let read_kmer_nodes =
                component_read_kmer_nodes(contigs, component, reads1, reads2, read_kmer_k);
            let read_kmer_edges =
                component_read_kmer_edges(contigs, component, reads1, reads2, read_kmer_k);
            let read_kmer_edges: Vec<RaptorReadKmerEdge> = read_kmer_edges
                .iter()
                .map(|((from, to), support)| RaptorReadKmerEdge {
                    from: from.clone(),
                    to: to.clone(),
                    support: *support,
                })
                .collect();
            let read_kmer_edges_sample = read_kmer_edges.iter().take(32).cloned().collect();
            RaptorComponentGraph {
                component_id: component.id,
                node_count: nodes.len(),
                edge_count: edges.len(),
                read_kmer_k,
                read_kmer_node_count: read_kmer_nodes.len(),
                read_kmer_edge_count: read_kmer_edges.len(),
                read_kmer_nodes: read_kmer_nodes.into_iter().collect(),
                read_kmer_edges,
                read_kmer_edges_sample,
                nodes,
                edges,
                assigned_read_count: component.assigned_read_count,
                assigned_pair_count: component.assigned_pair_count,
            }
        })
        .collect()
}

fn component_read_kmer_k(
    contigs: &[Contig],
    component: &RaptorComponent,
    reads1: &[String],
    reads2: Option<&[String]>,
    preferred_k: usize,
) -> usize {
    assigned_component_reads(contigs, component, reads1, reads2)
        .iter()
        .map(|read| read.len())
        .min()
        .map(|min_len| preferred_k.min(min_len).max(1))
        .unwrap_or(preferred_k.max(1))
}

fn component_read_kmer_edges(
    contigs: &[Contig],
    component: &RaptorComponent,
    reads1: &[String],
    reads2: Option<&[String]>,
    k: usize,
) -> BTreeMap<(String, String), usize> {
    let assigned_reads = assigned_component_reads(contigs, component, reads1, reads2);
    let mut edges = BTreeMap::new();
    for read in &assigned_reads {
        collect_read_kmer_edges(read, k, &mut edges);
        let read_rc = reverse_complement(read);
        collect_read_kmer_edges(&read_rc, k, &mut edges);
    }
    edges
}

fn component_read_kmer_nodes(
    contigs: &[Contig],
    component: &RaptorComponent,
    reads1: &[String],
    reads2: Option<&[String]>,
    k: usize,
) -> BTreeSet<String> {
    let assigned_reads = assigned_component_reads(contigs, component, reads1, reads2);
    let mut nodes = BTreeSet::new();
    for read in &assigned_reads {
        collect_observed_kmers(std::slice::from_ref(read), k, &mut nodes);
        let read_rc = reverse_complement(read);
        collect_observed_kmers(std::slice::from_ref(&read_rc), k, &mut nodes);
    }
    nodes
}

fn assigned_component_reads(
    contigs: &[Contig],
    component: &RaptorComponent,
    reads1: &[String],
    reads2: Option<&[String]>,
) -> Vec<String> {
    let component_contigs: Vec<&Contig> = component
        .contig_ids
        .iter()
        .filter_map(|contig_id| contigs.iter().find(|contig| contig.id == *contig_id))
        .collect();
    let mut reads = Vec::new();
    for read in reads1 {
        if component_contigs
            .iter()
            .any(|contig| read_matches_sequence(read, &contig.sequence))
        {
            reads.push(read.clone());
        }
    }
    if let Some(reads2) = reads2 {
        for read in reads2 {
            if component_contigs
                .iter()
                .any(|contig| read_matches_sequence(read, &contig.sequence))
            {
                reads.push(read.clone());
            }
        }
    }
    reads
}

fn collect_read_kmer_edges(read: &str, k: usize, edges: &mut BTreeMap<(String, String), usize>) {
    if read.len() <= k {
        return;
    }
    let kmers: Vec<String> = read
        .as_bytes()
        .windows(k)
        .map(|window| String::from_utf8_lossy(window).to_ascii_uppercase())
        .collect();
    for pair in kmers.windows(2) {
        *edges.entry((pair[0].clone(), pair[1].clone())).or_default() += 1;
    }
}

fn shared_read_support(
    left: &Contig,
    right: &Contig,
    reads1: &[String],
    reads2: Option<&[String]>,
) -> usize {
    let left = left.sequence.as_str();
    let right = right.sequence.as_str();
    let mut count = reads1
        .iter()
        .filter(|read| read_matches_both(read, left, right))
        .count();
    if let Some(reads2) = reads2 {
        count += reads2
            .iter()
            .filter(|read| read_matches_both(read, left, right))
            .count();
    }
    count
}

fn shared_pair_support(
    left: &Contig,
    right: &Contig,
    reads1: &[String],
    reads2: Option<&[String]>,
) -> usize {
    let Some(reads2) = reads2 else {
        return 0;
    };
    let left = left.sequence.as_str();
    let right = right.sequence.as_str();
    reads1
        .iter()
        .zip(reads2.iter())
        .filter(|(read1, read2)| {
            let pair_hits_left =
                read_matches_sequence(read1, left) || read_matches_sequence(read2, left);
            let pair_hits_right =
                read_matches_sequence(read1, right) || read_matches_sequence(read2, right);
            pair_hits_left && pair_hits_right
        })
        .count()
}

fn shared_observed_kmer_count(
    left: &Contig,
    right: &Contig,
    reads1: &[String],
    reads2: Option<&[String]>,
    k: usize,
) -> usize {
    let k = k.min(left.sequence.len()).min(right.sequence.len()).max(1);
    let mut observed = BTreeSet::new();
    collect_observed_kmers(reads1, k, &mut observed);
    if let Some(reads2) = reads2 {
        collect_observed_kmers(reads2, k, &mut observed);
    }
    observed
        .into_iter()
        .filter(|kmer| {
            sequence_contains(&left.sequence, kmer) && sequence_contains(&right.sequence, kmer)
        })
        .count()
}

fn collect_observed_kmers(reads: &[String], k: usize, observed: &mut BTreeSet<String>) {
    for read in reads {
        if read.len() < k {
            continue;
        }
        for window in read.as_bytes().windows(k) {
            observed.insert(String::from_utf8_lossy(window).to_ascii_uppercase());
        }
        let read_rc = reverse_complement(read);
        for window in read_rc.as_bytes().windows(k) {
            observed.insert(String::from_utf8_lossy(window).to_ascii_uppercase());
        }
    }
}

fn read_matches_both(read: &str, left: &str, right: &str) -> bool {
    read_matches_sequence(read, left) && read_matches_sequence(read, right)
}

fn read_matches_sequence(read: &str, sequence: &str) -> bool {
    let read_rc = reverse_complement(read);
    sequence_contains(sequence, read)
        || (!read_rc.eq_ignore_ascii_case(read) && sequence_contains(sequence, &read_rc))
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
        cluster_contigs_by_sequence_or_read_kmers, cluster_contigs_by_shared_sequence,
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
    fn read_kmer_connectivity_can_define_component_boundary() {
        let contigs = vec![
            contig(10, "AAAACCCCGGGG"),
            contig(20, "CCCCAAAATTTT"),
            contig(30, "TATATATATATA"),
        ];
        let reads1 = vec!["CCCC".to_string()];

        let sequence_only = cluster_contigs_by_shared_sequence(&contigs, 8);
        assert_eq!(sequence_only.len(), 3);

        let read_linked =
            cluster_contigs_by_sequence_or_read_kmers(&contigs, 8, &reads1, None, 4, 1);
        assert_eq!(read_linked.len(), 2);
        assert_eq!(read_linked[0].contig_ids, vec![10, 20]);
        assert_eq!(read_linked[1].contig_ids, vec![30]);
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
        let reads1 = vec!["CCCCGGGG".to_string()];
        let reads2 = vec!["CCCCGGGG".to_string()];
        let graphs = build_component_graphs(&contigs, &components, 8, &reads1, Some(&reads2));

        assert_eq!(graphs.len(), 2);
        assert_eq!(graphs[0].component_id, 0);
        assert_eq!(graphs[0].node_count, 2);
        assert_eq!(graphs[0].edge_count, 1);
        assert_eq!(graphs[0].read_kmer_k, 8);
        assert_eq!(graphs[0].read_kmer_node_count, 1);
        assert_eq!(graphs[0].read_kmer_edge_count, 0);
        assert_eq!(graphs[0].read_kmer_nodes, vec!["CCCCGGGG".to_string()]);
        assert!(graphs[0].read_kmer_edges.is_empty());
        assert!(graphs[0].read_kmer_edges_sample.is_empty());
        assert_eq!(graphs[0].edges[0].left_contig_id, 10);
        assert_eq!(graphs[0].edges[0].right_contig_id, 20);
        assert_eq!(graphs[0].edges[0].shared_bases, 8);
        assert_eq!(graphs[0].edges[0].shared_read_count, 2);
        assert_eq!(graphs[0].edges[0].shared_pair_count, 1);
        assert_eq!(graphs[0].edges[0].shared_observed_kmer_count, 1);
        assert_eq!(graphs[1].node_count, 1);
        assert_eq!(graphs[1].edge_count, 0);
    }

    #[test]
    fn component_graphs_include_read_kmer_edge_sample() {
        let contigs = vec![contig(0, "AAAACCCCGGGGTTTT")];
        let components = cluster_contigs_by_shared_sequence(&contigs, 8);
        let reads1 = vec!["AAAACCCCG".to_string(), "AAAACCCCGG".to_string()];

        let graphs = build_component_graphs(&contigs, &components, 8, &reads1, None);

        assert_eq!(graphs[0].read_kmer_k, 9);
        assert!(graphs[0].read_kmer_node_count >= 2);
        assert!(graphs[0].read_kmer_edge_count >= 1);
        assert_eq!(
            graphs[0].read_kmer_nodes.len(),
            graphs[0].read_kmer_node_count
        );
        assert_eq!(
            graphs[0].read_kmer_edges.len(),
            graphs[0].read_kmer_edge_count
        );
        assert!(!graphs[0].read_kmer_edges_sample.is_empty());
        assert!(graphs[0]
            .read_kmer_edges_sample
            .iter()
            .all(|edge| edge.support >= 1));
        assert_eq!(
            graphs[0].read_kmer_edges_sample,
            graphs[0]
                .read_kmer_edges
                .iter()
                .take(graphs[0].read_kmer_edges_sample.len())
                .cloned()
                .collect::<Vec<_>>()
        );
    }
}
