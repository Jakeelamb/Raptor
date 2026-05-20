use crate::graph::assembler::Contig;
use crate::graph::partition::{
    assign_read_evidence_to_components, build_component_graphs,
    cluster_contigs_by_sequence_or_read_kmers, RaptorComponent, RaptorComponentGraph,
};
use crate::io::fasta::try_open_fasta;
use crate::io::fastq::{stream_fastq_records_checked, try_open_fastq, FastqWriter};
use crate::kmer::kmer::reverse_complement;
use crate::pipeline::assemble::assemble_reads_with_gpu;
use crate::pipeline::normalize::{
    normalize_paired_with_config, normalize_single_with_config, NormalizeConfig, NormalizeSummary,
};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::{self, BufRead};
use std::path::{Path, PathBuf};
use std::time::Instant;

const MIN_SELECTED_NOVEL_SEQUENCE_FRACTION: f64 = 0.5;
const COMPONENT_CONTIG_SELECTION_METHOD: &str = "component_contig_evidence_score_v2";

#[derive(Debug, Clone)]
pub struct TrinityWorkflowConfig {
    pub input1: Option<String>,
    pub input2: Option<String>,
    pub samples_file: Option<String>,
    pub ss_lib_type: Option<String>,
    pub output_dir: String,
    pub output_fasta: Option<String>,
    pub report_json: Option<String>,
    pub normalize: bool,
    pub normalize_config: NormalizeConfig,
    pub min_len: usize,
    pub use_gpu: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct TrinityWorkflowReport {
    pub workflow: &'static str,
    pub input1: String,
    pub input2: Option<String>,
    pub samples_file: Option<String>,
    pub ss_lib_type: Option<String>,
    pub sample_count: usize,
    pub output_dir: String,
    pub normalized: bool,
    pub normalized_input1: Option<String>,
    pub normalized_input2: Option<String>,
    pub normalization: Option<NormalizeSummary>,
    pub assembly_fasta: String,
    pub assembly_metrics_json: String,
    pub assembly_metrics_tsv: String,
    pub components_json: String,
    pub component_graphs_json: String,
    pub component_paths_fasta: String,
    pub component_selected_isoforms_fasta: String,
    pub component_selected_isoforms_json: String,
    pub component_isoform_candidates_json: String,
    pub component_clustering: &'static str,
    pub component_count: usize,
    pub component_graph_count: usize,
    pub report_json: String,
    pub elapsed_seconds: f64,
}

pub fn run_trinity_workflow(config: TrinityWorkflowConfig) -> io::Result<TrinityWorkflowReport> {
    let start = Instant::now();
    let output_dir = PathBuf::from(&config.output_dir);
    fs::create_dir_all(&output_dir)?;

    let assembly_fasta = config
        .output_fasta
        .as_ref()
        .map(PathBuf::from)
        .unwrap_or_else(|| output_dir.join("raptor_trinity.fasta.gz"));

    if let Some(parent) = assembly_fasta.parent() {
        fs::create_dir_all(parent)?;
    }

    let resolved_inputs = resolve_workflow_inputs(&config, &output_dir)?;
    let ss_lib_type = validate_ss_lib_type(config.ss_lib_type.as_deref())?;

    let (assembly_input1, assembly_input2, normalization) = if config.normalize {
        normalize_for_workflow(
            &resolved_inputs.input1,
            resolved_inputs.input2.as_deref(),
            &config,
            &output_dir,
        )?
    } else {
        (
            resolved_inputs.input1.clone(),
            resolved_inputs.input2.clone(),
            None,
        )
    };

    assemble_reads_with_gpu(
        &assembly_input1,
        assembly_input2.as_deref(),
        path_to_str(&assembly_fasta)?,
        config.min_len,
        false,
        false,
        false,
        false,
        false,
        20,
        false,
        25,
        false,
        false,
        None,
        None,
        false,
        None,
        None,
        20,
        0.9,
        false,
        false,
        None,
        0.1,
        None,
        false,
        config.use_gpu,
    )?;

    let assembly_fasta_string = path_to_str(&assembly_fasta)?.to_string();
    let component_path = output_dir.join("raptor_components.json");
    let component_graphs_path = output_dir.join("raptor_component_graphs.json");
    let component_paths_fasta = output_dir.join("raptor_component_paths.fasta");
    let component_selected_isoforms_fasta =
        output_dir.join("raptor_component_selected_isoforms.fasta");
    let component_selected_isoforms_json =
        output_dir.join("raptor_component_selected_isoforms.json");
    let component_isoform_candidates_json =
        output_dir.join("raptor_component_isoform_candidates.json");
    let component_artifacts = write_component_artifacts(
        &assembly_fasta_string,
        &component_path,
        &component_graphs_path,
        &component_paths_fasta,
        &component_selected_isoforms_fasta,
        &component_selected_isoforms_json,
        &component_isoform_candidates_json,
        60,
        &assembly_input1,
        assembly_input2.as_deref(),
    )?;
    let report_path = config
        .report_json
        .as_ref()
        .map(PathBuf::from)
        .unwrap_or_else(|| output_dir.join("raptor_trinity_report.json"));

    let report = TrinityWorkflowReport {
        workflow: "raptor_trinity_de_novo",
        input1: resolved_inputs.input1,
        input2: resolved_inputs.input2,
        samples_file: config.samples_file,
        ss_lib_type,
        sample_count: resolved_inputs.sample_count,
        output_dir: config.output_dir,
        normalized: config.normalize,
        normalized_input1: if config.normalize {
            Some(assembly_input1)
        } else {
            None
        },
        normalized_input2: if config.normalize {
            assembly_input2
        } else {
            None
        },
        normalization,
        assembly_metrics_json: sidecar_path(&assembly_fasta_string, "assembly_metrics.json"),
        assembly_metrics_tsv: sidecar_path(&assembly_fasta_string, "assembly_metrics.tsv"),
        components_json: path_to_str(&component_path)?.to_string(),
        component_graphs_json: path_to_str(&component_graphs_path)?.to_string(),
        component_paths_fasta: path_to_str(&component_paths_fasta)?.to_string(),
        component_selected_isoforms_fasta: path_to_str(&component_selected_isoforms_fasta)?
            .to_string(),
        component_selected_isoforms_json: path_to_str(&component_selected_isoforms_json)?
            .to_string(),
        component_isoform_candidates_json: path_to_str(&component_isoform_candidates_json)?
            .to_string(),
        component_clustering: component_artifacts.clustering,
        component_count: component_artifacts.components.len(),
        component_graph_count: component_artifacts.component_graphs.len(),
        assembly_fasta: assembly_fasta_string,
        report_json: path_to_str(&report_path)?.to_string(),
        elapsed_seconds: start.elapsed().as_secs_f64(),
    };

    write_report(&report_path, &report)?;
    Ok(report)
}

struct ComponentArtifacts {
    components: Vec<RaptorComponent>,
    component_graphs: Vec<RaptorComponentGraph>,
    clustering: &'static str,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct TrinitySampleInput {
    condition: String,
    replicate: String,
    left: String,
    right: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ResolvedWorkflowInputs {
    input1: String,
    input2: Option<String>,
    sample_count: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
struct SelectedComponentIsoform {
    id: String,
    component_id: usize,
    component_rank: usize,
    source_contig_id: usize,
    length: usize,
    selection_method: &'static str,
    evidence_score: usize,
    component_assigned_reads: usize,
    component_assigned_pairs: usize,
    direct_read_support: usize,
    direct_pair_support: usize,
    overlapping_read_kmer_path_count: usize,
    max_overlapping_read_kmer_path_support: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
struct ComponentIsoformCandidate {
    id: String,
    component_id: usize,
    candidate_rank: usize,
    source_kind: &'static str,
    source_contig_id: Option<usize>,
    source_path_index: Option<usize>,
    length: usize,
    selected: bool,
    selection_method: &'static str,
    evidence_score: usize,
    component_assigned_reads: usize,
    component_assigned_pairs: usize,
    direct_read_support: usize,
    direct_pair_support: usize,
    overlapping_read_kmer_path_count: usize,
    max_overlapping_read_kmer_path_support: usize,
}

fn write_component_artifacts(
    assembly_fasta: &str,
    component_path: &Path,
    component_graphs_path: &Path,
    component_paths_fasta: &Path,
    component_selected_isoforms_fasta: &Path,
    component_selected_isoforms_json: &Path,
    component_isoform_candidates_json: &Path,
    min_shared_bases: usize,
    reads1_path: &str,
    reads2_path: Option<&str>,
) -> io::Result<ComponentArtifacts> {
    let contigs = read_contigs_from_fasta(assembly_fasta)?;
    let reads1 = read_fastq_sequences(reads1_path)?;
    let reads2 = reads2_path.map(read_fastq_sequences).transpose()?;
    if let Some(reads2) = &reads2 {
        if reads1.len() != reads2.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "paired FASTQ inputs have different record counts after normalization: {} != {}",
                    reads1.len(),
                    reads2.len()
                ),
            ));
        }
    }
    let mut components = cluster_contigs_by_sequence_or_read_kmers(
        &contigs,
        min_shared_bases,
        &reads1,
        reads2.as_deref(),
        25,
        1,
    );
    assign_read_evidence_to_components(&contigs, &mut components, &reads1, reads2.as_deref());
    let component_graphs = build_component_graphs(
        &contigs,
        &components,
        min_shared_bases,
        &reads1,
        reads2.as_deref(),
    );
    write_json(component_path, &components, "component report")?;
    write_json(
        component_graphs_path,
        &component_graphs,
        "component graph report",
    )?;
    write_component_paths_fasta(component_paths_fasta, &contigs, &component_graphs)?;
    write_component_selected_isoforms_fasta(
        component_selected_isoforms_fasta,
        &contigs,
        &component_graphs,
        &reads1,
        reads2.as_deref(),
    )?;
    let selected_isoforms =
        select_component_isoforms(&contigs, &component_graphs, &reads1, reads2.as_deref());
    write_json(
        component_selected_isoforms_json,
        &selected_isoforms,
        "selected component isoform report",
    )?;
    let isoform_candidates =
        component_isoform_candidates(&contigs, &component_graphs, &reads1, reads2.as_deref());
    write_json(
        component_isoform_candidates_json,
        &isoform_candidates,
        "component isoform candidate report",
    )?;
    Ok(ComponentArtifacts {
        components,
        component_graphs,
        clustering: "sequence_or_read_kmer",
    })
}

fn write_component_paths_fasta(
    path: &Path,
    contigs: &[Contig],
    component_graphs: &[RaptorComponentGraph],
) -> io::Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let contigs_by_id: std::collections::BTreeMap<usize, &Contig> =
        contigs.iter().map(|contig| (contig.id, contig)).collect();
    let mut output = String::new();
    for graph in component_graphs {
        for node in &graph.nodes {
            if let Some(contig) = contigs_by_id.get(&node.contig_id) {
                output.push_str(&format!(
                    ">component_{}_contig_{} length={}\n",
                    graph.component_id,
                    node.contig_id,
                    contig.sequence.len()
                ));
                write_wrapped_fasta_sequence(&mut output, &contig.sequence);
            }
        }
        for (path_idx, read_path) in graph.read_kmer_paths.iter().enumerate() {
            output.push_str(&format!(
                ">component_{}_path_{} edges={} min_support={}\n",
                graph.component_id, path_idx, read_path.edge_count, read_path.min_support
            ));
            write_wrapped_fasta_sequence(&mut output, &read_path.sequence);
        }
    }
    fs::write(path, output)
}

fn write_component_selected_isoforms_fasta(
    path: &Path,
    contigs: &[Contig],
    component_graphs: &[RaptorComponentGraph],
    reads1: &[String],
    reads2: Option<&[String]>,
) -> io::Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let selected_isoforms = select_component_isoforms(contigs, component_graphs, reads1, reads2);
    let contigs_by_id: std::collections::BTreeMap<usize, &Contig> =
        contigs.iter().map(|contig| (contig.id, contig)).collect();
    let mut output = String::new();
    for isoform in selected_isoforms {
        if let Some(contig) = contigs_by_id.get(&isoform.source_contig_id) {
            output.push_str(&format!(
                ">{} component={} rank={} source=contig_{} length={} score={} method={} direct_reads={} direct_pairs={} read_kmer_paths={} max_path_support={}\n",
                isoform.id,
                isoform.component_id,
                isoform.component_rank,
                isoform.source_contig_id,
                isoform.length,
                isoform.evidence_score,
                isoform.selection_method,
                isoform.direct_read_support,
                isoform.direct_pair_support,
                isoform.overlapping_read_kmer_path_count,
                isoform.max_overlapping_read_kmer_path_support
            ));
            write_wrapped_fasta_sequence(&mut output, &contig.sequence);
        }
    }
    fs::write(path, output)
}

fn select_component_isoforms(
    contigs: &[Contig],
    component_graphs: &[RaptorComponentGraph],
    reads1: &[String],
    reads2: Option<&[String]>,
) -> Vec<SelectedComponentIsoform> {
    let contigs_by_id: BTreeMap<usize, &Contig> =
        contigs.iter().map(|contig| (contig.id, contig)).collect();
    let mut selected = Vec::new();
    for graph in component_graphs {
        let component_isoforms =
            ranked_component_contig_isoforms(graph, &contigs_by_id, reads1, reads2);
        let component_isoforms = select_non_redundant_component_isoforms(
            component_isoforms,
            &contigs_by_id,
            reads2.is_some(),
        );
        for (idx, mut isoform) in component_isoforms.into_iter().enumerate() {
            isoform.component_rank = idx + 1;
            isoform.id = format!("component_{}_isoform_{}", graph.component_id, idx);
            selected.push(isoform);
        }
    }
    selected
}

fn component_isoform_candidates(
    contigs: &[Contig],
    component_graphs: &[RaptorComponentGraph],
    reads1: &[String],
    reads2: Option<&[String]>,
) -> Vec<ComponentIsoformCandidate> {
    let contigs_by_id: BTreeMap<usize, &Contig> =
        contigs.iter().map(|contig| (contig.id, contig)).collect();
    let mut candidates = Vec::new();
    for graph in component_graphs {
        let mut component_candidates = Vec::new();
        let contig_isoforms =
            ranked_component_contig_isoforms(graph, &contigs_by_id, reads1, reads2);
        let selected_contig_ids: BTreeSet<usize> = select_non_redundant_component_isoforms(
            contig_isoforms.clone(),
            &contigs_by_id,
            reads2.is_some(),
        )
        .into_iter()
        .map(|isoform| isoform.source_contig_id)
        .collect();
        for isoform in contig_isoforms {
            component_candidates.push(ComponentIsoformCandidate {
                id: String::new(),
                component_id: graph.component_id,
                candidate_rank: 0,
                source_kind: "contig",
                source_contig_id: Some(isoform.source_contig_id),
                source_path_index: None,
                length: isoform.length,
                selected: selected_contig_ids.contains(&isoform.source_contig_id),
                selection_method: COMPONENT_CONTIG_SELECTION_METHOD,
                evidence_score: isoform.evidence_score,
                component_assigned_reads: graph.assigned_read_count,
                component_assigned_pairs: graph.assigned_pair_count,
                direct_read_support: isoform.direct_read_support,
                direct_pair_support: isoform.direct_pair_support,
                overlapping_read_kmer_path_count: isoform.overlapping_read_kmer_path_count,
                max_overlapping_read_kmer_path_support: isoform
                    .max_overlapping_read_kmer_path_support,
            });
        }
        for (path_idx, read_path) in graph.read_kmer_paths.iter().enumerate() {
            component_candidates.push(ComponentIsoformCandidate {
                id: String::new(),
                component_id: graph.component_id,
                candidate_rank: 0,
                source_kind: "read_kmer_path",
                source_contig_id: None,
                source_path_index: Some(path_idx),
                length: read_path.sequence.len(),
                selected: false,
                selection_method: COMPONENT_CONTIG_SELECTION_METHOD,
                evidence_score: read_path.min_support * 2 + read_path.edge_count,
                component_assigned_reads: graph.assigned_read_count,
                component_assigned_pairs: graph.assigned_pair_count,
                direct_read_support: 0,
                direct_pair_support: 0,
                overlapping_read_kmer_path_count: 1,
                max_overlapping_read_kmer_path_support: read_path.min_support,
            });
        }
        component_candidates.sort_by(|left, right| {
            right
                .evidence_score
                .cmp(&left.evidence_score)
                .then_with(|| right.direct_pair_support.cmp(&left.direct_pair_support))
                .then_with(|| right.direct_read_support.cmp(&left.direct_read_support))
                .then_with(|| right.length.cmp(&left.length))
                .then_with(|| left.source_kind.cmp(right.source_kind))
                .then_with(|| left.source_contig_id.cmp(&right.source_contig_id))
                .then_with(|| left.source_path_index.cmp(&right.source_path_index))
        });
        for (idx, mut candidate) in component_candidates.into_iter().enumerate() {
            candidate.candidate_rank = idx + 1;
            candidate.id = format!("component_{}_candidate_{}", graph.component_id, idx);
            candidates.push(candidate);
        }
    }
    candidates
}

fn ranked_component_contig_isoforms(
    graph: &RaptorComponentGraph,
    contigs_by_id: &BTreeMap<usize, &Contig>,
    reads1: &[String],
    reads2: Option<&[String]>,
) -> Vec<SelectedComponentIsoform> {
    let mut component_isoforms = Vec::new();
    for node in &graph.nodes {
        let Some(contig) = contigs_by_id.get(&node.contig_id) else {
            continue;
        };
        let (direct_read_support, direct_pair_support) =
            contig_read_pair_support(&contig.sequence, reads1, reads2);
        let overlapping_paths: Vec<&crate::graph::partition::RaptorReadKmerPath> = graph
            .read_kmer_paths
            .iter()
            .filter(|path| sequences_overlap(&contig.sequence, &path.sequence))
            .collect();
        let max_path_support = overlapping_paths
            .iter()
            .map(|path| path.min_support)
            .max()
            .unwrap_or(0);
        let evidence_score = selected_isoform_evidence_score(
            direct_read_support,
            direct_pair_support,
            overlapping_paths.len(),
            max_path_support,
        );
        component_isoforms.push(SelectedComponentIsoform {
            id: String::new(),
            component_id: graph.component_id,
            component_rank: 0,
            source_contig_id: node.contig_id,
            length: contig.sequence.len(),
            selection_method: COMPONENT_CONTIG_SELECTION_METHOD,
            evidence_score,
            component_assigned_reads: graph.assigned_read_count,
            component_assigned_pairs: graph.assigned_pair_count,
            direct_read_support,
            direct_pair_support,
            overlapping_read_kmer_path_count: overlapping_paths.len(),
            max_overlapping_read_kmer_path_support: max_path_support,
        });
    }
    component_isoforms.sort_by(|left, right| {
        right
            .direct_read_support
            .cmp(&left.direct_read_support)
            .then_with(|| right.direct_pair_support.cmp(&left.direct_pair_support))
            .then_with(|| right.evidence_score.cmp(&left.evidence_score))
            .then_with(|| right.length.cmp(&left.length))
            .then_with(|| left.source_contig_id.cmp(&right.source_contig_id))
    });
    component_isoforms
}

fn select_non_redundant_component_isoforms(
    ranked_isoforms: Vec<SelectedComponentIsoform>,
    contigs_by_id: &BTreeMap<usize, &Contig>,
    require_pair_support: bool,
) -> Vec<SelectedComponentIsoform> {
    let mut selected = Vec::new();
    let mut selected_sequences: Vec<&str> = Vec::new();
    for isoform in ranked_isoforms {
        let Some(contig) = contigs_by_id.get(&isoform.source_contig_id) else {
            continue;
        };
        if isoform.direct_read_support == 0
            || (require_pair_support && isoform.direct_pair_support == 0)
        {
            continue;
        }
        if novel_sequence_fraction(&contig.sequence, &selected_sequences)
            >= MIN_SELECTED_NOVEL_SEQUENCE_FRACTION
        {
            selected_sequences.push(&contig.sequence);
            selected.push(isoform);
        }
    }
    selected
}

fn novel_sequence_fraction(candidate: &str, selected_sequences: &[&str]) -> f64 {
    if candidate.is_empty() {
        return 0.0;
    }
    let best_overlap = selected_sequences
        .iter()
        .map(|selected| {
            longest_common_substring_len(candidate.as_bytes(), selected.as_bytes()).max(
                longest_common_substring_len(
                    candidate.as_bytes(),
                    reverse_complement(selected).as_bytes(),
                ),
            )
        })
        .max()
        .unwrap_or(0);
    (candidate.len().saturating_sub(best_overlap)) as f64 / candidate.len() as f64
}

fn longest_common_substring_len(left: &[u8], right: &[u8]) -> usize {
    if left.is_empty() || right.is_empty() {
        return 0;
    }
    let mut previous = vec![0usize; right.len() + 1];
    let mut best = 0usize;
    for left_base in left {
        let mut current = vec![0usize; right.len() + 1];
        for (idx, right_base) in right.iter().enumerate() {
            if left_base == right_base {
                let len = previous[idx] + 1;
                current[idx + 1] = len;
                best = best.max(len);
            }
        }
        previous = current;
    }
    best
}

fn selected_isoform_evidence_score(
    direct_read_support: usize,
    direct_pair_support: usize,
    overlapping_read_kmer_path_count: usize,
    max_overlapping_read_kmer_path_support: usize,
) -> usize {
    direct_read_support
        + direct_pair_support * 5
        + overlapping_read_kmer_path_count * 3
        + max_overlapping_read_kmer_path_support * 2
}

fn contig_read_pair_support(
    contig: &str,
    reads1: &[String],
    reads2: Option<&[String]>,
) -> (usize, usize) {
    let read1_support = reads1
        .iter()
        .filter(|read| read_matches_contig(read, contig))
        .count();
    let Some(reads2) = reads2 else {
        return (read1_support, 0);
    };
    let mut read_support = read1_support;
    let mut pair_support = 0;
    for (read1, read2) in reads1.iter().zip(reads2.iter()) {
        let left_matches = read_matches_contig(read1, contig);
        let right_matches = read_matches_contig(read2, contig);
        if right_matches {
            read_support += 1;
        }
        if left_matches || right_matches {
            pair_support += 1;
        }
    }
    (read_support, pair_support)
}

fn read_matches_contig(read: &str, contig: &str) -> bool {
    contig.contains(read) || contig.contains(&reverse_complement(read))
}

fn sequences_overlap(left: &str, right: &str) -> bool {
    left.contains(right)
        || right.contains(left)
        || left.contains(&reverse_complement(right))
        || right.contains(&reverse_complement(left))
}

fn write_wrapped_fasta_sequence(output: &mut String, sequence: &str) {
    for chunk in sequence.as_bytes().chunks(80) {
        output.push_str(&String::from_utf8_lossy(chunk));
        output.push('\n');
    }
}

fn write_json<T: Serialize>(path: &Path, value: &T, label: &str) -> io::Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let file = fs::File::create(path)?;
    serde_json::to_writer_pretty(file, value).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("failed to serialize {label}: {err}"),
        )
    })
}

fn read_fastq_sequences(path: &str) -> io::Result<Vec<String>> {
    let reader = try_open_fastq(path)?;
    stream_fastq_records_checked(reader)
        .map(|record| record.map(|record| record.sequence))
        .collect()
}

fn read_contigs_from_fasta(path: &str) -> io::Result<Vec<Contig>> {
    let mut reader = try_open_fasta(path)?;
    let mut line = String::new();
    let mut contigs = Vec::new();
    let mut current_name: Option<String> = None;
    let mut sequence = String::new();

    loop {
        line.clear();
        let bytes = reader.read_line(&mut line)?;
        if bytes == 0 {
            break;
        }
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(header) = trimmed.strip_prefix('>') {
            if current_name.take().is_some() {
                contigs.push(Contig {
                    id: contigs.len(),
                    sequence: std::mem::take(&mut sequence),
                    kmer_path: Vec::new(),
                });
            }
            current_name = Some(header.to_string());
        } else {
            sequence.push_str(trimmed);
        }
    }

    if current_name.is_some() {
        contigs.push(Contig {
            id: contigs.len(),
            sequence,
            kmer_path: Vec::new(),
        });
    }

    Ok(contigs)
}

fn resolve_workflow_inputs(
    config: &TrinityWorkflowConfig,
    output_dir: &Path,
) -> io::Result<ResolvedWorkflowInputs> {
    if let Some(samples_file) = &config.samples_file {
        if config.input1.is_some() || config.input2.is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--samples-file cannot be combined with --input1/--input2",
            ));
        }
        let samples = read_samples_file(samples_file)?;
        return materialize_sample_inputs(&samples, output_dir);
    }

    let Some(input1) = &config.input1 else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "raptor trinity requires --input1 or --samples-file",
        ));
    };
    materialize_direct_inputs(input1, config.input2.as_deref(), output_dir)
}

fn split_input_list(input: &str) -> Vec<String> {
    input
        .split(',')
        .map(str::trim)
        .filter(|path| !path.is_empty())
        .map(str::to_string)
        .collect()
}

fn validate_ss_lib_type(ss_lib_type: Option<&str>) -> io::Result<Option<String>> {
    let Some(raw_value) = ss_lib_type else {
        return Ok(None);
    };
    let normalized = raw_value.trim().to_ascii_uppercase();
    match normalized.as_str() {
        "F" | "R" | "FR" | "RF" => Ok(Some(normalized)),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "--SS_lib_type must be one of F, R, FR, or RF; got '{}'",
                raw_value
            ),
        )),
    }
}

fn read_samples_file(path: &str) -> io::Result<Vec<TrinitySampleInput>> {
    let file = fs::File::open(path)?;
    let reader = io::BufReader::new(file);
    let mut samples = Vec::new();
    for (line_idx, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = if trimmed.contains('\t') {
            trimmed.split('\t').collect()
        } else {
            trimmed.split(',').collect()
        };
        let fields: Vec<&str> = fields
            .into_iter()
            .map(str::trim)
            .filter(|field| !field.is_empty())
            .collect();
        if fields.len() < 3 || fields.len() > 4 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "samples file line {}: expected condition replicate left [right], got {} field(s)",
                    line_idx + 1,
                    fields.len()
                ),
            ));
        }
        if samples.is_empty()
            && fields[0].eq_ignore_ascii_case("condition")
            && fields[1].eq_ignore_ascii_case("replicate")
        {
            continue;
        }
        samples.push(TrinitySampleInput {
            condition: fields[0].to_string(),
            replicate: fields[1].to_string(),
            left: fields[2].to_string(),
            right: fields.get(3).map(|field| field.to_string()),
        });
    }
    if samples.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "samples file did not contain any input rows",
        ));
    }
    let paired = samples[0].right.is_some();
    if samples
        .iter()
        .any(|sample| sample.right.is_some() != paired)
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "samples file cannot mix single-end and paired-end rows",
        ));
    }
    Ok(samples)
}

fn materialize_sample_inputs(
    samples: &[TrinitySampleInput],
    output_dir: &Path,
) -> io::Result<ResolvedWorkflowInputs> {
    if samples.len() == 1 {
        return Ok(ResolvedWorkflowInputs {
            input1: samples[0].left.clone(),
            input2: samples[0].right.clone(),
            sample_count: 1,
        });
    }

    let merged_left = output_dir.join("samples_file_left.fastq.gz");
    let merged_left_str = path_to_str(&merged_left)?.to_string();
    let mut left_writer = FastqWriter::try_new(&merged_left_str)?;
    let mut right_writer = if samples[0].right.is_some() {
        let merged_right = output_dir.join("samples_file_right.fastq.gz");
        let merged_right_str = path_to_str(&merged_right)?.to_string();
        Some((
            merged_right_str,
            FastqWriter::try_new(path_to_str(&merged_right)?)?,
        ))
    } else {
        None
    };

    for sample in samples {
        copy_fastq_records(&sample.left, &mut left_writer)?;
        match (&sample.right, &mut right_writer) {
            (Some(right), Some((_, writer))) => copy_fastq_records(right, writer)?,
            (None, None) => {}
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "samples file cannot mix single-end and paired-end rows",
                ));
            }
        }
    }

    Ok(ResolvedWorkflowInputs {
        input1: merged_left_str,
        input2: right_writer.map(|(path, _)| path),
        sample_count: samples.len(),
    })
}

fn materialize_direct_inputs(
    input1: &str,
    input2: Option<&str>,
    output_dir: &Path,
) -> io::Result<ResolvedWorkflowInputs> {
    let left_paths = split_input_list(input1);
    if left_paths.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--input1 did not contain any read paths",
        ));
    }
    let right_paths = input2.map(split_input_list);
    if let Some(right_paths) = &right_paths {
        if right_paths.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--input2 did not contain any read paths",
            ));
        }
        if right_paths.len() != left_paths.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "--input1 and --input2 comma-separated lists must contain the same number of paths ({} != {})",
                    left_paths.len(),
                    right_paths.len()
                ),
            ));
        }
    }

    if left_paths.len() == 1 {
        return Ok(ResolvedWorkflowInputs {
            input1: left_paths[0].clone(),
            input2: right_paths.map(|paths| paths[0].clone()),
            sample_count: 1,
        });
    }

    let merged_left = output_dir.join("comma_input_left.fastq.gz");
    let merged_left_str = path_to_str(&merged_left)?.to_string();
    {
        let mut left_writer = FastqWriter::try_new(&merged_left_str)?;
        for left_path in &left_paths {
            copy_fastq_records(left_path, &mut left_writer)?;
        }
    }

    let merged_right_str = if let Some(right_paths) = right_paths {
        let merged_right = output_dir.join("comma_input_right.fastq.gz");
        let merged_right_str = path_to_str(&merged_right)?.to_string();
        {
            let mut right_writer = FastqWriter::try_new(&merged_right_str)?;
            for right_path in &right_paths {
                copy_fastq_records(right_path, &mut right_writer)?;
            }
        }
        Some(merged_right_str)
    } else {
        None
    };

    Ok(ResolvedWorkflowInputs {
        input1: merged_left_str,
        input2: merged_right_str,
        sample_count: left_paths.len(),
    })
}

fn copy_fastq_records(path: &str, writer: &mut FastqWriter) -> io::Result<()> {
    let reader = try_open_fastq(path)?;
    for record in stream_fastq_records_checked(reader) {
        writer.write_record(&record?)?;
    }
    Ok(())
}

fn normalize_for_workflow(
    input1: &str,
    input2: Option<&str>,
    config: &TrinityWorkflowConfig,
    output_dir: &Path,
) -> io::Result<(String, Option<String>, Option<NormalizeSummary>)> {
    let prefix = output_dir.join("normalized_reads");
    let prefix_str = path_to_str(&prefix)?;
    if let Some(input2) = input2 {
        let summary = normalize_paired_with_config(
            input1,
            input2,
            prefix_str,
            config.normalize_config,
            true,
        )?;
        Ok((
            format!("{prefix_str}_R1.fastq.gz"),
            Some(format!("{prefix_str}_R2.fastq.gz")),
            Some(summary),
        ))
    } else {
        let summary =
            normalize_single_with_config(input1, prefix_str, config.normalize_config, true)?;
        Ok((format!("{prefix_str}.fastq.gz"), None, Some(summary)))
    }
}

fn write_report(path: &Path, report: &TrinityWorkflowReport) -> io::Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let file = fs::File::create(path)?;
    serde_json::to_writer_pretty(file, report).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("failed to serialize Trinity workflow report: {err}"),
        )
    })
}

fn path_to_str(path: &Path) -> io::Result<&str> {
    path.to_str()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "path is not valid UTF-8"))
}

fn sidecar_path(output_path: &str, suffix: &str) -> String {
    if output_path.ends_with(".fa.gz") {
        output_path.replace(".fa.gz", &format!(".{suffix}"))
    } else if output_path.ends_with(".fasta.gz") {
        output_path.replace(".fasta.gz", &format!(".{suffix}"))
    } else if output_path.ends_with(".fa") {
        output_path.replace(".fa", &format!(".{suffix}"))
    } else if output_path.ends_with(".fasta") {
        output_path.replace(".fasta", &format!(".{suffix}"))
    } else {
        format!("{output_path}.{suffix}")
    }
}

#[cfg(test)]
mod tests {
    use super::{run_trinity_workflow, TrinityWorkflowConfig, COMPONENT_CONTIG_SELECTION_METHOD};
    use crate::pipeline::normalize::NormalizeConfig;
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn trinity_workflow_runs_normalize_then_assemble_on_tiny_paired_reads() {
        let temp_dir = TempDir::new().expect("temp dir");
        let r1 = temp_dir.path().join("r1.fastq");
        let r2 = temp_dir.path().join("r2.fastq");
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        fs::write(&r1, format!("@r1/1\n{seq}\n+\n{}\n", "I".repeat(seq.len()))).expect("write r1");
        fs::write(&r2, format!("@r1/2\n{seq}\n+\n{}\n", "I".repeat(seq.len()))).expect("write r2");

        let output_dir = temp_dir.path().join("workflow");
        let report = run_trinity_workflow(TrinityWorkflowConfig {
            input1: Some(r1.to_string_lossy().into_owned()),
            input2: Some(r2.to_string_lossy().into_owned()),
            ss_lib_type: None,
            samples_file: None,
            output_dir: output_dir.to_string_lossy().into_owned(),
            output_fasta: None,
            report_json: None,
            normalize: true,
            normalize_config: NormalizeConfig {
                k: 5,
                target_coverage: u16::MAX,
                min_abundance: 1,
                max_reads: None,
                use_gpu: false,
            },
            min_len: 10,
            use_gpu: false,
        })
        .expect("workflow should run");

        assert!(report.normalized);
        assert_eq!(
            report
                .normalization
                .expect("normalization summary")
                .total_reads,
            1
        );
        assert_eq!(report.component_count, 1);
        assert_eq!(report.component_clustering, "sequence_or_read_kmer");
        assert!(std::path::Path::new(&report.assembly_fasta).exists());
        assert!(std::path::Path::new(&report.components_json).exists());
        assert!(std::path::Path::new(&report.component_graphs_json).exists());
        assert!(std::path::Path::new(&report.component_paths_fasta).exists());
        assert!(std::path::Path::new(&report.component_selected_isoforms_fasta).exists());
        assert!(std::path::Path::new(&report.component_selected_isoforms_json).exists());
        assert!(std::path::Path::new(&report.component_isoform_candidates_json).exists());
        let components: Vec<serde_json::Value> = serde_json::from_str(
            &fs::read_to_string(&report.components_json).expect("component json"),
        )
        .expect("parse components");
        assert_eq!(components[0]["assigned_read_count"], 2);
        assert_eq!(components[0]["assigned_pair_count"], 1);
        let component_graphs: Vec<serde_json::Value> = serde_json::from_str(
            &fs::read_to_string(&report.component_graphs_json).expect("component graph json"),
        )
        .expect("parse component graphs");
        assert_eq!(component_graphs[0]["node_count"], 1);
        assert_eq!(component_graphs[0]["assigned_pair_count"], 1);
        let selected_isoforms: Vec<serde_json::Value> = serde_json::from_str(
            &fs::read_to_string(&report.component_selected_isoforms_json)
                .expect("selected isoform json"),
        )
        .expect("parse selected isoforms");
        assert_eq!(selected_isoforms.len(), 1);
        assert_eq!(selected_isoforms[0]["component_rank"], 1);
        assert_eq!(
            selected_isoforms[0]["selection_method"],
            COMPONENT_CONTIG_SELECTION_METHOD
        );
        assert!(
            selected_isoforms[0]["evidence_score"]
                .as_u64()
                .expect("evidence score")
                >= 1
        );
        assert_eq!(selected_isoforms[0]["component_assigned_pairs"], 1);
        assert!(
            selected_isoforms[0]["direct_read_support"]
                .as_u64()
                .expect("direct read support")
                >= 1
        );
        let isoform_candidates: Vec<serde_json::Value> = serde_json::from_str(
            &fs::read_to_string(&report.component_isoform_candidates_json)
                .expect("isoform candidate json"),
        )
        .expect("parse isoform candidates");
        assert!(isoform_candidates.iter().any(|candidate| {
            candidate["source_kind"] == "contig" && candidate["selected"] == true
        }));
        assert!(output_dir.join("raptor_trinity_report.json").exists());
    }

    #[test]
    fn trinity_workflow_selects_isoforms_for_single_end_reads_without_pair_support() {
        let temp_dir = TempDir::new().expect("temp dir");
        let r1 = temp_dir.path().join("single.fastq");
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        fs::write(&r1, format!("@r1\n{seq}\n+\n{}\n", "I".repeat(seq.len())))
            .expect("write single-end reads");

        let output_dir = temp_dir.path().join("workflow_single");
        let report = run_trinity_workflow(TrinityWorkflowConfig {
            input1: Some(r1.to_string_lossy().into_owned()),
            input2: None,
            ss_lib_type: None,
            samples_file: None,
            output_dir: output_dir.to_string_lossy().into_owned(),
            output_fasta: None,
            report_json: None,
            normalize: true,
            normalize_config: NormalizeConfig {
                k: 5,
                target_coverage: u16::MAX,
                min_abundance: 1,
                max_reads: None,
                use_gpu: false,
            },
            min_len: 10,
            use_gpu: false,
        })
        .expect("single-end workflow should run");

        let selected_isoforms: Vec<serde_json::Value> = serde_json::from_str(
            &fs::read_to_string(&report.component_selected_isoforms_json)
                .expect("selected isoform json"),
        )
        .expect("parse selected isoforms");
        assert_eq!(selected_isoforms.len(), 1);
        assert_eq!(selected_isoforms[0]["component_assigned_pairs"], 0);
        assert_eq!(selected_isoforms[0]["direct_pair_support"], 0);
        assert!(
            selected_isoforms[0]["direct_read_support"]
                .as_u64()
                .expect("direct read support")
                >= 1
        );
    }

    #[test]
    fn trinity_workflow_accepts_trinity_style_samples_file() {
        let temp_dir = TempDir::new().expect("temp dir");
        let r1 = temp_dir.path().join("rep1_R1.fastq");
        let r2 = temp_dir.path().join("rep1_R2.fastq");
        let samples_file = temp_dir.path().join("samples.tsv");
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        fs::write(&r1, format!("@r1/1\n{seq}\n+\n{}\n", "I".repeat(seq.len()))).expect("write r1");
        fs::write(&r2, format!("@r1/2\n{seq}\n+\n{}\n", "I".repeat(seq.len()))).expect("write r2");
        fs::write(
            &samples_file,
            format!(
                "condition\treplicate\tleft\tright\ncondA\trep1\t{}\t{}\n",
                r1.display(),
                r2.display()
            ),
        )
        .expect("write samples file");

        let output_dir = temp_dir.path().join("workflow_samples");
        let report = run_trinity_workflow(TrinityWorkflowConfig {
            input1: None,
            input2: None,
            ss_lib_type: None,
            samples_file: Some(samples_file.to_string_lossy().into_owned()),
            output_dir: output_dir.to_string_lossy().into_owned(),
            output_fasta: None,
            report_json: None,
            normalize: true,
            normalize_config: NormalizeConfig {
                k: 5,
                target_coverage: u16::MAX,
                min_abundance: 1,
                max_reads: None,
                use_gpu: false,
            },
            min_len: 10,
            use_gpu: false,
        })
        .expect("samples-file workflow should run");

        assert_eq!(report.sample_count, 1);
        assert_eq!(
            report.samples_file.as_deref(),
            Some(samples_file.to_string_lossy().as_ref())
        );
        assert_eq!(report.input1, r1.to_string_lossy());
        assert_eq!(
            report.input2.as_deref(),
            Some(r2.to_string_lossy().as_ref())
        );
        assert!(std::path::Path::new(&report.assembly_fasta).exists());
    }

    #[test]
    fn trinity_workflow_merges_multi_row_samples_file() {
        let temp_dir = TempDir::new().expect("temp dir");
        let r1 = temp_dir.path().join("rep_R1.fastq");
        let r2 = temp_dir.path().join("rep_R2.fastq");
        let samples_file = temp_dir.path().join("samples.tsv");
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        fs::write(&r1, format!("@r1/1\n{seq}\n+\n{}\n", "I".repeat(seq.len()))).expect("write r1");
        fs::write(&r2, format!("@r1/2\n{seq}\n+\n{}\n", "I".repeat(seq.len()))).expect("write r2");
        fs::write(
            &samples_file,
            format!(
                "condition\treplicate\tleft\tright\ncondA\trep1\t{}\t{}\ncondA\trep2\t{}\t{}\n",
                r1.display(),
                r2.display(),
                r1.display(),
                r2.display()
            ),
        )
        .expect("write samples file");

        let output_dir = temp_dir.path().join("workflow_samples_multi");
        let report = run_trinity_workflow(TrinityWorkflowConfig {
            input1: None,
            input2: None,
            ss_lib_type: None,
            samples_file: Some(samples_file.to_string_lossy().into_owned()),
            output_dir: output_dir.to_string_lossy().into_owned(),
            output_fasta: None,
            report_json: None,
            normalize: true,
            normalize_config: NormalizeConfig {
                k: 5,
                target_coverage: u16::MAX,
                min_abundance: 1,
                max_reads: None,
                use_gpu: false,
            },
            min_len: 10,
            use_gpu: false,
        })
        .expect("multi-row samples-file workflow should run");

        assert_eq!(report.sample_count, 2);
        assert!(report.input1.ends_with("samples_file_left.fastq.gz"));
        assert!(report
            .input2
            .as_deref()
            .expect("merged right FASTQ")
            .ends_with("samples_file_right.fastq.gz"));
        assert!(std::path::Path::new(&report.input1).exists());
        assert!(std::path::Path::new(report.input2.as_ref().unwrap()).exists());
        assert!(std::path::Path::new(&report.assembly_fasta).exists());
    }

    #[test]
    fn trinity_workflow_merges_comma_separated_direct_inputs() {
        let temp_dir = TempDir::new().expect("temp dir");
        let r1 = temp_dir.path().join("rep_R1.fastq");
        let r2 = temp_dir.path().join("rep_R2.fastq");
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        fs::write(&r1, format!("@r1/1\n{seq}\n+\n{}\n", "I".repeat(seq.len()))).expect("write r1");
        fs::write(&r2, format!("@r1/2\n{seq}\n+\n{}\n", "I".repeat(seq.len()))).expect("write r2");

        let output_dir = temp_dir.path().join("workflow_comma_inputs");
        let report = run_trinity_workflow(TrinityWorkflowConfig {
            input1: Some(format!("{},{}", r1.display(), r1.display())),
            input2: Some(format!("{},{}", r2.display(), r2.display())),
            ss_lib_type: None,
            samples_file: None,
            output_dir: output_dir.to_string_lossy().into_owned(),
            output_fasta: None,
            report_json: None,
            normalize: true,
            normalize_config: NormalizeConfig {
                k: 5,
                target_coverage: u16::MAX,
                min_abundance: 1,
                max_reads: None,
                use_gpu: false,
            },
            min_len: 10,
            use_gpu: false,
        })
        .expect("comma-separated direct input workflow should run");

        assert_eq!(report.sample_count, 2);
        assert!(report.input1.ends_with("comma_input_left.fastq.gz"));
        assert!(report
            .input2
            .as_deref()
            .expect("merged right FASTQ")
            .ends_with("comma_input_right.fastq.gz"));
        assert!(std::path::Path::new(&report.input1).exists());
        assert!(std::path::Path::new(report.input2.as_ref().unwrap()).exists());
        assert!(std::path::Path::new(&report.assembly_fasta).exists());
    }

    #[test]
    fn trinity_workflow_rejects_mismatched_comma_input_lists() {
        let temp_dir = TempDir::new().expect("temp dir");
        let r1 = temp_dir.path().join("rep_R1.fastq");
        let r2 = temp_dir.path().join("rep_R2.fastq");
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        fs::write(&r1, format!("@r1/1\n{seq}\n+\n{}\n", "I".repeat(seq.len()))).expect("write r1");
        fs::write(&r2, format!("@r1/2\n{seq}\n+\n{}\n", "I".repeat(seq.len()))).expect("write r2");

        let output_dir = temp_dir.path().join("workflow_bad_comma_inputs");
        let err = run_trinity_workflow(TrinityWorkflowConfig {
            input1: Some(format!("{},{}", r1.display(), r1.display())),
            input2: Some(r2.to_string_lossy().into_owned()),
            ss_lib_type: None,
            samples_file: None,
            output_dir: output_dir.to_string_lossy().into_owned(),
            output_fasta: None,
            report_json: None,
            normalize: true,
            normalize_config: NormalizeConfig {
                k: 5,
                target_coverage: u16::MAX,
                min_abundance: 1,
                max_reads: None,
                use_gpu: false,
            },
            min_len: 10,
            use_gpu: false,
        })
        .expect_err("mismatched comma input lists should fail");

        assert_eq!(err.kind(), std::io::ErrorKind::InvalidInput);
        assert!(err.to_string().contains("same number of paths"));
    }

    #[test]
    fn trinity_workflow_reports_valid_strand_specific_library_type() {
        let temp_dir = TempDir::new().expect("temp dir");
        let r1 = temp_dir.path().join("r1.fastq");
        let r2 = temp_dir.path().join("r2.fastq");
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        fs::write(&r1, format!("@r1/1\n{seq}\n+\n{}\n", "I".repeat(seq.len()))).expect("write r1");
        fs::write(&r2, format!("@r1/2\n{seq}\n+\n{}\n", "I".repeat(seq.len()))).expect("write r2");

        let output_dir = temp_dir.path().join("workflow_stranded");
        let report = run_trinity_workflow(TrinityWorkflowConfig {
            input1: Some(r1.to_string_lossy().into_owned()),
            input2: Some(r2.to_string_lossy().into_owned()),
            ss_lib_type: Some("rf".to_string()),
            samples_file: None,
            output_dir: output_dir.to_string_lossy().into_owned(),
            output_fasta: None,
            report_json: None,
            normalize: true,
            normalize_config: NormalizeConfig {
                k: 5,
                target_coverage: u16::MAX,
                min_abundance: 1,
                max_reads: None,
                use_gpu: false,
            },
            min_len: 10,
            use_gpu: false,
        })
        .expect("stranded workflow should run");

        assert_eq!(report.ss_lib_type.as_deref(), Some("RF"));
        let report_payload: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(&report.report_json).expect("report json"))
                .expect("parse workflow report");
        assert_eq!(report_payload["ss_lib_type"], "RF");
    }

    #[test]
    fn trinity_workflow_rejects_invalid_strand_specific_library_type() {
        let temp_dir = TempDir::new().expect("temp dir");
        let r1 = temp_dir.path().join("r1.fastq");
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        fs::write(&r1, format!("@r1\n{seq}\n+\n{}\n", "I".repeat(seq.len()))).expect("write r1");

        let output_dir = temp_dir.path().join("workflow_bad_stranded");
        let err = run_trinity_workflow(TrinityWorkflowConfig {
            input1: Some(r1.to_string_lossy().into_owned()),
            input2: None,
            ss_lib_type: Some("bad".to_string()),
            samples_file: None,
            output_dir: output_dir.to_string_lossy().into_owned(),
            output_fasta: None,
            report_json: None,
            normalize: true,
            normalize_config: NormalizeConfig {
                k: 5,
                target_coverage: u16::MAX,
                min_abundance: 1,
                max_reads: None,
                use_gpu: false,
            },
            min_len: 10,
            use_gpu: false,
        })
        .expect_err("invalid strand type should fail");

        assert_eq!(err.kind(), std::io::ErrorKind::InvalidInput);
        assert!(err.to_string().contains("--SS_lib_type"));
    }
}
