use crate::graph::assembler::Contig;
use crate::graph::partition::{cluster_contigs_by_shared_sequence, RaptorComponent};
use crate::io::fasta::try_open_fasta;
use crate::pipeline::assemble::assemble_reads_with_gpu;
use crate::pipeline::normalize::{
    normalize_paired_with_config, normalize_single_with_config, NormalizeConfig, NormalizeSummary,
};
use serde::Serialize;
use std::fs;
use std::io::{self, BufRead};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Debug, Clone)]
pub struct TrinityWorkflowConfig {
    pub input1: String,
    pub input2: Option<String>,
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
    pub output_dir: String,
    pub normalized: bool,
    pub normalized_input1: Option<String>,
    pub normalized_input2: Option<String>,
    pub normalization: Option<NormalizeSummary>,
    pub assembly_fasta: String,
    pub assembly_metrics_json: String,
    pub assembly_metrics_tsv: String,
    pub components_json: String,
    pub component_count: usize,
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

    let (assembly_input1, assembly_input2, normalization) = if config.normalize {
        normalize_for_workflow(&config, &output_dir)?
    } else {
        (config.input1.clone(), config.input2.clone(), None)
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
    let components = write_component_report(&assembly_fasta_string, &component_path, 60)?;
    let report_path = config
        .report_json
        .as_ref()
        .map(PathBuf::from)
        .unwrap_or_else(|| output_dir.join("raptor_trinity_report.json"));

    let report = TrinityWorkflowReport {
        workflow: "raptor_trinity_de_novo",
        input1: config.input1,
        input2: config.input2,
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
        component_count: components.len(),
        assembly_fasta: assembly_fasta_string,
        report_json: path_to_str(&report_path)?.to_string(),
        elapsed_seconds: start.elapsed().as_secs_f64(),
    };

    write_report(&report_path, &report)?;
    Ok(report)
}

fn write_component_report(
    assembly_fasta: &str,
    component_path: &Path,
    min_shared_bases: usize,
) -> io::Result<Vec<RaptorComponent>> {
    let contigs = read_contigs_from_fasta(assembly_fasta)?;
    let components = cluster_contigs_by_shared_sequence(&contigs, min_shared_bases);
    if let Some(parent) = component_path.parent() {
        fs::create_dir_all(parent)?;
    }
    let file = fs::File::create(component_path)?;
    serde_json::to_writer_pretty(file, &components).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("failed to serialize component report: {err}"),
        )
    })?;
    Ok(components)
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

fn normalize_for_workflow(
    config: &TrinityWorkflowConfig,
    output_dir: &Path,
) -> io::Result<(String, Option<String>, Option<NormalizeSummary>)> {
    let prefix = output_dir.join("normalized_reads");
    let prefix_str = path_to_str(&prefix)?;
    if let Some(input2) = &config.input2 {
        let summary = normalize_paired_with_config(
            &config.input1,
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
        let summary = normalize_single_with_config(
            &config.input1,
            prefix_str,
            config.normalize_config,
            true,
        )?;
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
    use super::{run_trinity_workflow, TrinityWorkflowConfig};
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
            input1: r1.to_string_lossy().into_owned(),
            input2: Some(r2.to_string_lossy().into_owned()),
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
        assert!(std::path::Path::new(&report.assembly_fasta).exists());
        assert!(std::path::Path::new(&report.components_json).exists());
        assert!(output_dir.join("raptor_trinity_report.json").exists());
    }
}
