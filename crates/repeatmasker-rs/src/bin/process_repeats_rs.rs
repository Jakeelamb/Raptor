use clap::Parser;
use repeatmasker_rs::{
    adjudicate_annotations, adjudicate_annotations_with_alignments, catalog_stats_for_annotations,
    catalog_stats_for_annotations_with_alignments, compute_repeatmasker_mask_stats_from_fasta,
    encode_repeatmasker_alignment, load_repeatmasker_alignments_path,
    load_repeatmasker_catalog_path, write_alignment_tsv_path, write_annotations_tsv_path,
    write_chain_tsv_path, write_chain_tsv_with_alignments_path, write_gff3_path,
    write_out_annotations_path, write_repeat_summary_tsv_path, write_repeatmasker_tbl_path,
    write_repeatmasker_tbl_with_alignments_path, MaskMode, RepeatMaskerCatalog, RepeatMaskerError,
};
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Debug, Parser)]
#[command(
    name = "process_repeats_rs",
    about = "Rust RepeatMasker ProcessRepeats prototype for .out/.cat parsing, summary export, and masking"
)]
struct Cli {
    /// RepeatMasker annotation file to parse (.out, .cat, .gz)
    #[arg(long)]
    annotations: PathBuf,

    /// Optional normalized .out-style output path
    #[arg(long)]
    out: Option<PathBuf>,

    /// Optional GFF3 output path
    #[arg(long)]
    gff: Option<PathBuf>,

    /// Optional tabular TSV output path
    #[arg(long)]
    tsv: Option<PathBuf>,

    /// Optional repeat family/class summary TSV path
    #[arg(long)]
    summary_tsv: Option<PathBuf>,

    /// Optional RepeatMasker-style .tbl output path
    #[arg(long)]
    tbl: Option<PathBuf>,

    /// Optional source FASTA for exact stats and optional masking
    #[arg(long)]
    fasta: Option<PathBuf>,

    /// Optional RepeatMasker .align file to refine report-level adjudication
    #[arg(long)]
    align: Option<PathBuf>,

    /// Optional alignment TSV export derived from --align
    #[arg(long)]
    align_tsv: Option<PathBuf>,

    /// Optional fragment-chain TSV export derived from adjudicated annotations
    #[arg(long)]
    chain_tsv: Option<PathBuf>,

    /// Optional masked FASTA output path; requires --fasta
    #[arg(long)]
    masked_output: Option<PathBuf>,

    /// Use lowercase instead of N masking when writing --masked-output
    #[arg(long)]
    xsmall: bool,

    /// Use X instead of N masking when writing --masked-output
    #[arg(long)]
    x: bool,

    /// Optional JSON summary path
    #[arg(long)]
    stats_json: Option<PathBuf>,

    /// FASTA line width for --masked-output
    #[arg(long, default_value_t = 50)]
    wrap_width: usize,
}

fn selected_mask_mode(cli: &Cli) -> MaskMode {
    if cli.xsmall {
        MaskMode::Lowercase
    } else if cli.x {
        MaskMode::X
    } else {
        MaskMode::N
    }
}

fn source_name_from_annotation_path(path: &Path) -> String {
    let file_name = path
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("annotations");
    let stripped = file_name.strip_suffix(".gz").unwrap_or(file_name);
    stripped
        .strip_suffix(".cat")
        .or_else(|| stripped.strip_suffix(".out"))
        .unwrap_or(stripped)
        .to_string()
}

fn main() -> Result<(), RepeatMaskerError> {
    let cli = Cli::parse();
    let catalog = load_repeatmasker_catalog_path(&cli.annotations)?;
    let annotations = &catalog.annotations;
    let alignments = cli
        .align
        .as_deref()
        .map(load_repeatmasker_alignments_path)
        .transpose()?;
    let adjudicated_annotations = alignments.as_deref().map_or_else(
        || adjudicate_annotations(annotations),
        |alignments| adjudicate_annotations_with_alignments(annotations, alignments),
    );
    let report_catalog = RepeatMaskerCatalog {
        annotations: adjudicated_annotations.clone(),
        metadata: catalog.metadata.clone(),
    };
    let computed_mask_stats = if let Some(fasta_path) = &cli.fasta {
        Some(compute_repeatmasker_mask_stats_from_fasta(
            &adjudicated_annotations,
            fasta_path,
            selected_mask_mode(&cli),
            cli.masked_output
                .as_deref()
                .map(|path| (path, cli.wrap_width.max(1))),
        )?)
    } else if cli.masked_output.is_some() {
        return Err(RepeatMaskerError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "--masked-output requires --fasta",
        )));
    } else {
        None
    };

    if let Some(path) = &cli.out {
        write_out_annotations_path(path, annotations)?;
    }
    if let Some(path) = &cli.gff {
        write_gff3_path(path, annotations)?;
    }
    if let Some(path) = &cli.tsv {
        write_annotations_tsv_path(path, annotations)?;
    }
    if let Some(path) = &cli.summary_tsv {
        write_repeat_summary_tsv_path(path, &report_catalog)?;
    }
    if let Some(path) = &cli.tbl {
        if let Some(alignments) = &alignments {
            write_repeatmasker_tbl_with_alignments_path(
                path,
                &source_name_from_annotation_path(&cli.annotations),
                &catalog,
                computed_mask_stats.as_ref(),
                alignments,
            )?;
        } else {
            write_repeatmasker_tbl_path(
                path,
                &source_name_from_annotation_path(&cli.annotations),
                &catalog,
                computed_mask_stats.as_ref(),
            )?;
        }
    }
    if let Some(path) = &cli.align_tsv {
        let Some(alignments) = &alignments else {
            return Err(RepeatMaskerError::Io(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "--align-tsv requires --align",
            )));
        };
        write_alignment_tsv_path(path, alignments)?;
    }
    if let Some(path) = &cli.chain_tsv {
        if let Some(alignments) = &alignments {
            write_chain_tsv_with_alignments_path(path, &adjudicated_annotations, alignments)?;
        } else {
            write_chain_tsv_path(path, &adjudicated_annotations)?;
        }
    }

    let mut summary = serde_json::to_value(alignments.as_deref().map_or_else(
        || {
            catalog_stats_for_annotations(
                &adjudicated_annotations,
                Some(annotations.len()),
                catalog.metadata.clone(),
            )
        },
        |alignments| {
            catalog_stats_for_annotations_with_alignments(
                &adjudicated_annotations,
                Some(annotations.len()),
                catalog.metadata.clone(),
                alignments,
            )
        },
    ))
    .map_err(std::io::Error::other)?;
    if let Some(alignments) = &alignments {
        let aggregate_match_bases: usize = alignments
            .iter()
            .map(|record| encode_repeatmasker_alignment(record).match_bases)
            .sum();
        let aggregate_mismatch_bases: usize = alignments
            .iter()
            .map(|record| encode_repeatmasker_alignment(record).mismatch_bases)
            .sum();
        summary["alignments"] = serde_json::json!({
            "record_count": alignments.len(),
            "linked_record_count": alignments
                .iter()
                .filter(|record| !record.annotation_id.is_empty())
                .count(),
            "match_bases": aggregate_match_bases,
            "mismatch_bases": aggregate_mismatch_bases,
        });
    }
    if let Some(mask_stats) = &computed_mask_stats {
        summary["masking"] = serde_json::to_value(mask_stats).map_err(std::io::Error::other)?;
    }
    summary["chains"] = serde_json::json!({
        "count": summary["element_chain_count"].as_u64().unwrap_or(0),
    });

    if let Some(path) = cli.stats_json {
        let json = serde_json::to_string_pretty(&summary).map_err(std::io::Error::other)?;
        fs::write(path, json + "\n")?;
    }

    println!(
        "Parsed {} raw annotations; adjudicated to {} across {} sequences in {} element chains",
        annotations.len(),
        adjudicated_annotations.len(),
        summary["sequence_count"].as_u64().unwrap_or(0),
        summary["chains"]["count"].as_u64().unwrap_or(0)
    );
    Ok(())
}
