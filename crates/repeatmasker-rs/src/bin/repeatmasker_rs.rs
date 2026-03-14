use clap::{Parser, ValueEnum};
use repeatmasker_rs::{mask_repeatmasker_out_to_fasta, MaskMode, RepeatMaskerError};
use std::fs;
use std::path::PathBuf;

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
enum MaskModeArg {
    N,
    X,
    Lowercase,
}

impl From<MaskModeArg> for MaskMode {
    fn from(value: MaskModeArg) -> Self {
        match value {
            MaskModeArg::N => MaskMode::N,
            MaskModeArg::X => MaskMode::X,
            MaskModeArg::Lowercase => MaskMode::Lowercase,
        }
    }
}

#[derive(Debug, Parser)]
#[command(
    name = "repeatmasker-rs",
    about = "Rust RepeatMasker postprocessor prototype for masking FASTA from .out annotations"
)]
struct Cli {
    /// RepeatMasker .out annotation file
    #[arg(long)]
    annotations: PathBuf,

    /// Source FASTA to mask
    #[arg(long)]
    fasta: PathBuf,

    /// Output FASTA path for the masked sequence
    #[arg(long)]
    output: PathBuf,

    /// Masking mode: N, X, or lowercase
    #[arg(long, value_enum, default_value = "n")]
    mask: MaskModeArg,

    /// Optional path to write machine-readable masking stats as JSON
    #[arg(long)]
    stats_json: Option<PathBuf>,

    /// FASTA line width in the output file
    #[arg(long, default_value_t = 50)]
    wrap_width: usize,
}

fn main() -> Result<(), RepeatMaskerError> {
    let cli = Cli::parse();
    let stats = mask_repeatmasker_out_to_fasta(
        &cli.annotations,
        &cli.fasta,
        &cli.output,
        cli.mask.into(),
        cli.wrap_width.max(1),
    )?;

    if let Some(stats_path) = cli.stats_json {
        let json = serde_json::to_string_pretty(&stats).map_err(std::io::Error::other)?;
        fs::write(stats_path, json + "\n")?;
    }

    println!(
        "Masked {} sequences, {} bp total, {} bp masked, GC {:.4}",
        stats.sequence_count, stats.total_bases, stats.masked_bases, stats.gc_fraction
    );
    Ok(())
}
