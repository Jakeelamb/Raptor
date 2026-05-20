use raptor::pipeline::normalize::{normalize_paired_with_config, NormalizeConfig};
use std::fs;

fn main() {
    let args: Vec<String> = std::env::args().collect();

    let help_requested = args.iter().any(|arg| arg == "--help" || arg == "-h");
    if args.len() < 4 || help_requested {
        print_usage(&args[0]);
        std::process::exit(if help_requested { 0 } else { 1 });
    }

    let defaults = NormalizeConfig::default();
    let k = parse_or_default(args.get(4), defaults.k);
    let target = parse_or_default(args.get(5), defaults.target_coverage);
    let min_abundance = parse_or_default(args.get(6), defaults.min_abundance);
    let use_gpu = args.iter().any(|arg| arg == "--gpu");
    let streaming = true;

    println!(
        "Running paired-end normalization with parameters: k={}, target={}, min_abund={}",
        k, target, min_abundance
    );

    let pipeline_output_prefix = format!("{}_norm", args[3]);
    match normalize_paired_with_config(
        &args[1],
        &args[2],
        &pipeline_output_prefix,
        NormalizeConfig {
            k,
            target_coverage: target,
            min_abundance,
            use_gpu,
            ..NormalizeConfig::default()
        },
        streaming,
    ) {
        Ok(summary) => {
            if let Err(err) =
                move_pipeline_outputs_to_legacy_paths(&args[3], &pipeline_output_prefix)
            {
                eprintln!("Failed to move normalized output files: {}", err);
                std::process::exit(1);
            }
            println!(
                "Paired normalization complete. Kept {}/{} read pairs ({:.1}%)",
                summary.kept_reads,
                summary.total_reads,
                percentage(summary.kept_reads, summary.total_reads)
            );
        }
        Err(err) => {
            eprintln!("Paired normalization failed: {}", err);
            std::process::exit(1);
        }
    }
}

fn move_pipeline_outputs_to_legacy_paths(
    output_prefix: &str,
    pipeline_output_prefix: &str,
) -> std::io::Result<()> {
    fs::rename(
        format!("{}_R1.fastq.gz", pipeline_output_prefix),
        format!("{}_R1.norm.fastq.gz", output_prefix),
    )?;
    fs::rename(
        format!("{}_R2.fastq.gz", pipeline_output_prefix),
        format!("{}_R2.norm.fastq.gz", output_prefix),
    )
}

fn parse_or_default<T>(value: Option<&String>, default: T) -> T
where
    T: std::str::FromStr,
{
    value
        .and_then(|value| value.parse().ok())
        .unwrap_or(default)
}

#[inline]
fn percentage(kept: usize, total: usize) -> f64 {
    if total == 0 {
        0.0
    } else {
        (kept as f64 / total as f64) * 100.0
    }
}

fn print_usage(program_name: &str) {
    eprintln!("Paired-end FASTQ Read Normalizer");
    eprintln!("--------------------------------");
    eprintln!(
        "Usage: {} <R1.fastq(.gz)> <R2.fastq(.gz)> <output_prefix> [k-size] [target] [min_abund] [--gpu]",
        program_name
    );
    eprintln!();
    eprintln!("Parameters:");
    eprintln!("  <R1.fastq(.gz)>    Input FASTQ file for read 1");
    eprintln!("  <R2.fastq(.gz)>    Input FASTQ file for read 2");
    eprintln!("  <output_prefix>    Prefix for output files");
    eprintln!("  [k-size]           K-mer size (default: 25)");
    eprintln!("  [target]           Target coverage (default: 50)");
    eprintln!("  [min_abund]        Minimum k-mer abundance to keep (default: 2)");
    eprintln!("  [--gpu]            Request GPU acceleration where supported");
    eprintln!();
    eprintln!("Output files:");
    eprintln!("  <output_prefix>_R1.norm.fastq.gz");
    eprintln!("  <output_prefix>_R2.norm.fastq.gz");
}
