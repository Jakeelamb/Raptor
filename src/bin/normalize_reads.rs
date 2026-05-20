use raptor::pipeline::normalize::{
    normalize_single_with_config, NormalizeConfig, NormalizeSummary,
};
use std::io;

fn run(
    input_path: &str,
    output_prefix: &str,
    k: usize,
    target: u16,
    min_abund: u16,
) -> io::Result<NormalizeSummary> {
    normalize_single_with_config(
        input_path,
        &format!("{}_norm", output_prefix),
        NormalizeConfig {
            k,
            target_coverage: target,
            min_abundance: min_abund,
            ..NormalizeConfig::default()
        },
        true,
    )
}

#[inline]
fn percentage(kept: usize, total: usize) -> f64 {
    if total == 0 {
        0.0
    } else {
        (kept as f64 / total as f64) * 100.0
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: {} <input.fastq(.gz)> <output_prefix> [k-size] [target] [min_abund]",
            args[0]
        );
        std::process::exit(1);
    }

    let input_path = &args[1];
    let output_prefix = &args[2];
    let k = if args.len() > 3 {
        args[3].parse().unwrap_or(NormalizeConfig::default().k)
    } else {
        NormalizeConfig::default().k
    };
    let target = if args.len() > 4 {
        args[4]
            .parse()
            .unwrap_or(NormalizeConfig::default().target_coverage)
    } else {
        NormalizeConfig::default().target_coverage
    };
    let min_abund = if args.len() > 5 {
        args[5]
            .parse()
            .unwrap_or(NormalizeConfig::default().min_abundance)
    } else {
        NormalizeConfig::default().min_abundance
    };

    println!(
        "Running in normalization mode with parameters: k={}, target={}, min_abund={}",
        k, target, min_abund
    );

    match run(input_path, output_prefix, k, target, min_abund) {
        Ok(summary) => {
            println!(
                "Normalization complete. Kept {}/{} reads ({:.1}%)",
                summary.kept_reads,
                summary.total_reads,
                percentage(summary.kept_reads, summary.total_reads)
            );
        }
        Err(err) => {
            eprintln!("Normalization failed: {}", err);
            std::process::exit(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::run;
    use raptor::io::fastq::{stream_fastq_records_checked, try_open_fastq};
    use std::io::{self, Write};
    use tempfile::{NamedTempFile, TempDir};

    #[test]
    fn run_rejects_truncated_fastq_input() {
        let mut input = NamedTempFile::new().expect("create input fastq");
        writeln!(input, "@read_1").expect("write header");
        writeln!(input, "ACGT").expect("write sequence");
        writeln!(input, "+").expect("write plus line");
        input.flush().expect("flush input");

        let temp_dir = TempDir::new().expect("create temp dir");
        let output_prefix = temp_dir.path().join("normalized");
        let err = run(
            input.path().to_str().expect("utf8 input path"),
            output_prefix.to_str().expect("utf8 output prefix"),
            3,
            5,
            1,
        )
        .expect_err("truncated FASTQ must return an error");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn run_streams_output_without_in_memory_record_buffer() {
        let mut input = NamedTempFile::new().expect("create input fastq");
        writeln!(input, "@read_1").expect("write header");
        writeln!(input, "ACGTAC").expect("write sequence");
        writeln!(input, "+").expect("write plus line");
        writeln!(input, "IIIIII").expect("write quality");
        writeln!(input, "@read_2").expect("write header");
        writeln!(input, "TGCATG").expect("write sequence");
        writeln!(input, "+").expect("write plus line");
        writeln!(input, "IIIIII").expect("write quality");
        input.flush().expect("flush input");

        let temp_dir = TempDir::new().expect("create temp dir");
        let output_prefix = temp_dir.path().join("normalized");
        let summary = run(
            input.path().to_str().expect("utf8 input path"),
            output_prefix.to_str().expect("utf8 output prefix"),
            3,
            5,
            1,
        )
        .expect("normalization should succeed");

        assert_eq!(summary.total_reads, 2);
        assert!(summary.kept_reads <= summary.total_reads);

        let output_path = format!(
            "{}_norm.fastq.gz",
            output_prefix.to_str().expect("utf8 output prefix")
        );
        let reader = try_open_fastq(&output_path).expect("open normalized output");
        let mut observed = 0usize;
        for record in stream_fastq_records_checked(reader) {
            record.expect("valid normalized record");
            observed += 1;
        }
        assert_eq!(observed, summary.kept_reads);
    }
}
