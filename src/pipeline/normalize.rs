use crate::io::fastq::{
    stream_fastq_records_checked, stream_paired_fastq_records_checked, try_open_fastq, FastqWriter,
};
use crate::kmer::cms::CountMinSketch;
use crate::kmer::normalize::{should_keep_read_pair_with_scratch, should_keep_read_with_scratch};
use crate::kmer::nthash::NtHashIterator;
use serde::Serialize;
use std::io;
use std::time::Instant;
use tracing::info;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NormalizeConfig {
    pub k: usize,
    pub target_coverage: u16,
    pub min_abundance: u16,
    pub max_reads: Option<usize>,
    pub use_gpu: bool,
}

impl Default for NormalizeConfig {
    fn default() -> Self {
        Self {
            k: 25,
            target_coverage: 50,
            min_abundance: 2,
            max_reads: None,
            use_gpu: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub struct NormalizeSummary {
    pub total_reads: usize,
    pub kept_reads: usize,
}

#[inline]
fn should_stop_at_limit(processed: usize, max_reads: Option<usize>) -> bool {
    max_reads.is_some_and(|limit| processed >= limit)
}

/// Normalizes single-end reads using streaming API for better memory efficiency
pub fn normalize_single(
    input_path: &str,
    output_prefix: &str,
    _use_gpu: bool,
    streaming: bool,
) -> io::Result<()> {
    normalize_single_with_config(
        input_path,
        output_prefix,
        NormalizeConfig {
            use_gpu: _use_gpu,
            ..NormalizeConfig::default()
        },
        streaming,
    )
    .map(|_| ())
}

pub fn normalize_single_with_config(
    input_path: &str,
    output_prefix: &str,
    config: NormalizeConfig,
    streaming: bool,
) -> io::Result<NormalizeSummary> {
    info!("Normalizing single-end reads: {}", input_path);
    let start_time = Instant::now();

    let k = config.k;
    let target_coverage = config.target_coverage;
    let min_abundance = config.min_abundance;
    if config.use_gpu {
        info!("GPU normalization requested; current normalization path uses CPU ntHash/CMS");
    }

    info!("First pass: counting k-mers with k={}", k);
    let reader = try_open_fastq(input_path)?;

    let mut cms = CountMinSketch::new(4, 1 << 20);
    let mut total_records = 0;

    // Use streaming for k-mer counting (always more efficient)
    info!("Using streaming mode for k-mer counting with ntHash");
    for record in stream_fastq_records_checked(reader) {
        if should_stop_at_limit(total_records, config.max_reads) {
            break;
        }
        let record = record?;
        total_records += 1;

        // Use ntHash for O(1) rolling hash per k-mer
        for (_, hash) in NtHashIterator::new(record.sequence.as_bytes(), k) {
            cms.insert_hash(hash);
        }

        // Progress update
        if total_records % 100_000 == 0 {
            info!("Processed {} records in first pass...", total_records);
        }
    }
    let _ = streaming; // Suppress unused variable warning

    info!(
        "Completed first pass: processed {} records in {:.2?}",
        total_records,
        start_time.elapsed()
    );

    // Second pass: filter reads using streaming
    info!(
        "Second pass: filtering reads with target coverage {} and min abundance {}",
        target_coverage, min_abundance
    );

    let reader = try_open_fastq(input_path)?;
    let mut writer = FastqWriter::try_new(&format!("{}.fastq.gz", output_prefix))?;
    let mut kept_count = 0;
    let second_pass_start = Instant::now();
    let mut current_record = 0;
    let mut abundance_scratch = Vec::new();

    // Stream records for filtering
    for record in stream_fastq_records_checked(reader) {
        if should_stop_at_limit(current_record, config.max_reads) {
            break;
        }
        let record = record?;
        current_record += 1;

        if should_keep_read_with_scratch(
            &record,
            &cms,
            k,
            target_coverage,
            min_abundance,
            &mut abundance_scratch,
        ) {
            writer.write_record(&record)?;
            kept_count += 1;
        }

        // Progress update
        if current_record % 100_000 == 0 {
            info!(
                "Processed {}/{} records in second pass...",
                current_record, total_records
            );
        }
    }

    info!(
        "Normalization complete: {}/{} reads kept ({:.1}%) in {:.2?}",
        kept_count,
        total_records,
        percentage(kept_count, total_records),
        second_pass_start.elapsed()
    );
    Ok(NormalizeSummary {
        total_reads: total_records,
        kept_reads: kept_count,
    })
}

/// Normalizes paired-end reads using streaming API for better memory efficiency
pub fn normalize_paired(
    input_r1: &str,
    input_r2: &str,
    output_prefix: &str,
    _use_gpu: bool,
    streaming: bool,
) -> io::Result<()> {
    normalize_paired_with_config(
        input_r1,
        input_r2,
        output_prefix,
        NormalizeConfig {
            use_gpu: _use_gpu,
            ..NormalizeConfig::default()
        },
        streaming,
    )
    .map(|_| ())
}

pub fn normalize_paired_with_config(
    input_r1: &str,
    input_r2: &str,
    output_prefix: &str,
    config: NormalizeConfig,
    streaming: bool,
) -> io::Result<NormalizeSummary> {
    info!(
        "Normalizing paired-end reads: {} and {}",
        input_r1, input_r2
    );
    let start_time = Instant::now();

    let k = config.k;
    let target_coverage = config.target_coverage;
    let min_abundance = config.min_abundance;
    if config.use_gpu {
        info!("GPU normalization requested; current normalization path uses CPU ntHash/CMS");
    }

    info!("First pass: counting k-mers with k={}", k);
    let reader1 = try_open_fastq(input_r1)?;
    let reader2 = try_open_fastq(input_r2)?;

    let mut cms = CountMinSketch::new(4, 1 << 20);
    let mut total_pairs = 0;

    // Use streaming for k-mer counting (always more efficient)
    info!("Using streaming mode for paired k-mer counting with ntHash");
    for pair in stream_paired_fastq_records_checked(reader1, reader2) {
        if should_stop_at_limit(total_pairs, config.max_reads) {
            break;
        }
        let (r1, r2) = pair?;
        total_pairs += 1;

        // Process R1 with ntHash
        for (_, hash) in NtHashIterator::new(r1.sequence.as_bytes(), k) {
            cms.insert_hash(hash);
        }

        // Process R2 with ntHash
        for (_, hash) in NtHashIterator::new(r2.sequence.as_bytes(), k) {
            cms.insert_hash(hash);
        }

        // Progress update
        if total_pairs % 100_000 == 0 {
            info!("Processed {} read pairs in first pass...", total_pairs);
        }
    }
    let _ = streaming; // Suppress unused variable warning

    info!(
        "Completed first pass: processed {} read pairs in {:.2?}",
        total_pairs,
        start_time.elapsed()
    );

    // Second pass: filter read pairs using streaming
    info!(
        "Second pass: filtering reads with target coverage {} and min abundance {}",
        target_coverage, min_abundance
    );

    let reader1 = try_open_fastq(input_r1)?;
    let reader2 = try_open_fastq(input_r2)?;
    let mut writer1 = FastqWriter::try_new(&format!("{}_R1.fastq.gz", output_prefix))?;
    let mut writer2 = FastqWriter::try_new(&format!("{}_R2.fastq.gz", output_prefix))?;
    let mut kept_count = 0;
    let second_pass_start = Instant::now();
    let mut current_pair = 0;
    let mut scratch_r1 = Vec::new();
    let mut scratch_r2 = Vec::new();

    // Stream paired records for filtering
    for pair in stream_paired_fastq_records_checked(reader1, reader2) {
        if should_stop_at_limit(current_pair, config.max_reads) {
            break;
        }
        let (r1, r2) = pair?;
        current_pair += 1;

        if should_keep_read_pair_with_scratch(
            &r1,
            &r2,
            &cms,
            k,
            target_coverage,
            min_abundance,
            &mut scratch_r1,
            &mut scratch_r2,
        ) {
            writer1.write_record(&r1)?;
            writer2.write_record(&r2)?;
            kept_count += 1;
        }

        // Progress update
        if current_pair % 100_000 == 0 {
            info!(
                "Processed {}/{} read pairs in second pass...",
                current_pair, total_pairs
            );
        }
    }

    info!(
        "Paired normalization complete: {}/{} pairs kept ({:.1}%) in {:.2?}",
        kept_count,
        total_pairs,
        percentage(kept_count, total_pairs),
        second_pass_start.elapsed()
    );
    Ok(NormalizeSummary {
        total_reads: total_pairs,
        kept_reads: kept_count,
    })
}

#[inline]
fn percentage(kept: usize, total: usize) -> f64 {
    if total == 0 {
        0.0
    } else {
        (kept as f64 / total as f64) * 100.0
    }
}

#[cfg(test)]
mod tests {
    use super::{
        normalize_paired, normalize_paired_with_config, normalize_single,
        normalize_single_with_config, NormalizeConfig,
    };
    use crate::io::fastq::{
        stream_fastq_records_checked, stream_paired_fastq_records_checked, try_open_fastq,
    };
    use std::fs;
    use std::io::{self, Write};
    use tempfile::TempDir;

    #[test]
    fn normalize_single_rejects_truncated_fastq() {
        let temp_dir = TempDir::new().expect("temp dir");
        let input = temp_dir.path().join("input.fastq");
        let mut file = fs::File::create(&input).expect("create fastq");
        writeln!(file, "@r1").expect("header");
        writeln!(file, "ACGT").expect("sequence");
        writeln!(file, "+").expect("plus");
        file.flush().expect("flush");

        let output_prefix = temp_dir.path().join("normalized");
        let err = normalize_single(
            input.to_str().expect("utf8 path"),
            output_prefix.to_str().expect("utf8 path"),
            false,
            true,
        )
        .expect_err("expected truncated input error");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn normalize_paired_rejects_mismatched_record_counts() {
        let temp_dir = TempDir::new().expect("temp dir");
        let r1 = temp_dir.path().join("r1.fastq");
        let r2 = temp_dir.path().join("r2.fastq");

        fs::write(&r1, "@r1/1\nACGT\n+\nIIII\n@r2/1\nTGCA\n+\nIIII\n").expect("write r1");
        fs::write(&r2, "@r1/2\nACGT\n+\nIIII\n").expect("write r2");

        let output_prefix = temp_dir.path().join("normalized");
        let err = normalize_paired(
            r1.to_str().expect("utf8 path"),
            r2.to_str().expect("utf8 path"),
            output_prefix.to_str().expect("utf8 path"),
            false,
            true,
        )
        .expect_err("expected mismatched pair error");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("different record counts"));
    }

    #[test]
    fn normalize_single_handles_empty_input() {
        let temp_dir = TempDir::new().expect("temp dir");
        let input = temp_dir.path().join("empty.fastq");
        fs::File::create(&input).expect("create empty input");

        let output_prefix = temp_dir.path().join("normalized");
        normalize_single(
            input.to_str().expect("utf8 path"),
            output_prefix.to_str().expect("utf8 path"),
            false,
            true,
        )
        .expect("empty input should not fail");

        let output = temp_dir.path().join("normalized.fastq.gz");
        assert!(output.exists());
    }

    #[test]
    fn normalize_single_honors_max_reads_config() {
        let temp_dir = TempDir::new().expect("temp dir");
        let input = temp_dir.path().join("input.fastq");
        fs::write(
            &input,
            "@r1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n@r2\nACGTACGTACGT\n+\nIIIIIIIIIIII\n",
        )
        .expect("write input");

        let output_prefix = temp_dir.path().join("normalized");
        let summary = normalize_single_with_config(
            input.to_str().expect("utf8 path"),
            output_prefix.to_str().expect("utf8 path"),
            NormalizeConfig {
                k: 4,
                target_coverage: u16::MAX,
                min_abundance: 1,
                max_reads: Some(1),
                use_gpu: false,
            },
            true,
        )
        .expect("normalization should succeed");

        assert_eq!(summary.total_reads, 1);
        let output = temp_dir.path().join("normalized.fastq.gz");
        let reader = try_open_fastq(output.to_str().expect("utf8 path")).expect("open output");
        let observed = stream_fastq_records_checked(reader)
            .collect::<io::Result<Vec<_>>>()
            .expect("valid output");
        assert_eq!(observed.len(), summary.kept_reads);
        assert!(observed.len() <= 1);
    }

    #[test]
    fn normalize_paired_honors_max_reads_config() {
        let temp_dir = TempDir::new().expect("temp dir");
        let r1 = temp_dir.path().join("r1.fastq");
        let r2 = temp_dir.path().join("r2.fastq");

        fs::write(
            &r1,
            "@r1/1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n@r2/1\nTGCATGCATGCA\n+\nIIIIIIIIIIII\n",
        )
        .expect("write r1");
        fs::write(
            &r2,
            "@r1/2\nACGTACGTACGT\n+\nIIIIIIIIIIII\n@r2/2\nTGCATGCATGCA\n+\nIIIIIIIIIIII\n",
        )
        .expect("write r2");

        let output_prefix = temp_dir.path().join("normalized");
        let summary = normalize_paired_with_config(
            r1.to_str().expect("utf8 path"),
            r2.to_str().expect("utf8 path"),
            output_prefix.to_str().expect("utf8 path"),
            NormalizeConfig {
                k: 4,
                target_coverage: u16::MAX,
                min_abundance: 1,
                max_reads: Some(1),
                use_gpu: false,
            },
            true,
        )
        .expect("normalization should succeed");

        assert_eq!(summary.total_reads, 1);
        let output_r1 = temp_dir.path().join("normalized_R1.fastq.gz");
        let output_r2 = temp_dir.path().join("normalized_R2.fastq.gz");
        let reader1 = try_open_fastq(output_r1.to_str().expect("utf8 path")).expect("open r1");
        let reader2 = try_open_fastq(output_r2.to_str().expect("utf8 path")).expect("open r2");
        let observed = stream_paired_fastq_records_checked(reader1, reader2)
            .collect::<io::Result<Vec<_>>>()
            .expect("valid output");
        assert_eq!(observed.len(), summary.kept_reads);
        assert!(observed.len() <= 1);
    }
}
