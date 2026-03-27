use criterion::{
    black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion, Throughput,
};
use memmap2::Mmap;
use raptor::io::fastq::{for_each_fastq_sequence_checked, stream_fastq_records_checked};
use raptor::kmer::disk_counting_v2::{DiskCounterConfig, DiskKmerCounterV2};
use raptor::kmer::kmer::KmerU64;
use raptor::pipeline::large_genome_assembler::LargeGenomeStageBenchFixture;
use rayon::prelude::*;
use std::fs::{self, OpenOptions};
use std::io::{BufWriter, Cursor, Write};
use std::path::Path;
use tempfile::TempDir;

fn deterministic_dna_sequence(seed: usize, len: usize) -> Vec<u8> {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    (0..len)
        .map(|idx| BASES[(seed.wrapping_mul(17).wrapping_add(idx.wrapping_mul(13))) & 0b11])
        .collect()
}

fn generate_sequences(num_sequences: usize, seq_len: usize) -> Vec<Vec<u8>> {
    (0..num_sequences)
        .map(|idx| deterministic_dna_sequence(idx, seq_len))
        .collect()
}

fn build_fastq_fixture(num_sequences: usize, seq_len: usize) -> Vec<u8> {
    let mut fastq = Vec::with_capacity(num_sequences * (seq_len * 2 + 32));
    for idx in 0..num_sequences {
        let sequence = deterministic_dna_sequence(idx, seq_len);
        fastq.extend_from_slice(format!("@read_{idx}\n").as_bytes());
        fastq.extend_from_slice(&sequence);
        fastq.extend_from_slice(b"\n+\n");
        fastq.extend(std::iter::repeat_n(b'I', seq_len));
        fastq.push(b'\n');
    }
    fastq
}

fn encode_base_2bit(base: u8) -> Option<u64> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn kmer_mask(k: usize) -> u64 {
    if k >= 32 {
        u64::MAX
    } else {
        (1u64 << (k * 2)) - 1
    }
}

fn legacy_reopen_distribution(
    sequences: &[Vec<u8>],
    k: usize,
    num_buckets: usize,
    batch_size: usize,
    temp_dir: &Path,
) -> std::io::Result<u64> {
    fs::create_dir_all(temp_dir)?;
    let use_mask_bucket = num_buckets.is_power_of_two();
    let bucket_mask = (num_buckets.saturating_sub(1)) as u64;
    let rolling_mask = kmer_mask(k);
    let mut total_kmers = 0u64;

    for chunk in sequences.chunks(batch_size) {
        let mut writers = (0..num_buckets)
            .map(|bucket_idx| {
                let path = temp_dir.join(format!("bucket_{bucket_idx:05}.bin"));
                let file = OpenOptions::new().create(true).append(true).open(path)?;
                Ok(BufWriter::with_capacity(1 << 20, file))
            })
            .collect::<std::io::Result<Vec<_>>>()?;

        for sequence in chunk {
            if sequence.len() < k {
                continue;
            }

            let mut rolling = 0u64;
            let mut valid_run = 0usize;
            for &base in sequence {
                if let Some(base_bits) = encode_base_2bit(base) {
                    rolling = ((rolling << 2) | base_bits) & rolling_mask;
                    valid_run += 1;
                    if valid_run >= k {
                        let canonical = KmerU64 {
                            encoded: rolling,
                            len: k as u8,
                        }
                        .canonical()
                        .encoded;
                        let bucket_id = if use_mask_bucket {
                            (canonical & bucket_mask) as usize
                        } else {
                            (canonical % num_buckets as u64) as usize
                        };
                        writers[bucket_id].write_all(&canonical.to_le_bytes())?;
                        total_kmers = total_kmers.saturating_add(1);
                    }
                } else {
                    rolling = 0;
                    valid_run = 0;
                }
            }
        }

        for writer in &mut writers {
            writer.flush()?;
        }
    }

    for bucket_idx in 0..num_buckets {
        let path = temp_dir.join(format!("bucket_{bucket_idx:05}.bin"));
        let _ = fs::remove_file(path);
    }
    Ok(total_kmers)
}

fn persistent_distribution(
    sequences: &[Vec<u8>],
    k: usize,
    num_buckets: usize,
    batch_size: usize,
    temp_dir: &Path,
) -> std::io::Result<u64> {
    let config = DiskCounterConfig {
        k,
        num_buckets,
        min_count: 1,
        temp_dir: temp_dir.to_path_buf(),
        write_buffer_size: 1 << 20,
    };
    let mut counter = DiskKmerCounterV2::new(config)?;
    for chunk in sequences.chunks(batch_size) {
        counter.distribute(chunk.iter().map(Vec::as_slice))?;
    }
    let total_kmers = counter.stats().0;
    counter.cleanup()?;
    Ok(total_kmers)
}

fn prepare_bucket_count_fixture(
    sequences: &[Vec<u8>],
    k: usize,
    num_buckets: usize,
    batch_size: usize,
) -> std::io::Result<(TempDir, u64)> {
    let temp_dir = TempDir::new()?;
    let config = DiskCounterConfig {
        k,
        num_buckets,
        min_count: 1,
        temp_dir: temp_dir.path().to_path_buf(),
        write_buffer_size: 1 << 20,
    };
    let mut counter = DiskKmerCounterV2::new(config)?;
    for chunk in sequences.chunks(batch_size) {
        counter.distribute(chunk.iter().map(Vec::as_slice))?;
    }
    let total_kmers = counter.stats().0;
    std::mem::forget(counter);
    Ok((temp_dir, total_kmers))
}

fn legacy_count_bucket(path: &Path, min_count: u32) -> std::io::Result<usize> {
    let file = std::fs::File::open(path)?;
    let file_size = file.metadata()?.len() as usize;
    if !file_size.is_multiple_of(8) {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "bucket file '{}' has invalid byte length {} (not divisible by 8)",
                path.display(),
                file_size
            ),
        ));
    }
    let num_kmers = file_size / 8;
    if num_kmers == 0 {
        return Ok(0);
    }

    let mmap = unsafe { Mmap::map(&file)? };
    let mut kmers = Vec::with_capacity(num_kmers);
    for chunk in mmap.chunks_exact(8) {
        let bytes: [u8; 8] = chunk.try_into().expect("chunk size must be 8");
        kmers.push(u64::from_le_bytes(bytes));
    }
    kmers.sort_unstable();

    let mut unique = 0usize;
    let mut idx = 0usize;
    while idx < kmers.len() {
        let current = kmers[idx];
        idx += 1;
        let mut count = 1u32;
        while idx < kmers.len() && kmers[idx] == current {
            count = count.saturating_add(1);
            idx += 1;
        }
        if count >= min_count {
            unique += 1;
        }
    }

    Ok(unique)
}

fn legacy_count_all_buckets(
    temp_dir: &Path,
    num_buckets: usize,
    min_count: u32,
) -> std::io::Result<usize> {
    let bucket_counts: Vec<usize> = (0..num_buckets)
        .into_par_iter()
        .map(|bucket_idx| {
            let path = temp_dir.join(format!("bucket_{bucket_idx:05}.bin"));
            legacy_count_bucket(&path, min_count)
        })
        .collect::<std::io::Result<Vec<_>>>()?;

    Ok(bucket_counts.into_iter().sum())
}

fn radix_count_bucket(path: &Path, min_count: u32) -> std::io::Result<usize> {
    let file = std::fs::File::open(path)?;
    let file_size = file.metadata()?.len() as usize;
    if !file_size.is_multiple_of(8) {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "bucket file '{}' has invalid byte length {} (not divisible by 8)",
                path.display(),
                file_size
            ),
        ));
    }

    let num_kmers = file_size / 8;
    if num_kmers == 0 {
        return Ok(0);
    }

    let mmap = unsafe { Mmap::map(&file)? };
    let mut kmers = Vec::with_capacity(num_kmers);
    for chunk in mmap.chunks_exact(8) {
        let bytes: [u8; 8] = chunk.try_into().expect("chunk size must be 8");
        kmers.push(u64::from_le_bytes(bytes));
    }

    let mut scratch = vec![0u64; kmers.len()];
    let mut counts = vec![0usize; 256];
    let mut src: &mut [u64] = kmers.as_mut_slice();
    let mut dst: &mut [u64] = scratch.as_mut_slice();

    for shift in (0..64).step_by(8) {
        counts.fill(0);

        for &value in src.iter() {
            counts[((value >> shift) & 0xff) as usize] += 1;
        }

        let mut offset = 0usize;
        for count in counts.iter_mut() {
            let current = *count;
            *count = offset;
            offset += current;
        }

        for &value in src.iter() {
            let bucket = ((value >> shift) & 0xff) as usize;
            dst[counts[bucket]] = value;
            counts[bucket] += 1;
        }

        std::mem::swap(&mut src, &mut dst);
    }

    let mut unique = 0usize;
    let mut idx = 0usize;
    while idx < kmers.len() {
        let current = kmers[idx];
        idx += 1;
        let mut count = 1u32;
        while idx < kmers.len() && kmers[idx] == current {
            count = count.saturating_add(1);
            idx += 1;
        }
        if count >= min_count {
            unique += 1;
        }
    }

    Ok(unique)
}

fn radix_count_all_buckets(
    temp_dir: &Path,
    num_buckets: usize,
    min_count: u32,
) -> std::io::Result<usize> {
    let bucket_counts: Vec<usize> = (0..num_buckets)
        .into_par_iter()
        .map(|bucket_idx| {
            let path = temp_dir.join(format!("bucket_{bucket_idx:05}.bin"));
            radix_count_bucket(&path, min_count)
        })
        .collect::<std::io::Result<Vec<_>>>()?;

    Ok(bucket_counts.into_iter().sum())
}

fn bench_fastq_parsing(c: &mut Criterion) {
    let mut group = c.benchmark_group("fastq_parsing");

    for read_count in [10_000usize, 50_000] {
        let fastq = build_fastq_fixture(read_count, 150);
        group.throughput(Throughput::Bytes(fastq.len() as u64));

        group.bench_with_input(
            BenchmarkId::new("full_record_checked", read_count),
            &fastq,
            |b, fastq| {
                b.iter(|| {
                    let mut total_bases = 0usize;
                    for record in stream_fastq_records_checked(Cursor::new(fastq.as_slice())) {
                        let record = record.expect("benchmark FASTQ should parse");
                        total_bases += record.sequence.len();
                    }
                    black_box(total_bases)
                });
            },
        );

        group.bench_with_input(
            BenchmarkId::new("sequence_only_checked", read_count),
            &fastq,
            |b, fastq| {
                b.iter(|| {
                    let mut total_bases = 0usize;
                    for_each_fastq_sequence_checked(Cursor::new(fastq.as_slice()), |sequence| {
                        total_bases += sequence.len();
                        Ok(())
                    })
                    .expect("benchmark FASTQ should parse");
                    black_box(total_bases)
                });
            },
        );
    }

    group.finish();
}

fn bench_disk_distribution(c: &mut Criterion) {
    let mut group = c.benchmark_group("disk_distribution");
    let k = 31;
    let num_buckets = 256;
    let batch_size = 1_000;

    for sequence_count in [10_000usize, 40_000] {
        let sequences = generate_sequences(sequence_count, 150);
        let total_bases: usize = sequences.iter().map(Vec::len).sum();
        group.throughput(Throughput::Bytes(total_bases as u64));

        group.bench_with_input(
            BenchmarkId::new("legacy_reopen_per_batch", sequence_count),
            &sequences,
            |b, sequences| {
                b.iter_batched(
                    || TempDir::new().expect("temp dir"),
                    |temp_dir| {
                        let total = legacy_reopen_distribution(
                            sequences,
                            k,
                            num_buckets,
                            batch_size,
                            temp_dir.path(),
                        )
                        .expect("legacy distribution should succeed");
                        black_box(total)
                    },
                    BatchSize::SmallInput,
                );
            },
        );

        group.bench_with_input(
            BenchmarkId::new("persistent_writers", sequence_count),
            &sequences,
            |b, sequences| {
                b.iter_batched(
                    || TempDir::new().expect("temp dir"),
                    |temp_dir| {
                        let total = persistent_distribution(
                            sequences,
                            k,
                            num_buckets,
                            batch_size,
                            temp_dir.path(),
                        )
                        .expect("persistent distribution should succeed");
                        black_box(total)
                    },
                    BatchSize::SmallInput,
                );
            },
        );
    }

    group.finish();
}

fn bench_disk_bucket_count(c: &mut Criterion) {
    let mut group = c.benchmark_group("disk_bucket_count");
    let k = 31;
    let num_buckets = 64;
    let batch_size = 1_000;

    for sequence_count in [10_000usize, 40_000] {
        let sequences = generate_sequences(sequence_count, 150);
        let total_kmers = sequence_count.saturating_mul(150usize.saturating_sub(k) + 1);
        group.throughput(Throughput::Elements(total_kmers as u64));

        group.bench_with_input(
            BenchmarkId::new("legacy_sort_unstable", sequence_count),
            &sequences,
            |b, sequences| {
                b.iter_batched(
                    || {
                        prepare_bucket_count_fixture(sequences, k, num_buckets, batch_size)
                            .expect("bucket fixture should build")
                    },
                    |(temp_dir, total_kmers)| {
                        let unique = legacy_count_all_buckets(temp_dir.path(), num_buckets, 1)
                            .expect("legacy count phase should succeed");
                        black_box((unique, total_kmers));
                        drop(temp_dir);
                    },
                    BatchSize::SmallInput,
                );
            },
        );

        group.bench_with_input(
            BenchmarkId::new("radix_sort_u64", sequence_count),
            &sequences,
            |b, sequences| {
                b.iter_batched(
                    || {
                        prepare_bucket_count_fixture(sequences, k, num_buckets, batch_size)
                            .expect("bucket fixture should build")
                    },
                    |(temp_dir, total_kmers)| {
                        let unique = radix_count_all_buckets(temp_dir.path(), num_buckets, 1)
                            .expect("radix count phase should succeed");
                        black_box((unique, total_kmers));
                        drop(temp_dir);
                    },
                    BatchSize::SmallInput,
                );
            },
        );
    }

    group.finish();
}

fn bench_graph_analysis(c: &mut Criterion) {
    let mut group = c.benchmark_group("graph_analysis");

    for component_count in [128usize, 512, 2048] {
        let fixture = LargeGenomeStageBenchFixture::synthetic_branching(11, component_count, 6, 2);
        group.throughput(Throughput::Elements(fixture.graph_node_count() as u64));

        group.bench_with_input(
            BenchmarkId::new("legacy_weighted_graph", component_count),
            &fixture,
            |b, fixture| {
                b.iter(|| {
                    let stats = fixture.run_weighted_graph_analysis_legacy();
                    black_box(stats)
                });
            },
        );

        group.bench_with_input(
            BenchmarkId::new("direct_scan", component_count),
            &fixture,
            |b, fixture| {
                b.iter(|| {
                    let stats = fixture.run_weighted_graph_analysis_direct();
                    black_box(stats)
                });
            },
        );
    }

    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(20);
    targets = bench_fastq_parsing, bench_disk_distribution, bench_disk_bucket_count, bench_graph_analysis
}
criterion_main!(benches);
