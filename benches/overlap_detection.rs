use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId, Throughput};
use raptor::accel::cpu_backend::CpuBackend;
use raptor::accel::backend::ComputeBackend;
use raptor::kmer::minimizer::{get_minimizers, MinimizerIndex};
use rand::Rng;

/// Generate random DNA contigs for benchmarking
fn generate_contigs(num_contigs: usize, contig_len: usize) -> Vec<String> {
    let mut rng = rand::thread_rng();
    let bases = ['A', 'C', 'G', 'T'];

    (0..num_contigs)
        .map(|_| {
            (0..contig_len)
                .map(|_| bases[rng.gen_range(0..4)])
                .collect()
        })
        .collect()
}

/// Generate contigs with deliberate overlaps for testing
fn generate_overlapping_contigs(num_contigs: usize, contig_len: usize, overlap_len: usize) -> Vec<String> {
    let mut rng = rand::thread_rng();
    let bases = ['A', 'C', 'G', 'T'];

    let mut contigs = Vec::with_capacity(num_contigs);

    // First contig is random
    let first: String = (0..contig_len)
        .map(|_| bases[rng.gen_range(0..4)])
        .collect();
    contigs.push(first);

    // Subsequent contigs share suffix/prefix overlaps
    for i in 1..num_contigs {
        let prev = &contigs[i - 1];
        let suffix = &prev[prev.len() - overlap_len..];

        // Start with overlap, then random bases
        let new_bases: String = (0..contig_len - overlap_len)
            .map(|_| bases[rng.gen_range(0..4)])
            .collect();

        contigs.push(format!("{}{}", suffix, new_bases));
    }

    contigs
}

/// Benchmark overlap detection methods
fn bench_overlap_detection(c: &mut Criterion) {
    let mut group = c.benchmark_group("overlap_detection");

    // Test with different numbers of contigs
    for num_contigs in [10, 50, 100] {
        let contigs = generate_overlapping_contigs(num_contigs, 500, 50);
        let total_bases: usize = contigs.iter().map(|s| s.len()).sum();

        group.throughput(Throughput::Bytes(total_bases as u64));

        // Minimizer-based overlap (current optimized implementation)
        group.bench_with_input(
            BenchmarkId::new("minimizer_index", num_contigs),
            &contigs,
            |b, ctgs| {
                let backend = CpuBackend::new();
                b.iter(|| {
                    black_box(backend.find_overlaps(ctgs, 20, 2))
                });
            },
        );
    }

    group.finish();
}

/// Benchmark minimizer extraction
fn bench_minimizer_extraction(c: &mut Criterion) {
    let mut group = c.benchmark_group("minimizer_extraction");

    let sequence = generate_contigs(1, 10000)[0].clone();
    let bytes = sequence.as_bytes();

    group.throughput(Throughput::Bytes(sequence.len() as u64));

    // Different window sizes
    for w in [3, 5, 10] {
        group.bench_with_input(
            BenchmarkId::new("window_size", w),
            &w,
            |b, &w| {
                b.iter(|| {
                    black_box(get_minimizers(bytes, 15, w))
                });
            },
        );
    }

    group.finish();
}

/// Benchmark minimizer index building and querying
fn bench_minimizer_index(c: &mut Criterion) {
    let mut group = c.benchmark_group("minimizer_index");

    for num_seqs in [100, 500, 1000] {
        let sequences: Vec<Vec<u8>> = generate_contigs(num_seqs, 500)
            .into_iter()
            .map(|s| s.into_bytes())
            .collect();

        let total_bases: usize = sequences.iter().map(|s| s.len()).sum();
        group.throughput(Throughput::Bytes(total_bases as u64));

        // Index building
        group.bench_with_input(
            BenchmarkId::new("build_index", num_seqs),
            &sequences,
            |b, seqs| {
                b.iter(|| {
                    black_box(MinimizerIndex::build(seqs, 15, 5))
                });
            },
        );

        // Index querying
        let index = MinimizerIndex::build(&sequences, 15, 5);
        let query = &sequences[0];

        group.bench_with_input(
            BenchmarkId::new("query_index", num_seqs),
            &(index, query),
            |b, (idx, q)| {
                b.iter(|| {
                    black_box(idx.find_overlap_candidates(q, 2))
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_overlap_detection, bench_minimizer_extraction, bench_minimizer_index);
criterion_main!(benches);
