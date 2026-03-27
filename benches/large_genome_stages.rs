use criterion::{
    black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion, Throughput,
};
use raptor::pipeline::large_genome_assembler::{
    LargeGenomeErrorCorrectionBenchFixture, LargeGenomeStageBenchFixture,
};

fn bench_branch_threading(c: &mut Criterion) {
    let mut group = c.benchmark_group("large_genome_branch_threading");

    for component_count in [32usize, 128, 512] {
        let fixture = LargeGenomeStageBenchFixture::synthetic_branching(11, component_count, 6, 2);
        group.throughput(Throughput::Bytes(fixture.total_branch_read_bases() as u64));
        group.bench_with_input(
            BenchmarkId::new("reread", component_count),
            &fixture,
            |b, fixture| {
                b.iter(|| {
                    let stats = fixture.run_branch_threading();
                    black_box((stats.supported_edges, stats.edge_observations))
                });
            },
        );
        group.bench_with_input(
            BenchmarkId::new("spool", component_count),
            &fixture,
            |b, fixture| {
                b.iter(|| {
                    let stats = fixture.run_branch_threading_from_spool();
                    black_box((stats.supported_edges, stats.edge_observations))
                });
            },
        );
    }

    group.finish();
}

fn bench_contig_extraction(c: &mut Criterion) {
    let mut group = c.benchmark_group("large_genome_contig_extraction");

    for component_count in [32usize, 128, 512] {
        let fixture = LargeGenomeStageBenchFixture::synthetic_branching(11, component_count, 6, 2);
        group.throughput(Throughput::Elements(fixture.graph_node_count() as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(component_count),
            &fixture,
            |b, fixture| {
                b.iter(|| {
                    let contigs = fixture.run_contig_extraction();
                    let total_bases: usize = contigs.iter().map(String::len).sum();
                    black_box((contigs.len(), total_bases))
                });
            },
        );
    }

    group.finish();
}

fn bench_error_correction(c: &mut Criterion) {
    let mut group = c.benchmark_group("large_genome_error_correction");

    for trusted_kmer_count in [4_096usize, 16_384, 65_536] {
        let fixture = LargeGenomeErrorCorrectionBenchFixture::synthetic_singleton_correction(
            31,
            trusted_kmer_count,
            2,
        );
        group.throughput(Throughput::Elements(fixture.total_kmer_count() as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(trusted_kmer_count),
            &fixture,
            |b, fixture| {
                b.iter_batched(
                    || fixture.cloned_kmer_counts(),
                    |counts| {
                        let stats = fixture.run_error_correction(counts);
                        black_box((stats.corrected_kmers, stats.remaining_kmers))
                    },
                    BatchSize::LargeInput,
                );
            },
        );
    }

    for trusted_kmer_count in [4_096usize, 16_384, 65_536] {
        let fixture = LargeGenomeErrorCorrectionBenchFixture::synthetic_rejection_heavy(
            31,
            trusted_kmer_count,
            2,
        );
        group.throughput(Throughput::Elements(fixture.total_kmer_count() as u64));
        group.bench_with_input(
            BenchmarkId::new("rejection-heavy", trusted_kmer_count),
            &fixture,
            |b, fixture| {
                b.iter_batched(
                    || fixture.cloned_kmer_counts(),
                    |counts| {
                        let stats = fixture.run_error_correction(counts);
                        black_box((stats.corrected_kmers, stats.remaining_kmers))
                    },
                    BatchSize::LargeInput,
                );
            },
        );
    }

    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = bench_branch_threading, bench_contig_extraction, bench_error_correction
}
criterion_main!(benches);
