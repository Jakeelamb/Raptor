use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId, Throughput};
use raptor::accel::cpu_backend::CpuBackend;
use raptor::accel::backend::ComputeBackend;
use raptor::kmer::nthash::{nthash, NtHashIterator};
use raptor::kmer::kmer::KmerU64;
use ahash::AHashMap;
use rand::Rng;

/// Generate random DNA sequences for benchmarking
fn generate_sequences(num_seqs: usize, seq_len: usize) -> Vec<String> {
    let mut rng = rand::thread_rng();
    let bases = ['A', 'C', 'G', 'T'];

    (0..num_seqs)
        .map(|_| {
            (0..seq_len)
                .map(|_| bases[rng.gen_range(0..4)])
                .collect()
        })
        .collect()
}

/// Benchmark k-mer counting methods
fn bench_kmer_counting(c: &mut Criterion) {
    let mut group = c.benchmark_group("kmer_counting");

    // Test with different dataset sizes
    for num_seqs in [100, 1000, 10000] {
        let sequences = generate_sequences(num_seqs, 150); // 150bp reads
        let total_bases: usize = sequences.iter().map(|s| s.len()).sum();

        group.throughput(Throughput::Bytes(total_bases as u64));

        // u64-encoded counting (optimized)
        group.bench_with_input(
            BenchmarkId::new("u64_encoded", num_seqs),
            &sequences,
            |b, seqs| {
                let backend = CpuBackend::new();
                b.iter(|| {
                    black_box(backend.count_kmers_u64(seqs, 31))
                });
            },
        );

        // String-based counting (legacy)
        group.bench_with_input(
            BenchmarkId::new("string_based", num_seqs),
            &sequences,
            |b, seqs| {
                let backend = CpuBackend::new();
                b.iter(|| {
                    black_box(backend.count_kmers(seqs, 31))
                });
            },
        );

        // Bloom-filtered counting
        group.bench_with_input(
            BenchmarkId::new("bloom_filtered", num_seqs),
            &sequences,
            |b, seqs| {
                let backend = CpuBackend::new();
                b.iter(|| {
                    black_box(backend.count_kmers_u64_filtered(seqs, 31, 2))
                });
            },
        );
    }

    group.finish();
}

/// Benchmark ntHash vs naive hashing
fn bench_hashing(c: &mut Criterion) {
    let mut group = c.benchmark_group("hashing");

    let sequence = generate_sequences(1, 10000)[0].clone();
    let bytes = sequence.as_bytes();
    let k = 31;

    group.throughput(Throughput::Bytes(sequence.len() as u64));

    // ntHash rolling (O(1) per k-mer)
    group.bench_function("nthash_rolling", |b| {
        b.iter(|| {
            let count: usize = NtHashIterator::new(bytes, k).count();
            black_box(count)
        });
    });

    // KmerU64 sliding window
    group.bench_function("kmeru64_sliding", |b| {
        b.iter(|| {
            let mut count = 0usize;
            if bytes.len() >= k {
                if let Some(mut kmer) = KmerU64::from_slice(&bytes[0..k]) {
                    count += 1;
                    for i in k..bytes.len() {
                        if let Some(next) = kmer.extend(bytes[i]) {
                            kmer = next;
                            count += 1;
                        }
                    }
                }
            }
            black_box(count)
        });
    });

    // Naive recomputation (O(k) per k-mer)
    group.bench_function("naive_recompute", |b| {
        b.iter(|| {
            let mut count = 0usize;
            for i in 0..=bytes.len().saturating_sub(k) {
                if let Some(_) = nthash(&bytes[i..i+k]) {
                    count += 1;
                }
            }
            black_box(count)
        });
    });

    group.finish();
}

/// Benchmark adjacency building
fn bench_adjacency(c: &mut Criterion) {
    let mut group = c.benchmark_group("adjacency");

    for num_seqs in [100, 1000] {
        let sequences = generate_sequences(num_seqs, 150);
        let backend = CpuBackend::new();
        let kmer_counts = backend.count_kmers_u64(&sequences, 31);

        group.bench_with_input(
            BenchmarkId::new("build_adjacency_u64", num_seqs),
            &kmer_counts,
            |b, counts| {
                b.iter(|| {
                    black_box(backend.build_adjacency_u64(counts, 31))
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_kmer_counting, bench_hashing, bench_adjacency);
criterion_main!(benches);
