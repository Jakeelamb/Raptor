use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use raptor::graph::isoform_filter::kmer_jaccard_similarity;
use raptor::kmer::bloom::{BloomFilter, CountingBloomFilter};
use rand::Rng;

/// Generate random DNA sequences for benchmarking
fn generate_sequence(len: usize) -> String {
    let mut rng = rand::thread_rng();
    let bases = ['A', 'C', 'G', 'T'];
    (0..len).map(|_| bases[rng.gen_range(0..4)]).collect()
}

/// Generate similar sequences (for similarity testing)
fn generate_similar_sequences(len: usize, similarity: f64) -> (String, String) {
    let mut rng = rand::thread_rng();
    let bases = ['A', 'C', 'G', 'T'];

    let seq1: String = (0..len).map(|_| bases[rng.gen_range(0..4)]).collect();

    let seq2: String = seq1
        .chars()
        .map(|c| {
            if rng.gen::<f64>() < similarity {
                c
            } else {
                bases[rng.gen_range(0..4)]
            }
        })
        .collect();

    (seq1, seq2)
}

/// Benchmark similarity computation methods
fn bench_similarity(c: &mut Criterion) {
    let mut group = c.benchmark_group("similarity");

    // Test different sequence lengths
    for len in [100, 500, 1000] {
        let (seq1, seq2) = generate_similar_sequences(len, 0.9);

        // K-mer Jaccard similarity (optimized)
        group.bench_with_input(
            BenchmarkId::new("kmer_jaccard", len),
            &(&seq1, &seq2),
            |b, (s1, s2)| {
                b.iter(|| {
                    black_box(kmer_jaccard_similarity(s1, s2, 11))
                });
            },
        );
    }

    group.finish();
}

/// Benchmark Bloom filter operations
fn bench_bloom_filter(c: &mut Criterion) {
    let mut group = c.benchmark_group("bloom_filter");

    // Test different expected item counts
    for num_items in [10_000, 100_000, 1_000_000] {
        // Bloom filter creation
        group.bench_with_input(
            BenchmarkId::new("create", num_items),
            &num_items,
            |b, &n| {
                b.iter(|| {
                    black_box(BloomFilter::with_fp_rate(n, 0.01))
                });
            },
        );

        // Bloom filter insertion
        let mut bloom = BloomFilter::with_fp_rate(num_items, 0.01);
        group.bench_with_input(
            BenchmarkId::new("insert", num_items),
            &num_items,
            |b, &n| {
                b.iter(|| {
                    for i in 0..1000 {
                        bloom.insert(i as u64 * 12345);
                    }
                    black_box(())
                });
            },
        );

        // Bloom filter lookup
        for i in 0..10000 {
            bloom.insert(i as u64 * 12345);
        }
        group.bench_with_input(
            BenchmarkId::new("lookup", num_items),
            &num_items,
            |b, _| {
                b.iter(|| {
                    let mut found = 0;
                    for i in 0..1000 {
                        if bloom.may_contain(i as u64 * 12345) {
                            found += 1;
                        }
                    }
                    black_box(found)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark counting Bloom filter
fn bench_counting_bloom(c: &mut Criterion) {
    let mut group = c.benchmark_group("counting_bloom");

    for num_items in [10_000, 100_000] {
        let mut cbf = CountingBloomFilter::with_fp_rate(num_items, 0.01);

        // Insert and count
        group.bench_with_input(
            BenchmarkId::new("insert_count", num_items),
            &num_items,
            |b, _| {
                b.iter(|| {
                    for i in 0..1000 {
                        cbf.insert(i as u64 * 12345);
                    }
                    black_box(())
                });
            },
        );

        // Populate for threshold checking
        for i in 0..10000 {
            cbf.insert(i as u64 * 12345);
            cbf.insert(i as u64 * 12345); // Insert twice
        }

        // Check threshold
        group.bench_with_input(
            BenchmarkId::new("count_at_least", num_items),
            &num_items,
            |b, _| {
                b.iter(|| {
                    let mut above_threshold = 0;
                    for i in 0..1000 {
                        if cbf.count_at_least(i as u64 * 12345, 2) {
                            above_threshold += 1;
                        }
                    }
                    black_box(above_threshold)
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_similarity, bench_bloom_filter, bench_counting_bloom);
criterion_main!(benches);
