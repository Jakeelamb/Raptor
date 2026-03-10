use ahash::AHashMap;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use raptor::accel::backend::AdjacencyTableU64;
use raptor::graph::assembler::greedy_assembly_u64;
use raptor::kmer::kmer::encode_kmer;

fn index_to_kmer(mut idx: usize, k: usize) -> u64 {
    let mut chars = vec!['A'; k];
    for pos in (0..k).rev() {
        chars[pos] = match idx & 0b11 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        };
        idx >>= 2;
    }
    let seq: String = chars.into_iter().collect();
    encode_kmer(&seq).expect("generated k-mer is valid DNA")
}

fn build_linear_graph(k: usize, num_kmers: usize) -> (AHashMap<u64, u32>, AdjacencyTableU64) {
    let mut counts = AHashMap::with_capacity(num_kmers);
    let mut adjacency = AdjacencyTableU64::with_capacity(k as u8, num_kmers);
    let mut nodes = Vec::with_capacity(num_kmers);

    for i in 0..num_kmers {
        let kmer = index_to_kmer(i, k);
        nodes.push(kmer);
        counts.insert(kmer, (num_kmers - i) as u32 + 1);
    }

    for window in nodes.windows(2) {
        let from = window[0];
        let to = window[1];
        let count = *counts.get(&to).unwrap_or(&1);
        adjacency.add_edge(from, to, count);
    }

    (counts, adjacency)
}

fn bench_greedy_assembly(c: &mut Criterion) {
    let mut group = c.benchmark_group("greedy_assembly_u64");
    let k = 15;

    for num_kmers in [1_000usize, 10_000, 50_000] {
        let (counts, adjacency) = build_linear_graph(k, num_kmers);
        group.throughput(Throughput::Elements(num_kmers as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(num_kmers),
            &(counts, adjacency),
            |b, (counts, adjacency)| {
                b.iter(|| {
                    let contigs =
                        greedy_assembly_u64(k, black_box(counts), black_box(adjacency), k);
                    black_box(contigs.len());
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_greedy_assembly);
criterion_main!(benches);
