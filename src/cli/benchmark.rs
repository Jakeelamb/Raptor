use std::time::Instant;
use crate::io::fastq::{open_fastq, stream_fastq_records};
use crate::kmer::nthash::NtHashIterator;
use ahash::AHashMap;

/// Benchmark k-mer counting using ntHash for fast rolling updates.
pub fn benchmark_kmer_counting(input: &str, k: usize) {
    let start = Instant::now();
    let mut counts: AHashMap<u64, u32> = AHashMap::new();

    let reader = open_fastq(input);
    let records = stream_fastq_records(reader);

    println!("Starting k-mer counting benchmark with k={} (using ntHash)", k);

    for record in records {
        // Use ntHash for O(1) rolling hash per k-mer
        for (_, hash) in NtHashIterator::new(record.sequence.as_bytes(), k) {
            *counts.entry(hash).or_insert(0) += 1;
        }
    }

    let duration = start.elapsed();
    println!(
        "K-mer counting complete: {} unique kmers in {:.2?} (avg {:.2} k/sec)",
        counts.len(),
        duration,
        counts.len() as f64 / duration.as_secs_f64() / 1_000.0
    );

    // Show some basic stats about the k-mer distribution
    let mut total_kmers: u64 = 0;
    let mut max_count: u32 = 0;

    for &count in counts.values() {
        total_kmers += count as u64;
        if count > max_count {
            max_count = count;
        }
    }

    println!("Total k-mers processed: {}", total_kmers);
    println!("Max k-mer frequency: {}", max_count);
    println!("Average k-mer frequency: {:.2}", total_kmers as f64 / counts.len() as f64);
} 