use std::time::Instant;
use crate::io::fastq::{open_fastq, read_fastq_records};
use crate::kmer::kmer::canonical_kmer;
use std::collections::HashMap;

pub fn benchmark_kmer_counting(input: &str, k: usize) {
    let start = Instant::now();
    let mut counts = HashMap::new();

    let reader = open_fastq(input);
    let records = read_fastq_records(reader);

    println!("Starting k-mer counting benchmark with k={}", k);
    
    for record in records {
        for i in 0..=record.sequence.len().saturating_sub(k) {
            if let Some(kmer) = canonical_kmer(&record.sequence[i..i + k]) {
                *counts.entry(kmer).or_insert(0) += 1;
            }
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
    let mut total_kmers = 0;
    let mut max_count = 0;
    
    for &count in counts.values() {
        total_kmers += count;
        if count > max_count {
            max_count = count;
        }
    }
    
    println!("Total k-mers processed: {}", total_kmers);
    println!("Max k-mer frequency: {}", max_count);
    println!("Average k-mer frequency: {:.2}", total_kmers as f64 / counts.len() as f64);
} 