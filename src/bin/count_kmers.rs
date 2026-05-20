use raptor::io::fastq::{stream_fastq_records, try_open_fastq};
use raptor::kmer::cms::CountMinSketch;
use raptor::kmer::kmer::KmerU64;
use raptor::kmer::nthash::NtHashIterator;
use std::collections::HashMap;
use std::time::Instant;

fn main() -> std::process::ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <input.fastq(.gz)> [k-mer length]", args[0]);
        return std::process::ExitCode::FAILURE;
    }

    let input_path = &args[1];
    let k = if args.len() > 2 {
        match args[2].parse::<usize>() {
            Ok(k) => k,
            Err(err) => {
                eprintln!("Invalid k-mer length '{}': {}", args[2], err);
                return std::process::ExitCode::FAILURE;
            }
        }
    } else {
        21 // Default k-mer size
    };

    println!("Processing {} with k-mer size: {}", input_path, k);
    let start = Instant::now();

    // Process the FASTQ file using streaming for memory efficiency
    let reader = match try_open_fastq(input_path) {
        Ok(reader) => reader,
        Err(err) => {
            eprintln!("Unable to open FASTQ file '{}': {}", input_path, err);
            return std::process::ExitCode::FAILURE;
        }
    };
    let records = stream_fastq_records(reader);

    // Set up k-mer counting
    let mut cms = CountMinSketch::new(4, 1 << 20);
    let mut total_kmers = 0usize;
    let mut total_reads = 0usize;
    let mut kmer_histogram: HashMap<u16, usize> = HashMap::new();
    // Track exact counts using u64 encoding for memory efficiency
    let mut exact_counts: HashMap<u64, u16> = HashMap::new();

    // Process each read
    for record in records {
        total_reads += 1;

        let bytes = record.sequence.as_bytes();
        if bytes.len() < k {
            continue;
        }

        // Process each k-mer using ntHash for fast hashing
        for (pos, hash) in NtHashIterator::new(bytes, k) {
            // Count in CMS using hash directly
            cms.insert_hash(hash);

            if let Some(kmer) = KmerU64::from_slice(&bytes[pos..pos + k]) {
                let canonical = kmer.canonical();
                let count = exact_counts.entry(canonical.encoded).or_insert(0);
                *count = count.saturating_add(1);
            }

            total_kmers += 1;
        }
    }

    // Record the elapsed time
    let elapsed = start.elapsed();

    // Extract and update the histogram from exact counts
    for &count in exact_counts.values() {
        *kmer_histogram.entry(count).or_insert(0) += 1;
    }

    // Print statistics
    println!("\nStatistics:");
    println!("  Total reads processed: {}", total_reads);
    println!("  Total k-mers processed: {}", total_kmers);
    println!("  Unique k-mers tracked: {}", exact_counts.len());
    println!("  Processing time: {:.2?}", elapsed);

    // Print histogram (top 10 frequencies)
    println!("\nK-mer frequency histogram (from exact counts):");
    let mut hist_entries: Vec<_> = kmer_histogram.into_iter().collect();
    hist_entries.sort_by_key(|&(count, _)| count);

    for (count, occurrences) in hist_entries.iter().take(10) {
        println!("  Count {}: {} k-mers", count, occurrences);
    }

    // Print top k-mers by frequency
    println!("\nTop 5 most frequent k-mers:");
    let mut top_kmers: Vec<_> = exact_counts.into_iter().collect();
    top_kmers.sort_by(|a, b| b.1.cmp(&a.1));

    for (kmer_encoded, count) in top_kmers.iter().take(5) {
        // Decode u64 back to string for display
        let kmer_str = raptor::kmer::kmer::decode_kmer(*kmer_encoded, k);
        println!("  {}: {} occurrences", kmer_str, count);
    }

    println!("\nComplete.");
    std::process::ExitCode::SUCCESS
}
