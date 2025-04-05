use raptor::io::fastq::{open_fastq, read_paired_fastq_records, FastqWriter};
use raptor::kmer::cms::CountMinSketch;
use raptor::kmer::normalize::{should_keep_read_pair, estimate_read_abundance};
use raptor::kmer::kmer::canonical_kmer;
use raptor::accel::gpu::kmer_gpu::GpuKmerCounter;

// Process reads in chunks to reduce memory usage
const CHUNK_SIZE: usize = 100_000; // Process this many read pairs at a time

// Default values for modern RNA-Seq datasets
const DEFAULT_MAX_READS: usize = 5_000_000;
const DEFAULT_COVERAGE_TARGET: usize = 500;
const DEFAULT_MIN_ABUNDANCE: usize = 1;
const DEFAULT_K: usize = 25;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    
    if args.len() < 4 || args.contains(&"--help".to_string()) || args.contains(&"-h".to_string()) {
        print_usage(&args[0]);
        std::process::exit(if args.contains(&"--help".to_string()) || args.contains(&"-h".to_string()) { 0 } else { 1 });
    }

    let input_r1 = &args[1];
    let input_r2 = &args[2];
    let output_prefix = &args[3];
    
    // Allow command-line overrides for parameters
    let k = if args.len() > 4 { args[4].parse().unwrap_or(DEFAULT_K) } else { DEFAULT_K };
    let target = if args.len() > 5 { args[5].parse().unwrap_or(DEFAULT_COVERAGE_TARGET) } else { DEFAULT_COVERAGE_TARGET };
    let min_abund = if args.len() > 6 { args[6].parse().unwrap_or(DEFAULT_MIN_ABUNDANCE) } else { DEFAULT_MIN_ABUNDANCE };
    let gpu_enabled = args.contains(&"--gpu".to_string());
    let streaming_mode = args.contains(&"--streaming".to_string()); // Fixed: don't default to true
    
    println!("Running paired-end normalization with parameters: k={}, target={}, min_abund={}", k, target, min_abund);
    
    let mut writer1 = FastqWriter::new(&format!("{}_R1.norm.fastq.gz", output_prefix));
    let mut writer2 = FastqWriter::new(&format!("{}_R2.norm.fastq.gz", output_prefix));
    let mut total_pairs = 0;
    let mut kept_pairs = 0;
    
    if gpu_enabled {
        // GPU mode - process in chunks to avoid memory issues
        println!("Running GPU-accelerated normalization for paired-end reads");
        println!("Memory efficient mode: Processing in chunks of {} read pairs", CHUNK_SIZE);
        
        let reader1 = open_fastq(input_r1);
        let reader2 = open_fastq(input_r2);
        let mut pair_iterator = read_paired_fastq_records(reader1, reader2);
        
        // Process in chunks
        loop {
            let mut chunk = Vec::with_capacity(CHUNK_SIZE);
            
            // Collect a chunk of read pairs
            for _ in 0..CHUNK_SIZE {
                if let Some(pair) = pair_iterator.next() {
                    chunk.push(pair);
                } else {
                    break;
                }
            }
            
            // If no more reads, we're done
            if chunk.is_empty() {
                break;
            }
            
            total_pairs += chunk.len();
            
            // Flatten sequences from both read pairs
            let sequences: Vec<String> = chunk.iter()
                .flat_map(|(r1, r2)| vec![r1.sequence.clone(), r2.sequence.clone()])
                .collect();
            
            // Initialize GPU counter with a limited batch size
            let counter = GpuKmerCounter::new(k, sequences.len());
            let gpu_counts = counter.count(&sequences, k);
            
            // Filter pairs based on GPU k-mer counts
            for (r1, r2) in &chunk {
                let est_r1 = estimate_read_abundance(&r1.sequence, k, &gpu_counts);
                let est_r2 = estimate_read_abundance(&r2.sequence, k, &gpu_counts);
                
                let keep = est_r1.median <= 50 && est_r1.min >= 2 &&
                           est_r2.median <= 50 && est_r2.min >= 2;
                
                if keep {
                    writer1.write_record(r1).expect("Failed to write R1 record");
                    writer2.write_record(r2).expect("Failed to write R2 record");
                    kept_pairs += 1;
                }
            }
            
            println!("Processed {} read pairs, kept {} so far", total_pairs, kept_pairs);
            
            // Clear chunk to free memory
            chunk.clear();
        }
        
        println!("Paired normalization complete. Kept {}/{} read pairs ({:.1}%)", 
            kept_pairs, total_pairs, (kept_pairs as f64 / total_pairs as f64) * 100.0);
    } else if streaming_mode {
        println!("Using streaming mode for memory efficiency");
        println!("Memory efficient mode: Processing in chunks of {} read pairs", CHUNK_SIZE);
        
        // Create a much smaller CMS to reduce memory usage
        let mut cms = CountMinSketch::new(2, 1 << 17); // 2 hash functions, even smaller table
        
        let reader1 = open_fastq(input_r1);
        let reader2 = open_fastq(input_r2);
        let mut pair_iterator = read_paired_fastq_records(reader1, reader2);
        
        // Process reads in chunks
        loop {
            let mut chunk = Vec::with_capacity(CHUNK_SIZE);
            
            // Collect a chunk of read pairs
            for _ in 0..CHUNK_SIZE {
                if let Some(pair) = pair_iterator.next() {
                    chunk.push(pair);
                } else {
                    break;
                }
            }
            
            // If no more reads, we're done
            if chunk.is_empty() {
                break;
            }
            
            // Count k-mers for this chunk
            for (r1, r2) in &chunk {
                // Process R1
                if r1.sequence.len() >= k {
                    for i in 0..=r1.sequence.len() - k {
                        if let Some(kmer) = canonical_kmer(&r1.sequence[i..i + k]) {
                            cms.insert(&kmer);
                        }
                    }
                }
                
                // Process R2
                if r2.sequence.len() >= k {
                    for i in 0..=r2.sequence.len() - k {
                        if let Some(kmer) = canonical_kmer(&r2.sequence[i..i + k]) {
                            cms.insert(&kmer);
                        }
                    }
                }
            }
            
            // Filter and write reads from this chunk
            for (r1, r2) in &chunk {
                total_pairs += 1;
                
                if should_keep_read_pair(r1, r2, &cms, k, target as u16, min_abund as u16) {
                    writer1.write_record(r1).expect("Failed to write R1 record");
                    writer2.write_record(r2).expect("Failed to write R2 record");
                    kept_pairs += 1;
                }
            }
            
            println!("Processed {} read pairs, kept {} so far", total_pairs, kept_pairs);
            
            // Clear chunk to free memory
            chunk.clear();
        }
        
        println!("Paired normalization complete. Kept {}/{} read pairs ({:.1}%)", 
            kept_pairs, total_pairs, (kept_pairs as f64 / total_pairs as f64) * 100.0);
    } else {
        // CPU mode without streaming
        println!("WARNING: Using non-streaming mode can consume significant memory.");
        println!("Consider using --streaming flag for large datasets.");
        
        let reader1 = open_fastq(input_r1);
        let reader2 = open_fastq(input_r2);
        let pairs = read_paired_fastq_records(reader1, reader2);
        let all_pairs: Vec<(_, _)> = pairs.collect();
        total_pairs = all_pairs.len();
        
        // CPU-based normalization
        let mut cms = CountMinSketch::new(4, 1 << 20);
        
        // First pass - count k-mers from both mates
        for (r1, r2) in &all_pairs {
            // Process R1
            if r1.sequence.len() >= k {
                for i in 0..=r1.sequence.len() - k {
                    if let Some(kmer) = canonical_kmer(&r1.sequence[i..i + k]) {
                        cms.insert(&kmer);
                    }
                }
            }
            
            // Process R2
            if r2.sequence.len() >= k {
                for i in 0..=r2.sequence.len() - k {
                    if let Some(kmer) = canonical_kmer(&r2.sequence[i..i + k]) {
                        cms.insert(&kmer);
                    }
                }
            }
        }
        
        // Second pass - filter read pairs
        for (r1, r2) in &all_pairs {
            if should_keep_read_pair(r1, r2, &cms, k, target as u16, min_abund as u16) {
                writer1.write_record(r1).expect("Failed to write R1 record");
                writer2.write_record(r2).expect("Failed to write R2 record");
                kept_pairs += 1;
            }
        }
        
        println!("Paired normalization complete. Kept {}/{} read pairs ({:.1}%)", 
            kept_pairs, total_pairs, (kept_pairs as f64 / total_pairs as f64) * 100.0);
    }
}

fn print_usage(program_name: &str) {
    eprintln!("Paired-end FASTQ Read Normalizer");
    eprintln!("--------------------------------");
    eprintln!("Usage: {} <R1.fastq(.gz)> <R2.fastq(.gz)> <output_prefix> [k-size] [target] [min_abund] [--gpu] [--streaming]", program_name);
    eprintln!();
    eprintln!("Parameters:");
    eprintln!("  <R1.fastq(.gz)>    Input FASTQ file for read 1 (can be gzipped)");
    eprintln!("  <R2.fastq(.gz)>    Input FASTQ file for read 2 (can be gzipped)");
    eprintln!("  <output_prefix>    Prefix for output files");
    eprintln!("  [k-size]           K-mer size (default: 25)");
    eprintln!("  [target]           Target coverage (default: 500)");
    eprintln!("  [min_abund]        Minimum k-mer abundance to keep (default: 1)");
    eprintln!("  [--gpu]            Use GPU acceleration for k-mer counting if available");
    eprintln!("  [--streaming]      Process reads in chunks to reduce memory usage");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  {} R1.fastq.gz R2.fastq.gz output", program_name);
    eprintln!("  {} R1.fastq.gz R2.fastq.gz output 21 50 2 --gpu", program_name);
    eprintln!();
    eprintln!("Output files:");
    eprintln!("  <output_prefix>_R1.norm.fastq.gz  - Normalized reads for mate 1");
    eprintln!("  <output_prefix>_R2.norm.fastq.gz  - Normalized reads for mate 2");
    eprintln!();
} 