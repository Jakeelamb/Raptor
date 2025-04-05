use raptor::io::fastq::{open_fastq, read_paired_fastq_records, FastqWriter};
use raptor::kmer::cms::CountMinSketch;
use raptor::kmer::normalize::{should_keep_read_pair, estimate_read_abundance, AbundanceStats};
use raptor::kmer::kmer::canonical_kmer;
use raptor::accel::gpu::kmer_gpu::GpuKmerCounter;

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
    let k = if args.len() > 4 { args[4].parse().unwrap_or(15) } else { 15 };
    let target = if args.len() > 5 { args[5].parse().unwrap_or(5) } else { 5 };
    let min_abund = if args.len() > 6 { args[6].parse().unwrap_or(1) } else { 1 };
    let gpu_enabled = args.contains(&"--gpu".to_string());
    
    println!("Running paired-end normalization with parameters: k={}, target={}, min_abund={}", k, target, min_abund);
    
    let reader1 = open_fastq(input_r1);
    let reader2 = open_fastq(input_r2);
    let pairs = read_paired_fastq_records(reader1, reader2);
    let mut all_pairs: Vec<(_, _)> = pairs.collect();
    
    if gpu_enabled {
        println!("Running GPU-accelerated normalization for paired-end reads");
        
        // Flatten sequences from both read pairs
        let sequences: Vec<String> = all_pairs.iter()
            .flat_map(|(r1, r2)| vec![r1.sequence.clone(), r2.sequence.clone()])
            .collect();
        
        // Initialize GPU counter
        let counter = GpuKmerCounter::new(k, sequences.len());
        let gpu_counts = counter.count(&sequences, k);
        
        // Second pass - filter pairs based on GPU k-mer counts
        let mut writer1 = FastqWriter::new(&format!("{}_R1.norm.fastq.gz", output_prefix));
        let mut writer2 = FastqWriter::new(&format!("{}_R2.norm.fastq.gz", output_prefix));
        let mut kept_count = 0;
        let total_count = all_pairs.len();
        
        for (r1, r2) in &all_pairs {
            let est_r1 = estimate_read_abundance(&r1.sequence, k, &gpu_counts);
            let est_r2 = estimate_read_abundance(&r2.sequence, k, &gpu_counts);
            
            let keep = est_r1.median <= 50 && est_r1.min >= 2 &&
                       est_r2.median <= 50 && est_r2.min >= 2;
            
            if keep {
                writer1.write_record(r1).expect("Failed to write R1 record");
                writer2.write_record(r2).expect("Failed to write R2 record");
                kept_count += 1;
            }
        }
        
        println!("Paired normalization complete. Kept {}/{} read pairs ({:.1}%)", 
            kept_count, total_count, (kept_count as f64 / total_count as f64) * 100.0);
    } else {
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
        let mut writer1 = FastqWriter::new(&format!("{}_R1.norm.fastq.gz", output_prefix));
        let mut writer2 = FastqWriter::new(&format!("{}_R2.norm.fastq.gz", output_prefix));
        let mut kept_count = 0;
        let total_count = all_pairs.len();
        
        for (r1, r2) in &all_pairs {
            if should_keep_read_pair(r1, r2, &cms, k, target as u16, min_abund as u16) {
                writer1.write_record(r1).expect("Failed to write R1 record");
                writer2.write_record(r2).expect("Failed to write R2 record");
                kept_count += 1;
            }
        }
        
        println!("Paired normalization complete. Kept {}/{} read pairs ({:.1}%)", 
            kept_count, total_count, (kept_count as f64 / total_count as f64) * 100.0);
    }
}

fn print_usage(program_name: &str) {
    eprintln!("Paired-end FASTQ Read Normalizer");
    eprintln!("--------------------------------");
    eprintln!("Usage: {} <R1.fastq(.gz)> <R2.fastq(.gz)> <output_prefix> [k-size] [target] [min_abund] [--gpu]", program_name);
    eprintln!();
    eprintln!("Parameters:");
    eprintln!("  <R1.fastq(.gz)>    Input FASTQ file for read 1 (can be gzipped)");
    eprintln!("  <R2.fastq(.gz)>    Input FASTQ file for read 2 (can be gzipped)");
    eprintln!("  <output_prefix>    Prefix for output files");
    eprintln!("  [k-size]           K-mer size (default: 15)");
    eprintln!("  [target]           Target coverage (default: 5)");
    eprintln!("  [min_abund]        Minimum k-mer abundance to keep (default: 1)");
    eprintln!("  [--gpu]            Use GPU acceleration for k-mer counting if available");
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