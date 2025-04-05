use crate::io::fastq::{open_fastq, read_fastq_records, read_paired_fastq_records, FastqWriter,
                             stream_fastq_records, stream_paired_fastq_records, FastqRecord};
use crate::kmer::cms::CountMinSketch;
use crate::kmer::normalize::{should_keep_read, should_keep_read_pair};
use crate::kmer::kmer::canonical_kmer;
use tracing::info;
use std::time::Instant;

/// Normalizes single-end reads using streaming API for better memory efficiency
pub fn normalize_single(input_path: &str, output_prefix: &str, _use_gpu: bool, streaming: bool) {
    info!("Normalizing single-end reads: {}", input_path);
    let start_time = Instant::now();
    
    // First pass: count k-mers using streaming
    let k = 25;
    let target_coverage = 50;
    let min_abundance = 2;
    
    info!("First pass: counting k-mers with k={}", k);
    let reader = open_fastq(input_path);
    
    let mut cms = CountMinSketch::new(4, 1 << 20);
    let mut total_records = 0;
    
    if streaming {
        info!("Using streaming mode for k-mer counting");
        // Stream records for k-mer counting
        for record in stream_fastq_records(reader) {
            total_records += 1;
            
            if record.sequence.len() >= k {
                for i in 0..=record.sequence.len() - k {
                    if let Some(kmer) = canonical_kmer(&record.sequence[i..i + k]) {
                        cms.insert(&kmer);
                    }
                }
            }
            
            // Progress update
            if total_records % 100_000 == 0 {
                info!("Processed {} records in first pass...", total_records);
            }
        }
    } else {
        info!("Using in-memory mode for k-mer counting");
        // Read all records into memory for k-mer counting
        for record in read_fastq_records(reader) {
            total_records += 1;
            
            if record.sequence.len() >= k {
                for i in 0..=record.sequence.len() - k {
                    if let Some(kmer) = canonical_kmer(&record.sequence[i..i + k]) {
                        cms.insert(&kmer);
                    }
                }
            }
            
            // Progress update
            if total_records % 100_000 == 0 {
                info!("Processed {} records in first pass...", total_records);
            }
        }
    }
    
    info!("Completed first pass: processed {} records in {:.2?}", 
         total_records, start_time.elapsed());
    
    // Second pass: filter reads using streaming
    info!("Second pass: filtering reads with target coverage {} and min abundance {}", 
          target_coverage, min_abundance);
    
    let reader = open_fastq(input_path);
    let mut writer = FastqWriter::new(&format!("{}.fastq.gz", output_prefix));
    let mut kept_count = 0;
    let second_pass_start = Instant::now();
    let mut current_record = 0;
    
    if streaming {
        info!("Using streaming mode for filtering");
        // Stream records for filtering
        for record in stream_fastq_records(reader) {
            current_record += 1;
            
            if should_keep_read(&record, &cms, k, target_coverage, min_abundance) {
                writer.write_record(&record).expect("Failed to write record");
                kept_count += 1;
            }
            
            // Progress update
            if current_record % 100_000 == 0 {
                info!("Processed {}/{} records in second pass...", current_record, total_records);
            }
        }
    } else {
        info!("Using in-memory mode for filtering");
        // Load all records into memory for filtering
        for record in read_fastq_records(reader) {
            current_record += 1;
            
            if should_keep_read(&record, &cms, k, target_coverage, min_abundance) {
                writer.write_record(&record).expect("Failed to write record");
                kept_count += 1;
            }
            
            // Progress update
            if current_record % 100_000 == 0 {
                info!("Processed {}/{} records in second pass...", current_record, total_records);
            }
        }
    }
    
    info!("Normalization complete: {}/{} reads kept ({:.1}%) in {:.2?}", 
        kept_count, total_records, (kept_count as f64 / total_records as f64) * 100.0,
        second_pass_start.elapsed());
}

/// Normalizes paired-end reads using streaming API for better memory efficiency
pub fn normalize_paired(input_r1: &str, input_r2: &str, output_prefix: &str, _use_gpu: bool, streaming: bool) {
    info!("Normalizing paired-end reads: {} and {}", input_r1, input_r2);
    let start_time = Instant::now();
    
    // First pass: count k-mers from both mates using streaming
    let k = 25;
    let target_coverage = 50;
    let min_abundance = 2;
    
    info!("First pass: counting k-mers with k={}", k);
    let reader1 = open_fastq(input_r1);
    let reader2 = open_fastq(input_r2);
    
    let mut cms = CountMinSketch::new(4, 1 << 20);
    let mut total_pairs = 0;
    
    if streaming {
        info!("Using streaming mode for paired k-mer counting");
        // Stream paired records for k-mer counting
        for (r1, r2) in stream_paired_fastq_records(reader1, reader2) {
            total_pairs += 1;
            
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
            
            // Progress update
            if total_pairs % 100_000 == 0 {
                info!("Processed {} read pairs in first pass...", total_pairs);
            }
        }
    } else {
        info!("Using in-memory mode for paired k-mer counting");
        // Load all paired records into memory for k-mer counting
        for (r1, r2) in read_paired_fastq_records(reader1, reader2) {
            total_pairs += 1;
            
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
            
            // Progress update
            if total_pairs % 100_000 == 0 {
                info!("Processed {} read pairs in first pass...", total_pairs);
            }
        }
    }
    
    info!("Completed first pass: processed {} read pairs in {:.2?}", 
          total_pairs, start_time.elapsed());
    
    // Second pass: filter read pairs using streaming
    info!("Second pass: filtering reads with target coverage {} and min abundance {}", 
          target_coverage, min_abundance);
    
    let reader1 = open_fastq(input_r1);
    let reader2 = open_fastq(input_r2);
    let mut writer1 = FastqWriter::new(&format!("{}_R1.fastq.gz", output_prefix));
    let mut writer2 = FastqWriter::new(&format!("{}_R2.fastq.gz", output_prefix));
    let mut kept_count = 0;
    let second_pass_start = Instant::now();
    let mut current_pair = 0;
    
    if streaming {
        info!("Using streaming mode for paired filtering");
        // Stream paired records for filtering
        for (r1, r2) in stream_paired_fastq_records(reader1, reader2) {
            current_pair += 1;
            
            if should_keep_read_pair(&r1, &r2, &cms, k, target_coverage, min_abundance) {
                writer1.write_record(&r1).expect("Failed to write R1");
                writer2.write_record(&r2).expect("Failed to write R2");
                kept_count += 1;
            }
            
            // Progress update
            if current_pair % 100_000 == 0 {
                info!("Processed {}/{} read pairs in second pass...", current_pair, total_pairs);
            }
        }
    } else {
        info!("Using in-memory mode for paired filtering");
        // Load all paired records into memory for filtering
        for (r1, r2) in read_paired_fastq_records(reader1, reader2) {
            current_pair += 1;
            
            if should_keep_read_pair(&r1, &r2, &cms, k, target_coverage, min_abundance) {
                writer1.write_record(&r1).expect("Failed to write R1");
                writer2.write_record(&r2).expect("Failed to write R2");
                kept_count += 1;
            }
            
            // Progress update
            if current_pair % 100_000 == 0 {
                info!("Processed {}/{} read pairs in second pass...", current_pair, total_pairs);
            }
        }
    }
    
    info!("Paired normalization complete: {}/{} pairs kept ({:.1}%) in {:.2?}", 
        kept_count, total_pairs, (kept_count as f64 / total_pairs as f64) * 100.0,
        second_pass_start.elapsed());
} 