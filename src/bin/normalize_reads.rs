use raptor::io::fastq::{open_fastq, stream_fastq_records, FastqWriter};
use raptor::kmer::cms::CountMinSketch;
use raptor::kmer::normalize::should_keep_read;
#[allow(deprecated)]
use raptor::kmer::kmer::canonical_kmer;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <input.fastq(.gz)> <output_prefix> [k-size] [target] [min_abund]", args[0]);
        std::process::exit(1);
    }

    let input_path = &args[1];
    let output_prefix = &args[2];

    // Allow command-line overrides for parameters
    let k = if args.len() > 3 { args[3].parse().unwrap_or(15) } else { 15 }; // Shorter k size
    let target = if args.len() > 4 { args[4].parse().unwrap_or(5) } else { 5 }; // Lower target coverage
    let min_abund = if args.len() > 5 { args[5].parse().unwrap_or(1) } else { 1 }; // Accept even low-abundance k-mers

    println!("Running in normalization mode with parameters: k={}, target={}, min_abund={}", k, target, min_abund);

    // First pass: count k-mers using streaming for memory efficiency
    let reader = open_fastq(input_path);
    let records = stream_fastq_records(reader);
    let mut cms = CountMinSketch::new(4, 1 << 20);
    let mut all_records = vec![];
    
    for record in records {
        // Only consider k-mers if the sequence is long enough
        if record.sequence.len() >= k {
            for i in 0..=record.sequence.len() - k {
                if let Some(kmer) = canonical_kmer(&record.sequence[i..i + k]) {
                    cms.insert(&kmer);
                }
            }
        }
        all_records.push(record);
    }
    
    // Second pass: filter reads
    let mut writer = FastqWriter::new(&format!("{}_norm.fastq.gz", output_prefix));
    let mut kept_count = 0;
    let total_count = all_records.len();
    
    for record in &all_records {
        if should_keep_read(record, &cms, k, target, min_abund) {
            writer.write_record(record).expect("Failed to write record");
            kept_count += 1;
        }
    }
    
    println!("Normalization complete. Kept {}/{} reads ({:.1}%)", 
        kept_count, total_count, (kept_count as f64 / total_count as f64) * 100.0);
}
