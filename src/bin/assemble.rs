use assembler::io::fastq::{open_fastq, read_fastq_records};
use assembler::io::fasta::FastaWriter;
use assembler::kmer::kmer::canonical_kmer;
use assembler::graph::assembler::{greedy_assembly};
use rayon::ThreadPoolBuilder;

use std::collections::HashMap;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <input.fastq(.gz)> <output.fasta(.gz)> [min_len]", args[0]);
        std::process::exit(1);
    }
    let threads_arg = args.iter().position(|x| x == "--threads")
    .and_then(|i| args.get(i + 1))
    .and_then(|s| s.parse::<usize>().ok());

    let num_threads = threads_arg.unwrap_or_else(num_cpus::get);

    ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .expect("Failed to build thread pool");

    let input_path = &args[1];
    let output_path = &args[2];
    let min_len = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(100);

    let k = 7;
    let reader = open_fastq(input_path);
    let records = read_fastq_records(reader);

    let mut kmer_counts = HashMap::new();
    for record in records {
        for i in 0..=record.sequence.len().saturating_sub(k) {
            if let Some(kmer) = canonical_kmer(&record.sequence[i..i + k]) {
                *kmer_counts.entry(kmer).or_insert(0) += 1;
            }
        }
    }

    let contigs = greedy_assembly(k, &kmer_counts, min_len);

    let mut writer = FastaWriter::new(output_path);
    for (i, contig) in contigs.iter().enumerate() {
        writer.write_contig(contig, i + 1);
    }

    println!("Wrote {} contigs to {}", contigs.len(), output_path);
}
