use raptor::gpu::kmer_gpu::GpuKmerCounter;
use std::path::Path;
use std::time::Instant;

fn main() {
    let input_path_str = "sample_large.fastq";
    let k = 25;

    let reader = raptor::io::fastq::open_fastq(input_path_str);
    let records: Vec<_> = raptor::io::fastq::read_fastq_records(reader).collect();
    let sequences: Vec<String> = records.iter().map(|r| r.sequence.clone()).collect();

    let counter = GpuKmerCounter::new(k, sequences.len());
    
    let start = Instant::now();
    let _counts = counter.count(&sequences, k);
    let elapsed = start.elapsed();
    
    println!("Time to count k-mers on GPU: {:?}", elapsed);
}
