use raptor::gpu::kmer_gpu::GpuKmerCounter;
use std::path::Path;
use std::time::Instant;

fn main() {
    let input_path = Path::new("sample_large.fastq");
    let k = 25;

    let reader = raptor::io::fastq::open_fastq(&input_path);
    let records: Vec<_> = raptor::io::fastq::read_fastq_records(reader).collect();
    let sequences: Vec<String> = records.iter().map(|r| r.sequence.clone()).collect();

    let counter = GpuKmerCounter::new(k, sequences.len());
    let sketch = counter.count(&sequences, k);

    for (i, c) in sketch.iter().enumerate().filter(|(_, &c)| c > 0) {
        println!("Hash {} → {}", i, c);
    }
}
