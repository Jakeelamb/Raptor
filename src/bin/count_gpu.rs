use assembler::gpu::kmer_gpu::GpuKmerCounter;

fn main() {
    let input_path = std::env::args().nth(1).unwrap();
    let k = 25;

    let reader = assembler::io::fastq::open_fastq(&input_path);
    let records: Vec<_> = assembler::io::fastq::read_fastq_records(reader).collect();
    let sequences: Vec<String> = records.iter().map(|r| r.sequence.clone()).collect();

    let counter = GpuKmerCounter::new(k, sequences.len());
    let sketch = counter.count(&sequences, k);

    for (i, c) in sketch.iter().enumerate().filter(|(_, &c)| c > 0) {
        println!("Hash {} → {}", i, c);
    }
}
