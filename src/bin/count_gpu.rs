#[cfg(feature = "gpu")]
use raptor::gpu::kmer_gpu::GpuKmerCounter;
#[cfg(feature = "gpu")]
use std::time::Instant;

fn main() {
    let input_path_str = "sample_large.fastq";
    let k = 25;

    let reader = raptor::io::fastq::open_fastq(input_path_str);
    let records: Vec<_> = raptor::io::fastq::read_fastq_records(reader).collect();
    let sequences: Vec<String> = records.iter().map(|r| r.sequence.clone()).collect();

    #[cfg(feature = "gpu")]
    {
        // Create counter with k-mer size, max reads, and hash table size (24 bits = 16M entries)
        match GpuKmerCounter::new(k, sequences.len(), 24) {
            Ok(counter) => {
                let start = Instant::now();
                match counter.count(&sequences) {
                    Ok(counts) => {
                        let elapsed = start.elapsed();
                        println!("Time to count k-mers on GPU: {:?}", elapsed);
                        println!("Found {} unique k-mers", counts.len());
                    }
                    Err(e) => {
                        eprintln!("Error counting k-mers: {}", e);
                    }
                }
            }
            Err(e) => {
                eprintln!("Error creating GPU counter: {}", e);
            }
        }
    }

    #[cfg(not(feature = "gpu"))]
    {
        let _ = (k, sequences); // Suppress unused warnings
        println!("GPU support is not enabled. Compile with '--features gpu' to enable GPU support.");
    }
}
