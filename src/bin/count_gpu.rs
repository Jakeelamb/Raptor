#[cfg(feature = "gpu")]
use raptor::gpu::kmer_gpu::GpuKmerCounter;
#[cfg(feature = "gpu")]
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
        25
    };

    // Use streaming for memory efficiency
    let reader = match raptor::io::fastq::try_open_fastq(input_path) {
        Ok(reader) => reader,
        Err(err) => {
            eprintln!("Unable to open FASTQ file '{}': {}", input_path, err);
            return std::process::ExitCode::FAILURE;
        }
    };
    let records: Vec<_> = raptor::io::fastq::stream_fastq_records(reader).collect();
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
                        return std::process::ExitCode::FAILURE;
                    }
                }
            }
            Err(e) => {
                eprintln!("Error creating GPU counter: {}", e);
                return std::process::ExitCode::FAILURE;
            }
        }
        std::process::ExitCode::SUCCESS
    }

    #[cfg(not(feature = "gpu"))]
    {
        let _ = (k, sequences); // Suppress unused warnings
        println!(
            "GPU support is not enabled. Compile with '--features gpu' to enable GPU support."
        );
        std::process::ExitCode::FAILURE
    }
}
