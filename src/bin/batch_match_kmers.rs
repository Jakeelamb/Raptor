use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use clap::Parser;
use raptor::accel::simd::batch_match_kmers_simd;

/// Batch k-mer matching using SIMD and parallel processing acceleration
#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    /// Input file with query sequences (one per line)
    #[clap(short, long)]
    queries: PathBuf,
    
    /// Input file with target sequences to search in (one per line)
    #[clap(short, long)]
    targets: PathBuf,
    
    /// Output file for results
    #[clap(short, long)]
    output: PathBuf,
    
    /// Maximum mismatches allowed
    #[clap(short, long, default_value = "2")]
    mismatches: usize,
    
    /// Minimum required overlap length
    #[clap(short = 'n', long, default_value = "15")]
    min_overlap: usize,
}

fn main() -> io::Result<()> {
    let args = Args::parse();
    
    // Validate parameters
    if args.min_overlap < 5 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Minimum overlap length must be at least 5"
        ));
    }
    
    if args.mismatches > args.min_overlap / 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Maximum mismatches ({}) should not exceed half the minimum overlap length ({})",
                   args.mismatches, args.min_overlap / 2)
        ));
    }
    
    // Read queries
    let query_file = File::open(&args.queries).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("Failed to open query file '{}': {}", args.queries.display(), e)
        )
    })?;
    
    let reader = BufReader::new(query_file);
    let queries: Vec<String> = reader.lines()
        .filter_map(Result::ok)
        .collect();
        
    if queries.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Query file contains no valid sequences"
        ));
    }
        
    println!("Loaded {} query sequences", queries.len());
    
    // Read targets
    let target_file = File::open(&args.targets).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("Failed to open target file '{}': {}", args.targets.display(), e)
        )
    })?;
    
    let reader = BufReader::new(target_file);
    let targets: Vec<String> = reader.lines()
        .filter_map(Result::ok)
        .collect();
        
    if targets.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Target file contains no valid sequences"
        ));
    }
        
    println!("Loaded {} target sequences", targets.len());
    
    // Open output file
    let output_file = File::create(&args.output).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("Failed to create output file '{}': {}", args.output.display(), e)
        )
    })?;
    
    let mut writer = BufWriter::new(output_file);
    
    writeln!(writer, "query_idx,target_idx,overlap_size,shift,mismatches")?;
    
    // Process each query against all targets
    let mut total_matches = 0;
    let target_bytes: Vec<_> = targets.iter().map(|s| s.as_bytes()).collect();
    
    for (q_idx, query) in queries.iter().enumerate() {
        let query_bytes = query.as_bytes();
        
        // Use SIMD batch matching
        let matches = batch_match_kmers_simd(
            query_bytes,
            &target_bytes,
            args.min_overlap,
            args.mismatches
        );
        
        // Write results
        for (t_idx, (overlap, mismatch_count)) in matches.iter().enumerate() {
            if let Some((overlap_size, shift)) = overlap {
                writeln!(writer, "{},{},{},{},{}", q_idx, t_idx, overlap_size, shift, mismatch_count)?;
                total_matches += 1;
            }
        }
        
        // Progress update
        if (q_idx + 1) % 100 == 0 || q_idx + 1 == queries.len() {
            println!("Processed {}/{} queries, found {} matches so far", 
                q_idx + 1, queries.len(), total_matches);
        }
    }
    
    println!("Complete! Found {} total matches", total_matches);
    Ok(())
}