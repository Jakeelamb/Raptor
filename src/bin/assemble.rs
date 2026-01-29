use raptor::io::fastq::{open_fastq, read_fastq_records, FastqRecord};
use raptor::io::fasta::FastaWriter;
#[allow(deprecated)]
use raptor::kmer::kmer::canonical_kmer;
use raptor::kmer::kmer::encode_kmer;
use raptor::graph::assembler::Contig;
use rayon::ThreadPoolBuilder;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::fs::OpenOptions;
use std::io::Write;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} input.fastq output.fasta [--dev-mode]", args[0]);
        std::process::exit(1);
    }
    
    let input_path = &args[1];
    let output_path = &args[2];
    
    // Default parameters - modified for better results with low-coverage data
    let k = 25;  // K-mer size set to 25
    let min_coverage = 2;  // Lower coverage threshold
    let min_length = 50;   // Shorter minimum contig length
    let threads = num_cpus::get();
    let dev_mode = args.contains(&"--dev-mode".to_string());
    
    if dev_mode {
        println!("DEV MODE ENABLED: Will output additional debug information and skip filtering thresholds");
    }
    
    println!("Assembling with k={}, min_coverage={}, threads={}, min_length={}", 
             k, min_coverage, threads, min_length);
    
    // Initialize thread pool
    ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    
    // Read input FASTQ file
    println!("Reading input file {}...", input_path);
    let reader = open_fastq(input_path);
    let records: Vec<FastqRecord> = read_fastq_records(reader).collect();
    println!("Read {} records", records.len());
    
    // Build k-mer map from reads
    println!("Building k-mer map...");
    let mut kmer_map = HashMap::new();
    for record in &records {
        let sequence = &record.sequence;
        for i in 0..=sequence.len().saturating_sub(k) {
            if i + k <= sequence.len() {
                let kmer = &sequence[i..i+k];
                if let Some(canonical) = canonical_kmer(kmer) {
                    *kmer_map.entry(canonical).or_insert(0u32) += 1;
                }
            }
        }
    }
    
    // Filter low-coverage k-mers
    println!("Filtering k-mers with coverage < {}...", min_coverage);
    
    let original_kmer_count = kmer_map.len();
    
    if dev_mode {
        let debug_dir = Path::new(output_path).parent().unwrap_or(Path::new("."));
        let debug_file = debug_dir.join("debug_kmers.tsv");
        
        let mut file = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(debug_file)
            .unwrap_or_else(|e| {
                eprintln!("Warning: Could not create debug file: {}", e);
                OpenOptions::new().write(true).open("/dev/null").unwrap()
            });
            
        // Write all k-mers with their counts for analysis
        writeln!(&mut file, "kmer\tcount").unwrap();
        for (kmer, count) in &kmer_map {
            writeln!(&mut file, "{}\t{}", kmer, count).unwrap();
        }
    }
    
    // Only filter in non-dev mode
    if !dev_mode {
        kmer_map.retain(|_, &mut count| count >= min_coverage);
    }
    
    println!("Retained {} k-mers ({}%)", kmer_map.len(), 
             (kmer_map.len() as f32 / original_kmer_count as f32 * 100.0).round());
    
    // Assemble contigs using greedy algorithm
    println!("Assembling contigs...");
    let mut contigs = assemble_contigs(k, &kmer_map, min_length);
    
    // Filter short contigs (unless in dev mode)
    if !dev_mode {
        let pre_filter_count = contigs.len();
        contigs.retain(|c| c.sequence.len() >= min_length);
        println!("Filtered {} short contigs, retained {}", pre_filter_count - contigs.len(), contigs.len());
    } else if contigs.len() > 0 {
        // In dev mode, write all contigs regardless of length
        let debug_dir = Path::new(output_path).parent().unwrap_or(Path::new("."));
        let debug_file = debug_dir.join("debug_contigs.fasta");
        
        let mut debug_writer = FastaWriter::new(debug_file.to_str().unwrap());
        for (i, contig) in contigs.iter().enumerate() {
            if let Err(e) = debug_writer.write_contig(contig, i + 1) {
                eprintln!("Error writing debug contig {}: {}", i + 1, e);
            }
        }
    }
    
    println!("Assembled {} contigs", contigs.len());
    
    // Write output FASTA file
    println!("Writing output file {}...", output_path);
    let mut writer = FastaWriter::new(output_path);
    
    for (i, contig) in contigs.iter().enumerate() {
        if let Err(e) = writer.write_contig(contig, i + 1) {
            eprintln!("Error writing contig {}: {}", i + 1, e);
        }
    }
    
    println!("Done!");
}

// Improved version of greedy assembly with extension
#[allow(deprecated)]
fn assemble_contigs(k: usize, kmer_counts: &HashMap<String, u32>, min_len: usize) -> Vec<Contig> {
    let mut contigs = Vec::new();
    let mut used = HashSet::new();

    // Sort k-mers by abundance
    let mut sorted_kmers: Vec<(&String, &u32)> = kmer_counts.iter().collect();
    sorted_kmers.sort_by(|a, b| b.1.cmp(a.1));

    for (kmer, &_count) in sorted_kmers {
        if used.contains(kmer) {
            continue;
        }

        // Try to extend this k-mer in both directions
        let mut contig = kmer.clone();
        let mut path: Vec<u64> = vec![encode_kmer(kmer).unwrap_or(0)];
        used.insert(kmer.clone());

        // Forward extension
        let mut current = kmer.clone();
        loop {
            let suffix = &current[1..];
            let mut best_next = None;
            let mut best_count = 0;

            // Find the best extension
            for (ext_kmer, &ext_count) in kmer_counts.iter() {
                if used.contains(ext_kmer) {
                    continue;
                }

                if ext_kmer.starts_with(suffix) && ext_count > best_count {
                    best_next = Some(ext_kmer);
                    best_count = ext_count;
                }
            }

            if let Some(next) = best_next {
                used.insert(next.clone());
                path.push(encode_kmer(next).unwrap_or(0));
                let extension = &next[k-1..];
                contig.push_str(extension);
                current = next.clone();
            } else {
                break;
            }
        }

        if contig.len() >= min_len {
            contigs.push(Contig {
                id: contigs.len(),
                sequence: contig,
                kmer_path: path,
            });
        }
    }

    contigs
}
