use raptor::io::fastq::{open_fastq, stream_fastq_records, FastqRecord};
use raptor::io::fasta::FastaWriter;
use raptor::kmer::kmer::{KmerU64, decode_kmer};
use raptor::graph::assembler::Contig;
use rayon::ThreadPoolBuilder;
use ahash::{AHashMap, AHashSet};
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

    // Read input FASTQ file using streaming for memory efficiency
    println!("Reading input file {}...", input_path);
    let reader = open_fastq(input_path);
    let records: Vec<FastqRecord> = stream_fastq_records(reader).collect();
    println!("Read {} records", records.len());
    
    // Build k-mer map from reads using u64 encoding for memory efficiency
    println!("Building k-mer map with ntHash...");
    let mut kmer_map: AHashMap<u64, u32> = AHashMap::new();
    for record in &records {
        let bytes = record.sequence.as_bytes();
        // Use KmerU64 for canonical encoding
        if bytes.len() >= k {
            if let Some(mut kmer) = KmerU64::from_slice(&bytes[0..k]) {
                let canonical = kmer.canonical();
                *kmer_map.entry(canonical.encoded).or_insert(0) += 1;

                // Sliding window
                for i in k..bytes.len() {
                    if let Some(next) = kmer.extend(bytes[i]) {
                        kmer = next;
                        let canonical = kmer.canonical();
                        *kmer_map.entry(canonical.encoded).or_insert(0) += 1;
                    } else if let Some(fresh) = KmerU64::from_slice(&bytes[i + 1 - k..i + 1]) {
                        kmer = fresh;
                        let canonical = kmer.canonical();
                        *kmer_map.entry(canonical.encoded).or_insert(0) += 1;
                    }
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
        for (&kmer_encoded, &count) in &kmer_map {
            let kmer_str = decode_kmer(kmer_encoded, k);
            writeln!(&mut file, "{}\t{}", kmer_str, count).unwrap();
        }
    }
    
    // Only filter in non-dev mode
    if !dev_mode {
        kmer_map.retain(|_, &mut count| count >= min_coverage);
    }
    
    println!("Retained {} k-mers ({}%)", kmer_map.len(), 
             (kmer_map.len() as f32 / original_kmer_count as f32 * 100.0).round());
    
    // Assemble contigs using greedy algorithm with u64 encoding
    println!("Assembling contigs...");
    let mut contigs = assemble_contigs_u64(k, &kmer_map, min_length);
    
    // Filter short contigs (unless in dev mode)
    if !dev_mode {
        let pre_filter_count = contigs.len();
        contigs.retain(|c| c.sequence.len() >= min_length);
        println!("Filtered {} short contigs, retained {}", pre_filter_count - contigs.len(), contigs.len());
    } else if !contigs.is_empty() {
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

/// Greedy assembly using u64-encoded k-mers for memory efficiency.
/// Uses suffix matching via bit operations.
fn assemble_contigs_u64(k: usize, kmer_counts: &AHashMap<u64, u32>, min_len: usize) -> Vec<Contig> {
    const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
    let mut contigs = Vec::new();
    let mut used: AHashSet<u64> = AHashSet::new();

    // Sort k-mers by abundance (descending)
    let mut sorted_kmers: Vec<(&u64, &u32)> = kmer_counts.iter().collect();
    sorted_kmers.sort_by(|a, b| b.1.cmp(a.1));

    // Mask for suffix (k-1 bases)
    let suffix_mask: u64 = (1u64 << ((k - 1) * 2)) - 1;

    for (&kmer, &_count) in sorted_kmers {
        if used.contains(&kmer) {
            continue;
        }

        // Start with this k-mer
        let mut contig = decode_kmer(kmer, k);
        let mut path: Vec<u64> = vec![kmer];
        used.insert(kmer);

        // Forward extension
        let mut current = kmer;
        loop {
            // Get suffix (last k-1 bases)
            let suffix = current & suffix_mask;

            let mut best_next: Option<u64> = None;
            let mut best_count = 0u32;

            // Try all 4 possible extensions
            for base in 0u64..4 {
                let next = (suffix << 2) | base;
                // Get canonical form
                let next_canonical = raptor::kmer::kmer::canonical_kmer_u64(next, k);

                if let Some(&ext_count) = kmer_counts.get(&next_canonical) {
                    if !used.contains(&next_canonical) && ext_count > best_count {
                        best_next = Some(next_canonical);
                        best_count = ext_count;
                    }
                }
            }

            if let Some(next) = best_next {
                used.insert(next);
                path.push(next);
                // Extract last base and append
                let last_base = BASES[(next & 0b11) as usize];
                contig.push(last_base);
                current = next;
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
