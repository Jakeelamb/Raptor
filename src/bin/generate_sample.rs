use rand::{Rng, thread_rng};
use rand::distributions::Uniform;
use std::fs::File;
use std::io::{BufWriter, Write};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let output_path = if args.len() > 1 {
        &args[1]
    } else {
        "sample_large.fastq"
    };
    
    let num_reads = 300; // More reads for better statistics
    let read_length = 50;
    
    println!("Generating {} reads of {}bp each to {}", num_reads, read_length, output_path);
    
    let file = File::create(output_path).expect("Failed to create output file");
    let mut writer = BufWriter::new(file);
    let mut rng = thread_rng();
    
    // Distribution for random base selection: 0=A, 1=C, 2=G, 3=T
    let base_dist = Uniform::from(0..4);
    
    // Create some motifs that will repeat across reads with different frequencies
    let common_motifs = vec![
        "GATTACA",       // Will appear in ~40% of reads
        "ACGTACGT",      // Will appear in ~30% of reads
    ];
    
    let rare_motifs = vec![
        "TTAGGGTTAGGG",  // Will appear in ~5% of reads
        "AAAAAAAAA",     // Will appear in ~5% of reads
        "CTGACTGA",      // Will appear in ~5% of reads
    ];
    
    // Also create some "high coverage regions" - identical reads that appear multiple times
    let num_high_coverage_reads = 50;
    let coverage_factor = 4; // Each high-coverage read will be repeated this many times
    
    // First, generate the high-coverage reads
    let mut high_coverage_reads = Vec::new();
    for _ in 0..num_high_coverage_reads {
        let mut sequence = String::with_capacity(read_length);
        // Guaranteed to include one of the common motifs to make them "special"
        let motif_idx = rng.gen_range(0..common_motifs.len());
        let motif = &common_motifs[motif_idx];
        let pos = rng.gen_range(0..(read_length - motif.len()));
        
        // Generate random sequence before motif
        for _ in 0..pos {
            sequence.push(base_nucleotide(&mut rng, base_dist));
        }
        
        // Add the motif
        sequence.push_str(motif);
        
        // Generate random sequence after motif
        for _ in (pos + motif.len())..read_length {
            sequence.push(base_nucleotide(&mut rng, base_dist));
        }
        
        high_coverage_reads.push(sequence);
    }
    
    let mut read_id = 1;
    
    // Then output each high-coverage read multiple times
    for sequence in &high_coverage_reads {
        for _ in 0..coverage_factor {
            // Generate quality scores
            let quality = "I".repeat(read_length);
            
            // Write FASTQ record
            writeln!(writer, "@SEQ_{}", read_id).unwrap();
            writeln!(writer, "{}", sequence).unwrap();
            writeln!(writer, "+").unwrap();
            writeln!(writer, "{}", quality).unwrap();
            
            read_id += 1;
        }
    }
    
    // Generate the remaining normal reads
    let remaining_reads = num_reads - (num_high_coverage_reads * coverage_factor);
    for _ in 0..remaining_reads {
        let sequence = if rng.gen_bool(0.4) {
            // 40% chance to include a common motif
            let motif_idx = rng.gen_range(0..common_motifs.len());
            generate_sequence_with_motif(&mut rng, read_length, &common_motifs[motif_idx], base_dist)
        } else if rng.gen_bool(0.15) {
            // ~15% chance (of the remaining 60%) to include a rare motif
            let motif_idx = rng.gen_range(0..rare_motifs.len());
            generate_sequence_with_motif(&mut rng, read_length, &rare_motifs[motif_idx], base_dist)
        } else {
            // Remaining ~45% are completely random
            (0..read_length)
                .map(|_| base_nucleotide(&mut rng, base_dist))
                .collect()
        };
        
        // Generate quality scores
        let quality = "I".repeat(read_length);
        
        // Write FASTQ record
        writeln!(writer, "@SEQ_{}", read_id).unwrap();
        writeln!(writer, "{}", sequence).unwrap();
        writeln!(writer, "+").unwrap();
        writeln!(writer, "{}", quality).unwrap();
        
        read_id += 1;
    }
    
    println!("Sample file generated successfully!");
}

// Helper function to generate a random base
fn base_nucleotide(rng: &mut impl Rng, dist: Uniform<u8>) -> char {
    match rng.sample(dist) {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => unreachable!(),
    }
}

// Helper function to generate a sequence with a motif
fn generate_sequence_with_motif(rng: &mut impl Rng, length: usize, motif: &str, dist: Uniform<u8>) -> String {
    let mut sequence = String::with_capacity(length);
    
    // If motif is longer than the sequence, truncate it
    if motif.len() >= length {
        return motif[0..length].to_string();
    }
    
    // Choose position to insert motif
    let pos = rng.gen_range(0..(length - motif.len()));
    
    // Generate random sequence before motif
    for _ in 0..pos {
        sequence.push(base_nucleotide(rng, dist));
    }
    
    // Add the motif
    sequence.push_str(motif);
    
    // Generate random sequence after motif
    for _ in (pos + motif.len())..length {
        sequence.push(base_nucleotide(rng, dist));
    }
    
    sequence
} 