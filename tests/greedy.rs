#[allow(deprecated)]
use raptor::graph::assembler::greedy_assembly;
use raptor::kmer::kmer::{canonical_kmer, decode_kmer};
use std::collections::HashMap;

#[test]
#[allow(deprecated)]
fn test_greedy_extension_builds_correct_contig() {
    let mut kmers = HashMap::new();
    kmers.insert("ATGCG".to_string(), 10);
    kmers.insert("TGCGA".to_string(), 8);
    kmers.insert("GCGAT".to_string(), 6);
    kmers.insert("CGATT".to_string(), 4);

    // Debug canonical k-mers
    for kmer in kmers.keys() {
        println!("Kmer: {}, Canonical: {:?}", kmer, canonical_kmer(kmer));
    }

    let contigs = greedy_assembly(5, &kmers, 5);

    println!("Generated {} contigs:", contigs.len());
    for (i, contig) in contigs.iter().enumerate() {
        println!("Contig {}: {}", i, contig.sequence);
        let path_decoded: Vec<String> = contig.kmer_path.iter()
            .map(|&k| decode_kmer(k, 5))
            .collect();
        println!("  Path: {:?}", path_decoded);
    }

    assert_eq!(contigs.len(), 1);
    assert_eq!(contigs[0].sequence, "ATGCGATT");
}
