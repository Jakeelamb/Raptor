use crate::accel::simd::match_kmers_with_overlap;
use crate::graph::assembler::Contig;
use rayon::prelude::*;

/// Find overlaps between contigs with a minimum length and maximum number of mismatches
pub fn find_overlaps(contigs: &[Contig], min_overlap: usize, max_mismatches: usize) -> Vec<(usize, usize, usize)> {
    let mut links = Vec::new();
    
    // Parallel iteration over all pairs of contigs
    let contig_pairs: Vec<_> = (0..contigs.len())
        .flat_map(|i| (0..contigs.len()).map(move |j| (i, j)))
        .filter(|&(i, j)| i != j) // Don't compare a contig with itself
        .collect();
    
    let overlaps: Vec<_> = contig_pairs.par_iter()
        .filter_map(|&(i, j)| {
            let from = &contigs[i].sequence;
            let to = &contigs[j].sequence;
            
            // Find suffix-prefix overlap
            match_kmers_with_overlap(from, to, min_overlap, max_mismatches)
                .map(|(shift, _)| {
                    let overlap_len = from.len() - shift;
                    (i, j, overlap_len)
                })
        })
        .collect();
    
    links.extend(overlaps);
    links
} 