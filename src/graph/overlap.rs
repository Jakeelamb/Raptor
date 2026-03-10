use crate::accel::simd::match_kmers_with_overlap;
use crate::graph::assembler::Contig;
use rayon::prelude::*;

/// Find overlaps between contigs with a minimum length and maximum number of mismatches
pub fn find_overlaps(
    contigs: &[Contig],
    min_overlap: usize,
    max_mismatches: usize,
) -> Vec<(usize, usize, usize)> {
    // Process each source contig in parallel and avoid building an O(n²) pair list.
    let mut links: Vec<(usize, usize, usize)> = (0..contigs.len())
        .into_par_iter()
        .map(|i| {
            let mut local = Vec::new();
            let from = &contigs[i].sequence;
            for j in 0..contigs.len() {
                if i == j {
                    continue;
                }
                let to = &contigs[j].sequence;
                if let Some((shift, _)) =
                    match_kmers_with_overlap(from, to, min_overlap, max_mismatches)
                {
                    let overlap_len = from.len() - shift;
                    local.push((i, j, overlap_len));
                }
            }
            local
        })
        .flatten()
        .collect();

    links.sort_unstable_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
    links
}

#[cfg(test)]
mod tests {
    use super::find_overlaps;
    use crate::graph::assembler::Contig;

    #[test]
    fn find_overlaps_has_stable_pair_order() {
        let contigs = vec![
            Contig {
                id: 0,
                sequence: "AAACCC".to_string(),
                kmer_path: Vec::new(),
            },
            Contig {
                id: 1,
                sequence: "CCCGGG".to_string(),
                kmer_path: Vec::new(),
            },
            Contig {
                id: 2,
                sequence: "GGGTTT".to_string(),
                kmer_path: Vec::new(),
            },
        ];

        let overlaps = find_overlaps(&contigs, 3, 0);
        assert!(!overlaps.is_empty());
        assert!(overlaps
            .windows(2)
            .all(|w| (w[0].0, w[0].1) <= (w[1].0, w[1].1)));
    }
}
