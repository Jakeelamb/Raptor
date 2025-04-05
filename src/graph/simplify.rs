use crate::graph::assembler::Contig;
use crate::kmer::rle::rle_encode;

/// Collapse contigs with identical or highly similar RLE-encoded sequences
pub fn collapse_repeats(contigs: Vec<Contig>, min_repeat_len: usize) -> Vec<Contig> {
    let mut collapsed: Vec<Contig> = Vec::new();
    let mut seen: Vec<Vec<(u8, u8)>> = Vec::new();

    'outer: for contig in contigs {
        let rle = rle_encode(&contig.sequence);

        for existing in &seen {
            let match_len = rle.iter()
                .zip(existing.iter())
                .take_while(|&(&(a, ac), &(b, bc))| a == b && ac == bc)
                .count();

            if match_len >= min_repeat_len {
                // Collapse repeat — keep one representative
                continue 'outer;
            }
        }

        seen.push(rle);
        collapsed.push(contig);
    }

    collapsed
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_collapse_repeats() {
        // Create sample contigs with similar RLE patterns
        let contig1 = Contig {
            sequence: "AAAACCCCCGGGGTTTT".to_string(),
            kmer_path: vec!["AAAAC".to_string()],
        };
        
        let contig2 = Contig {
            sequence: "AAAACCCCGGGGTTTT".to_string(), // Similar but one less C
            kmer_path: vec!["AAAAC".to_string()],
        };
        
        let contig3 = Contig {
            sequence: "TTTTGGGGCCCCAAAA".to_string(), // Different pattern
            kmer_path: vec!["TTTTG".to_string()],
        };
        
        // Test with min_repeat_len = 2 (should collapse contigs with at least 2 matching RLE tuples)
        let result = collapse_repeats(vec![contig1.clone(), contig2.clone(), contig3.clone()], 2);
        // All 3 are kept because the RLE patterns are different enough with default check
        assert_eq!(result.len(), 3);
        
        // Test with min_repeat_len = 4 (this test uses different contigs with more distinct patterns)
        let contigs = vec![
            Contig {
                sequence: "AAAAAAGGGGGG".to_string(),
                kmer_path: vec![],
            },
            Contig {
                sequence: "TTTTTTCCCCCC".to_string(),
                kmer_path: vec![],
            },
            Contig {
                sequence: "GGGGGGAAAAAA".to_string(),
                kmer_path: vec![],
            }
        ];
        let result = collapse_repeats(contigs, 4);
        // All 3 should be preserved as they have distinct RLE patterns
        assert_eq!(result.len(), 3);
    }
} 