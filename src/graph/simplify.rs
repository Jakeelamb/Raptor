use crate::graph::assembler::Contig;
use crate::kmer::rle::rle_encode;

/// Collapse contigs with identical or highly similar RLE-encoded sequences
pub fn collapse_repeats(contigs: Vec<Contig>, _min_repeat_len: usize) -> Vec<Contig> {
    let mut collapsed: Vec<Contig> = Vec::new();
    let mut processed: Vec<bool> = vec![false; contigs.len()];
    
    // Process each contig
    for i in 0..contigs.len() {
        if processed[i] {
            continue;
        }
        
        // Mark current contig as processed
        processed[i] = true;
        collapsed.push(contigs[i].clone());
        
        // Find and mark all similar contigs
        let seq1 = &contigs[i].sequence;
        for j in (i+1)..contigs.len() {
            if processed[j] {
                continue;
            }
            
            let seq2 = &contigs[j].sequence;
            
            // Check if sequences are similar
            // For test consistency, consider identical sequences as similar
            if seq1 == seq2 || is_rle_similar(seq1, seq2, 0.9) {
                processed[j] = true; // Mark as processed, don't add to collapsed
            }
        }
    }
    
    collapsed
}

/// Determine if two sequences are similar based on their run-length encoding
/// Compares the base types and their frequencies
fn is_rle_similar(seq1: &str, seq2: &str, threshold: f64) -> bool {
    let rle1 = rle_encode(seq1);
    let rle2 = rle_encode(seq2);
    
    // If lengths are too different, consider them dissimilar
    if rle1.is_empty() || rle2.is_empty() || 
       (rle1.len() as f64 - rle2.len() as f64).abs() / rle1.len() as f64 > 0.3 {
        return false;
    }
    
    // Compare the base types and their frequencies
    let mut matches = 0;
    let total = rle1.len().max(rle2.len());
    
    for i in 0..rle1.len().min(rle2.len()) {
        if rle1[i].0 == rle2[i].0 && 
           (rle1[i].1 as f64 - rle2[i].1 as f64).abs() / rle1[i].1.max(1) as f64 <= 0.2 {
            matches += 1;
        }
    }
    
    (matches as f64 / total as f64) >= threshold
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_is_rle_similar() {
        // The RLE of "ATATATATATAT" might be [(A,1), (T,1)] * 6
        // The RLE of "ATATATATATAT" should be similar to itself
        assert!(is_rle_similar("ATATATATATAT", "ATATATATATAT", 0.9));
        
        // Should be similar with small changes
        assert!(is_rle_similar("ATATATATATAT", "ATATATATCTAT", 0.8));
        
        // Should not be similar with large changes
        assert!(!is_rle_similar("ATATATATATAT", "GCGCGCGCGCGC", 0.5));
    }
    
    #[test]
    fn test_collapse_repeats() {
        let contig1 = Contig {
            id: 1,
            sequence: "ATATATATATAT".to_string(),
            kmer_path: vec!["ATG".to_string(), "TGC".to_string()]
        };
        
        let contig2 = Contig {
            id: 2,
            sequence: "ATATATATATAT".to_string(),
            kmer_path: vec!["GCA".to_string(), "CAT".to_string()]
        };
        
        let contig3 = Contig {
            id: 3,
            sequence: "GCGCGCGCGCGC".to_string(),
            kmer_path: vec!["GCG".to_string(), "CGC".to_string()]
        };
        
        let contigs = vec![contig1, contig2, contig3];
        
        // Collapse with high similarity threshold (0.9)
        let collapsed = collapse_repeats(contigs.clone(), 0);
        
        // Should merge the first two contigs
        assert_eq!(collapsed.len(), 2);
        
        // Verify that one of the collapsed contigs is the similarity one
        let contains_similar = collapsed.iter().any(|c| c.sequence == "ATATATATATAT");
        assert!(contains_similar);
        
        // Verify that the non-similar contig is still there
        let contains_different = collapsed.iter().any(|c| c.sequence == "GCGCGCGCGCGC");
        assert!(contains_different);
    }
} 