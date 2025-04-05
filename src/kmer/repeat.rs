use crate::kmer::rle::rle_encode;

/// Detect repeats in a DNA sequence using run-length encoding
pub fn detect_rle_repeats(seq: &str, min_len: usize, min_unit: usize) -> Vec<String> {
    let rle = rle_encode(seq);
    
    let mut repeats = Vec::new();
    
    for i in 0..rle.len() {
        let (base, count) = rle[i];
        if count as usize >= min_len {
            repeats.push(format!("{}{}", base as char, count));
        }
    }
    
    // Filter by minimum number of units
    if min_unit > 1 {
        repeats.retain(|r| {
            let count = r.chars().skip(1).collect::<String>().parse::<usize>().unwrap_or(0);
            count >= min_unit
        });
    }
    
    repeats
}

/// Find positions of a repeat in a DNA sequence
pub fn find_repeat_positions(seq: &str, pattern: &str, min_len: usize) -> Vec<usize> {
    let mut positions = Vec::new();
    let mut i = 0;
    
    while i + pattern.len() <= seq.len() {
        if &seq[i..i + pattern.len()] == pattern {
            let mut j = i + pattern.len();
            while j + pattern.len() <= seq.len() && &seq[j..j + pattern.len()] == pattern {
                j += pattern.len();
            }
            
            if j - i >= min_len {
                positions.push(i);
            }
            
            i = j;
        } else {
            i += 1;
        }
    }
    
    positions
}

/// Collapse repeats in a sequence using RLE for display
pub fn collapse_repeats(seq: &str, min_length: usize) -> String {
    let rle = rle_encode(seq);
    let mut output = String::new();
    
    for (base, count) in rle {
        if count as usize >= min_length {
            output.push_str(&format!("[{}]{}", base as char, count));
        } else {
            for _ in 0..count {
                output.push(base as char);
            }
        }
    }
    
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_detect_repeats() {
        // Test sequence with simple repeats
        let seq = "AAAACCCCCGGGGGTTTTT";
        let repeats = detect_rle_repeats(seq, 4, 1);
        
        // Should find A4, C5, G5, T5 - all 4 homopolymers
        assert_eq!(repeats.len(), 4);
        assert!(repeats.contains(&"A4".to_string()));
        assert!(repeats.contains(&"C5".to_string()));
        assert!(repeats.contains(&"G5".to_string()));
        assert!(repeats.contains(&"T5".to_string()));
        
        // Test with higher min_unit (should find fewer repeats)
        let repeats = detect_rle_repeats(seq, 4, 5);
        // With min_unit=5, we only expect C5, G5, T5 (not A4)
        assert_eq!(repeats.len(), 3);
    }
    
    #[test]
    fn test_collapse_repeats() {
        // Test homopolymer compression
        let seq = "AAAACGGGGGT";
        let collapsed = collapse_repeats(seq, 4);
        assert_eq!(collapsed, "[A]4C[G]5T");
        
        // Test with mixed repeats and non-repeats
        let seq = "AAACCCGGTT";
        let collapsed = collapse_repeats(seq, 3);
        // Expect: AAA (3 As) -> [A]3, CCC (3 Cs) -> [C]3, GG (2 Gs) -> GG, TT (2 Ts) -> TT
        let expected = "[A]3[C]3GGTT";
        println!("Expected: {}, Actual: {}", expected, collapsed);
        assert_eq!(collapsed, expected);
    }
} 