pub type Kmer = String;

/// Encodes a DNA k-mer to a 64-bit integer (2 bits per nucleotide, max 32-mer)
pub fn encode_kmer(seq: &str) -> Option<u64> {
    let mut val: u64 = 0;
    for &b in seq.as_bytes() {
        val <<= 2;
        val |= match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => return None,
        };
    }
    Some(val)
}

/// Returns the canonical form of a k-mer (lexicographically smaller of forward and reverse complement)
pub fn canonical_kmer(seq: &str) -> Option<String> {
    // Bail early for empty sequences
    if seq.is_empty() {
        return Some(String::new());
    }

    // Check all characters are valid nucleotides
    if !seq.chars().all(|c| matches!(c, 'A' | 'C' | 'G' | 'T')) {
        return None;
    }

    // Get the reverse complement
    let rc = reverse_complement(seq);
    
    // Return the lexicographically smaller one
    // Important: We preserve the original string format without converting to uppercase
    if seq.to_string() <= rc {
        Some(seq.to_string())
    } else {
        Some(rc)
    }
}

/// Returns the reverse complement of a DNA sequence
pub fn reverse_complement(seq: &str) -> String {
    // For assembler purposes, we can directly use the SIMD version as it's now correctly implemented
    crate::accel::simd::reverse_complement_simd(seq)
}
