pub type Kmer = String;

/// Compact k-mer representation using 2 bits per nucleotide.
/// Supports k-mers up to 32bp (64 bits / 2 bits per base).
/// Memory: 9 bytes vs ~55 bytes for String-based k-mers.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct KmerU64 {
    /// 2-bit encoded nucleotides: A=00, C=01, G=10, T=11
    pub encoded: u64,
    /// Length of the k-mer (1-32)
    pub len: u8,
}

impl KmerU64 {
    /// Create a new KmerU64 from a byte slice. Zero-copy encoding.
    /// Returns None if sequence contains invalid nucleotides or is > 32bp.
    #[inline]
    pub fn from_slice(seq: &[u8]) -> Option<Self> {
        if seq.len() > 32 || seq.is_empty() {
            return None;
        }
        let mut encoded: u64 = 0;
        for &b in seq {
            encoded <<= 2;
            encoded |= match b {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => return None,
            };
        }
        Some(Self {
            encoded,
            len: seq.len() as u8,
        })
    }

    /// Create from a string slice (convenience wrapper).
    #[inline]
    pub fn from_str(seq: &str) -> Option<Self> {
        Self::from_slice(seq.as_bytes())
    }

    /// Compute reverse complement using bit manipulation only.
    /// Complement: A<->T (00<->11), C<->G (01<->10) = XOR with 11
    /// Then reverse the 2-bit pairs.
    #[inline]
    pub fn reverse_complement(&self) -> Self {
        let k = self.len as u32;
        // Complement all bases (XOR with 0b11 for each 2-bit pair)
        let mask = (1u64 << (k * 2)) - 1;
        let complemented = self.encoded ^ mask;

        // Reverse the 2-bit pairs
        let mut reversed: u64 = 0;
        let mut temp = complemented;
        for _ in 0..k {
            reversed = (reversed << 2) | (temp & 0b11);
            temp >>= 2;
        }

        Self {
            encoded: reversed,
            len: self.len,
        }
    }

    /// Return the canonical form (lexicographically smaller of self and reverse complement).
    /// Pure bit manipulation, no allocations.
    #[inline]
    pub fn canonical(&self) -> Self {
        let rc = self.reverse_complement();
        if self.encoded <= rc.encoded {
            *self
        } else {
            rc
        }
    }

    /// Decode back to a String (for output/debugging).
    pub fn decode(&self) -> String {
        let mut result = String::with_capacity(self.len as usize);
        let temp = self.encoded;
        let bases = [b'A', b'C', b'G', b'T'];

        // Extract bases from most significant to least significant
        for i in (0..self.len).rev() {
            let shift = i as u32 * 2;
            let base_idx = ((temp >> shift) & 0b11) as usize;
            result.push(bases[base_idx] as char);
        }
        result
    }

    /// Get the suffix of length k-1 (drop first base).
    #[inline]
    pub fn suffix(&self) -> u64 {
        let mask = (1u64 << ((self.len as u32 - 1) * 2)) - 1;
        self.encoded & mask
    }

    /// Get the prefix of length k-1 (drop last base).
    #[inline]
    pub fn prefix(&self) -> u64 {
        self.encoded >> 2
    }

    /// Extend the k-mer by appending a base (for sliding window).
    /// Returns the new encoded value with the same length (drops first base).
    #[inline]
    pub fn extend(&self, base: u8) -> Option<Self> {
        let base_bits = match base {
            b'A' | b'a' => 0u64,
            b'C' | b'c' => 1u64,
            b'G' | b'g' => 2u64,
            b'T' | b't' => 3u64,
            _ => return None,
        };
        let mask = (1u64 << (self.len as u32 * 2)) - 1;
        let new_encoded = ((self.encoded << 2) | base_bits) & mask;
        Some(Self {
            encoded: new_encoded,
            len: self.len,
        })
    }
}

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

/// Decode a u64-encoded k-mer back to a String.
pub fn decode_kmer(encoded: u64, k: usize) -> String {
    let mut result = String::with_capacity(k);
    let bases = ['A', 'C', 'G', 'T'];
    for i in (0..k).rev() {
        let shift = i * 2;
        let base_idx = ((encoded >> shift) & 0b11) as usize;
        result.push(bases[base_idx]);
    }
    result
}

/// Compute canonical form of u64-encoded k-mer (returns smaller of fwd/rc).
#[inline]
pub fn canonical_kmer_u64(encoded: u64, k: usize) -> u64 {
    let rc = reverse_complement_u64(encoded, k);
    if encoded <= rc { encoded } else { rc }
}

/// Compute reverse complement of u64-encoded k-mer using bit manipulation.
#[inline]
pub fn reverse_complement_u64(encoded: u64, k: usize) -> u64 {
    let mask = (1u64 << (k * 2)) - 1;
    let complemented = encoded ^ mask;

    let mut reversed: u64 = 0;
    let mut temp = complemented;
    for _ in 0..k {
        reversed = (reversed << 2) | (temp & 0b11);
        temp >>= 2;
    }
    reversed
}

/// Returns the canonical form of a k-mer (lexicographically smaller of forward and reverse complement)
///
/// DEPRECATED: Use KmerU64::canonical() for better performance (9 bytes vs ~55 bytes per k-mer).
#[deprecated(note = "Use KmerU64::canonical() for 6x memory reduction")]
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
