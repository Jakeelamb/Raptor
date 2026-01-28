/// 2-bit k-mer encoding utilities for GPU processing
///
/// DNA bases are encoded as:
/// - A = 00 (0)
/// - C = 01 (1)
/// - G = 10 (2)
/// - T = 11 (3)
///
/// This allows k-mers up to k=32 to fit in a single u64

/// Encode a single nucleotide to 2 bits
#[inline]
pub fn encode_base(base: u8) -> Option<u64> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None, // N or other invalid bases
    }
}

/// Decode 2 bits back to a nucleotide
#[inline]
pub fn decode_base(bits: u64) -> char {
    match bits & 0b11 {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => unreachable!(),
    }
}

/// Get the complement of an encoded base
#[inline]
pub fn complement_base(bits: u64) -> u64 {
    // A(00) <-> T(11), C(01) <-> G(10)
    bits ^ 0b11
}

/// Encode a k-mer string to a u64
///
/// Returns None if the k-mer contains invalid bases (N) or is too long (k > 32)
pub fn encode_kmer(kmer: &str) -> Option<u64> {
    if kmer.len() > 32 {
        return None;
    }

    let mut encoded = 0u64;
    for base in kmer.bytes() {
        encoded = (encoded << 2) | encode_base(base)?;
    }
    Some(encoded)
}

/// Encode a k-mer from bytes to a u64
pub fn encode_kmer_bytes(kmer: &[u8]) -> Option<u64> {
    if kmer.len() > 32 {
        return None;
    }

    let mut encoded = 0u64;
    for &base in kmer {
        encoded = (encoded << 2) | encode_base(base)?;
    }
    Some(encoded)
}

/// Decode a u64 back to a k-mer string
pub fn decode_kmer(encoded: u64, k: usize) -> String {
    let mut result = String::with_capacity(k);
    for i in (0..k).rev() {
        let bits = (encoded >> (i * 2)) & 0b11;
        result.push(decode_base(bits));
    }
    result
}

/// Compute the reverse complement of an encoded k-mer
///
/// For a k-mer of length k encoded in 2*k bits:
/// 1. Complement each base (XOR with 11)
/// 2. Reverse the order of bases
pub fn reverse_complement_encoded(encoded: u64, k: usize) -> u64 {
    let mut result = 0u64;
    let mut temp = encoded;

    for _ in 0..k {
        // Get the last base, complement it, and add to result
        let base = temp & 0b11;
        result = (result << 2) | complement_base(base);
        temp >>= 2;
    }

    result
}

/// Get the canonical form of an encoded k-mer (lexicographically smaller of kmer and its reverse complement)
pub fn canonical_encoded(encoded: u64, k: usize) -> u64 {
    let rc = reverse_complement_encoded(encoded, k);
    std::cmp::min(encoded, rc)
}

/// Encode a k-mer and return its canonical form
pub fn canonical_encode(kmer: &str) -> Option<u64> {
    let encoded = encode_kmer(kmer)?;
    Some(canonical_encoded(encoded, kmer.len()))
}

/// MurmurHash3 finalizer for better hash distribution
#[inline]
pub fn murmur_hash3_finalize(mut h: u64) -> u64 {
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51afd7ed558ccd);
    h ^= h >> 33;
    h = h.wrapping_mul(0xc4ceb9fe1a85ec53);
    h ^= h >> 33;
    h
}

/// Hash an encoded k-mer using MurmurHash3 finalizer
#[inline]
pub fn hash_kmer(encoded: u64) -> u64 {
    murmur_hash3_finalize(encoded)
}

/// GPU-friendly k-mer table structure
#[derive(Debug, Clone)]
pub struct GpuKmerTable {
    /// 2-bit encoded k-mers
    pub kmers_encoded: Vec<u64>,
    /// Counts for each k-mer
    pub counts: Vec<u32>,
    /// K-mer size
    pub k: usize,
}

impl GpuKmerTable {
    pub fn new(k: usize) -> Self {
        Self {
            kmers_encoded: Vec::new(),
            counts: Vec::new(),
            k,
        }
    }

    pub fn with_capacity(k: usize, capacity: usize) -> Self {
        Self {
            kmers_encoded: Vec::with_capacity(capacity),
            counts: Vec::with_capacity(capacity),
            k,
        }
    }

    pub fn insert(&mut self, encoded: u64, count: u32) {
        self.kmers_encoded.push(encoded);
        self.counts.push(count);
    }

    pub fn len(&self) -> usize {
        self.kmers_encoded.len()
    }

    pub fn is_empty(&self) -> bool {
        self.kmers_encoded.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_base() {
        assert_eq!(encode_base(b'A'), Some(0));
        assert_eq!(encode_base(b'C'), Some(1));
        assert_eq!(encode_base(b'G'), Some(2));
        assert_eq!(encode_base(b'T'), Some(3));
        assert_eq!(encode_base(b'N'), None);

        assert_eq!(decode_base(0), 'A');
        assert_eq!(decode_base(1), 'C');
        assert_eq!(decode_base(2), 'G');
        assert_eq!(decode_base(3), 'T');
    }

    #[test]
    fn test_encode_decode_kmer() {
        let kmer = "ACGT";
        let encoded = encode_kmer(kmer).unwrap();
        let decoded = decode_kmer(encoded, 4);
        assert_eq!(decoded, kmer);

        // Test longer k-mer
        let kmer = "ACGTACGTACGTACGT";
        let encoded = encode_kmer(kmer).unwrap();
        let decoded = decode_kmer(encoded, 16);
        assert_eq!(decoded, kmer);
    }

    #[test]
    fn test_reverse_complement() {
        // ACGT -> ACGT (palindrome)
        let kmer = "ACGT";
        let encoded = encode_kmer(kmer).unwrap();
        let rc = reverse_complement_encoded(encoded, 4);
        let rc_decoded = decode_kmer(rc, 4);
        assert_eq!(rc_decoded, "ACGT");

        // AAAA -> TTTT
        let kmer = "AAAA";
        let encoded = encode_kmer(kmer).unwrap();
        let rc = reverse_complement_encoded(encoded, 4);
        let rc_decoded = decode_kmer(rc, 4);
        assert_eq!(rc_decoded, "TTTT");

        // AACG -> CGTT
        let kmer = "AACG";
        let encoded = encode_kmer(kmer).unwrap();
        let rc = reverse_complement_encoded(encoded, 4);
        let rc_decoded = decode_kmer(rc, 4);
        assert_eq!(rc_decoded, "CGTT");
    }

    #[test]
    fn test_canonical() {
        // AACG vs CGTT -> AACG is smaller
        let kmer = "AACG";
        let canonical = canonical_encode(kmer).unwrap();
        let canonical_str = decode_kmer(canonical, 4);
        assert_eq!(canonical_str, "AACG");

        // TTAA vs TTAA (reverse complement) -> TTAA
        let kmer = "TTAA";
        let canonical = canonical_encode(kmer).unwrap();
        let canonical_str = decode_kmer(canonical, 4);
        // TTAA reverse complement is TTAA, so it should be TTAA
        assert_eq!(canonical_str, "TTAA");
    }

    #[test]
    fn test_max_kmer_length() {
        // 32-mers should work
        let kmer = "ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32 bases
        assert!(encode_kmer(kmer).is_some());

        // 33-mers should fail
        let kmer = "ACGTACGTACGTACGTACGTACGTACGTACGTA"; // 33 bases
        assert!(encode_kmer(kmer).is_none());
    }

    #[test]
    fn test_invalid_bases() {
        // K-mers with N should fail
        assert!(encode_kmer("ACNT").is_none());
        assert!(encode_kmer("NACGT").is_none());
    }
}
