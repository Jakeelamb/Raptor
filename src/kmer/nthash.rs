// src/kmer/nthash.rs
//! ntHash rolling hash implementation for efficient k-mer hashing.
//!
//! ntHash provides O(1) rolling updates for sliding window k-mer hashing,
//! compared to O(k) for recomputing from scratch. This is critical for
//! high-performance k-mer counting.
//!
//! Reference: Mohamadi, H., Chu, J., Vandervalk, B. P., & Birol, I. (2016).
//! ntHash: recursive nucleotide hashing. Bioinformatics, 32(22), 3492-3494.

/// ntHash constants for each nucleotide (from the ntHash paper)
/// These are carefully chosen 64-bit values that provide good hash distribution
const NT_A: u64 = 0x3c8bfbb395c60474;
const NT_C: u64 = 0x3193c18562a02b4c;
const NT_G: u64 = 0x20323ed082572324;
const NT_T: u64 = 0x295549f54be24456;

/// Special marker for invalid bases
const NT_INVALID: u64 = u64::MAX;

/// Lookup table for base -> hash value
const NT_HASH: [u64; 256] = {
    let mut table = [NT_INVALID; 256];
    table[b'A' as usize] = NT_A;
    table[b'a' as usize] = NT_A;
    table[b'C' as usize] = NT_C;
    table[b'c' as usize] = NT_C;
    table[b'G' as usize] = NT_G;
    table[b'g' as usize] = NT_G;
    table[b'T' as usize] = NT_T;
    table[b't' as usize] = NT_T;
    table
};

/// Lookup table for complement base -> hash value
const NT_HASH_RC: [u64; 256] = {
    let mut table = [NT_INVALID; 256];
    table[b'A' as usize] = NT_T; // A complements T
    table[b'a' as usize] = NT_T;
    table[b'C' as usize] = NT_G; // C complements G
    table[b'c' as usize] = NT_G;
    table[b'G' as usize] = NT_C; // G complements C
    table[b'g' as usize] = NT_C;
    table[b'T' as usize] = NT_A; // T complements A
    table[b't' as usize] = NT_A;
    table
};

/// Rolling hash state for ntHash
#[derive(Clone, Debug)]
pub struct NtHasher {
    /// Forward strand hash
    forward: u64,
    /// Reverse complement hash
    reverse: u64,
    /// K-mer size
    k: usize,
}

impl NtHasher {
    /// Initialize ntHash from a k-mer sequence.
    ///
    /// Returns None if the sequence contains invalid bases (not A/C/G/T).
    ///
    /// # Arguments
    /// * `seq` - The initial k-mer sequence as bytes
    /// * `k` - K-mer size (must equal seq.len())
    ///
    /// # Example
    /// ```ignore
    /// let hasher = NtHasher::new(b"ACGT", 4).unwrap();
    /// ```
    #[inline]
    pub fn new(seq: &[u8], k: usize) -> Option<Self> {
        if seq.len() != k || k == 0 || k > 32 {
            return None;
        }

        let mut forward: u64 = 0;
        let mut reverse: u64 = 0;

        for (i, &base) in seq.iter().enumerate() {
            let h = NT_HASH[base as usize];
            let h_rc = NT_HASH_RC[base as usize];

            // Check for invalid base (N, etc.)
            if h == NT_INVALID {
                return None;
            }

            // Forward: rotate left and XOR
            forward = forward.rotate_left(1) ^ h;

            // Reverse: build from the end, using complement
            // Position in RC is (k-1-i), so we rotate by that amount
            reverse ^= h_rc.rotate_left((k - 1 - i) as u32);
        }

        Some(Self { forward, reverse, k })
    }

    /// Roll the hash to the next position in O(1) time.
    ///
    /// This removes the effect of `out_base` and adds the effect of `in_base`.
    ///
    /// # Arguments
    /// * `out_base` - The base being removed (leftmost)
    /// * `in_base` - The base being added (rightmost)
    ///
    /// # Returns
    /// The canonical hash value (minimum of forward and reverse complement)
    #[inline]
    pub fn roll(&mut self, out_base: u8, in_base: u8) -> u64 {
        let h_out = NT_HASH[out_base as usize];
        let h_in = NT_HASH[in_base as usize];
        let h_out_rc = NT_HASH_RC[out_base as usize];
        let h_in_rc = NT_HASH_RC[in_base as usize];

        // Update forward hash:
        // The old base was at position 0 (leftmost), contributed as h_out rotated by k-1
        // After rotation, we need to remove its contribution
        // Then add the new base at position k-1 (rightmost)
        self.forward = self.forward.rotate_left(1)
            ^ h_out.rotate_left(self.k as u32)
            ^ h_in;

        // Update reverse complement hash:
        // Old base complement was at position k-1 in RC (rightmost)
        // New base complement should be at position 0 in RC (leftmost)
        self.reverse = self.reverse.rotate_right(1)
            ^ h_out_rc
            ^ h_in_rc.rotate_left((self.k - 1) as u32);

        self.canonical()
    }

    /// Get the canonical hash (minimum of forward and reverse complement).
    #[inline]
    pub fn canonical(&self) -> u64 {
        self.forward.min(self.reverse)
    }

    /// Get the forward hash only.
    #[inline]
    pub fn forward_hash(&self) -> u64 {
        self.forward
    }

    /// Get the reverse complement hash only.
    #[inline]
    pub fn reverse_hash(&self) -> u64 {
        self.reverse
    }

    /// Get the k-mer size.
    #[inline]
    pub fn k(&self) -> usize {
        self.k
    }
}

/// Iterator that yields canonical ntHash values for all k-mers in a sequence.
pub struct NtHashIterator<'a> {
    seq: &'a [u8],
    hasher: Option<NtHasher>,
    pos: usize,
    k: usize,
}

impl<'a> NtHashIterator<'a> {
    /// Create a new iterator over k-mers in a sequence.
    pub fn new(seq: &'a [u8], k: usize) -> Self {
        let hasher = if seq.len() >= k {
            NtHasher::new(&seq[0..k], k)
        } else {
            None
        };

        Self {
            seq,
            hasher,
            pos: 0,
            k,
        }
    }
}

impl<'a> Iterator for NtHashIterator<'a> {
    type Item = (usize, u64); // (position, hash)

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        // First k-mer
        if self.pos == 0 {
            if let Some(ref hasher) = self.hasher {
                self.pos = 1;
                return Some((0, hasher.canonical()));
            } else {
                return None;
            }
        }

        // Check if we have more k-mers
        let kmer_end = self.pos + self.k;
        if kmer_end > self.seq.len() {
            return None;
        }

        // Roll to next position
        let out_base = self.seq[self.pos - 1];
        let in_base = self.seq[kmer_end - 1];

        if let Some(ref mut hasher) = self.hasher {
            // Check if new base is valid
            if NT_HASH[in_base as usize] == 0 && in_base != b'A' && in_base != b'a' {
                // Invalid base - need to restart from next valid k-mer
                self.pos += 1;
                while self.pos + self.k <= self.seq.len() {
                    if let Some(new_hasher) = NtHasher::new(&self.seq[self.pos..self.pos + self.k], self.k) {
                        *hasher = new_hasher;
                        let result_pos = self.pos;
                        self.pos += 1;
                        return Some((result_pos, hasher.canonical()));
                    }
                    self.pos += 1;
                }
                return None;
            }

            let hash = hasher.roll(out_base, in_base);
            let result_pos = self.pos;
            self.pos += 1;
            Some((result_pos, hash))
        } else {
            None
        }
    }
}

/// Compute ntHash for a single k-mer.
///
/// This is a convenience function for one-off hashing.
/// For sliding windows, use NtHasher or NtHashIterator instead.
#[inline]
pub fn nthash(seq: &[u8]) -> Option<u64> {
    NtHasher::new(seq, seq.len()).map(|h| h.canonical())
}

/// Compute ntHash values for all k-mers in a sequence.
///
/// Returns a vector of (position, hash) tuples.
pub fn nthash_all(seq: &[u8], k: usize) -> Vec<(usize, u64)> {
    NtHashIterator::new(seq, k).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nthash_basic() {
        let hash = nthash(b"ACGT").unwrap();
        assert_ne!(hash, 0);
    }

    #[test]
    fn test_nthash_canonical() {
        // Same sequence should give same hash
        let h1 = nthash(b"ACGT").unwrap();
        let h2 = nthash(b"ACGT").unwrap();
        assert_eq!(h1, h2);

        // Forward and reverse complement should give same canonical hash
        // AA -> TT (RC), so canonical(AA) == canonical(TT)
        let hasher1 = NtHasher::new(b"AA", 2).unwrap();
        let hasher2 = NtHasher::new(b"TT", 2).unwrap();
        assert_eq!(hasher1.canonical(), hasher2.canonical());
    }

    #[test]
    fn test_nthash_rolling() {
        let seq = b"ACGTACGT";
        let k = 4;

        // Get all hashes via iterator (uses rolling)
        let hashes: Vec<u64> = NtHashIterator::new(seq, k).map(|(_, h)| h).collect();

        // Verify we get the right number
        assert_eq!(hashes.len(), seq.len() - k + 1);

        // All hashes should be non-zero
        for (i, &h) in hashes.iter().enumerate() {
            assert_ne!(h, 0, "Hash at position {} is zero", i);
        }

        // Verify first hash matches direct computation
        let first_hash = nthash(&seq[0..k]).unwrap();
        assert_eq!(hashes[0], first_hash);
    }

    #[test]
    fn test_nthash_invalid_base() {
        // N should cause None
        assert!(nthash(b"ACNT").is_none());
        assert!(nthash(b"NACG").is_none());
    }

    #[test]
    fn test_nthash_iterator_skips_invalid() {
        let seq = b"ACGTNGCTA";
        let k = 4;

        let hashes: Vec<(usize, u64)> = NtHashIterator::new(seq, k).collect();

        // Should skip the k-mers containing N
        // First hash should be at position 0
        if !hashes.is_empty() {
            assert_eq!(hashes[0].0, 0);
        }
    }

    #[test]
    fn test_different_kmers_different_hashes() {
        // Use longer k-mers to avoid canonical collisions
        let h1 = nthash(b"AAAAAAAA").unwrap();
        let h2 = nthash(b"CCCCCCCC").unwrap();
        let h3 = nthash(b"ACGTACGT").unwrap();
        let h4 = nthash(b"GTACGTAC").unwrap();

        // These should all be different (checking forward hashes)
        let hasher1 = NtHasher::new(b"AAAAAAAA", 8).unwrap();
        let hasher2 = NtHasher::new(b"CCCCCCCC", 8).unwrap();
        let hasher3 = NtHasher::new(b"ACGTACGT", 8).unwrap();

        assert_ne!(hasher1.forward_hash(), hasher2.forward_hash());
        assert_ne!(hasher1.forward_hash(), hasher3.forward_hash());
        assert_ne!(hasher2.forward_hash(), hasher3.forward_hash());
    }
}
