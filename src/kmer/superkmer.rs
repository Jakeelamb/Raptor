//! Super-k-mer compression for memory-efficient k-mer storage
//!
//! A super-k-mer is a maximal sequence of consecutive k-mers that share
//! the same minimizer. Instead of storing each k-mer individually, we
//! store the super-k-mer sequence and can reconstruct all k-mers from it.
//!
//! Example with k=5, w=3:
//! Sequence: ACGTACGTACG
//! K-mers:   ACGTA, CGTAC, GTACG, TACGT, ACGTA, CGTAC, GTACG
//!
//! If ACGTA, CGTAC, GTACG share minimizer "ACG", they form a super-k-mer.
//! Storage: "ACGTACG" (7 bases) instead of 3 × 5 = 15 bases
//! Compression ratio: ~2-3x typical
//!
//! For large genomes, this reduces disk I/O by 50-70%.

use crate::kmer::nthash::NtHashIterator;

/// A super-k-mer: consecutive k-mers sharing the same minimizer
#[derive(Debug, Clone)]
pub struct SuperKmer {
    /// The sequence containing all k-mers (length = k + num_kmers - 1)
    pub sequence: Vec<u8>,
    /// The minimizer hash (used for bucket assignment)
    pub minimizer_hash: u64,
    /// Number of k-mers in this super-k-mer
    pub num_kmers: usize,
}

impl SuperKmer {
    /// Create a new super-k-mer starting with a single k-mer
    pub fn new(kmer: &[u8], minimizer_hash: u64) -> Self {
        Self {
            sequence: kmer.to_vec(),
            minimizer_hash,
            num_kmers: 1,
        }
    }

    /// Try to extend the super-k-mer with the next k-mer
    /// Returns true if extension was successful (same minimizer)
    pub fn try_extend(&mut self, next_base: u8, next_minimizer: u64) -> bool {
        if next_minimizer == self.minimizer_hash {
            self.sequence.push(next_base);
            self.num_kmers += 1;
            true
        } else {
            false
        }
    }

    /// Get the length of the underlying sequence
    pub fn sequence_len(&self) -> usize {
        self.sequence.len()
    }

    /// Iterate over all k-mers in this super-k-mer
    pub fn iter_kmers(&self, k: usize) -> impl Iterator<Item = &[u8]> {
        (0..self.num_kmers).map(move |i| &self.sequence[i..i + k])
    }

    /// Encode to compact binary format for disk storage
    /// Format: [minimizer_hash: 8 bytes][len: 2 bytes][sequence: len bytes]
    pub fn encode(&self) -> Vec<u8> {
        let mut encoded = Vec::with_capacity(10 + self.sequence.len());
        encoded.extend_from_slice(&self.minimizer_hash.to_le_bytes());
        encoded.extend_from_slice(&(self.sequence.len() as u16).to_le_bytes());
        encoded.extend_from_slice(&self.sequence);
        encoded
    }

    /// Decode from binary format
    pub fn decode(data: &[u8], k: usize) -> Option<Self> {
        if data.len() < 10 {
            return None;
        }

        let minimizer_hash = u64::from_le_bytes(data[0..8].try_into().ok()?);
        let seq_len = u16::from_le_bytes(data[8..10].try_into().ok()?) as usize;

        if data.len() < 10 + seq_len {
            return None;
        }

        let sequence = data[10..10 + seq_len].to_vec();
        let num_kmers = seq_len.saturating_sub(k - 1);

        Some(Self {
            sequence,
            minimizer_hash,
            num_kmers,
        })
    }

    /// Get encoded size in bytes
    pub fn encoded_size(&self) -> usize {
        10 + self.sequence.len()
    }
}

/// Extract super-k-mers from a sequence
pub struct SuperKmerExtractor {
    k: usize,
    w: usize,  // Minimizer window size
}

impl SuperKmerExtractor {
    pub fn new(k: usize, w: usize) -> Self {
        Self { k, w }
    }

    /// Extract all super-k-mers from a sequence
    pub fn extract(&self, seq: &[u8]) -> Vec<SuperKmer> {
        if seq.len() < self.k {
            return Vec::new();
        }

        let mut super_kmers = Vec::new();
        let mut current: Option<SuperKmer> = None;

        // Use ntHash for fast minimizer computation
        let hashes: Vec<u64> = NtHashIterator::new(seq, self.w)
            .map(|(_, h)| h)
            .collect();

        if hashes.is_empty() {
            return Vec::new();
        }

        // Find minimizer for each k-mer position
        for i in 0..=seq.len().saturating_sub(self.k) {
            // Get minimizer for this k-mer's window
            let window_start = i;
            let window_end = (i + self.k - self.w + 1).min(hashes.len());

            if window_start >= window_end {
                continue;
            }

            let minimizer = hashes[window_start..window_end]
                .iter()
                .min()
                .copied()
                .unwrap_or(0);

            match &mut current {
                None => {
                    // Start new super-k-mer
                    current = Some(SuperKmer::new(&seq[i..i + self.k], minimizer));
                }
                Some(sk) => {
                    // Try to extend existing super-k-mer
                    if minimizer == sk.minimizer_hash {
                        sk.sequence.push(seq[i + self.k - 1]);
                        sk.num_kmers += 1;
                    } else {
                        // Minimizer changed, save current and start new
                        super_kmers.push(std::mem::replace(
                            sk,
                            SuperKmer::new(&seq[i..i + self.k], minimizer),
                        ));
                    }
                }
            }
        }

        // Don't forget the last super-k-mer
        if let Some(sk) = current {
            super_kmers.push(sk);
        }

        super_kmers
    }

    /// Calculate compression ratio for a sequence
    pub fn compression_ratio(&self, seq: &[u8]) -> f64 {
        let super_kmers = self.extract(seq);

        let original_size: usize = super_kmers.iter().map(|sk| sk.num_kmers * self.k).sum();
        let compressed_size: usize = super_kmers.iter().map(|sk| sk.sequence_len()).sum();

        if compressed_size == 0 {
            1.0
        } else {
            original_size as f64 / compressed_size as f64
        }
    }
}

/// Efficient 2-bit encoding for super-k-mer sequences
pub struct TwoBitEncoder;

impl TwoBitEncoder {
    /// Encode sequence to 2-bit format (4 bases per byte)
    /// Returns None if sequence contains invalid bases
    pub fn encode(seq: &[u8]) -> Option<Vec<u8>> {
        let num_bytes = (seq.len() + 3) / 4;
        let mut encoded = vec![0u8; num_bytes];

        for (i, &base) in seq.iter().enumerate() {
            let bits = match base {
                b'A' | b'a' => 0b00,
                b'C' | b'c' => 0b01,
                b'G' | b'g' => 0b10,
                b'T' | b't' => 0b11,
                _ => return None,  // Invalid base
            };

            let byte_idx = i / 4;
            let bit_offset = (3 - (i % 4)) * 2;  // MSB first
            encoded[byte_idx] |= bits << bit_offset;
        }

        Some(encoded)
    }

    /// Decode 2-bit format back to sequence
    pub fn decode(encoded: &[u8], len: usize) -> Vec<u8> {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        let mut seq = Vec::with_capacity(len);

        for i in 0..len {
            let byte_idx = i / 4;
            let bit_offset = (3 - (i % 4)) * 2;
            let bits = (encoded[byte_idx] >> bit_offset) & 0b11;
            seq.push(BASES[bits as usize]);
        }

        seq
    }

    /// Calculate encoded size for a sequence of given length
    pub fn encoded_size(seq_len: usize) -> usize {
        (seq_len + 3) / 4
    }
}

/// Compact super-k-mer with 2-bit encoding for maximum compression
#[derive(Debug, Clone)]
pub struct CompactSuperKmer {
    /// 2-bit encoded sequence
    pub encoded_seq: Vec<u8>,
    /// Original sequence length
    pub seq_len: u16,
    /// Minimizer hash for bucket assignment
    pub minimizer_hash: u64,
}

impl CompactSuperKmer {
    pub fn from_super_kmer(sk: &SuperKmer) -> Option<Self> {
        let encoded_seq = TwoBitEncoder::encode(&sk.sequence)?;
        Some(Self {
            encoded_seq,
            seq_len: sk.sequence.len() as u16,
            minimizer_hash: sk.minimizer_hash,
        })
    }

    /// Get the sequence
    pub fn sequence(&self) -> Vec<u8> {
        TwoBitEncoder::decode(&self.encoded_seq, self.seq_len as usize)
    }

    /// Encode to binary format for disk storage
    /// Format: [minimizer: 8][seq_len: 2][encoded_seq: variable]
    pub fn encode_to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(10 + self.encoded_seq.len());
        bytes.extend_from_slice(&self.minimizer_hash.to_le_bytes());
        bytes.extend_from_slice(&self.seq_len.to_le_bytes());
        bytes.extend_from_slice(&self.encoded_seq);
        bytes
    }

    /// Calculate storage size
    pub fn storage_size(&self) -> usize {
        10 + self.encoded_seq.len()
    }

    /// Compression ratio compared to storing raw k-mers
    pub fn compression_ratio(&self, k: usize) -> f64 {
        let num_kmers = (self.seq_len as usize).saturating_sub(k - 1);
        let raw_size = num_kmers * k;
        if self.storage_size() == 0 {
            1.0
        } else {
            raw_size as f64 / self.storage_size() as f64
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_super_kmer_extraction() {
        let extractor = SuperKmerExtractor::new(5, 3);
        let seq = b"ACGTACGTACGTACGTACGT";

        let super_kmers = extractor.extract(seq);
        assert!(!super_kmers.is_empty());

        // Verify all k-mers can be reconstructed
        let mut reconstructed_count = 0;
        for sk in &super_kmers {
            reconstructed_count += sk.num_kmers;
        }
        assert_eq!(reconstructed_count, seq.len() - 5 + 1);
    }

    #[test]
    fn test_two_bit_encoding() {
        let seq = b"ACGTACGT";
        let encoded = TwoBitEncoder::encode(seq).unwrap();
        let decoded = TwoBitEncoder::decode(&encoded, seq.len());
        assert_eq!(seq.to_vec(), decoded);
    }

    #[test]
    fn test_compression_ratio() {
        let extractor = SuperKmerExtractor::new(31, 15);
        // Repetitive sequence should compress well
        let repetitive = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let ratio = extractor.compression_ratio(repetitive);
        assert!(ratio > 1.5, "Expected compression ratio > 1.5, got {}", ratio);
    }

    #[test]
    fn test_super_kmer_encode_decode() {
        let sk = SuperKmer::new(b"ACGTACGT", 12345);
        let encoded = sk.encode();
        let decoded = SuperKmer::decode(&encoded, 8).unwrap();
        assert_eq!(sk.sequence, decoded.sequence);
        assert_eq!(sk.minimizer_hash, decoded.minimizer_hash);
    }
}
