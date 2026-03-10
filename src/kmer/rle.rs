/// Run-length encode a DNA string: AATT → [(A,2), (T,2)]
pub fn rle_encode(seq: &str) -> Vec<(u8, u8)> {
    if seq.is_empty() {
        return Vec::new();
    }

    let mut result = Vec::new();
    let mut bytes = seq.bytes();
    let mut current = bytes.next().expect("empty sequences are handled above");
    let mut count: u8 = 1;

    for next in bytes {
        if next == current {
            if count == u8::MAX {
                // Split very long homopolymers into multiple chunks to avoid overflow.
                result.push((current, u8::MAX));
                count = 1;
            } else {
                count += 1;
            }
        } else {
            result.push((current, count));
            current = next;
            count = 1;
        }
    }
    result.push((current, count));

    result
}

/// Return the number of RLE entries that `rle_encode(seq)` would produce,
/// without allocating the encoded vector.
pub fn rle_encoded_len(seq: &str) -> usize {
    let bytes = seq.as_bytes();
    if bytes.is_empty() {
        return 0;
    }

    let mut runs = 0usize;
    let mut current = bytes[0];
    let mut run_len = 1usize;

    for &base in &bytes[1..] {
        if base == current {
            run_len += 1;
        } else {
            runs += run_len.div_ceil(u8::MAX as usize);
            current = base;
            run_len = 1;
        }
    }
    runs + run_len.div_ceil(u8::MAX as usize)
}

/// Decode RLE back to a DNA string
pub fn rle_decode(encoded: &[(u8, u8)]) -> String {
    encoded
        .iter()
        .flat_map(|(b, c)| std::iter::repeat(*b).take(*c as usize))
        .map(|b| b as char)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rle_encoding() {
        let seq = "AAAACCCGGGGTTT";
        let encoded = rle_encode(seq);
        assert_eq!(encoded, vec![(b'A', 4), (b'C', 3), (b'G', 4), (b'T', 3)]);

        let decoded = rle_decode(&encoded);
        assert_eq!(decoded, seq);
    }

    #[test]
    fn test_rle_single_bases() {
        let seq = "ACGT";
        let encoded = rle_encode(seq);
        assert_eq!(encoded, vec![(b'A', 1), (b'C', 1), (b'G', 1), (b'T', 1)]);

        let decoded = rle_decode(&encoded);
        assert_eq!(decoded, seq);
    }

    #[test]
    fn test_rle_empty() {
        let seq = "";
        let encoded = rle_encode(seq);
        assert_eq!(encoded, vec![]);

        let decoded = rle_decode(&encoded);
        assert_eq!(decoded, seq);
    }

    #[test]
    fn test_rle_long_run_splits_at_u8_max_without_overflow() {
        let seq = "A".repeat(256 + 44);
        let encoded = rle_encode(&seq);

        assert_eq!(encoded, vec![(b'A', 255), (b'A', 45)]);
        assert_eq!(rle_encoded_len(&seq), encoded.len());
        let decoded = rle_decode(&encoded);
        assert_eq!(decoded, seq);
    }

    #[test]
    fn test_rle_encoded_len_matches_materialized_encoding() {
        let cases = [
            "",
            "A",
            "ACGT",
            "AAAACCCGGGGTTT",
            "AAATTTCCCAAAGGG",
            "AACCGGTTNNNN",
        ];

        for seq in cases {
            assert_eq!(rle_encoded_len(seq), rle_encode(seq).len());
        }

        let long = format!("{}{}{}", "A".repeat(700), "C".repeat(2), "G".repeat(512));
        assert_eq!(rle_encoded_len(&long), rle_encode(&long).len());
    }
}
