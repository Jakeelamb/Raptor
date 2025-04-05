/// Run-length encode a DNA string: AATT → [(A,2), (T,2)]
pub fn rle_encode(seq: &str) -> Vec<(u8, u8)> {
    let mut result = Vec::new();
    let mut chars = seq.bytes().peekable();

    while let Some(b) = chars.next() {
        let mut count = 1;
        while let Some(&next) = chars.peek() {
            if next == b {
                count += 1;
                chars.next();
            } else {
                break;
            }
        }
        result.push((b, count));
    }

    result
}

/// Decode RLE back to a DNA string
pub fn rle_decode(encoded: &[(u8, u8)]) -> String {
    encoded.iter().flat_map(|(b, c)| std::iter::repeat(*b).take(*c as usize)).map(|b| b as char).collect()
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
}
