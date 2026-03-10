// src/kmer/normalize.rs
use crate::io::fastq::FastqRecord;
use crate::kmer::cms::CountMinSketch;
use crate::kmer::nthash::NtHashIterator;

/// Statistics about a read's k-mer abundance
#[derive(Debug, Clone, Copy)]
pub struct AbundanceStats {
    pub median: u32,
    pub min: u32,
}

/// Estimates the abundance statistics for a read sequence.
/// Uses ntHash for fast O(1) rolling hash computation.
pub fn estimate_read_abundance(seq: &str, k: usize, sketch: &[u32]) -> AbundanceStats {
    let bytes = seq.as_bytes();
    let sketch_mask = (sketch.len() - 1) as u64;

    let mut abund: Vec<u32> = NtHashIterator::new(bytes, k)
        .map(|(_, hash)| sketch[(hash & sketch_mask) as usize])
        .collect();

    if abund.is_empty() {
        return AbundanceStats { median: 0, min: 0 };
    }

    let mid = abund.len() / 2;
    abund.select_nth_unstable(mid);
    AbundanceStats {
        median: abund[mid],
        min: *abund.iter().min().unwrap(),
    }
}

/// Determines whether a read should be kept based on its k-mer coverage.
/// Uses ntHash for fast O(1) rolling hash computation.
pub fn should_keep_read(
    record: &FastqRecord,
    cms: &CountMinSketch,
    k: usize,
    target: u16,
    min_abund: u16,
) -> bool {
    if record.sequence.len() < k {
        return false; // Skip reads shorter than k
    }

    let bytes = record.sequence.as_bytes();

    // Get median abundance of k-mers in the read using ntHash
    let mut abundances: Vec<u16> = NtHashIterator::new(bytes, k)
        .map(|(_, hash)| cms.estimate_hash(hash))
        .collect();

    // If no valid k-mers, skip this read
    if abundances.is_empty() {
        return false;
    }

    let median_idx = abundances.len() / 2;
    let (_, median_abund, _) = abundances.select_nth_unstable(median_idx);
    let median_abund = *median_abund;

    // Skip reads with low abundance (potential errors)
    if median_abund < min_abund {
        return false;
    }

    // Keep high-abundance reads with probability target/abundance
    if median_abund > target {
        let keep_prob = target as f64 / median_abund as f64;
        return deterministic_roll(record) < keep_prob;
    }

    // Always keep reads at or below target abundance
    true
}

/// Determines whether a read pair should be kept based on both reads' k-mer coverage
pub fn should_keep_read_pair(
    r1: &FastqRecord,
    r2: &FastqRecord,
    cms: &CountMinSketch,
    k: usize,
    target: u16,
    min_abund: u16,
) -> bool {
    // Keep a pair only if both reads should be kept
    should_keep_read(r1, cms, k, target, min_abund)
        && should_keep_read(r2, cms, k, target, min_abund)
}

#[inline]
fn deterministic_roll(record: &FastqRecord) -> f64 {
    const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    const FNV_PRIME: u64 = 0x0000_0001_0000_01b3;

    fn mix(mut hash: u64, bytes: &[u8]) -> u64 {
        for &b in bytes {
            hash ^= b as u64;
            hash = hash.wrapping_mul(FNV_PRIME);
        }
        hash
    }

    let mut hash = FNV_OFFSET;
    hash = mix(hash, record.header.as_bytes());
    hash = hash.wrapping_mul(FNV_PRIME) ^ 0xff;
    hash = mix(hash, record.sequence.as_bytes());
    hash = hash.wrapping_mul(FNV_PRIME) ^ 0xfe;
    hash = mix(hash, record.quality.as_bytes());

    // Map top 53 bits to [0, 1) for stable float comparisons.
    const INV_2POW53: f64 = 1.0 / ((1u64 << 53) as f64);
    ((hash >> 11) as f64) * INV_2POW53
}

#[cfg(test)]
mod tests {
    use super::should_keep_read;
    use crate::io::fastq::FastqRecord;
    use crate::kmer::cms::CountMinSketch;
    use crate::kmer::nthash::NtHashIterator;

    fn build_record(header: &str, sequence: &str) -> FastqRecord {
        FastqRecord {
            header: header.to_string(),
            sequence: sequence.to_string(),
            plus: "+".to_string(),
            quality: "I".repeat(sequence.len()),
        }
    }

    fn build_cms_with_depth(
        sequence: &str,
        k: usize,
        depth: usize,
        repeats: usize,
    ) -> CountMinSketch {
        let mut cms = CountMinSketch::new(depth, 1 << 15);
        for _ in 0..repeats {
            for (_, hash) in NtHashIterator::new(sequence.as_bytes(), k) {
                cms.insert_hash(hash);
            }
        }
        cms
    }

    #[test]
    fn downsampling_decision_is_deterministic_per_read() {
        let record = build_record("@read_42", "ACGTACGTACGTACGTACGTACGTACGT");
        let cms = build_cms_with_depth(&record.sequence, 7, 4, 120);

        let baseline = should_keep_read(&record, &cms, 7, 20, 2);
        for _ in 0..64 {
            assert_eq!(should_keep_read(&record, &cms, 7, 20, 2), baseline);
        }
    }

    #[test]
    fn deterministic_downsampling_preserves_expected_fraction_across_many_reads() {
        let sequence = "TGCATGCATGCATGCATGCATGCATGCA";
        let cms = build_cms_with_depth(sequence, 7, 4, 100);
        let target = 10;
        let min_abund = 2;
        let mut abundances: Vec<u16> = NtHashIterator::new(sequence.as_bytes(), 7)
            .map(|(_, hash)| cms.estimate_hash(hash))
            .collect();
        let mid = abundances.len() / 2;
        let (_, median_abund, _) = abundances.select_nth_unstable(mid);
        let expected = target as f64 / *median_abund as f64;

        let mut kept = 0usize;
        let mut total = 0usize;
        for i in 0..2000 {
            let header = format!("@read_{i}");
            let record = build_record(&header, sequence);
            if should_keep_read(&record, &cms, 7, target, min_abund) {
                kept += 1;
            }
            total += 1;
        }

        let observed = kept as f64 / total as f64;
        assert!((observed - expected).abs() < 0.05);
    }

    #[test]
    fn short_reads_are_rejected() {
        let record = build_record("@short", "ACGT");
        let cms = CountMinSketch::new(4, 1024);
        assert!(!should_keep_read(&record, &cms, 7, 20, 2));
    }
}
