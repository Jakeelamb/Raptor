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
    if k == 0 || seq.len() < k || sketch.is_empty() {
        return AbundanceStats { median: 0, min: 0 };
    }

    let bytes = seq.as_bytes();
    let mut abund: Vec<u32> = Vec::with_capacity(bytes.len() - k + 1);
    let mut min_abundance = u32::MAX;

    if sketch.len().is_power_of_two() {
        let sketch_mask = (sketch.len() - 1) as u64;
        for (_, hash) in NtHashIterator::new(bytes, k) {
            let count = sketch[(hash & sketch_mask) as usize];
            min_abundance = min_abundance.min(count);
            abund.push(count);
        }
    } else {
        let sketch_len = sketch.len() as u64;
        for (_, hash) in NtHashIterator::new(bytes, k) {
            let count = sketch[(hash % sketch_len) as usize];
            min_abundance = min_abundance.min(count);
            abund.push(count);
        }
    }

    if abund.is_empty() {
        return AbundanceStats { median: 0, min: 0 };
    }

    let mid = abund.len() / 2;
    abund.select_nth_unstable(mid);
    AbundanceStats {
        median: abund[mid],
        min: min_abundance,
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
    let mut scratch = Vec::new();
    should_keep_read_with_scratch(record, cms, k, target, min_abund, &mut scratch)
}

/// Determines whether a read should be kept based on its k-mer coverage
/// using a caller-provided scratch buffer to avoid per-read allocations.
pub fn should_keep_read_with_scratch(
    record: &FastqRecord,
    cms: &CountMinSketch,
    k: usize,
    target: u16,
    min_abund: u16,
    scratch: &mut Vec<u16>,
) -> bool {
    if record.sequence.len() < k {
        return false; // Skip reads shorter than k
    }

    let bytes = record.sequence.as_bytes();

    // Get median abundance of k-mers in the read using ntHash
    scratch.clear();
    scratch.extend(NtHashIterator::new(bytes, k).map(|(_, hash)| cms.estimate_hash(hash)));

    // If no valid k-mers, skip this read
    if scratch.is_empty() {
        return false;
    }

    let median_idx = scratch.len() / 2;
    let (_, median_abund, _) = scratch.select_nth_unstable(median_idx);
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
    let mut scratch_r1 = Vec::new();
    let mut scratch_r2 = Vec::new();
    should_keep_read_pair_with_scratch(
        r1,
        r2,
        cms,
        k,
        target,
        min_abund,
        &mut scratch_r1,
        &mut scratch_r2,
    )
}

/// Determines whether a read pair should be kept based on both reads' k-mer
/// coverage while reusing caller-provided scratch buffers.
pub fn should_keep_read_pair_with_scratch(
    r1: &FastqRecord,
    r2: &FastqRecord,
    cms: &CountMinSketch,
    k: usize,
    target: u16,
    min_abund: u16,
    scratch_r1: &mut Vec<u16>,
    scratch_r2: &mut Vec<u16>,
) -> bool {
    // Keep a pair only if both reads should be kept
    should_keep_read_with_scratch(r1, cms, k, target, min_abund, scratch_r1)
        && should_keep_read_with_scratch(r2, cms, k, target, min_abund, scratch_r2)
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
    use super::{
        estimate_read_abundance, should_keep_read, should_keep_read_pair,
        should_keep_read_pair_with_scratch,
    };
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

    fn reference_abundance(seq: &str, k: usize, sketch: &[u32]) -> (u32, u32) {
        if k == 0 || seq.len() < k || sketch.is_empty() {
            return (0, 0);
        }

        let sketch_len = sketch.len() as u64;
        let mut values: Vec<u32> = NtHashIterator::new(seq.as_bytes(), k)
            .map(|(_, hash)| sketch[(hash % sketch_len) as usize])
            .collect();

        if values.is_empty() {
            return (0, 0);
        }

        let min = values.iter().copied().min().unwrap_or(0);
        let mid = values.len() / 2;
        values.select_nth_unstable(mid);
        (values[mid], min)
    }

    #[test]
    fn estimate_read_abundance_handles_empty_sketch_and_invalid_k() {
        let seq = "ACGTACGT";
        let empty: Vec<u32> = Vec::new();

        let zero_k = estimate_read_abundance(seq, 0, &[1, 2, 3, 4]);
        assert_eq!(zero_k.median, 0);
        assert_eq!(zero_k.min, 0);

        let short_seq = estimate_read_abundance("ACG", 5, &[1, 2, 3, 4]);
        assert_eq!(short_seq.median, 0);
        assert_eq!(short_seq.min, 0);

        let empty_sketch = estimate_read_abundance(seq, 3, &empty);
        assert_eq!(empty_sketch.median, 0);
        assert_eq!(empty_sketch.min, 0);
    }

    #[test]
    fn estimate_read_abundance_matches_reference_on_non_power_of_two_sketch() {
        let seq = "ACGTACGTACGT";
        let sketch = vec![9, 2, 5, 1, 4, 7, 3];

        let observed = estimate_read_abundance(seq, 5, &sketch);
        let (expected_median, expected_min) = reference_abundance(seq, 5, &sketch);

        assert_eq!(observed.median, expected_median);
        assert_eq!(observed.min, expected_min);
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

    #[test]
    fn scratch_pair_decision_matches_non_scratch_path() {
        let r1 = build_record("@r1", "ACGTACGTACGTACGTACGTACGTACGT");
        let r2 = build_record("@r2", "TGCATGCATGCATGCATGCATGCATGCA");
        let mut cms = CountMinSketch::new(4, 1 << 15);
        for _ in 0..64 {
            for (_, hash) in NtHashIterator::new(r1.sequence.as_bytes(), 7) {
                cms.insert_hash(hash);
            }
            for (_, hash) in NtHashIterator::new(r2.sequence.as_bytes(), 7) {
                cms.insert_hash(hash);
            }
        }

        let baseline = should_keep_read_pair(&r1, &r2, &cms, 7, 20, 2);
        let mut scratch_r1 = Vec::new();
        let mut scratch_r2 = Vec::new();
        let observed = should_keep_read_pair_with_scratch(
            &r1,
            &r2,
            &cms,
            7,
            20,
            2,
            &mut scratch_r1,
            &mut scratch_r2,
        );
        assert_eq!(observed, baseline);
    }
}
