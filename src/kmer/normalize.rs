// src/kmer/normalize.rs
use crate::io::fastq::FastqRecord;
use crate::kmer::cms::CountMinSketch;
use crate::kmer::nthash::NtHashIterator;
use rand::Rng;

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

    // Sort abundances and find median
    abundances.sort_unstable();
    let median_idx = abundances.len() / 2;
    let median_abund = abundances[median_idx];

    // Skip reads with low abundance (potential errors)
    if median_abund < min_abund {
        return false;
    }

    // Keep high-abundance reads with probability target/abundance
    if median_abund > target {
        let mut rng = rand::thread_rng();
        let keep_prob = target as f64 / median_abund as f64;
        return rng.gen::<f64>() < keep_prob;
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
    should_keep_read(r1, cms, k, target, min_abund) && 
    should_keep_read(r2, cms, k, target, min_abund)
}

fn median(values: &mut [u16]) -> u16 {
    values.sort_unstable();
    let mid = values.len() / 2;
    if values.len() % 2 == 0 {
        // Even number of elements, take average of middle two
        (values[mid - 1] + values[mid]) / 2
    } else {
        // Odd number of elements, take middle one
        values[mid]
    }
}
