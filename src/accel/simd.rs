//! SIMD-accelerated k-mer operations
#![allow(dead_code)]

use std::arch::x86_64::*;
use std::cmp;
use rayon::prelude::*;

/// Check if AVX2 is available at runtime
#[inline]
fn is_avx2_available() -> bool {
    #[cfg(target_feature = "avx2")]
    {
        true
    }
    #[cfg(not(target_feature = "avx2"))]
    {
        is_x86_feature_detected!("avx2")
    }
}

/// Calculate Hamming distance using AVX2 (32 bytes at a time).
///
/// This provides 2x speedup over SSE2 for sequences >= 32 bytes.
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn hamming_distance_avx2_inner(a: &[u8], b: &[u8]) -> usize {
    let common_length = std::cmp::min(a.len(), b.len());
    let mut dist = 0;
    let mut i = 0;

    // Process 32 bytes at a time with AVX2
    while i + 32 <= common_length {
        let a_chunk = _mm256_loadu_si256(a.as_ptr().add(i) as *const __m256i);
        let b_chunk = _mm256_loadu_si256(b.as_ptr().add(i) as *const __m256i);
        let cmp = _mm256_cmpeq_epi8(a_chunk, b_chunk);
        let mask = _mm256_movemask_epi8(cmp) as u32;
        // Count mismatches (bits that are 0 in the mask)
        dist += 32 - mask.count_ones() as usize;
        i += 32;
    }

    // Process remaining 16 bytes with SSE2
    while i + 16 <= common_length {
        let a_chunk = _mm_loadu_si128(a.as_ptr().add(i) as *const __m128i);
        let b_chunk = _mm_loadu_si128(b.as_ptr().add(i) as *const __m128i);
        let cmp = _mm_cmpeq_epi8(a_chunk, b_chunk);
        let mask = _mm_movemask_epi8(cmp);
        dist += 16 - mask.count_ones() as usize;
        i += 16;
    }

    // Scalar fallback for remaining bytes
    for j in i..common_length {
        if a[j] != b[j] {
            dist += 1;
        }
    }

    // Add length difference as mismatches
    dist += a.len().abs_diff(b.len());
    dist
}

/// Calculate Hamming distance using SSE2 (16 bytes at a time).
#[inline]
unsafe fn hamming_distance_sse2_inner(a: &[u8], b: &[u8]) -> usize {
    let common_length = std::cmp::min(a.len(), b.len());
    let mut dist = 0;
    let mut i = 0;

    // Process 16 bytes at a time
    while i + 16 <= common_length {
        let a_chunk = _mm_loadu_si128(a.as_ptr().add(i) as *const __m128i);
        let b_chunk = _mm_loadu_si128(b.as_ptr().add(i) as *const __m128i);
        let cmp = _mm_cmpeq_epi8(a_chunk, b_chunk);
        let mask = _mm_movemask_epi8(cmp);
        dist += 16 - mask.count_ones() as usize;
        i += 16;
    }

    // Scalar fallback
    for j in i..common_length {
        if a[j] != b[j] {
            dist += 1;
        }
    }

    dist += a.len().abs_diff(b.len());
    dist
}

/// Calculate Hamming distance (number of mismatching positions) between two byte slices.
///
/// Automatically selects the best SIMD implementation:
/// - AVX2 (32 bytes/iteration) if available
/// - SSE2 (16 bytes/iteration) as fallback
/// - Scalar for remaining bytes
#[inline]
pub fn hamming_distance_simd(a: &[u8], b: &[u8]) -> usize {
    // For short sequences, scalar is faster (no SIMD setup overhead)
    if a.len() < 16 && b.len() < 16 {
        let common_length = std::cmp::min(a.len(), b.len());
        let mut dist = 0;
        for i in 0..common_length {
            if a[i] != b[i] {
                dist += 1;
            }
        }
        return dist + a.len().abs_diff(b.len());
    }

    unsafe {
        if is_avx2_available() {
            hamming_distance_avx2_inner(a, b)
        } else {
            hamming_distance_sse2_inner(a, b)
        }
    }
}

/// Returns the reverse complement of a DNA sequence
pub fn reverse_complement_simd(seq: &str) -> String {
    let bytes = seq.as_bytes();
    let mut result = Vec::with_capacity(bytes.len());
    
    // Process the sequence in reverse order
    for i in (0..bytes.len()).rev() {
        // Get the complement of the current nucleotide
        let complement = match bytes[i] {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N', // Non-standard bases are replaced with N
        };
        
        result.push(complement);
    }
    
    // Convert back to a UTF-8 string
    String::from_utf8(result).unwrap()
}

/// Match k-mers with potential mismatches using SIMD acceleration
pub fn match_kmers_simd(a: &[u8], b: &[u8], max_mismatches: usize) -> bool {
    if a.len() != b.len() {
        return false;
    }
    
    let dist = hamming_distance_simd(a, b);
    dist <= max_mismatches
}

/// Find overlap between sequences with a maximum allowed distance
pub fn match_kmers_with_overlap(query: &str, target: &str, min_overlap: usize, max_distance: usize) -> Option<(usize, usize)> {
    let q = query.as_bytes();
    let t = target.as_bytes();

    let max_shift = cmp::min(q.len(), t.len()) - min_overlap;

    for shift in 0..=max_shift {
        let q_sub = &q[shift..];
        let t_sub = &t[..cmp::min(q_sub.len(), t.len())];
        
        if t_sub.len() < min_overlap {
            continue;
        }
        
        let overlap_len = cmp::min(q_sub.len(), t_sub.len());
        let q_overlap = &q_sub[..overlap_len];
        let t_overlap = &t_sub[..overlap_len];

        let dist = hamming_distance_simd(q_overlap, t_overlap);
        if dist <= max_distance {
            return Some((shift, dist));
        }
    }

    None
}

/// Process a batch of sequences looking for overlaps with a query sequence
/// Returns a vector of matches with (candidate, shift, distance)
pub fn batch_match_kmers_simd<'a>(
    query: &[u8],
    candidates: &[&[u8]],
    min_overlap: usize,
    max_mismatch: usize,
) -> Vec<(Option<(usize, usize)>, usize)> {
    // Process candidates in parallel for better performance
    candidates.par_iter()
        .map(|candidate| {
            // Find the best overlap between query and this candidate
            find_best_overlap(query, candidate, min_overlap, max_mismatch)
        })
        .collect()
}

/// Find the best overlap between two sequences (as byte arrays)
/// Returns the best overlap position, length, and mismatch count
#[inline]
fn find_best_overlap(query: &[u8], target: &[u8], min_overlap: usize, max_distance: usize) -> (Option<(usize, usize)>, usize) {
    // If either sequence is shorter than min_overlap, no match is possible
    if query.len() < min_overlap || target.len() < min_overlap {
        return (None, 0);
    }

    let max_shift = cmp::min(query.len(), target.len()) - min_overlap;
    let mut best_overlap: Option<(usize, usize)> = None;
    let mut best_distance = max_distance + 1; // Start with worse than max allowed

    // Try overlaps with query at the beginning (target shifted right)
    for shift in 0..=max_shift {
        let query_sub = &query[..cmp::min(query.len(), target.len() - shift)];
        let target_sub = &target[shift..shift + query_sub.len()];
        
        if query_sub.len() < min_overlap {
            continue;
        }
        
        let dist = hamming_distance_simd(query_sub, target_sub);
        if dist <= max_distance && (best_overlap.is_none() || dist < best_distance) {
            best_overlap = Some((query_sub.len(), shift));
            best_distance = dist;
        }
    }

    // Try overlaps with target at the beginning (query shifted right)
    for shift in 1..=max_shift {
        let query_sub = &query[shift..];
        let target_sub = &target[..cmp::min(target.len(), query_sub.len())];
        
        if target_sub.len() < min_overlap {
            continue;
        }
        
        let dist = hamming_distance_simd(query_sub, target_sub);
        if dist <= max_distance && (best_overlap.is_none() || dist < best_distance) {
            best_overlap = Some((target_sub.len(), shift + target.len()));  // Encode that this is a different type of overlap
            best_distance = dist;
        }
    }

    (best_overlap, best_distance)
}

/// Bit-parallel edit distance (Myers' algorithm)
pub fn edit_distance_bp(a: &str, b: &str, max_dist: usize) -> Option<usize> {
    let m = a.len();
    if m > 64 { return None; } // Limited to 64-bit

    let mut pv = !0u64;
    let mut mv = 0u64;
    let mut score = m;

    let mut peq = [0u64; 256];
    for (i, &b) in a.as_bytes().iter().enumerate() {
        peq[b as usize] |= 1 << i;
    }

    for &b in b.as_bytes() {
        let eq = peq[b as usize];
        let _xv = eq | mv;
        let xh = (((eq & pv).wrapping_add(pv)) ^ pv) | eq;

        let ph = mv | !(xh | pv);
        let mh = pv & xh;

        if (ph & (1 << (m - 1))) != 0 { score += 1; }
        if (mh & (1 << (m - 1))) != 0 { score -= 1; }

        pv = (mh << 1) | 1;
        mv = ph << 1;
    }

    if score <= max_dist { Some(score) } else { None }
}
