// GPU kernel for parallel overlap detection between contigs
//
// This kernel computes overlaps between all pairs of contigs in parallel,
// transforming O(n²) sequential comparisons to massively parallel execution.

// Calculate Hamming distance between two sequences
inline uint hamming_distance(
    __global const uchar* seq1, int start1,
    __global const uchar* seq2, int start2,
    int len
) {
    uint dist = 0;
    for (int i = 0; i < len; i++) {
        if (seq1[start1 + i] != seq2[start2 + i]) {
            dist++;
        }
    }
    return dist;
}

// Main 2D parallel overlap detection kernel
// Each work item (i, j) checks if contig i overlaps with contig j
__kernel void find_overlaps_2d(
    __global const uchar* sequences,     // Flattened contig sequences
    __global const int* offsets,         // Start offset for each contig
    __global const int* lengths,         // Length of each contig
    __global int* overlap_from,          // Output: source contig index
    __global int* overlap_to,            // Output: target contig index
    __global int* overlap_len,           // Output: overlap length
    __global int* overlap_count,         // Atomic counter for results
    const uint num_contigs,              // Number of contigs
    const uint min_overlap,              // Minimum overlap length
    const uint max_mismatches,           // Maximum allowed mismatches
    const uint max_results               // Maximum number of results to store
) {
    uint i = get_global_id(0);  // From contig index
    uint j = get_global_id(1);  // To contig index

    // Skip self-comparison and out-of-bounds
    if (i >= num_contigs || j >= num_contigs || i == j) {
        return;
    }

    int len_i = lengths[i];
    int len_j = lengths[j];
    int start_i = offsets[i];
    int start_j = offsets[j];

    // Try to find suffix-prefix overlap: suffix of i with prefix of j
    int max_possible_overlap = min(len_i, len_j);
    if (max_possible_overlap < (int)min_overlap) {
        return;
    }

    // Try different overlap lengths, starting from the largest
    for (int overlap = max_possible_overlap; overlap >= (int)min_overlap; overlap--) {
        int suffix_start = start_i + len_i - overlap;  // Start of suffix in contig i
        int prefix_start = start_j;                     // Start of prefix in contig j

        uint dist = hamming_distance(
            sequences, suffix_start,
            sequences, prefix_start,
            overlap
        );

        if (dist <= max_mismatches) {
            // Found a valid overlap
            int result_idx = atomic_inc(overlap_count);

            if (result_idx < (int)max_results) {
                overlap_from[result_idx] = i;
                overlap_to[result_idx] = j;
                overlap_len[result_idx] = overlap;
            }

            // Found the best (longest) overlap for this pair
            return;
        }
    }
}

// Alternative: 1D linearized kernel for smaller datasets
// Each work item processes one contig pair
__kernel void find_overlaps_1d(
    __global const uchar* sequences,
    __global const int* offsets,
    __global const int* lengths,
    __global int* overlap_from,
    __global int* overlap_to,
    __global int* overlap_len,
    __global int* overlap_count,
    const uint num_contigs,
    const uint min_overlap,
    const uint max_mismatches,
    const uint max_results
) {
    uint gid = get_global_id(0);

    // Convert linear index to pair indices
    // For n contigs, we have n*(n-1) pairs (excluding self-pairs)
    uint total_pairs = num_contigs * (num_contigs - 1);
    if (gid >= total_pairs) return;

    // Decode pair indices from linear index
    uint i = gid / (num_contigs - 1);
    uint j_offset = gid % (num_contigs - 1);
    uint j = (j_offset < i) ? j_offset : j_offset + 1;  // Skip self-comparison

    int len_i = lengths[i];
    int len_j = lengths[j];
    int start_i = offsets[i];
    int start_j = offsets[j];

    int max_possible_overlap = min(len_i, len_j);
    if (max_possible_overlap < (int)min_overlap) {
        return;
    }

    for (int overlap = max_possible_overlap; overlap >= (int)min_overlap; overlap--) {
        int suffix_start = start_i + len_i - overlap;
        int prefix_start = start_j;

        uint dist = hamming_distance(
            sequences, suffix_start,
            sequences, prefix_start,
            overlap
        );

        if (dist <= max_mismatches) {
            int result_idx = atomic_inc(overlap_count);

            if (result_idx < (int)max_results) {
                overlap_from[result_idx] = i;
                overlap_to[result_idx] = j;
                overlap_len[result_idx] = overlap;
            }
            return;
        }
    }
}

// Batch overlap detection with SIMD-style vectorization hint
// Processes 4 overlap lengths in parallel per work item
__kernel void find_overlaps_batched(
    __global const uchar* sequences,
    __global const int* offsets,
    __global const int* lengths,
    __global int* overlap_from,
    __global int* overlap_to,
    __global int* overlap_len,
    __global int* overlap_count,
    const uint num_contigs,
    const uint min_overlap,
    const uint max_mismatches,
    const uint max_results
) {
    uint pair_idx = get_global_id(0);
    uint batch_idx = get_global_id(1);  // Which batch of overlap lengths to try

    uint total_pairs = num_contigs * (num_contigs - 1);
    if (pair_idx >= total_pairs) return;

    uint i = pair_idx / (num_contigs - 1);
    uint j_offset = pair_idx % (num_contigs - 1);
    uint j = (j_offset < i) ? j_offset : j_offset + 1;

    int len_i = lengths[i];
    int len_j = lengths[j];
    int start_i = offsets[i];
    int start_j = offsets[j];

    int max_possible_overlap = min(len_i, len_j);

    // Each batch handles 4 overlap lengths
    int batch_size = 4;
    int overlap_start = max_possible_overlap - batch_idx * batch_size;
    int overlap_end = max(overlap_start - batch_size, (int)min_overlap - 1);

    for (int overlap = overlap_start; overlap > overlap_end; overlap--) {
        if (overlap < (int)min_overlap) break;

        int suffix_start = start_i + len_i - overlap;
        int prefix_start = start_j;

        uint dist = hamming_distance(
            sequences, suffix_start,
            sequences, prefix_start,
            overlap
        );

        if (dist <= max_mismatches) {
            int result_idx = atomic_inc(overlap_count);

            if (result_idx < (int)max_results) {
                overlap_from[result_idx] = i;
                overlap_to[result_idx] = j;
                overlap_len[result_idx] = overlap;
            }
            return;
        }
    }
}
