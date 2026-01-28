// Improved GPU k-mer counting kernel with:
// - MurmurHash for better distribution
// - Canonical k-mer handling (reverse complement)
// - Proper N-base handling
// - Configurable hash table size

// Base encoding: A=0, C=1, G=2, T=3, N=255 (invalid)
inline uchar encode_base(uchar base) {
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 255; // Invalid/N base
    }
}

// Complement: A<->T, C<->G
inline uchar complement_base(uchar encoded) {
    return encoded ^ 3;
}

// MurmurHash3 finalizer for 64-bit
inline ulong murmur_finalize(ulong h) {
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdUL;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53UL;
    h ^= h >> 33;
    return h;
}

// Encode k-mer to 64-bit value (2 bits per base)
// Returns 0xFFFFFFFFFFFFFFFF if k-mer contains invalid bases
inline ulong encode_kmer(__global const uchar* seq, int start, int k) {
    ulong encoded = 0;
    for (int i = 0; i < k; i++) {
        uchar base = encode_base(seq[start + i]);
        if (base == 255) {
            return 0xFFFFFFFFFFFFFFFFUL; // Invalid k-mer
        }
        encoded = (encoded << 2) | base;
    }
    return encoded;
}

// Compute reverse complement of encoded k-mer
inline ulong reverse_complement(ulong encoded, int k) {
    ulong result = 0;
    for (int i = 0; i < k; i++) {
        uchar base = encoded & 3;
        result = (result << 2) | complement_base(base);
        encoded >>= 2;
    }
    return result;
}

// Get canonical k-mer (min of kmer and reverse complement)
inline ulong canonical_kmer(ulong encoded, int k) {
    ulong rc = reverse_complement(encoded, k);
    return (encoded < rc) ? encoded : rc;
}

// Main kernel: count k-mers from sequences
// Each work item processes one read
__kernel void count_kmers_v2(
    __global const uchar* sequences,    // Flattened sequence data
    __global const int* offsets,        // Start offset for each read
    __global uint* hash_table_counts,   // Hash table for counts
    __global ulong* hash_table_kmers,   // Hash table for k-mer values (for collision resolution)
    const uint k,                       // K-mer size
    const uint num_reads,               // Number of reads
    const uint table_size,              // Hash table size (should be power of 2)
    const uint table_mask               // table_size - 1 for fast modulo
) {
    uint gid = get_global_id(0);
    if (gid >= num_reads) return;

    int start = offsets[gid];
    int end = offsets[gid + 1];
    int len = end - start;

    if (len < (int)k) return;

    // Process each k-mer position
    for (int i = 0; i <= len - (int)k; i++) {
        ulong encoded = encode_kmer(sequences, start + i, k);

        // Skip k-mers with invalid bases
        if (encoded == 0xFFFFFFFFFFFFFFFFUL) continue;

        // Get canonical form
        ulong canon = canonical_kmer(encoded, k);

        // Hash the canonical k-mer
        ulong hash = murmur_finalize(canon);
        uint slot = (uint)(hash & table_mask);

        // Linear probing for collision resolution
        // Note: This is a simple approach; more sophisticated methods
        // could be used for better performance
        for (uint probe = 0; probe < 32; probe++) {
            uint idx = (slot + probe) & table_mask;

            // Try to claim this slot
            ulong expected = 0;
            ulong old = atom_cmpxchg(&hash_table_kmers[idx], expected, canon);

            if (old == 0 || old == canon) {
                // Successfully claimed or already our k-mer
                atomic_inc(&hash_table_counts[idx]);
                break;
            }
            // Collision with different k-mer, try next slot
        }
    }
}

// Alternative kernel for small datasets: direct atomic increment
// Uses simpler hash without collision resolution
__kernel void count_kmers_simple(
    __global const uchar* sequences,
    __global const int* offsets,
    __global uint* counts,
    const uint k,
    const uint num_reads,
    const uint table_size,
    const uint table_mask
) {
    uint gid = get_global_id(0);
    if (gid >= num_reads) return;

    int start = offsets[gid];
    int end = offsets[gid + 1];
    int len = end - start;

    if (len < (int)k) return;

    for (int i = 0; i <= len - (int)k; i++) {
        ulong encoded = encode_kmer(sequences, start + i, k);
        if (encoded == 0xFFFFFFFFFFFFFFFFUL) continue;

        ulong canon = canonical_kmer(encoded, k);
        ulong hash = murmur_finalize(canon);
        uint slot = (uint)(hash & table_mask);

        atomic_inc(&counts[slot]);
    }
}

// Kernel for batch k-mer extraction (returns encoded k-mers for CPU post-processing)
__kernel void extract_kmers(
    __global const uchar* sequences,
    __global const int* offsets,
    __global ulong* output_kmers,       // Output: encoded canonical k-mers
    __global int* output_counts,        // Output: count per k-mer position (always 1, for reduction)
    const uint k,
    const uint num_reads,
    const uint max_kmers_per_read       // Maximum k-mers to extract per read
) {
    uint gid = get_global_id(0);
    if (gid >= num_reads) return;

    int start = offsets[gid];
    int end = offsets[gid + 1];
    int len = end - start;

    if (len < (int)k) return;

    int num_kmers = len - k + 1;
    if (num_kmers > (int)max_kmers_per_read) {
        num_kmers = max_kmers_per_read;
    }

    int output_base = gid * max_kmers_per_read;

    for (int i = 0; i < num_kmers; i++) {
        ulong encoded = encode_kmer(sequences, start + i, k);
        ulong canon = 0xFFFFFFFFFFFFFFFFUL;

        if (encoded != 0xFFFFFFFFFFFFFFFFUL) {
            canon = canonical_kmer(encoded, k);
        }

        output_kmers[output_base + i] = canon;
        output_counts[output_base + i] = (canon != 0xFFFFFFFFFFFFFFFFUL) ? 1 : 0;
    }

    // Fill remaining slots with invalid markers
    for (int i = num_kmers; i < (int)max_kmers_per_read; i++) {
        output_kmers[output_base + i] = 0xFFFFFFFFFFFFFFFFUL;
        output_counts[output_base + i] = 0;
    }
}
