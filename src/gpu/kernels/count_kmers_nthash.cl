// GPU k-mer counting kernel using ntHash for O(1) rolling updates.
//
// ntHash provides constant-time hash updates when sliding the k-mer window,
// making it 3-4x faster than recomputing hashes from scratch.
//
// Reference: Mohamadi et al. (2016) "ntHash: recursive nucleotide hashing"

// ntHash constants from the paper
#define NT_A 0x3c8bfbb395c60474UL
#define NT_C 0x3193c18562a02b4cUL
#define NT_G 0x20323ed082572324UL
#define NT_T 0x295549f54be24456UL
#define NT_INVALID 0xFFFFFFFFFFFFFFFFUL

// Lookup table for base -> ntHash value (in constant memory for fast access)
__constant ulong NT_HASH[4] = { NT_A, NT_C, NT_G, NT_T };
__constant ulong NT_HASH_RC[4] = { NT_T, NT_G, NT_C, NT_A };  // Complement values

// Encode base to 0-3, return 255 for invalid
inline uchar encode_base(uchar base) {
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 255;
    }
}

// Get ntHash value for a base
inline ulong get_nt_hash(uchar base) {
    uchar enc = encode_base(base);
    return (enc < 4) ? NT_HASH[enc] : NT_INVALID;
}

// Get reverse complement ntHash value
inline ulong get_nt_hash_rc(uchar base) {
    uchar enc = encode_base(base);
    return (enc < 4) ? NT_HASH_RC[enc] : NT_INVALID;
}

// Rotate left 64-bit
inline ulong rotate_left_64(ulong x, uint n) {
    return (x << n) | (x >> (64 - n));
}

// Rotate right 64-bit
inline ulong rotate_right_64(ulong x, uint n) {
    return (x >> n) | (x << (64 - n));
}

// Initialize ntHash for first k-mer
// Returns NT_INVALID if any base is invalid
inline ulong nthash_init(
    __global const uchar* seq,
    int start,
    int k,
    ulong* forward_out,
    ulong* reverse_out
) {
    ulong forward = 0;
    ulong reverse = 0;

    for (int i = 0; i < k; i++) {
        ulong h = get_nt_hash(seq[start + i]);
        ulong h_rc = get_nt_hash_rc(seq[start + i]);

        if (h == NT_INVALID) {
            return NT_INVALID;
        }

        // Forward: rotate left and XOR
        forward = rotate_left_64(forward, 1) ^ h;

        // Reverse: build from end using complement
        reverse ^= rotate_left_64(h_rc, k - 1 - i);
    }

    *forward_out = forward;
    *reverse_out = reverse;

    // Return canonical (minimum of forward and reverse)
    return (forward < reverse) ? forward : reverse;
}

// Roll ntHash to next position in O(1)
inline ulong nthash_roll(
    uchar out_base,
    uchar in_base,
    int k,
    ulong* forward,
    ulong* reverse
) {
    ulong h_out = get_nt_hash(out_base);
    ulong h_in = get_nt_hash(in_base);
    ulong h_out_rc = get_nt_hash_rc(out_base);
    ulong h_in_rc = get_nt_hash_rc(in_base);

    if (h_in == NT_INVALID) {
        return NT_INVALID;
    }

    // Update forward hash
    *forward = rotate_left_64(*forward, 1)
             ^ rotate_left_64(h_out, k)
             ^ h_in;

    // Update reverse complement hash
    *reverse = rotate_right_64(*reverse, 1)
             ^ h_out_rc
             ^ rotate_left_64(h_in_rc, k - 1);

    return (*forward < *reverse) ? *forward : *reverse;
}

// Main kernel: Count k-mers using ntHash rolling updates
// Each work item processes one read with O(1) per k-mer
__kernel void count_kmers_nthash(
    __global const uchar* sequences,    // Flattened sequence data
    __global const int* offsets,        // Start offset for each read
    __global uint* hash_table_counts,   // Hash table for counts
    __global ulong* hash_table_kmers,   // Hash table for k-mer hashes
    const uint k,                       // K-mer size
    const uint num_reads,               // Number of reads
    const uint table_size,              // Hash table size
    const uint table_mask               // table_size - 1
) {
    uint gid = get_global_id(0);
    if (gid >= num_reads) return;

    int start = offsets[gid];
    int end = offsets[gid + 1];
    int len = end - start;

    if (len < (int)k) return;

    ulong forward, reverse;

    // Initialize first k-mer
    ulong hash = nthash_init(sequences, start, k, &forward, &reverse);

    if (hash != NT_INVALID) {
        // Insert first k-mer
        uint slot = (uint)(hash & table_mask);
        for (uint probe = 0; probe < 32; probe++) {
            uint idx = (slot + probe) & table_mask;
            ulong expected = 0;
            ulong old = atom_cmpxchg(&hash_table_kmers[idx], expected, hash);
            if (old == 0 || old == hash) {
                atomic_inc(&hash_table_counts[idx]);
                break;
            }
        }
    }

    // Roll through remaining k-mers with O(1) updates
    for (int i = 1; i <= len - (int)k; i++) {
        uchar out_base = sequences[start + i - 1];
        uchar in_base = sequences[start + i + k - 1];

        hash = nthash_roll(out_base, in_base, k, &forward, &reverse);

        if (hash == NT_INVALID) {
            // Invalid base - need to reinitialize
            hash = nthash_init(sequences, start + i, k, &forward, &reverse);
            if (hash == NT_INVALID) continue;
        }

        // Insert k-mer with linear probing
        uint slot = (uint)(hash & table_mask);
        for (uint probe = 0; probe < 32; probe++) {
            uint idx = (slot + probe) & table_mask;
            ulong expected = 0;
            ulong old = atom_cmpxchg(&hash_table_kmers[idx], expected, hash);
            if (old == 0 || old == hash) {
                atomic_inc(&hash_table_counts[idx]);
                break;
            }
        }
    }
}

// Simplified kernel for sketch-based counting (no collision resolution)
// Faster but approximate - good for initial filtering
__kernel void count_kmers_nthash_sketch(
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

    ulong forward, reverse;
    ulong hash = nthash_init(sequences, start, k, &forward, &reverse);

    if (hash != NT_INVALID) {
        uint slot = (uint)(hash & table_mask);
        atomic_inc(&counts[slot]);
    }

    for (int i = 1; i <= len - (int)k; i++) {
        uchar out_base = sequences[start + i - 1];
        uchar in_base = sequences[start + i + k - 1];

        hash = nthash_roll(out_base, in_base, k, &forward, &reverse);

        if (hash == NT_INVALID) {
            hash = nthash_init(sequences, start + i, k, &forward, &reverse);
            if (hash == NT_INVALID) continue;
        }

        uint slot = (uint)(hash & table_mask);
        atomic_inc(&counts[slot]);
    }
}

// Two-pass counting with Bloom filter pre-filtering
// Pass 1: Populate 4-bit counting Bloom filter
// Pass 2: Only count k-mers that appear 2+ times
__kernel void bloom_filter_pass1(
    __global const uchar* sequences,
    __global const int* offsets,
    __global ulong* bloom_counters,     // 4-bit counters packed in u64
    const uint k,
    const uint num_reads,
    const uint num_counters,
    const uint counter_mask
) {
    uint gid = get_global_id(0);
    if (gid >= num_reads) return;

    int start = offsets[gid];
    int end = offsets[gid + 1];
    int len = end - start;

    if (len < (int)k) return;

    ulong forward, reverse;
    ulong hash = nthash_init(sequences, start, k, &forward, &reverse);

    // Insert into Bloom filter (increment 4-bit counter)
    if (hash != NT_INVALID) {
        uint counter_idx = (uint)(hash & counter_mask);
        uint word_idx = counter_idx / 16;
        uint nibble_offset = (counter_idx % 16) * 4;

        // Atomic increment of 4-bit counter (saturate at 15)
        ulong old_val, new_val;
        do {
            old_val = bloom_counters[word_idx];
            ulong current = (old_val >> nibble_offset) & 0xF;
            if (current >= 15) break;  // Already saturated
            new_val = (old_val & ~(0xFUL << nibble_offset))
                    | ((current + 1) << nibble_offset);
        } while (atom_cmpxchg(&bloom_counters[word_idx], old_val, new_val) != old_val);
    }

    for (int i = 1; i <= len - (int)k; i++) {
        uchar out_base = sequences[start + i - 1];
        uchar in_base = sequences[start + i + k - 1];

        hash = nthash_roll(out_base, in_base, k, &forward, &reverse);

        if (hash == NT_INVALID) {
            hash = nthash_init(sequences, start + i, k, &forward, &reverse);
            if (hash == NT_INVALID) continue;
        }

        uint counter_idx = (uint)(hash & counter_mask);
        uint word_idx = counter_idx / 16;
        uint nibble_offset = (counter_idx % 16) * 4;

        ulong old_val, new_val;
        do {
            old_val = bloom_counters[word_idx];
            ulong current = (old_val >> nibble_offset) & 0xF;
            if (current >= 15) break;
            new_val = (old_val & ~(0xFUL << nibble_offset))
                    | ((current + 1) << nibble_offset);
        } while (atom_cmpxchg(&bloom_counters[word_idx], old_val, new_val) != old_val);
    }
}

// Pass 2: Count only k-mers with Bloom count >= threshold
__kernel void bloom_filter_pass2(
    __global const uchar* sequences,
    __global const int* offsets,
    __global const ulong* bloom_counters,  // Read-only from pass 1
    __global uint* hash_table_counts,
    __global ulong* hash_table_kmers,
    const uint k,
    const uint num_reads,
    const uint num_counters,
    const uint counter_mask,
    const uint table_size,
    const uint table_mask,
    const uint min_bloom_count            // Threshold (typically 2)
) {
    uint gid = get_global_id(0);
    if (gid >= num_reads) return;

    int start = offsets[gid];
    int end = offsets[gid + 1];
    int len = end - start;

    if (len < (int)k) return;

    ulong forward, reverse;
    ulong hash = nthash_init(sequences, start, k, &forward, &reverse);

    if (hash != NT_INVALID) {
        // Check Bloom filter
        uint counter_idx = (uint)(hash & counter_mask);
        uint word_idx = counter_idx / 16;
        uint nibble_offset = (counter_idx % 16) * 4;
        ulong count = (bloom_counters[word_idx] >> nibble_offset) & 0xF;

        if (count >= min_bloom_count) {
            // This k-mer appears multiple times - count it
            uint slot = (uint)(hash & table_mask);
            for (uint probe = 0; probe < 32; probe++) {
                uint idx = (slot + probe) & table_mask;
                ulong expected = 0;
                ulong old = atom_cmpxchg(&hash_table_kmers[idx], expected, hash);
                if (old == 0 || old == hash) {
                    atomic_inc(&hash_table_counts[idx]);
                    break;
                }
            }
        }
    }

    for (int i = 1; i <= len - (int)k; i++) {
        uchar out_base = sequences[start + i - 1];
        uchar in_base = sequences[start + i + k - 1];

        hash = nthash_roll(out_base, in_base, k, &forward, &reverse);

        if (hash == NT_INVALID) {
            hash = nthash_init(sequences, start + i, k, &forward, &reverse);
            if (hash == NT_INVALID) continue;
        }

        uint counter_idx = (uint)(hash & counter_mask);
        uint word_idx = counter_idx / 16;
        uint nibble_offset = (counter_idx % 16) * 4;
        ulong count = (bloom_counters[word_idx] >> nibble_offset) & 0xF;

        if (count >= min_bloom_count) {
            uint slot = (uint)(hash & table_mask);
            for (uint probe = 0; probe < 32; probe++) {
                uint idx = (slot + probe) & table_mask;
                ulong expected = 0;
                ulong old = atom_cmpxchg(&hash_table_kmers[idx], expected, hash);
                if (old == 0 || old == hash) {
                    atomic_inc(&hash_table_counts[idx]);
                    break;
                }
            }
        }
    }
}
