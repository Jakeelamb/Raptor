__kernel void count_kmers(
    __global const uchar* sequences,
    __global const int* offsets,
    __global uint* counts,
    const uint k,
    const uint max_reads
) {
    uint gid = get_global_id(0);
    if (gid >= max_reads) return;

    int offset = offsets[gid];
    int next_offset = offsets[gid + 1];
    int len = next_offset - offset;

    for (int i = 0; i <= len - k; i++) {
        uint hash = 0;
        for (int j = 0; j < k; j++) {
            uchar base = sequences[offset + i + j];
            hash = hash * 4 + (base == 'A' ? 0 : base == 'C' ? 1 : base == 'G' ? 2 : 3);
        }
        atomic_inc(&counts[hash % 65536]); // naive sketch
    }
}
