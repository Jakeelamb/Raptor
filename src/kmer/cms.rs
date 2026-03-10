use ahash::RandomState;
use std::hash::Hash;

/// A Count-Min Sketch implementation for k-mer counting
pub struct CountMinSketch {
    matrix: Vec<Vec<u16>>,
    width: usize,
    depth: usize,
    hashers: Vec<RandomState>,
}

impl CountMinSketch {
    /// Create a new Count-Min Sketch with the given depth and width
    pub fn new(depth: usize, width: usize) -> Self {
        // Initialize matrix with zeros
        let matrix = vec![vec![0; width]; depth];

        // Initialize hashers with deterministic row-specific seeds.
        // This makes sketch behavior reproducible across runs.
        let mut hashers = Vec::with_capacity(depth);
        for row in 0..depth {
            hashers.push(Self::hasher_for_row(row));
        }

        CountMinSketch {
            matrix,
            width,
            depth,
            hashers,
        }
    }

    /// Insert a k-mer into the sketch
    pub fn insert<T: Hash>(&mut self, item: &T) {
        for i in 0..self.depth {
            let idx = self.hash_index(item, i);
            self.matrix[i][idx] = self.matrix[i][idx].saturating_add(1);
        }
    }

    /// Get the estimated count for a k-mer
    pub fn estimate<T: Hash>(&self, item: &T) -> u16 {
        let mut min_count = u16::MAX;
        for i in 0..self.depth {
            let idx = self.hash_index(item, i);
            min_count = min_count.min(self.matrix[i][idx]);
        }
        min_count
    }

    /// Insert a pre-computed hash into the sketch.
    /// Uses double hashing to generate multiple indices from a single hash.
    #[inline]
    pub fn insert_hash(&mut self, hash: u64) {
        let h1 = hash as usize;
        let h2 = (hash >> 32) as usize;
        for i in 0..self.depth {
            let idx = (h1.wrapping_add(i.wrapping_mul(h2))) % self.width;
            self.matrix[i][idx] = self.matrix[i][idx].saturating_add(1);
        }
    }

    /// Get the estimated count for a pre-computed hash.
    /// Uses double hashing to generate multiple indices from a single hash.
    #[inline]
    pub fn estimate_hash(&self, hash: u64) -> u16 {
        let h1 = hash as usize;
        let h2 = (hash >> 32) as usize;
        let mut min_count = u16::MAX;
        for i in 0..self.depth {
            let idx = (h1.wrapping_add(i.wrapping_mul(h2))) % self.width;
            min_count = min_count.min(self.matrix[i][idx]);
        }
        min_count
    }

    /// Calculate a hash index for the given item and row
    fn hash_index<T: Hash>(&self, item: &T, row: usize) -> usize {
        let hash = self.hashers[row].hash_one(item) as usize;
        hash % self.width
    }

    #[inline]
    fn hasher_for_row(row: usize) -> RandomState {
        const S0: u64 = 0x243f_6a88_85a3_08d3;
        const S1: u64 = 0x1319_8a2e_0370_7344;
        const S2: u64 = 0xa409_3822_299f_31d0;
        const S3: u64 = 0x082e_fa98_ec4e_6c89;
        let row = row as u64;
        RandomState::with_seeds(
            S0 ^ row.wrapping_mul(0x9e37_79b9_7f4a_7c15),
            S1 ^ row.wrapping_mul(0xc2b2_ae3d_27d4_eb4f),
            S2 ^ row.wrapping_mul(0x1656_67b1_9e37_79f9),
            S3 ^ row.wrapping_mul(0x27d4_eb2f_1656_67c5),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::CountMinSketch;

    #[test]
    fn sketch_estimates_are_reproducible_across_instances() {
        let mut a = CountMinSketch::new(4, 1024);
        let mut b = CountMinSketch::new(4, 1024);

        for _ in 0..11 {
            a.insert(&"ACGT");
            b.insert(&"ACGT");
        }
        for _ in 0..7 {
            a.insert(&"TGCA");
            b.insert(&"TGCA");
        }

        assert_eq!(a.estimate(&"ACGT"), b.estimate(&"ACGT"));
        assert_eq!(a.estimate(&"TGCA"), b.estimate(&"TGCA"));
        assert_eq!(a.estimate(&"GGGG"), b.estimate(&"GGGG"));
    }
}
