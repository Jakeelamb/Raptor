use ahash::RandomState;
use std::hash::{Hash, Hasher, BuildHasher};

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
        
        // Initialize hashers with different seeds
        let mut hashers = Vec::with_capacity(depth);
        for _ in 0..depth {
            hashers.push(RandomState::new());
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
            if self.matrix[i][idx] < u16::MAX {
                self.matrix[i][idx] += 1;
            }
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
    
    /// Calculate a hash index for the given item and row
    fn hash_index<T: Hash>(&self, item: &T, row: usize) -> usize {
        let mut hasher = self.hashers[row].build_hasher();
        item.hash(&mut hasher);
        let hash = hasher.finish() as usize;
        hash % self.width
    }
}
