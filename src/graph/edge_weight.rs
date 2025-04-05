use std::collections::HashMap;

/// Compute edge weight considering both coverage and paired-end support
///
/// # Arguments
///
/// * `from` - Source node ID
/// * `to` - Target node ID
/// * `coverage` - Vector of coverage values per node
/// * `pair_support` - Number of read pairs linking these two contigs
///
/// # Returns
///
/// Enhanced edge weight as f32
pub fn compute_edge_weight(from: usize, to: usize, coverage: &[f32], pair_support: usize) -> f32 {
    // Base weight is the average coverage of the two nodes
    let base_weight = if from < coverage.len() && to < coverage.len() {
        (coverage[from] + coverage[to]) / 2.0
    } else {
        0.0
    };
    
    // Add bonus for paired-end support (5.0 per pair)
    // This makes read pair evidence strongly influence the path selection
    let pair_bonus = pair_support as f32 * 5.0;
    
    base_weight + pair_bonus
}

/// Track and record paired-end links between contigs
pub struct PairLinkTracker {
    /// Maps (from_node, to_node) to count of supporting pairs
    pair_links: HashMap<(usize, usize), usize>,
}

impl PairLinkTracker {
    /// Create a new PairLinkTracker
    pub fn new() -> Self {
        Self {
            pair_links: HashMap::new(),
        }
    }
    
    /// Add a pair link between two contigs
    pub fn add_pair_link(&mut self, from: usize, to: usize) {
        *self.pair_links.entry((from, to)).or_insert(0) += 1;
    }
    
    /// Get the number of pairs supporting an edge
    pub fn get_pair_support(&self, from: usize, to: usize) -> usize {
        *self.pair_links.get(&(from, to)).unwrap_or(&0)
    }
    
    /// Get all pair links
    pub fn get_all_pair_links(&self) -> &HashMap<(usize, usize), usize> {
        &self.pair_links
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_compute_edge_weight() {
        let coverage = vec![10.0, 20.0, 30.0];
        
        // Test with no pair support
        let weight1 = compute_edge_weight(0, 1, &coverage, 0);
        assert_eq!(weight1, 15.0); // (10 + 20) / 2
        
        // Test with pair support
        let weight2 = compute_edge_weight(1, 2, &coverage, 2);
        assert_eq!(weight2, 35.0); // (20 + 30) / 2 + 2 * 5.0
    }
    
    #[test]
    fn test_pair_link_tracker() {
        let mut tracker = PairLinkTracker::new();
        
        // Add some links
        tracker.add_pair_link(0, 1);
        tracker.add_pair_link(1, 2);
        tracker.add_pair_link(0, 1); // Add another link for the same pair
        
        // Check pair support
        assert_eq!(tracker.get_pair_support(0, 1), 2);
        assert_eq!(tracker.get_pair_support(1, 2), 1);
        assert_eq!(tracker.get_pair_support(2, 3), 0); // Non-existent pair
    }
} 