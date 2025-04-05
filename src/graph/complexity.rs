use crate::graph::stitch::Path;
use std::collections::HashMap;

pub struct PathStats {
    pub total_paths: usize,
    pub average_length: f64,
    pub branch_count: usize,
}

pub fn compute_path_stats(paths: &[Path]) -> PathStats {
    let total = paths.len();
    let total_len: usize = paths.iter().map(|p| p.segments.len()).sum();
    let avg_len = if total > 0 { total_len as f64 / total as f64 } else { 0.0 };

    // crude: branch = segment reused across ≥2 paths
    let mut seen = HashMap::new();
    for p in paths {
        for &seg in &p.segments {
            *seen.entry(seg).or_insert(0) += 1;
        }
    }

    let branches = seen.values().filter(|&&v| v > 1).count();

    PathStats {
        total_paths: total,
        average_length: avg_len,
        branch_count: branches,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_compute_path_stats() {
        // Create test paths
        let path1 = Path {
            id: 0,
            segments: vec![0, 1, 2],
            overlaps: vec![5, 5],
        };
        
        let path2 = Path {
            id: 1,
            segments: vec![0, 3, 4],
            overlaps: vec![5, 5],
        };
        
        let path3 = Path {
            id: 2,
            segments: vec![5, 6, 7],
            overlaps: vec![5, 5],
        };
        
        let paths = vec![path1, path2, path3];
        
        let stats = compute_path_stats(&paths);
        
        assert_eq!(stats.total_paths, 3);
        assert_eq!(stats.average_length, 3.0);
        // Only segment 0 is shared between paths
        assert_eq!(stats.branch_count, 1);
    }
} 