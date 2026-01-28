use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use std::io::{BufWriter, Write};
use fasthash::xx;

/// Partition contigs/sequences by minimizer hash bucket
pub fn partition_by_minimizer(sequences: &[String], k: usize, buckets: usize) -> HashMap<usize, Vec<String>> {
    let mut map: HashMap<usize, Vec<String>> = HashMap::new();

    for seq in sequences {
        if seq.len() < k { continue; }
        
        // Find the lexicographically minimum k-mer
        let min_kmer = seq
            .as_bytes()
            .windows(k)
            .map(|w| w.iter().map(|&b| b as char).collect::<String>())
            .min()
            .unwrap();

        // Use xxHash for fast hashing
        let hash = xx::hash64(min_kmer.as_bytes()) as usize % buckets;
        map.entry(hash).or_default().push(seq.clone());
    }

    map
}

/// Save partitioned sequences to separate files
pub fn save_partitions(partitions: &HashMap<usize, Vec<String>>, output_dir: &str, prefix: &str) -> Vec<String> {
    let mut file_paths = Vec::new();
    
    // Create output directory if it doesn't exist
    std::fs::create_dir_all(output_dir).expect("Failed to create output directory");
    
    for (bucket_id, sequences) in partitions {
        let file_name = format!("{}_{:02}.fasta", prefix, bucket_id);
        let file_path = Path::new(output_dir).join(&file_name);
        
        let file = File::create(&file_path).expect("Failed to create partition file");
        let mut writer = BufWriter::new(file);
        
        for (i, seq) in sequences.iter().enumerate() {
            writeln!(writer, ">seq_{}_{}", bucket_id, i).expect("Failed to write header");
            writeln!(writer, "{}", seq).expect("Failed to write sequence");
        }
        
        file_paths.push(file_path.to_string_lossy().to_string());
    }
    
    file_paths
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_partition_by_minimizer() {
        let sequences = vec![
            "ATCGATCGATCG".to_string(),
            "GATCGATCGATC".to_string(),
            "CGATCGATCGAT".to_string(),
            "TCGATCGATCGA".to_string(),
        ];
        
        let partitions = partition_by_minimizer(&sequences, 5, 2);
        
        // Should distribute sequences into buckets
        assert!(partitions.len() <= 2);
        assert!(partitions.values().all(|v| !v.is_empty()));
        
        // Sum of sequences in all buckets should equal original count
        let total_seqs: usize = partitions.values().map(|v| v.len()).sum();
        assert_eq!(total_seqs, sequences.len());
    }
}
