use fasthash::xx;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Partition contigs/sequences by minimizer hash bucket
pub fn partition_by_minimizer(
    sequences: &[String],
    k: usize,
    buckets: usize,
) -> HashMap<usize, Vec<String>> {
    if k == 0 || buckets == 0 {
        return HashMap::new();
    }

    let mut map: HashMap<usize, Vec<String>> = HashMap::new();

    for seq in sequences {
        if seq.len() < k {
            continue;
        }

        // Find the lexicographically minimum k-mer without allocating per window.
        let mut windows = seq.as_bytes().windows(k);
        let mut min_kmer = windows
            .next()
            .expect("sequence length is checked to be >= k above");
        for window in windows {
            if window < min_kmer {
                min_kmer = window;
            }
        }

        // Use xxHash for fast hashing
        let hash = xx::hash64(min_kmer) as usize % buckets;
        map.entry(hash).or_default().push(seq.clone());
    }

    map
}

/// Save partitioned sequences to separate files
pub fn save_partitions(
    partitions: &HashMap<usize, Vec<String>>,
    output_dir: &str,
    prefix: &str,
) -> Vec<String> {
    let mut file_paths = Vec::new();

    // Create output directory if it doesn't exist
    std::fs::create_dir_all(output_dir).expect("Failed to create output directory");

    let mut bucket_ids: Vec<usize> = partitions.keys().copied().collect();
    bucket_ids.sort_unstable();

    for bucket_id in bucket_ids {
        let sequences = &partitions[&bucket_id];
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
    use std::path::Path;
    use tempfile::TempDir;

    fn partition_by_minimizer_reference(
        sequences: &[String],
        k: usize,
        buckets: usize,
    ) -> HashMap<usize, Vec<String>> {
        if k == 0 || buckets == 0 {
            return HashMap::new();
        }

        let mut map: HashMap<usize, Vec<String>> = HashMap::new();
        for seq in sequences {
            if seq.len() < k {
                continue;
            }

            let min_kmer = seq
                .as_bytes()
                .windows(k)
                .map(|w| w.iter().map(|&b| b as char).collect::<String>())
                .min()
                .unwrap();
            let hash = xx::hash64(min_kmer.as_bytes()) as usize % buckets;
            map.entry(hash).or_default().push(seq.clone());
        }

        map
    }

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

    #[test]
    fn partition_by_minimizer_handles_invalid_parameters() {
        let sequences = vec!["ACGT".to_string()];
        assert!(partition_by_minimizer(&sequences, 0, 4).is_empty());
        assert!(partition_by_minimizer(&sequences, 3, 0).is_empty());
    }

    #[test]
    fn partition_by_minimizer_matches_reference_implementation() {
        let sequences = vec![
            "ATCGATCGATCG".to_string(),
            "TTTTACGTTTTT".to_string(),
            "GCGCGCGCGCGA".to_string(),
            "AACCAACCAACC".to_string(),
            "CCCCAAAAGGGG".to_string(),
        ];

        for &k in &[1usize, 3, 5, 7] {
            let observed = partition_by_minimizer(&sequences, k, 11);
            let expected = partition_by_minimizer_reference(&sequences, k, 11);
            assert_eq!(observed, expected);
        }
    }

    #[test]
    fn save_partitions_outputs_sorted_bucket_paths() {
        let mut partitions = HashMap::new();
        partitions.insert(12usize, vec!["AAAA".to_string()]);
        partitions.insert(2usize, vec!["CCCC".to_string()]);
        partitions.insert(7usize, vec!["GGGG".to_string()]);

        let tmp = TempDir::new().unwrap();
        let paths = save_partitions(&partitions, tmp.path().to_str().unwrap(), "partition");

        let names: Vec<String> = paths
            .iter()
            .map(|p| {
                Path::new(p)
                    .file_name()
                    .unwrap()
                    .to_string_lossy()
                    .to_string()
            })
            .collect();
        assert_eq!(
            names,
            vec![
                "partition_02.fasta".to_string(),
                "partition_07.fasta".to_string(),
                "partition_12.fasta".to_string(),
            ]
        );
    }
}
