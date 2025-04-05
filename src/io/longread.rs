use crate::io::fastq::FastqRecord;
use crate::accel::simd::match_kmers_simd;

/// Filter nanopore reads by length
pub fn filter_nanopore_reads(reads: &[FastqRecord], min_len: usize) -> Vec<FastqRecord> {
    reads.iter()
        .filter(|r| r.sequence.len() >= min_len)
        .cloned()
        .collect()
}

/// Align long reads to a contig using basic exact matching
/// Returns a vector of (start, end) positions
pub fn align_long_reads(contig: &str, long_reads: &[FastqRecord]) -> Vec<(usize, usize)> {
    let contig_bytes = contig.as_bytes();
    
    long_reads.iter()
        .filter_map(|r| {
            let read_bytes = r.sequence.as_bytes();
            
            // For large contigs, we use a sliding window approach
            if contig.len() > 1000 {
                for window_start in (0..contig.len()).step_by(500) {
                    let window_end = (window_start + 1000).min(contig.len());
                    let window = &contig_bytes[window_start..window_end];
                    
                    // Try to find an approximate match using SIMD acceleration
                    if match_kmers_simd(window, read_bytes, 5) {
                        // Found an approximate match, now find the exact position
                        if let Some(pos) = find_exact_position(contig, &r.sequence, window_start) {
                            return Some((pos, pos + r.sequence.len()));
                        }
                    }
                }
                None
            } else {
                // For small contigs, just do a direct search
                contig.find(&r.sequence)
                    .map(|pos| (pos, pos + r.sequence.len()))
            }
        })
        .collect()
}

/// Find the exact position of a read in a contig given an approximate starting position
fn find_exact_position(contig: &str, read: &str, approximate_start: usize) -> Option<usize> {
    let search_start = approximate_start.saturating_sub(100);
    let search_end = (approximate_start + 100).min(contig.len());
    let search_region = &contig[search_start..search_end];
    
    search_region.find(read).map(|pos| search_start + pos)
}

/// Integrate long read alignments into the polishing workflow
pub fn integrate_long_reads_for_polishing(
    contig: &str, 
    short_reads: &[String], 
    long_reads: &[FastqRecord],
    _window_size: usize
) -> String {
    use crate::graph::polish::polish_contig_string;
    
    // First, align long reads to get coverage data
    let long_read_alignments = align_long_reads(contig, long_reads);
    
    // Convert long read alignments to sequences for polishing
    let mut all_reads = short_reads.to_vec();
    
    for (start, end) in long_read_alignments {
        if start < contig.len() && end <= contig.len() {
            // Add the aligned portion of the long read
            all_reads.push(contig[start..end].to_string());
        }
    }
    
    // Use the string-based polishing function with both short and long reads
    polish_contig_string(contig, &all_reads, 0.7)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_filter_nanopore_reads() {
        let reads = vec![
            FastqRecord {
                header: "@read1".to_string(),
                sequence: "ATCGATCG".to_string(), // 8 bp
                plus: "+".to_string(),
                quality: "IIIIIIII".to_string(),
            },
            FastqRecord {
                header: "@read2".to_string(),
                sequence: "ATCG".to_string(), // 4 bp
                plus: "+".to_string(),
                quality: "IIII".to_string(),
            },
            FastqRecord {
                header: "@read3".to_string(),
                sequence: "ATCGATCGATCGATCG".to_string(), // 16 bp
                plus: "+".to_string(),
                quality: "IIIIIIIIIIIIIIII".to_string(),
            },
        ];
        
        let filtered = filter_nanopore_reads(&reads, 10);
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].header, "@read3");
    }
    
    #[test]
    fn test_align_long_reads() {
        let contig = "ATCGATCGATCGATCGATCGATCGATCG";
        let reads = vec![
            FastqRecord {
                header: "@read1".to_string(),
                sequence: "ATCGATCG".to_string(),
                plus: "+".to_string(),
                quality: "IIIIIIII".to_string(),
            },
            FastqRecord {
                header: "@read2".to_string(),
                sequence: "GATCGATCG".to_string(),
                plus: "+".to_string(),
                quality: "IIIIIIIII".to_string(),
            },
        ];
        
        let alignments = align_long_reads(contig, &reads);
        assert_eq!(alignments.len(), 2);
        assert_eq!(alignments[0], (0, 8)); // First read aligns at the beginning
        assert_eq!(alignments[1], (3, 12)); // Second read aligns after 3 bases
    }
} 