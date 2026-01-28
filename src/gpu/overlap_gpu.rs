#[cfg(feature = "gpu")]
use ocl::{Buffer, ProQue, flags};

/// GPU-accelerated overlap detection between contigs
#[cfg(feature = "gpu")]
pub struct GpuOverlapFinder {
    pro_que: ProQue,
    max_results: usize,
}

#[cfg(feature = "gpu")]
impl GpuOverlapFinder {
    /// Create a new GPU overlap finder
    ///
    /// # Arguments
    /// * `max_contigs` - Maximum number of contigs to process
    /// * `max_results` - Maximum number of overlaps to return
    pub fn new(max_contigs: usize, max_results: usize) -> Result<Self, String> {
        let kernel_src = include_str!("kernels/overlap.cl");

        // Use 2D work size for the 2D kernel
        let pro_que = ProQue::builder()
            .src(kernel_src)
            .dims([max_contigs, max_contigs])
            .build()
            .map_err(|e| format!("Failed to build OpenCL program: {}", e))?;

        Ok(GpuOverlapFinder { pro_que, max_results })
    }

    /// Find overlaps between contigs
    ///
    /// # Arguments
    /// * `contigs` - Contig sequences
    /// * `min_overlap` - Minimum overlap length
    /// * `max_mismatches` - Maximum allowed mismatches
    ///
    /// # Returns
    /// Vector of (from_idx, to_idx, overlap_len) tuples
    pub fn find_overlaps(
        &self,
        contigs: &[String],
        min_overlap: usize,
        max_mismatches: usize,
    ) -> Result<Vec<(usize, usize, usize)>, String> {
        if contigs.is_empty() || contigs.len() < 2 {
            return Ok(Vec::new());
        }

        let num_contigs = contigs.len();

        // Flatten sequences
        let sequence_data: Vec<u8> = contigs.iter().flat_map(|s| s.bytes()).collect();

        // Build offset and length arrays
        let mut offsets = Vec::with_capacity(num_contigs);
        let mut lengths = Vec::with_capacity(num_contigs);
        let mut pos = 0i32;

        for contig in contigs {
            offsets.push(pos);
            lengths.push(contig.len() as i32);
            pos += contig.len() as i32;
        }

        // Create GPU buffers
        let sequences_buf = Buffer::<u8>::builder()
            .queue(self.pro_que.queue().clone())
            .flags(flags::MEM_READ_ONLY)
            .len(sequence_data.len())
            .copy_host_slice(&sequence_data)
            .build()
            .map_err(|e| format!("Failed to create sequence buffer: {}", e))?;

        let offsets_buf = Buffer::<i32>::builder()
            .queue(self.pro_que.queue().clone())
            .flags(flags::MEM_READ_ONLY)
            .len(offsets.len())
            .copy_host_slice(&offsets)
            .build()
            .map_err(|e| format!("Failed to create offsets buffer: {}", e))?;

        let lengths_buf = Buffer::<i32>::builder()
            .queue(self.pro_que.queue().clone())
            .flags(flags::MEM_READ_ONLY)
            .len(lengths.len())
            .copy_host_slice(&lengths)
            .build()
            .map_err(|e| format!("Failed to create lengths buffer: {}", e))?;

        // Output buffers
        let from_buf = Buffer::<i32>::builder()
            .queue(self.pro_que.queue().clone())
            .flags(flags::MEM_WRITE_ONLY)
            .len(self.max_results)
            .build()
            .map_err(|e| format!("Failed to create from buffer: {}", e))?;

        let to_buf = Buffer::<i32>::builder()
            .queue(self.pro_que.queue().clone())
            .flags(flags::MEM_WRITE_ONLY)
            .len(self.max_results)
            .build()
            .map_err(|e| format!("Failed to create to buffer: {}", e))?;

        let len_buf = Buffer::<i32>::builder()
            .queue(self.pro_que.queue().clone())
            .flags(flags::MEM_WRITE_ONLY)
            .len(self.max_results)
            .build()
            .map_err(|e| format!("Failed to create len buffer: {}", e))?;

        // Atomic counter for results
        let count_buf = Buffer::<i32>::builder()
            .queue(self.pro_que.queue().clone())
            .flags(flags::MEM_READ_WRITE)
            .len(1)
            .fill_val(0i32)
            .build()
            .map_err(|e| format!("Failed to create count buffer: {}", e))?;

        // Choose kernel based on number of contigs
        let kernel = if num_contigs <= 256 {
            // Use 2D kernel for smaller datasets
            self.pro_que
                .kernel_builder("find_overlaps_2d")
                .global_work_size([num_contigs, num_contigs])
                .arg(&sequences_buf)
                .arg(&offsets_buf)
                .arg(&lengths_buf)
                .arg(&from_buf)
                .arg(&to_buf)
                .arg(&len_buf)
                .arg(&count_buf)
                .arg(num_contigs as u32)
                .arg(min_overlap as u32)
                .arg(max_mismatches as u32)
                .arg(self.max_results as u32)
                .build()
                .map_err(|e| format!("Failed to build 2D kernel: {}", e))?
        } else {
            // Use 1D linearized kernel for larger datasets
            let total_pairs = num_contigs * (num_contigs - 1);
            self.pro_que
                .kernel_builder("find_overlaps_1d")
                .global_work_size(total_pairs)
                .arg(&sequences_buf)
                .arg(&offsets_buf)
                .arg(&lengths_buf)
                .arg(&from_buf)
                .arg(&to_buf)
                .arg(&len_buf)
                .arg(&count_buf)
                .arg(num_contigs as u32)
                .arg(min_overlap as u32)
                .arg(max_mismatches as u32)
                .arg(self.max_results as u32)
                .build()
                .map_err(|e| format!("Failed to build 1D kernel: {}", e))?
        };

        // Execute kernel
        unsafe {
            kernel.enq().map_err(|e| format!("Failed to execute kernel: {}", e))?;
        }

        // Read back result count
        let mut count = vec![0i32; 1];
        count_buf.read(&mut count)
            .enq()
            .map_err(|e| format!("Failed to read count: {}", e))?;

        let result_count = (count[0] as usize).min(self.max_results);

        if result_count == 0 {
            return Ok(Vec::new());
        }

        // Read back results
        let mut from_indices = vec![0i32; result_count];
        let mut to_indices = vec![0i32; result_count];
        let mut overlap_lens = vec![0i32; result_count];

        from_buf.read(&mut from_indices)
            .enq()
            .map_err(|e| format!("Failed to read from indices: {}", e))?;

        to_buf.read(&mut to_indices)
            .enq()
            .map_err(|e| format!("Failed to read to indices: {}", e))?;

        len_buf.read(&mut overlap_lens)
            .enq()
            .map_err(|e| format!("Failed to read overlap lengths: {}", e))?;

        // Convert to result format
        let results: Vec<(usize, usize, usize)> = from_indices
            .iter()
            .zip(to_indices.iter())
            .zip(overlap_lens.iter())
            .map(|((&from, &to), &len)| (from as usize, to as usize, len as usize))
            .collect();

        Ok(results)
    }
}

/// CPU fallback for overlap finding
pub fn find_overlaps_cpu(
    contigs: &[String],
    min_overlap: usize,
    max_mismatches: usize,
) -> Vec<(usize, usize, usize)> {
    use rayon::prelude::*;

    if contigs.len() < 2 {
        return Vec::new();
    }

    // Generate all pairs
    let pairs: Vec<(usize, usize)> = (0..contigs.len())
        .flat_map(|i| (0..contigs.len()).map(move |j| (i, j)))
        .filter(|&(i, j)| i != j)
        .collect();

    // Parallel overlap detection
    pairs
        .par_iter()
        .filter_map(|&(i, j)| {
            find_suffix_prefix_overlap(&contigs[i], &contigs[j], min_overlap, max_mismatches)
                .map(|overlap_len| (i, j, overlap_len))
        })
        .collect()
}

/// Find the longest suffix-prefix overlap between two sequences
fn find_suffix_prefix_overlap(
    from: &str,
    to: &str,
    min_overlap: usize,
    max_mismatches: usize,
) -> Option<usize> {
    let max_possible = from.len().min(to.len());
    if max_possible < min_overlap {
        return None;
    }

    let from_bytes = from.as_bytes();
    let to_bytes = to.as_bytes();

    // Try from longest to shortest overlap
    for overlap in (min_overlap..=max_possible).rev() {
        let suffix_start = from.len() - overlap;
        let mut mismatches = 0;

        for i in 0..overlap {
            if from_bytes[suffix_start + i] != to_bytes[i] {
                mismatches += 1;
                if mismatches > max_mismatches {
                    break;
                }
            }
        }

        if mismatches <= max_mismatches {
            return Some(overlap);
        }
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_suffix_prefix_overlap() {
        // Perfect overlap
        let overlap = find_suffix_prefix_overlap("ACGTACGT", "ACGTAAAA", 4, 0);
        assert_eq!(overlap, Some(4));

        // Overlap with mismatch
        let overlap = find_suffix_prefix_overlap("ACGTACGT", "ACGTTAAA", 4, 1);
        assert_eq!(overlap, Some(4));

        // No overlap
        let overlap = find_suffix_prefix_overlap("AAAA", "TTTT", 2, 0);
        assert!(overlap.is_none());
    }

    #[test]
    fn test_cpu_overlap_finding() {
        let contigs = vec![
            "ACGTACGTACGT".to_string(),
            "ACGTAAAAAAA".to_string(),
            "TTTTTTTACGT".to_string(),
        ];

        let overlaps = find_overlaps_cpu(&contigs, 4, 0);

        // Should find overlap from contig 0 to contig 1 (ACGT)
        assert!(overlaps.iter().any(|&(from, to, _)| from == 0 && to == 1));
    }
}
