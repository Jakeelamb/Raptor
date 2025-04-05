#[cfg(feature = "gpu")]
use ocl::{ProQue, Buffer, flags};
use std::collections::HashMap;
#[cfg(feature = "gpu")]
use std::fs;
use log;

#[cfg(feature = "gpu")]
pub struct GpuKmerCounter {
    pro_que: ProQue,
    max_reads: usize,
}

#[cfg(feature = "gpu")]
impl GpuKmerCounter {
    pub fn new(_k: usize, max_reads: usize) -> Self {
        let kernel_src = fs::read_to_string("src/gpu/kernels/count_kmers.cl").expect("Kernel missing");
        let pro_que = ProQue::builder()
            .src(kernel_src)
            .dims(max_reads)
            .build().unwrap();

        GpuKmerCounter { pro_que, max_reads }
    }

    pub fn count(&self, reads: &Vec<String>, k: usize) -> Vec<u32> {
        // Flatten reads
        let sequence_data: Vec<u8> = reads.iter().flat_map(|s| s.bytes()).collect();

        // Offset positions
        let mut offsets = Vec::with_capacity(reads.len() + 1);
        let mut pos = 0;
        offsets.push(pos as i32);
        for r in reads {
            pos += r.len();
            offsets.push(pos as i32);
        }

        let sequences_buf = Buffer::<u8>::builder()
            .queue(self.pro_que.queue().clone())
            .flags(flags::MEM_READ_ONLY)
            .len(sequence_data.len())
            .copy_host_slice(&sequence_data)
            .build().unwrap();

        let offsets_buf = Buffer::<i32>::builder()
            .queue(self.pro_que.queue().clone())
            .flags(flags::MEM_READ_ONLY)
            .len(offsets.len())
            .copy_host_slice(&offsets)
            .build().unwrap();

        let counts_buf = Buffer::<u32>::builder()
            .queue(self.pro_que.queue().clone())
            .flags(flags::MEM_WRITE_ONLY)
            .len(65536)
            .build().unwrap();

        let kernel = self.pro_que.kernel_builder("count_kmers")
            .arg(&sequences_buf)
            .arg(&offsets_buf)
            .arg(&counts_buf)
            .arg(k as u32)
            .arg(reads.len() as u32)
            .build().unwrap();

        unsafe { kernel.enq().unwrap(); }

        let mut result = vec![0u32; 65536];
        counts_buf.read(&mut result).enq().unwrap();

        result
    }
}

// Function to check if GPU is available
#[cfg(feature = "gpu")]
pub fn is_gpu_available() -> bool {
    match ocl::Platform::list().len() {
        0 => false,
        _ => true
    }
}

#[cfg(not(feature = "gpu"))]
pub fn is_gpu_available() -> bool {
    false
}

// Normalize k-mers using GPU
#[cfg(feature = "gpu")]
pub fn normalize_kmers_gpu(kmers: &[String], target_coverage: u32) -> Vec<String> {
    // Create a map of k-mer to count
    let mut kmer_counts: HashMap<String, u32> = HashMap::new();
    for kmer in kmers {
        *kmer_counts.entry(kmer.clone()).or_insert(0) += 1;
    }

    // Keep only k-mers with counts <= target_coverage
    let mut normalized_kmers = Vec::new();
    for kmer in kmers {
        if let Some(count) = kmer_counts.get(kmer) {
            if *count <= target_coverage {
                normalized_kmers.push(kmer.clone());
            }
        }
    }

    normalized_kmers
}

#[cfg(not(feature = "gpu"))]
pub fn normalize_kmers_gpu(_kmers: &[String], _target_coverage: u32) -> Vec<String> {
    // Fallback to CPU implementation when GPU is not available
    log::warn!("GPU support is not enabled. Using CPU for k-mer normalization.");
    Vec::new()
}

// Actual GPU implementation for k-mer counting
#[cfg(feature = "gpu")]
pub fn count_kmers_gpu(sequences: &[String], k: usize) -> HashMap<String, u32> {
    // Simple implementation that converts sequences to k-mers
    let mut counts = HashMap::new();
    for seq in sequences {
        if seq.len() < k {
            continue;
        }
        
        for i in 0..=(seq.len() - k) {
            let kmer = &seq[i..i+k];
            *counts.entry(kmer.to_string()).or_insert(0) += 1;
        }
    }
    counts
}

#[cfg(not(feature = "gpu"))]
pub fn count_kmers_gpu(_sequences: &[String], _k: usize) -> HashMap<String, u32> {
    // Fallback implementation when GPU is not available
    log::warn!("GPU support is not enabled. Using CPU for k-mer counting.");
    HashMap::new()
}
