#[cfg(feature = "gpu")]
use ocl::{Buffer, ProQue, flags};
use std::collections::HashMap;

#[cfg(feature = "gpu")]
use crate::gpu::kmer_encode::{decode_kmer, hash_kmer, canonical_encoded};

/// GPU k-mer counter using OpenCL
#[cfg(feature = "gpu")]
pub struct GpuKmerCounter {
    pro_que: ProQue,
    table_size: usize,
    k: usize,
}

#[cfg(feature = "gpu")]
impl GpuKmerCounter {
    /// Create a new GPU k-mer counter
    ///
    /// # Arguments
    /// * `k` - K-mer size (must be <= 32)
    /// * `max_reads` - Maximum number of reads to process
    /// * `table_size_bits` - Log2 of hash table size (e.g., 24 for 16M entries)
    pub fn new(k: usize, max_reads: usize, table_size_bits: usize) -> Result<Self, String> {
        if k > 32 {
            return Err("K-mer size must be <= 32 for 64-bit encoding".to_string());
        }

        let table_size = 1usize << table_size_bits;

        // Load and compile the kernel
        let kernel_src = include_str!("kernels/count_kmers_v2.cl");

        let pro_que = ProQue::builder()
            .src(kernel_src)
            .dims(max_reads)
            .build()
            .map_err(|e| format!("Failed to build OpenCL program: {}", e))?;

        Ok(GpuKmerCounter { pro_que, table_size, k })
    }

    /// Count k-mers in the given reads
    ///
    /// Returns a HashMap of k-mer strings to counts
    pub fn count(&self, reads: &[String]) -> Result<HashMap<String, u32>, String> {
        if reads.is_empty() {
            return Ok(HashMap::new());
        }

        // Flatten reads into a single byte array
        let sequence_data: Vec<u8> = reads.iter().flat_map(|s| s.bytes()).collect();

        // Build offset array
        let mut offsets = Vec::with_capacity(reads.len() + 1);
        let mut pos = 0i32;
        offsets.push(pos);
        for r in reads {
            pos += r.len() as i32;
            offsets.push(pos);
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

        // Hash table buffers
        let counts_buf = Buffer::<u32>::builder()
            .queue(self.pro_que.queue().clone())
            .flags(flags::MEM_READ_WRITE)
            .len(self.table_size)
            .fill_val(0u32)
            .build()
            .map_err(|e| format!("Failed to create counts buffer: {}", e))?;

        let kmers_buf = Buffer::<u64>::builder()
            .queue(self.pro_que.queue().clone())
            .flags(flags::MEM_READ_WRITE)
            .len(self.table_size)
            .fill_val(0u64)
            .build()
            .map_err(|e| format!("Failed to create k-mers buffer: {}", e))?;

        let table_mask = (self.table_size - 1) as u32;

        // Build and execute kernel
        let kernel = self.pro_que
            .kernel_builder("count_kmers_v2")
            .arg(&sequences_buf)
            .arg(&offsets_buf)
            .arg(&counts_buf)
            .arg(&kmers_buf)
            .arg(self.k as u32)
            .arg(reads.len() as u32)
            .arg(self.table_size as u32)
            .arg(table_mask)
            .build()
            .map_err(|e| format!("Failed to build kernel: {}", e))?;

        unsafe {
            kernel.enq().map_err(|e| format!("Failed to execute kernel: {}", e))?;
        }

        // Read back results
        let mut counts = vec![0u32; self.table_size];
        let mut kmers = vec![0u64; self.table_size];

        counts_buf.read(&mut counts)
            .enq()
            .map_err(|e| format!("Failed to read counts: {}", e))?;

        kmers_buf.read(&mut kmers)
            .enq()
            .map_err(|e| format!("Failed to read k-mers: {}", e))?;

        // Convert to HashMap
        let mut result = HashMap::new();
        for (kmer_encoded, &count) in kmers.iter().zip(counts.iter()) {
            if count > 0 && *kmer_encoded != 0 {
                let kmer_str = decode_kmer(*kmer_encoded, self.k);
                result.insert(kmer_str, count);
            }
        }

        Ok(result)
    }

    /// Count k-mers using the simpler kernel (for smaller datasets)
    pub fn count_simple(&self, reads: &[String]) -> Result<Vec<u32>, String> {
        if reads.is_empty() {
            return Ok(vec![0u32; self.table_size]);
        }

        let sequence_data: Vec<u8> = reads.iter().flat_map(|s| s.bytes()).collect();

        let mut offsets = Vec::with_capacity(reads.len() + 1);
        let mut pos = 0i32;
        offsets.push(pos);
        for r in reads {
            pos += r.len() as i32;
            offsets.push(pos);
        }

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

        let counts_buf = Buffer::<u32>::builder()
            .queue(self.pro_que.queue().clone())
            .flags(flags::MEM_READ_WRITE)
            .len(self.table_size)
            .fill_val(0u32)
            .build()
            .map_err(|e| format!("Failed to create counts buffer: {}", e))?;

        let table_mask = (self.table_size - 1) as u32;

        let kernel = self.pro_que
            .kernel_builder("count_kmers_simple")
            .arg(&sequences_buf)
            .arg(&offsets_buf)
            .arg(&counts_buf)
            .arg(self.k as u32)
            .arg(reads.len() as u32)
            .arg(self.table_size as u32)
            .arg(table_mask)
            .build()
            .map_err(|e| format!("Failed to build kernel: {}", e))?;

        unsafe {
            kernel.enq().map_err(|e| format!("Failed to execute kernel: {}", e))?;
        }

        let mut result = vec![0u32; self.table_size];
        counts_buf.read(&mut result)
            .enq()
            .map_err(|e| format!("Failed to read counts: {}", e))?;

        Ok(result)
    }
}

/// Check if GPU is available
#[cfg(feature = "gpu")]
pub fn is_gpu_available() -> bool {
    match ocl::Platform::list() {
        platforms if !platforms.is_empty() => {
            // Check if any platform has devices
            platforms.iter().any(|p| {
                ocl::Device::list_all(p).map(|d| !d.is_empty()).unwrap_or(false)
            })
        }
        _ => false,
    }
}

#[cfg(not(feature = "gpu"))]
pub fn is_gpu_available() -> bool {
    false
}

/// Get GPU device information
#[cfg(feature = "gpu")]
pub fn get_gpu_info() -> Option<String> {
    let platforms = ocl::Platform::list();
    if platforms.is_empty() {
        return None;
    }

    let mut info = String::new();
    for platform in platforms {
        if let Ok(name) = platform.name() {
            info.push_str(&format!("Platform: {}\n", name));
        }
        if let Ok(devices) = ocl::Device::list_all(&platform) {
            for device in devices {
                if let Ok(name) = device.name() {
                    info.push_str(&format!("  Device: {}\n", name));
                }
            }
        }
    }

    if info.is_empty() {
        None
    } else {
        Some(info)
    }
}

#[cfg(not(feature = "gpu"))]
pub fn get_gpu_info() -> Option<String> {
    None
}

/// CPU fallback for k-mer counting (used when GPU is not available)
pub fn count_kmers_cpu(sequences: &[String], k: usize) -> HashMap<String, u32> {
    use crate::kmer::kmer::canonical_kmer;

    let mut counts = HashMap::new();
    for seq in sequences {
        if seq.len() < k {
            continue;
        }

        for i in 0..=seq.len() - k {
            if let Some(kmer) = canonical_kmer(&seq[i..i + k]) {
                *counts.entry(kmer).or_insert(0) += 1;
            }
        }
    }
    counts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cpu_fallback() {
        let sequences = vec![
            "ACGTACGT".to_string(),
            "CGTACGTA".to_string(),
        ];

        let counts = count_kmers_cpu(&sequences, 4);
        assert!(!counts.is_empty());
        // ACGT should appear multiple times
        assert!(counts.contains_key("ACGT"));
    }

    #[test]
    fn test_gpu_availability() {
        // This test just checks that the function doesn't crash
        let _available = is_gpu_available();
    }
}
