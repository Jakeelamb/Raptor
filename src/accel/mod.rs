pub mod simd;
pub mod gpu;
pub mod backend;
pub mod cpu_backend;
pub mod gpu_backend;

// Public API exports - some may be unused internally but available for library users
#[allow(unused_imports)]
pub use backend::{ComputeBackend, AdjacencyTable, GPU_THRESHOLD_READS};
pub use cpu_backend::CpuBackend;
#[allow(unused_imports)]
pub use gpu_backend::GpuBackend;

/// Create the appropriate backend based on configuration
///
/// # Arguments
/// * `use_gpu` - Whether to prefer GPU acceleration
/// * `num_sequences` - Number of sequences to process (for threshold check)
/// * `max_contigs` - Expected maximum number of contigs
///
/// # Returns
/// Box containing the selected backend
pub fn create_backend(
    use_gpu: bool,
    num_sequences: usize,
    _max_contigs: usize,
) -> Box<dyn ComputeBackend> {
    if use_gpu && num_sequences >= GPU_THRESHOLD_READS {
        #[cfg(feature = "gpu")]
        {
            if let Some(gpu) = gpu_backend::try_create_gpu_backend(num_sequences, max_contigs) {
                tracing::info!(
                    "Using GPU backend for {} sequences (threshold: {})",
                    num_sequences,
                    GPU_THRESHOLD_READS
                );
                return Box::new(gpu);
            }
            tracing::warn!("GPU requested but not available, using CPU backend");
        }

        #[cfg(not(feature = "gpu"))]
        {
            tracing::warn!("GPU requested but not compiled in (--features gpu), using CPU backend");
        }
    } else if use_gpu && num_sequences < GPU_THRESHOLD_READS {
        tracing::info!(
            "Dataset too small ({} < {} threshold), using CPU backend instead of GPU",
            num_sequences,
            GPU_THRESHOLD_READS
        );
    }

    Box::new(CpuBackend::new())
}

/// Create a CPU backend (convenience function)
#[allow(dead_code)]
pub fn create_cpu_backend() -> Box<dyn ComputeBackend> {
    Box::new(CpuBackend::new())
}

/// Create a CPU backend with specific thread count
#[allow(dead_code)]
pub fn create_cpu_backend_with_threads(threads: usize) -> Box<dyn ComputeBackend> {
    Box::new(CpuBackend::with_threads(threads))
}
