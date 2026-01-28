pub mod kmer_gpu;
pub mod kmer_encode;
pub mod overlap_gpu;

pub use kmer_gpu::{is_gpu_available, get_gpu_info, count_kmers_cpu};
pub use kmer_encode::{encode_kmer, decode_kmer, canonical_encode, GpuKmerTable};
pub use overlap_gpu::find_overlaps_cpu;

#[cfg(feature = "gpu")]
pub use kmer_gpu::GpuKmerCounter;
#[cfg(feature = "gpu")]
pub use overlap_gpu::GpuOverlapFinder;
