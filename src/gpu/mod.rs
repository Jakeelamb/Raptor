pub mod kmer_gpu;
pub mod kmer_encode;
pub mod overlap_gpu;

// Public API exports - some may be unused internally but available for library users
#[allow(unused_imports)]
pub use kmer_gpu::{is_gpu_available, get_gpu_info, count_kmers_cpu, count_kmers_nthash_cpu};
#[allow(unused_imports)]
pub use kmer_encode::{encode_kmer, decode_kmer, canonical_encode, GpuKmerTable};
#[allow(unused_imports)]
pub use overlap_gpu::find_overlaps_cpu;

#[cfg(feature = "gpu")]
pub use kmer_gpu::{GpuKmerCounter, GpuKmerCounterNtHash};
#[cfg(feature = "gpu")]
pub use overlap_gpu::GpuOverlapFinder;
