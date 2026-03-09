pub mod kmer_encode;
pub mod kmer_gpu;
pub mod overlap_gpu;

// Public API exports - some may be unused internally but available for library users
#[allow(unused_imports)]
pub use kmer_encode::{canonical_encode, decode_kmer, encode_kmer, GpuKmerTable};
#[allow(unused_imports)]
pub use kmer_gpu::{count_kmers_cpu, count_kmers_nthash_cpu, get_gpu_info, is_gpu_available};
#[allow(unused_imports)]
pub use overlap_gpu::find_overlaps_cpu;

#[cfg(feature = "gpu")]
#[allow(unused_imports)]
pub use kmer_gpu::{GpuKmerCounter, GpuKmerCounterNtHash};
#[cfg(feature = "gpu")]
#[allow(unused_imports)]
pub use overlap_gpu::GpuOverlapFinder;
