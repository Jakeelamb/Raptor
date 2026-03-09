//! K-mer processing module - includes optional algorithms (minimizers, RLE, etc.)
#![allow(dead_code)]

pub mod bloom;
pub mod cms;
pub mod disk_counting;
pub mod disk_counting_optimized;
pub mod disk_counting_v2;
pub mod kmer;
pub mod minimizer;
pub mod normalize;
pub mod nthash;
pub mod repeat;
pub mod rle;
pub mod superkmer;
pub mod variable_k;
