//! K-mer processing module - includes optional algorithms (minimizers, RLE, etc.)
#![allow(dead_code)]

pub mod cms;
pub mod variable_k;
pub mod rle;
pub mod kmer;
pub mod normalize;
pub mod repeat;
pub mod nthash;
pub mod bloom;
pub mod minimizer;
pub mod disk_counting;
pub mod disk_counting_v2;
pub mod superkmer;
