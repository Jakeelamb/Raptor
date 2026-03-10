//! Shared stats helpers for CLI usage.
#![allow(dead_code)]

pub type Stats = crate::stats::Stats;

#[inline]
pub fn calculate_stats(path: &str) -> std::io::Result<Stats> {
    crate::stats::calculate_stats(path)
}
