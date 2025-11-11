//! ARM NEON SIMD optimizations for CAF operations.
//!
//! # Evidence
//!
//! Entry 020-025 (apple-silicon-bio-bench):
//! - **Base counting**: 25× speedup (Cohen's d = 5.87)
//! - **Quality filtering**: 25× speedup (Cohen's d = 5.87)
//! - **MAPQ filtering**: 16× speedup (Cohen's d = 4.82)
//!
//! # Architecture
//!
//! This module provides NEON-optimized operations on CAF columnar data:
//! - Process 16 elements at a time with SIMD instructions
//! - Automatic fallback to scalar on non-ARM platforms
//! - Property-based testing ensures NEON == scalar

pub mod base_counting;
pub mod quality_filter;
pub mod mapq_filter;

// Re-export main functions
pub use base_counting::{count_bases, count_bases_scalar, BaseCounts};
pub use quality_filter::{mean_quality, mean_quality_scalar, filter_records_by_quality};
pub use mapq_filter::{filter_records_by_mapq, count_high_mapq, count_high_mapq_scalar};

// Re-export NEON-specific functions (ARM only)
#[cfg(target_arch = "aarch64")]
pub use base_counting::count_bases_neon;

#[cfg(target_arch = "aarch64")]
pub use quality_filter::mean_quality_neon;

#[cfg(target_arch = "aarch64")]
pub use mapq_filter::count_high_mapq_neon;
