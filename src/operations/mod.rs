//! ARM NEON-optimized operations
//!
//! This module implements Rule 1 (ARM NEON SIMD) from OPTIMIZATION_RULES.md,
//! providing 16-25× speedup on ARM platforms with automatic scalar fallback.

pub mod base_counting;

pub use base_counting::count_bases;

// Week 1-2: Additional operations
// mod gc_content;
// mod quality_filter;
// pub use gc_content::gc_content;
// pub use quality_filter::filter_by_quality;
