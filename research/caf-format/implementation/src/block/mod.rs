//! Block-level operations for CAF format.
//!
//! This module provides block builder and reader for converting between
//! row-oriented records and columnar blocks.
//!
//! # Block Size
//!
//! Default block size is 10,000 records (Rule 2 from OPTIMIZATION_RULES.md).
//! This balances:
//! - SIMD efficiency (amortizes overhead)
//! - Memory usage (~5 MB per block)
//! - Random access granularity

pub mod builder;
pub mod reader;

pub use builder::{BlockBuilder, AlignmentRecord};
pub use reader::BlockReader;

/// Default block size: 10,000 records (from OPTIMIZATION_RULES.md Rule 2).
pub const DEFAULT_BLOCK_SIZE: u32 = 10_000;
