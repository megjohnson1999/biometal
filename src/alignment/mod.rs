//! Sequence alignment algorithms
//!
//! This module provides implementations of sequence alignment algorithms
//! including Smith-Waterman (local alignment) with CPU, NEON, and Metal GPU
//! implementations.
//!
//! # Evidence Base
//!
//! - GPU Smith-Waterman: Expected 10-50× speedup (CUDA literature)
//! - NEON Smith-Waterman: Expected 2-4× speedup (Rule 1, but limited by dependencies)
//! - Complexity >0.70: Dynamic programming dominates (exceeds ASBB GPU threshold)
//!
//! # Architecture
//!
//! Three implementations with automatic dispatch:
//! - **Naive CPU**: Reference implementation (correctness baseline)
//! - **NEON CPU**: ARM SIMD optimization (2-4× speedup, portable to Graviton)
//! - **Metal GPU**: Apple Silicon exclusive (10-50× speedup, batch processing)
//!
//! # Examples
//!
//! ```
//! use biometal::alignment::{smith_waterman, ScoringMatrix};
//!
//! let query = b"ACGTACGT";
//! let reference = b"ACGTACGT";
//! let scoring = ScoringMatrix::default();
//!
//! let alignment = smith_waterman(query, reference, &scoring);
//! assert_eq!(alignment.score, 16); // 8 matches × 2 = 16
//! ```

pub mod cigar;
pub mod scoring;
pub mod smith_waterman;

#[cfg(feature = "gpu")]
pub mod gpu;

// Re-export public API
pub use cigar::{CigarOp, compress_cigar};
pub use scoring::ScoringMatrix;
pub use smith_waterman::{smith_waterman, smith_waterman_naive, Alignment};

#[cfg(target_arch = "aarch64")]
pub use smith_waterman::smith_waterman_neon;

#[cfg(feature = "gpu")]
pub use gpu::{smith_waterman_batch_gpu, GpuAlignmentBatch};

#[cfg(feature = "gpu")]
pub use smith_waterman::smith_waterman_gpu;
