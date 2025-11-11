//! # CAF (Columnar Alignment Format) Library
//!
//! CAF is a columnar binary format for DNA/RNA sequence alignments optimized for
//! ARM NEON SIMD operations and analytical bioinformatics workloads.
//!
//! ## Design Goals
//!
//! - **Performance**: 5-10× speedup over BAM for analytical operations
//! - **Correctness**: 100% lossless BAM ↔ CAF conversion
//! - **Modern**: zstd/lz4 compression, ARM NEON optimization
//! - **Evidence-based**: Design from biometal's OPTIMIZATION_RULES.md
//!
//! ## Quick Start
//!
//! ```rust,no_run
//! use caf::{CafReader, CafWriter};
//! use std::path::Path;
//!
//! # fn main() -> Result<(), caf::error::CafError> {
//! // Convert BAM → CAF
//! let bam_path = Path::new("input.bam");
//! let caf_path = Path::new("output.caf");
//! // bam_to_caf(bam_path, caf_path)?;
//! # Ok(())
//! # }
//! ```
//!
//! ## Module Structure
//!
//! - [`error`]: Error types and Result alias
//! - [`types`]: Core data structures (CafHeader, CafBlock, etc.)
//! - [`format`]: Binary format parsing/writing
//! - [`block`]: Columnar block operations
//! - [`column`]: Column-specific encoding
//! - [`compression`]: Compression strategies (zstd, lz4, RLE)
//! - [`conversion`]: BAM ↔ CAF conversion
//! - [`io`]: File I/O (CafReader, CafWriter traits)
//! - [`query`]: Region queries and filtering
//! - [`neon`]: ARM NEON optimizations (feature-gated)
//! - [`validation`]: Checksums and round-trip testing
//!
//! ## Features
//!
//! - `neon`: Enable ARM NEON SIMD optimizations (default on aarch64)
//! - `parallel`: Enable parallel block decompression
//!
//! ## Performance
//!
//! CAF achieves significant speedups for analytical operations:
//!
//! | Operation | BAM | CAF (NEON) | Speedup |
//! |-----------|-----|------------|---------|
//! | Quality filter Q30 | 2.00s | 0.08s | **25×** |
//! | Base counting | 1.50s | 0.06s | **25×** |
//! | MAPQ > 30 filter | 0.50s | 0.03s | **16×** |
//! | Parse 100K records | 2.56s | 0.25s | **10×** |
//!
//! (Targets, validation pending - see benchmarks/)
//!
//! ## References
//!
//! - Specification: `SPECIFICATION.md`
//! - Literature Review: `LITERATURE_REVIEW.md` (28 citations)
//! - Research Plan: `RESEARCH_PLAN.md`

// Public API
pub mod error;
pub mod types;

// Core modules
pub mod format;
pub mod block;
pub mod column;
pub mod compression;
pub mod conversion;
pub mod io;
pub mod query;
pub mod validation;

// Platform-specific optimizations
#[cfg(target_arch = "aarch64")]
pub mod neon;

// Re-exports for convenience
pub use error::{CafError, Result};
pub use types::{CafHeader, CafBlock, CafIndex, CafFooter, BlockMeta};
pub use block::{BlockBuilder, BlockReader, AlignmentRecord};
pub use io::{CafReader, CafWriter, CafFileWriter, CafFileReader};

// Conversion utilities
pub use conversion::{bam_to_caf, caf_to_sam};

/// CAF format version
pub const CAF_VERSION_MAJOR: u8 = 1;
pub const CAF_VERSION_MINOR: u8 = 0;

/// Default block size (10,000 records from OPTIMIZATION_RULES.md Rule 2)
pub const DEFAULT_BLOCK_SIZE: u32 = 10_000;

/// Magic number: "CAF\x01"
pub const CAF_MAGIC: [u8; 4] = [b'C', b'A', b'F', 0x01];

/// Footer magic: "CAFE"
pub const CAF_FOOTER_MAGIC: [u8; 4] = [b'C', b'A', b'F', b'E'];
