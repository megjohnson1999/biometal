//! I/O module: Streaming parsers and compression
//!
//! This module implements streaming architecture with constant memory (~5 MB)
//! regardless of dataset size, following Rules 2-6 from OPTIMIZATION_RULES.md.

pub mod compression;
mod fastq;

pub use compression::{decompress_bgzip_parallel, CompressedReader, DataSource, MMAP_THRESHOLD};
pub use fastq::FastqStream;

// Week 3-4: Implement network streaming (Rule 6)
// #[cfg(feature = "network")]
// mod network;
// #[cfg(feature = "network")]
// pub use network::NetworkStream;
