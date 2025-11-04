//! I/O module: Streaming parsers and compression
//!
//! This module implements streaming architecture with constant memory (~5 MB)
//! regardless of dataset size, following Rules 2-6 from OPTIMIZATION_RULES.md.

mod fastq;

pub use fastq::FastqStream;

// Week 1-2: Implement compression module (Rules 3-4)
// mod compression;
// pub use compression::{BgzipReader, parallel_decompress};

// Week 3-4: Implement network streaming (Rule 6)
// #[cfg(feature = "network")]
// mod network;
// #[cfg(feature = "network")]
// pub use network::NetworkStream;
