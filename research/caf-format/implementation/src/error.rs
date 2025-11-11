//! Error types for CAF operations.
//!
//! This module defines all error conditions that can occur during CAF file
//! parsing, writing, conversion, and validation.

use std::io;
use thiserror::Error;

/// CAF-specific error types.
///
/// All errors provide detailed context including file offsets, block IDs,
/// and specific failure information to aid debugging.
#[derive(Debug, Error)]
pub enum CafError {
    /// Invalid magic number at start of file
    #[error("Invalid magic number: expected CAF\\x01, got {0:?}")]
    InvalidMagic([u8; 4]),

    /// Unsupported CAF version
    #[error("Unsupported CAF version {major}.{minor} (this reader supports 1.x)")]
    UnsupportedVersion {
        /// Major version number
        major: u8,
        /// Minor version number
        minor: u8,
    },

    /// CRC32 checksum validation failure
    #[error(
        "Checksum mismatch in block {block_id}: expected 0x{expected:08x}, got 0x{actual:08x}"
    )]
    ChecksumMismatch {
        /// Block ID where mismatch occurred
        block_id: u32,
        /// Expected checksum value
        expected: u32,
        /// Actual computed checksum
        actual: u32,
    },

    /// Decompression failure
    #[error("Decompression failed for column '{column}' in block {block_id}: {source}")]
    DecompressionError {
        /// Column name (e.g., "sequences", "qualities")
        column: String,
        /// Block ID
        block_id: u32,
        /// Underlying compression error
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },

    /// Compression failure
    #[error("Compression failed for column '{column}' in block {block_id}: {source}")]
    CompressionError {
        /// Column name
        column: String,
        /// Block ID
        block_id: u32,
        /// Underlying compression error
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },

    /// Index structure inconsistency
    #[error(
        "Index inconsistency: header claims {header_blocks} blocks, index has {index_blocks}"
    )]
    IndexInconsistency {
        /// Number of blocks claimed in header
        header_blocks: u32,
        /// Number of blocks in index
        index_blocks: u32,
    },

    /// Block size mismatch
    #[error("Block {block_id} size mismatch: expected {expected} bytes, got {actual} bytes")]
    BlockSizeMismatch {
        /// Block ID
        block_id: u32,
        /// Expected size (from metadata)
        expected: u32,
        /// Actual size (decompressed)
        actual: u32,
    },

    /// Invalid block ID
    #[error("Invalid block ID {block_id}: file only has {num_blocks} blocks")]
    InvalidBlockId {
        /// Requested block ID
        block_id: u32,
        /// Total blocks in file
        num_blocks: u32,
    },

    /// Invalid region query
    #[error("Invalid region query: {message}")]
    InvalidRegion {
        /// Error description
        message: String,
    },

    /// Serialization/deserialization error
    #[error("Serialization error: {0}")]
    Serialization(#[from] bincode::Error),

    /// I/O error
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),

    /// BAM conversion error
    #[error("BAM conversion error: {message}")]
    BamConversion {
        /// Error description
        message: String,
    },

    /// Lossless conversion validation failure
    #[error("Lossless conversion failed: record {record_id} differs after round-trip")]
    LosslessValidationFailed {
        /// Record index that differs
        record_id: usize,
    },

    /// Column encoding error
    #[error("Column encoding error for '{column}': {message}")]
    ColumnEncoding {
        /// Column name
        column: String,
        /// Error description
        message: String,
    },

    /// Generic error for unexpected conditions
    #[error("{0}")]
    Other(String),
}

/// Specialized Result type for CAF operations.
pub type Result<T> = std::result::Result<T, CafError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_messages() {
        let err = CafError::InvalidMagic([b'B', b'A', b'M', 0x01]);
        assert!(err.to_string().contains("Invalid magic number"));

        let err = CafError::UnsupportedVersion { major: 2, minor: 0 };
        assert!(err.to_string().contains("Unsupported CAF version 2.0"));

        let err = CafError::ChecksumMismatch {
            block_id: 42,
            expected: 0xDEADBEEF,
            actual: 0xBADC0FFE,
        };
        assert!(err.to_string().contains("block 42"));
        assert!(err.to_string().contains("0xdeadbeef"));
    }
}
