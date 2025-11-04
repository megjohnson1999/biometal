//! Error types for biometal

use std::fmt;

/// Result type alias for biometal operations
pub type Result<T> = std::result::Result<T, BiometalError>;

/// Error types that can occur in biometal
#[derive(Debug)]
pub enum BiometalError {
    /// I/O error
    Io(std::io::Error),
    
    /// Invalid FASTQ format
    InvalidFastqFormat {
        /// Line number where error occurred
        line: usize,
        /// Error message
        msg: String,
    },
    
    /// Invalid FASTA format
    InvalidFastaFormat {
        /// Line number where error occurred
        line: usize,
        /// Error message
        msg: String,
    },
    
    /// Compression/decompression error
    Compression(String),
    
    /// Network error (Week 3-4)
    #[cfg(feature = "network")]
    Network(String),
}

impl fmt::Display for BiometalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BiometalError::Io(e) => write!(f, "I/O error: {}", e),
            BiometalError::InvalidFastqFormat { line, msg } => {
                write!(f, "Invalid FASTQ format at line {}: {}", line, msg)
            }
            BiometalError::InvalidFastaFormat { line, msg } => {
                write!(f, "Invalid FASTA format at line {}: {}", line, msg)
            }
            BiometalError::Compression(msg) => write!(f, "Compression error: {}", msg),
            #[cfg(feature = "network")]
            BiometalError::Network(msg) => write!(f, "Network error: {}", msg),
        }
    }
}

impl std::error::Error for BiometalError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            BiometalError::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for BiometalError {
    fn from(error: std::io::Error) -> Self {
        BiometalError::Io(error)
    }
}
