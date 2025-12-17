//! CLI command modules for biometal toolkit
//!
//! This module organizes biometal's command-line interface into logical groups:
//!
//! - `statistics`: NEON-optimized sequence statistics (16.7-25.1× speedup)
//! - `sequence`: Core sequence transformations (reverse complement, etc.)
//! - `quality`: Quality-based operations (trimming, masking, filtering)
//! - `conversion`: Format conversion utilities (FASTQ/FASTA, counting)
//! - `pattern`: Pattern matching operations (find, count, multi-pattern)
//! - `kmer`: K-mer operations (minimizers, spectrum analysis)
//!
//! # Design Principles
//!
//! - Unix philosophy: One command = one primitive function
//! - Comprehensive I/O: Support files, stdin/stdout, HTTP/HTTPS, SRA
//! - Performance-first: Prioritize NEON-optimized operations
//! - Evidence-based: Commands organized by proven performance characteristics

pub mod conversion;
pub mod kmer;
pub mod pattern;
pub mod quality;
pub mod sequence;
pub mod statistics;