//! I/O module: Streaming parsers and compression
//!
//! This module implements streaming architecture with constant memory (~5 MB)
//! regardless of dataset size, following Rules 2-6 from OPTIMIZATION_RULES.md.

pub mod compression;
pub mod fasta;
mod fastq;
mod paired;
pub mod sink;

pub use compression::{decompress_bgzip_parallel, CompressedReader, CompressedWriter, DataSource, STREAMING_THRESHOLD};
pub use fasta::{FaiIndex, FastaStream, FastaWriter};
pub use fastq::{FastqStream, FastqWriter};
pub use paired::PairedFastqStream;
pub use sink::DataSink;

// Week 3-4: Network streaming (Rule 6)
#[cfg(feature = "network")]
pub mod network;
#[cfg(feature = "network")]
pub use network::{HttpClient, HttpReader};

// Week 3-4: SRA integration (Rule 6)
#[cfg(feature = "network")]
pub mod sra;
#[cfg(feature = "network")]
pub use sra::{is_sra_accession, sra_to_url};

// Native BAM implementation (Phase 1-6)
// ARM-optimized BAM parser with parallel BGZF decompression
// See experiments/native-bam-implementation/ for design and profiling
pub mod bam;
pub use bam::BamReader;

// CRAM implementation (v1.11.0+)
// Reference-based compressed alignment format
// Uses noodles-cram for CRAM 3.0/3.1 support
pub mod cram;
pub use cram::CramReader;

// BCF implementation (v1.12.0+)
// Binary Variant Call Format - compressed VCF
// Native implementation with typed value system and BGZF compression
pub mod bcf;
pub use bcf::BcfReader;
