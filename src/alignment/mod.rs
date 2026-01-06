//! Sequence alignment algorithms and pattern matching workflows
//!
//! This module provides implementations of sequence alignment algorithms
//! including Smith-Waterman (local alignment) with CPU, NEON, and Metal GPU
//! implementations, plus streaming pattern matching workflows for bioinformatics analysis.
//!
//! # Core Components
//!
//! ## Smith-Waterman Alignment
//! Three implementations with automatic dispatch:
//! - **Naive CPU**: Reference implementation (correctness baseline)
//! - **NEON CPU**: ARM SIMD optimization (2-4× speedup, portable to Graviton)
//! - **Metal GPU**: Apple Silicon exclusive (5-10× speedup, batch processing)
//!
//! ## Streaming Read Mapping
//! Novel windowed reference processing maintaining ~5MB constant memory:
//! - Processes reference genome in overlapping windows
//! - Enables read mapping without large memory indices
//! - Constant memory regardless of reference/dataset size
//!
//! ## Pattern Matching Workflows
//! Streaming pattern matching for practical bioinformatics:
//! - **Motif Finding**: Regulatory sequences, transcription factor binding sites
//! - **Primer Detection**: PCR primer sites, qPCR assays
//! - **Adapter Detection**: Illumina sequencing adapters, contamination removal
//!
//! # Evidence Base
//!
//! - GPU Smith-Waterman: Expected 10-50× speedup (CUDA literature)
//! - NEON Smith-Waterman: Expected 2-4× speedup (Rule 1, but limited by dependencies)
//! - Complexity >0.70: Dynamic programming dominates (exceeds ASBB GPU threshold)
//! - Streaming architecture: 99.5% memory reduction vs traditional aligners
//!
//! # Examples
//!
//! ## Basic Alignment
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
//!
//! ## Streaming Read Mapping
//! ```no_run
//! use biometal::alignment::{StreamingMapper, StreamingMapperConfig};
//! use std::path::Path;
//!
//! let config = StreamingMapperConfig::default();
//! let mut mapper = StreamingMapper::new(config);
//!
//! // Map reads with constant memory usage
//! for mapping in mapper.map_reads_streaming(
//!     Path::new("reference.fasta"),
//!     Path::new("reads.fastq")
//! )? {
//!     let result = mapping?;
//!     println!("Mapped {} to position {}", result.query_id, result.global_ref_start);
//! }
//! ```
//!
//! ## Pattern Matching
//! ```no_run
//! use biometal::alignment::{MotifFinder, MotifPattern, ScoringMatrix};
//! use biometal::FastqStream;
//! use std::path::Path;
//!
//! let motifs = vec![
//!     MotifPattern::new("TATAAA", "TATA box", 40),
//!     MotifPattern::new("CAAT", "CAAT box", 30),
//! ];
//!
//! let mut finder = MotifFinder::new(motifs, ScoringMatrix::default());
//! let stream = FastqStream::from_path("sequences.fastq")?;
//!
//! for motif_match in finder.find_motifs_streaming(stream) {
//!     let m = motif_match?;
//!     println!("Found {} in {}: score {}", m.motif_name, m.sequence_id, m.alignment.score);
//! }
//! ```

pub mod cigar;
pub mod scoring;
pub mod smith_waterman;
pub mod streaming_mapper;
pub mod pattern_matcher;

#[cfg(feature = "gpu")]
pub mod gpu;

// Re-export public API
pub use cigar::{CigarOp, compress_cigar};
pub use scoring::ScoringMatrix;
pub use smith_waterman::{smith_waterman, smith_waterman_naive, Alignment};
pub use streaming_mapper::{
    StreamingMapper, StreamingMapperConfig, MappingResult, StreamingMappingIterator,
};
pub use pattern_matcher::{
    MotifFinder, PrimerFinder, AdapterDetector, MotifPattern, MotifMatch,
};

#[cfg(target_arch = "aarch64")]
pub use smith_waterman::smith_waterman_neon;

#[cfg(feature = "gpu")]
pub use gpu::{smith_waterman_batch_gpu, GpuAlignmentBatch};

#[cfg(feature = "gpu")]
pub use smith_waterman::smith_waterman_gpu;
