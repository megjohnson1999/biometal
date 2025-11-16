//! Index formats for genomic data files
//!
//! This module provides support for index formats that enable random access
//! to genomic data files:
//!
//! - **TBI (Tabix)**: Index for tab-delimited files (BED, VCF, GFF3)
//! - **CSI**: Coordinate-sorted index (successor to BAI/TBI, supports larger chromosomes)
//!
//! # Overview
//!
//! Index formats enable O(log n) region queries on sorted, BGZF-compressed
//! genomic files. They work with:
//! - BED (genomic intervals)
//! - VCF (variant calls)
//! - GFF3 (gene features)
//! - BAM (alignments)
//! - Generic tab-delimited files
//!
//! ## TBI vs CSI
//!
//! - **TBI**: Standard tabix index, limited to 512 Mbp chromosomes
//! - **CSI**: Extended index with configurable binning, supports larger references
//!
//! # Example - TBI
//!
//! ```no_run
//! use biometal::formats::index::TbiIndex;
//!
//! # fn main() -> biometal::Result<()> {
//! // Load tabix index
//! let index = TbiIndex::from_path("variants.vcf.gz.tbi")?;
//!
//! // Query region
//! let chunks = index.query("chr1", 1000000, 2000000)?;
//! println!("Found {} chunks for region", chunks.len());
//! # Ok(())
//! # }
//! ```
//!
//! # Example - CSI
//!
//! ```no_run
//! use biometal::formats::index::CsiIndex;
//!
//! # fn main() -> biometal::Result<()> {
//! // Load CSI index
//! let index = CsiIndex::from_path("alignments.bam.csi")?;
//!
//! // Query region by index (if names not available)
//! if let Some(chunks) = index.query_by_index(0, 1000000, 2000000)? {
//!     println!("Found {} chunks", chunks.len());
//! }
//! # Ok(())
//! # }
//! ```

pub mod csi;
pub mod tbi;

pub use csi::CsiIndex;
pub use tbi::TbiIndex;
