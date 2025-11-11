//! File I/O (CafReader, CafWriter implementations).
//!
//! This module provides high-level interfaces for reading and writing CAF files:
//! - `CafFileWriter`: Create CAF files with automatic block management
//! - `CafFileReader`: Read CAF files with random access and region queries
//!
//! # Example
//!
//! ```no_run
//! use caf::io::CafFileWriter;
//! use caf::block::AlignmentRecord;
//! use std::path::Path;
//!
//! # fn main() -> Result<(), caf::CafError> {
//! // Write a CAF file
//! let mut writer = CafFileWriter::create(Path::new("output.caf"))?;
//! writer.add_record(AlignmentRecord::default())?;
//! writer.finalize()?;
//! # Ok(())
//! # }
//! ```

pub mod writer;
pub mod reader;

pub use writer::CafFileWriter;
pub use reader::CafFileReader;

use crate::{CafBlock, CafHeader, Result};

/// CAF file reader trait (deprecated - use CafFileReader).
pub trait CafReader {
    /// Read the file header.
    fn read_header(&mut self) -> Result<CafHeader>;

    /// Read a specific block.
    fn read_block(&mut self, block_id: usize) -> Result<CafBlock>;

    /// Query blocks overlapping a genomic region.
    fn query_region(&self, chr: &str, start: i32, end: i32) -> Result<Vec<CafBlock>>;
}

/// CAF file writer trait (deprecated - use CafFileWriter).
pub trait CafWriter {
    /// Write the file header.
    fn write_header(&mut self, header: &CafHeader) -> Result<()>;

    /// Write a block.
    fn write_block(&mut self, block: &CafBlock) -> Result<()>;

    /// Finalize the file (write index + footer).
    fn finalize(&mut self) -> Result<()>;
}
