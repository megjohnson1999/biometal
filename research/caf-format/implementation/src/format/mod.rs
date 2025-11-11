//! Binary format handling for CAF files.
//!
//! This module handles low-level binary format parsing and writing:
//! - Magic number validation (`magic`)
//! - Header parsing/writing (`header`)
//! - Index structure (`index`)
//! - Footer structure (`footer`)
//!
//! # Example
//!
//! ```rust,no_run
//! use std::fs::File;
//! use caf::format::{read_magic, read_header, read_index, read_footer};
//!
//! # fn main() -> Result<(), caf::error::CafError> {
//! let mut file = File::open("data.caf")?;
//!
//! // Read in order: magic, header, blocks (not shown), index, footer
//! let version = read_magic(&mut file)?;
//! let header = read_header(&mut file)?;
//! // ... read blocks ...
//! let index = read_index(&mut file)?;
//! let footer = read_footer(&mut file)?;
//! # Ok(())
//! # }
//! ```

pub mod header;
pub mod index;
pub mod footer;
pub mod magic;

// Re-exports
pub use header::{read_header, write_header};
pub use index::{read_index, write_index};
pub use footer::{read_footer, write_footer};
pub use magic::{read_magic, write_magic, validate_magic};
