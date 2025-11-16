//! BCF (Binary Call Format) reader.
//!
//! BCF is the binary representation of VCF (Variant Call Format), providing:
//! - **Compact storage**: ~5-10x smaller than text VCF
//! - **Fast parsing**: Binary encoding with typed values
//! - **Random access**: BGZF compression enables indexing
//! - **Lossless conversion**: Full fidelity with VCF format
//!
//! # Format Structure
//!
//! BCF files consist of:
//! - **Magic bytes**: `BCF\2\2` (version 2.2)
//! - **Text header**: VCF header with length prefix
//! - **Binary records**: BGZF-compressed variant records
//!
//! # Key Differences from VCF
//!
//! - **Binary encoding**: Typed values (int8/16/32, float, string)
//! - **Dictionary-based**: Field names referenced by integer indices
//! - **0-based coordinates**: Unlike VCF's 1-based positions
//! - **BGZF compressed**: Enables random access with CSI/TBI indices
//!
//! # Typed Value System
//!
//! BCF uses a type byte followed by data:
//! - **Lower 4 bits**: Type (1=int8, 2=int16, 3=int32, 5=float, 7=char)
//! - **Upper 4 bits**: Array length (15 = length follows)
//!
//! # Examples
//!
//! ## Basic usage
//!
//! ```no_run
//! use biometal::io::bcf::BcfReader;
//! use std::fs::File;
//!
//! # fn main() -> biometal::Result<()> {
//! let file = File::open("variants.bcf")?;
//! let mut reader = BcfReader::new(file)?;
//!
//! // Access header
//! println!("Samples: {:?}", reader.header().samples());
//!
//! // Parse records
//! for result in reader.records() {
//!     let record = result?;
//!     println!("{}\t{}\t{} -> {}",
//!              record.chrom(),
//!              record.pos(),
//!              record.reference(),
//!              record.alternate().join(","));
//! }
//! # Ok(())
//! # }
//! ```
//!
//! ## INFO field access
//!
//! ```no_run
//! # use biometal::io::bcf::BcfReader;
//! # use std::fs::File;
//! # fn main() -> biometal::Result<()> {
//! # let file = File::open("variants.bcf")?;
//! # let mut reader = BcfReader::new(file)?;
//! for result in reader.records() {
//!     let record = result?;
//!
//!     // Access INFO fields
//!     if let Some(dp) = record.info("DP")? {
//!         println!("Depth: {:?}", dp);
//!     }
//!     if let Some(af) = record.info("AF")? {
//!         println!("Allele frequency: {:?}", af);
//!     }
//! }
//! # Ok(())
//! # }
//! ```

use crate::error::{BiometalError, Result};
use crate::io::compression::{CompressedReader, DataSource};
use std::collections::HashMap;
use std::io::Read;
use std::path::Path;

mod header;
mod record;
mod typed_value;

pub use header::BcfHeader;
pub use record::BcfRecord;
pub use typed_value::{TypedValue, ValueType};

/// BCF magic bytes for version 2.2
const BCF_MAGIC: &[u8; 5] = b"BCF\x02\x02";

/// BCF reader with streaming API.
///
/// Provides constant-memory iteration over BCF records using BGZF decompression.
pub struct BcfReader {
    /// BGZF-compressed reader
    reader: CompressedReader,
    /// Parsed header with field dictionaries
    header: BcfHeader,
}

impl BcfReader {
    /// Create a new BCF reader from a file path.
    ///
    /// Reads and validates the magic bytes and header.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::bcf::BcfReader;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let reader = BcfReader::from_path("variants.bcf")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let source = DataSource::from_path(path.as_ref());
        Self::new(source)
    }

    /// Create a new BCF reader from a DataSource.
    pub fn new(source: DataSource) -> Result<Self> {
        let mut reader = CompressedReader::new(source)?;

        // Read and verify magic
        let mut magic = [0u8; 5];
        reader.read_exact(&mut magic)?;
        if &magic != BCF_MAGIC {
            return Err(BiometalError::InvalidInput {
                msg: format!(
                    "Invalid BCF magic: expected {:?}, got {:?}",
                    BCF_MAGIC, magic
                ),
            });
        }

        // Parse header
        let header = BcfHeader::parse(&mut reader)?;

        Ok(BcfReader { reader, header })
    }

    /// Get reference to the header.
    pub fn header(&self) -> &BcfHeader {
        &self.header
    }

    /// Create an iterator over BCF records.
    ///
    /// Returns a streaming iterator that yields records one at a time
    /// with constant memory usage.
    pub fn records(&mut self) -> BcfRecords<'_> {
        BcfRecords { reader: self }
    }
}

/// Iterator over BCF records.
///
/// Provides streaming access to BCF records with constant memory usage.
pub struct BcfRecords<'a> {
    reader: &'a mut BcfReader,
}

impl<'a> Iterator for BcfRecords<'a> {
    type Item = Result<BcfRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match BcfRecord::parse(&mut self.reader.reader, &self.reader.header) {
            Ok(record) => Some(Ok(record)),
            Err(BiometalError::Io(ref e)) if e.kind() == std::io::ErrorKind::UnexpectedEof => None,
            Err(e) => Some(Err(e)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bcf_magic_constant() {
        assert_eq!(BCF_MAGIC, b"BCF\x02\x02");
        assert_eq!(BCF_MAGIC.len(), 5);
    }
}
