//! CAF footer parsing and writing.
//!
//! The footer is stored at the end of the CAF file and contains:
//! - Index offset (for fast seeking)
//! - Total blocks and records
//! - Full-file checksum
//! - Magic number "CAFE"

use crate::{types::CafFooter, CafError, Result};
use std::io::{Read, Write};

/// Read CAF footer from a reader.
///
/// The footer is always the last 32 bytes of a CAF file.
///
/// # Errors
///
/// - `CafError::InvalidMagic` if footer magic doesn't match "CAFE"
/// - `CafError::Io` for I/O errors
///
/// # Example
///
/// ```rust,no_run
/// use std::fs::File;
/// use std::io::Seek;
/// # use caf::format::read_footer;
/// # fn main() -> Result<(), caf::error::CafError> {
/// let mut file = File::open("data.caf")?;
/// // Seek to last 32 bytes
/// let footer = read_footer(&mut file)?;
/// println!("Index offset: {}", footer.index_offset);
/// # Ok(())
/// # }
/// ```
pub fn read_footer<R: Read>(reader: &mut R) -> Result<CafFooter> {
    let footer: CafFooter = bincode::deserialize_from(reader)?;

    // Validate footer magic
    if !footer.validate_magic() {
        return Err(CafError::InvalidMagic(footer.magic));
    }

    Ok(footer)
}

/// Write CAF footer to a writer.
///
/// # Example
///
/// ```rust,no_run
/// use std::fs::File;
/// use caf::types::CafFooter;
/// # use caf::format::write_footer;
/// # fn main() -> Result<(), caf::error::CafError> {
/// let mut file = File::create("data.caf")?;
/// let footer = CafFooter::new(1024, 10, 100000, 0xDEADBEEF);
/// write_footer(&mut file, &footer)?;
/// # Ok(())
/// # }
/// ```
pub fn write_footer<W: Write>(writer: &mut W, footer: &CafFooter) -> Result<()> {
    bincode::serialize_into(writer, footer)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::CAF_FOOTER_MAGIC;
    use std::io::Cursor;

    #[test]
    fn test_write_read_footer_roundtrip() {
        let footer = CafFooter::new(1024, 10, 100000, 0xDEADBEEF);

        // Write
        let mut buffer = Vec::new();
        write_footer(&mut buffer, &footer).unwrap();

        // Read
        let mut reader = Cursor::new(buffer);
        let recovered = read_footer(&mut reader).unwrap();

        assert_eq!(recovered.index_offset, 1024);
        assert_eq!(recovered.num_blocks, 10);
        assert_eq!(recovered.total_records, 100000);
        assert_eq!(recovered.checksum, 0xDEADBEEF);
        assert!(recovered.validate_magic());
    }

    #[test]
    fn test_footer_magic_validation() {
        let footer = CafFooter::new(0, 0, 0, 0);
        assert_eq!(footer.magic, CAF_FOOTER_MAGIC);
        assert!(footer.validate_magic());
    }

    #[test]
    fn test_invalid_footer_magic() {
        let mut footer = CafFooter::new(0, 0, 0, 0);
        footer.magic = [b'B', b'A', b'M', b'!'];

        let mut buffer = Vec::new();
        write_footer(&mut buffer, &footer).unwrap();

        let mut reader = Cursor::new(buffer);
        let result = read_footer(&mut reader);

        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), CafError::InvalidMagic(_)));
    }
}
