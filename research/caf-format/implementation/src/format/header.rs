//! CAF header parsing and writing.
//!
//! The CAF header contains:
//! - Version information (u16)
//! - Block configuration (block_size, num_blocks)
//! - Compressed SAM header
//! - Reference sequences (names + lengths)
//! - Column schema

use crate::{types::CafHeader, CafError, Result};
use std::io::{Read, Write};

/// Read CAF header from a reader.
///
/// The header is serialized using bincode and contains all metadata
/// needed to interpret the columnar blocks.
///
/// # Errors
///
/// - `CafError::Serialization` if header deserialization fails
/// - `CafError::UnsupportedVersion` if version is not supported
/// - `CafError::Io` for I/O errors
///
/// # Example
///
/// ```rust,no_run
/// use std::fs::File;
/// # use caf::format::read_header;
/// # fn main() -> Result<(), caf::error::CafError> {
/// let mut file = File::open("data.caf")?;
/// // Skip magic number first
/// let header = read_header(&mut file)?;
/// println!("CAF version: {}.{}", header.version_major(), header.version_minor());
/// # Ok(())
/// # }
/// ```
pub fn read_header<R: Read>(reader: &mut R) -> Result<CafHeader> {
    // Deserialize header using bincode
    let header: CafHeader = bincode::deserialize_from(reader)?;

    // Validate version
    let major = header.version_major();
    let minor = header.version_minor();

    if major != 1 {
        return Err(CafError::UnsupportedVersion { major, minor });
    }

    // Validate block size
    if header.block_size == 0 {
        return Err(CafError::Other("Header block_size must be > 0".to_string()));
    }

    // Validate reference sequence consistency
    if header.ref_names.len() != header.num_refs as usize {
        return Err(CafError::Other(format!(
            "Header reference inconsistency: num_refs={} but {} names provided",
            header.num_refs,
            header.ref_names.len()
        )));
    }

    if header.ref_lengths.len() != header.num_refs as usize {
        return Err(CafError::Other(format!(
            "Header reference inconsistency: num_refs={} but {} lengths provided",
            header.num_refs,
            header.ref_lengths.len()
        )));
    }

    Ok(header)
}

/// Write CAF header to a writer.
///
/// # Errors
///
/// - `CafError::Serialization` if header serialization fails
/// - `CafError::Io` for I/O errors
///
/// # Example
///
/// ```rust,no_run
/// use std::fs::File;
/// use caf::types::CafHeader;
/// # use caf::format::write_header;
/// # fn main() -> Result<(), caf::error::CafError> {
/// let mut file = File::create("data.caf")?;
/// let header = CafHeader::new(10_000, vec![]);
/// write_header(&mut file, &header)?;
/// # Ok(())
/// # }
/// ```
pub fn write_header<W: Write>(writer: &mut W, header: &CafHeader) -> Result<()> {
    // Serialize header using bincode
    bincode::serialize_into(writer, header)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn create_test_header() -> CafHeader {
        let mut header = CafHeader::new(10_000, b"@HD\tVN:1.0\n".to_vec());
        header.num_refs = 2;
        header.ref_names = vec!["chr1".to_string(), "chr2".to_string()];
        header.ref_lengths = vec![248956422, 242193529];
        header
    }

    #[test]
    fn test_write_read_header_roundtrip() {
        let original = create_test_header();

        // Write to buffer
        let mut buffer = Vec::new();
        write_header(&mut buffer, &original).unwrap();

        // Read back
        let mut reader = Cursor::new(buffer);
        let recovered = read_header(&mut reader).unwrap();

        // Verify
        assert_eq!(recovered.version, original.version);
        assert_eq!(recovered.block_size, original.block_size);
        assert_eq!(recovered.num_refs, original.num_refs);
        assert_eq!(recovered.ref_names, original.ref_names);
        assert_eq!(recovered.ref_lengths, original.ref_lengths);
    }

    #[test]
    fn test_header_version_validation() {
        let header = create_test_header();
        assert_eq!(header.version_major(), 1);
        assert_eq!(header.version_minor(), 0);
    }

    #[test]
    fn test_unsupported_version() {
        // Create header with version 2.0
        let mut header = create_test_header();
        header.version = 0x0200; // Version 2.0

        let mut buffer = Vec::new();
        write_header(&mut buffer, &header).unwrap();

        let mut reader = Cursor::new(buffer);
        let result = read_header(&mut reader);

        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            CafError::UnsupportedVersion { major: 2, minor: 0 }
        ));
    }

    #[test]
    fn test_empty_sam_header() {
        let header = CafHeader::new(10_000, vec![]);

        let mut buffer = Vec::new();
        write_header(&mut buffer, &header).unwrap();

        let mut reader = Cursor::new(buffer);
        let recovered = read_header(&mut reader).unwrap();

        assert_eq!(recovered.sam_header, Vec::<u8>::new());
    }

    #[test]
    fn test_header_with_references() {
        let mut header = CafHeader::new(10_000, vec![]);
        header.num_refs = 3;
        header.ref_names = vec![
            "chr1".to_string(),
            "chr2".to_string(),
            "chrM".to_string(),
        ];
        header.ref_lengths = vec![248956422, 242193529, 16569];

        let mut buffer = Vec::new();
        write_header(&mut buffer, &header).unwrap();

        let mut reader = Cursor::new(buffer);
        let recovered = read_header(&mut reader).unwrap();

        assert_eq!(recovered.num_refs, 3);
        assert_eq!(recovered.ref_names.len(), 3);
        assert_eq!(recovered.ref_lengths.len(), 3);
        assert_eq!(recovered.ref_names[2], "chrM");
        assert_eq!(recovered.ref_lengths[2], 16569);
    }
}
