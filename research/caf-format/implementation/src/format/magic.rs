//! Magic number validation for CAF files.
//!
//! CAF files begin with a 4-byte magic number: "CAF\x01"
//! - Bytes 0-2: "CAF" (format identifier)
//! - Byte 3: Major version number (0x01 for v1.x)

use crate::{error::CafError, Result, CAF_MAGIC};
use std::io::{Read, Write};

/// Read and validate CAF magic number from a reader.
///
/// # Errors
///
/// Returns `CafError::InvalidMagic` if the magic number doesn't match "CAF\x01".
///
/// # Example
///
/// ```rust,no_run
/// use std::io::Cursor;
/// # use caf::format::read_magic;
/// # fn main() -> Result<(), caf::error::CafError> {
/// let data = b"CAF\x01...";
/// let mut reader = Cursor::new(data);
/// let version = read_magic(&mut reader)?;
/// assert_eq!(version, 1);
/// # Ok(())
/// # }
/// ```
pub fn read_magic<R: Read>(reader: &mut R) -> Result<u8> {
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;

    // Validate magic bytes
    if magic[0..3] != CAF_MAGIC[0..3] {
        return Err(CafError::InvalidMagic(magic));
    }

    // Extract and return major version
    let major_version = magic[3];
    Ok(major_version)
}

/// Write CAF magic number to a writer.
///
/// # Example
///
/// ```rust,no_run
/// use std::io::Cursor;
/// # use caf::format::write_magic;
/// # fn main() -> Result<(), caf::error::CafError> {
/// let mut buffer = Vec::new();
/// write_magic(&mut buffer)?;
/// assert_eq!(&buffer[..], b"CAF\x01");
/// # Ok(())
/// # }
/// ```
pub fn write_magic<W: Write>(writer: &mut W) -> Result<()> {
    writer.write_all(&CAF_MAGIC)?;
    Ok(())
}

/// Validate that a byte array contains a valid CAF magic number.
///
/// Returns the major version number if valid.
pub fn validate_magic(magic: &[u8; 4]) -> Result<u8> {
    if magic[0..3] != CAF_MAGIC[0..3] {
        return Err(CafError::InvalidMagic(*magic));
    }
    Ok(magic[3])
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_read_magic_valid() {
        let data = b"CAF\x01rest of file";
        let mut reader = Cursor::new(data);
        let version = read_magic(&mut reader).unwrap();
        assert_eq!(version, 1);
        assert_eq!(reader.position(), 4); // Should have read exactly 4 bytes
    }

    #[test]
    fn test_read_magic_invalid() {
        let data = b"BAM\x01rest of file";
        let mut reader = Cursor::new(data);
        let result = read_magic(&mut reader);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), CafError::InvalidMagic(_)));
    }

    #[test]
    fn test_read_magic_wrong_prefix() {
        let data = b"CA\x00\x01rest";
        let mut reader = Cursor::new(data);
        let result = read_magic(&mut reader);
        assert!(result.is_err());
    }

    #[test]
    fn test_read_magic_too_short() {
        let data = b"CAF";
        let mut reader = Cursor::new(data);
        let result = read_magic(&mut reader);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), CafError::Io(_)));
    }

    #[test]
    fn test_write_magic() {
        let mut buffer = Vec::new();
        write_magic(&mut buffer).unwrap();
        assert_eq!(buffer, vec![b'C', b'A', b'F', 0x01]);
    }

    #[test]
    fn test_write_read_roundtrip() {
        let mut buffer = Vec::new();
        write_magic(&mut buffer).unwrap();

        let mut reader = Cursor::new(buffer);
        let version = read_magic(&mut reader).unwrap();
        assert_eq!(version, 1);
    }

    #[test]
    fn test_validate_magic_valid() {
        let magic = [b'C', b'A', b'F', 0x01];
        let version = validate_magic(&magic).unwrap();
        assert_eq!(version, 1);
    }

    #[test]
    fn test_validate_magic_invalid() {
        let magic = [b'B', b'A', b'M', 0x01];
        let result = validate_magic(&magic);
        assert!(result.is_err());
    }

    #[test]
    fn test_validate_magic_future_version() {
        let magic = [b'C', b'A', b'F', 0x02];
        let version = validate_magic(&magic).unwrap();
        assert_eq!(version, 2); // Future version should validate but return version
    }
}
