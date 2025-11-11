//! CAF index parsing and writing.
//!
//! The index provides block-level random access to CAF files.
//! It contains file offsets and metadata for each block.

use crate::{types::CafIndex, Result};
use std::io::{Read, Write};

/// Read CAF index from a reader.
///
/// The index is stored at the end of the file and contains
/// offsets and metadata for all blocks.
///
/// # Example
///
/// ```rust,no_run
/// use std::fs::File;
/// use std::io::Seek;
/// # use caf::format::read_index;
/// # fn main() -> Result<(), caf::error::CafError> {
/// let mut file = File::open("data.caf")?;
/// // Seek to index offset (from footer)
/// let index = read_index(&mut file)?;
/// println!("Blocks: {}", index.num_blocks);
/// # Ok(())
/// # }
/// ```
pub fn read_index<R: Read>(reader: &mut R) -> Result<CafIndex> {
    let index: CafIndex = bincode::deserialize_from(reader)?;
    Ok(index)
}

/// Write CAF index to a writer.
///
/// # Example
///
/// ```rust,no_run
/// use std::fs::File;
/// use caf::types::CafIndex;
/// # use caf::format::write_index;
/// # fn main() -> Result<(), caf::error::CafError> {
/// let mut file = File::create("data.caf")?;
/// let index = CafIndex::new();
/// write_index(&mut file, &index)?;
/// # Ok(())
/// # }
/// ```
pub fn write_index<W: Write>(writer: &mut W, index: &CafIndex) -> Result<()> {
    bincode::serialize_into(writer, index)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::BlockMeta;
    use std::io::Cursor;

    #[test]
    fn test_write_read_index_roundtrip() {
        let mut index = CafIndex::new();
        index.add_block(
            1024,
            BlockMeta {
                num_records: 10000,
                ref_id: 0,
                start_pos: 1000,
                end_pos: 2000,
                compressed_size: 50000,
                uncompressed_size: 150000,
                checksum: 0xDEADBEEF,
            },
        );

        // Write
        let mut buffer = Vec::new();
        write_index(&mut buffer, &index).unwrap();

        // Read
        let mut reader = Cursor::new(buffer);
        let recovered = read_index(&mut reader).unwrap();

        assert_eq!(recovered.num_blocks, 1);
        assert_eq!(recovered.block_offsets[0], 1024);
        assert_eq!(recovered.block_metadata[0].num_records, 10000);
    }

    #[test]
    fn test_empty_index() {
        let index = CafIndex::new();

        let mut buffer = Vec::new();
        write_index(&mut buffer, &index).unwrap();

        let mut reader = Cursor::new(buffer);
        let recovered = read_index(&mut reader).unwrap();

        assert_eq!(recovered.num_blocks, 0);
        assert_eq!(recovered.block_offsets.len(), 0);
        assert_eq!(recovered.block_metadata.len(), 0);
    }
}
