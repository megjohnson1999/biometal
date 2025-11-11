//! CAF file reader implementation.
//!
//! The CafFileReader provides a high-level interface for reading CAF files:
//! 1. Read header and index
//! 2. Random access to blocks
//! 3. Region queries
//! 4. Streaming record iteration
//!
//! # Example
//!
//! ```no_run
//! use caf::io::CafFileReader;
//! use std::path::Path;
//!
//! # fn main() -> Result<(), caf::CafError> {
//! let reader = CafFileReader::open(Path::new("input.caf"))?;
//!
//! // Iterate over all records
//! for record in reader.records()? {
//!     let record = record?;
//!     println!("Position: {}", record.position);
//! }
//! # Ok(())
//! # }
//! ```

use crate::{
    block::{AlignmentRecord, BlockReader},
    format::{read_header, read_magic, read_index, read_footer},
    types::{CafHeader, CafIndex, CafFooter, CafBlock},
    CafError, Result,
};
use std::fs::File;
use std::io::{BufReader, Seek, SeekFrom};
use std::path::Path;

/// CAF file reader with index-based random access.
pub struct CafFileReader {
    /// Buffered file reader
    reader: BufReader<File>,

    /// File header
    header: CafHeader,

    /// Block index for random access
    index: CafIndex,

    /// Footer metadata
    footer: CafFooter,
}

impl CafFileReader {
    /// Open a CAF file for reading.
    ///
    /// This reads the magic number, header, footer, and index to prepare for
    /// random access and streaming.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use caf::io::CafFileReader;
    /// use std::path::Path;
    ///
    /// # fn main() -> Result<(), caf::CafError> {
    /// let reader = CafFileReader::open(Path::new("input.caf"))?;
    /// println!("Total blocks: {}", reader.num_blocks());
    /// # Ok(())
    /// # }
    /// ```
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        // Read and validate magic number
        read_magic(&mut reader)?;

        // Read header
        let header = read_header(&mut reader)?;

        // Read footer (at end of file - 32 bytes)
        reader.seek(SeekFrom::End(-32))?;
        let footer = read_footer(&mut reader)?;

        // Read index (using offset from footer)
        reader.seek(SeekFrom::Start(footer.index_offset))?;
        let index = read_index(&mut reader)?;

        // Validate consistency
        if index.num_blocks != footer.num_blocks {
            return Err(CafError::Other(format!(
                "Index/footer mismatch: index has {} blocks, footer says {}",
                index.num_blocks, footer.num_blocks
            )));
        }

        Ok(Self {
            reader,
            header,
            index,
            footer,
        })
    }

    /// Get the file header.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use caf::io::CafFileReader;
    /// use std::path::Path;
    ///
    /// # fn main() -> Result<(), caf::CafError> {
    /// let reader = CafFileReader::open(Path::new("input.caf"))?;
    /// let header = reader.header();
    /// println!("Block size: {}", header.block_size);
    /// # Ok(())
    /// # }
    /// ```
    pub fn header(&self) -> &CafHeader {
        &self.header
    }

    /// Get the block index.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use caf::io::CafFileReader;
    /// use std::path::Path;
    ///
    /// # fn main() -> Result<(), caf::CafError> {
    /// let reader = CafFileReader::open(Path::new("input.caf"))?;
    /// let index = reader.index();
    /// println!("Total blocks: {}", index.num_blocks);
    /// # Ok(())
    /// # }
    /// ```
    pub fn index(&self) -> &CafIndex {
        &self.index
    }

    /// Get the total number of blocks.
    pub fn num_blocks(&self) -> u32 {
        self.index.num_blocks
    }

    /// Get the total number of records (from footer).
    pub fn total_records(&self) -> u64 {
        self.footer.total_records
    }

    /// Read a specific block by block ID.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use caf::io::CafFileReader;
    /// use std::path::Path;
    ///
    /// # fn main() -> Result<(), caf::CafError> {
    /// let mut reader = CafFileReader::open(Path::new("input.caf"))?;
    /// let block_reader = reader.read_block(0)?;
    ///
    /// // Access records in this block
    /// let record = block_reader.get_record(0)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn read_block(&mut self, block_id: u32) -> Result<BlockReader> {
        if block_id >= self.index.num_blocks {
            return Err(CafError::Other(format!(
                "Block ID {} out of bounds (max {})",
                block_id,
                self.index.num_blocks - 1
            )));
        }

        // Get block offset from index
        let offset = self.index.block_offsets[block_id as usize];

        // Seek to block position
        self.reader.seek(SeekFrom::Start(offset))?;

        // Deserialize block
        let block: CafBlock = bincode::deserialize_from(&mut self.reader)?;

        // Create reader with dictionary from header (validates checksum)
        let quality_dict = self.header.quality_dict.as_deref();
        BlockReader::with_dictionary(block, quality_dict)
    }

    /// Iterate over all records in the file.
    ///
    /// This provides a streaming interface that reads blocks on demand.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use caf::io::CafFileReader;
    /// use std::path::Path;
    ///
    /// # fn main() -> Result<(), caf::CafError> {
    /// let reader = CafFileReader::open(Path::new("input.caf"))?;
    ///
    /// for record in reader.records()? {
    ///     let record = record?;
    ///     println!("Position: {}", record.position);
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn records(self) -> Result<CafRecordIterator> {
        Ok(CafRecordIterator {
            reader: self,
            current_block: 0,
            current_block_reader: None,
            current_record_in_block: 0,
        })
    }

    /// Query records in a specific genomic region.
    ///
    /// Returns an iterator over records overlapping the region [start, end).
    ///
    /// # Example
    ///
    /// ```no_run
    /// use caf::io::CafFileReader;
    /// use std::path::Path;
    ///
    /// # fn main() -> Result<(), caf::CafError> {
    /// let reader = CafFileReader::open(Path::new("input.caf"))?;
    ///
    /// // Query chr1:1000-2000
    /// for record in reader.query_region(0, 1000, 2000)? {
    ///     let record = record?;
    ///     println!("Position: {}", record.position);
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn query_region(
        self,
        ref_id: i32,
        start: i32,
        end: i32,
    ) -> Result<CafRegionIterator> {
        // Find blocks overlapping the region
        let overlapping_blocks = self.find_overlapping_blocks(ref_id, start, end);

        Ok(CafRegionIterator {
            reader: self,
            ref_id,
            start,
            end,
            block_ids: overlapping_blocks,
            current_block_index: 0,
            current_block_reader: None,
            current_record_in_block: 0,
        })
    }

    /// Find blocks that may overlap a genomic region.
    ///
    /// This is a simple linear scan through block metadata. A production
    /// implementation would use a spatial index (R-tree, interval tree, etc.).
    fn find_overlapping_blocks(&self, ref_id: i32, start: i32, end: i32) -> Vec<u32> {
        let mut overlapping = Vec::new();

        for block_id in 0..self.index.num_blocks {
            let meta = &self.index.block_metadata[block_id as usize];

            // Check if block overlaps region
            // (In current implementation, ref_id/positions aren't tracked in metadata,
            //  so we conservatively include all blocks. This will be optimized later.)
            if meta.ref_id == -1 || meta.ref_id == ref_id {
                // Conservative: include if we don't have position info
                if meta.start_pos == 0 && meta.end_pos == 0 {
                    overlapping.push(block_id);
                } else if meta.start_pos < end && meta.end_pos >= start {
                    overlapping.push(block_id);
                }
            }
        }

        overlapping
    }
}

/// Iterator over all records in a CAF file.
pub struct CafRecordIterator {
    reader: CafFileReader,
    current_block: u32,
    current_block_reader: Option<BlockReader>,
    current_record_in_block: usize,
}

impl Iterator for CafRecordIterator {
    type Item = Result<AlignmentRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // If we don't have a block reader, load the next block
            if self.current_block_reader.is_none() {
                if self.current_block >= self.reader.num_blocks() {
                    return None; // No more blocks
                }

                // Load block
                match self.reader.read_block(self.current_block) {
                    Ok(block_reader) => {
                        self.current_block_reader = Some(block_reader);
                        self.current_record_in_block = 0;
                        self.current_block += 1;
                    }
                    Err(e) => return Some(Err(e)),
                }
            }

            // Try to get next record from current block
            if let Some(ref block_reader) = self.current_block_reader {
                if self.current_record_in_block < block_reader.num_records() {
                    let record = block_reader.get_record(self.current_record_in_block);
                    self.current_record_in_block += 1;
                    return Some(record);
                } else {
                    // Exhausted current block, move to next
                    self.current_block_reader = None;
                }
            }
        }
    }
}

/// Iterator over records in a genomic region.
pub struct CafRegionIterator {
    reader: CafFileReader,
    ref_id: i32,
    start: i32,
    end: i32,
    block_ids: Vec<u32>,
    current_block_index: usize,
    current_block_reader: Option<BlockReader>,
    current_record_in_block: usize,
}

impl Iterator for CafRegionIterator {
    type Item = Result<AlignmentRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // If we don't have a block reader, load the next block
            if self.current_block_reader.is_none() {
                if self.current_block_index >= self.block_ids.len() {
                    return None; // No more blocks in region
                }

                // Load block
                let block_id = self.block_ids[self.current_block_index];
                match self.reader.read_block(block_id) {
                    Ok(block_reader) => {
                        self.current_block_reader = Some(block_reader);
                        self.current_record_in_block = 0;
                        self.current_block_index += 1;
                    }
                    Err(e) => return Some(Err(e)),
                }
            }

            // Try to get next record from current block
            if let Some(ref block_reader) = self.current_block_reader {
                if self.current_record_in_block < block_reader.num_records() {
                    match block_reader.get_record(self.current_record_in_block) {
                        Ok(record) => {
                            self.current_record_in_block += 1;

                            // Filter by region
                            if record.ref_id == self.ref_id
                                && record.position >= self.start
                                && record.position < self.end
                            {
                                return Some(Ok(record));
                            }
                            // Otherwise continue to next record
                        }
                        Err(e) => return Some(Err(e)),
                    }
                } else {
                    // Exhausted current block, move to next
                    self.current_block_reader = None;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::CafFileWriter;
    use crate::block::AlignmentRecord;
    use tempfile::NamedTempFile;

    fn create_test_record(ref_id: i32, position: i32) -> AlignmentRecord {
        AlignmentRecord {
            ref_id,
            position,
            mapq: 60,
            flags: 99,
            sequence: b"ACGT".to_vec(),
            qualities: b"IIII".to_vec(),
            cigar: vec![(0, 4)],
            read_name: b"read1".to_vec(),
            mate_ref_id: ref_id,
            mate_position: position + 100,
            template_length: 200,
        }
    }

    #[test]
    fn test_reader_open() {
        // Create a test file
        let temp = NamedTempFile::new().unwrap();
        let mut writer = CafFileWriter::create(temp.path()).unwrap();
        writer.add_record(create_test_record(0, 1000)).unwrap();
        writer.finalize().unwrap();

        // Open for reading
        let reader = CafFileReader::open(temp.path()).unwrap();
        assert_eq!(reader.num_blocks(), 1);
        assert_eq!(reader.total_records(), 1);
    }

    #[test]
    fn test_reader_read_block() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = CafFileWriter::create(temp.path()).unwrap();

        for i in 0..100 {
            writer.add_record(create_test_record(0, 1000 + i)).unwrap();
        }
        writer.finalize().unwrap();

        let mut reader = CafFileReader::open(temp.path()).unwrap();
        let block_reader = reader.read_block(0).unwrap();

        assert_eq!(block_reader.num_records(), 100);

        let record = block_reader.get_record(0).unwrap();
        assert_eq!(record.position, 1000);
    }

    #[test]
    fn test_reader_records_iterator() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = CafFileWriter::create(temp.path()).unwrap();

        for i in 0..250 {
            writer.add_record(create_test_record(0, 1000 + i)).unwrap();
        }
        writer.finalize().unwrap();

        let reader = CafFileReader::open(temp.path()).unwrap();
        let records: Vec<_> = reader.records().unwrap().collect();

        assert_eq!(records.len(), 250);

        // Verify first and last
        assert_eq!(records[0].as_ref().unwrap().position, 1000);
        assert_eq!(records[249].as_ref().unwrap().position, 1249);
    }

    #[test]
    fn test_reader_multiple_blocks() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = CafFileWriter::create(temp.path()).unwrap();

        // Create 25,000 records (2.5 blocks)
        for i in 0..25_000 {
            writer.add_record(create_test_record(0, 1000 + (i % 10000))).unwrap();
        }
        writer.finalize().unwrap();

        let reader = CafFileReader::open(temp.path()).unwrap();
        assert_eq!(reader.num_blocks(), 3);
        assert_eq!(reader.total_records(), 25_000);

        let records: Vec<_> = reader.records().unwrap().collect();
        assert_eq!(records.len(), 25_000);
    }

    #[test]
    fn test_reader_region_query() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = CafFileWriter::create(temp.path()).unwrap();

        // Create records at different positions
        for i in 0..100 {
            writer.add_record(create_test_record(0, 1000 + i * 10)).unwrap();
        }
        writer.finalize().unwrap();

        let reader = CafFileReader::open(temp.path()).unwrap();

        // Query region 1000-1500 (should get ~50 records)
        let records: Vec<_> = reader
            .query_region(0, 1000, 1500)
            .unwrap()
            .collect();

        // Verify all returned records are in range
        for record in &records {
            let record = record.as_ref().unwrap();
            assert!(record.position >= 1000);
            assert!(record.position < 1500);
        }

        assert!(records.len() > 0);
        assert!(records.len() <= 50);
    }

    #[test]
    fn test_reader_invalid_block_id() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = CafFileWriter::create(temp.path()).unwrap();
        writer.add_record(create_test_record(0, 1000)).unwrap();
        writer.finalize().unwrap();

        let mut reader = CafFileReader::open(temp.path()).unwrap();
        let result = reader.read_block(999);

        assert!(result.is_err());
    }
}
