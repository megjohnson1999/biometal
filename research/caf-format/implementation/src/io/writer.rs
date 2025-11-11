//! CAF file writer implementation.
//!
//! The CafFileWriter provides a high-level interface for creating CAF files:
//! 1. Write header with SAM metadata
//! 2. Add records (automatically accumulated into 10K blocks)
//! 3. Finalize (writes index + footer)
//!
//! # Example
//!
//! ```no_run
//! use caf::io::CafFileWriter;
//! use caf::block::AlignmentRecord;
//! use std::path::Path;
//!
//! # fn main() -> Result<(), caf::CafError> {
//! let mut writer = CafFileWriter::create(Path::new("output.caf"))?;
//!
//! // Set references
//! writer.set_references(vec!["chr1".to_string()], vec![248956422])?;
//!
//! // Add records
//! writer.add_record(AlignmentRecord {
//!     ref_id: 0,
//!     position: 1000,
//!     mapq: 60,
//!     sequence: b"ACGT".to_vec(),
//!     ..Default::default()
//! })?;
//!
//! // Finalize (writes index + footer)
//! writer.finalize()?;
//! # Ok(())
//! # }
//! ```

use crate::{
    block::{AlignmentRecord, BlockBuilder, DEFAULT_BLOCK_SIZE},
    format::{write_header, write_magic, write_index, write_footer},
    types::{CafHeader, CafIndex, CafFooter, BlockMeta},
    CafError, Result,
};
use std::fs::File;
use std::io::{BufWriter, Write, Seek, SeekFrom};
use std::path::Path;

/// CAF file writer with automatic block management.
pub struct CafFileWriter {
    /// Buffered file writer
    writer: BufWriter<File>,

    /// File header (updated on finalize)
    header: CafHeader,

    /// Current block builder
    block_builder: BlockBuilder,

    /// Block index (built as we write)
    index: CafIndex,

    /// Current file position (for tracking offsets)
    current_offset: u64,

    /// Block size (default 10,000)
    block_size: u32,

    /// Total records written
    total_records: u64,

    /// Whether file has been finalized
    finalized: bool,
}

impl CafFileWriter {
    /// Create a new CAF file for writing.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use caf::io::CafFileWriter;
    /// use std::path::Path;
    ///
    /// # fn main() -> Result<(), caf::CafError> {
    /// let writer = CafFileWriter::create(Path::new("output.caf"))?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Write magic number
        write_magic(&mut writer)?;

        // Create empty header and write it immediately
        let header = CafHeader::new(DEFAULT_BLOCK_SIZE, vec![]);

        // Serialize and write header
        let header_bytes = bincode::serialize(&header)?;
        writer.write_all(&header_bytes)?;

        // Track current offset (magic + header)
        let current_offset = 4 + header_bytes.len() as u64;

        Ok(Self {
            writer,
            header,
            block_builder: BlockBuilder::new(0, DEFAULT_BLOCK_SIZE),
            index: CafIndex::new(),
            current_offset,
            block_size: DEFAULT_BLOCK_SIZE,
            total_records: 0,
            finalized: false,
        })
    }

    /// Set reference sequences.
    ///
    /// Must be called before adding records.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use caf::io::CafFileWriter;
    /// use std::path::Path;
    ///
    /// # fn main() -> Result<(), caf::CafError> {
    /// let mut writer = CafFileWriter::create(Path::new("output.caf"))?;
    /// writer.set_references(
    ///     vec!["chr1".to_string(), "chr2".to_string()],
    ///     vec![248956422, 242193529]
    /// )?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn set_references(&mut self, names: Vec<String>, lengths: Vec<i32>) -> Result<()> {
        if names.len() != lengths.len() {
            return Err(CafError::Other(
                "Reference names and lengths must have same length".to_string(),
            ));
        }

        // Don't allow changing references after records have been added
        if self.total_records > 0 {
            return Err(CafError::Other(
                "Cannot set references after adding records. \
                 Call set_references() immediately after create()".to_string(),
            ));
        }

        // Calculate old header size
        let old_header_bytes = bincode::serialize(&self.header)?;
        let old_header_size = old_header_bytes.len() as u64;

        // Update header
        self.header.num_refs = names.len() as u32;
        self.header.ref_names = names;
        self.header.ref_lengths = lengths;

        // Serialize new header
        let new_header_bytes = bincode::serialize(&self.header)?;
        let new_header_size = new_header_bytes.len() as u64;

        // Adjust current_offset if header size changed
        let size_diff = (new_header_size as i64) - (old_header_size as i64);
        self.current_offset = (self.current_offset as i64 + size_diff) as u64;

        // Flush buffered writer before direct file access
        self.writer.flush()?;

        // Rewrite header at position 4 (using underlying file)
        let file = self.writer.get_mut();
        file.seek(SeekFrom::Start(4))?;
        file.write_all(&new_header_bytes)?;
        file.seek(SeekFrom::Start(self.current_offset))?; // Restore position
        file.flush()?;

        Ok(())
    }

    /// Set SAM header text.
    ///
    /// Must be called before adding records.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use caf::io::CafFileWriter;
    /// use std::path::Path;
    ///
    /// # fn main() -> Result<(), caf::CafError> {
    /// let mut writer = CafFileWriter::create(Path::new("output.caf"))?;
    /// writer.set_sam_header(b"@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:248956422\n".to_vec())?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn set_sam_header(&mut self, sam_header: Vec<u8>) -> Result<()> {
        // Don't allow changing SAM header after records have been added
        if self.total_records > 0 {
            return Err(CafError::Other(
                "Cannot set SAM header after adding records. \
                 Call set_sam_header() immediately after create()".to_string(),
            ));
        }

        // Calculate old header size
        let old_header_bytes = bincode::serialize(&self.header)?;
        let old_header_size = old_header_bytes.len() as u64;

        // Update header
        self.header.sam_header = sam_header;

        // Serialize new header
        let new_header_bytes = bincode::serialize(&self.header)?;
        let new_header_size = new_header_bytes.len() as u64;

        // Adjust current_offset if header size changed
        let size_diff = (new_header_size as i64) - (old_header_size as i64);
        self.current_offset = (self.current_offset as i64 + size_diff) as u64;

        // Flush buffered writer before direct file access
        self.writer.flush()?;

        // Rewrite header at position 4 (using underlying file)
        let file = self.writer.get_mut();
        file.seek(SeekFrom::Start(4))?;
        file.write_all(&new_header_bytes)?;
        file.seek(SeekFrom::Start(self.current_offset))?; // Restore position
        file.flush()?;

        Ok(())
    }

    /// Set the quality score dictionary for compression.
    ///
    /// When provided, quality scores will use zstd dictionary compression
    /// instead of raw storage, improving compression ratio by 1.5-3×.
    ///
    /// Must be called before adding any records.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use caf::io::CafFileWriter;
    /// use caf::compression::{train_dictionary, DICTIONARY_SIZE};
    /// use std::path::Path;
    ///
    /// # fn main() -> Result<(), caf::CafError> {
    /// let mut writer = CafFileWriter::create(Path::new("output.caf"))?;
    ///
    /// // Train dictionary from quality score samples
    /// let quality_samples = vec![/* ... */];
    /// let dict = train_dictionary(&quality_samples, DICTIONARY_SIZE)?;
    ///
    /// writer.set_quality_dict(dict)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn set_quality_dict(&mut self, dictionary: Vec<u8>) -> Result<()> {
        // Don't allow changing dictionary after records have been added
        if self.total_records > 0 {
            return Err(CafError::Other(
                "Cannot set quality dictionary after adding records. \
                 Call set_quality_dict() immediately after create()".to_string(),
            ));
        }

        // Calculate old header size
        let old_header_bytes = bincode::serialize(&self.header)?;
        let old_header_size = old_header_bytes.len() as u64;

        // Update header with dictionary
        self.header.quality_dict = Some(dictionary.clone());

        // Serialize new header
        let new_header_bytes = bincode::serialize(&self.header)?;
        let new_header_size = new_header_bytes.len() as u64;

        // Adjust current_offset if header size changed
        let size_diff = (new_header_size as i64) - (old_header_size as i64);
        self.current_offset = (self.current_offset as i64 + size_diff) as u64;

        // Set dictionary on current block builder
        self.block_builder.set_quality_dict(dictionary);

        // Flush buffered writer before direct file access
        self.writer.flush()?;

        // Rewrite header at position 4 (using underlying file)
        let file = self.writer.get_mut();
        file.seek(SeekFrom::Start(4))?;
        file.write_all(&new_header_bytes)?;
        file.seek(SeekFrom::Start(self.current_offset))?; // Restore position
        file.flush()?;

        Ok(())
    }

    /// Add a record to the file.
    ///
    /// Records are automatically accumulated into blocks. When a block reaches
    /// the configured size (default 10,000), it is automatically written to disk.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use caf::io::CafFileWriter;
    /// use caf::block::AlignmentRecord;
    /// use std::path::Path;
    ///
    /// # fn main() -> Result<(), caf::CafError> {
    /// let mut writer = CafFileWriter::create(Path::new("output.caf"))?;
    ///
    /// writer.add_record(AlignmentRecord {
    ///     ref_id: 0,
    ///     position: 1000,
    ///     mapq: 60,
    ///     sequence: b"ACGT".to_vec(),
    ///     ..Default::default()
    /// })?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn add_record(&mut self, record: AlignmentRecord) -> Result<()> {
        if self.finalized {
            return Err(CafError::Other(
                "Cannot add records after finalization".to_string(),
            ));
        }

        // Add to current block
        self.block_builder.add_record(record)?;
        self.total_records += 1;

        // If block is full, write it
        if self.block_builder.is_full() {
            self.flush_block()?;
        }

        Ok(())
    }

    /// Flush the current block to disk.
    ///
    /// This is called automatically when a block reaches capacity,
    /// but can be called manually to force a block write.
    fn flush_block(&mut self) -> Result<()> {
        if self.block_builder.is_empty() {
            return Ok(());
        }

        // Build the block
        let block = self.block_builder.build()?;

        // Track block offset
        let block_offset = self.current_offset;

        // Create block metadata for index
        let meta = BlockMeta {
            num_records: block.num_records,
            ref_id: -1, // TODO: Track actual ref_id range
            start_pos: 0, // TODO: Track actual position range
            end_pos: 0,
            compressed_size: 0, // Will be calculated from serialized size
            uncompressed_size: block.uncompressed_size,
            checksum: block.checksum,
        };

        // Serialize and write block
        let block_bytes = bincode::serialize(&block)?;
        self.writer.write_all(&block_bytes)?;
        self.current_offset += block_bytes.len() as u64;

        // Update metadata with actual compressed size
        let meta = BlockMeta {
            compressed_size: block_bytes.len() as u32,
            ..meta
        };

        // Add to index
        self.index.add_block(block_offset, meta);

        // Reset builder for next block
        self.block_builder = BlockBuilder::new(self.index.num_blocks, self.block_size);

        // Set dictionary on new builder if we have one
        if let Some(ref dict) = self.header.quality_dict {
            self.block_builder.set_quality_dict(dict.clone());
        }

        Ok(())
    }

    /// Finalize the file by writing index and footer.
    ///
    /// This MUST be called to produce a valid CAF file. It:
    /// 1. Flushes any remaining records
    /// 2. Writes the block index
    /// 3. Writes the footer
    /// 4. Updates the header with final block count
    ///
    /// # Example
    ///
    /// ```no_run
    /// use caf::io::CafFileWriter;
    /// use caf::block::AlignmentRecord;
    /// use std::path::Path;
    ///
    /// # fn main() -> Result<(), caf::CafError> {
    /// let mut writer = CafFileWriter::create(Path::new("output.caf"))?;
    /// writer.add_record(AlignmentRecord::default())?;
    /// writer.finalize()?; // Required!
    /// # Ok(())
    /// # }
    /// ```
    pub fn finalize(&mut self) -> Result<()> {
        if self.finalized {
            return Err(CafError::Other("File already finalized".to_string()));
        }

        // Flush any remaining records
        if !self.block_builder.is_empty() {
            self.flush_block()?;
        }

        // Record index offset
        let index_offset = self.current_offset;
        self.index.index_offset = index_offset;

        // Write index
        write_index(&mut self.writer, &self.index)?;
        let index_bytes = bincode::serialize(&self.index)?;
        self.current_offset += index_bytes.len() as u64;

        // Write footer
        let footer = CafFooter::new(
            index_offset,
            self.index.num_blocks,
            self.total_records,
            0, // TODO: Calculate file checksum
        );
        write_footer(&mut self.writer, &footer)?;

        // Flush writer
        self.writer.flush()?;

        // Update header with final block count and rewrite
        self.header.num_blocks = self.index.num_blocks;
        let header_bytes = bincode::serialize(&self.header)?;

        let file = self.writer.get_mut();
        file.seek(SeekFrom::Start(4))?; // After magic
        file.write_all(&header_bytes)?;
        file.flush()?;

        self.finalized = true;

        Ok(())
    }

    /// Get total number of records written.
    pub fn total_records(&self) -> u64 {
        self.total_records
    }

    /// Get number of blocks written.
    pub fn num_blocks(&self) -> u32 {
        self.index.num_blocks
    }
}

impl Drop for CafFileWriter {
    fn drop(&mut self) {
        if !self.finalized {
            eprintln!("Warning: CafFileWriter dropped without calling finalize()");
            eprintln!("The output file may be incomplete or invalid.");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
    fn test_writer_create() {
        let temp = NamedTempFile::new().unwrap();
        let writer = CafFileWriter::create(temp.path()).unwrap();
        assert_eq!(writer.total_records(), 0);
        assert_eq!(writer.num_blocks(), 0);
    }

    #[test]
    fn test_writer_add_single_record() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = CafFileWriter::create(temp.path()).unwrap();

        writer.add_record(create_test_record(0, 1000)).unwrap();
        assert_eq!(writer.total_records(), 1);

        writer.finalize().unwrap();
        assert_eq!(writer.num_blocks(), 1);
    }

    #[test]
    fn test_writer_multiple_records() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = CafFileWriter::create(temp.path()).unwrap();

        for i in 0..100 {
            writer.add_record(create_test_record(0, 1000 + i)).unwrap();
        }

        assert_eq!(writer.total_records(), 100);
        writer.finalize().unwrap();
        assert_eq!(writer.num_blocks(), 1);
    }

    #[test]
    fn test_writer_multiple_blocks() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = CafFileWriter::create(temp.path()).unwrap();

        // Add enough records for 2+ blocks (10K each)
        for i in 0..25_000 {
            writer.add_record(create_test_record(0, 1000 + (i % 10000))).unwrap();
        }

        assert_eq!(writer.total_records(), 25_000);
        writer.finalize().unwrap();
        assert_eq!(writer.num_blocks(), 3); // 10K + 10K + 5K
    }

    #[test]
    fn test_writer_set_references() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = CafFileWriter::create(temp.path()).unwrap();

        writer.set_references(
            vec!["chr1".to_string(), "chr2".to_string()],
            vec![248956422, 242193529],
        ).unwrap();

        writer.finalize().unwrap();
    }

    #[test]
    fn test_writer_finalize_required() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = CafFileWriter::create(temp.path()).unwrap();

        writer.add_record(create_test_record(0, 1000)).unwrap();

        // File should be incomplete without finalize
        drop(writer);
        // (Warning printed to stderr)
    }

    #[test]
    fn test_writer_cannot_add_after_finalize() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = CafFileWriter::create(temp.path()).unwrap();

        writer.add_record(create_test_record(0, 1000)).unwrap();
        writer.finalize().unwrap();

        let result = writer.add_record(create_test_record(0, 2000));
        assert!(result.is_err());
    }

    #[test]
    fn test_write_and_read_back() {
        use crate::io::CafFileReader;

        let temp = NamedTempFile::new().unwrap();

        // Write
        let mut writer = CafFileWriter::create(temp.path()).unwrap();
        writer.add_record(create_test_record(0, 1000)).unwrap();

        println!("Before finalize - blocks: {}, records: {}",
                 writer.num_blocks(), writer.total_records());

        writer.finalize().unwrap();

        println!("After finalize - blocks: {}, records: {}",
                 writer.num_blocks(), writer.total_records());

        // Read back
        let reader = CafFileReader::open(temp.path()).unwrap();
        println!("Reader - blocks: {}, records: {}",
                 reader.num_blocks(), reader.total_records());

        assert_eq!(reader.num_blocks(), 1);
        assert_eq!(reader.total_records(), 1);
    }
}
