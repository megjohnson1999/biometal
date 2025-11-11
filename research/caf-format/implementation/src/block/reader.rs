//! Block reader for decoding columnar blocks to records.
//!
//! The BlockReader takes a CafBlock in columnar format and reconstructs
//! the original alignment records by:
//! 1. Decompressing each column
//! 2. Decoding column-specific encodings (delta, zigzag)
//! 3. Reconstructing records from columns using offset arrays
//! 4. Validating checksums
//!
//! # Example
//!
//! ```
//! use caf::block::{BlockBuilder, BlockReader, AlignmentRecord};
//!
//! // Build a block
//! let mut builder = BlockBuilder::new(0, 10_000);
//! builder.add_record(AlignmentRecord {
//!     ref_id: 0,
//!     position: 1000,
//!     mapq: 60,
//!     sequence: b"ACGT".to_vec(),
//!     ..Default::default()
//! })?;
//! let block = builder.build()?;
//!
//! // Read it back
//! let reader = BlockReader::new(block)?;
//! let records: Vec<_> = reader.into_iter().collect::<Result<Vec<_>, _>>()?;
//! assert_eq!(records.len(), 1);
//! assert_eq!(records[0].position, 1000);
//! # Ok::<(), caf::CafError>(())
//! ```

use crate::{
    column::{decode_integers_delta, decode_sequence_ascii, decode_qualities_raw, decode_cigar_op},
    compression::decompress,
    types::{CafBlock, ColumnData, CompressedColumn},
    CafError, Result,
};
use super::builder::AlignmentRecord;

/// Decompressed columns ready for record reconstruction.
#[derive(Debug)]
struct DecompressedColumns {
    ref_ids: Vec<i32>,
    positions: Vec<i32>,
    mapq: Vec<u8>,
    flags: Vec<u16>,
    sequences: Vec<u8>,
    seq_offsets: Vec<u32>,
    qualities: Vec<u8>,
    qual_offsets: Vec<u32>,
    cigar_ops: Vec<u32>,
    cigar_offsets: Vec<u32>,
    read_names: Vec<u8>,
    read_name_offsets: Vec<u32>,
    mate_ref_ids: Vec<i32>,
    mate_positions: Vec<i32>,
    template_lengths: Vec<i32>,
}

/// Block reader for decoding columnar blocks.
#[derive(Debug)]
pub struct BlockReader {
    /// Number of records in block
    num_records: usize,

    /// Decompressed columns
    columns: DecompressedColumns,
}

impl BlockReader {
    /// Create a new block reader.
    ///
    /// This decompresses and decodes all columns upfront for efficient
    /// random access to records.
    ///
    /// # Errors
    ///
    /// Returns error if:
    /// - Column decompression fails
    /// - Column decoding fails (e.g., invalid sequence bases)
    /// - Checksum validation fails
    ///
    /// # Example
    ///
    /// ```
    /// use caf::block::{BlockBuilder, BlockReader, AlignmentRecord};
    ///
    /// let mut builder = BlockBuilder::new(0, 10_000);
    /// builder.add_record(AlignmentRecord::default())?;
    /// let block = builder.build()?;
    ///
    /// let reader = BlockReader::new(block)?;
    /// assert_eq!(reader.len(), 1);
    /// # Ok::<(), caf::CafError>(())
    /// ```
    pub fn new(block: CafBlock) -> Result<Self> {
        Self::with_dictionary(block, None)
    }

    /// Create a new block reader with optional quality score dictionary.
    ///
    /// When a dictionary is provided, quality scores compressed with zstd
    /// will be decompressed using the dictionary.
    pub fn with_dictionary(block: CafBlock, quality_dict: Option<&[u8]>) -> Result<Self> {
        // Validate checksum
        Self::validate_checksum(&block)?;

        // Decompress and decode all columns
        let columns = Self::decompress_columns(&block.columns, quality_dict)?;

        Ok(Self {
            num_records: block.num_records as usize,
            columns,
        })
    }

    /// Get number of records in block.
    pub fn len(&self) -> usize {
        self.num_records
    }

    /// Check if block is empty.
    pub fn is_empty(&self) -> bool {
        self.num_records == 0
    }

    /// Get number of records (alias for len()).
    pub fn num_records(&self) -> usize {
        self.num_records
    }

    /// Get quality scores column (column-selective reading).
    ///
    /// Provides direct access to quality scores without decompressing other columns.
    /// This is the key to CAF's streaming analytics efficiency.
    ///
    /// Returns quality scores as concatenated bytes with offset array.
    pub fn qualities(&self) -> (&[u8], &[u32]) {
        (&self.columns.qualities, &self.columns.qual_offsets)
    }

    /// Get MAPQ column (column-selective reading).
    ///
    /// Provides direct access to MAPQ values without decompressing other columns.
    pub fn mapq_values(&self) -> &[u8] {
        &self.columns.mapq
    }

    /// Get sequence column (column-selective reading).
    ///
    /// Provides direct access to sequences without decompressing other columns.
    pub fn sequences(&self) -> (&[u8], &[u32]) {
        (&self.columns.sequences, &self.columns.seq_offsets)
    }

    /// Get a specific record by index.
    ///
    /// # Errors
    ///
    /// Returns error if index is out of bounds.
    ///
    /// # Example
    ///
    /// ```
    /// use caf::block::{BlockBuilder, BlockReader, AlignmentRecord};
    ///
    /// let mut builder = BlockBuilder::new(0, 10_000);
    /// builder.add_record(AlignmentRecord {
    ///     position: 1000,
    ///     ..Default::default()
    /// })?;
    /// let block = builder.build()?;
    ///
    /// let reader = BlockReader::new(block)?;
    /// let record = reader.get_record(0)?;
    /// assert_eq!(record.position, 1000);
    /// # Ok::<(), caf::CafError>(())
    /// ```
    pub fn get_record(&self, index: usize) -> Result<AlignmentRecord> {
        if index >= self.num_records {
            return Err(CafError::Other(format!(
                "Record index {} out of bounds (max {})",
                index, self.num_records
            )));
        }

        self.reconstruct_record(index)
    }

    /// Validate block checksum.
    fn validate_checksum(block: &CafBlock) -> Result<()> {
        let mut hasher = crc32fast::Hasher::new();

        // Hash all compressed column data (same as builder)
        hasher.update(&block.columns.ref_ids.data);
        hasher.update(&block.columns.positions.data);
        hasher.update(&block.columns.mapq.data);
        hasher.update(&block.columns.flags.data);
        hasher.update(&block.columns.sequences.data);
        hasher.update(&block.columns.seq_offsets.data);
        hasher.update(&block.columns.qualities.data);
        hasher.update(&block.columns.qual_offsets.data);
        hasher.update(&block.columns.cigar_ops.data);
        hasher.update(&block.columns.cigar_offsets.data);
        hasher.update(&block.columns.read_names.data);
        hasher.update(&block.columns.read_name_offsets.data);
        hasher.update(&block.columns.mate_ref_ids.data);
        hasher.update(&block.columns.mate_positions.data);
        hasher.update(&block.columns.template_lengths.data);

        let calculated = hasher.finalize();

        if calculated != block.checksum {
            return Err(CafError::ChecksumMismatch {
                block_id: block.block_id,
                expected: block.checksum,
                actual: calculated,
            });
        }

        Ok(())
    }

    /// Decompress and decode all columns.
    fn decompress_columns(columns: &ColumnData, quality_dict: Option<&[u8]>) -> Result<DecompressedColumns> {
        // Decompress integer columns
        let ref_ids = Self::decompress_i32_column(&columns.ref_ids)?;
        let positions_delta = Self::decompress_i32_column(&columns.positions)?;
        let mapq = Self::decompress_u8_column(&columns.mapq)?;
        let flags = Self::decompress_u16_column(&columns.flags)?;

        // Decompress variable-length columns
        let sequences_encoded = Self::decompress_u8_column(&columns.sequences)?;
        let seq_offsets = Self::decompress_u32_column(&columns.seq_offsets)?;

        // Decompress qualities with dictionary support
        let qualities_encoded = Self::decompress_u8_column_with_dict(&columns.qualities, quality_dict)?;
        let qual_offsets = Self::decompress_u32_column(&columns.qual_offsets)?;
        let cigar_ops = Self::decompress_u32_column(&columns.cigar_ops)?;
        let cigar_offsets = Self::decompress_u32_column(&columns.cigar_offsets)?;
        let read_names = Self::decompress_u8_column(&columns.read_names)?;
        let read_name_offsets = Self::decompress_u32_column(&columns.read_name_offsets)?;

        // Decompress mate columns
        let mate_ref_ids = Self::decompress_i32_column(&columns.mate_ref_ids)?;
        let mate_positions_delta = Self::decompress_i32_column(&columns.mate_positions)?;
        let template_lengths = Self::decompress_i32_column(&columns.template_lengths)?;

        // Decode columns
        let positions = decode_integers_delta(&positions_delta);
        let sequences = decode_sequence_ascii(&sequences_encoded)?;
        let qualities = decode_qualities_raw(&qualities_encoded);
        let mate_positions = decode_integers_delta(&mate_positions_delta);

        Ok(DecompressedColumns {
            ref_ids,
            positions,
            mapq,
            flags,
            sequences,
            seq_offsets,
            qualities,
            qual_offsets,
            cigar_ops,
            cigar_offsets,
            read_names,
            read_name_offsets,
            mate_ref_ids,
            mate_positions,
            template_lengths,
        })
    }

    /// Decompress i32 column.
    fn decompress_i32_column(column: &CompressedColumn<i32>) -> Result<Vec<i32>> {
        let bytes = decompress(&column.data, column.compression_type, column.uncompressed_len as usize)?;

        // Convert bytes back to i32 slice
        let count = bytes.len() / std::mem::size_of::<i32>();
        let mut values = Vec::with_capacity(count);

        for chunk in bytes.chunks_exact(std::mem::size_of::<i32>()) {
            let value = i32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
            values.push(value);
        }

        Ok(values)
    }

    /// Decompress u8 column.
    fn decompress_u8_column(column: &CompressedColumn<u8>) -> Result<Vec<u8>> {
        decompress(&column.data, column.compression_type, column.uncompressed_len as usize)
    }

    /// Decompress u8 column with optional dictionary support.
    ///
    /// If the column is zstd-compressed and a dictionary is provided,
    /// uses dictionary decompression. Otherwise falls back to standard decompression.
    fn decompress_u8_column_with_dict(
        column: &CompressedColumn<u8>,
        dict: Option<&[u8]>,
    ) -> Result<Vec<u8>> {
        use crate::types::CompressionType;

        // If column is zstd-compressed and we have a dictionary, use it
        if column.compression_type == CompressionType::Zstd {
            if let Some(dictionary) = dict {
                return crate::compression::decompress_zstd_dict(&column.data, dictionary);
            }
        }

        // Otherwise use standard decompression
        decompress(&column.data, column.compression_type, column.uncompressed_len as usize)
    }

    /// Decompress u16 column.
    fn decompress_u16_column(column: &CompressedColumn<u16>) -> Result<Vec<u16>> {
        let bytes = decompress(&column.data, column.compression_type, column.uncompressed_len as usize)?;

        let count = bytes.len() / std::mem::size_of::<u16>();
        let mut values = Vec::with_capacity(count);

        for chunk in bytes.chunks_exact(std::mem::size_of::<u16>()) {
            let value = u16::from_le_bytes([chunk[0], chunk[1]]);
            values.push(value);
        }

        Ok(values)
    }

    /// Decompress u32 column.
    fn decompress_u32_column(column: &CompressedColumn<u32>) -> Result<Vec<u32>> {
        let bytes = decompress(&column.data, column.compression_type, column.uncompressed_len as usize)?;

        let count = bytes.len() / std::mem::size_of::<u32>();
        let mut values = Vec::with_capacity(count);

        for chunk in bytes.chunks_exact(std::mem::size_of::<u32>()) {
            let value = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
            values.push(value);
        }

        Ok(values)
    }

    /// Reconstruct a record from columns.
    fn reconstruct_record(&self, index: usize) -> Result<AlignmentRecord> {
        // Extract scalar fields
        let ref_id = self.columns.ref_ids[index];
        let position = self.columns.positions[index];
        let mapq = self.columns.mapq[index];
        let flags = self.columns.flags[index];
        let mate_ref_id = self.columns.mate_ref_ids[index];
        let mate_position = self.columns.mate_positions[index];
        let template_length = self.columns.template_lengths[index];

        // Extract variable-length fields using offsets
        let sequence = self.extract_slice(
            &self.columns.sequences,
            &self.columns.seq_offsets,
            index,
        )?;

        let qualities = self.extract_slice(
            &self.columns.qualities,
            &self.columns.qual_offsets,
            index,
        )?;

        let cigar = self.extract_cigar(index)?;

        let read_name = self.extract_slice(
            &self.columns.read_names,
            &self.columns.read_name_offsets,
            index,
        )?;

        Ok(AlignmentRecord {
            ref_id,
            position,
            mapq,
            flags,
            sequence,
            qualities,
            cigar,
            read_name,
            mate_ref_id,
            mate_position,
            template_length,
        })
    }

    /// Extract a slice using offset array.
    fn extract_slice(&self, data: &[u8], offsets: &[u32], index: usize) -> Result<Vec<u8>> {
        let start = offsets[index] as usize;
        let end = offsets[index + 1] as usize;

        if end > data.len() {
            return Err(CafError::Other(format!(
                "Invalid offset: end {} exceeds data length {}",
                end, data.len()
            )));
        }

        // For read names, remove null terminator if present
        let slice = &data[start..end];
        if !slice.is_empty() && slice[slice.len() - 1] == 0 {
            Ok(slice[..slice.len() - 1].to_vec())
        } else {
            Ok(slice.to_vec())
        }
    }

    /// Extract CIGAR operations for a record.
    fn extract_cigar(&self, index: usize) -> Result<Vec<(u8, u32)>> {
        let start = self.columns.cigar_offsets[index] as usize;
        let end = self.columns.cigar_offsets[index + 1] as usize;

        if end > self.columns.cigar_ops.len() {
            return Err(CafError::Other(format!(
                "Invalid CIGAR offset: end {} exceeds ops length {}",
                end, self.columns.cigar_ops.len()
            )));
        }

        let mut cigar = Vec::with_capacity(end - start);
        for &encoded in &self.columns.cigar_ops[start..end] {
            let (op, len) = decode_cigar_op(encoded);
            cigar.push((op, len));
        }

        Ok(cigar)
    }

    /// Convert to iterator over records.
    pub fn into_iter(self) -> BlockRecordIterator {
        BlockRecordIterator {
            reader: self,
            index: 0,
        }
    }
}

/// Iterator over records in a block.
pub struct BlockRecordIterator {
    reader: BlockReader,
    index: usize,
}

impl Iterator for BlockRecordIterator {
    type Item = Result<AlignmentRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.reader.num_records {
            return None;
        }

        let result = self.reader.reconstruct_record(self.index);
        self.index += 1;
        Some(result)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.reader.num_records - self.index;
        (remaining, Some(remaining))
    }
}

impl ExactSizeIterator for BlockRecordIterator {
    fn len(&self) -> usize {
        self.reader.num_records - self.index
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::block::BlockBuilder;

    fn create_test_record(ref_id: i32, position: i32) -> AlignmentRecord {
        AlignmentRecord {
            ref_id,
            position,
            mapq: 60,
            flags: 99,
            sequence: b"ACGT".to_vec(),
            qualities: b"IIII".to_vec(),
            cigar: vec![(0, 4)], // 4M
            read_name: b"read1".to_vec(),
            mate_ref_id: ref_id,
            mate_position: position + 100,
            template_length: 200,
        }
    }

    #[test]
    fn test_reader_single_record() {
        let mut builder = BlockBuilder::new(0, 10_000);
        builder.add_record(create_test_record(0, 1000)).unwrap();
        let block = builder.build().unwrap();

        let reader = BlockReader::new(block).unwrap();
        assert_eq!(reader.len(), 1);
        assert!(!reader.is_empty());

        let record = reader.get_record(0).unwrap();
        assert_eq!(record.ref_id, 0);
        assert_eq!(record.position, 1000);
        assert_eq!(record.mapq, 60);
        assert_eq!(record.sequence, b"ACGT");
    }

    #[test]
    fn test_reader_multiple_records() {
        let mut builder = BlockBuilder::new(0, 10_000);
        for i in 0..100 {
            builder.add_record(create_test_record(0, 1000 + i)).unwrap();
        }
        let block = builder.build().unwrap();

        let reader = BlockReader::new(block).unwrap();
        assert_eq!(reader.len(), 100);

        // Check first and last records
        let first = reader.get_record(0).unwrap();
        assert_eq!(first.position, 1000);

        let last = reader.get_record(99).unwrap();
        assert_eq!(last.position, 1099);
    }

    #[test]
    fn test_reader_iterator() {
        let mut builder = BlockBuilder::new(0, 10_000);
        for i in 0..10 {
            builder.add_record(create_test_record(0, 1000 + i)).unwrap();
        }
        let block = builder.build().unwrap();

        let reader = BlockReader::new(block).unwrap();
        let records: Vec<_> = reader.into_iter().collect::<Result<Vec<_>>>().unwrap();

        assert_eq!(records.len(), 10);
        for (i, record) in records.iter().enumerate() {
            assert_eq!(record.position, 1000 + i as i32);
        }
    }

    #[test]
    fn test_reader_out_of_bounds() {
        let mut builder = BlockBuilder::new(0, 10_000);
        builder.add_record(create_test_record(0, 1000)).unwrap();
        let block = builder.build().unwrap();

        let reader = BlockReader::new(block).unwrap();
        let result = reader.get_record(1);
        assert!(result.is_err());
    }

    #[test]
    fn test_checksum_validation() {
        let mut builder = BlockBuilder::new(0, 10_000);
        builder.add_record(create_test_record(0, 1000)).unwrap();
        let mut block = builder.build().unwrap();

        // Corrupt the checksum
        block.checksum = 0xDEADBEEF;

        let result = BlockReader::new(block);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), CafError::ChecksumMismatch { .. }));
    }

    #[test]
    fn test_round_trip() {
        // Build block with diverse data
        let mut builder = BlockBuilder::new(0, 10_000);

        builder.add_record(AlignmentRecord {
            ref_id: 0,
            position: 1000,
            mapq: 60,
            flags: 99,
            sequence: b"ACGTACGT".to_vec(),
            qualities: b"IIIIIIII".to_vec(),
            cigar: vec![(0, 8)], // 8M
            read_name: b"read1".to_vec(),
            mate_ref_id: 0,
            mate_position: 1100,
            template_length: 200,
        }).unwrap();

        builder.add_record(AlignmentRecord {
            ref_id: 1,
            position: 2000,
            mapq: 30,
            flags: 147,
            sequence: b"GGCCAA".to_vec(),
            qualities: b"HHHHHH".to_vec(),
            cigar: vec![(0, 4), (1, 2)], // 4M2I
            read_name: b"read2".to_vec(),
            mate_ref_id: 1,
            mate_position: 2200,
            template_length: 300,
        }).unwrap();

        let block = builder.build().unwrap();
        let reader = BlockReader::new(block).unwrap();

        // Verify first record
        let r1 = reader.get_record(0).unwrap();
        assert_eq!(r1.ref_id, 0);
        assert_eq!(r1.position, 1000);
        assert_eq!(r1.mapq, 60);
        assert_eq!(r1.flags, 99);
        assert_eq!(r1.sequence, b"ACGTACGT");
        assert_eq!(r1.qualities, b"IIIIIIII");
        assert_eq!(r1.cigar, vec![(0, 8)]);
        assert_eq!(r1.read_name, b"read1");
        assert_eq!(r1.mate_ref_id, 0);
        assert_eq!(r1.mate_position, 1100);
        assert_eq!(r1.template_length, 200);

        // Verify second record
        let r2 = reader.get_record(1).unwrap();
        assert_eq!(r2.ref_id, 1);
        assert_eq!(r2.position, 2000);
        assert_eq!(r2.mapq, 30);
        assert_eq!(r2.flags, 147);
        assert_eq!(r2.sequence, b"GGCCAA");
        assert_eq!(r2.qualities, b"HHHHHH");
        assert_eq!(r2.cigar, vec![(0, 4), (1, 2)]);
        assert_eq!(r2.read_name, b"read2");
        assert_eq!(r2.mate_ref_id, 1);
        assert_eq!(r2.mate_position, 2200);
        assert_eq!(r2.template_length, 300);
    }
}

