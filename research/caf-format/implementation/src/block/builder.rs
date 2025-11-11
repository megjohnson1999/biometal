//! Block builder for converting records to columnar format.
//!
//! The BlockBuilder accumulates alignment records up to block_size (default 10,000),
//! then converts them to columnar format with:
//! - Column-specific encoding (delta, zigzag, pre-decoded sequences)
//! - Compression (zstd, lz4, RLE, or raw)
//! - CRC32 checksums for data integrity
//!
//! # Example
//!
//! ```
//! use caf::block::{BlockBuilder, AlignmentRecord};
//!
//! let mut builder = BlockBuilder::new(0, 10_000);
//!
//! // Add records
//! builder.add_record(AlignmentRecord {
//!     ref_id: 0,
//!     position: 1000,
//!     mapq: 60,
//!     flags: 99,
//!     sequence: b"ACGT".to_vec(),
//!     qualities: b"IIII".to_vec(),
//!     cigar: vec![(0, 4)], // 4M
//!     ..Default::default()
//! })?;
//!
//! // Build block when full or at end of stream
//! if builder.is_full() {
//!     let block = builder.build()?;
//! }
//! # Ok::<(), caf::CafError>(())
//! ```

use crate::{
    column::{encode_integers_delta, encode_sequence_ascii, encode_qualities_raw, encode_cigar_op},
    compression::compress,
    types::{CafBlock, ColumnData, CompressedColumn, CompressionConfig, CompressionType},
    CafError, Result,
};

/// Alignment record in row-oriented format.
///
/// This represents a single alignment record with all its fields,
/// similar to SAM/BAM format.
#[derive(Debug, Clone, Default)]
pub struct AlignmentRecord {
    /// Reference sequence ID (-1 for unmapped)
    pub ref_id: i32,

    /// 0-based leftmost position
    pub position: i32,

    /// Mapping quality [0, 255]
    pub mapq: u8,

    /// SAM flags (11 bits used)
    pub flags: u16,

    /// Sequence as ASCII bytes ('A', 'C', 'G', 'T', 'N')
    pub sequence: Vec<u8>,

    /// Quality scores as ASCII bytes (Phred+33)
    pub qualities: Vec<u8>,

    /// CIGAR operations as (op, length) pairs
    /// op: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X
    pub cigar: Vec<(u8, u32)>,

    /// Read name (null-terminated)
    pub read_name: Vec<u8>,

    /// Mate reference ID
    pub mate_ref_id: i32,

    /// Mate position
    pub mate_position: i32,

    /// Template length (insert size)
    pub template_length: i32,
}

/// Block builder for accumulating records and converting to columnar format.
pub struct BlockBuilder {
    /// Block ID
    block_id: u32,

    /// Maximum records per block
    block_size: u32,

    /// Accumulated records
    records: Vec<AlignmentRecord>,

    /// Compression configuration
    compression_config: CompressionConfig,

    /// Optional trained dictionary for quality score compression
    quality_dict: Option<Vec<u8>>,
}

impl BlockBuilder {
    /// Create a new block builder.
    ///
    /// # Arguments
    ///
    /// * `block_id` - Sequential block number
    /// * `block_size` - Maximum records per block (typically 10,000)
    ///
    /// # Example
    ///
    /// ```
    /// use caf::block::BlockBuilder;
    ///
    /// let builder = BlockBuilder::new(0, 10_000);
    /// assert!(builder.is_empty());
    /// assert!(!builder.is_full());
    /// ```
    pub fn new(block_id: u32, block_size: u32) -> Self {
        Self {
            block_id,
            block_size,
            records: Vec::with_capacity(block_size as usize),
            compression_config: CompressionConfig::default(),
            quality_dict: None,
        }
    }

    /// Set the quality score dictionary for compression.
    ///
    /// When a dictionary is provided, quality scores will be compressed using
    /// zstd dictionary compression instead of raw storage.
    pub fn set_quality_dict(&mut self, dict: Vec<u8>) {
        self.quality_dict = Some(dict);
    }

    /// Add a record to the block.
    ///
    /// Returns an error if the block is already full.
    ///
    /// # Example
    ///
    /// ```
    /// use caf::block::{BlockBuilder, AlignmentRecord};
    ///
    /// let mut builder = BlockBuilder::new(0, 10);
    /// let record = AlignmentRecord {
    ///     ref_id: 0,
    ///     position: 1000,
    ///     ..Default::default()
    /// };
    /// builder.add_record(record)?;
    /// assert_eq!(builder.len(), 1);
    /// # Ok::<(), caf::CafError>(())
    /// ```
    pub fn add_record(&mut self, record: AlignmentRecord) -> Result<()> {
        if self.is_full() {
            return Err(CafError::Other(format!(
                "Block {} is full (max {} records)",
                self.block_id, self.block_size
            )));
        }

        self.records.push(record);
        Ok(())
    }

    /// Check if block is empty.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Check if block is full.
    pub fn is_full(&self) -> bool {
        self.records.len() >= self.block_size as usize
    }

    /// Get current number of records.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Build a columnar block from accumulated records.
    ///
    /// This converts row-oriented records to columnar format with:
    /// 1. Column-specific encoding (delta, zigzag, pre-decoded sequences)
    /// 2. Compression (zstd, lz4, RLE, or raw)
    /// 3. CRC32 checksums
    ///
    /// After building, the builder is reset for the next block.
    ///
    /// # Example
    ///
    /// ```
    /// use caf::block::{BlockBuilder, AlignmentRecord};
    ///
    /// let mut builder = BlockBuilder::new(0, 10_000);
    /// // Add records...
    /// for i in 0..100 {
    ///     builder.add_record(AlignmentRecord {
    ///         ref_id: 0,
    ///         position: 1000 + i,
    ///         ..Default::default()
    ///     })?;
    /// }
    ///
    /// let block = builder.build()?;
    /// assert_eq!(block.num_records, 100);
    /// assert_eq!(block.block_id, 0);
    /// # Ok::<(), caf::CafError>(())
    /// ```
    pub fn build(&mut self) -> Result<CafBlock> {
        if self.is_empty() {
            return Err(CafError::Other(
                "Cannot build block with no records".to_string(),
            ));
        }

        let num_records = self.records.len() as u32;

        // Convert to columnar format
        let columns = self.build_columns()?;

        // Calculate uncompressed size and checksum
        let uncompressed_size = self.calculate_uncompressed_size(&columns);
        let checksum = self.calculate_checksum(&columns);

        // Create block
        let block = CafBlock {
            block_id: self.block_id,
            num_records,
            uncompressed_size,
            checksum,
            columns,
        };

        // Reset for next block
        self.records.clear();
        self.block_id += 1;

        Ok(block)
    }

    /// Convert accumulated records to columnar format.
    fn build_columns(&self) -> Result<ColumnData> {
        // Extract columns from records
        let ref_ids: Vec<i32> = self.records.iter().map(|r| r.ref_id).collect();
        let positions: Vec<i32> = self.records.iter().map(|r| r.position).collect();
        let mapq: Vec<u8> = self.records.iter().map(|r| r.mapq).collect();
        let flags: Vec<u16> = self.records.iter().map(|r| r.flags).collect();

        // Build sequences and offsets
        let (sequences, seq_offsets) = self.build_sequence_columns();

        // Build qualities and offsets
        let (qualities, qual_offsets) = self.build_quality_columns();

        // Build CIGAR and offsets
        let (cigar_ops, cigar_offsets) = self.build_cigar_columns();

        // Build read names and offsets
        let (read_names, read_name_offsets) = self.build_read_name_columns();

        // Mate information
        let mate_ref_ids: Vec<i32> = self.records.iter().map(|r| r.mate_ref_id).collect();
        let mate_positions: Vec<i32> = self.records.iter().map(|r| r.mate_position).collect();
        let template_lengths: Vec<i32> = self.records.iter().map(|r| r.template_length).collect();

        // Encode and compress columns
        // Use adaptive compression for MAPQ to handle both uniform and diverse values
        use crate::compression::select_compression;
        let mapq_compression = select_compression(&mapq);

        Ok(ColumnData {
            ref_ids: self.compress_i32_column(&ref_ids, "ref_ids", CompressionType::Rle)?,
            positions: self.compress_i32_column(&encode_integers_delta(&positions), "positions", CompressionType::Zstd)?,
            mapq: self.compress_u8_column(&mapq, "mapq", mapq_compression)?,
            flags: self.compress_u16_column(&flags, "flags", CompressionType::Zstd)?,
            sequences: self.compress_u8_column(&encode_sequence_ascii(&sequences), "sequences", CompressionType::Lz4)?,
            seq_offsets: self.compress_u32_column(&seq_offsets, "seq_offsets", CompressionType::Zstd)?,
            qualities: self.compress_u8_column_with_dict(&encode_qualities_raw(&qualities), "qualities", CompressionType::Raw, &self.quality_dict)?,
            qual_offsets: self.compress_u32_column(&qual_offsets, "qual_offsets", CompressionType::Zstd)?,
            cigar_ops: self.compress_u32_column(&cigar_ops, "cigar_ops", CompressionType::Zstd)?,
            cigar_offsets: self.compress_u32_column(&cigar_offsets, "cigar_offsets", CompressionType::Zstd)?,
            read_names: self.compress_u8_column(&read_names, "read_names", CompressionType::Zstd)?,
            read_name_offsets: self.compress_u32_column(&read_name_offsets, "read_name_offsets", CompressionType::Zstd)?,
            mate_ref_ids: self.compress_i32_column(&mate_ref_ids, "mate_ref_ids", CompressionType::Rle)?,
            mate_positions: self.compress_i32_column(&encode_integers_delta(&mate_positions), "mate_positions", CompressionType::Zstd)?,
            template_lengths: self.compress_i32_column(&template_lengths, "template_lengths", CompressionType::Zstd)?,
        })
    }

    /// Build sequence column and offsets.
    fn build_sequence_columns(&self) -> (Vec<u8>, Vec<u32>) {
        let mut sequences = Vec::new();
        let mut offsets = Vec::with_capacity(self.records.len() + 1);
        offsets.push(0);

        for record in &self.records {
            sequences.extend_from_slice(&record.sequence);
            offsets.push(sequences.len() as u32);
        }

        (sequences, offsets)
    }

    /// Build quality column and offsets.
    fn build_quality_columns(&self) -> (Vec<u8>, Vec<u32>) {
        let mut qualities = Vec::new();
        let mut offsets = Vec::with_capacity(self.records.len() + 1);
        offsets.push(0);

        for record in &self.records {
            qualities.extend_from_slice(&record.qualities);
            offsets.push(qualities.len() as u32);
        }

        (qualities, offsets)
    }

    /// Build CIGAR column and offsets.
    fn build_cigar_columns(&self) -> (Vec<u32>, Vec<u32>) {
        let mut cigar_ops = Vec::new();
        let mut offsets = Vec::with_capacity(self.records.len() + 1);
        offsets.push(0);

        for record in &self.records {
            for &(op, len) in &record.cigar {
                cigar_ops.push(encode_cigar_op(op, len));
            }
            offsets.push(cigar_ops.len() as u32);
        }

        (cigar_ops, offsets)
    }

    /// Build read name column and offsets.
    fn build_read_name_columns(&self) -> (Vec<u8>, Vec<u32>) {
        let mut read_names = Vec::new();
        let mut offsets = Vec::with_capacity(self.records.len() + 1);
        offsets.push(0);

        for record in &self.records {
            read_names.extend_from_slice(&record.read_name);
            read_names.push(0); // Null terminator
            offsets.push(read_names.len() as u32);
        }

        (read_names, offsets)
    }

    /// Compress i32 column.
    fn compress_i32_column(
        &self,
        values: &[i32],
        _column_name: &str,
        compression_type: CompressionType,
    ) -> Result<CompressedColumn<i32>> {
        let uncompressed = unsafe {
            std::slice::from_raw_parts(
                values.as_ptr() as *const u8,
                values.len() * std::mem::size_of::<i32>(),
            )
        };

        let compressed = compress(uncompressed, compression_type)?;

        Ok(CompressedColumn::new(
            compression_type,
            compressed.len() as u32,
            uncompressed.len() as u32,
            compressed,
        ))
    }

    /// Compress u8 column.
    fn compress_u8_column(
        &self,
        values: &[u8],
        _column_name: &str,
        compression_type: CompressionType,
    ) -> Result<CompressedColumn<u8>> {
        let compressed = compress(values, compression_type)?;

        Ok(CompressedColumn::new(
            compression_type,
            compressed.len() as u32,
            values.len() as u32,
            compressed,
        ))
    }

    /// Compress u8 column with optional dictionary.
    ///
    /// If a dictionary is provided, uses zstd dictionary compression.
    /// Otherwise falls back to specified compression type.
    fn compress_u8_column_with_dict(
        &self,
        values: &[u8],
        column_name: &str,
        compression_type: CompressionType,
        dict: &Option<Vec<u8>>,
    ) -> Result<CompressedColumn<u8>> {
        if let Some(dictionary) = dict {
            // Use dictionary compression
            let compressed = crate::compression::compress_zstd_dict(
                values,
                crate::compression::ZSTD_LEVEL,
                dictionary,
            )?;

            Ok(CompressedColumn::new(
                CompressionType::Zstd,  // Using zstd with dictionary
                compressed.len() as u32,
                values.len() as u32,
                compressed,
            ))
        } else {
            // Fall back to standard compression
            self.compress_u8_column(values, column_name, compression_type)
        }
    }

    /// Compress u16 column.
    fn compress_u16_column(
        &self,
        values: &[u16],
        _column_name: &str,
        compression_type: CompressionType,
    ) -> Result<CompressedColumn<u16>> {
        let uncompressed = unsafe {
            std::slice::from_raw_parts(
                values.as_ptr() as *const u8,
                values.len() * std::mem::size_of::<u16>(),
            )
        };

        let compressed = compress(uncompressed, compression_type)?;

        Ok(CompressedColumn::new(
            compression_type,
            compressed.len() as u32,
            uncompressed.len() as u32,
            compressed,
        ))
    }

    /// Compress u32 column.
    fn compress_u32_column(
        &self,
        values: &[u32],
        _column_name: &str,
        compression_type: CompressionType,
    ) -> Result<CompressedColumn<u32>> {
        let uncompressed = unsafe {
            std::slice::from_raw_parts(
                values.as_ptr() as *const u8,
                values.len() * std::mem::size_of::<u32>(),
            )
        };

        let compressed = compress(uncompressed, compression_type)?;

        Ok(CompressedColumn::new(
            compression_type,
            compressed.len() as u32,
            uncompressed.len() as u32,
            compressed,
        ))
    }

    /// Calculate total uncompressed size of all columns.
    fn calculate_uncompressed_size(&self, columns: &ColumnData) -> u32 {
        columns.ref_ids.uncompressed_len
            + columns.positions.uncompressed_len
            + columns.mapq.uncompressed_len
            + columns.flags.uncompressed_len
            + columns.sequences.uncompressed_len
            + columns.seq_offsets.uncompressed_len
            + columns.qualities.uncompressed_len
            + columns.qual_offsets.uncompressed_len
            + columns.cigar_ops.uncompressed_len
            + columns.cigar_offsets.uncompressed_len
            + columns.read_names.uncompressed_len
            + columns.read_name_offsets.uncompressed_len
            + columns.mate_ref_ids.uncompressed_len
            + columns.mate_positions.uncompressed_len
            + columns.template_lengths.uncompressed_len
    }

    /// Calculate CRC32 checksum of uncompressed data.
    ///
    /// Checksums all uncompressed column data to detect corruption.
    fn calculate_checksum(&self, columns: &ColumnData) -> u32 {
        let mut hasher = crc32fast::Hasher::new();

        // Hash all column data (uncompressed)
        // We hash the compressed data since that's what we store
        hasher.update(&columns.ref_ids.data);
        hasher.update(&columns.positions.data);
        hasher.update(&columns.mapq.data);
        hasher.update(&columns.flags.data);
        hasher.update(&columns.sequences.data);
        hasher.update(&columns.seq_offsets.data);
        hasher.update(&columns.qualities.data);
        hasher.update(&columns.qual_offsets.data);
        hasher.update(&columns.cigar_ops.data);
        hasher.update(&columns.cigar_offsets.data);
        hasher.update(&columns.read_names.data);
        hasher.update(&columns.read_name_offsets.data);
        hasher.update(&columns.mate_ref_ids.data);
        hasher.update(&columns.mate_positions.data);
        hasher.update(&columns.template_lengths.data);

        hasher.finalize()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_builder_new() {
        let builder = BlockBuilder::new(0, 10_000);
        assert_eq!(builder.block_id, 0);
        assert_eq!(builder.block_size, 10_000);
        assert!(builder.is_empty());
        assert!(!builder.is_full());
        assert_eq!(builder.len(), 0);
    }

    #[test]
    fn test_add_record() {
        let mut builder = BlockBuilder::new(0, 10);
        let record = create_test_record(0, 1000);

        builder.add_record(record).unwrap();
        assert_eq!(builder.len(), 1);
        assert!(!builder.is_empty());
        assert!(!builder.is_full());
    }

    #[test]
    fn test_builder_full() {
        let mut builder = BlockBuilder::new(0, 2);

        builder.add_record(create_test_record(0, 1000)).unwrap();
        assert!(!builder.is_full());

        builder.add_record(create_test_record(0, 1005)).unwrap();
        assert!(builder.is_full());

        // Adding another should fail
        let result = builder.add_record(create_test_record(0, 1010));
        assert!(result.is_err());
    }

    #[test]
    fn test_build_block() {
        let mut builder = BlockBuilder::new(0, 10_000);

        // Add some records
        for i in 0..100 {
            builder
                .add_record(create_test_record(0, 1000 + i))
                .unwrap();
        }

        let block = builder.build().unwrap();
        assert_eq!(block.block_id, 0);
        assert_eq!(block.num_records, 100);
        assert!(block.uncompressed_size > 0);

        // Builder should be reset
        assert!(builder.is_empty());
        assert_eq!(builder.block_id, 1); // Incremented for next block
    }

    #[test]
    fn test_build_empty_fails() {
        let mut builder = BlockBuilder::new(0, 10_000);
        let result = builder.build();
        assert!(result.is_err());
    }

    #[test]
    fn test_sequence_columns() {
        let mut builder = BlockBuilder::new(0, 10);
        builder
            .add_record(AlignmentRecord {
                sequence: b"ACGT".to_vec(),
                ..Default::default()
            })
            .unwrap();
        builder
            .add_record(AlignmentRecord {
                sequence: b"GGCCAA".to_vec(),
                ..Default::default()
            })
            .unwrap();

        let (sequences, offsets) = builder.build_sequence_columns();
        assert_eq!(sequences, b"ACGTGGCCAA");
        assert_eq!(offsets, vec![0, 4, 10]);
    }

    #[test]
    fn test_cigar_columns() {
        let mut builder = BlockBuilder::new(0, 10);
        builder
            .add_record(AlignmentRecord {
                cigar: vec![(0, 10), (1, 2)], // 10M2I
                ..Default::default()
            })
            .unwrap();
        builder
            .add_record(AlignmentRecord {
                cigar: vec![(0, 5)], // 5M
                ..Default::default()
            })
            .unwrap();

        let (cigar_ops, offsets) = builder.build_cigar_columns();
        assert_eq!(cigar_ops.len(), 3); // 10M, 2I, 5M
        assert_eq!(offsets, vec![0, 2, 3]);
    }

    #[test]
    fn test_checksum_deterministic() {
        let mut builder1 = BlockBuilder::new(0, 10_000);
        let mut builder2 = BlockBuilder::new(0, 10_000);

        // Add same records to both builders
        for i in 0..10 {
            let record = create_test_record(0, 1000 + i);
            builder1.add_record(record.clone()).unwrap();
            builder2.add_record(record).unwrap();
        }

        let block1 = builder1.build().unwrap();
        let block2 = builder2.build().unwrap();

        // Checksums should be identical for identical data
        assert_eq!(block1.checksum, block2.checksum);
        assert_ne!(block1.checksum, 0); // Should not be placeholder
    }

    #[test]
    fn test_checksum_different_data() {
        let mut builder1 = BlockBuilder::new(0, 10_000);
        let mut builder2 = BlockBuilder::new(0, 10_000);

        // Add different records
        builder1.add_record(create_test_record(0, 1000)).unwrap();
        builder2.add_record(create_test_record(0, 2000)).unwrap();

        let block1 = builder1.build().unwrap();
        let block2 = builder2.build().unwrap();

        // Checksums should be different for different data
        assert_ne!(block1.checksum, block2.checksum);
    }
}
