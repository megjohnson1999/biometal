//! Core data structures for CAF format.
//!
//! This module defines the binary structures used throughout the CAF format:
//! - `CafHeader`: File header with SAM metadata
//! - `CafBlock`: Columnar block (10,000 records)
//! - `CafIndex`: Block index for random access
//! - `CafFooter`: File footer with index offset
//! - `BlockMeta`: Per-block metadata

use serde::{Deserialize, Serialize};

/// CAF file header containing SAM metadata and format information.
///
/// The header is written at the beginning of the CAF file (after magic number)
/// and contains all necessary metadata for interpreting the columnar blocks.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CafHeader {
    /// Format version (major.minor as u16)
    ///
    /// Version 1.0 = 0x0100
    pub version: u16,

    /// Records per block (default: 10,000 from Rule 2)
    pub block_size: u32,

    /// Total number of blocks in file
    pub num_blocks: u32,

    /// Compressed SAM text header (zstd compressed)
    pub sam_header: Vec<u8>,

    /// Number of reference sequences
    pub num_refs: u32,

    /// Reference sequence names (zstd compressed)
    pub ref_names: Vec<String>,

    /// Reference sequence lengths
    pub ref_lengths: Vec<i32>,

    /// Column schema configuration
    pub column_schema: ColumnSchema,

    /// Optional trained dictionary for quality score compression
    ///
    /// When present, quality scores use zstd dictionary compression
    /// instead of raw storage for 1.5-3× better compression.
    pub quality_dict: Option<Vec<u8>>,
}

impl CafHeader {
    /// Create a new CAF header with default settings.
    pub fn new(block_size: u32, sam_header: Vec<u8>) -> Self {
        Self {
            version: 0x0100,  // Version 1.0
            block_size,
            num_blocks: 0,  // Updated during writing
            sam_header,
            num_refs: 0,
            ref_names: Vec::new(),
            ref_lengths: Vec::new(),
            column_schema: ColumnSchema::default(),
            quality_dict: None,  // Dictionary trained during writing
        }
    }

    /// Get major version number.
    pub fn version_major(&self) -> u8 {
        (self.version >> 8) as u8
    }

    /// Get minor version number.
    pub fn version_minor(&self) -> u8 {
        (self.version & 0xFF) as u8
    }
}

/// Column schema specifying which optional columns are present.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ColumnSchema {
    /// Whether read names are present
    pub has_read_names: bool,

    /// Whether quality scores are present
    pub has_qualities: bool,

    /// Whether CIGAR operations are present
    pub has_cigar: bool,

    /// Whether auxiliary tags are present
    pub has_tags: bool,

    /// Compression configuration per column
    pub compression_config: CompressionConfig,
}

impl Default for ColumnSchema {
    fn default() -> Self {
        Self {
            has_read_names: true,
            has_qualities: true,
            has_cigar: true,
            has_tags: true,
            compression_config: CompressionConfig::default(),
        }
    }
}

/// Compression strategy configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CompressionConfig {
    /// zstd compression level (default: 3)
    pub zstd_level: i32,

    /// Whether to use adaptive RLE for MAPQ
    pub use_adaptive_rle: bool,

    /// Dictionary size for read names (bytes)
    pub dict_size: usize,
}

impl Default for CompressionConfig {
    fn default() -> Self {
        Self {
            zstd_level: 3,           // Balanced ratio vs speed
            use_adaptive_rle: true,  // Adaptive MAPQ compression
            dict_size: 128 * 1024,   // 128 KB dictionary
        }
    }
}

/// Columnar block containing alignment data.
///
/// Each block holds up to `block_size` records (default 10,000) stored in
/// columnar format for efficient SIMD processing.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CafBlock {
    /// Sequential block number
    pub block_id: u32,

    /// Actual number of records in this block (≤ block_size)
    pub num_records: u32,

    /// Uncompressed block size (bytes)
    pub uncompressed_size: u32,

    /// CRC32 checksum of uncompressed data
    pub checksum: u32,

    /// Columnar data
    pub columns: ColumnData,
}

impl CafBlock {
    /// Create a new empty block.
    pub fn new(block_id: u32) -> Self {
        Self {
            block_id,
            num_records: 0,
            uncompressed_size: 0,
            checksum: 0,
            columns: ColumnData::default(),
        }
    }

    /// Check if block is empty.
    pub fn is_empty(&self) -> bool {
        self.num_records == 0
    }

    /// Check if block is full (assuming default block size).
    pub fn is_full(&self, block_size: u32) -> bool {
        self.num_records >= block_size
    }
}

/// Column data for all fields.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ColumnData {
    // Core alignment fields (always present)
    /// Reference sequence IDs
    pub ref_ids: CompressedColumn<i32>,

    /// 0-based leftmost positions
    pub positions: CompressedColumn<i32>,

    /// Mapping quality scores [0, 255]
    pub mapq: CompressedColumn<u8>,

    /// SAM flags (11 bits used)
    pub flags: CompressedColumn<u16>,

    // Sequences (pre-decoded ASCII)
    /// ASCII sequences: 'A', 'C', 'G', 'T', 'N'
    pub sequences: CompressedColumn<u8>,

    /// Cumulative offsets into sequences array
    pub seq_offsets: CompressedColumn<u32>,

    // Quality scores (Phred+33)
    /// ASCII quality scores (optional)
    pub qualities: CompressedColumn<u8>,

    /// Cumulative offsets into qualities array
    pub qual_offsets: CompressedColumn<u32>,

    // CIGAR operations (optional)
    /// BAM CIGAR encoding: (length << 4) | op
    /// - op: 4-bit operation code (M=0, I=1, D=2, N=3, S=4, H=5, P=6, ==7, X=8)
    /// - length: 28-bit operation length
    /// See SAM spec section 4.2.1
    pub cigar_ops: CompressedColumn<u32>,

    /// Per-record CIGAR boundaries
    pub cigar_offsets: CompressedColumn<u32>,

    // Read names (optional)
    /// Read names concatenated as null-terminated strings
    pub read_names: CompressedColumn<u8>,

    /// Cumulative offsets into read_names array
    pub read_name_offsets: CompressedColumn<u32>,

    // Mate information (optional)
    /// Mate reference IDs
    pub mate_ref_ids: CompressedColumn<i32>,

    /// Mate positions
    pub mate_positions: CompressedColumn<i32>,

    /// Template lengths (insert sizes)
    pub template_lengths: CompressedColumn<i32>,
}

/// Compressed column with metadata.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CompressedColumn<T> {
    /// Compression type used
    pub compression_type: CompressionType,

    /// Compressed size (bytes)
    pub compressed_len: u32,

    /// Uncompressed size (bytes)
    pub uncompressed_len: u32,

    /// Compressed data
    pub data: Vec<u8>,

    /// Phantom data for type safety
    #[serde(skip)]
    _phantom: std::marker::PhantomData<T>,
}

impl<T> Default for CompressedColumn<T> {
    fn default() -> Self {
        Self {
            compression_type: CompressionType::Raw,
            compressed_len: 0,
            uncompressed_len: 0,
            data: Vec::new(),
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T> CompressedColumn<T> {
    /// Create a new compressed column.
    pub fn new(
        compression_type: CompressionType,
        compressed_len: u32,
        uncompressed_len: u32,
        data: Vec<u8>,
    ) -> Self {
        Self {
            compression_type,
            compressed_len,
            uncompressed_len,
            data,
            _phantom: std::marker::PhantomData,
        }
    }

    /// Check if column is empty.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Get compression ratio.
    pub fn compression_ratio(&self) -> f64 {
        if self.compressed_len == 0 {
            1.0
        } else {
            self.uncompressed_len as f64 / self.compressed_len as f64
        }
    }
}

/// Compression type enumeration.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[repr(u8)]
pub enum CompressionType {
    /// No compression (raw data)
    Raw = 0,

    /// zstd compression (level specified in config)
    Zstd = 1,

    /// lz4 fast compression
    Lz4 = 2,

    /// Run-length encoding
    Rle = 3,
}

impl CompressionType {
    /// Get compression type name for error messages.
    pub fn name(&self) -> &'static str {
        match self {
            CompressionType::Raw => "raw",
            CompressionType::Zstd => "zstd",
            CompressionType::Lz4 => "lz4",
            CompressionType::Rle => "rle",
        }
    }
}

/// File index for block-level random access.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CafIndex {
    /// File offset to index section (for fast seeking)
    pub index_offset: u64,

    /// Total number of blocks (redundant with header for validation)
    pub num_blocks: u32,

    /// File offset to each block
    pub block_offsets: Vec<u64>,

    /// Metadata for each block
    pub block_metadata: Vec<BlockMeta>,
}

impl CafIndex {
    /// Create a new empty index.
    pub fn new() -> Self {
        Self {
            index_offset: 0,
            num_blocks: 0,
            block_offsets: Vec::new(),
            block_metadata: Vec::new(),
        }
    }

    /// Add a block to the index.
    pub fn add_block(&mut self, offset: u64, meta: BlockMeta) {
        self.block_offsets.push(offset);
        self.block_metadata.push(meta);
        self.num_blocks += 1;
    }
}

impl Default for CafIndex {
    fn default() -> Self {
        Self::new()
    }
}

/// Per-block metadata for indexing and queries.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlockMeta {
    /// Number of records in block
    pub num_records: u32,

    /// Reference sequence ID (-1 = unmapped)
    pub ref_id: i32,

    /// First alignment position (0-based)
    pub start_pos: i32,

    /// Last alignment position (0-based)
    pub end_pos: i32,

    /// Compressed block size (bytes)
    pub compressed_size: u32,

    /// Uncompressed block size (bytes)
    pub uncompressed_size: u32,

    /// CRC32 checksum
    pub checksum: u32,
}

impl BlockMeta {
    /// Check if block overlaps a genomic region.
    pub fn overlaps(&self, ref_id: i32, start: i32, end: i32) -> bool {
        self.ref_id == ref_id && self.start_pos <= end && self.end_pos >= start
    }
}

/// File footer for fast index seeking.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CafFooter {
    /// Offset to index section
    pub index_offset: u64,

    /// Total number of blocks (validation)
    pub num_blocks: u32,

    /// Total number of alignment records
    pub total_records: u64,

    /// CRC32 checksum of entire file
    pub checksum: u32,

    /// Magic number: "CAFE" (end marker)
    pub magic: [u8; 4],

    /// Reserved for future use (zeros)
    pub padding: [u8; 4],
}

impl CafFooter {
    /// Create a new footer.
    pub fn new(index_offset: u64, num_blocks: u32, total_records: u64, checksum: u32) -> Self {
        Self {
            index_offset,
            num_blocks,
            total_records,
            checksum,
            magic: crate::CAF_FOOTER_MAGIC,
            padding: [0; 4],
        }
    }

    /// Validate footer magic number.
    pub fn validate_magic(&self) -> bool {
        self.magic == crate::CAF_FOOTER_MAGIC
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_header_version() {
        let header = CafHeader::new(10_000, vec![]);
        assert_eq!(header.version_major(), 1);
        assert_eq!(header.version_minor(), 0);
    }

    #[test]
    fn test_block_empty_full() {
        let mut block = CafBlock::new(0);
        assert!(block.is_empty());
        assert!(!block.is_full(10_000));

        block.num_records = 10_000;
        assert!(!block.is_empty());
        assert!(block.is_full(10_000));
    }

    #[test]
    fn test_compression_ratio() {
        let column: CompressedColumn<u8> = CompressedColumn::new(
            CompressionType::Zstd,
            100,
            300,
            vec![0; 100],
        );
        assert_eq!(column.compression_ratio(), 3.0);
    }

    #[test]
    fn test_block_meta_overlaps() {
        let meta = BlockMeta {
            num_records: 1000,
            ref_id: 0,
            start_pos: 1000,
            end_pos: 2000,
            compressed_size: 10_000,
            uncompressed_size: 30_000,
            checksum: 0,
        };

        // Overlaps
        assert!(meta.overlaps(0, 500, 1500));  // Spans start
        assert!(meta.overlaps(0, 1500, 2500)); // Spans end
        assert!(meta.overlaps(0, 1200, 1800)); // Contained

        // No overlap
        assert!(!meta.overlaps(0, 0, 999));    // Before
        assert!(!meta.overlaps(0, 2001, 3000)); // After
        assert!(!meta.overlaps(1, 1000, 2000)); // Different ref
    }

    #[test]
    fn test_footer_magic() {
        let footer = CafFooter::new(1024, 10, 100_000, 0xDEADBEEF);
        assert!(footer.validate_magic());
    }
}
