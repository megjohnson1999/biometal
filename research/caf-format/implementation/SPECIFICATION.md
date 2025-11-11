# CAF Format Specification

**Version**: 1.0.0
**Date**: November 10, 2025
**Status**: Production

---

## Table of Contents

1. [Overview](#overview)
2. [Design Principles](#design-principles)
3. [File Structure](#file-structure)
4. [Header Format](#header-format)
5. [Block Format](#block-format)
6. [Column Encoding](#column-encoding)
7. [Compression](#compression)
8. [Dictionary Compression](#dictionary-compression)
9. [Index Structure](#index-structure)
10. [Checksums](#checksums)
11. [Backward Compatibility](#backward-compatibility)
12. [Implementation Notes](#implementation-notes)

---

## Overview

CAF (Columnar Alignment Format) is a columnar storage format for DNA/RNA sequence alignments optimized for ARM NEON processors. It provides:

- **Columnar storage**: Enables vectorized operations on Apple Silicon
- **Block-based organization**: 10,000 records per block for optimal NEON performance
- **Zstandard compression**: Modern compression with dictionary training
- **Fast random access**: Block index for O(1) seeking
- **Lossless conversion**: Full BAM/SAM compatibility

**Target Performance**: 5-10× speedup vs BAM for common operations (base counting, quality filtering)

**Target File Size**: 1.5-2× BAM file size (achieved via dictionary compression)

---

## Design Principles

### 1. ARM-First Architecture

CAF is designed from the ground up for ARM NEON SIMD:
- Columnar layout enables vectorized processing
- 10,000 record blocks maximize NEON throughput (Rule 2)
- Cache-friendly memory access patterns

### 2. Evidence-Based Optimization

All design choices backed by experimental validation:
- **Rule 1**: NEON SIMD for element-wise operations (16-25× speedup)
- **Rule 2**: Block size of 10K records (avoids 82-86% overhead)
- **Rule 5**: Streaming architecture (~5 MB constant memory)

See: OPTIMIZATION_RULES.md in parent biometal project

### 3. Production Quality

- CRC32 checksums on all blocks
- Structured error handling
- Backward compatibility guarantees
- Comprehensive validation

---

## File Structure

```
┌──────────────────┐
│  Magic Number    │  8 bytes: 0x43414601 ("CAF\x01")
├──────────────────┤
│  File Header     │  Variable (bincode serialized)
├──────────────────┤
│  Block 0         │  Up to 10,000 records
├──────────────────┤
│  Block 1         │  Up to 10,000 records
├──────────────────┤
│  ...             │
├──────────────────┤
│  Block N-1       │  Up to 10,000 records
├──────────────────┤
│  Block Index     │  Offset table for random access
├──────────────────┤
│  Footer          │  Metadata + validation
└──────────────────┘
```

### File Layout

1. **Magic Number** (8 bytes): `0x43414601` ("CAF\x01")
2. **Header**: File metadata, reference sequences, dictionary (if present)
3. **Blocks**: Columnar data blocks (10K records each)
4. **Index**: Block offset table for random access
5. **Footer**: File validation metadata

---

## Header Format

The header is bincode-serialized and contains:

```rust
pub struct CafHeader {
    /// Format version (currently 1.0.0)
    pub version: String,

    /// Reference sequence information
    pub references: Vec<ReferenceInfo>,

    /// Number of blocks in file
    pub num_blocks: u32,

    /// Records per block (typically 10,000)
    pub block_size: u32,

    /// Optional trained dictionary for quality score compression
    pub quality_dict: Option<Vec<u8>>,
}

pub struct ReferenceInfo {
    /// Reference sequence name (e.g., "chr1")
    pub name: String,

    /// Reference sequence length in bases
    pub length: u32,
}
```

### Dictionary Compression (v1.0+)

Starting in version 1.0, the header may contain a trained Zstandard dictionary for quality score compression:

- **Field**: `quality_dict: Option<Vec<u8>>`
- **Purpose**: Compress quality scores using patterns learned from the dataset
- **Training**: Dictionary trained on ~30,000 quality score samples
- **Size**: Typically 110 KB (industry standard)
- **Benefit**: 86% file size reduction vs raw quality storage
- **Compatibility**: If absent, reader falls back to standard decompression

See: [Dictionary Compression](#dictionary-compression) section for details.

---

## Block Format

Each block contains up to 10,000 alignment records in columnar form:

```rust
pub struct CafBlock {
    /// Block identifier (0-indexed)
    pub block_id: u32,

    /// Number of records in this block (≤10,000)
    pub num_records: u32,

    /// CRC32 checksum of all column data
    pub checksum: u32,

    /// Compressed columnar data
    pub columns: ColumnData,
}
```

### Column Organization

```rust
pub struct ColumnData {
    /// Reference sequence IDs (compressed i32)
    pub ref_ids: CompressedColumn<i32>,

    /// Alignment positions (compressed i32)
    pub positions: CompressedColumn<i32>,

    /// Mapping quality scores (compressed u8)
    pub mapqs: CompressedColumn<u8>,

    /// SAM flags (compressed u16)
    pub flags: CompressedColumn<u16>,

    /// Mate reference IDs (compressed i32)
    pub mate_ref_ids: CompressedColumn<i32>,

    /// Mate positions (compressed i32)
    pub mate_positions: CompressedColumn<i32>,

    /// Template lengths (compressed i32)
    pub template_lengths: CompressedColumn<i32>,

    /// Read names (compressed string)
    pub read_names: CompressedColumn<String>,

    /// DNA sequences (4-bit encoded, compressed)
    pub sequences: CompressedColumn<u8>,

    /// Quality scores (compressed u8, may use dictionary)
    pub qualities: CompressedColumn<u8>,

    /// CIGAR strings (RLE encoded, compressed)
    pub cigars: CompressedColumn<u8>,

    /// Auxiliary tags (compressed bytes)
    pub tags: CompressedColumn<u8>,
}
```

---

## Column Encoding

### Integer Columns

**Columns**: `ref_ids`, `positions`, `mapqs`, `flags`, etc.

**Encoding**: Raw binary (native endianness)

**Compression**:
- Adaptive based on value distribution
- Small values (0-255): Raw storage
- Larger values: Zstandard level 3

**Rationale**: Simple, fast decompression. Compression ratio improves with value clustering.

### Sequence Column

**Format**: 4-bit encoding (2 bases per byte)

```
Base encoding:
  A = 0x1
  C = 0x2
  G = 0x4
  T = 0x8
  N = 0xF
```

**Packing**: High nibble = first base, low nibble = second base

Example: "ACGT" → `[0x12, 0x48]`

**Compression**: Zstandard level 3 after 4-bit packing

**Rationale**: 50% space reduction vs ASCII, maintains fast decoding

### Quality Scores Column

**Format**: ASCII-33 encoded (Phred+33)

**Compression Strategy**:

1. **With Dictionary** (v1.0+):
   - Use trained Zstandard dictionary from header
   - Achieves ~10× compression ratio
   - 86% file size reduction vs raw storage

2. **Without Dictionary** (backward compatible):
   - Raw storage (no compression)
   - Ensures compatibility with older readers

**Rationale**: Quality scores occupy 70% of compressed file size. Dictionary compression dramatically reduces this without losing information.

### CIGAR Column

**Encoding**: Run-length encoded (RLE)

```rust
pub struct RleCigar {
    /// Operation type (M/I/D/S/H/N/P/X/=)
    pub op: u8,

    /// Operation length
    pub len: u32,
}
```

**Compression**: Zstandard level 3 after RLE encoding

**Rationale**: CIGARs often have long runs (e.g., "100M"). RLE + Zstd achieves ~20× compression.

### Read Names Column

**Encoding**: Length-prefixed UTF-8 strings

```
Format: [len: u32][utf8_bytes: len]
```

**Compression**: Zstandard level 3

**Rationale**: Read names often share common prefixes (e.g., instrument IDs). Zstandard's dictionary learns these patterns.

### Tags Column

**Encoding**: Raw bytes (preserves BAM tag format)

**Compression**: Zstandard level 3

**Rationale**: Tags vary widely. Raw storage + Zstd provides good compression without complex encoding.

---

## Compression

### Compression Types

```rust
pub enum CompressionType {
    /// No compression (fast, large)
    Raw = 0,

    /// Zstandard level 3 (balanced)
    Zstd = 1,

    /// LZ4 (fastest decompression)
    Lz4 = 2,
}
```

### Compression Strategy

**Default**: Zstandard level 3
- Balanced speed vs compression
- ~2× faster than gzip, better ratio
- Well-suited for columnar data

**Adaptive Compression**:

1. **Small columns** (<100 bytes): Raw storage
   - Avoids compression overhead
   - Example: MAPQ column with all values = 60

2. **Medium columns** (100 bytes - 10 KB): Zstandard
   - Good compression ratio
   - Fast decompression

3. **Large columns** (>10 KB): Zstandard
   - Compression essential for file size
   - Amortized decompression cost

### Compressed Column Format

```rust
pub struct CompressedColumn<T> {
    /// Compression type used
    pub compression_type: CompressionType,

    /// Compressed data size in bytes
    pub compressed_len: u32,

    /// Uncompressed data size in bytes
    pub uncompressed_len: u32,

    /// Compressed data bytes
    pub data: Vec<u8>,
}
```

---

## Dictionary Compression

**Version**: 1.0+
**Purpose**: Dramatically reduce file size (86% reduction achieved)
**Target**: Quality score columns (70% of file size)

### Motivation

Quality scores in genomic files:
- Occupy 70% of lossless compressed file size
- Show repetitive patterns within datasets
- Are highly compressible with trained dictionaries

**Problem**: CAF v0.1 stored quality scores raw → 11.6 MB file (11.9× vs BAM)
**Solution**: Zstandard dictionary training → 1.6 MB file (1.6× vs BAM)
**Result**: **86% file size reduction**, achieved 1.5-2× target

### Dictionary Training Process

#### Phase 1: Sample Collection

Collect quality scores from first 30,000 records:

```rust
const SAMPLE_SIZE: usize = 30_000;

let mut quality_samples: Vec<Vec<u8>> = Vec::new();
for record in bam_records.take(SAMPLE_SIZE) {
    if !record.quality.is_empty() {
        quality_samples.push(record.quality.clone());
    }
}
```

**Rationale**: 30K samples provide sufficient pattern coverage while keeping training fast (<20ms)

#### Phase 2: Dictionary Training

Train 110 KB dictionary using Zstandard:

```rust
pub const DICTIONARY_SIZE: usize = 110_000;

pub fn train_dictionary(samples: &[Vec<u8>], dict_size: usize) -> Result<Vec<u8>> {
    // Flatten samples into continuous buffer
    let sample_sizes: Vec<usize> = samples.iter().map(|s| s.len()).collect();
    let total_size: usize = sample_sizes.iter().sum();
    let mut continuous = Vec::with_capacity(total_size);
    for sample in samples {
        continuous.extend_from_slice(sample);
    }

    // Train dictionary
    zstd::dict::from_continuous(&continuous, &sample_sizes, dict_size)
}
```

**Dictionary Size**: 110 KB (industry standard, balances compression vs memory)

#### Phase 3: Compression with Dictionary

Store dictionary in header, use for all quality columns:

```rust
// Writer: Set dictionary before adding records
writer.set_quality_dict(trained_dict)?;

// Compress quality column with dictionary
let compressed = compress_zstd_dict(
    quality_data,
    ZSTD_LEVEL,
    &dictionary,
)?;
```

### Dictionary Storage

**Location**: File header (`CafHeader.quality_dict`)

**Format**: Raw bytes (Zstandard-trained dictionary)

**Size**: Typically 110-112 KB

**Overhead**: One-time cost per file (amortized across all blocks)

### Dictionary Decompression

Reader extracts dictionary from header and uses for all blocks:

```rust
// Reader: Extract dictionary from header
let quality_dict = header.quality_dict.as_deref();

// Decompress quality column with dictionary
fn decompress_quality_column(
    column: &CompressedColumn<u8>,
    dict: Option<&[u8]>,
) -> Result<Vec<u8>> {
    if column.compression_type == CompressionType::Zstd {
        if let Some(dictionary) = dict {
            // Use dictionary decompression
            return decompress_zstd_dict(&column.data, dictionary);
        }
    }

    // Fallback: standard decompression
    decompress(&column.data, column.compression_type, column.uncompressed_len)
}
```

### Performance Impact

**Training Overhead**:
- Collection: ~100 ms (30K records)
- Training: ~10-20 ms
- Total: <150 ms (amortized across entire file)

**Compression Speed**:
- Dictionary overhead: <1% CPU time
- Minimal impact on conversion time (702 ms for 100K records)

**Decompression Speed**:
- Dictionary overhead: ~5% slower than raw
- Still faster than gzip
- 86% size reduction outweighs speed cost

**File Size**:
- Before: 11.6 MB (11.9× vs BAM)
- After: 1.6 MB (1.6× vs BAM)
- Reduction: **86%**
- Target: 1.5-2× (achieved)

### Compatibility

**Forward Compatibility**: Readers that don't support dictionaries:
- Detect `CompressionType::Zstd` on quality column
- Attempt standard zstd decompression (will fail with "Dictionary mismatch")
- Should check header version and warn user to upgrade

**Backward Compatibility**: Files without dictionary:
- `quality_dict` field is `None`
- Reader falls back to standard decompression
- No performance degradation

**Version Detection**:

```rust
// Check if file uses dictionary compression
if header.quality_dict.is_some() {
    // Version 1.0+: Dictionary compression
    // Use decompress_zstd_dict()
} else {
    // Version 0.1: Raw quality scores
    // Use standard decompress()
}
```

---

## Index Structure

The block index enables fast random access:

```rust
pub struct BlockIndex {
    /// Number of blocks in file
    pub num_blocks: u32,

    /// Byte offset of each block from file start
    pub block_offsets: Vec<u64>,
}
```

**Location**: End of file, before footer

**Purpose**: O(1) seeking to any block without scanning

**Example**:
```
Block 0: offset = 1024
Block 1: offset = 124512
Block 2: offset = 248901
...
```

**Usage**:
```rust
// Seek to block 42
let offset = index.block_offsets[42];
reader.seek(SeekFrom::Start(offset))?;
let block = reader.read_block()?;
```

---

## Checksums

### Block Checksums

**Algorithm**: CRC32C (hardware-accelerated on ARM)

**Coverage**: All column data in block

**Validation**: On block read, before decompression

```rust
pub fn validate_checksum(block: &CafBlock) -> Result<()> {
    let computed = compute_checksum(&block.columns);
    if computed != block.checksum {
        return Err(CafError::ChecksumMismatch {
            block_id: block.block_id,
            expected: block.checksum,
            actual: computed,
        });
    }
    Ok(())
}
```

**Rationale**: Detects corruption during I/O or storage. CRC32C chosen for ARM hardware acceleration (Rule 1).

### File Checksums

**Future**: Footer may include file-level checksum (not yet implemented)

---

## Backward Compatibility

### Version Policy

CAF follows semantic versioning:

- **Major version**: Breaking changes (e.g., column format change)
- **Minor version**: Backward-compatible additions (e.g., new compression type)
- **Patch version**: Bug fixes

Current version: **1.0.0**

### Reading Old Files

Version 1.0 readers MUST support:

- **v0.1**: Files without `quality_dict` field
  - Fallback to raw quality decompression
  - No performance degradation

### Writing Compatible Files

Writers SHOULD:

- Include version in header (`version: "1.0.0"`)
- Use dictionary compression by default (for file size)
- Allow disabling dictionary for v0.1 compatibility (if needed)

---

## Implementation Notes

### Memory Management

**Streaming Architecture**:
- Read one block at a time (~5 MB uncompressed)
- Decompress columns on demand
- Avoid accumulating records in memory

**Target**: Constant ~5 MB memory usage (Rule 5)

### Error Handling

All operations return `Result<T, CafError>`:

```rust
pub enum CafError {
    /// I/O errors
    Io(std::io::Error),

    /// Compression errors
    CompressionError {
        column: String,
        block_id: u32,
        source: Box<dyn std::error::Error>,
    },

    /// Decompression errors
    DecompressionError {
        column: String,
        block_id: u32,
        source: Box<dyn std::error::Error>,
    },

    /// Checksum validation failures
    ChecksumMismatch {
        block_id: u32,
        expected: u32,
        actual: u32,
    },

    /// Invalid format
    InvalidFormat(String),

    /// Other errors
    Other(String),
}
```

**No panics in library code** - always use `Result`.

### NEON Optimization

CAF's columnar layout enables ARM NEON vectorization:

**Example**: Base counting on sequence column
```rust
#[cfg(target_arch = "aarch64")]
pub unsafe fn count_bases_neon(sequences: &[u8]) -> BaseCounts {
    // Process 16 bases at once with NEON
    // Expected: 25× speedup vs scalar
}
```

**Target Operations**:
- Base counting (25× speedup)
- Quality filtering (25× speedup)
- MAPQ filtering (16× speedup)

See: biometal OPTIMIZATION_RULES.md for evidence base.

### Testing

**Unit Tests**: Column encoding/decoding, compression, checksums

**Integration Tests**:
- Round-trip BAM → CAF → SAM (lossless)
- Dictionary compression validation
- Checksum validation

**Property-Based Tests**:
- NEON == scalar correctness
- Round-trip preservation
- Arbitrary record generation

**Benchmarks** (N=30):
- Quality filtering (Q30)
- Base counting
- MAPQ filtering
- Full parsing (100K records)

---

## Appendix: Design Trade-offs

### Why 10K Records Per Block?

**Evidence**: Entry 027 from apple-silicon-bio-bench (1,440 measurements)

**Finding**: Single-record processing incurs 82-86% overhead on NEON operations

**Solution**: Block-based processing amortizes overhead, maintains NEON speedup

**Trade-off**: Larger blocks = more memory, but 10K records = ~5 MB (acceptable)

### Why Zstandard Over Gzip?

**Evidence**: Entry 029 (CPU parallel prototype)

**Finding**: Zstd decompresses 2× faster than gzip with better compression ratio

**Benefit**: Enables faster file parsing, smaller files

**Trade-off**: Requires zstd library, but widely available (Rust crate, system libs)

### Why Dictionary Compression?

**Evidence**: Achieved 86% file size reduction in testing

**Finding**: Quality scores are highly repetitive within datasets

**Benefit**: CAF achieves target file size (1.6× vs BAM) vs 11.9× without dictionary

**Trade-off**: Training overhead (~150 ms), but amortized across entire file

### Why Columnar Format?

**Evidence**: Enables ARM NEON vectorization (16-25× speedup)

**Finding**: Row-based formats (BAM) require deinterleaving before SIMD

**Benefit**: Direct SIMD operations on columns (quality filtering, base counting)

**Trade-off**: Random record access slower, but block index mitigates this

---

## Future Extensions

**Potential Additions** (not yet implemented):

1. **Read Name Dictionary**: Similar to quality dictionary, for read name compression
2. **Parallel Block Compression**: Multi-threaded compression during conversion
3. **Streaming Decompression**: Decompress next block while processing current
4. **Region Queries**: Index by reference position for fast region extraction
5. **GPU Acceleration**: Metal/Vulkan for sequence operations

---

## References

1. BAM Format Specification v1.6
2. SAM Format Specification v1.6
3. Zstandard Compression (RFC 8878)
4. biometal OPTIMIZATION_RULES.md (evidence base)
5. apple-silicon-bio-bench (1,357 experiments)

---

## Changelog

### Version 1.0.0 (November 10, 2025)

**Added**:
- Dictionary compression for quality scores (86% file size reduction)
- `quality_dict` field in `CafHeader`
- Dictionary training from 30K samples
- Dictionary-aware compression/decompression
- Backward compatibility with v0.1

**Performance**:
- File size: 11.6 MB → 1.6 MB (86% reduction)
- Target achieved: 1.6× vs BAM (within 1.5-2× target)
- Round-trip verified: BAM → CAF → SAM lossless

**Implementation**:
- ~305 lines of code (writer + reader)
- Training overhead: <150 ms
- Conversion time: 702 ms for 100K records

---

**Document Version**: 1.0.0
**Last Updated**: November 10, 2025
