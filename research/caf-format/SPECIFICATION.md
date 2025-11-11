# CAF Specification v1.0

**Format Name**: Columnar Alignment Format (CAF)
**Version**: 1.0.0
**Status**: Research Specification (Finalized)
**Authors**: Scott Handley
**Last Updated**: November 10, 2025
**License**: MIT

---

## Document History

| Version | Date | Changes |
|---------|------|---------|
| 1.0.0 | 2025-11-10 | Initial finalized specification for publication |
| 0.9.0 | 2025-11-10 | Draft from experiments/columnar-alignment-format/ |

---

## Abstract

CAF is a **columnar binary format** for DNA/RNA sequence alignments optimized for analytical bioinformatics operations on modern ARM hardware. Unlike row-oriented formats (BAM [1], CRAM [3]), CAF stores alignment data in columnar blocks enabling:

1. **ARM NEON SIMD vectorization** [10, 11]: 16-25× speedup for element-wise operations
2. **Modern compression** [17, 19]: zstd/lz4 (2-8× faster than gzip)
3. **Analytical query optimization** [6, 7]: Predicate pushdown, column skipping
4. **Hardware-aware design** [21]: Optimized for 2025+ systems (ARM, NVMe, multi-core)

**Design Goal**: 5-10× performance improvement over BAM for analytical operations while maintaining lossless conversion compatibility.

**File Extension**: `.caf`
**Index Extension**: `.caf.idx`
**MIME Type**: `application/x-caf` (proposed)

---

## Design Rationale

CAF addresses the gap between BAM's 2009 design constraints [1] and 2025 hardware capabilities [21, 22, 23]:

| Constraint | 2009 (BAM era) | 2025 (CAF era) | CAF Design |
|------------|----------------|----------------|------------|
| **Storage** | HDD: $100/TB, slow | NVMe: $10/TB, fast | Optimize for speed over compression |
| **CPU** | Single-core, cheap | 8-16 cores, expensive | SIMD + parallel decompression |
| **Memory** | 4-8 GB, scarce | 16-64 GB, abundant | Larger working sets (10K blocks) |
| **SIMD** | 128-bit SSE2, limited | 128-bit NEON, standard | Columnar layout for vectorization |
| **Compression** | gzip (1992) | zstd/lz4 (2015+) | Modern algorithms |

**Evidence base**: 6 optimization rules from 1,357 experiments (biometal OPTIMIZATION_RULES.md)

---

## File Structure

### High-Level Layout

```
Offset  | Size      | Description
--------|-----------|------------------------------------------
0       | 4 bytes   | Magic number: "CAF\x01"
4       | Variable  | Header (compressed SAM header + metadata)
?       | Variable  | Index (block offsets + metadata)
?       | Variable  | Data blocks (columnar, compressed)
```

**Binary encoding**: Little-endian (consistent with BAM [1])
**Alignment**: All integer fields 4-byte aligned
**Serialization**: bincode v1.3 (Rust) or equivalent

### Magic Number (4 bytes)

```
Offset | Type | Value      | Description
-------|------|------------|-------------------------
0      | u8   | 0x43 ('C') | Format identifier
1      | u8   | 0x41 ('A') | Format identifier
2      | u8   | 0x46 ('F') | Format identifier
3      | u8   | 0x01       | Major version (allows evolution)
```

**Version encoding**: Byte 3 = major version (1.0.x → 0x01, 2.0.x → 0x02)
**Future compatibility**: Readers must check magic number and version

### Header Structure

```rust
struct CafHeader {
    // Version (2 bytes)
    version: u16,                  // 1.0 = 0x0100 (major.minor)

    // Block configuration (8 bytes)
    block_size: u32,               // Records/block (default: 10,000)
    num_blocks: u32,               // Total blocks in file

    // SAM header (variable)
    sam_header_len: u32,           // Compressed length
    sam_header: Vec<u8>,           // zstd compressed SAM text header

    // Reference sequences (variable)
    num_refs: u32,                 // Number of reference sequences
    ref_names: Vec<String>,        // Reference names (zstd compressed)
    ref_lengths: Vec<i32>,         // Reference lengths

    // Column schema (variable)
    column_schema: ColumnSchema,   // Column types and compression
}

struct ColumnSchema {
    // Which optional columns are present
    has_read_names: bool,
    has_qualities: bool,
    has_cigar: bool,
    has_tags: bool,

    // Compression strategy per column
    compression_config: CompressionConfig,
}
```

**Serialization order**:
1. version (u16, 2 bytes)
2. block_size (u32, 4 bytes)
3. num_blocks (u32, 4 bytes)
4. sam_header_len (u32, 4 bytes)
5. sam_header (Vec<u8>, variable)
6. num_refs (u32, 4 bytes)
7. ref_names + ref_lengths (variable)
8. column_schema (variable)

### Index Structure

**Stored in main file** (not separate like BAI, for simplicity)

```rust
struct CafIndex {
    // Index metadata (12 bytes)
    index_offset: u64,             // File offset to index (for seeking)
    num_blocks: u32,               // Redundant with header (validation)

    // Block offsets (8 bytes × num_blocks)
    block_offsets: Vec<u64>,       // File offset to each block

    // Block metadata (28 bytes × num_blocks)
    block_metadata: Vec<BlockMeta>,
}

struct BlockMeta {
    num_records: u32,              // Actual records in block (≤ block_size)
    ref_id: i32,                   // Reference sequence ID (-1 = unmapped)
    start_pos: i32,                // First alignment position (0-based)
    end_pos: i32,                  // Last alignment position (0-based)
    compressed_size: u32,          // Block size on disk
    uncompressed_size: u32,        // Block size in memory
    checksum: u32,                 // CRC32 of uncompressed data
}
```

**Index placement**: Written **after** data blocks (similar to BAM/BAI)
**Index offset**: Stored in footer at end of file for fast seeking

### File Footer (32 bytes)

```rust
struct CafFooter {
    index_offset: u64,             // Offset to index section
    num_blocks: u32,               // Total blocks (validation)
    total_records: u64,            // Total alignment records
    checksum: u32,                 // CRC32 of entire file
    magic: [u8; 4],                // "CAFE" (end marker)
    padding: [u8; 4],              // Reserved (zeros)
}
```

**Purpose**: Fast seeking to index without full file scan

---

## Columnar Block Structure

**Block size**: 10,000 records (from OPTIMIZATION_RULES.md Rule 2, Entry 027)
**Rationale**: Validated through 1,440 measurements, balances SIMD benefits vs overhead

### Block Layout

Each block is independently compressed and can be processed in parallel.

```rust
struct CafBlock {
    // Block header (16 bytes, uncompressed)
    block_id: u32,                 // Sequential block number
    num_records: u32,              // Actual records (≤ 10,000)
    uncompressed_size: u32,        // Total size before compression
    checksum: u32,                 // CRC32 for validation

    // Column data (compressed, order matters for efficient decoding)
    columns: ColumnData,
}

struct ColumnData {
    // 1. Core alignment fields (always present)
    ref_ids: CompressedColumn<i32>,      // Reference sequence IDs
    positions: CompressedColumn<i32>,    // 0-based leftmost positions
    mapq: CompressedColumn<u8>,          // Mapping quality [0, 255]
    flags: CompressedColumn<u16>,        // SAM flags (11 bits used)

    // 2. Sequences (pre-decoded for NEON)
    sequences: CompressedColumn<u8>,     // ASCII: 'A','C','G','T','N'
    seq_offsets: CompressedColumn<u32>,  // Cumulative offsets into sequences

    // 3. Quality scores (Phred+33)
    qualities: CompressedColumn<u8>,     // ASCII quality scores (optional)
    qual_offsets: CompressedColumn<u32>, // Same structure as seq_offsets

    // 4. CIGAR operations (optional)
    cigar_ops: CompressedColumn<u32>,    // Standard BAM encoding (op << 28 | len)
    cigar_offsets: CompressedColumn<u32>, // Per-record CIGAR boundaries

    // 5. Read names (optional)
    read_names: CompressedColumn<String>, // Dictionary-compressed

    // 6. Mate information (optional)
    mate_ref_ids: CompressedColumn<i32>,  // Mate reference IDs
    mate_positions: CompressedColumn<i32>, // Mate positions
    template_lengths: CompressedColumn<i32>, // Insert sizes

    // 7. Optional tags (flexible schema)
    tags: TagStorage,                     // Nested columnar
}
```

### Compressed Column Format

Each column is independently compressed:

```rust
struct CompressedColumn<T> {
    compression_type: CompressionType,  // zstd, lz4, rle, or raw
    compressed_len: u32,                // Bytes on disk
    uncompressed_len: u32,              // Bytes in memory
    data: Vec<u8>,                      // Compressed bytes
}

enum CompressionType {
    Raw = 0,         // No compression (e.g., quality scores)
    Zstd = 1,        // zstd level 3 (integers, CIGAR)
    Lz4 = 2,         // lz4 (sequences)
    Rle = 3,         // Run-length encoding (MAPQ)
}
```

**Compression decision logic** (per column, per block):
```rust
fn choose_compression(column_type: ColumnType, data: &[u8]) -> CompressionType {
    match column_type {
        ColumnType::Positions | ColumnType::RefIds | ColumnType::Flags => {
            CompressionType::Zstd  // Sorted integers, good compression
        },
        ColumnType::Sequences => {
            CompressionType::Lz4   // Speed priority for NEON
        },
        ColumnType::Qualities => {
            CompressionType::Raw   // High entropy, skip overhead
        },
        ColumnType::Mapq => {
            // Adaptive: RLE if >50% same value, else raw
            if has_high_repetition(data) {
                CompressionType::Rle
            } else {
                CompressionType::Raw
            }
        },
        _ => CompressionType::Zstd  // Default
    }
}
```

---

## Column Encodings & Compression

### Compression Parameters (Finalized)

| Column | Type | Compression | Level/Params | Rationale |
|--------|------|-------------|--------------|-----------|
| **ref_ids** | `i32` | zstd | level 3 | Sorted, high redundancy [17] |
| **positions** | `i32` | zstd | level 3 | Sorted, delta encoding [20] |
| **mapq** | `u8` | Adaptive RLE | threshold 50% | Often repetitive |
| **flags** | `u16` | zstd | level 3 | Low cardinality patterns |
| **sequences** | `u8` (ASCII) | lz4 | default | Speed priority (>1 GB/s) [19] |
| **seq_offsets** | `u32` | zstd | level 3 | Monotonic increasing |
| **qualities** | `u8` | Raw | none | High entropy, incompressible |
| **cigar_ops** | `u32` | zstd | level 3 | Moderate redundancy |
| **read_names** | `String` | Dictionary | + zstd 3 | Common prefixes |

**zstd parameters**:
- Level 3: Compression ratio ~3-4×, speed ~300-400 MB/s decompression [17]
- Dictionary size: 128 KB (for read names)
- Checksum: CRC32 (builtin)

**lz4 parameters**:
- Mode: Fast (not HC)
- Block size: 64 KB (same as BGZF for comparison)
- Expected ratio: 2-3× on DNA sequences
- Decompression: >1 GB/s single-threaded [19]

### Detailed Column Specifications

### 1. Integers (ref_ids, positions, flags)

**Encoding**: Native little-endian i32/u16
**Compression**: zstd level 3 [17]
**Rationale**:
- Sorted integers enable delta encoding
- zstd level 3 balances ratio vs speed
- Decompression 2-3× faster than gzip [17, 18]

### MAPQ

**Encoding**: u8 [0, 255]
**Compression**: RLE (run-length) or raw
**Decision logic**:
```rust
// If >50% same value, use RLE
if most_common_count > block_size / 2 {
    RleEncoded { value, count }
} else {
    Raw(Vec<u8>)
}
```

### 2. Sequences (Critical for NEON Performance)

**Encoding**: ASCII pre-decoded (uppercase 'A', 'C', 'G', 'T', 'N')
**Compression**: lz4 fast mode [19]
**Storage**: Variable-length array with offset index

**Rationale** (key CAF innovation):
- **BAM approach**: 4-bit encoding (2 bases/byte) requires unpacking
  ```rust
  // BAM decoding: CPU-intensive bit manipulation
  let base = match (packed >> (4 * i)) & 0xF {
      1 => b'A', 2 => b'C', 4 => b'G', 8 => b'T', _ => b'N'
  };
  ```
- **CAF approach**: Pre-decoded ASCII enables direct NEON operations [10, 11]
  ```rust
  // CAF: Zero-copy NEON processing
  let bases = vld1q_u8(sequences.as_ptr());  // Load 16 bases directly
  ```

**Trade-off analysis**:
- Raw storage: 2× larger than 4-bit (1 byte vs 0.5 bytes per base)
- With lz4 compression: 2-3× ratio → **final size ~0.33-0.5 bytes/base**
- **Net effect**: Similar file size, 10-25× faster NEON operations [biometal benchmarks]

**Storage layout**:
```rust
// Sequences stored contiguously
sequences: Vec<u8> = [
    b'A', b'C', b'G', ...,  // Record 0 (150 bases)
    b'T', b'G', b'C', ...,  // Record 1 (150 bases)
    ...
]

// Offset index (one extra element for last boundary)
seq_offsets: Vec<u32> = [0, 150, 300, 450, ...]

// Access record i:
let start = seq_offsets[i] as usize;
let end = seq_offsets[i+1] as usize;
let sequence = &sequences[start..end];
```

**NEON optimization opportunity**:
```rust
#[cfg(target_arch = "aarch64")]
unsafe fn count_bases_neon(sequences: &[u8]) -> BaseCounts {
    // Process 16 bases in parallel (proven 25× speedup in biometal)
    let a_vec = vdupq_n_u8(b'A');
    let c_vec = vdupq_n_u8(b'C');
    let g_vec = vdupq_n_u8(b'G');
    let t_vec = vdupq_n_u8(b'T');

    for chunk in sequences.chunks_exact(16) {
        let bases = vld1q_u8(chunk.as_ptr());
        // Count all 4 bases in parallel...
    }
}
```

### Quality Scores

**Encoding**: ASCII Phred+33
**Compression**: Raw (no compression)
**Rationale**:
- Quality scores are high-entropy (random)
- Compression ratio <1.1×, not worth overhead
- Instant access for NEON filtering

### CIGAR

**Encoding**: Standard BAM u32 (4 bits op, 28 bits len)
**Compression**: RLE for repetitive CIGARs
**Alternative**: Raw if heterogeneous

### Read Names

**Encoding**: Dictionary + indices
**Rationale**: Many reads share common prefixes

```rust
// Example
struct ReadNameDict {
    prefixes: Vec<String>,   // ["SRR123.", "HWI:", ...]
    indices: Vec<u16>,       // [0, 0, 0, 1, ...]
    suffixes: Vec<String>,   // ["1", "2", "3", "1000", ...]
}
```

### Tags (Optional Fields)

**Encoding**: Nested columnar
**Schema**: Flexible, defined in header

```rust
struct TagColumn {
    tag_name: [u8; 2],           // e.g., "NM", "MD"
    tag_type: TagType,           // i32, String, etc.
    values: Vec<TagValue>,        // Type-specific storage
    present_mask: BitVec,         // Which records have this tag
}
```

---

## NEON Optimization Opportunities

### 1. Quality Filtering (Target: 20× speedup)

```rust
#[cfg(target_arch = "aarch64")]
unsafe fn filter_quality_neon(block: &CafBlock, min_q: u8) -> Vec<bool> {
    let threshold = vdupq_n_u8(min_q);
    let mut masks = vec![false; block.qualities.len()];

    for (i, chunk) in block.qualities.chunks_exact(16).enumerate() {
        let quals = vld1q_u8(chunk.as_ptr());
        let cmp = vcgeq_u8(quals, threshold);  // 16 compares in parallel

        // Store result
        vst1q_u8(masks[i * 16..].as_mut_ptr(), cmp);
    }

    masks
}
```

**Expected**: 16-25× faster than scalar (proven in biometal)

### 2. Base Counting (Target: 20× speedup)

```rust
// Reuse biometal's proven NEON kernel
// Process 16 ASCII bases in parallel
// Count A/C/G/T/N frequencies
```

### 3. MAPQ Filtering (Target: 16× speedup)

```rust
// Direct access to mapq array (no unpacking)
// Parallel threshold comparison
for chunk in block.mapq.chunks(16) {
    let mapqs = vld1q_u8(chunk.as_ptr());
    let mask = vcgtq_u8(mapqs, threshold);
    // 16 comparisons in one instruction
}
```

---

## Index Format (.caf.idx)

**Purpose**: Block-level random access (chr:start-end queries)

**Structure**:
```
[Magic: "CAFI" (4 bytes)]
[Version: u16]
[Num blocks: u32]
[Block offsets: [u64; num_blocks]]
[Block metadata: [BlockMeta; num_blocks]]
```

**Query algorithm**:
```rust
fn query_region(idx: &CafIndex, chr: &str, start: i32, end: i32) -> Vec<usize> {
    idx.block_metadata
        .iter()
        .enumerate()
        .filter(|(_, meta)| {
            meta.ref_id == chr_id &&
            meta.start_pos <= end &&
            meta.end_pos >= start
        })
        .map(|(i, _)| i)
        .collect()
}
```

**Granularity**: Block-level (10,000 records)
- Coarser than BAM's BAI (64 KB chunks)
- Sufficient for analytical queries
- Not suitable for single-base genome browsers

---

## Lossless Conversion: BAM ↔ CAF

### BAM → CAF

1. Read BAM header → Store in CAF header (compressed)
2. Accumulate 10,000 records
3. Convert row-oriented → columnar:
   - Extract positions, mapq, flags to arrays
   - Decode 4-bit sequences → ASCII
   - Extract quality scores
   - Parse CIGAR
   - Extract tags
4. Compress columns
5. Write CAF block
6. Build index

### CAF → BAM

1. Read CAF header → Reconstruct SAM header
2. For each block:
   - Decompress columns
   - For each record i:
     - Reconstruct BAM record from columns[i]
     - Encode sequences ASCII → 4-bit
     - Write BAM record
3. Write using standard BAM writer

### Validation

**Property**: `bam_to_caf(caf_to_bam(x)) == x`

**Testing**:
- Differential testing (1,000+ diverse BAM files)
- Property-based testing (proptest)
- Byte-level comparison

---

## Storage Trade-offs

### Expected Sizes (100K records)

| Format | Size | Ratio |
|--------|------|-------|
| SAM (uncompressed) | 500 MB | 10.0× |
| BAM (gzip) | 50 MB | 1.0× (baseline) |
| CAF (zstd/lz4) | 75-100 MB | 1.5-2.0× |

**Breakdown by column**:
- Positions: ~1-2 bytes (compressed) ✅ Similar to BAM
- Sequences: ~0.5 bytes (lz4) ✅ Similar to BAM (4-bit)
- Qualities: ~1 byte (raw) ❌ Slightly worse than BAM (0.9)
- Metadata: ~0.5 bytes ✅ Similar to BAM

**Trade-off acceptance**: 1.5-2× larger for 5-10× faster analytical operations

---

## Performance Targets

| Operation | BAM | CAF (NEON) | Target Speedup |
|-----------|-----|------------|----------------|
| Quality filter Q30 | 2.00 sec | 0.08 sec | **25×** |
| Base counting | 1.50 sec | 0.06 sec | **25×** |
| MAPQ > 30 filter | 0.50 sec | 0.03 sec | **16×** |
| Parse 100K records | 2.56 sec | 0.25 sec | **10×** |
| **Overall** | **1.0×** | **5-10×** | **Target** |

**Validation**: Statistical tests (N=30, p < 0.05)

---

## Use Cases

### ✅ Recommended For

1. **Batch processing pipelines**
   - Filter, aggregate, transform
   - Constant memory streaming
   - Parallel block processing

2. **ML data preparation**
   - Extract features from alignments
   - Quality filtering
   - Fast iteration

3. **Quality control**
   - Coverage analysis
   - Mapping quality assessment
   - Base composition checks

4. **Research workflows**
   - Experimental analysis
   - Custom operations
   - Performance-critical paths

### ❌ Not Recommended For

1. **Genome browsers** (IGV, etc.)
   - Requires base-level random access
   - CAF is block-level (10K records)

2. **Long-term archival**
   - BAM is universal standard
   - CAF larger (1.5-2×)

3. **Clinical pipelines**
   - Regulatory compliance requires BAM
   - CAF is experimental

4. **Immediate sharing**
   - Recipients may not have CAF readers
   - Convert to BAM for distribution

---

## Future Extensions

### Version 1.1 (Potential)

**GPU optimization**:
- Metal/CUDA-friendly memory layout
- Batch kernel launches
- Unified memory for sequences

**Neural Engine**:
- Quantized integer operations
- Base counting on Neural Engine
- Quality filtering acceleration

**Advanced indexing**:
- Finer-grained than blocks
- Spatial indexing (R-tree)
- Predicate pushdown

---

## Open Research Questions

1. **Optimal block size**: 10K validated, but test 1K, 100K?
2. **Compression tuning**: Which zstd level? When to use lz4 vs raw?
3. **Tag encoding**: Fixed schema vs flexible?
4. **GPU integration**: Phase 1 or deferred?

---

## Error Handling and Validation

### File Validation Requirements

**Mandatory checks on read**:
1. **Magic number**: First 4 bytes must be `CAF\x01`
2. **Version compatibility**: Major version must be supported (v1.x)
3. **Checksums**: Validate CRC32 for each block and full file
4. **Index consistency**: `num_blocks` in header must match index
5. **Block metadata**: `uncompressed_size` must match actual decompressed size

### Error Types

```rust
#[derive(Debug, Error)]
pub enum CafError {
    #[error("Invalid magic number: expected CAF\\x01, got {0:?}")]
    InvalidMagic([u8; 4]),

    #[error("Unsupported version: {major}.{minor}")]
    UnsupportedVersion { major: u8, minor: u8 },

    #[error("Checksum mismatch in block {block_id}: expected {expected:08x}, got {actual:08x}")]
    ChecksumMismatch {
        block_id: u32,
        expected: u32,
        actual: u32,
    },

    #[error("Decompression failed for {column} in block {block_id}: {source}")]
    DecompressionError {
        column: String,
        block_id: u32,
        source: Box<dyn std::error::Error>,
    },

    #[error("Index inconsistency: header claims {header_blocks} blocks, index has {index_blocks}")]
    IndexInconsistency {
        header_blocks: u32,
        index_blocks: u32,
    },

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}
```

### Corruption Detection

**CRC32 checksums** at multiple levels:
1. **Per-block**: Checksum of uncompressed column data
2. **Full-file**: Checksum in footer covers entire file
3. **Validation**: Fail fast on first checksum mismatch

**Recovery strategies**:
- Block-level corruption: Skip corrupted block, continue reading
- Header corruption: Fatal error, file unreadable
- Index corruption: Rebuild index from data blocks (slower but recoverable)

---

## Conformance and Testing

### Validation Suite Requirements

**Round-trip testing**:
```rust
fn test_lossless_conversion(bam_path: &Path) {
    let bam = BamReader::from_path(bam_path).unwrap();
    let caf = bam_to_caf(&bam).unwrap();
    let bam2 = caf_to_bam(&caf).unwrap();

    // Byte-level comparison
    assert_eq!(bam.records(), bam2.records());
}
```

**Property-based testing** [27]:
```rust
proptest! {
    #[test]
    fn test_block_roundtrip(records in vec(arbitrary_bam_record(), 1..10000)) {
        let block = records_to_caf_block(&records);
        let recovered = caf_block_to_records(&block);
        prop_assert_eq!(records, recovered);
    }
}
```

**Benchmark validation** [26]:
- N = 30 iterations per operation
- Statistical significance: p < 0.05 (paired t-test)
- Effect size: Cohen's d > 0.8 (large effect)
- Reproducibility: Fixed datasets, public code

---

## Implementation Recommendations

### Suggested Development Order

**Phase 1: Core functionality** (Weeks 2-3)
1. CAF writer: BAM → CAF conversion
2. CAF reader: CAF → BAM conversion
3. Lossless round-trip validation

**Phase 2: NEON optimization** (Weeks 4-5)
4. NEON quality filtering
5. NEON base counting
6. NEON MAPQ filtering
7. Benchmark vs BAM (N=30)

**Phase 3: Robustness** (Week 6)
8. Error handling (checksums, validation)
9. Property-based testing
10. Cross-platform validation

**Phase 4: Documentation** (Weeks 7-8)
11. API documentation
12. Usage examples
13. Manuscript writing

### Reference Implementation

**Language**: Rust 1.70+ (2021 edition)
**Dependencies**:
- `zstd = "0.13"` - Modern compression [17]
- `lz4 = "1.24"` - Fast decompression [19]
- `bincode = "1.3"` - Binary serialization
- `crc32fast = "1.3"` - Checksum validation

**Traits to implement**:
```rust
pub trait CafReader {
    fn read_block(&mut self, block_id: usize) -> Result<CafBlock, CafError>;
    fn query_region(&self, chr: &str, start: i32, end: i32) -> Result<Vec<CafBlock>, CafError>;
}

pub trait CafWriter {
    fn write_block(&mut self, block: &CafBlock) -> Result<(), CafError>;
    fn finalize(&mut self) -> Result<(), CafError>;  // Write index + footer
}
```

---

## References

Complete references available in `research/caf-format/LITERATURE_REVIEW.md` (28 citations).

**Key citations**:
- **[1]** Li H, et al. (2009). BAM format specification. *Bioinformatics*, 25(16), 2078-2079
- **[3]** Fritz MH, et al. (2011). CRAM reference-based compression. *Genome Research*, 21(5), 734-740
- **[6-7]** Melnik et al. (2010). Columnar analytics (Dremel). *VLDB*, 3(1-2), 330-339
- **[10-11]** ARM NEON programming guides
- **[12]** Zhao et al. (2013). SSW SIMD alignment (16-25× speedup). *PLOS ONE*, 8(12), e82138
- **[17]** Collet & Kucherawy (2021). zstd RFC 8878. *IETF*
- **[19]** Collet (2011-2023). LZ4 specification
- **[20]** Abadi et al. (2008). Column-stores vs row-stores. *SIGMOD*, 967-980
- **[21]** Hennessy & Patterson (2019). Computer Architecture (6th ed.)
- **[26-27]** Hoefler & Belli (2015). Scientific benchmarking. Mytkowicz et al. (2009). Measurement pitfalls

See `LITERATURE_REVIEW.md` for complete bibliography.

---

## Appendices

### Appendix A: Size Estimation Calculator

```python
def estimate_caf_size(num_records, avg_seq_len=150, avg_qual=True):
    """Estimate CAF file size given dataset characteristics"""
    # Positions, flags, mapq, mate info: ~12 bytes compressed
    metadata = num_records * 12

    # Sequences: ASCII (1 byte) compressed 2.5× with lz4
    sequences = (num_records * avg_seq_len) / 2.5

    # Qualities: raw (incompressible)
    qualities = num_records * avg_seq_len if avg_qual else 0

    # CIGAR, read names: ~6 bytes compressed
    other = num_records * 6

    total = metadata + sequences + qualities + other
    return total
```

**Example**: 1M records, 150bp average:
- Metadata: 1M × 12 = 12 MB
- Sequences: 1M × 150 / 2.5 = 60 MB
- Qualities: 1M × 150 = 150 MB
- Other: 1M × 6 = 6 MB
- **Total: ~228 MB** (vs BAM ~150 MB = 1.5× larger)

### Appendix B: Comparison Matrix

| Feature | BAM | CRAM | CAF |
|---------|-----|------|-----|
| **Year** | 2009 | 2011 | 2025 |
| **Layout** | Row | Row | Columnar |
| **Sequence encoding** | 4-bit | Reference | ASCII |
| **Compression** | gzip | rANS + ref | zstd/lz4 |
| **SIMD support** | Limited | Limited | Extensive (NEON) |
| **Random access** | BAI (fine) | CRAI | Block-level |
| **Size vs BAM** | 1.0× | 0.4-0.7× | 1.5-2.0× |
| **Speed (analytical)** | 1.0× | 0.3-0.5× | 5-10× (target) |
| **Best for** | General | Archival | Analytics, ARM |
| **Lossless** | Yes | Yes | Yes (from BAM) |

---

**Specification Status**: ✅ **Finalized v1.0.0** (November 10, 2025)
**Implementation Status**: Not started (begins Week 2)
**Validation**: Pending experimental results (Weeks 4-6)
**Next Review**: After Phase 1 implementation

---

**For questions or suggestions, see**: `research/caf-format/README.md`
**Project repository**: biometal @ github.com/scotthandley/biometal
