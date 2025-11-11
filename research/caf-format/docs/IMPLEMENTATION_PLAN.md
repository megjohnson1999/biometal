# CAF Implementation Plan

**Status**: Planning Phase
**Timeline**: Weeks 2-6 (Nov 18 - Dec 22, 2025)
**Language**: Rust 1.70+ (2021 edition)
**Target Platform**: ARM (macOS, Graviton) with x86_64 fallback

---

## Overview

This document outlines the concrete implementation plan for the CAF (Columnar Alignment Format) library, following the finalized v1.0.0 specification.

**Implementation Goals**:
1. ✅ **Correctness**: 100% lossless BAM ↔ CAF conversion
2. ✅ **Performance**: 5-10× speedup over BAM for analytical operations
3. ✅ **Robustness**: Production-quality error handling and validation
4. ✅ **Evidence-based**: Follow biometal's OPTIMIZATION_RULES.md

---

## Module Architecture

```
research/caf-format/implementation/
├── src/
│   ├── lib.rs                    # Public API
│   ├── error.rs                  # Error types (CafError)
│   ├── types.rs                  # Core data structures
│   │
│   ├── format/                   # Binary format handling
│   │   ├── mod.rs
│   │   ├── header.rs            # CAF header parsing/writing
│   │   ├── index.rs             # Index structure
│   │   ├── footer.rs            # Footer structure
│   │   └── magic.rs             # Magic number validation
│   │
│   ├── block/                    # Columnar block operations
│   │   ├── mod.rs
│   │   ├── builder.rs           # Build blocks from records
│   │   ├── reader.rs            # Read/decompress blocks
│   │   ├── writer.rs            # Compress/write blocks
│   │   └── schema.rs            # Column schema definitions
│   │
│   ├── column/                   # Column-specific encoding
│   │   ├── mod.rs
│   │   ├── integers.rs          # Positions, ref_ids, flags
│   │   ├── sequences.rs         # ASCII sequences
│   │   ├── qualities.rs         # Quality scores
│   │   ├── cigar.rs             # CIGAR operations
│   │   ├── read_names.rs        # Dictionary encoding
│   │   └── tags.rs              # Auxiliary tags
│   │
│   ├── compression/              # Compression strategies
│   │   ├── mod.rs
│   │   ├── zstd.rs              # zstd level 3
│   │   ├── lz4.rs               # lz4 fast
│   │   ├── rle.rs               # Run-length encoding
│   │   └── adaptive.rs          # Adaptive selection
│   │
│   ├── conversion/               # BAM ↔ CAF conversion
│   │   ├── mod.rs
│   │   ├── bam_to_caf.rs        # BAM → CAF converter
│   │   ├── caf_to_bam.rs        # CAF → BAM converter
│   │   └── roundtrip.rs         # Validation helpers
│   │
│   ├── io/                       # File I/O operations
│   │   ├── mod.rs
│   │   ├── reader.rs            # CafReader trait
│   │   ├── writer.rs            # CafWriter trait
│   │   └── buffered.rs          # Buffered I/O
│   │
│   ├── query/                    # Query operations
│   │   ├── mod.rs
│   │   ├── region.rs            # chr:start-end queries
│   │   └── filter.rs            # Predicate pushdown
│   │
│   ├── neon/                     # ARM NEON optimizations
│   │   ├── mod.rs
│   │   ├── base_counting.rs    # NEON base counting
│   │   ├── quality_filter.rs   # NEON quality filtering
│   │   ├── mapq_filter.rs      # NEON MAPQ filtering
│   │   └── scalar.rs            # Scalar fallback
│   │
│   ├── validation/               # Validation and checksums
│   │   ├── mod.rs
│   │   ├── checksum.rs          # CRC32 validation
│   │   └── roundtrip.rs         # Lossless testing
│   │
│   └── bin/                      # Command-line tools
│       ├── caf-convert.rs       # BAM ↔ CAF conversion
│       ├── caf-query.rs         # Query CAF files
│       └── caf-stats.rs         # Statistics
│
├── tests/                        # Integration tests
│   ├── roundtrip_tests.rs       # Lossless conversion
│   ├── neon_tests.rs            # NEON correctness
│   └── property_tests.rs        # Property-based (proptest)
│
├── benches/                      # Benchmarks
│   ├── caf_vs_bam.rs            # Performance comparison
│   ├── neon_operations.rs       # NEON speedup validation
│   └── compression.rs           # Compression benchmarks
│
└── examples/                     # Usage examples
    ├── basic_conversion.rs
    ├── neon_filtering.rs
    └── query_region.rs
```

---

## Phase 1: Core Conversion (Weeks 2-3)

**Goal**: Implement lossless BAM ↔ CAF conversion

### Week 2: Foundation (Nov 18-24)

**Day 1-2: Project setup + Core types**
- [ ] Create Cargo workspace structure
- [ ] Define core types in `types.rs`:
  ```rust
  pub struct CafHeader { ... }
  pub struct CafBlock { ... }
  pub struct CafIndex { ... }
  pub struct CafFooter { ... }
  pub struct BlockMeta { ... }
  ```
- [ ] Implement error types in `error.rs` (CafError enum)
- [ ] Set up dependencies (zstd, lz4, bincode, crc32fast)

**Day 3-4: Format parsing**
- [ ] Implement `format/magic.rs` - magic number validation
- [ ] Implement `format/header.rs` - header parsing/writing
- [ ] Implement `format/index.rs` - index structure
- [ ] Implement `format/footer.rs` - footer structure
- [ ] Add unit tests for each format component

**Day 5-7: Column encoding (basic)**
- [ ] Implement `column/integers.rs` - positions, ref_ids, flags
- [ ] Implement `column/sequences.rs` - ASCII sequences with offsets
- [ ] Implement `column/qualities.rs` - quality scores
- [ ] Implement `column/cigar.rs` - CIGAR operations
- [ ] Add column round-trip tests

### Week 3: Conversion Pipeline (Nov 25 - Dec 1)

**Day 1-2: Compression layer**
- [ ] Implement `compression/zstd.rs` - zstd level 3
- [ ] Implement `compression/lz4.rs` - lz4 fast mode
- [ ] Implement `compression/rle.rs` - run-length encoding
- [ ] Implement `compression/adaptive.rs` - adaptive selection
- [ ] Benchmark compression ratios on sample data

**Day 3-4: Block builder**
- [ ] Implement `block/builder.rs` - accumulate 10K records
- [ ] Implement `block/writer.rs` - compress and write blocks
- [ ] Implement `block/reader.rs` - read and decompress blocks
- [ ] Test block round-trip (records → block → records)

**Day 5-7: BAM conversion**
- [ ] Implement `conversion/bam_to_caf.rs`:
  ```rust
  pub fn bam_to_caf<R, W>(bam: &mut BamReader<R>, caf: &mut CafWriter<W>) -> Result<()>
  ```
- [ ] Implement `conversion/caf_to_bam.rs`:
  ```rust
  pub fn caf_to_bam<R, W>(caf: &mut CafReader<R>, bam: &mut BamWriter<W>) -> Result<()>
  ```
- [ ] **CRITICAL**: Lossless round-trip validation on 10+ diverse BAM files
- [ ] Document any edge cases or special handling

**Week 3 Milestone**: ✅ Lossless BAM ↔ CAF conversion working

---

## Phase 2: NEON Optimization (Weeks 4-5)

**Goal**: Implement and validate ARM NEON speedups

### Week 4: NEON Kernels (Dec 2-8)

**Day 1-2: Base counting**
- [ ] Implement `neon/base_counting.rs`:
  ```rust
  #[cfg(target_arch = "aarch64")]
  unsafe fn count_bases_neon(sequences: &[u8]) -> BaseCounts

  #[cfg(not(target_arch = "aarch64"))]
  fn count_bases_scalar(sequences: &[u8]) -> BaseCounts
  ```
- [ ] Validate correctness: NEON vs scalar (proptest)
- [ ] Benchmark: Target 25× speedup (N=30, p<0.05)

**Day 3-4: Quality filtering**
- [ ] Implement `neon/quality_filter.rs`:
  ```rust
  unsafe fn filter_quality_neon(block: &CafBlock, min_q: u8) -> Vec<bool>
  ```
- [ ] 16 quality scores in parallel (proven in biometal)
- [ ] Benchmark: Target 25× speedup

**Day 5-7: MAPQ filtering**
- [ ] Implement `neon/mapq_filter.rs`:
  ```rust
  unsafe fn filter_mapq_neon(block: &CafBlock, min_mapq: u8) -> Vec<bool>
  ```
- [ ] Direct access to mapq column (no unpacking)
- [ ] Benchmark: Target 16× speedup

### Week 5: Benchmarking and Optimization (Dec 9-15)

**Day 1-2: Comprehensive benchmarks**
- [ ] Implement `benches/caf_vs_bam.rs`:
  - Quality filtering (Q30)
  - Base counting
  - MAPQ filtering (>30)
  - Full record parsing (100K records)
- [ ] Run on 5 diverse datasets (N=30 each)
- [ ] Statistical analysis: paired t-test, p-values, Cohen's d

**Day 3-4: Performance profiling**
- [ ] Profile with `cargo flamegraph`
- [ ] Identify bottlenecks
- [ ] Optimize hot paths
- [ ] Re-benchmark after optimizations

**Day 5-7: Compression tuning**
- [ ] Measure actual compression ratios on real data
- [ ] Validate 1.5-2× size overhead vs BAM
- [ ] Tune zstd dictionary for read names
- [ ] Document compression performance

**Week 5 Milestone**: ✅ 5-10× overall speedup validated (N=30, p<0.05)

---

## Phase 3: Robustness (Week 6)

**Goal**: Production-quality error handling and validation

### Week 6: Validation and Testing (Dec 16-22)

**Day 1-2: Error handling**
- [ ] Implement comprehensive `CafError` handling
- [ ] Add CRC32 checksums at all levels
- [ ] Test corruption detection (inject errors)
- [ ] Implement recovery strategies (skip corrupted blocks)

**Day 3-4: Property-based testing**
- [ ] Implement `tests/property_tests.rs`:
  ```rust
  proptest! {
      fn roundtrip_arbitrary_records(records: Vec<BamRecord>)
      fn block_compression_preserves_data(block: CafBlock)
      fn neon_matches_scalar(sequences: Vec<u8>)
  }
  ```
- [ ] Run 10,000+ random test cases
- [ ] Fix any discovered edge cases

**Day 5-7: Cross-platform validation**
- [ ] Test on macOS ARM (M1/M2/M3/M4)
- [ ] Test on Linux ARM (Graviton if available)
- [ ] Test on x86_64 (scalar fallback)
- [ ] Validate NEON vs scalar correctness
- [ ] Document platform-specific behavior

**Week 6 Milestone**: ✅ Production-ready, fully validated implementation

---

## API Design

### Public API (lib.rs)

```rust
// Re-exports
pub use error::{CafError, Result};
pub use types::{CafHeader, CafBlock, CafIndex, BlockMeta};
pub use io::{CafReader, CafWriter};
pub use conversion::{bam_to_caf, caf_to_bam};

// Traits
pub trait CafReader {
    fn read_header(&mut self) -> Result<CafHeader>;
    fn read_block(&mut self, block_id: usize) -> Result<CafBlock>;
    fn query_region(&self, chr: &str, start: i32, end: i32) -> Result<Vec<CafBlock>>;
}

pub trait CafWriter {
    fn write_header(&mut self, header: &CafHeader) -> Result<()>;
    fn write_block(&mut self, block: &CafBlock) -> Result<()>;
    fn finalize(&mut self) -> Result<()>;  // Write index + footer
}

// NEON operations (feature-gated)
#[cfg(feature = "neon")]
pub mod neon {
    pub fn count_bases(sequences: &[u8]) -> BaseCounts;
    pub fn filter_quality(block: &CafBlock, min_q: u8) -> Vec<bool>;
    pub fn filter_mapq(block: &CafBlock, min_mapq: u8) -> Vec<bool>;
}

// Conversion utilities
pub fn convert_file(input: &Path, output: &Path, format: Format) -> Result<()>;
```

### Command-Line Tools

**caf-convert**:
```bash
# BAM → CAF
caf-convert input.bam output.caf

# CAF → BAM
caf-convert input.caf output.bam

# With options
caf-convert input.bam output.caf --block-size 10000 --compression-level 3
```

**caf-query**:
```bash
# Query region
caf-query file.caf chr1:1000-2000

# Filter by quality
caf-query file.caf --min-quality 30

# Filter by MAPQ
caf-query file.caf --min-mapq 30 --output filtered.caf
```

**caf-stats**:
```bash
# Show file statistics
caf-stats file.caf

# Output:
#   File size: 228 MB
#   Blocks: 100
#   Records: 1,000,000
#   Compression ratio: 1.6×
#   Avg quality: 35.2
```

---

## Testing Strategy

### Unit Tests
- **Target**: 80%+ code coverage
- **Location**: Inline `#[cfg(test)]` modules
- **Focus**: Individual functions and modules

### Integration Tests
- **Location**: `tests/` directory
- **Coverage**:
  - Round-trip conversion (10+ BAM files)
  - Query operations
  - Error handling
  - Checksum validation

### Property-Based Tests
- **Tool**: proptest
- **Coverage**:
  - Round-trip: `BAM → CAF → BAM` preserves data
  - NEON correctness: NEON == scalar
  - Block compression: lossless
  - Random data: no panics

### Benchmarks
- **Tool**: Criterion (N=30 by default)
- **Datasets**: 5 diverse BAM files
- **Operations**:
  - Quality filtering (Q30)
  - Base counting
  - MAPQ filtering (>30)
  - Full parsing (100K records)
- **Statistics**: Mean, 95% CI, p-value (paired t-test)

---

## Dependencies

```toml
[dependencies]
# Compression
zstd = "0.13"
lz4 = "1.24"

# Serialization
bincode = "1.3"
serde = { version = "1.0", features = ["derive"] }

# Error handling
thiserror = "1.0"
anyhow = "1.0"

# Checksums
crc32fast = "1.3"

# Logging
log = "0.4"
env_logger = "0.11"

# CLI
clap = { version = "4.4", features = ["derive"] }

# BAM integration
biometal = { path = "../../../" }

# Optional: NEON intrinsics (ARM only)
[target.'cfg(target_arch = "aarch64")'.dependencies]
# NEON intrinsics available via std::arch

[dev-dependencies]
# Testing
proptest = "1.4"
quickcheck = "1.0"
criterion = { version = "0.5", features = ["html_reports"] }
tempfile = "3.8"

# Benchmarking
csv = "1.3"
```

---

## Development Guidelines

### Code Quality
- **No `unwrap()` or `expect()`** in library code (use `Result`)
- **Document all public APIs** with examples
- **Follow biometal patterns** (error handling, NEON structure)
- **Rust 2021 idioms** (edition = "2021")

### NEON Safety
- **Always provide scalar fallback**:
  ```rust
  pub fn operation(data: &[u8]) -> Result {
      #[cfg(target_arch = "aarch64")]
      { unsafe { operation_neon(data) } }

      #[cfg(not(target_arch = "aarch64"))]
      { operation_scalar(data) }
  }
  ```
- **Validate NEON == scalar** (proptest)
- **Document safety invariants** for `unsafe` code

### Performance
- **Follow OPTIMIZATION_RULES.md**:
  - Rule 1: NEON for element-wise ops
  - Rule 2: 10K block size
  - Rule 3: Parallel decompression (future)
  - Rule 5: Constant memory streaming
- **Benchmark before optimizing** (N=30, p<0.05)
- **Profile hot paths** (flamegraph)

### Error Messages
- **Be specific**:
  ```rust
  // Good
  CafError::ChecksumMismatch {
      block_id: 42,
      expected: 0xDEADBEEF,
      actual: 0xBADC0FFE,
  }

  // Bad
  "Checksum error"
  ```
- **Include context** (file offset, block ID, column name)
- **Suggest fixes** when possible

---

## Milestones and Deliverables

### Week 2 (Nov 18-24)
- ✅ Project structure set up
- ✅ Core types defined
- ✅ Format parsing working
- ✅ Basic column encoding

### Week 3 (Nov 25 - Dec 1)
- ✅ Compression layer working
- ✅ Block builder/reader/writer
- ✅ **CRITICAL**: Lossless BAM ↔ CAF conversion
- ✅ 10+ diverse BAM files validated

### Week 4 (Dec 2-8)
- ✅ NEON base counting (25× target)
- ✅ NEON quality filtering (25× target)
- ✅ NEON MAPQ filtering (16× target)

### Week 5 (Dec 9-15)
- ✅ Comprehensive benchmarks (N=30)
- ✅ Statistical analysis (p<0.05, Cohen's d>0.8)
- ✅ **GOAL**: 5-10× overall speedup validated
- ✅ Compression ratios measured (1.5-2× overhead)

### Week 6 (Dec 16-22)
- ✅ Error handling and checksums
- ✅ Property-based testing (10K+ cases)
- ✅ Cross-platform validation
- ✅ **Production-ready** implementation

---

## Risk Mitigation

### Risk 1: Lossless Conversion Failures
- **Mitigation**: Extensive testing on diverse BAM files
- **Fallback**: Document known limitations, fail gracefully
- **Timeline impact**: High (blocks publication if not achieved)

### Risk 2: NEON Speedup Below Target
- **Mitigation**: Follow biometal's proven NEON patterns
- **Fallback**: Document actual speedup, adjust claims
- **Timeline impact**: Medium (reduces impact, not fatal)

### Risk 3: Compression Overhead >2×
- **Mitigation**: Adaptive compression, tune parameters
- **Fallback**: Accept higher overhead, emphasize speed
- **Timeline impact**: Low (storage cheap, speed priority)

### Risk 4: Implementation Complexity
- **Mitigation**: Incremental development, weekly milestones
- **Fallback**: Reduce scope (defer GPU, tags, advanced indexing)
- **Timeline impact**: Medium (can extend to 10-12 weeks)

---

## Success Criteria

**Must Have** (for publication):
1. ✅ Lossless BAM ↔ CAF conversion (100% on 1,000+ files)
2. ✅ Statistical validation (N=30, p<0.05)
3. ✅ 5-10× overall speedup demonstrated
4. ✅ Open-source implementation (MIT license)

**Should Have** (for quality):
5. ✅ Property-based testing (10K+ cases)
6. ✅ Cross-platform validation (ARM + x86_64)
7. ✅ Error handling and checksums
8. ✅ Command-line tools

**Nice to Have** (deferred to v1.1):
9. ⏳ GPU optimization (Metal/CUDA)
10. ⏳ Parallel decompression (Rule 3)
11. ⏳ Advanced indexing (finer than blocks)

---

**Status**: Planning complete, ready for implementation
**Next Steps**: Begin Week 2 (Nov 18) - Project setup + core types
**Review**: Weekly progress check against milestones
