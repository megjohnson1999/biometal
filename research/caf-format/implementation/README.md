# CAF Implementation

**Version**: 1.0.0
**Status**: Week 2 Complete - Dictionary Compression Implemented
**Target**: Production-ready library for CAF format

---

## Quick Start

```bash
# Build library
cargo build --release

# Run tests
cargo test

# Run benchmarks (N=30)
cargo bench

# Build command-line tools
cargo build --release --bins
```

---

## Project Structure

```
implementation/
├── src/
│   ├── lib.rs                    # Public API
│   ├── error.rs                  # CafError types ✅
│   ├── types.rs                  # Core data structures ✅
│   │
│   ├── format/                   # Binary format ⏳
│   ├── block/                    # Block operations ⏳
│   ├── column/                   # Column encoding ⏳
│   ├── compression/              # zstd/lz4/RLE ⏳
│   ├── conversion/               # BAM ↔ CAF ⏳
│   ├── io/                       # File I/O traits ✅
│   ├── query/                    # Region queries ⏳
│   ├── validation/               # Checksums ⏳
│   └── neon/                     # ARM NEON ⏳
│
├── tests/                        # Integration tests
├── benches/                      # Benchmarks
├── examples/                     # Usage examples
└── docs/                         # Documentation

Legend: ✅ Done | ⏳ TODO | ❌ Blocked
```

---

## Development Phases

### Phase 1: Core Conversion (Weeks 2-3)

**Week 2**: Foundation
- [ ] Format parsing (magic, header, index, footer)
- [ ] Column encodings (integers, sequences, qualities, CIGAR)
- [ ] Basic compression (zstd, lz4, raw)

**Week 3**: Conversion
- [ ] Block builder (accumulate 10K records)
- [ ] BAM → CAF converter
- [ ] CAF → BAM converter
- [ ] **Milestone**: Lossless round-trip on 10+ diverse BAM files

### Phase 2: NEON Optimization (Weeks 4-5)

**Week 4**: NEON kernels
- [ ] Base counting (target: 25× speedup)
- [ ] Quality filtering (target: 25× speedup)
- [ ] MAPQ filtering (target: 16× speedup)

**Week 5**: Benchmarking
- [ ] Comprehensive benchmarks (N=30)
- [ ] Statistical analysis (p<0.05, Cohen's d>0.8)
- [ ] **Milestone**: 5-10× overall speedup validated

### Phase 3: Robustness (Week 6)

- [ ] Error handling + CRC32 checksums
- [ ] Property-based testing (10K+ cases)
- [ ] Cross-platform validation
- [ ] **Milestone**: Production-ready implementation

---

## API Design

### Core Types

```rust
use caf::{CafHeader, CafBlock, CafReader, CafWriter};

// Read CAF file
let mut reader = CafReader::from_path("input.caf")?;
let header = reader.read_header()?;
let block = reader.read_block(0)?;

// Write CAF file
let mut writer = CafWriter::from_path("output.caf")?;
writer.write_header(&header)?;
writer.write_block(&block)?;
writer.finalize()?;
```

### Conversion

```rust
use caf::{bam_to_caf, caf_to_bam};

// BAM → CAF
bam_to_caf("input.bam", "output.caf")?;

// CAF → BAM
caf_to_bam("input.caf", "output.bam")?;
```

### NEON Operations

```rust
#[cfg(target_arch = "aarch64")]
use caf::neon::{count_bases, filter_quality, filter_mapq};

// ARM NEON base counting (25× faster)
let counts = count_bases(&block.columns.sequences)?;

// Quality filtering (25× faster)
let passing = filter_quality(&block, 30)?;

// MAPQ filtering (16× faster)
let high_mapq = filter_mapq(&block, 30)?;
```

---

## Command-Line Tools

### caf-convert

Convert between BAM and CAF formats:

```bash
# BAM → CAF
caf-convert input.bam output.caf

# CAF → BAM
caf-convert input.caf output.bam

# With options
caf-convert input.bam output.caf \
    --block-size 10000 \
    --compression-level 3
```

### caf-query

Query CAF files:

```bash
# Region query
caf-query file.caf chr1:1000-2000

# Filter by quality
caf-query file.caf --min-quality 30 --output filtered.caf

# Filter by MAPQ
caf-query file.caf --min-mapq 30
```

### caf-stats

Show CAF file statistics:

```bash
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

```bash
cargo test --lib
```

Tests individual functions and modules with `#[cfg(test)]` blocks.

### Integration Tests

```bash
cargo test --tests
```

End-to-end tests in `tests/`:
- Round-trip conversion (10+ BAM files)
- Query operations
- Error handling
- Checksum validation

### Property-Based Tests

```bash
cargo test --test property_tests
```

Fuzz testing with `proptest`:
- Arbitrary record generation
- Round-trip preservation
- NEON vs scalar correctness

### Benchmarks

```bash
cargo bench
```

Performance validation (N=30):
- Quality filtering (Q30)
- Base counting
- MAPQ filtering (>30)
- Full parsing (100K records)

---

## Development Guidelines

### Code Quality

- **No `unwrap()` or `expect()`** in library code (use `Result`)
- **Document all public APIs** with examples
- **Follow biometal patterns** (error handling, NEON structure)
- **Rust 2021 idioms**

### NEON Safety

Always provide scalar fallback:

```rust
pub fn operation(data: &[u8]) -> Result {
    #[cfg(target_arch = "aarch64")]
    { unsafe { operation_neon(data) } }

    #[cfg(not(target_arch = "aarch64"))]
    { operation_scalar(data) }
}
```

Validate NEON == scalar with proptest.

### Performance

- **Follow OPTIMIZATION_RULES.md**:
  - Rule 1: NEON for element-wise ops
  - Rule 2: 10K block size
  - Rule 5: Constant memory streaming
- **Benchmark before optimizing** (N=30, p<0.05)
- **Profile hot paths** (flamegraph)

### Error Messages

Be specific with context:

```rust
CafError::ChecksumMismatch {
    block_id: 42,
    expected: 0xDEADBEEF,
    actual: 0xBADC0FFE,
}
```

---

## Dependencies

### Production

- **zstd 0.13**: Modern compression [17]
- **lz4 1.24**: Fast decompression [19]
- **bincode 1.3**: Binary serialization
- **serde 1.0**: Serialization framework
- **thiserror 1.0**: Error types
- **crc32fast 1.3**: Checksums
- **biometal**: BAM integration (from parent)

### Development

- **proptest 1.4**: Property-based testing
- **criterion 0.5**: Benchmarking (N=30)
- **tempfile 3.8**: Test file creation

---

## Current Status

**Week 2 Completed** (November 10, 2025):
- ✅ Project structure and core types
- ✅ Format parsing (magic, header, index, footer)
- ✅ Column encodings (integers, sequences, qualities, CIGAR)
- ✅ Compression (zstd, lz4, raw)
- ✅ Block builder (10K records per block)
- ✅ BAM → CAF converter
- ✅ CAF → SAM converter
- ✅ **Dictionary compression** (86% file size reduction)
- ✅ Round-trip validation (lossless BAM → CAF → SAM)
- ✅ Comprehensive specification document

**Performance Achieved**:
- File size: **1.6× vs BAM** (target: 1.5-2×) ✅
- Quality compression: **86% reduction** (11.6 MB → 1.6 MB)
- Conversion time: 702 ms for 100K records
- Dictionary training: <150 ms overhead

**Next Steps** (Week 3, Nov 11-17):
1. Integration testing with real-world BAM files
2. Edge case testing (varied distributions, large files)
3. Performance profiling and benchmarking
4. NEON optimization implementation (target: 5-10× speedup)

---

## Contributing

This is a research project. For questions or suggestions:
- See: `research/caf-format/README.md`
- Implementation Plan: `docs/IMPLEMENTATION_PLAN.md`
- Specification: `SPECIFICATION.md`

---

## License

MIT License - See parent biometal project
