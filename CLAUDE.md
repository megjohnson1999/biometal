# biometal: Claude Development Guide

**Project**: biometal - ARM-native bioinformatics library
**Status**: v0.1.0 (Early Development, Week 1-2 starting Nov 4, 2025)
**Timeline**: 6 weeks to v1.0.0 (Nov 4 - Dec 15, 2025)

---

## Mission

Democratize bioinformatics by enabling 5TB dataset analysis on consumer hardware through:
- **Streaming architecture**: Constant ~5 MB memory (not load-all)
- **ARM-native performance**: 16-25× NEON speedup
- **Network streaming**: Analyze without downloading
- **Evidence-based optimization**: Every rule validated experimentally

**Target audiences**: LMIC researchers, small labs, students, field researchers, ML practitioners

---

## Core Principles

### 1. Evidence-Based Design

Every optimization in biometal comes from validated experimental results (ASBB project):
- **Follow OPTIMIZATION_RULES.md**: 6 rules distilled from 1,357 experiments
- **Don't guess**: If unsure about an optimization, refer to ASBB evidence
- **Document rationale**: Link implementations to specific rules/lab notebook entries

**Example**:
```rust
// Rule 2: Block size from Entry 027 (1,440 measurements)
const BLOCK_SIZE: usize = 10_000; // Evidence-based, not arbitrary
```

### 2. Streaming-First Architecture

**Always design for constant memory**, not batch processing:
- ❌ `Vec<FastqRecord>` - accumulates in memory
- ✅ `Iterator<Item = FastqRecord>` - constant memory

**Memory target**: ~5 MB regardless of dataset size (Rule 5)

### 3. ARM-Native with Portable Fallback

**Always provide both ARM and fallback implementations**:
```rust
#[cfg(target_arch = "aarch64")]
pub fn operation_neon(input: &[u8]) -> Result {
    // ARM NEON implementation (16-25× faster)
}

#[cfg(not(target_arch = "aarch64"))]
pub fn operation_scalar(input: &[u8]) -> Result {
    // Scalar fallback for x86_64
}

pub fn operation(input: &[u8]) -> Result {
    #[cfg(target_arch = "aarch64")]
    { operation_neon(input) }
    
    #[cfg(not(target_arch = "aarch64"))]
    { operation_scalar(input) }
}
```

**Platform support priority**: Mac → Linux ARM (Graviton) → x86_64 fallback

### 4. Production Quality from Day 1

This is a production library, not research code:
- ✅ Comprehensive error handling (`Result<T, BiometalError>`)
- ✅ Documentation with examples (every public API)
- ✅ Property-based testing (proptest)
- ✅ Benchmarks (criterion)
- ❌ Panics in library code (use `Result` instead)
- ❌ `unwrap()` or `expect()` in library code

---

## Project Structure

```
biometal/
├── src/
│   ├── lib.rs              # Public API, re-exports
│   ├── io/                 # Streaming parsers (Rules 3-5)
│   │   ├── mod.rs
│   │   ├── fastq.rs        # FASTQ streaming parser
│   │   ├── fasta.rs        # FASTA streaming parser
│   │   ├── compression.rs  # Bgzip + mmap (Rules 3-4)
│   │   └── network.rs      # HTTP/SRA streaming (Rule 6, Week 3-4)
│   ├── operations/         # NEON-optimized operations (Rule 1)
│   │   ├── mod.rs
│   │   ├── base_counting.rs
│   │   ├── gc_content.rs
│   │   ├── quality_filter.rs
│   │   └── ...
│   ├── optimization/       # Auto-detection and tuning
│   │   ├── mod.rs
│   │   ├── platform.rs     # Platform detection
│   │   └── thresholds.rs   # Evidence-based thresholds
│   ├── error.rs            # Error types
│   └── types.rs            # Common types (FastqRecord, etc.)
├── benches/
│   └── operations.rs       # Criterion benchmarks
├── examples/
│   ├── basic_streaming.rs
│   ├── network_analysis.rs
│   └── neon_operations.rs
├── docs/
│   └── architecture.md
├── OPTIMIZATION_RULES.md   # Evidence base (from ASBB)
├── README.md               # User-facing documentation
├── CLAUDE.md               # This file
└── Cargo.toml
```

---

## Development Phases

### Week 1-2: Core Infrastructure + I/O Optimization (Nov 4-15)

**Goals**:
- Streaming FASTQ/FASTA parser
- ARM NEON operations (base_counting, gc_content, quality_filter)
- Parallel bgzip + smart mmap (Rules 3-4)
- Block-based processing (Rule 2)

**Key files to create**:
1. `src/io/fastq.rs` - Streaming FASTQ parser (Rule 5)
2. `src/io/compression.rs` - Parallel bgzip + mmap (Rules 3-4)
3. `src/operations/base_counting.rs` - NEON implementation (Rule 1)
4. `src/operations/gc_content.rs` - NEON implementation (Rule 1)
5. `src/types.rs` - FastqRecord, FastaRecord, etc.

**Deliverable**: biometal v0.1.0 (local file streaming)

### Week 3-4: Network Streaming (Nov 18-29)

**Goals**:
- HTTP/HTTPS streaming with range requests
- Smart LRU caching (configurable size)
- Background prefetching
- SRA toolkit integration

**Key files to create**:
1. `src/io/network.rs` - HTTP streaming + caching (Rule 6)
2. `src/io/sra.rs` - SRA toolkit wrapper

**Deliverable**: biometal v0.2.0 (network streaming)

### Week 5-6: Python Bindings + Polish (Dec 2-13)

**Goals**:
- PyO3 wrappers for Python ecosystem
- K-mer utilities (for BERT preprocessing)
- Example notebooks
- Cross-platform testing

**Key files to create**:
1. `src/python/` - PyO3 bindings
2. `examples/*.ipynb` - Jupyter notebooks

**Deliverable**: biometal v0.3.0 (ML-ready)

### Week 7+: Production Release (Dec 16+)

**Goals**:
- Extended operation coverage
- Comprehensive documentation
- Cross-platform testing (Mac, Graviton, RPi)
- Publish to crates.io

**Deliverable**: biometal v1.0.0 (production)

---

## Implementation Guidelines

### Rule 1: ARM NEON SIMD (16-25× speedup)

**When**: Element-wise operations (complexity 0.30-0.40)

**Pattern**:
```rust
#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::*;

#[cfg(target_arch = "aarch64")]
pub unsafe fn count_bases_neon(seq: &[u8]) -> [u32; 4] {
    // Process 16 bytes at a time with NEON
    // See OPTIMIZATION_RULES.md Rule 1 for full example
}

#[cfg(not(target_arch = "aarch64"))]
pub fn count_bases_scalar(seq: &[u8]) -> [u32; 4] {
    // Scalar fallback
}
```

**Evidence**: Entry 020-025 (307 experiments, 9,210 measurements)

### Rule 2: Block-Based Processing (10K records)

**Why**: Preserves NEON speedup (avoids 82-86% overhead)

**Pattern**:
```rust
const BLOCK_SIZE: usize = 10_000; // From Entry 027

pub struct FastqStream<R: BufRead> {
    reader: R,
    block_buffer: Vec<FastqRecord>,
}

impl<R: BufRead> FastqStream<R> {
    fn process_block(&mut self) -> Result<ProcessedBlock> {
        self.block_buffer.clear();
        
        // Fill block (up to 10K records)
        while self.block_buffer.len() < BLOCK_SIZE {
            match self.read_record()? {
                Some(record) => self.block_buffer.push(record),
                None => break,
            }
        }
        
        // Process entire block with NEON
        let results = unsafe { process_block_neon(&self.block_buffer) };
        Ok(ProcessedBlock::new(results))
    }
}
```

**Evidence**: Entry 027 (1,440 measurements)

### Rule 3: Parallel Bgzip (6.5× speedup)

**When**: All bgzip-compressed files

**Pattern**:
```rust
use rayon::prelude::*;

pub fn decompress_bgzip_parallel(compressed: &[u8]) -> io::Result<Vec<u8>> {
    let blocks = parse_bgzip_blocks(compressed)?;
    
    let decompressed: Vec<_> = blocks
        .par_iter()
        .map(|block| decompress_block(block))
        .collect::<io::Result<Vec<_>>>()?;
    
    Ok(decompressed.concat())
}
```

**Evidence**: Entry 029 (CPU parallel prototype)

### Rule 4: Smart mmap (2.5× additional, threshold-based)

**When**: Files ≥50 MB on macOS

**Pattern**:
```rust
const MMAP_THRESHOLD: u64 = 50 * 1024 * 1024; // 50 MB

enum DataSource {
    StandardIo(Vec<u8>),
    MemoryMapped(Mmap),
}

impl DataSource {
    pub fn open(path: &Path) -> io::Result<Self> {
        let size = std::fs::metadata(path)?.len();
        
        if size >= MMAP_THRESHOLD {
            Self::open_mmap(path) // 2.5× faster for large files
        } else {
            Ok(Self::StandardIo(std::fs::read(path)?)) // Avoid overhead
        }
    }
}
```

**Evidence**: Entry 032 (scale validation, 0.54-544 MB)

### Rule 5: Constant-Memory Streaming (~5 MB)

**Always**: Design for streaming, not batch

**Pattern**:
```rust
pub struct FastqStream<R: BufRead> {
    reader: R,
    line_buffer: String,
    // NO Vec<FastqRecord> accumulation!
}

impl<R: BufRead> Iterator for FastqStream<R> {
    type Item = io::Result<FastqRecord>;
    
    fn next(&mut self) -> Option<Self::Item> {
        // Read one record, return, discard
        // Memory stays constant
    }
}
```

**Evidence**: Entry 026 (720 measurements, 99.5% reduction)

### Rule 6: Network Streaming (Week 3-4)

**Why**: I/O dominates 264-352× (makes network streaming critical)

**Pattern** (Week 3-4 implementation):
```rust
pub enum DataSource {
    Local(PathBuf),
    Http(Url),
    Sra(String),
}

pub struct StreamingReader {
    source: DataSource,
    cache: LruCache<BlockId, Vec<u8>>,
    prefetch: Prefetcher,
}
```

**Evidence**: Entry 028 (360 measurements)

---

## Error Handling

**Always use `Result` types**, never panic in library code:

```rust
#[derive(Debug, thiserror::Error)]
pub enum BiometalError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    
    #[error("Invalid FASTQ format at line {line}: {msg}")]
    InvalidFormat { line: usize, msg: String },
    
    #[error("Network error: {0}")]
    Network(String),
    
    #[error("Compression error: {0}")]
    Compression(String),
}

pub type Result<T> = std::result::Result<T, BiometalError>;
```

---

## Testing Strategy

### 1. Property-Based Testing (proptest)

```rust
use proptest::prelude::*;

proptest! {
    #[test]
    fn test_base_counting_matches_naive(seq in "[ACGT]{1,1000}") {
        let neon_result = count_bases_neon(seq.as_bytes());
        let naive_result = count_bases_naive(seq.as_bytes());
        prop_assert_eq!(neon_result, naive_result);
    }
}
```

### 2. Benchmarking (criterion)

```rust
use criterion::{criterion_group, criterion_main, Criterion};

fn bench_base_counting(c: &mut Criterion) {
    let seq = generate_sequence(100_000);
    
    c.bench_function("base_counting_neon", |b| {
        b.iter(|| count_bases_neon(&seq))
    });
}

criterion_group!(benches, bench_base_counting);
criterion_main!(benches);
```

---

## Documentation

**Every public API must have**:
1. Doc comment explaining what it does
2. Example showing how to use it
3. Link to evidence (when implementing optimizations)

**Example**:
```rust
/// Stream FASTQ records from a file with constant memory.
///
/// Uses block-based processing (10K records) to preserve ARM NEON speedup
/// while maintaining streaming benefits. Memory footprint remains constant
/// at ~5 MB regardless of file size.
///
/// # Evidence
///
/// - Rule 2 (Block-based): Entry 027, 1,440 measurements
/// - Rule 5 (Streaming): Entry 026, 99.5% memory reduction
///
/// # Example
///
/// ```
/// use biometal::FastqStream;
///
/// let stream = FastqStream::from_path("large.fq.gz")?;
/// for record in stream {
///     let record = record?;
///     // Process one record at a time
/// }
/// # Ok::<(), biometal::Error>(())
/// ```
pub struct FastqStream<R: BufRead> { /* ... */ }
```

---

## Common Pitfalls to Avoid

### 1. Accumulating Records in Memory

**Bad**:
```rust
let mut records = Vec::new();
for line in reader.lines() {
    records.push(parse_record(line)?);
}
// Memory grows linearly with file size!
```

**Good**:
```rust
for record in FastqStream::from_path(path)? {
    let record = record?;
    // Process immediately, no accumulation
}
// Memory stays constant at ~5 MB
```

### 2. Using `unwrap()` in Library Code

**Bad**:
```rust
pub fn operation(input: &[u8]) -> Output {
    parse(input).unwrap() // Panics on invalid input!
}
```

**Good**:
```rust
pub fn operation(input: &[u8]) -> Result<Output> {
    parse(input).map_err(|e| BiometalError::InvalidFormat {
        line: 0,
        msg: e.to_string(),
    })
}
```

### 3. Implementing Optimizations Without Evidence

**Bad**:
```rust
const BLOCK_SIZE: usize = 8_192; // Arbitrary choice
```

**Good**:
```rust
// Rule 2: Block size from Entry 027 (1,440 measurements)
// Evidence: 10K records balances SIMD efficiency and memory
const BLOCK_SIZE: usize = 10_000;
```

### 4. Platform-Specific Code Without Fallback

**Bad**:
```rust
pub fn operation(input: &[u8]) -> Result {
    // Only works on ARM!
    unsafe { operation_neon(input) }
}
```

**Good**:
```rust
pub fn operation(input: &[u8]) -> Result {
    #[cfg(target_arch = "aarch64")]
    { unsafe { operation_neon(input) } }
    
    #[cfg(not(target_arch = "aarch64"))]
    { operation_scalar(input) }
}
```

---

## For Claude: Session Guidelines

### What to Emphasize

- Evidence-based design (follow OPTIMIZATION_RULES.md)
- Streaming-first architecture (constant memory)
- ARM-native with portable fallback
- Production quality (error handling, docs, tests)

### What NOT to Do

- Don't make up optimization parameters (refer to evidence)
- Don't accumulate records in memory (streaming only)
- Don't panic in library code (use Result)
- Don't implement ARM-only code without scalar fallback

### Decision Framework

**When user asks "how should I implement X?"**:
1. Check OPTIMIZATION_RULES.md for relevant rule
2. Follow the implementation pattern for that rule
3. Link to evidence (lab notebook entry)

**When user proposes optimization**:
1. Is this validated in ASBB experiments?
2. If yes: Which rule/entry documents it?
3. If no: Suggest validating first or using proven approach

---

## Quick Reference

**Evidence base**: 1,357 experiments, 40,710 measurements (N=30)  
**Full methodology**: [apple-silicon-bio-bench](https://github.com/shandley/apple-silicon-bio-bench)  
**Optimization rules**: See [OPTIMIZATION_RULES.md](OPTIMIZATION_RULES.md)  
**Timeline**: 6 weeks to v1.0.0 (Nov 4 - Dec 15, 2025)

---

**Last Updated**: November 4, 2025  
**Phase**: Week 1-2 (Core Infrastructure + I/O Optimization)  
**Next Milestone**: biometal v0.1.0 (Nov 15, 2025)
