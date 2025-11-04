# biometal: Claude Development Guide

**Project**: biometal - ARM-native bioinformatics library with Apple Silicon breakthroughs
**Status**: v0.1.0 (Early Development, Week 1-2 starting Nov 4, 2025)
**Timeline**: 10 weeks to v1.1.0 (Nov 4 - Jan 12, 2026)

---

## Mission: Democratize Bioinformatics on Modern Laptops

Enable world-class genomics analysis on consumer hardware, making bioinformatics accessible to:
- **LMIC researchers** (limited compute budgets)
- **Small labs and students** (laptop-only workflows)
- **Field researchers** (portable analysis)
- **ML practitioners** (local preprocessing pipelines)

### Core Innovations

1. **Streaming Architecture** (Rule 5)
   - Constant ~5 MB memory regardless of dataset size
   - Analyze 5TB datasets without downloading
   - 99.5% memory reduction vs. traditional tools

2. **ARM NEON Performance** (Rule 1)
   - 16-25× speedup for sequence operations
   - Works across Mac (M-series), AWS Graviton, Raspberry Pi
   - Automatic scalar fallback for x86_64

3. **Apple Silicon Breakthroughs** ⭐ NEW
   - **Metal GPU + Unified Memory**: 10-50× pileup generation (zero-copy architecture)
   - **NEON CIGAR parsing**: 10-20× SAM/BAM operations (unpublished optimization)
   - **First bioinformatics library** to exploit Apple's true unified memory

4. **Network Streaming** (Rule 6)
   - Analyze remote datasets without downloading (264-352× I/O speedup)
   - Smart LRU caching + background prefetching
   - HTTP/HTTPS and SRA toolkit integration

5. **Evidence-Based Design**
   - Every optimization validated with N=30 statistical rigor
   - 1,357 experiments, 40,710 measurements (ASBB project)
   - Published methodology: [apple-silicon-bio-bench](https://github.com/shandley/apple-silicon-bio-bench)

**Vision**: Make M-series MacBooks as capable for genomics as $50K HPC clusters

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

**Platform support priority**: Mac (with Metal) → Linux ARM (Graviton) → x86_64 fallback

### 4. Apple Silicon: Exploit Unique Hardware

**Apple's unified memory architecture enables breakthrough optimizations impossible on CUDA/x86**:

**Metal GPU + Zero-Copy UMA**:
```rust
// Traditional GPU (CUDA): Memory copy bottleneck
CPU RAM ←[PCIe 16 GB/s]→ GPU RAM  // Slow!

// Apple Silicon: True unified memory
CPU ←[400 GB/s]→ Shared RAM ←[400 GB/s]→ GPU  // 25× bandwidth!
```

**When to use Metal**:
- ✅ Compute-bound operations (pileup, coverage, depth)
- ✅ Parallel accumulation across genomic positions
- ✅ Operations benefiting from 1000s of parallel threads
- ❌ I/O-bound operations (Metal won't help)

**Metal + NEON synergy**:
```rust
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
pub fn pileup_metal(bam: &[u8]) -> Result<Pileup> {
    // Parse with NEON, process with Metal GPU
}
```

**Hardware tiers**:
1. **Mac with Metal**: Full acceleration (NEON + Metal + UMA)
2. **Linux ARM**: NEON only (still 16-25× faster)
3. **x86_64**: Scalar fallback (portable)

### 5. Delegate Parsing, Own Acceleration

**Don't reinvent battle-tested parsers** - focus on unique optimizations:

**Strategy**:
- ✅ Use **noodles** for complex formats (BAM, CRAM, VCF)
- ✅ Implement simple formats ourselves (FASTQ, FASTA, SAM text)
- ✅ Add Apple Silicon acceleration layer on top
- ✅ Maintain consistent biometal API

**Example**:
```rust
// Delegate parsing to noodles
pub struct BamStream {
    reader: noodles::bam::Reader,  // Spec-compliant
}

// Add our unique Metal acceleration
pub fn generate_pileup_metal(bam: &BamStream) -> Result<Pileup> {
    // World-first: Metal GPU + UMA for pileup
}
```

### 6. Production Quality from Day 1

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
│   │   ├── fastq.rs        # FASTQ streaming parser (Week 1-2)
│   │   ├── fasta.rs        # FASTA streaming parser (Week 1-2)
│   │   ├── paired_end.rs   # Paired-end support (R1/R2, interleaved) (Week 2)
│   │   ├── compression.rs  # Parallel bgzip + mmap (Rules 3-4) (Week 1-2)
│   │   ├── network.rs      # HTTP/SRA streaming (Rule 6) (Week 3-4)
│   │   ├── sam.rs          # Text SAM parser (Week 7)
│   │   ├── bam.rs          # BAM wrapper (noodles + Metal) (Week 7-8)
│   │   └── cram.rs         # CRAM wrapper (noodles) (Week 8)
│   │
│   ├── operations/         # NEON-optimized operations (Rule 1)
│   │   ├── mod.rs
│   │   ├── base_counting.rs    # Week 1-2
│   │   ├── gc_content.rs       # Week 1-2
│   │   ├── quality_filter.rs   # Week 1-2
│   │   ├── complexity.rs       # Week 2
│   │   ├── cigar.rs            # NEON CIGAR parsing (Week 9) ⭐
│   │   └── kmer.rs             # K-mer generation (Week 5-6)
│   │
│   ├── bam_ops/            # BAM-specific operations ⭐ NEW
│   │   ├── mod.rs
│   │   ├── pileup.rs       # Metal GPU pileup (Week 9-10) ⭐
│   │   ├── coverage.rs     # Metal GPU coverage (Week 10) ⭐
│   │   └── depth.rs        # Depth calculation
│   │
│   ├── metal/              # Metal GPU support ⭐ NEW (macOS only)
│   │   ├── mod.rs
│   │   ├── shaders.metal   # Metal compute shaders
│   │   └── buffer.rs       # Zero-copy buffer management
│   │
│   ├── optimization/       # Auto-detection and tuning
│   │   ├── mod.rs
│   │   ├── platform.rs     # Platform detection (NEON, Metal, AMX)
│   │   └── thresholds.rs   # Evidence-based thresholds
│   │
│   ├── error.rs            # Error types ✅ DONE
│   └── types.rs            # Common types (FastqRecord, BamRecord, etc.)
│
├── benches/
│   ├── fastq_parsing.rs    # FASTQ benchmarks
│   ├── neon_operations.rs  # NEON benchmarks
│   ├── metal_pileup.rs     # Metal GPU benchmarks ⭐
│   └── compression.rs      # Bgzip benchmarks
│
├── examples/
│   ├── basic_fastq_streaming.rs
│   ├── paired_end_reads.rs
│   ├── network_streaming.rs
│   ├── bam_pileup_metal.rs  # Metal GPU example ⭐
│   └── format_conversion.rs
│
├── tests/
│   ├── integration/
│   │   ├── fastq_tests.rs
│   │   ├── bam_tests.rs
│   │   └── paired_end_tests.rs
│   └── fixtures/           # Test data
│
├── docs/
│   ├── architecture.md
│   ├── apple_silicon.md     # Apple Silicon optimizations ⭐
│   └── metal_programming.md # Metal GPU guide ⭐
│
├── OPTIMIZATION_RULES.md    # Evidence base (ASBB + new BAM experiments)
├── README.md                # User-facing documentation
├── CLAUDE.md                # This file
└── Cargo.toml
```

**Key**: ⭐ = Apple Silicon breakthrough feature, ✅ = Complete

---

## Development Phases: Phased Releases

### 🎯 v1.0: FASTQ/FASTA Core (6 weeks, Nov 4 - Dec 15)

**Mission**: Production-ready streaming FASTQ/FASTA with evidence-based optimization

#### Week 1-2: Core Infrastructure + I/O (Nov 4-15)

**Goals**:
- Streaming FASTQ/FASTA parser with constant memory
- Paired-end support (R1/R2 separate + interleaved)
- ARM NEON operations (base_counting, gc_content, quality_filter)
- Parallel bgzip + smart mmap (Rules 3-4)
- Block-based processing (Rule 2)

**Critical path**:
1. `src/io/fastq.rs` - FASTQ streaming parser ⚡ PRIORITY 1
2. `src/io/compression.rs` - Parallel bgzip + mmap ⚡ PRIORITY 1
3. `src/io/paired_end.rs` - Paired-end support
4. `src/operations/base_counting.rs` - First NEON operation (establishes pattern)
5. `src/operations/gc_content.rs` - Second NEON operation
6. Property tests + benchmarks

**Deliverable**: biometal v0.1.0 (local FASTQ/FASTA streaming)

#### Week 3-4: Network Streaming (Nov 18-29)

**Goals**:
- HTTP/HTTPS streaming with range requests
- Smart LRU caching (configurable size)
- Background prefetching
- Paired-end network support
- SRA toolkit integration

**Key files**:
1. `src/io/network.rs` - HTTP streaming + caching (Rule 6)
2. `src/io/sra.rs` - SRA toolkit wrapper
3. Integration tests with real network sources

**Deliverable**: biometal v0.2.0 (network streaming)

#### Week 5-6: Python Bindings + Polish (Dec 2-13)

**Goals**:
- PyO3 wrappers for Python ecosystem
- K-mer utilities (for ML BERT preprocessing)
- Example notebooks (Jupyter)
- Cross-platform testing (Mac, Graviton, x86_64)
- Comprehensive documentation

**Key files**:
1. `src/python/` - PyO3 bindings
2. `examples/*.ipynb` - Jupyter notebooks
3. Integration examples

**Deliverable**: biometal v1.0.0 (production FASTQ/FASTA library)

---

### 🚀 v1.1: BAM/SAM + Apple Silicon Breakthrough (4 weeks, Dec 16 - Jan 12)

**Mission**: World-first Metal GPU + UMA optimization for BAM operations

#### Week 7-8: BAM/SAM Foundation (Dec 16-29)

**Goals**:
- Delegate BAM/CRAM parsing to noodles (spec-compliant)
- Text-based SAM parser (simpler, we implement)
- Consistent biometal API wrapper
- Basic BAM streaming (no Metal yet)

**Key files**:
1. `src/io/sam.rs` - Text SAM parser with NEON optimization
2. `src/io/bam.rs` - noodles wrapper with biometal API
3. `src/io/cram.rs` - noodles wrapper
4. Integration tests

**Deliverable**: biometal v1.1.0-alpha (BAM support, no Metal yet)

#### Week 9: NEON CIGAR Optimization (Dec 30 - Jan 5)

**Goals**:
- NEON-optimized CIGAR string parsing (unpublished)
- Benchmark vs. scalar (expect 10-20× like Rule 1)
- Extend ASBB with Entry 034 (N=30 experiments)

**Key files**:
1. `src/operations/cigar.rs` - NEON CIGAR parser ⭐
2. Benchmarks + validation experiments

**Expected**: 10-20× speedup (analogous to Rule 1)

**Deliverable**: biometal v1.1.0-beta (NEON BAM operations)

#### Week 10: Metal GPU Pileup (Jan 6-12)

**Goals**:
- Metal compute shader for pileup generation
- Zero-copy UMA buffer management
- Benchmark vs. SAMtools/PaCBAM
- Extend ASBB with Entry 035 (N=30 experiments)

**Key files**:
1. `src/metal/shaders.metal` - Metal compute shaders ⭐
2. `src/bam_ops/pileup.rs` - Metal pileup generation ⭐
3. `src/metal/buffer.rs` - Zero-copy buffers

**Expected**: 10-50× speedup (vs. SAMtools mpileup)

**Deliverable**: biometal v1.1.0 (BAM with Metal breakthrough) 🎉

---

### 🌟 v1.2+: Extended BAM Operations (2+ weeks, Jan 13+)

**Mission**: Comprehensive BAM toolkit with full Metal acceleration

#### Week 11+: Coverage, Depth, and Beyond

**Goals**:
- Metal GPU coverage calculation (Entry 036)
- Additional NEON optimizations (quality score operations)
- BED, VCF, GFF/GTF parsers (optional)
- Neural Engine exploration (future: ML variant calling)

**Key files**:
1. `src/bam_ops/coverage.rs` - Metal coverage ⭐
2. `src/bam_ops/depth.rs` - Depth calculation
3. Annotation format parsers (BED, VCF, GTF)

**Deliverable**: biometal v1.2.0+ (comprehensive genomics toolkit)

---

## Apple Silicon Breakthroughs: Metal + UMA

### Why Apple Silicon is Different

**Traditional GPU Computing (CUDA)**:
```
Problem: PCIe bottleneck limits genomics speedup

CPU RAM ─┬─[PCIe: 16 GB/s]─→ GPU RAM
         │                   ↓
         │              GPU Compute
         │                   ↓
         └─[PCIe: 16 GB/s]←─ Results

Bottleneck: Memory copies eat 50-80% of potential speedup
```

**Apple Silicon Unified Memory**:
```
Breakthrough: Zero-copy shared memory

         ┌─→ CPU Cores
         │
Shared RAM ──┬── [400 GB/s internal fabric]
         │
         ├─→ GPU Cores (Metal)
         │
         ├─→ Neural Engine
         │
         └─→ AMX Matrix Coprocessor

Advantage: No memory copies, 25× bandwidth vs. PCIe
```

### Rule 7: Metal GPU for Compute-Bound Operations (NEW)

**When to use Metal** (Week 9-10):
- ✅ Pileup generation (accumulate 50-100 reads per position)
- ✅ Coverage calculation (parallel position counting)
- ✅ Depth computation (statistical accumulation)
- ❌ I/O-bound operations (parsing, filtering - use NEON instead)

**Pattern**:
```rust
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
pub fn generate_pileup_metal(
    bam_records: &[BamRecord],
    region: GenomicRegion,
) -> Result<Pileup> {
    use metal::{Device, CommandQueue, MTLResourceOptions};

    let device = Device::system_default()
        .ok_or(BiometalError::MetalNotAvailable)?;

    // Zero-copy: Share CPU buffer with GPU
    let shared_buffer = device.new_buffer_with_data(
        bam_records.as_ptr() as *const _,
        (bam_records.len() * std::mem::size_of::<BamRecord>()) as u64,
        MTLResourceOptions::StorageModeShared, // KEY: Shared mode = zero-copy
    );

    // Metal compute shader processes in parallel
    let pileup = compute_pileup_gpu(device, shared_buffer, region)?;

    Ok(pileup)
}

// Fallback for non-macOS ARM (Graviton, RPi)
#[cfg(all(target_arch = "aarch64", not(target_os = "macos")))]
pub fn generate_pileup_metal(
    bam_records: &[BamRecord],
    region: GenomicRegion,
) -> Result<Pileup> {
    // Use NEON CPU parallel instead
    generate_pileup_neon(bam_records, region)
}
```

**Evidence**: Entry 035 (Week 10, to be validated with N=30)

**Expected speedup**: 10-50× vs. SAMtools mpileup (based on CUDA analogues + UMA advantage)

### Metal Programming Guidelines

**Shader structure** (`src/metal/shaders.metal`):
```metal
#include <metal_stdlib>
using namespace metal;

// Accumulate bases at each genomic position
kernel void pileup_accumulate(
    constant BamRecord* records [[buffer(0)]],
    device atomic_uint* pileup [[buffer(1)]],
    constant uint& num_records [[buffer(2)]],
    constant uint& start_position [[buffer(3)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= num_records) return;

    const BamRecord& record = records[gid];
    uint pos = record.position;

    // Each thread processes one read
    for (uint i = 0; i < record.length; i++) {
        uint genome_pos = pos + i - start_position;
        atomic_fetch_add_explicit(&pileup[genome_pos], 1, memory_order_relaxed);
    }
}
```

**When to use Metal vs. NEON**:

| Operation | Use Metal? | Use NEON? | Rationale |
|-----------|-----------|-----------|-----------|
| Pileup generation | ✅ Yes | Also yes (fallback) | Embarrassingly parallel, 1000s of positions |
| Coverage calculation | ✅ Yes | Also yes (fallback) | Parallel accumulation |
| CIGAR parsing | ❌ No | ✅ Yes | String parsing (SIMD efficient) |
| Base counting | ❌ No | ✅ Yes | NEON already 16-25× (Rule 1) |
| Quality filtering | ❌ No | ✅ Yes | NEON already 25× (Rule 1) |
| File I/O | ❌ No | Use rayon | I/O bound, not compute |

**Key insight**: Metal shines for operations with massive parallelism (1000s of independent positions). NEON better for sequence-level operations (already proven 16-25×).

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

### Evidence Base
- **Current**: 1,357 experiments, 40,710 measurements (N=30)
- **Source**: [apple-silicon-bio-bench](https://github.com/shandley/apple-silicon-bio-bench)
- **Rules**: 6 FASTQ/FASTA optimization rules (OPTIMIZATION_RULES.md)
- **Future**: +90 BAM experiments (Week 9-10) → 1,447 total experiments

### Timeline & Milestones

**v1.0 (6 weeks): FASTQ/FASTA Core** (Nov 4 - Dec 15, 2025)
- ✅ Streaming architecture (Rule 5)
- ✅ ARM NEON operations (Rule 1, 16-25×)
- ✅ Parallel bgzip + mmap (Rules 3-4, 16.3×)
- ✅ Paired-end support (R1/R2 + interleaved)
- ✅ Network streaming (Rule 6)
- ✅ Python bindings (PyO3)

**v1.1 (4 weeks): BAM + Apple Silicon Breakthrough** (Dec 16 - Jan 12, 2026)
- ⭐ Metal GPU pileup (10-50× expected, Entry 035)
- ⭐ NEON CIGAR parsing (10-20× expected, Entry 034)
- ✅ noodles BAM/CRAM wrappers
- ✅ Text SAM parser
- **World-first**: Metal GPU + UMA for genomics

**v1.2+ (2+ weeks): Extended Operations** (Jan 13+)
- Metal coverage/depth (Entry 036)
- Annotation formats (BED, VCF, GTF)
- Neural Engine exploration

### Platform Tiers
1. **Mac (Apple Silicon)**: Full acceleration (NEON + Metal + UMA)
2. **Linux ARM** (Graviton, RPi): NEON only (16-25×)
3. **x86_64**: Scalar fallback (portable)

### Unique Advantages
- Zero-copy unified memory (25× bandwidth vs. CUDA PCIe)
- First bioinformatics library with Metal GPU
- Evidence-based design (every optimization validated)
- Democratize genomics on consumer laptops

---

**Last Updated**: November 4, 2025
**Current Phase**: Week 1-2 (Core Infrastructure + I/O Optimization)
**Current Task**: Implementing FASTQ parser + compression module
**Next Milestone**: biometal v0.1.0 (Nov 15, 2025)
