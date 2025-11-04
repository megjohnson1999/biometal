# biometal Architecture

This document describes the internal architecture of biometal, focusing on the streaming-first design and network streaming implementation.

## Table of Contents

1. [Design Principles](#design-principles)
2. [Streaming Architecture](#streaming-architecture)
3. [Network Streaming](#network-streaming)
4. [Memory Management](#memory-management)
5. [ARM NEON Integration](#arm-neon-integration)
6. [I/O Optimization](#io-optimization)
7. [Error Handling](#error-handling)

---

## Design Principles

biometal is built on three core principles:

### 1. Evidence-Based Design

Every optimization in biometal is validated through experimental measurement:
- **1,357 experiments** with 40,710 measurements (N=30)
- **Statistical rigor**: 95% confidence intervals, Cohen's d effect sizes
- **Lab notebook**: 33 entries documenting full experimental log
- **See**: [OPTIMIZATION_RULES.md](../OPTIMIZATION_RULES.md) for complete evidence base

### 2. Streaming-First Architecture

All data structures designed for constant memory, not batch processing:
- **Iterator pattern**: Process one record at a time
- **No accumulation**: Never collect all records in memory
- **Block-based processing**: 10K records per block (Rule 2)
- **Memory target**: ~5 MB regardless of dataset size (Rule 5)

### 3. Platform Portability

Code works across ARM and x86_64 platforms:
- **ARM NEON**: 16-25× speedup on ARM platforms (Rule 1)
- **Scalar fallback**: Automatic on x86_64
- **No vendor lock-in**: Works on Mac, Graviton, Ampere, Raspberry Pi
- **Feature flags**: Optional network streaming, Python bindings

---

## Streaming Architecture

### Core Abstraction: DataSource

All input sources unified through a single enum:

```rust
pub enum DataSource {
    Local(PathBuf),              // Local file path
    Http(String),                // HTTP/HTTPS URL
    Sra(String),                 // SRA accession (auto-converted to HTTP)
}
```

**Design rationale**:
- User code identical regardless of source
- Network streaming transparent to application
- Easy to add new sources (S3, GCS, etc.)

### Streaming Parser: FastqStream

```rust
pub struct FastqStream<R: BufRead> {
    reader: R,
    line_buffer: String,
    // No Vec<FastqRecord> - constant memory!
}

impl<R: BufRead> Iterator for FastqStream<R> {
    type Item = io::Result<FastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        // Read one record, return, discard
        // Memory footprint constant
    }
}
```

**Key insight**: By implementing `Iterator`, we force streaming semantics:
- User cannot accidentally accumulate all records
- Compiler optimizations work better (iterator fusion)
- Composable with standard library (`map`, `filter`, `take`, etc.)

### Block-Based Processing

**Problem**: NEON speedup can be lost with naive streaming (Entry 027)

**Solution**: Process in blocks of 10K records (Evidence: Entry 027, 1,440 measurements)

```rust
const BLOCK_SIZE: usize = 10_000;

pub struct BlockProcessor<R: BufRead> {
    stream: FastqStream<R>,
    block_buffer: Vec<FastqRecord>,  // Reused buffer
}

impl<R: BufRead> BlockProcessor<R> {
    fn next_block(&mut self) -> Option<Vec<ProcessedResult>> {
        self.block_buffer.clear();

        // Fill block (up to 10K records)
        while self.block_buffer.len() < BLOCK_SIZE {
            match self.stream.next() {
                Some(Ok(record)) => self.block_buffer.push(record),
                Some(Err(e)) => return Some(Err(e)),
                None => break,
            }
        }

        if self.block_buffer.is_empty() {
            return None;
        }

        // Process entire block with NEON
        Some(process_block_neon(&self.block_buffer))
    }
}
```

**Tradeoff analysis**:
- **Pro**: Preserves NEON speedup (avoids 82-86% overhead)
- **Pro**: Better cache locality (sequential block access)
- **Con**: 10K records in memory (~1-2 MB)
- **Verdict**: Worth it - 16-25× NEON speedup dominates memory cost

---

## Network Streaming

Network streaming is the most complex part of biometal. It addresses the I/O bottleneck that dominates bioinformatics workloads (Entry 028: 264-352× slower than compute).

### Architecture Overview

```
┌─────────────────────────────────────────────────────────────┐
│                        FastqStream                          │
│                    (Application Layer)                      │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│                       HttpReader                            │
│              (Streaming + Prefetch Layer)                   │
│  • Range requests (64 KB chunks)                            │
│  • Background prefetching (4 blocks ahead)                  │
│  • Chunk size configuration                                 │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│                      HttpClient                             │
│                  (Caching + Retry Layer)                    │
│  • LRU cache (50 MB byte-bounded)                           │
│  • Retry logic (3 attempts, exponential backoff)            │
│  • Concurrent request handling                              │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│                        reqwest                              │
│                   (HTTP Client Library)                     │
└─────────────────────────────────────────────────────────────┘
```

### Component 1: HttpClient (Caching + Retry)

**Responsibilities**:
1. Byte-bounded LRU cache (50 MB default)
2. Retry logic with exponential backoff
3. Thread-safe operation (Arc + Mutex)

**Implementation**:

```rust
pub struct HttpClient {
    client: reqwest::blocking::Client,
    cache: Arc<Mutex<LruCache<CacheKey, Vec<u8>>>>,
}

impl HttpClient {
    pub fn fetch_range(&self, url: &str, start: u64, end: u64)
        -> io::Result<Vec<u8>>
    {
        let key = CacheKey { url: url.to_string(), start, end };

        // Check cache first
        {
            let mut cache = self.cache.lock().unwrap();
            if let Some(data) = cache.get(&key) {
                return Ok(data.clone());  // Cache hit!
            }
        }

        // Cache miss - fetch from network with retry
        let data = self.fetch_with_retry(url, start, end)?;

        // Store in cache
        {
            let mut cache = self.cache.lock().unwrap();
            cache.put(key, data.clone());
        }

        Ok(data)
    }

    fn fetch_with_retry(&self, url: &str, start: u64, end: u64)
        -> io::Result<Vec<u8>>
    {
        let mut delay = Duration::from_millis(1000);

        for attempt in 1..=3 {
            match self.do_fetch(url, start, end) {
                Ok(data) => return Ok(data),
                Err(e) if attempt < 3 => {
                    eprintln!("Attempt {} failed, retrying in {:?}",
                              attempt, delay);
                    thread::sleep(delay);
                    delay *= 2;  // Exponential backoff
                }
                Err(e) => return Err(e),  // Final attempt failed
            }
        }

        unreachable!()
    }
}
```

**Cache design choices**:

| Choice | Rationale |
|--------|-----------|
| LRU eviction | Sequential streaming has good locality |
| Byte-bounded (50 MB) | Prevents unbounded growth with large blocks |
| Thread-safe (Mutex) | Enables background prefetching |
| Cache entire ranges | Avoids partial block complications |

### Component 2: HttpReader (Streaming + Prefetch)

**Responsibilities**:
1. Implement `Read` trait for streaming
2. Background prefetching to hide latency
3. Sequential block reading with range requests

**Implementation**:

```rust
pub struct HttpReader {
    client: HttpClient,
    url: String,
    position: u64,
    total_size: Option<u64>,
    chunk_size: usize,        // Default: 64 KB
    prefetch_count: usize,    // Default: 4 blocks
}

impl Read for HttpReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        // Calculate range for current read
        let start = self.position;
        let end = start + self.chunk_size as u64;

        // Fetch current block (may be cached)
        let data = self.client.fetch_range(&self.url, start, end)?;

        // Copy to user buffer
        let n = buf.len().min(data.len());
        buf[..n].copy_from_slice(&data[..n]);

        self.position += n as u64;

        // Trigger background prefetching for upcoming blocks
        if self.prefetch_count > 0 {
            self.trigger_prefetch();
        }

        Ok(n)
    }
}

impl HttpReader {
    fn trigger_prefetch(&self) {
        let mut ranges = Vec::with_capacity(self.prefetch_count);
        let mut prefetch_position = self.position;

        for _ in 0..self.prefetch_count {
            let start = prefetch_position;
            let end = start + self.chunk_size as u64;

            // Check EOF
            if let Some(total) = self.total_size {
                if start >= total {
                    break;
                }
            }

            ranges.push((start, end));
            prefetch_position = end;
        }

        if !ranges.is_empty() {
            self.client.prefetch(&self.url, &ranges);
        }
    }
}
```

**Prefetch design**:

```
Current read position: 128 KB
Prefetch count: 4

┌──────────┬──────────┬──────────┬──────────┬──────────┬──────────┐
│  64 KB   │  64 KB   │  64 KB   │  64 KB   │  64 KB   │  64 KB   │
└──────────┴──────────┴──────────┴──────────┴──────────┴──────────┘
     ▲          │          │          │          │
     │          │          │          │          │
  Current    Prefetch  Prefetch  Prefetch  Prefetch
   Read        #1        #2        #3        #4
  (from     (background threads fetch these)
  cache)
```

**Latency hiding analysis**:

Without prefetching:
```
Read 1: 50ms network → Process 1ms
Read 2: 50ms network → Process 1ms
Read 3: 50ms network → Process 1ms
Total: 153ms for 3 reads
```

With prefetching (4 blocks ahead):
```
Read 1: 50ms network → Process 1ms → (Prefetch 2-5 starts)
Read 2: 0ms (cached!) → Process 1ms
Read 3: 0ms (cached!) → Process 1ms
Total: 53ms for 3 reads (2.9× faster)
```

**Evidence**: Entry 028 - I/O bottleneck 264-352× slower than compute. Prefetching essential.

### Component 3: SRA Integration

**Problem**: SRA accessions need to be converted to HTTP URLs

**Solution**: Automatic conversion using NCBI's S3 public access pattern

```rust
const SRA_BASE_URL: &str = "https://sra-pub-run-odp.s3.amazonaws.com/sra";

pub fn sra_to_url(accession: &str) -> Result<String> {
    // Validate: SRR/SRX/SRS/SRP + digits
    validate_accession(accession)?;

    // Extract directory prefix (first 6 chars)
    let dir_prefix = &accession[..6];

    // Build URL: base/prefix/accession/accession
    Ok(format!("{}/{}/{}/{}",
        SRA_BASE_URL, dir_prefix, accession, accession))
}
```

**URL pattern examples**:
```
SRR000001 → https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR000/SRR000001/SRR000001
SRR390728 → https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR390/SRR390728/SRR390728
SRX123456 → https://sra-pub-run-odp.s3.amazonaws.com/sra/SRX123/SRX123456/SRX123456
```

**Integration with DataSource**:

```rust
impl DataSource {
    pub fn open(&self) -> io::Result<Box<dyn BufRead>> {
        match self {
            DataSource::Local(path) => {
                // Local file with smart mmap + parallel bgzip
                open_local(path)
            }
            DataSource::Http(url) => {
                // HTTP streaming with caching + prefetch
                let reader = HttpReader::new(url)?;
                Ok(Box::new(BufReader::new(reader)))
            }
            DataSource::Sra(accession) => {
                // Convert SRA → HTTP, then stream
                let url = sra_to_url(accession)?;
                let reader = HttpReader::new(&url)?;
                Ok(Box::new(BufReader::new(reader)))
            }
        }
    }
}
```

**Design benefit**: SRA streaming is completely transparent to application code.

---

## Memory Management

### Memory Budget

```
Total Memory = Streaming + Cache + Prefetch + Operations

Component           Default     Range        Notes
─────────────────────────────────────────────────────────────
Streaming buffer    5 MB        5 MB         Constant (Rule 5)
LRU cache          50 MB       10-200 MB     Byte-bounded
Prefetch buffer   256 KB       128 KB-4 MB   (count × chunk_size)
Block buffer      1-2 MB       1-2 MB        10K records (Rule 2)
NEON operations   <100 KB      <100 KB       In-place operations
─────────────────────────────────────────────────────────────
TOTAL             ~55 MB       ~16-206 MB    Regardless of dataset
```

### Constant Memory Guarantee

**Key insight**: Streaming + bounded cache = constant memory

```rust
// ❌ BAD: Accumulates in memory (grows with dataset)
let records: Vec<_> = stream.collect();

// ✅ GOOD: Constant memory (streaming)
for record in stream {
    process(record?);
}
```

**Evidence**: Entry 026 - 99.5% memory reduction (1,344 MB → 5 MB for 1M sequences)

---

## ARM NEON Integration

### SIMD Operations

ARM NEON provides 16-25× speedup for element-wise operations (Rule 1).

**Architecture pattern**:

```rust
#[cfg(target_arch = "aarch64")]
pub fn operation_neon(input: &[u8]) -> Result {
    use std::arch::aarch64::*;

    unsafe {
        // NEON intrinsics (16 bytes at a time)
        // See src/operations/*.rs for implementations
    }
}

#[cfg(not(target_arch = "aarch64"))]
pub fn operation_scalar(input: &[u8]) -> Result {
    // Portable scalar fallback
}

pub fn operation(input: &[u8]) -> Result {
    #[cfg(target_arch = "aarch64")]
    { operation_neon(input) }

    #[cfg(not(target_arch = "aarch64"))]
    { operation_scalar(input) }
}
```

**Speedup validation** (Evidence: Entry 020-025, 9,210 measurements):

| Operation | Speedup | Cohen's d | CI (95%) |
|-----------|---------|-----------|----------|
| base_counting | 16.7× | 4.82 | ±0.14 |
| gc_content | 20.3× | 5.12 | ±0.18 |
| mean_quality | 25.1× | 5.87 | ±0.21 |

---

## I/O Optimization

### Parallel Bgzip Decompression (Rule 3)

**Evidence**: Entry 029 - 6.5× speedup

```rust
pub fn decompress_bgzip_parallel(compressed: &[u8]) -> io::Result<Vec<u8>> {
    use rayon::prelude::*;

    // Parse BGZF blocks
    let blocks = parse_bgzip_blocks(compressed)?;

    // Decompress in parallel
    let decompressed: Vec<_> = blocks
        .par_iter()
        .map(|block| decompress_block(block))
        .collect::<io::Result<Vec<_>>>()?;

    Ok(decompressed.concat())
}
```

### Smart mmap (Rule 4)

**Evidence**: Entry 032 - Additional 2.5× speedup for large files on macOS

```rust
const MMAP_THRESHOLD: u64 = 50 * 1024 * 1024;  // 50 MB

pub fn open_with_mmap(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let size = std::fs::metadata(path)?.len();

    if size >= MMAP_THRESHOLD {
        // Use mmap for large files (2.5× faster)
        let mmap = unsafe { Mmap::map(&File::open(path)?)? };

        #[cfg(target_os = "macos")]
        {
            // Hint sequential access
            madvise(mmap.as_ptr(), mmap.len(), MADV_SEQUENTIAL)?;
        }

        Ok(Box::new(BufReader::new(Cursor::new(mmap))))
    } else {
        // Standard I/O for small files (avoid overhead)
        Ok(Box::new(BufReader::new(File::open(path)?)))
    }
}
```

**Threshold rationale**: 50 MB balances mmap benefits vs overhead

---

## Error Handling

### Error Types

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

### Retry Logic

Network failures handled automatically with exponential backoff:

```
Attempt 1: Fail → Wait 1s
Attempt 2: Fail → Wait 2s
Attempt 3: Fail → Return error
```

**Rationale**: Transient network issues are common, but persistent failures should bubble up.

---

## Performance Summary

### Memory

- **Streaming**: ~5 MB regardless of dataset size
- **Total**: ~55 MB with caching + prefetching
- **Comparison**: 99.5%+ reduction vs batch processing

### Compute

- **ARM NEON**: 16-25× speedup on element-wise operations
- **Scalar fallback**: Automatic on x86_64
- **Block processing**: Preserves NEON efficiency

### I/O

- **Parallel bgzip**: 6.5× speedup
- **Smart mmap**: Additional 2.5× for large files (macOS)
- **Combined**: 16.3× I/O speedup

### Network

- **Prefetching**: Hides network latency (2-3× effective speedup)
- **Caching**: Reduces redundant requests (50 MB LRU)
- **Retry**: Handles transient failures automatically

---

## Future Directions

### Planned Enhancements

1. **Configurable cache size** (currently hardcoded 50 MB)
2. **Compression-aware prefetching** (compressed blocks are smaller)
3. **Multi-connection streaming** (parallel range requests)
4. **GPU acceleration** (for highly parallel operations)

### Research Questions

1. **Optimal prefetch count**: May vary by network characteristics
2. **Adaptive chunk sizing**: Tune based on observed latency
3. **Cache replacement policy**: LRU vs LFU vs adaptive
4. **Compression ratio prediction**: Estimate block sizes for better prefetching

---

## References

- **Evidence base**: [OPTIMIZATION_RULES.md](../OPTIMIZATION_RULES.md)
- **Full methodology**: [apple-silicon-bio-bench](https://github.com/shandley/apple-silicon-bio-bench)
- **Lab notebook**: 33 entries, 1,357 experiments, 40,710 measurements

---

**Last Updated**: November 4, 2025
**Version**: v0.2.2 (Network Streaming Complete)
