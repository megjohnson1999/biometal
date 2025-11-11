# biometal Performance Optimization Guide

**Version**: 1.6.0
**Last Updated**: November 10, 2025

This guide helps you maximize biometal's performance for your specific use case.

---

## Table of Contents

1. [Quick Wins](#quick-wins)
2. [Understanding the 6 Optimization Rules](#understanding-the-6-optimization-rules)
3. [Platform-Specific Optimization](#platform-specific-optimization)
4. [Workflow-Specific Tips](#workflow-specific-tips)
5. [Performance Profiling](#performance-profiling)
6. [Common Anti-Patterns](#common-anti-patterns)
7. [Benchmarking Your Code](#benchmarking-your-code)

---

## Quick Wins

These are the fastest ways to improve performance:

### 1. Use ARM Hardware
**Impact**: 16-25× speedup
**Cost**: Hardware upgrade

- **Apple Silicon** (M1/M2/M3/M4): 16-25× speedup
- **AWS Graviton**: 6-10× speedup
- **Intel/AMD**: 1× (scalar fallback)

### 2. Stream, Don't Accumulate
**Impact**: 99.5% memory reduction, enables processing of arbitrarily large files
**Cost**: Free (code pattern change)

```python
# ❌ BAD: Accumulates in memory
records = []
for record in stream:
    records.append(record)
process_all(records)  # Requires 100× more memory

# ✅ GOOD: Constant memory
for record in stream:
    process(record)  # Immediate processing, ~5 MB total
```

### 3. Use Indexed Queries for BAM
**Impact**: 1.68-500× speedup (increases with file size)
**Cost**: Free (requires `.bai` index file)

```python
# ❌ SLOW: Full scan (O(n))
for record in biometal.BamReader.from_path("file.bam"):
    if record.position and 1000000 <= record.position < 2000000:
        process(record)

# ✅ FAST: Indexed query (O(log n))
index = biometal.BaiIndex.from_path("file.bam.bai")
for record in biometal.BamReader.query_region(
    "file.bam", index, "chr1", 1000000, 2000000
):
    process(record)
```

### 4. Let biometal Handle Decompression
**Impact**: Automatic parallel decompression for BGZF files
**Cost**: Free

```python
# ✅ GOOD: biometal handles .gz automatically
stream = biometal.FastqStream.from_path("reads.fastq.gz")

# ❌ BAD: Manual decompression adds overhead
import gzip
with gzip.open("reads.fastq.gz") as f:
    # Slower, more memory
```

### 5. Use Batch Processing for Network Streaming
**Impact**: Better throughput, reduced latency
**Cost**: Minor code complexity

```python
# ✅ Better throughput
batch = []
for record in biometal.FastqStream.from_url(url):
    batch.append(record)
    if len(batch) >= 1000:
        process_batch(batch)
        batch = []
```

---

## Understanding the 6 Optimization Rules

biometal's design is based on 6 evidence-based optimization rules validated through 1,357 experiments (40,710 measurements). Understanding these rules helps you write optimal code.

### Rule 1: ARM NEON SIMD (16-25× speedup)

**What**: Single Instruction, Multiple Data parallelism on ARM processors
**When**: Element-wise operations (base counting, GC content, quality scores)
**Evidence**: Entry 020-025 (307 experiments, 9,210 measurements)

**How it Works**:
```
Scalar (x86_64):          NEON (ARM):
Process 1 byte at a time  Process 16 bytes at once
A → check → count         AAAAAAAAAAAAAAAA → check all → count all
C → check → count         CCCCCCCCCCCCCCCC → check all → count all
G → check → count         GGGGGGGGGGGGGGGG → check all → count all
(slow)                    (16× faster)
```

**When Applied**:
- `count_bases()`: 16.7× speedup
- `gc_content()`: 20.3× speedup
- `quality_filter()`: 25.1× speedup
- BAM sequence decoding: 4.62× speedup

**User Action**: None required! biometal detects ARM hardware and uses NEON automatically.

### Rule 2: Block-Based Processing (10K records)

**What**: Process records in blocks instead of one-at-a-time
**When**: Streaming contexts where you want to preserve NEON speedup
**Evidence**: Entry 027 (1,440 measurements)

**Why It Matters**:
- Record-by-record NEON overhead: 82-86%
- Block size 10K overhead: 4-8%
- Sweet spot: ~1.5 MB for 150bp reads

**Implementation**: biometal's `FastqStream` uses block-based processing internally. You benefit automatically when iterating.

**User Action**: None required for FASTQ/FASTA streaming. Already implemented.

### Rule 3: Parallel Bgzip (6.5× speedup)

**What**: Decompress multiple BGZF blocks in parallel
**When**: All bgzip-compressed files (.bam, .fastq.gz with BGZF format)
**Evidence**: Entry 029 (CPU parallel prototype)
**Status**: ⏳ Planned for Phase 2 (not yet implemented)

**Expected Impact**:
- BGZF decompression: 6.5× faster
- BAM parsing: 5× → 32.5× (combined with existing optimizations)

**User Action**: None required once implemented (automatic).

### Rule 4: Smart mmap (2.5× additional speedup)

**What**: Memory-map files ≥50 MB on macOS for faster I/O
**When**: Large files on Mac ARM (≥50 MB threshold)
**Evidence**: Entry 032 (scale validation, 0.54-544 MB)
**Status**: ⏳ Planned for Phase 2 (not yet implemented)

**Expected Impact**:
- File I/O: 2.5× faster for large files
- Threshold-based: Only activates when beneficial

**User Action**: None required once implemented (automatic).

### Rule 5: Constant-Memory Streaming (~5 MB)

**What**: Design for streaming, not batch processing
**When**: Always! Core architectural principle
**Evidence**: Entry 026 (720 measurements, 99.5% reduction)

**Memory Comparison**:

| Dataset Size | Traditional | biometal | Reduction |
|--------------|-------------|----------|-----------|
| 100K seqs | 134 MB | 5 MB | 96.3% |
| 1M seqs | 1,344 MB | 5 MB | 99.5% |
| 5TB | 5,000 GB | 5 MB | 99.9999% |

**User Action**: **Critical!** You must use streaming patterns:

```python
# ❌ BAD: Accumulates in memory
records = list(biometal.FastqStream.from_path("file.fq.gz"))
# Memory: 1,344 MB for 1M sequences

# ✅ GOOD: Constant memory
for record in biometal.FastqStream.from_path("file.fq.gz"):
    process(record)
# Memory: ~5 MB regardless of file size
```

### Rule 6: Network Streaming

**What**: Stream data directly from HTTP/HTTPS without downloading
**Why**: I/O dominates 264-352× more than computation
**Evidence**: Entry 028 (360 measurements)

**Benefits**:
- Start analysis immediately (no download wait)
- No disk space required
- Sample large files without full download

**User Action**: Use `.from_url()` methods:

```python
# Stream from HTTP (no download)
stream = biometal.FastqStream.from_url("https://example.com/reads.fq.gz")
for record in stream:
    process(record)
```

**Performance Tips**:
- Network speed is the bottleneck (not CPU)
- Use batch processing for better throughput
- Consider local caching for repeated access

---

## Platform-Specific Optimization

### Apple Silicon (M1/M2/M3/M4)

**Baseline Performance**: 16-25× speedup (ARM NEON)

**Additional Optimizations**:

1. **Use latest hardware**: M3/M4 > M2 > M1
2. **High-Performance mode**: `sudo pmset -a forcechargeidle 0`
3. **Thermal management**: Ensure good cooling for sustained workloads
4. **Memory bandwidth**: Unified memory architecture is already optimal

**Expected Performance**:
- Base counting: 5,254 Kseq/s (16.7× vs scalar)
- GC content: 5,954 Kseq/s (20.3× vs scalar)
- Quality filter: 6,143 Kseq/s (25.1× vs scalar)
- BAM parsing: 55.1 MiB/s (5.0× vs scalar, will be 32.5× with parallel BGZF)

### AWS Graviton (ARM)

**Baseline Performance**: 6-10× speedup (ARM NEON, portable)

**Additional Optimizations**:

1. **Instance type**: Use Graviton 3 (c7g) > Graviton 2 (c6g)
2. **Network**: Use enhanced networking for streaming workloads
3. **Storage**: Use NVMe instance storage for local files
4. **Parallelism**: Run multiple instances for embarrassingly parallel workflows

**Expected Performance**:
- Roughly 50-60% of Apple Silicon performance
- Still significantly faster than x86_64 (6-10× vs 1×)

### Linux x86_64 / Intel Mac

**Baseline Performance**: 1× (scalar fallback, no SIMD)

**Optimization Strategy**: Minimize computation, maximize I/O efficiency

1. **Use indexed queries**: BAI indexing still provides 1.68-500× speedup
2. **Network streaming**: I/O bottleneck means network streaming is still viable
3. **Constant memory**: Streaming architecture still provides 99.5% memory reduction
4. **Consider ARM**: For compute-intensive workloads, consider migrating to ARM hardware

**When x86_64 is Fine**:
- I/O-bound workflows (network streaming, disk reads)
- Small files (<100 MB)
- Infrequent processing
- Already have x86_64 infrastructure

**When to Consider ARM**:
- CPU-bound workflows (base counting, GC content, quality filtering)
- Large files (>1 GB)
- Frequent processing (daily pipelines)
- New infrastructure investments

---

## Workflow-Specific Tips

### FASTQ Quality Control

**Bottleneck**: Base counting, GC content, quality score calculations

**Optimizations**:
1. ✅ Use ARM hardware (16-25× speedup)
2. ✅ Stream, don't accumulate
3. ✅ Use biometal's optimized operations

**Example**:
```python
import biometal

# Optimized QC workflow
total_reads = 0
high_quality = 0
total_gc = 0

for record in biometal.FastqStream.from_path("reads.fq.gz"):
    total_reads += 1

    # ARM-optimized operations (16-25× faster)
    if record.mean_quality() >= 30:
        high_quality += 1

    gc = biometal.gc_content(record.sequence)
    total_gc += gc

# Results
print(f"High quality: {high_quality/total_reads*100:.1f}%")
print(f"Mean GC: {total_gc/total_reads:.2f}%")
```

**Performance**: On Apple M1 Max, processes 1 GB FASTQ in ~2-3 seconds with constant ~5 MB memory.

### BAM Alignment Analysis

**Bottleneck**: File decompression, sequence decoding, CIGAR parsing

**Optimizations**:
1. ✅ Use indexed queries for targeted analysis (1.68-500× speedup)
2. ✅ Stream, don't accumulate
3. ✅ ARM NEON sequence decoding (4.62× speedup in v1.5.0)
4. ⏳ Parallel BGZF decompression (6.5× speedup, coming in Phase 2)

**Example (Targeted Analysis)**:
```python
import biometal

# Load BAI index once
index = biometal.BaiIndex.from_path("alignments.bam.bai")

# Query specific regions (O(log n))
regions = [("chr1", 1e6, 2e6), ("chr1", 5e6, 6e6), ...]

for ref, start, end in regions:
    count = 0
    for record in biometal.BamReader.query_region(
        "alignments.bam", index, ref, start, end
    ):
        if record.is_mapped and record.mapq >= 30:
            count += 1

    print(f"{ref}:{start}-{end}: {count} reads")
```

**Performance**:
- Current: 55.1 MiB/s (BGZF + NEON sequence decoding)
- Phase 2: ~350 MiB/s (+ parallel BGZF decompression)

**Example (Full Genome)**:
```python
import biometal

# Sequential read for genome-wide stats
total_mapped = 0
total_unmapped = 0

for record in biometal.BamReader.from_path("alignments.bam"):
    if record.is_mapped:
        total_mapped += 1
    else:
        total_unmapped += 1

print(f"Mapping rate: {total_mapped/(total_mapped+total_unmapped)*100:.2f}%")
```

### K-mer Analysis

**Bottleneck**: K-mer extraction and counting

**Optimizations**:
1. ✅ Stream, don't accumulate
2. ✅ Use Counter for efficient counting
3. ⚠️ K-mer operations not yet NEON-optimized (planned)

**Example**:
```python
import biometal
from collections import Counter

# Streaming k-mer extraction
k = 21
kmer_counts = Counter()

for record in biometal.FastqStream.from_path("reads.fq.gz"):
    kmers = biometal.extract_kmers(record.sequence, k)
    kmer_counts.update(kmers)

# Top 10 most common k-mers
for kmer, count in kmer_counts.most_common(10):
    print(f"{kmer.decode()}: {count:,}")
```

**Performance**: Constant ~5 MB memory, processes 1 GB FASTQ in ~5-10 seconds.

### Network Streaming

**Bottleneck**: Network I/O (not CPU)

**Optimizations**:
1. ✅ Use batch processing for better throughput
2. ✅ Enable prefetching (automatic in biometal)
3. ⚠️ Consider connection pooling for multiple URLs

**Example**:
```python
import biometal

# Stream from HTTP (no download)
url = "https://example.com/reads.fastq.gz"

batch = []
for record in biometal.FastqStream.from_url(url):
    batch.append(record)

    if len(batch) >= 1000:
        # Process batch
        process_batch(batch)
        batch = []

# Process remaining
if batch:
    process_batch(batch)
```

**Performance**: ~60-80 MB/s throughput on 1 Gbps connection, constant ~5 MB memory.

---

## Performance Profiling

### Measuring Your Performance

**1. Time Operations**:
```python
import time

start = time.time()
# Your code here
elapsed = time.time() - start
print(f"Elapsed: {elapsed:.2f} seconds")
```

**2. Memory Usage**:
```python
import tracemalloc

tracemalloc.start()
# Your code here
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()

print(f"Current: {current / 1024 / 1024:.1f} MB")
print(f"Peak: {peak / 1024 / 1024:.1f} MB")
```

**3. Throughput**:
```python
import time

start = time.time()
record_count = 0

for record in biometal.FastqStream.from_path("reads.fq.gz"):
    record_count += 1

elapsed = time.time() - start
throughput = record_count / elapsed

print(f"Throughput: {throughput:,.0f} records/sec")
```

### Expected Performance Targets

**FASTQ Processing** (Apple M1 Max):
- Base counting: 5,254,000 sequences/sec
- GC content: 5,954,000 sequences/sec
- Quality filter: 6,143,000 sequences/sec
- Memory: ~5 MB constant

**BAM Processing** (Apple M1 Max):
- Sequential read: 55.1 MiB/s (compressed file)
- Record throughput: 5.82 million records/sec
- Indexed query: 1.68-500× faster than full scan
- Memory: ~5 MB constant

**Network Streaming** (1 Gbps):
- FASTQ from URL: 60-80 MB/s
- Latency: <100 ms to first record
- Memory: ~5 MB constant

### Profiling Tools

**Python**:
```bash
# Time profiling
python -m cProfile -s cumulative your_script.py

# Line profiling
pip install line_profiler
kernprof -l -v your_script.py

# Memory profiling
pip install memory_profiler
python -m memory_profiler your_script.py
```

**Rust**:
```bash
# Criterion benchmarks
cargo bench

# Flamegraph profiling
cargo install flamegraph
cargo flamegraph --bin your_binary
```

---

## Common Anti-Patterns

### Anti-Pattern 1: Accumulating Records

**Problem**: Loads entire dataset into memory

```python
# ❌ BAD: Accumulates 1,344 MB for 1M sequences
records = []
for record in biometal.FastqStream.from_path("reads.fq.gz"):
    records.append(record)

# Then process all at once
for record in records:
    process(record)
```

**Solution**: Process immediately

```python
# ✅ GOOD: Constant 5 MB memory
for record in biometal.FastqStream.from_path("reads.fq.gz"):
    process(record)  # Immediate processing
```

**Impact**: 99.5% memory reduction, enables processing of arbitrarily large files

### Anti-Pattern 2: Manual Decompression

**Problem**: Adds overhead, disables optimizations

```python
# ❌ BAD: Manual gzip decompression
import gzip
with gzip.open("reads.fastq.gz", 'rt') as f:
    for line in f:
        # Parse manually, slower
        pass
```

**Solution**: Let biometal handle it

```python
# ✅ GOOD: Automatic decompression with optimizations
for record in biometal.FastqStream.from_path("reads.fastq.gz"):
    process(record)
```

**Impact**: Automatic parallel BGZF decompression (when available)

### Anti-Pattern 3: Full Scan Instead of Indexed Query

**Problem**: O(n) instead of O(log n) for BAM region queries

```python
# ❌ SLOW: Full scan (O(n))
for record in biometal.BamReader.from_path("file.bam"):
    if (record.reference_id == 0 and
        record.position and
        1000000 <= record.position < 2000000):
        process(record)
```

**Solution**: Use indexed query

```python
# ✅ FAST: Indexed query (O(log n))
index = biometal.BaiIndex.from_path("file.bam.bai")
for record in biometal.BamReader.query_region(
    "file.bam", index, "chr1", 1000000, 2000000
):
    process(record)
```

**Impact**: 1.68-500× speedup (increases with file size)

### Anti-Pattern 4: Not Using ARM-Optimized Operations

**Problem**: Reimplements operations that biometal optimizes

```python
# ❌ SLOW: Manual base counting (scalar, slow)
def count_bases_manual(sequence):
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for base in sequence:
        if base in counts:
            counts[base] += 1
    return counts

for record in stream:
    counts = count_bases_manual(record.sequence)
```

**Solution**: Use biometal's optimized operations

```python
# ✅ FAST: ARM NEON optimized (16.7× speedup)
for record in stream:
    counts = biometal.count_bases(record.sequence)
```

**Impact**: 16-25× speedup on ARM hardware

### Anti-Pattern 5: Loading Index Repeatedly

**Problem**: Wastes time loading index for each query

```python
# ❌ SLOW: Loads index N times
for ref, start, end in regions:
    index = biometal.BaiIndex.from_path("file.bam.bai")  # Reload!
    for record in biometal.BamReader.query_region(...):
        process(record)
```

**Solution**: Load once, reuse

```python
# ✅ FAST: Load index once
index = biometal.BaiIndex.from_path("file.bam.bai")
for ref, start, end in regions:
    for record in biometal.BamReader.query_region(...):
        process(record)
```

**Impact**: Eliminates N × 1ms overhead (negligible per query, but adds up)

---

## Benchmarking Your Code

### Creating a Benchmark

**Python**:
```python
import time
import biometal

def benchmark_operation(file_path, num_runs=5):
    """Benchmark an operation with multiple runs."""
    times = []

    for _ in range(num_runs):
        start = time.time()

        # Your operation here
        count = 0
        for record in biometal.FastqStream.from_path(file_path):
            count += 1

        elapsed = time.time() - start
        times.append(elapsed)

    # Statistics
    mean_time = sum(times) / len(times)
    min_time = min(times)
    max_time = max(times)

    print(f"Mean: {mean_time:.2f}s")
    print(f"Min: {min_time:.2f}s")
    print(f"Max: {max_time:.2f}s")
    print(f"Throughput: {count/mean_time:,.0f} records/sec")

# Run benchmark
benchmark_operation("reads.fastq.gz")
```

**Rust** (Criterion):
```rust
use criterion::{criterion_group, criterion_main, Criterion};
use biometal::FastqStream;

fn bench_fastq_parsing(c: &mut Criterion) {
    c.bench_function("fastq_parse_1m", |b| {
        b.iter(|| {
            let stream = FastqStream::from_path("reads.fq.gz").unwrap();
            let count: usize = stream.count();
            count
        })
    });
}

criterion_group!(benches, bench_fastq_parsing);
criterion_main!(benches);
```

### Comparing Implementations

```python
import time
import biometal

def benchmark_comparison(file_path):
    """Compare two implementations."""

    # Implementation 1: ARM-optimized
    start = time.time()
    count1 = 0
    gc_sum1 = 0
    for record in biometal.FastqStream.from_path(file_path):
        count1 += 1
        gc_sum1 += biometal.gc_content(record.sequence)
    time1 = time.time() - start

    # Implementation 2: Manual (scalar)
    start = time.time()
    count2 = 0
    gc_sum2 = 0
    for record in biometal.FastqStream.from_path(file_path):
        count2 += 1
        # Manual GC calculation (slower)
        seq = record.sequence.decode()
        gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
        gc_sum2 += gc
    time2 = time.time() - start

    # Results
    print(f"ARM-optimized: {time1:.2f}s ({count1/time1:,.0f} records/sec)")
    print(f"Manual scalar: {time2:.2f}s ({count2/time2:,.0f} records/sec)")
    print(f"Speedup: {time2/time1:.2f}×")

benchmark_comparison("reads.fastq.gz")
```

---

## Summary

### Performance Checklist

When optimizing biometal code, check:

- [ ] Using ARM hardware when possible (16-25× speedup)
- [ ] Streaming, not accumulating (99.5% memory reduction)
- [ ] Using indexed queries for BAM when appropriate (1.68-500× speedup)
- [ ] Letting biometal handle decompression (automatic optimizations)
- [ ] Using biometal's optimized operations (count_bases, gc_content, etc.)
- [ ] Loading indexes once and reusing (eliminates repeated overhead)
- [ ] Batch processing for network streaming (better throughput)
- [ ] Profiling to identify actual bottlenecks (measure, don't guess)

### Performance Priorities

1. **Architecture First**: Use streaming patterns (Rule 5)
2. **Hardware Next**: Use ARM if compute-bound (Rule 1)
3. **Algorithm Then**: Use indexed queries for targeted analysis
4. **Micro-optimizations Last**: Profile before optimizing details

### Getting Help

- **User Guide**: [docs/USER_GUIDE.md](USER_GUIDE.md)
- **Architecture**: [docs/ARCHITECTURE.md](ARCHITECTURE.md)
- **Optimization Rules**: [OPTIMIZATION_RULES.md](../OPTIMIZATION_RULES.md)
- **GitHub Issues**: Report performance issues with benchmarks

---

**Version**: 1.6.0
**Last Updated**: November 10, 2025
**Author**: biometal team

For questions or feedback, open an issue on [GitHub](https://github.com/scotthandley/biometal).
