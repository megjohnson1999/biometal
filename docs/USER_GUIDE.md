# biometal User Guide

**Version**: 1.6.0
**Last Updated**: November 10, 2025

Welcome to biometal! This guide will help you get started with high-performance bioinformatics on ARM and x86_64 platforms.

---

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Core Concepts](#core-concepts)
4. [Quick Start](#quick-start)
5. [Common Workflows](#common-workflows)
6. [Performance Guide](#performance-guide)
7. [Troubleshooting](#troubleshooting)
8. [Migration Guide](#migration-guide)
9. [API Reference](#api-reference)

---

## Introduction

### What is biometal?

biometal is a high-performance bioinformatics library designed for:
- **ARM-native performance**: 16-25× speedup on Apple Silicon (M1/M2/M3/M4) and AWS Graviton
- **Constant memory usage**: ~5 MB regardless of dataset size (streaming architecture)
- **Network streaming**: Analyze data without downloading
- **Production quality**: Robust error handling, comprehensive testing, Python and Rust APIs

### Why Use biometal?

**For Python Users**:
- Drop-in replacement for common pysam/samtools operations
- 16-25× faster on ARM processors (Mac M-series, AWS Graviton)
- Constant memory usage (no OOM errors on large files)
- Simple installation: `pip install biometal-rs`

**For Rust Users**:
- High-performance parsers for FASTQ, FASTA, BAM/SAM
- Zero-copy parsing where possible
- ARM NEON SIMD optimizations with x86_64 fallbacks
- Network streaming support built-in

**Use Cases**:
- Quality control on FASTQ files
- BAM file analysis and filtering
- Indexed region queries (O(log n) vs O(n))
- K-mer analysis and minimizer extraction
- Large-scale dataset processing
- Network streaming from public databases

### Supported Platforms

| Platform | SIMD Support | Speedup | Status |
|----------|-------------|---------|--------|
| **Mac ARM** (M1/M2/M3/M4) | ARM NEON | 16-25× | ✅ Optimized |
| **Linux ARM** (Graviton) | ARM NEON | 6-10× | ✅ Portable |
| **x86_64** (Intel/AMD) | None | 1× | ✅ Fallback |

---

## Installation

### Python Installation

**From PyPI** (Recommended):
```bash
pip install biometal-rs
```

**Verify Installation**:
```python
import biometal
print(biometal.__version__)  # Should print "1.6.0" or higher
```

**Platform-Specific Notes**:

- **macOS (ARM)**: Pre-built wheels available, no compilation needed
- **macOS (Intel)**: Pre-built wheels available, scalar fallback
- **Linux (ARM/x86_64)**: May compile from source (requires Rust toolchain)
- **Windows**: Limited support, use WSL2 recommended

**Troubleshooting Installation**:
```bash
# If pip install fails, try upgrading pip first
pip install --upgrade pip

# For Linux ARM/x86_64, install Rust if compilation is needed
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

### Rust Installation

Add to `Cargo.toml`:
```toml
[dependencies]
biometal = "1.6.0"
```

Or via command line:
```bash
cargo add biometal
```

**Feature Flags**:
```toml
[dependencies]
biometal = { version = "1.6.0", features = ["network", "sra"] }
```

Available features:
- `network`: HTTP streaming support
- `sra`: SRA toolkit integration (requires `fasterq-dump`)

---

## Core Concepts

### 1. Streaming Architecture

biometal uses **streaming iterators** instead of loading entire files into memory.

**Traditional Approach** (BAD):
```python
# ❌ Loads entire file into memory
records = []
for record in reader:
    records.append(record)  # Accumulates in RAM

# If file is 50 GB, you need >50 GB RAM
process_all(records)
```

**biometal Approach** (GOOD):
```python
# ✅ Constant ~5 MB memory
for record in biometal.FastqStream.from_path("reads.fastq.gz"):
    # Process each record immediately
    if record.mean_quality() >= 30:
        process(record)
    # Record is discarded after processing
```

**Benefits**:
- Constant memory (~5 MB) regardless of file size
- Can process 5 TB datasets on consumer hardware
- Enables network streaming (analyze without downloading)

### 2. ARM NEON Optimization

biometal automatically uses ARM NEON SIMD instructions when available.

**What is NEON?**
- ARM's SIMD (Single Instruction, Multiple Data) instruction set
- Processes 16 bytes in parallel per CPU cycle
- Available on Apple Silicon (M1/M2/M3/M4) and AWS Graviton

**Performance Impact**:
```python
# On Apple M1 Max:
# Base counting with NEON: 16-25× faster than scalar
# Quality filtering: 10-15× faster
# GC content calculation: 18-22× faster
```

**No Code Changes Required**:
```python
# This automatically uses NEON on ARM platforms
counts = biometal.count_bases(sequence)
# On ARM: Uses NEON (16× faster)
# On x86_64: Uses scalar fallback (1× baseline)
```

### 3. Zero-Copy Philosophy

biometal minimizes memory allocations and copies.

**Example: Sequence Access**:
```python
# Returns view into existing buffer (no copy)
sequence = record.sequence  # bytes object, zero-copy

# Only allocates when conversion is needed
sequence_str = sequence.decode('utf-8')  # Allocation here
```

**Impact**:
- Lower memory pressure
- Reduced GC overhead (Python)
- Better cache locality

### 4. Error Handling

biometal uses Rust's `Result` type for robust error handling.

**Python**: Exceptions with detailed error messages
```python
try:
    reader = biometal.FastqStream.from_path("missing.fastq")
except FileNotFoundError as e:
    print(f"Error: {e}")
```

**Rust**: Result types
```rust
use biometal::io::FastqStream;

match FastqStream::from_path("reads.fastq.gz") {
    Ok(reader) => { /* process */ },
    Err(e) => eprintln!("Error: {}", e),
}
```

---

## Quick Start

### Python Quick Start

**1. Read FASTQ File**:
```python
import biometal

# Streaming FASTQ reader (gzip-compressed or plain text)
for record in biometal.FastqStream.from_path("reads.fastq.gz"):
    print(f"Name: {record.name}")
    print(f"Sequence: {record.sequence.decode()}")
    print(f"Quality: {record.quality.decode()}")
    print(f"Mean Quality: {record.mean_quality():.1f}")
    break  # Just show first record
```

**2. Filter by Quality**:
```python
import biometal

high_quality_count = 0
for record in biometal.FastqStream.from_path("reads.fastq.gz"):
    if record.mean_quality() >= 30:
        high_quality_count += 1

print(f"High-quality reads: {high_quality_count}")
```

**3. Count Bases (ARM-Optimized)**:
```python
import biometal

total_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

for record in biometal.FastqStream.from_path("reads.fastq.gz"):
    counts = biometal.count_bases(record.sequence)
    total_counts['A'] += counts['A']
    total_counts['C'] += counts['C']
    total_counts['G'] += counts['G']
    total_counts['T'] += counts['T']

print(f"Base composition: {total_counts}")
```

**4. Read BAM File**:
```python
import biometal

# Sequential BAM reading (streaming)
for record in biometal.BamReader.from_path("alignments.bam"):
    if record.is_mapped and record.mapq >= 30:
        print(f"{record.name}: chr{record.reference_id} @ {record.position}")
```

**5. Indexed BAM Query** (NEW in v1.6.0):
```python
import biometal

# Load BAI index
index = biometal.BaiIndex.from_path("alignments.bam.bai")

# Query specific region (1.68-500× faster than full scan)
for record in biometal.BamReader.query_region(
    "alignments.bam",
    index,
    "chr1",
    1000000,  # start
    2000000   # end
):
    if record.is_mapped and record.mapq >= 30:
        print(f"{record.name}: {record.position}")
```

### Rust Quick Start

**1. Read FASTQ File**:
```rust
use biometal::io::FastqStream;

let mut reader = FastqStream::from_path("reads.fastq.gz")?;

for result in reader {
    let record = result?;
    println!("Name: {}", record.name);
    println!("Sequence: {}", std::str::from_utf8(record.sequence)?);
    println!("Quality: {}", std::str::from_utf8(record.quality)?);
    break; // Just show first record
}
```

**2. Count Bases (ARM-Optimized)**:
```rust
use biometal::operations::count_bases;

let sequence = b"ACGTACGT";
let counts = count_bases(sequence);
println!("A: {}, C: {}, G: {}, T: {}", counts.a, counts.c, counts.g, counts.t);
```

**3. Read BAM File**:
```rust
use biometal::io::bam::BamReader;

let mut reader = BamReader::from_path("alignments.bam")?;

for result in reader {
    let record = result?;
    if record.is_mapped() && record.mapq().unwrap_or(0) >= 30 {
        println!("{}: ref {} @ {}", record.name(), record.reference_id(), record.position());
    }
}
```

---

## Common Workflows

### Workflow 1: FASTQ Quality Control

**Goal**: Analyze FASTQ file quality metrics

```python
import biometal
from collections import defaultdict

# Initialize counters
total_reads = 0
total_bases = 0
quality_histogram = defaultdict(int)
length_histogram = defaultdict(int)
base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}

# Stream through file (constant memory)
for record in biometal.FastqStream.from_path("reads.fastq.gz"):
    total_reads += 1
    total_bases += len(record.sequence)

    # Quality metrics
    mean_q = record.mean_quality()
    quality_histogram[int(mean_q)] += 1

    # Length distribution
    length_histogram[len(record.sequence)] += 1

    # Base composition (ARM-optimized)
    counts = biometal.count_bases(record.sequence)
    base_counts['A'] += counts['A']
    base_counts['C'] += counts['C']
    base_counts['G'] += counts['G']
    base_counts['T'] += counts['T']

# Calculate statistics
gc_content = (base_counts['G'] + base_counts['C']) / total_bases * 100
mean_length = total_bases / total_reads

# Print report
print(f"Total reads: {total_reads:,}")
print(f"Total bases: {total_bases:,}")
print(f"Mean length: {mean_length:.1f}")
print(f"GC content: {gc_content:.2f}%")
print(f"\nBase composition:")
for base, count in sorted(base_counts.items()):
    pct = count / total_bases * 100
    print(f"  {base}: {count:>12,} ({pct:5.2f}%)")

# Quality distribution
print(f"\nQuality distribution:")
for q in sorted(quality_histogram.keys()):
    count = quality_histogram[q]
    print(f"  Q{q:2d}: {count:>8,} reads")
```

**Performance**: On Apple M1 Max, processes 1 GB FASTQ in ~2-3 seconds with constant ~5 MB memory.

### Workflow 2: BAM File Analysis

**Goal**: Analyze alignment statistics from BAM file

```python
import biometal
from collections import defaultdict

# Initialize counters
total_reads = 0
mapped_reads = 0
unmapped_reads = 0
mapq_histogram = defaultdict(int)
coverage = defaultdict(int)

# Stream through BAM file (constant memory)
for record in biometal.BamReader.from_path("alignments.bam"):
    total_reads += 1

    if record.is_mapped:
        mapped_reads += 1

        # MAPQ distribution
        if record.mapq is not None:
            mapq_histogram[record.mapq] += 1

        # Coverage calculation (chr1 only, for example)
        if record.reference_id == 0 and record.position is not None:
            pos = record.position
            for op in record.cigar:
                if op.consumes_reference():
                    for i in range(op.length):
                        coverage[pos] += 1
                        pos += 1
    else:
        unmapped_reads += 1

# Calculate statistics
mapping_rate = mapped_reads / total_reads * 100 if total_reads > 0 else 0
mean_coverage = sum(coverage.values()) / len(coverage) if coverage else 0

# Print report
print(f"Total reads: {total_reads:,}")
print(f"Mapped: {mapped_reads:,} ({mapping_rate:.2f}%)")
print(f"Unmapped: {unmapped_reads:,}")
print(f"\nMAPQ distribution:")
for q in sorted(mapq_histogram.keys())[-10:]:  # Top 10
    count = mapq_histogram[q]
    print(f"  Q{q:2d}: {count:>8,} reads")

if coverage:
    print(f"\nChr1 mean coverage: {mean_coverage:.2f}×")
```

### Workflow 3: Indexed BAM Region Query (v1.6.0)

**Goal**: Fast extraction of reads from specific genomic region

```python
import biometal

# Load BAI index once
index = biometal.BaiIndex.from_path("alignments.bam.bai")

# Define regions of interest
regions = [
    ("chr1", 1000000, 2000000),
    ("chr1", 5000000, 6000000),
    ("chr2", 100000, 200000),
]

# Query each region (O(log n) vs O(n) for full scan)
for ref_name, start, end in regions:
    high_quality_count = 0
    total_bases = 0
    coverage = {}

    for record in biometal.BamReader.query_region(
        "alignments.bam",
        index,
        ref_name,
        start,
        end
    ):
        if record.is_mapped and record.mapq >= 30:
            high_quality_count += 1
            total_bases += len(record.sequence)

            # Track coverage
            if record.position is not None:
                pos = record.position
                for op in record.cigar:
                    if op.consumes_reference():
                        for i in range(op.length):
                            coverage[pos] = coverage.get(pos, 0) + 1
                            pos += 1

    # Calculate statistics
    region_size = end - start
    mean_cov = sum(coverage.values()) / len(coverage) if coverage else 0

    print(f"\nRegion: {ref_name}:{start:,}-{end:,}")
    print(f"  High-quality reads: {high_quality_count:,}")
    print(f"  Total bases: {total_bases:,}")
    print(f"  Mean coverage: {mean_cov:.2f}×")
```

**Performance**: On 1 GB BAM file, indexed queries are 1.68× faster than full scan for small files, and 10-500× faster for large files (1-10 GB).

**When to Use Indexed Queries**:
- ✅ Extracting specific regions (exons, genes, intervals)
- ✅ Coverage analysis for targeted regions
- ✅ Variant calling on specific loci
- ✅ Multi-sample region comparison
- ❌ Full-genome analysis (use sequential read instead)

### Workflow 4: K-mer Analysis

**Goal**: Extract k-mers and count frequencies

```python
import biometal
from collections import Counter

# K-mer extraction (k=21)
k = 21
kmer_counts = Counter()

for record in biometal.FastqStream.from_path("reads.fastq.gz"):
    # Extract k-mers from sequence
    kmers = biometal.extract_kmers(record.sequence, k)
    kmer_counts.update(kmers)

# Top 10 most common k-mers
print(f"Top 10 most common {k}-mers:")
for kmer, count in kmer_counts.most_common(10):
    print(f"  {kmer.decode()}: {count:,}")

# K-mer spectrum analysis
histogram = Counter(kmer_counts.values())
print(f"\nK-mer frequency distribution:")
for freq in sorted(histogram.keys())[:20]:  # First 20 frequencies
    count = histogram[freq]
    print(f"  {freq:>3}×: {count:>8,} k-mers")
```

### Workflow 5: Network Streaming

**Goal**: Analyze data without downloading

```python
import biometal

# Stream directly from HTTP URL (no download)
url = "https://example.com/reads.fastq.gz"

total_reads = 0
high_quality_reads = 0

for record in biometal.FastqStream.from_url(url):
    total_reads += 1
    if record.mean_quality() >= 30:
        high_quality_reads += 1

    # Stop after analyzing first 10,000 reads
    if total_reads >= 10000:
        break

hq_rate = high_quality_reads / total_reads * 100
print(f"High-quality rate: {hq_rate:.2f}% (n={total_reads:,})")
```

**Benefits**:
- No disk space required
- Start analysis immediately (no download wait)
- Sample large files without downloading entirely
- Ideal for exploratory analysis

### Workflow 6: Quality Filtering and Output

**Goal**: Filter reads by quality and save to new file

```python
import biometal

# Create output writer
output = biometal.FastqSink.create("filtered_reads.fastq.gz")

total_reads = 0
kept_reads = 0

# Filter by quality threshold
for record in biometal.FastqStream.from_path("reads.fastq.gz"):
    total_reads += 1

    if record.mean_quality() >= 30 and len(record.sequence) >= 50:
        output.write(record)
        kept_reads += 1

# Close output
output.close()

retention_rate = kept_reads / total_reads * 100
print(f"Kept {kept_reads:,} / {total_reads:,} reads ({retention_rate:.2f}%)")
```

---

## Performance Guide

### Getting Best Performance

**1. Use ARM Hardware When Possible**
- Apple Silicon (M1/M2/M3/M4): 16-25× speedup
- AWS Graviton (ARM): 6-10× speedup
- Intel/AMD: 1× baseline (scalar fallback)

**2. Stream, Don't Accumulate**
```python
# ❌ BAD: Accumulates in memory
records = []
for record in reader:
    records.append(record)

# ✅ GOOD: Constant memory
for record in reader:
    process(record)  # Immediate processing
```

**3. Use Indexed Queries for BAM**
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

**4. Let biometal Handle Decompression**
```python
# ✅ GOOD: biometal handles .gz automatically
for record in biometal.FastqStream.from_path("reads.fastq.gz"):
    process(record)

# ❌ BAD: Manual decompression adds overhead
import gzip
with gzip.open("reads.fastq.gz") as f:
    for line in f:  # Slower, more memory
        ...
```

**5. Batch Processing for Network Streaming**
```python
# ✅ Better throughput with batching
batch = []
for record in biometal.FastqStream.from_url(url):
    batch.append(record)
    if len(batch) >= 1000:
        process_batch(batch)
        batch = []
```

### Performance Benchmarks

**FASTQ Processing** (1 GB file, Apple M1 Max):
```
Base counting:    16-25× faster than pure Python
Quality filter:   10-15× faster than pure Python
GC content:       18-22× faster than pure Python
Memory:           ~5 MB constant (vs 1+ GB for batch loading)
```

**BAM Processing** (100K records, Apple M1 Max):
```
Sequential read:  ~180 ms (full parsing)
Indexed query:    ~11 ms (small region)
Speedup:          1.68× (small file) to 500× (large file)
```

**Network Streaming** (HTTP, 1 Gbps connection):
```
FASTQ from URL:   ~60-80 MB/s throughput
No local storage: 0 GB disk usage
Start latency:    <100 ms
```

---

## Troubleshooting

### Common Issues

#### Issue 1: Import Error in Python

**Symptom**:
```python
import biometal
# ImportError: No module named 'biometal'
```

**Solution**:
```bash
# Check if installed
pip list | grep biometal

# If not installed
pip install biometal-rs

# If installed but not found, check Python version
python --version  # Should be 3.8+
pip --version     # Should match Python version
```

#### Issue 2: Performance Not as Expected

**Symptom**: No speedup on ARM hardware

**Diagnosis**:
```python
import biometal
print(biometal.get_optimization_info())
# Should show "ARM NEON: enabled" on ARM platforms
```

**Solution**:
- If NEON not detected, check you're on ARM platform
- Verify with: `uname -m` (should show `arm64` or `aarch64`)
- Reinstall with: `pip install --force-reinstall biometal-rs`

#### Issue 3: File Not Found

**Symptom**:
```python
reader = biometal.FastqStream.from_path("reads.fastq.gz")
# FileNotFoundError: File not found: reads.fastq.gz
```

**Solution**:
```python
# Use absolute path
import os
path = os.path.abspath("reads.fastq.gz")
reader = biometal.FastqStream.from_path(path)

# Or check current directory
print(os.getcwd())
print(os.listdir('.'))
```

#### Issue 4: Memory Issues (Large Files)

**Symptom**: OOM (Out of Memory) errors

**Diagnosis**:
```python
# ❌ BAD: Are you accumulating records?
records = []
for record in reader:
    records.append(record)  # DON'T DO THIS
```

**Solution**:
```python
# ✅ GOOD: Stream and process immediately
for record in reader:
    process(record)  # Process and discard
    # No accumulation = constant memory
```

#### Issue 5: BAI Index Not Found

**Symptom**:
```python
index = biometal.BaiIndex.from_path("file.bam.bai")
# FileNotFoundError
```

**Solution**:
```bash
# Generate BAI index using samtools
samtools index file.bam
# This creates file.bam.bai

# Or in Python
import biometal
index = biometal.BaiIndex.from_path("file.bam.bai")
```

#### Issue 6: Slow Network Streaming

**Symptom**: Network streaming slower than expected

**Diagnosis**:
- Check network speed: `curl -o /dev/null [url]`
- Check server limits: Some servers throttle connections

**Solution**:
```python
# Use connection pooling for multiple requests
# Enable keep-alive connections
# Consider downloading large files locally if repeated access needed
```

### Getting Help

**Bug Reports**: [GitHub Issues](https://github.com/scotthandley/biometal/issues)

**Questions**:
- GitHub Discussions
- Biostars (tag: `biometal`)

**Documentation**:
- API Reference: See inline documentation
- Examples: `examples/` directory
- Architecture: `docs/ARCHITECTURE.md`

---

## Migration Guide

### From pysam

**pysam** is a Python wrapper around HTSlib. biometal offers similar functionality with better performance on ARM.

#### Reading BAM Files

**pysam**:
```python
import pysam

# Sequential reading
with pysam.AlignmentFile("file.bam", "rb") as f:
    for read in f:
        if read.is_mapped and read.mapping_quality >= 30:
            print(f"{read.query_name}: {read.reference_start}")
```

**biometal**:
```python
import biometal

# Sequential reading
for record in biometal.BamReader.from_path("file.bam"):
    if record.is_mapped and record.mapq >= 30:
        print(f"{record.name}: {record.position}")
```

#### Indexed Region Queries

**pysam**:
```python
import pysam

with pysam.AlignmentFile("file.bam", "rb") as f:
    for read in f.fetch("chr1", 1000000, 2000000):
        if read.mapping_quality >= 30:
            process(read)
```

**biometal**:
```python
import biometal

index = biometal.BaiIndex.from_path("file.bam.bai")
for record in biometal.BamReader.query_region(
    "file.bam", index, "chr1", 1000000, 2000000
):
    if record.mapq >= 30:
        process(record)
```

#### Key Differences

| Feature | pysam | biometal |
|---------|-------|----------|
| **Platform** | x86_64 optimized | ARM optimized (16-25×) |
| **Memory** | Variable | Constant (~5 MB) |
| **API** | Context managers | Iterators |
| **Installation** | Requires compilation | Pre-built wheels (Python) |
| **Formats** | BAM/SAM/CRAM/VCF | BAM/SAM/FASTQ/FASTA |

### From samtools

**samtools** is a C-based command-line toolkit. biometal offers programmatic API with better ARM performance.

#### Viewing BAM Records

**samtools**:
```bash
samtools view file.bam | head -10
```

**biometal**:
```python
import biometal

for i, record in enumerate(biometal.BamReader.from_path("file.bam")):
    print(f"{record.name}\t{record.sequence.decode()}")
    if i >= 9:
        break
```

#### Region Queries

**samtools**:
```bash
samtools view file.bam chr1:1000000-2000000
```

**biometal**:
```python
import biometal

index = biometal.BaiIndex.from_path("file.bam.bai")
for record in biometal.BamReader.query_region(
    "file.bam", index, "chr1", 1000000, 2000000
):
    print(f"{record.name}\t{record.position}")
```

#### Flagstat

**samtools**:
```bash
samtools flagstat file.bam
```

**biometal**:
```python
import biometal

total = mapped = unmapped = 0
for record in biometal.BamReader.from_path("file.bam"):
    total += 1
    if record.is_mapped:
        mapped += 1
    else:
        unmapped += 1

print(f"{total} total")
print(f"{mapped} mapped ({mapped/total*100:.2f}%)")
print(f"{unmapped} unmapped ({unmapped/total*100:.2f}%)")
```

### Migration Checklist

- [ ] Install biometal: `pip install biometal-rs`
- [ ] Replace `pysam.AlignmentFile` with `biometal.BamReader`
- [ ] Replace `pysam.FastxFile` with `biometal.FastqStream`
- [ ] Replace `.fetch()` with `.query_region()`
- [ ] Remove context managers (`with`), use iterators
- [ ] Test on ARM hardware for performance gains
- [ ] Update documentation/comments

---

## API Reference

### Python API

**Core Modules**:
- `biometal.FastqStream` - FASTQ streaming reader
- `biometal.FastaStream` - FASTA streaming reader
- `biometal.BamReader` - BAM/SAM streaming reader
- `biometal.BaiIndex` - BAI index support (v1.6.0)
- `biometal.FastqSink` - FASTQ writer
- `biometal.count_bases()` - Base counting (ARM-optimized)
- `biometal.gc_content()` - GC content (ARM-optimized)
- `biometal.reverse_complement()` - Sequence reversal
- `biometal.extract_kmers()` - K-mer extraction

**Detailed Documentation**:
```python
import biometal
help(biometal.BamReader)  # Inline documentation
```

### Rust API

**Core Modules**:
- `biometal::io::FastqStream` - FASTQ parsing
- `biometal::io::FastaStream` - FASTA parsing
- `biometal::io::bam::BamReader` - BAM parsing
- `biometal::io::bam::BaiIndex` - BAI index support
- `biometal::operations::count_bases` - Base counting
- `biometal::operations::gc_content` - GC content
- `biometal::operations::reverse_complement` - Sequence ops

**Documentation**:
```bash
cargo doc --open
```

---

## Next Steps

### Learn More

- **Architecture**: See `docs/ARCHITECTURE.md` for design details
- **Performance**: See `docs/PERFORMANCE_TUNING.md` for optimization rules
- **Python**: See `docs/PYTHON.md` for Python-specific details
- **BAM**: See `docs/BAM_API.md` for BAM format details

### Examples

Explore the `examples/` directory:
- `examples/basic_fastq_streaming.rs` - FASTQ streaming basics
- `examples/bam_streaming.rs` - BAM parsing examples
- `examples/bam_advanced_filtering.py` - BAM filtering patterns
- `examples/python_bai_usage.py` - Indexed query examples (v1.6.0)
- `examples/http_streaming.rs` - Network streaming
- `examples/kmer_operations_full.rs` - K-mer analysis

### Contributing

biometal is open-source and welcomes contributions:
- Bug reports: GitHub Issues
- Feature requests: GitHub Discussions
- Code contributions: Pull requests welcome

### Stay Updated

- **GitHub**: Watch repository for updates
- **Changelog**: `CHANGELOG.md` for version history
- **Releases**: GitHub releases for new versions

---

**Version**: 1.6.0
**Last Updated**: November 10, 2025
**License**: MIT
**Author**: Scott Handley

For questions or feedback, open an issue on [GitHub](https://github.com/scotthandley/biometal).
