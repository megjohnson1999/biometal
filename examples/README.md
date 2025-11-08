# biometal Examples

Practical code examples demonstrating biometal's streaming architecture and ARM-native performance.

---

## Python Examples

### Basic Usage

#### [basic_usage.py](basic_usage.py)
Introductory examples covering core functionality:
- Streaming FASTQ files with constant memory
- GC content and base counting (16-25× ARM NEON speedup)
- Quality score analysis
- Network streaming basics

**Run**:
```bash
python basic_usage.py
```

### BAM Alignment Analysis

#### [bam_advanced_filtering.py](bam_advanced_filtering.py) ⭐ NEW
Advanced BAM filtering techniques for production workflows:

**Features**:
- **Region queries**: Extract alignments from chr:start-end
- **Complex filters**: Combine MAPQ, flags, strand, length criteria
- **Helper classes**: `Region`, `FilterCriteria` for reusable filters
- **Common workflows**: Variant calling, RNA-seq, coverage analysis
- **Statistics**: Flag distributions, MAPQ histograms, coverage

**Components**:
```python
# Region Query Helper
region = Region(reference_id=0, start=1000, end=2000)
for record in filter_by_region("alignments.bam", region, min_mapq=30):
    process(record)

# Complex Filter
criteria = FilterCriteria(
    min_mapq=30,
    require_primary=True,
    require_forward=True,
    min_length=50,
    max_length=150
)
for record in filter_bam("alignments.bam", criteria):
    process(record)

# Pre-built Workflows
for record in variant_calling_filter("alignments.bam"):
    # High-quality (MAPQ≥30), primary, mapped alignments
    call_variants(record)
```

**Common Use Cases**:
- Variant calling preprocessing (MAPQ ≥ 30 filter)
- Strand-specific RNA-seq analysis
- Coverage analysis and CNV detection
- Quality control statistics
- Region extraction for targeted analysis

**Run**:
```bash
python bam_advanced_filtering.py path/to/alignments.bam
```

**Performance**:
- 4.54 million records/sec throughput
- 43.0 MiB/s compressed file processing
- 4× speedup via automatic parallel BGZF
- Constant ~5 MB memory footprint

---

## Rust Examples

### FASTQ Streaming

#### [basic_fastq_streaming.rs](basic_fastq_streaming.rs)
Demonstrates streaming FASTQ parsing:
- Iterator-based API (constant memory)
- Error handling patterns
- Record processing

**Run**:
```bash
cargo run --example basic_fastq_streaming
```

### BAM Parsing

#### [bam_streaming.rs](bam_streaming.rs)
BAM file streaming with parallel BGZF:
- Automatic compression detection
- Header and record parsing
- Filtering by MAPQ and flags
- Streaming iterator pattern

**Run**:
```bash
cargo run --example bam_streaming
```

#### [bam_to_sam.rs](bam_to_sam.rs)
Convert BAM to SAM format:
- Streaming BAM reader
- SAM text format output
- Constant memory conversion
- Production-ready converter

**Run**:
```bash
cargo run --example bam_to_sam -- input.bam output.sam
```

#### [bam_reader_test.rs](bam_reader_test.rs)
Basic BAM reader demonstration:
- Opening BAM files
- Accessing header information
- Reading first N records
- Simple filtering

**Run**:
```bash
cargo run --example bam_reader_test
```

### Network Streaming

#### [http_streaming.rs](http_streaming.rs)
HTTP range request streaming:
- Stream from URLs without downloading
- Smart caching with LRU
- Background prefetching
- Network configuration

**Run**:
```bash
cargo run --example http_streaming
```

#### [sra_streaming.rs](sra_streaming.rs)
SRA (Sequence Read Archive) streaming:
- SRA toolkit integration
- Accession-based streaming
- Network error handling
- Production pipeline integration

**Run**:
```bash
cargo run --example sra_streaming -- SRR390728
```

### K-mer Operations

#### [kmer_operations_full.rs](kmer_operations_full.rs)
Complete k-mer analysis workflows:
- K-mer extraction (DNABert preprocessing)
- Minimizer computation (minimap2-style)
- K-mer spectrum analysis
- Parallel extraction (2.2× speedup)

**Run**:
```bash
cargo run --example kmer_operations_full
```

### Sequence Operations

#### [sequence_operations.rs](sequence_operations.rs)
Sequence manipulation operations:
- Reverse complement
- Trimming (fixed position, quality-based)
- Masking (quality-based)
- Validation (DNA, RNA)

**Run**:
```bash
cargo run --example sequence_operations
```

### Compression & Tuning

#### [analyze_bgzip_blocks.rs](analyze_bgzip_blocks.rs)
Analyze BGZF compression structure:
- Block boundary detection
- Compression ratio analysis
- Performance profiling

**Run**:
```bash
cargo run --example analyze_bgzip_blocks -- file.fq.gz
```

#### [prefetch_tuning.rs](prefetch_tuning.rs)
Network prefetch configuration:
- Cache size tuning
- Prefetch distance optimization
- Performance benchmarking

**Run**:
```bash
cargo run --example prefetch_tuning
```

---

## Quick Start

### Python Examples

1. **Install biometal**:
   ```bash
   pip install biometal-rs
   ```

2. **Run an example**:
   ```bash
   python basic_usage.py
   ```

3. **Advanced filtering**:
   ```bash
   python bam_advanced_filtering.py your_alignments.bam
   ```

### Rust Examples

1. **Clone repository**:
   ```bash
   git clone https://github.com/shandley/biometal
   cd biometal
   ```

2. **Run an example**:
   ```bash
   cargo run --example basic_fastq_streaming
   ```

3. **List all examples**:
   ```bash
   cargo run --example
   ```

---

## Example Categories

### By Format
- **FASTQ**: `basic_fastq_streaming.rs`, `basic_usage.py`
- **BAM/SAM**: `bam_streaming.rs`, `bam_to_sam.rs`, `bam_advanced_filtering.py`
- **FASTA**: Covered in basic examples
- **Network**: `http_streaming.rs`, `sra_streaming.rs`

### By Operation
- **Streaming**: All examples use streaming architecture
- **Filtering**: `bam_advanced_filtering.py`, `bam_streaming.rs`
- **K-mers**: `kmer_operations_full.rs`
- **Sequence ops**: `sequence_operations.rs`
- **Conversion**: `bam_to_sam.rs`

### By Level
- **Beginner**: `basic_usage.py`, `basic_fastq_streaming.rs`
- **Intermediate**: `bam_streaming.rs`, `kmer_operations_full.rs`
- **Advanced**: `bam_advanced_filtering.py`, `sra_streaming.rs`

---

## Performance Notes

All examples demonstrate biometal's performance characteristics:

### ARM NEON Speedup
- **Base counting**: 16.7× faster on Apple Silicon
- **GC content**: 20.3× faster on Apple Silicon
- **Quality filter**: 25.1× faster on Apple Silicon

### Parallel BGZF
- **BAM parsing**: 4× faster (6.5× on decompression)
- **Throughput**: 4.54 million records/sec
- **Bandwidth**: 43.0 MiB/s compressed

### Memory Efficiency
- **Constant ~5 MB**: Regardless of file size
- **Stream 5TB files**: On 24GB consumer laptops
- **99.5% reduction**: vs traditional load-all approaches

---

## Common Use Cases

### Quality Control
```bash
# Python: Advanced BAM QC
python bam_advanced_filtering.py alignments.bam

# Rust: FASTQ streaming QC
cargo run --example sequence_operations
```

### Variant Calling Pipeline
```python
# Extract high-quality alignments
from bam_advanced_filtering import variant_calling_filter

for record in variant_calling_filter("alignments.bam"):
    # Process with variant caller
    call_variants(record)
```

### ML Preprocessing
```bash
# K-mer extraction for DNABert
cargo run --example kmer_operations_full
```

### Network Analysis
```bash
# Analyze public data without downloading
cargo run --example sra_streaming -- SRR390728
```

---

## Contributing Examples

Have a useful workflow? Contribute an example!

1. **Python examples**: Add to `examples/` with clear docstrings
2. **Rust examples**: Add to `examples/` with `//!` module docs
3. **Documentation**: Update this README with your example
4. **Test**: Ensure example runs successfully

See [CLAUDE.md](../CLAUDE.md) for development guidelines.

---

## Resources

- **Documentation**: https://docs.rs/biometal
- **Tutorials**: [../notebooks/README.md](../notebooks/README.md)
- **API Reference**: https://docs.rs/biometal
- **GitHub**: https://github.com/shandley/biometal
- **Issues**: https://github.com/shandley/biometal/issues

---

**biometal v1.2.0+** - ARM-native bioinformatics with streaming architecture
