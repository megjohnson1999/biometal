# biometal Python Tutorials - Jupyter Notebooks

**Learn biometal through real bioinformatics workflows**

Interactive Jupyter notebooks teaching biometal's Python API from basics to production pipelines.

---

## 📚 Learning Path

### Beginner Level (Start Here!)

#### [01_getting_started.ipynb](01_getting_started.ipynb)
**Duration**: 15-20 minutes | **Status**: ✅ Complete

Learn the basics of biometal:
- Installing and importing biometal
- Streaming FASTQ files with constant memory
- Calculating GC content and base composition
- Analyzing quality scores
- Understanding ARM NEON performance benefits

**Perfect for**: First-time users, understanding streaming architecture

---

### Intermediate Level (Coming Soon)

#### [02_quality_control_pipeline.ipynb](02_quality_control_pipeline.ipynb)
**Duration**: 30-40 minutes | **Status**: ✅ Complete

Build complete QC pipelines:
- Fixed position trimming (adapter removal)
- Quality-based trimming (Trimmomatic-style sliding window)
- Length filtering (meets_length_requirement)
- Quality-based masking (variant calling pipelines)
- Complete trim → filter → mask workflow
- Production-ready reusable function

**Perfect for**: Pre-processing reads for downstream analysis, **showcases v1.2.0 Phase 4 features**

#### [03_kmer_analysis.ipynb](03_kmer_analysis.ipynb)
**Duration**: 30-40 minutes | **Status**: ✅ Complete

K-mer extraction for machine learning:
- K-mer extraction for DNABert/DNABERT-2 (k=3-6)
- Minimizers for indexing (minimap2-style, w=10)
- K-mer spectrum analysis (frequency distribution)
- Parallel extraction (2.2× speedup, ≥1000 sequences)
- Complete ML preprocessing pipeline
- Evidence-based design (Entry 034)

**Perfect for**: ML practitioners, transformer model preprocessing, **showcases v1.1.0 k-mer features**

#### [04_network_streaming.ipynb](04_network_streaming.ipynb)
**Duration**: 30-40 minutes | **Status**: ✅ Complete

Network streaming without downloading:
- HTTP streaming with range requests
- Constant memory (~5 MB for any file size)
- Public data access (ENA, cloud storage)
- Network configuration and tuning
- SRA concepts and limitations
- Complete remote QC pipeline

**Perfect for**: Working with public datasets, cloud analysis, **democratizing access to TB-scale data**

#### [05_bam_alignment_analysis.ipynb](05_bam_alignment_analysis.ipynb)
**Duration**: 30-40 minutes | **Status**: ✅ Complete (Nov 8, 2025)

BAM alignment file analysis with parallel BGZF:
- Stream BAM files with constant memory (~5 MB)
- 4× speedup via automatic parallel BGZF decompression
- Filter by MAPQ, flags, and alignment properties
- Calculate alignment statistics (mapping rate, coverage)
- Production QC workflows
- 4.54 million records/sec throughput

**Perfect for**: Alignment QC, variant calling pipelines, **showcases production BAM parser with 4× speedup**

---

### Advanced Level (Future)

#### 05_complete_pipeline.ipynb
**Status**: 📋 Planned

End-to-end production workflow:
- Stream from SRA
- Complete QC pipeline
- K-mer extraction
- Complexity filtering
- Output to formats
- Performance metrics

**Perfect for**: Production bioinformatics pipelines

---

## 🚀 Quick Start

### Installation

```bash
# Install biometal
pip install biometal-rs

# Install Jupyter
pip install jupyter matplotlib seaborn pandas

# Clone repo (for notebooks)
git clone https://github.com/shandley/biometal
cd biometal/notebooks

# Start Jupyter
jupyter notebook
```

### Run in Browser

Open `01_getting_started.ipynb` and run all cells (Cell → Run All).

---

## 📊 What You'll Learn

### Core Concepts
✅ **Streaming Architecture**: Analyze 5TB datasets with 5 MB memory
✅ **ARM NEON Performance**: 16-25× speedup on Apple Silicon
✅ **Network Streaming**: Analyze without downloading
✅ **Production Quality**: Real workflows, not toy examples

### Practical Skills
✅ Quality control pipelines (trim, filter, mask)
✅ K-mer extraction for ML (DNABert preprocessing)
✅ Network streaming (HTTP, analyze without downloading)
✅ BAM alignment analysis (QC, filtering, statistics)
✅ Memory-efficient processing (constant ~5 MB)

### biometal Features (v1.2.0+)
✅ **Core Operations**: GC content, base counting, quality scores (v1.0.0)
✅ **K-mer Operations**: Extraction, minimizers, spectrum (v1.1.0)
✅ **Phase 4 Operations**: Trimming, masking, sequence manipulation (v1.2.0)
✅ **Network Streaming**: HTTP range requests, constant memory (v1.0.0)
✅ **BAM/SAM Parser**: Parallel BGZF, 4× speedup, constant memory (v1.2.0+, Nov 8 2025)

---

## 💻 Requirements

### Software
- **Python**: 3.9+ (tested on 3.9-3.14)
- **biometal-rs**: 1.2.0+ (`pip install biometal-rs`)
- **Jupyter**: `pip install jupyter`
- **Visualization** (optional): `pip install matplotlib seaborn pandas`

### Hardware
- **Minimum**: Any modern laptop (x86_64 or ARM)
- **Recommended**: Apple Silicon (M1/M2/M3/M4) for NEON acceleration
- **Memory**: 8 GB RAM (biometal uses ~5 MB, but Jupyter needs more)

### Data
- Test data included in notebooks (small synthetic files)
- Real data streamed from NCBI SRA (no download needed!)

---

## 🎯 Learning Outcomes

After completing these notebooks, you will be able to:

1. **Stream and analyze** large FASTQ/FASTA files with constant memory
2. **Build QC pipelines** for read preprocessing (trim → filter → mask)
3. **Extract k-mers** for machine learning (DNABert, transformers)
4. **Stream from HTTP** to analyze public data without downloading
5. **Analyze BAM alignments** with 4× speedup (filter, QC, statistics)
6. **Write production pipelines** using biometal's streaming API
7. **Optimize performance** on Apple Silicon (ARM NEON + parallel BGZF)

---

## 📖 Tutorial Format

Each notebook follows this structure:

1. **Learning Objectives** - What you'll learn
2. **Prerequisites** - What you need to know
3. **Concepts** - Background and motivation
4. **Hands-On Code** - Real examples with explanations
5. **Visualizations** - Charts and plots showing results
6. **Best Practices** - Production tips
7. **Exercises** - Try it yourself
8. **What's Next** - Link to next tutorial

---

## 🔬 Real Data Examples

### Included Test Data
- Synthetic FASTQ files (generated in notebooks)
- Small E. coli samples (~1000 reads)
- Bacterial sequences for k-mer analysis

### Network Streaming Examples
- **Public FASTQ data**: ENA, cloud storage (S3, GCS, Azure)
  - Stream via HTTP (no download required!)
  - Constant memory (~5 MB)
  - Demonstrated in notebook 04
- **Future**: Direct SRA streaming (requires SRA Toolkit)

---

## 🆘 Getting Help

### Common Issues

**Q: Import error `ModuleNotFoundError: No module named 'biometal'`**
A: Install with `pip install biometal-rs` (note the `-rs`)

**Q: Notebook says "Kernel starting, please wait"**
A: Give it 10-30 seconds on first start

**Q: Can't find test data files**
A: Test data is generated in the notebooks (run all cells)

**Q: Performance seems slow**
A: Check your platform (run platform check cell). NEON acceleration only on ARM.

### Support

- **GitHub Issues**: https://github.com/shandley/biometal/issues
- **Discussions**: https://github.com/shandley/biometal/discussions
- **Documentation**: https://docs.rs/biometal

---

## 🤝 Contributing

Found an issue or want to improve a notebook?

1. Open an issue: https://github.com/shandley/biometal/issues
2. Suggest improvements
3. Share your own examples!

---

## 📜 License

These tutorials are part of biometal and licensed under MIT OR Apache-2.0 (your choice).

---

## 🙏 Acknowledgments

- **Evidence Base**: 1,357 experiments from apple-silicon-bio-bench
- **Community**: Feedback from bioinformaticians worldwide
- **Platform**: Optimized for Apple Silicon (M1/M2/M3/M4)

---

**biometal v1.2.0** - ARM-native bioinformatics with streaming architecture

**Mission**: Democratizing bioinformatics compute for LMIC researchers, small labs, students, and field researchers.
