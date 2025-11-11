# CAF Research Plan: Detailed Methodology

**Project**: Columnar Alignment Format (CAF)
**Duration**: 8 weeks (Nov 10, 2025 - Jan 10, 2026)
**Type**: Experimental research with publication target

---

## Research Objectives

### Primary Objective

**Design, implement, and validate a columnar alignment format (CAF) optimized for ARM NEON that achieves 5-10× performance improvement over BAM for analytical bioinformatics operations.**

### Secondary Objectives

1. Demonstrate lossless BAM ↔ CAF conversion
2. Characterize storage trade-offs (1.5-2× larger files)
3. Validate across multiple platforms (ARM, x86_64)
4. Publish open-source implementation
5. Submit manuscript to peer-reviewed journal

---

## Research Questions

### RQ1: Performance
**Can CAF achieve 5-10× speedup over BAM for analytical operations using ARM NEON?**

**Hypotheses**:
- H1.1: Quality filtering ≥20× faster (CAF vs BAM)
- H1.2: Base counting ≥20× faster (CAF vs BAM)
- H1.3: MAPQ filtering ≥16× faster (CAF vs BAM)
- H1.4: Overall workflow ≥5× faster (CAF vs BAM)

**Statistical test**: Paired t-test (N=30, p < 0.05)

### RQ2: Correctness
**Does CAF maintain 100% lossless conversion with BAM?**

**Validation**:
- Differential testing (1,000+ diverse BAM files)
- Property-based testing (proptest)
- Round-trip validation (BAM → CAF → BAM)

**Success criterion**: 100% bit-identical reconstruction

### RQ3: Storage Trade-offs
**What is the storage overhead of CAF compared to BAM?**

**Measurements**:
- File size comparison (N=30 diverse datasets)
- Compression ratio analysis
- Storage vs performance trade-off curves

**Acceptable range**: 1.5-2.0× larger than BAM

### RQ4: Generalizability
**Does CAF performance generalize across diverse datasets and platforms?**

**Validation**:
- Dataset diversity: WGS, exome, RNA-seq, targeted
- Platform diversity: M1/M2/M3/M4, Graviton, x86_64
- Workload diversity: Filter, aggregate, transform

---

## Methodology

### Phase 0: Preparation (Week 1, Nov 10-17)

#### Literature Review

**Scope**:
1. **Alignment formats**: SAM/BAM (2009), CRAM (2015), alternatives
2. **Columnar formats**: Parquet, Arrow, ORC
3. **SIMD optimization**: NEON, AVX2, prior art
4. **Compression**: zstd, lz4, benchmarks

**Deliverable**: `LITERATURE_REVIEW.md` with 20+ citations

#### Specification Finalization

**Tasks**:
- Finalize columnar block structure
- Define compression strategy
- Specify index format
- Document API design

**Deliverable**: `SPECIFICATION.md` v1.0

#### Infrastructure Setup

**Tasks**:
- Create research project structure
- Set up Rust workspace
- Configure benchmarking (Criterion, N=30)
- Prepare test datasets

**Deliverable**: Working build system

---

### Phase 1: Implementation (Weeks 2-3, Nov 18 - Dec 1)

#### Core Data Structures

**Week 2, Days 1-2**:
```rust
pub struct CafBlock {
    // Columnar arrays
    positions: Vec<i32>,         // zstd compressed
    mapq: Vec<u8>,               // raw/RLE
    flags: Vec<u16>,             // bitpacked
    sequences: Vec<u8>,          // ASCII, lz4
    seq_offsets: Vec<u32>,       // record boundaries
    qualities: Vec<u8>,          // raw
    qual_offsets: Vec<u32>,      // shared with seq
    cigar_ops: Vec<u32>,         // RLE
    cigar_offsets: Vec<u32>,     // per-record
    read_names: Vec<String>,     // dictionary
    tags: TagStorage,            // nested columnar
}
```

**Validation**: Unit tests for each field

#### BAM → CAF Converter

**Week 2, Days 3-5**:
1. Read BAM using biometal's BamReader
2. Accumulate 10,000 records
3. Convert to columnar layout
4. Compress columns (zstd/lz4)
5. Write CAF block

**Validation**: Differential testing (10 diverse BAM files)

#### CAF → BAM Converter

**Week 2, Days 6-7**:
1. Read CAF block
2. Decompress columns
3. Convert to row-oriented BAM records
4. Write using BAM writer

**Validation**: Round-trip test (BAM → CAF → BAM, bit-identical)

#### Comprehensive Testing

**Week 3**:
- Unit tests (100+ tests)
- Integration tests (real BAM files)
- Edge cases (empty files, large CIGARs, etc.)
- Property-based tests (proptest)

**Success criterion**: 100% test pass rate

---

### Phase 2: NEON Optimization (Weeks 4-5, Dec 2-15)

#### NEON Kernel Implementation

**Week 4, Days 1-3**: Quality Filtering
```rust
#[cfg(target_arch = "aarch64")]
unsafe fn filter_quality_neon(block: &CafBlock, min_q: u8) -> Vec<usize> {
    let threshold = vdupq_n_u8(min_q);
    let mut passing = Vec::new();

    // Process 16 quality scores at once
    for (i, chunk) in block.qualities.chunks(16).enumerate() {
        let quals = vld1q_u8(chunk.as_ptr());
        let mask = vcgeq_u8(quals, threshold);

        if vaddvq_u8(mask) > 0 {
            // At least one passes
            passing.push(i);
        }
    }

    passing
}
```

**Target**: ≥20× speedup vs scalar

**Week 4, Days 4-5**: Base Counting
- Reuse proven biometal NEON kernel
- Adapt for CAF's pre-decoded ASCII sequences
- Target: ≥20× speedup vs scalar

**Week 5, Days 1-2**: MAPQ Filtering
- Parallel comparison across mapq array
- Target: ≥16× speedup vs scalar

**Week 5, Days 3-5**: Benchmarking
- Criterion benchmarks (N=30 samples)
- CAF vs BAM comparison
- Statistical analysis (t-test, p < 0.05)

**Deliverable**: `benchmarks/BENCHMARK_RESULTS.md` with statistical validation

---

### Phase 3: Validation & Analysis (Week 6, Dec 16-22)

#### Correctness Validation

**Tasks**:
1. Differential testing (1,000+ BAM files)
   - Diverse sources: 1000 Genomes, GIAB, synthetic
   - Edge cases: Empty, large, unusual tags

2. Property-based testing
   ```rust
   proptest! {
       #[test]
       fn roundtrip_preserves_all_fields(bam in arb_bam()) {
           let caf = convert_to_caf(&bam)?;
           let reconstructed = convert_to_bam(&caf)?;
           prop_assert_eq!(bam, reconstructed);
       }
   }
   ```

3. Cross-platform validation
   - Mac ARM (M1 Max)
   - AWS Graviton (Linux ARM)
   - GitHub Actions (x86_64)

**Success criterion**: 100% lossless, all platforms

#### Performance Validation

**Benchmark Protocol** (N=30):
1. Select 30 diverse BAM files (100 KB - 100 MB)
2. For each file, run 30 iterations:
   - BAM quality filter (baseline)
   - CAF quality filter (NEON)
   - Measure time, memory
3. Statistical analysis:
   - Mean speedup
   - 95% confidence intervals
   - Paired t-test (p < 0.05)
4. Repeat for each operation (base count, MAPQ filter)

**Deliverable**: Publication-quality figures and tables

#### Storage Analysis

**Measurements**:
- File size: CAF vs BAM (N=30)
- Compression ratios by column
- Storage-performance trade-off curves

**Deliverable**: Trade-off analysis document

---

### Phase 4: Publication (Weeks 7-8, Dec 23 - Jan 10)

#### Manuscript Drafting

**Week 7**:
- Abstract (200 words)
- Introduction (300 words)
- Methods (400 words)
- Results (400 words)
- Discussion (300 words)
- Figures (3-4)
- Tables (2-3)

**Week 8**:
- Revisions
- Supplementary materials
- Code/data availability statements
- Submission preparation

#### Figures & Tables

**Figure 1**: CAF block structure (schematic)
**Figure 2**: Performance comparison (bar chart, N=30, error bars)
**Figure 3**: Storage trade-offs (scatter plot)
**Figure 4**: NEON optimization results (heatmap)

**Table 1**: Benchmark results (mean ± SD, p-values)
**Table 2**: Storage comparison
**Table 3**: Platform validation

#### Submission

**Target**: Bioinformatics (Oxford) - Application Note
**Deadline**: January 10, 2026
**Format**: 2,000 words, 3 figures, 2 tables

---

## Statistical Methodology

### Sample Size Justification

**N=30 per benchmark**:
- Criterion default for statistical power
- Sufficient for t-test (power > 0.8, α = 0.05, effect size > 1.0)
- Matches biometal's apple-silicon-bio-bench protocol

### Statistical Tests

**Primary analysis**: Paired t-test
- Null hypothesis: No difference between CAF and BAM
- Alternative: CAF faster than BAM (one-tailed)
- Significance level: α = 0.05
- Effect size: Cohen's d

**Secondary analysis**:
- Descriptive statistics (mean, median, SD)
- 95% confidence intervals
- Distribution plots (histograms, box plots)

### Multiple Testing Correction

**Bonferroni correction**: α / n comparisons
- 4 operations tested → α = 0.05 / 4 = 0.0125

---

## Data Management

### Test Datasets

**Public sources**:
1. **1000 Genomes** (WGS, high coverage)
2. **GIAB** (high-confidence variants)
3. **SRA** (diverse species, protocols)
4. **Synthetic** (edge cases, controlled)

**Diversity criteria**:
- File size: 1 KB - 1 GB
- Coverage: 1× - 100×
- Read length: 50 bp - 150 bp
- Platforms: Illumina, PacBio, ONT

### Data Sharing

**Code**: GitHub (MIT license)
- Full implementation
- Test suite
- Benchmarking scripts
- Analysis code

**Data**: Zenodo (DOI)
- Benchmark results (CSV)
- Test datasets (if not public)
- Analysis outputs

**Reproducibility**:
- Docker container
- Conda environment
- Random seeds documented

---

## Risk Management

### Technical Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| NEON speedup <5× | Low | High | BAM validation shows ≥2× floor, columnar should improve |
| Lossless conversion fails | Low | Critical | Comprehensive differential testing |
| Storage overhead >2× | Medium | Medium | Acceptable trade-off, clearly documented |
| Platform issues (x86_64) | Medium | Low | Scalar fallback, still useful on ARM |

### Schedule Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| Implementation delays | Medium | Medium | Phased approach, can cut Phase 3 (indexing) |
| Benchmark analysis slow | Low | Low | Automate with scripts |
| Publication delays | Medium | Low | Pre-print on bioRxiv while in review |

### Scientific Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| Results not significant | Low | High | Pilot data suggests >10× speedup |
| Reviewers reject approach | Medium | Medium | Thorough validation, open data |
| Low community interest | High | Low | Research contribution regardless |

---

## Quality Assurance

### Code Review
- Self-review with checklist
- Automated testing (CI)
- Code coverage ≥95%
- Clippy warnings = 0

### Data Validation
- All benchmarks automated
- Random seeds fixed
- Results reproducible
- Version control for all analysis

### Documentation
- Inline code comments
- API documentation (rustdoc)
- User guide
- Research diary (daily updates)

---

## Milestones & Deliverables

### Phase 0 (Week 1)
- [ ] Literature review complete (20+ citations)
- [ ] Specification finalized (v1.0)
- [ ] Research infrastructure ready

### Phase 1 (Weeks 2-3)
- [ ] CAF reader/writer implemented
- [ ] Lossless conversion validated
- [ ] 100+ tests passing

### Phase 2 (Weeks 4-5)
- [ ] NEON kernels implemented
- [ ] Benchmarks complete (N=30)
- [ ] 5-10× speedup demonstrated

### Phase 3 (Week 6)
- [ ] Correctness validated (1,000+ files)
- [ ] Statistical analysis complete
- [ ] Figures and tables ready

### Phase 4 (Weeks 7-8)
- [ ] Manuscript drafted
- [ ] Submitted to journal
- [ ] Code and data published

---

## Success Metrics

### Quantitative
- ✅ Performance: ≥5× speedup (p < 0.05)
- ✅ Correctness: 100% lossless
- ✅ Storage: 1.5-2.0× overhead
- ✅ Tests: ≥95% coverage, 100% pass rate

### Qualitative
- ✅ Publication accepted
- ✅ Code open-sourced
- ✅ Community interest (stars, citations)
- ✅ Implementation used in pipelines

---

## Timeline Summary

```
Week 1:  Preparation
Week 2:  Implementation (core)
Week 3:  Implementation (testing)
Week 4:  NEON optimization
Week 5:  Benchmarking
Week 6:  Validation & analysis
Week 7:  Manuscript drafting
Week 8:  Revision & submission
```

**Target submission**: January 10, 2026
**Estimated publication**: Q2 2026 (if accepted)

---

**Document Status**: Draft v1.0
**Last Updated**: November 10, 2025
**Next Review**: Weekly (every Monday)
