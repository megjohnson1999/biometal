# biometal: Claude Development Guide

**Project**: biometal - ARM-native bioinformatics library
**Latest Release**: v1.7.0 (November 13, 2025)
**Current Focus**: Core development, compression optimized
**Research Status**: CAF columnar format evaluated (Nov 4-11) - Archived

---

## Mission

Democratize bioinformatics by enabling 5TB dataset analysis on consumer hardware through:
- **Streaming architecture**: Constant ~5 MB memory (not load-all)
- **ARM-native performance**: 16-25× NEON speedup
- **Network streaming**: Analyze without downloading
- **Evidence-based optimization**: Every rule validated experimentally

**Target**: LMIC researchers, small labs, students, field researchers, ML practitioners

---

## Core Principles

### 1. Evidence-Based Design

Every optimization comes from validated experimental results (apple-silicon-bio-bench):
- Follow OPTIMIZATION_RULES.md: 6 rules from 1,357 experiments (N=30)
- Don't guess: Reference ASBB evidence for optimization decisions
- Document rationale: Link implementations to specific rules/entries

### 2. Streaming-First Architecture

Always design for constant memory:
- Bad: `Vec<FastqRecord>` (accumulates in memory)
- Good: `Iterator<Item = FastqRecord>` (constant memory)
- Target: ~5 MB regardless of dataset size (Rule 5)

### 3. ARM-Native with Portable Fallback

Always provide both ARM and fallback implementations:
```rust
#[cfg(target_arch = "aarch64")]
pub fn operation_neon(input: &[u8]) -> Result { /* 16-25× faster */ }

#[cfg(not(target_arch = "aarch64"))]
pub fn operation_scalar(input: &[u8]) -> Result { /* x86_64 fallback */ }

pub fn operation(input: &[u8]) -> Result {
    #[cfg(target_arch = "aarch64")]
    { operation_neon(input) }
    #[cfg(not(target_arch = "aarch64"))]
    { operation_scalar(input) }
}
```

Platform priority: Mac ARM → Linux ARM (Graviton) → x86_64 fallback

### 4. Production Quality

- Use `Result<T, BiometalError>` (no panics in library)
- Document every public API with examples
- Property-based testing (proptest)
- Benchmarks (criterion, N=30)
- No `unwrap()` or `expect()` in library code

---

## Recent Research: CAF Columnar Format (Archived)

**Period**: November 4-11, 2025
**Location**: research/caf-format/implementation/
**Status**: ✅ Research Complete - Archived

### What Was Built
- Production-quality columnar alignment format (~6,000 lines)
- BAM ↔ CAF conversion (lossless)
- Dictionary compression (86% quality score reduction)
- ARM NEON optimizations (16× base counting speedup)
- Comprehensive documentation and benchmarks (N=30)

### Key Findings
- **File size**: 1.6× LARGER than BAM (vs target 0.5-1.0×) ❌
- **Column-selective reading**: 1.4× speedup (vs predicted 3-5×) ⚠️
- **NEON effectiveness**: Varies by operation (1.3-16× speedup)
- **Compiler auto-vec**: Can beat explicit SIMD for simple operations
- **Research value**: Good methodology, valuable negative results ✅

### Verdict
CAF doesn't outperform BAM/CRAM enough to justify adoption (files 60% larger, modest speedups). Valuable for methodology demonstration and lessons learned. See `research/caf-format/implementation/CAF_FINAL_REPORT.md` for full analysis.

### Lessons Applied to biometal
1. Profile before assuming SIMD is faster
2. Data layout significantly affects performance
3. Benchmark with N=30 for statistical rigor
4. Document negative results honestly

**Status**: Archived for reference, returning to biometal core development

---

## Recent Optimization: Compression Backend (v1.7.0)

**Period**: November 12-13, 2025
**Location**: BACKEND_COMPARISON_FINDINGS.md, COMPRESSION_INVESTIGATION_FINDINGS.md
**Status**: ✅ COMPLETE - Production deployed

### What Was Achieved
- Comprehensive backend comparison (rust_backend, zlib-ng, cloudflare_zlib)
- cloudflare_zlib delivers best performance (1.67× decompression, 2.29× compression)
- Public compression API with configurable levels (fast/default/best)
- All 411 tests passing, production-ready

### Performance Impact
- **BAM parsing**: 55 MiB/s → 92 MiB/s (+67% improvement)
- **Decompression**: 1.67× faster vs rust_backend (miniz_oxide)
- **Compression (default)**: 64 MB/s (2.29× vs rust_backend)
- **Compression (fast)**: 358 MB/s (5.6× vs default)
- **Compression is now faster than decompression** (358 MB/s vs 290 MB/s in fast mode)

### Key Findings
- cloudflare_zlib 3-7% faster than zlib-ng (consistent across file sizes)
- Fast compression mode only 3-5% larger files (minimal quality penalty)
- Best compression mode (1.2× slower than default) provides 5-10% smaller files
- Performance scales linearly across file sizes (5MB → 544MB)

**Status**: Full details in BACKEND_COMPARISON_FINDINGS.md and COMPRESSION_INVESTIGATION_FINDINGS.md

---

## Current Status (v1.7.0)

### Released Features
- **FASTQ/FASTA** streaming parsers (constant memory)
- **ARM NEON** operations (base counting, GC content, quality filtering)
- **Sequence manipulation** (reverse_complement, trimming, masking)
- **K-mer operations** (extraction, minimizers, spectrum)
- **Network streaming** (HTTP, SRA)
- **BAM/SAM parser** (v1.4.0, November 8, 2025 - production-ready)
  - Full BAM parsing (header, records, CIGAR, tags, sequences)
  - ARM NEON sequence decoding (+27.5% BAM parsing speedup)
  - SAM writing for downstream tools
  - Robustness features (oversized CIGAR, malformed record handling)
  - 70 tests passing (integration complete)
- **Python bindings** (PyO3 0.27, 40+ functions)
  - Full BAM/FASTQ/FASTA support
  - CIGAR operations, SAM writing
- **Tests**: 347 passing (260 library + 87 doc)

### Optimization Rules Implemented

| Rule | Feature | Status | Impact |
|------|---------|--------|--------|
| **Rule 1** | ARM NEON SIMD | ✅ v1.0.0 | 16-25× speedup |
| **Rule 2** | Block-based processing | ✅ v1.0.0 | Preserves NEON gains |
| **Rule 3** | Parallel BGZF | ⏳ Phase 2 | 6.5× (planned) |
| **Rule 4** | Smart mmap | ⏳ Phase 2 | 2.5× (planned) |
| **Rule 5** | Constant-memory streaming | ✅ v1.0.0 | 99.5% memory reduction |
| **Rule 6** | Network streaming | ✅ v1.0.0 | Enables remote analysis |

**Current**: 4/6 rules (67%)
**Phase 2 Target**: 6/6 rules (100%), 27× combined speedup

### Distribution
- **PyPI**: biometal-rs v1.4.0 (pip install biometal-rs)
- **crates.io**: biometal v1.4.0 (cargo add biometal)

---

## Current Work: Phase 1 Consolidation

### Weeks 1-2: ✅ COMPLETE (November 10, 2025)

**Documentation Sprint (Week 1)**:
- ✅ User Guide (25,000+ words): docs/USER_GUIDE.md
- ✅ Performance Guide (10,000+ words): docs/PERFORMANCE_OPTIMIZATION_GUIDE.md
- ✅ BAI Tutorial (Jupyter): notebooks/07_bai_indexed_queries.ipynb
- ✅ Enhanced API docs: src/io/bam/index.rs

**Performance Benchmarking (Week 2)**:
- ✅ Comprehensive comparison vs samtools/pysam: benchmarks/comparison/BENCHMARK_COMPARISON.md
- ✅ Validated all v1.6.0 claims (1.68× indexed speedup, 10-200× memory advantage)
- ✅ Real-world scenario analysis (3 production use cases)
- ✅ Automated benchmark framework: benchmarks/comparison/samtools_vs_biometal.sh

**Status**: 50% of Phase 1 complete, strong foundation for community launch

### Weeks 3-4: 🔄 IN PROGRESS

**Community Building (Week 3)**:
- [ ] Blog post announcing v1.6.0
- [ ] Social media campaign (Twitter, Reddit, Biostars, LinkedIn)
- [ ] Engage with tool maintainers (samtools, pysam, HTSlib)
- [ ] Set up GitHub discussions and issue templates

**Quality Assurance (Week 4)**:
- [ ] Property-based testing expansion
- [ ] Fuzz testing for robustness
- [ ] Cross-platform validation (Graviton, x86_64)
- [ ] Memory safety audit (Valgrind, ASAN, Miri)

---

## Future Work

### Phase 2: High-ROI Performance (Weeks 5-9)

**Objective**: Implement remaining optimization rules for 27× combined speedup

**Rule 3: Parallel BGZF Decompression** (6.5× speedup)
- Current: 55 MiB/s
- Target: 358 MiB/s
- Effort: 40-60 hours
- Evidence: Entry 029 (CPU parallel prototype)

**Rule 4: Smart mmap** (2.5× additional)
- Combined with Rule 3: 895 MiB/s
- Total: **16× improvement** over current
- Effort: 40-60 hours
- Evidence: Entry 032 (scale validation)

**Expected Outcome**:
- Sequential BAM parsing: 55 MiB/s → **895 MiB/s**
- All 6 optimization rules implemented (100%)
- Validated against ASBB experiments

### Phase 3: Strategic Expansion (Weeks 10-14)

**Deferred pending Phase 1 community feedback**:
- Format expansion (CRAM, VCF, CSI index)
- Horizontal expansion (GFF/GTF, BED)
- Community-driven features

---

## Key Documentation

### For Users
- **📘 User Guide**: docs/USER_GUIDE.md - Comprehensive onboarding (installation → optimization)
- **📓 Tutorials**: notebooks/ - 7 Jupyter notebooks (including BAI indexed queries)
- **⚡ Performance Guide**: docs/PERFORMANCE_OPTIMIZATION_GUIDE.md - Maximize performance
- **📊 Benchmarks**: benchmarks/comparison/BENCHMARK_COMPARISON.md - vs samtools/pysam

### For Developers
- **📐 Architecture**: docs/ARCHITECTURE.md - Technical design
- **🔬 Optimization Rules**: OPTIMIZATION_RULES.md - Evidence-based optimization (6 rules)
- **🐍 Python Bindings**: docs/PYTHON.md - Python-specific details
- **🧬 BAM API**: docs/BAM_API.md - Complete BAM/SAM parser reference

### For Planning
- **📈 Phase 1 Progress**: PHASE1_PROGRESS_REPORT.md - Current consolidation status
- **🗺️ Strategic Analysis**: NEXT_STEPS_ANALYSIS.md - Long-term roadmap (3-phase, 14 weeks)
- **📝 Changelog**: CHANGELOG.md - Version history

---

## Development Workflow

### Session Guidelines

**What to Emphasize**:
- Evidence-based design (follow OPTIMIZATION_RULES.md)
- Streaming-first architecture (constant memory)
- ARM-native with portable fallback
- Production quality (error handling, docs, tests)

**What NOT to Do**:
- Don't make up optimization parameters (refer to evidence)
- Don't accumulate records in memory (streaming only)
- Don't panic in library code (use Result)
- Don't implement ARM-only code without scalar fallback

**Decision Framework**:
When implementing features:
1. Check OPTIMIZATION_RULES.md for relevant rule
2. Follow the implementation pattern for that rule
3. Link to evidence (lab notebook entry)

When evaluating optimizations:
1. Is this validated in ASBB experiments?
2. If yes: Which rule/entry documents it?
3. If no: Suggest validating first or using proven approach

### Testing Strategy

**Property-Based Testing**:
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

**Benchmarking**:
```rust
use criterion::{criterion_group, criterion_main, Criterion};

fn bench_operation(c: &mut Criterion) {
    let data = generate_test_data(100_000);
    c.bench_function("operation_neon", |b| {
        b.iter(|| operation_neon(&data))
    });
}

criterion_group!(benches, bench_operation);
criterion_main!(benches);
```

---

## Performance Expectations

### Current Performance (v1.7.0)

| Operation | Scalar | Optimized | Speedup |
|-----------|--------|-----------|---------|
| Base counting | 315 Kseq/s | 5,254 Kseq/s | **16.7× (NEON)** |
| GC content | 294 Kseq/s | 5,954 Kseq/s | **20.3× (NEON)** |
| Quality filter | 245 Kseq/s | 6,143 Kseq/s | **25.1× (NEON)** |
| BAM parsing | ~11 MiB/s | 92.0 MiB/s | **8.4× (BGZF + NEON + cloudflare_zlib)** |
| BAM indexed query | O(n) full scan | O(log n) indexed | **1.68-500× (scales with file size)** |

### Phase 2 Target Performance

| Operation | Current | Phase 2 Target | Improvement |
|-----------|---------|----------------|-------------|
| BAM parsing | 92 MiB/s | **895 MiB/s** | **9.7× (Rules 3+4)** |
| Memory usage | **5 MB** | **5 MB** | Constant (maintained) |
| Optimization rules | 4/6 (67%) | **6/6 (100%)** | Complete |

---

## Competitive Position (Validated Week 2)

### vs samtools

| Metric | biometal | samtools | Advantage |
|--------|----------|----------|-----------|
| Sequential parsing | 55.1 MiB/s | ~45-50 MiB/s | ✅ Competitive |
| Indexed queries | 1.68-500× | ~1.2-1.5× | ✅ **Superior** |
| Memory | **5 MB** | 20-50 MB | ✅ **10× lower** |
| ARM NEON | **4-25× speedup** | None | ✅ **Exclusive** |

### vs pysam

| Metric | biometal | pysam | Advantage |
|--------|----------|-------|-----------|
| Python performance | ~45 MiB/s | ~30-40 MiB/s | ✅ 1.5-2× faster |
| Memory | **5 MB** | 50 MB-1 GB | ✅ **10-200× lower** |
| API | Streaming | Context managers | ✅ Simpler |
| ARM NEON | **4-25× speedup** | None | ✅ **Exclusive** |

**Production Use Cases**:
- ✅ Large-file targeted analysis (indexed queries)
- ✅ Memory-constrained environments
- ✅ ARM infrastructure (Apple Silicon, Graviton)
- ✅ Terabyte-scale streaming

---

## Quick Reference

### Evidence Base
- **1,357 experiments**, 40,710 measurements (N=30)
- Source: apple-silicon-bio-bench
- Rules: 6 optimization rules (OPTIMIZATION_RULES.md)

### Platform Support
1. **Mac ARM** (M1/M2/M3/M4): 16-25× NEON speedup (optimized)
2. **Linux ARM** (Graviton): 6-10× NEON speedup (portable)
3. **x86_64**: 1× scalar fallback (portable)

### File Formats
- ✅ FASTQ, FASTA (v1.0.0)
- ✅ BAM, SAM (v1.4.0)
- ✅ BAI index (v1.6.0)
- ⏳ CSI index (partial)
- ❌ CRAM, VCF, BCF (future)

### Tests
- **582 passing** (100% pass rate)
  - 354 library tests
  - 81 BAM tests
  - 26 BAI Python tests
  - 121 documentation tests

---

## Session Checklist

When starting a new session:
- [ ] Review current phase status (PHASE1_PROGRESS_REPORT.md)
- [ ] Check CHANGELOG.md for recent changes
- [ ] Review relevant optimization rules (OPTIMIZATION_RULES.md)
- [ ] Check open issues/PRs if community-facing work

When implementing features:
- [ ] Follow evidence-based design (link to ASBB entry)
- [ ] Use streaming architecture (constant memory)
- [ ] Provide ARM + fallback implementations
- [ ] Add property-based tests
- [ ] Benchmark with criterion (N=30)
- [ ] Document with examples

When wrapping up:
- [ ] Update CHANGELOG.md
- [ ] Run full test suite
- [ ] Update relevant planning documents
- [ ] Document any decisions made

---

**Last Updated**: November 10, 2025 (Post v1.6.0 release, Phase 1 Weeks 1-2 complete)
**Next Milestone**: Phase 1 Weeks 3-4 (Community Building + Quality Assurance)
**Long-term Goal**: Phase 2 (Rules 3+4, 16× performance improvement)
