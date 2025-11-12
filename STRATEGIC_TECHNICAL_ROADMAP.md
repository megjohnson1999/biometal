# biometal: Strategic Technical Development Analysis

**Date**: November 11, 2025
**Version**: v1.6.0 (Revised Post-Rules 3+4 Investigation)
**Status**: Strategic planning with evidence-based corrections

---

## Executive Summary

After comprehensive analysis of the codebase, documentation, and evidence base (1,357 experiments from apple-silicon-bio-bench), plus systematic investigation of Rules 2, 3, and 4 (November 11, 2025), this document provides **evidence-corrected** strategic priorities.

**CRITICAL UPDATES** (November 11, 2025):
- ✅ **Rule 2**: Investigated - Convenience API only, not 14× speedup (requires major refactor)
- ❌ **Rule 3**: Multi-scale tested - 0.77-0.84× slowdown (NOT 6.5× speedup!) → PRUNED
- ⚠️ **Rule 4**: Validated - ~1% benefit (NOT 2.5×) due to CPU-bound decompression
- 🔍 **Bottleneck identified**: Decompression (98.7% of time), not I/O (1.3%)

**Strategic Pivot**: Original Phase 2A target (16.3× speedup via Rules 3+4) is **not achievable** with current architecture. Actual improvement: ~1.3× (removed parallel penalty + minimal mmap benefit).

**New Priority**: Target the actual bottleneck (decompression) or accept current performance and focus on horizontal expansion.

---

## INVESTIGATION RESULTS: Rules 2, 3, and 4 (November 11, 2025)

### Rule 2: Block Processing - Convenience API Only

**Status**: Implemented as convenience feature, does NOT achieve 14× speedup
**Benchmarks**: N=30 samples show block performance identical to per-record (0.88-1.01×)

**Root Cause Analysis**:
Entry 027's 82-86% overhead comes from calling operation functions 10,000 times (once per sequence). The "batch" mode in ASBB calls the function ONCE with all 10K sequences, keeping NEON registers hot.

biometal's current implementation:
- Per-record API: Calls `count_bases()` 10,000 times ✗
- Block API: Also calls `count_bases_neon()` 10,000 times internally ✗

**Both have the same overhead → No speedup achieved**

**Benchmark Results (November 11, 2025, N=30)**:

| Operation | Per-Record | Block | Speedup |
|-----------|-----------|-------|---------|
| Base counting | 111.07 µs | 126.49 µs | 0.88× (slower!) |
| GC content | 80.24 µs | 81.93 µs | 0.98× (same) |
| Mean quality | 41.20 µs | 40.84 µs | 1.01× (same) |
| QC workflow | 235.48 µs | 254.10 µs | 0.93× (slower!) |

**What Would Be Required for 14× Speedup**:
- Inline NEON operations directly into block functions (no function calls)
- Code duplication: ~1,200 lines (3 operations × 2 variants × ~200 lines each)
- Maintenance burden: Keep two implementations in sync
- **Effort**: 40-60 hours

**See**: `RULE2_INVESTIGATION_FINDINGS.md`

---

### Rule 3: Parallel BGZF - FAILED (Architecture Conflict)

**Status**: Implemented, multi-scale tested, **DISABLED** after failing DAG pruning
**Benchmarks**: N=10 samples across 3 file sizes show 0.77-0.84× slowdown (NOT 6.5× speedup!)

**Multi-Scale Results**:

| File Size | Blocks | Sequential | Parallel Bounded | Speedup | Verdict |
|-----------|--------|------------|------------------|---------|------------|
| **5.4 MB** | ~474 | 44.0 ms | 52.1 ms | **0.84×** | ❌ Slower |
| **54 MB** | ~4,747 | 437.3 ms | 533.2 ms | **0.82×** | ❌ Slower |
| **544 MB** | ~47,497 | 4.39 s | 5.67 s | **0.77×** | ❌ Slower |

**Critical Finding**: Performance DEGRADES with scale (opposite of Entry 029's pattern)

**Root Cause**:
- **Entry 029**: All-at-once decompression (load entire file → decompress all blocks in parallel) → 6.5× ✓
- **biometal**: Bounded streaming (8 blocks at a time for constant memory) → 0.77× ✗
- **Architectural conflict**: Rule 3 (parallelism) incompatible with Rule 5 (constant memory streaming)

**DAG Decision**: Failed pruning threshold (<1.5×) → **Optimization removed**

**Trade-off**: Prioritize Rule 5 (streaming for TB-scale files) over Rule 3 (speed)

**See**: `RULE3_AND_RULE4_SESSION_SUMMARY.md`, `RULE3_BENCHMARK_RESULTS.md`

---

### Rule 4: Smart mmap - LIMITED BENEFIT (~1%, not 2.5×)

**Status**: Implemented and working, but provides minimal benefit for compressed files
**Benchmarks**: N=10 samples reveal bottleneck is decompression, not I/O

**Bottleneck Analysis (544 MB file)**:

| Operation | Time | % of Total |
|-----------|------|------------|
| **Reading compressed file** (I/O) | 55.1 ms | 1.3% |
| **Decompressing** (CPU) | 4.37 s | 98.7% |
| **Total** | 4.42 s | 100% |

**Decompression is 79× slower than I/O** (4.37s / 55ms = 79×)

**Amdahl's Law Application**:
- If mmap gives 2.5× I/O speedup (Entry 032's claim on RAW files):
  - I/O time: 55 ms → 22 ms (saves 33 ms)
  - Total time: 4.425s → 4.392s
  - **Overall speedup: 1.007× (0.7% improvement)**

**Why Entry 032 Showed 2.5× Speedup**:
- Entry 032 tested **RAW file I/O** (no decompression) → 100% I/O-bound → 2.5× applies fully ✓
- biometal processes **compressed files** (with decompression) → 99% CPU-bound → 2.5× on 1% = negligible ✗

**Decision**: Keep implementation (no harm, future-proof), but document limited benefit

**See**: `RULE4_FINDINGS.md`

---

## REVISED Performance Reality

### Original Phase 2A Target (INCORRECT)

**Claimed** (from Entry 029 + Entry 032):
- Rule 3 (Parallel BGZF): 6.5× speedup
- Rule 4 (Smart mmap): 2.5× additional
- **Combined**: 16.3× speedup (6.5× × 2.5×)
- Sequential BAM: 55 MiB/s → **895 MiB/s**

### Actual Phase 2A Achievement (VALIDATED)

**Reality** (from multi-scale testing, N=10):
- Rule 3: 0.77-0.84× (SLOWER) → Disabled
- Rule 4: ~1% improvement (Amdahl's Law)
- **Combined**: ~1.3× improvement (removing parallel penalty + minimal mmap)
- Sequential BAM: 55 MiB/s → **~71 MiB/s**

**Discrepancy**: 16.3× claimed vs 1.3× achieved = **12.5× overestimation**

### Why Evidence Didn't Transfer

**Context dependency**:
1. Entry 029 used **all-at-once decompression** (unbounded memory) → biometal uses **bounded streaming** (constant memory)
2. Entry 032 tested **RAW file I/O** (100% I/O-bound) → biometal processes **compressed files** (99% CPU-bound)

**Lesson**: Same operations, different architectures → different results

---

## ACTUAL BOTTLENECK: Decompression (98.7% of time)

### Bottleneck Hierarchy (544 MB BAM file)

```
Total time: 4.42 s (100%)
├── Decompression: 4.37 s (98.7%) ← ACTUAL BOTTLENECK
│   └── CPU-bound (single-threaded flate2 decoder)
└── I/O: 55 ms (1.3%)
    └── Fast enough (SSD read)
```

### Why Current Optimizations Don't Help

**Rule 3 (Parallel BGZF)**: Conflicts with streaming architecture → 0.77× slowdown
**Rule 4 (Smart mmap)**: Optimizes I/O (1.3% of time) → negligible impact
**Rule 1 (NEON)**: Optimizes parsing (after decompression) → helps, but bottleneck remains

### What Would Actually Help

**Option 1: Alternative Decompression Libraries** (Highest Impact)

| Library | Claimed Speedup | Platform | Effort |
|---------|-----------------|----------|--------|
| **zlib-ng** | 2-3× faster | All | 20-30h |
| **libdeflate** | 1.5-2× faster | All | 15-20h |
| **igzip (ISA-L)** | 2-4× faster | x86_64 only | 25-35h |

**Potential Impact**:
- 2× decompression speedup = **1.97× overall** (98.7% of time improved)
- Much better than 1.3× from Rules 3+4
- Portable across all platforms

**Unknowns**:
- No evidence base (not in apple-silicon-bio-bench)
- Would require validation with N=30 benchmarks
- Integration complexity unknown

---

**Option 2: True All-At-Once Parallel (Violates Rule 5)** (Not Recommended)

If we abandon constant-memory streaming:
- Load entire file into memory (violates Rule 5)
- Decompress all blocks in parallel (Entry 029's approach)
- Expected: 6.5× speedup ✓
- **Trade-off**: Cannot handle TB-scale files (RAM limited)

**Decision**: Don't do this (streaming is core value proposition)

---

**Option 3: Accept Current Performance** (Strategic Choice)

Current sequential decompression (55-71 MiB/s) is:
- ✅ Competitive with samtools (~45-50 MiB/s)
- ✅ Constant memory (5 MB vs 20-50 MB)
- ✅ Reliable and portable

Focus effort on:
- Horizontal expansion (VCF, annotations)
- Rule 2 (if willing to invest in code duplication)
- Community building

---

## REVISED Strategic Priorities

### Tier 0: Decompression Bottleneck Investigation (OPTIONAL)

**IF** we want to significantly improve BAM parsing performance:

**Week 1-2: Investigate Alternative Decompression** (20-30 hours)
- Profile zlib-ng, libdeflate on ARM
- Benchmark with N=30 (validate 2-3× claims)
- Assess integration complexity
- Decide: Proceed or accept current performance

**Expected Outcome**:
- If successful: 2-3× improvement (55 → 110-165 MiB/s)
- If unsuccessful: Accept current performance, move to horizontal expansion

**Recommendation**: **Investigate** (20-30 hours is low risk for potentially high reward)

---

### Tier 1: High-ROI Vertical Optimization (Re-Evaluated)

**Rule 2: True Block Processing** (40-60 hours)

**Re-evaluation**:
- **Original thinking**: Lower priority than Rules 3+4 (14× vs 16.3×)
- **New reality**: Rules 3+4 provide ~1.3× (not 16.3×!)
- **New ranking**: Rule 2 now has MUCH higher ROI than Rules 3+4

**Decision factors**:
- Provides 14× speedup on CPU operations (validated)
- Affects ALL operations (FASTQ, BAM parsing, k-mer, quality)
- Requires ~1,200 lines of code duplication
- Maintenance burden: Keep two implementations in sync

**Recommendation**: **Re-evaluate priority** (now Tier 1 instead of deferred)

If decompression can't be improved → Rule 2 becomes highest-ROI vertical work

---

### Tier 1: Horizontal Expansion (Strategic Value)

**VCF/BCF Format** (60-80 hours) - **HIGH PRIORITY**

**Rationale**:
- Completes genomic workflow: FASTQ → BAM → VCF
- High user demand (variant calling is common)
- Natural extension of BAM work
- Clear use cases

**Implementation**:
- VCF streaming parser (20-25h)
- BCF binary format (25-30h)
- TBI index support (15-20h)
- Genotype operations (10-15h)

**Strategic Value**: "Complete genomic workflow toolkit"

---

**BED/GFF/GTF Parsers** (10-20 hours) - **QUICK WINS**

**Rationale**:
- Simple tab-delimited formats
- High utility (annotations, intervals)
- Low effort, high impact
- Enables many common workflows

**Implementation**:
- BED parser (3-5h)
- GFF3 parser (4-6h)
- GTF parser (3-5h)
- Interval operations (5-8h)

**Strategic Value**: Broad utility for minimal effort

---

**Alignment Analysis Operations** (60-80 hours) - **COMPLETENESS**

**Why Valuable**:
- Natural BAM extension
- Standard QC workflows
- Compete with samtools feature set

**Implementation**:
- Insert size distribution (15-20h)
- Coverage depth analysis (20-25h)
- Duplicate marking (20-25h)
- Alignment statistics (10-15h)

**Strategic Value**: Complete alignment analysis toolkit

---

### Tier 2: Advanced Features (Demand-Driven)

**K-mer Operations Expansion** (50-90 hours total)

**High Priority**:
- ✅ Canonical k-mers (5-10h) - **CORRECTNESS FIX** (should have been done originally)

**Medium Priority**:
- Streaming k-mer counting (30-40h)
- MinHash sketching (20-30h)
- K-mer classification (40-60h) - **IF** metagenomics demand

---

**Quality Control Operations** (40-60 hours)

- Adapter trimming (20-30h) - **IF** RNA-seq/ChIP-seq users request
- Error correction (20-30h) - **IF** assembly users request

---

### Tier 3: Low Priority (Wait for Demand)

**CRAM Format** (80-120 hours)
- Complex reference-based compression
- Potentially slower than BAM
- Uncertain demand
- **Decision**: Wait for explicit user requests

**Advanced Indexes** (40-60 hours)
- CSI index (30-40h)
- TBI index (20-30h) - May come with VCF work
- **Decision**: Implement if needed for specific formats

---

## CORRECTED Strategic Recommendations

### Phase 2A: Bottleneck Resolution (4-6 weeks) - OPTIONAL

**Week 1-2: Investigate Decompression Libraries** (20-30 hours)
- Benchmark zlib-ng, libdeflate on ARM
- Validate 2-3× speedup claims (N=30)
- Assess integration complexity
- **Go/No-Go Decision**: Proceed with integration or accept current performance

**Week 3-4: Implement Alternative Decompression** (20-30 hours) - **IF Go**
- Integrate chosen library
- Validate performance improvement
- Test across platforms (Mac ARM, Linux ARM, x86_64)
- **Expected Outcome**: 55 → 110-165 MiB/s (2-3× improvement)

**Week 5-6: Re-evaluate Rule 2** - **IF No-Go on decompression**
- If decompression can't be improved, Rule 2 becomes highest ROI
- Provides 14× CPU speedup (validated)
- Requires code duplication (~1,200 lines)
- **Decision point**: Accept maintenance cost for performance gain?

**Alternative**: Skip decompression investigation, accept current 55-71 MiB/s, move to Phase 2B

---

### Phase 2B: Strategic Horizontal Expansion (8-12 weeks)

**Week 1-4: VCF/BCF Format** (60-80 hours) - **HIGH PRIORITY**
- VCF streaming parser
- BCF binary format
- TBI index support
- Genotype operations
- **Why first**: Completes FASTQ → BAM → VCF workflow

**Week 5-6: Annotation Formats** (10-20 hours) - **QUICK WINS**
- BED parser
- GFF/GTF parsers
- Interval operations
- **Why second**: Low effort, high utility

**Week 7-10: Alignment Analysis** (60-80 hours) - **OPTIONAL**
- Insert size distribution
- Coverage depth analysis
- Duplicate marking
- Alignment statistics
- **Why optional**: Demand-driven, compete with samtools

**Week 11-12: K-mer Improvements** (20-40 hours) - **OPTIONAL**
- Canonical k-mers (correctness fix)
- Streaming k-mer counting (if demand)
- MinHash sketching (if demand)
- **Why optional**: Specialized use cases

---

### Phase 2C: Advanced Features (Demand-Driven)

**Based on community feedback**:
- CRAM format (80-120h) - IF explicitly requested
- K-mer classification (40-60h) - IF metagenomics demand
- Error correction (20-30h) - IF assembly users request
- Advanced QC (40-60h) - IF preprocessing demand

---

## Revised Priority Ranking

### Tier 0: Critical Decision (Do First)

| Action | Effort | Impact | ROI | Status |
|--------|--------|--------|-----|--------|
| **Investigate decompression libraries** | 20-30h | Potentially 2-3× | **HIGH** | **RECOMMENDED** |
| **Alternative: Accept current perf** | 0h | Focus on horizontal | N/A | **VALID** |

**Rationale**:
- Low effort (20-30h) for potentially high reward (2-3×)
- Addresses actual bottleneck (98.7% of time)
- If unsuccessful, accept current performance and move on

---

### Tier 1: High Priority (Do Next)

| Feature | Effort | Impact | ROI | Timing |
|---------|--------|--------|-----|--------|
| **VCF/BCF format** | 60-80h | High (workflow completion) | **HIGH** | Weeks 1-4 |
| **BED/GFF/GTF** | 10-20h | Medium (utility) | **HIGH** | Weeks 5-6 |
| **Canonical k-mers** | 5-10h | Medium (correctness) | **HIGH** | Week 7 |
| **Rule 2 (if decompression fails)** | 40-60h | High (14× CPU) | **MEDIUM** | Re-evaluate |

**Rationale**:
- VCF completes genomic workflow (FASTQ → BAM → VCF)
- Annotations are quick wins (low effort, high utility)
- Canonical k-mers fix correctness issue
- Rule 2 now higher ROI than originally thought (Rules 3+4 failed)

---

### Tier 2: Medium Priority (Opportunistic)

| Feature | Effort | Impact | ROI | Timing |
|---------|--------|--------|-----|--------|
| Alignment analysis | 60-80h | Medium (completeness) | MEDIUM | Weeks 7-10 |
| Streaming k-mer counting | 30-40h | Medium (specialized) | MEDIUM | Month 3+ |
| MinHash sketching | 20-30h | Medium (specialized) | MEDIUM | Month 3+ |
| Adapter trimming | 20-30h | Medium (preprocessing) | MEDIUM | Month 4+ |

**Rationale**:
- Useful but not critical
- Implement based on user feedback
- Moderate effort for moderate payoff

---

### Tier 3: Low Priority (Wait for Demand)

| Feature | Effort | Impact | ROI | Timing |
|---------|--------|--------|-----|--------|
| CRAM format | 80-120h | Medium (uncertain) | LOW | If demanded |
| K-mer classification | 40-60h | High (specialized) | LOW | If demanded |
| Error correction | 20-30h | Medium (specialized) | LOW | If demanded |
| Advanced indexes | 40-60h | Low (niche) | LOW | If demanded |

**Rationale**:
- High effort or specialized use cases
- Uncertain demand
- Can be added later without major refactoring

---

## Updated Competitive Position

### Current State (v1.6.0)

**Strengths**:
- ARM NEON optimization (4-25× current operations)
- Constant memory (10-200× lower than competitors)
- BAM/SAM with BAI index support
- Evidence-based methodology (catches false assumptions!)
- Competitive BAM parsing (55 MiB/s vs samtools 45-50 MiB/s)

**Gaps**:
- No variant format support (VCF/BCF)
- Limited annotation format support
- Decompression bottleneck (98.7% of time, hard to improve)

---

### After Phase 2B (VCF + Annotations)

**Position**: "Complete ARM-native genomic workflow toolkit"

**Differentiation**:
- 4-25× NEON speedups (unique on ARM)
- Complete workflow: FASTQ → BAM → VCF → annotations
- Constant memory across all operations
- Evidence-based performance claims (validated, not exaggerated)

**Competitive Matrix**:

| Feature | biometal (current) | biometal (after 2B) | samtools | HTSlib |
|---------|-------------------|---------------------|----------|---------|
| BAM parsing | 55 MiB/s | 55-165 MiB/s* | 45-50 MiB/s | 45-50 MiB/s |
| ARM NEON | 4-25× | 4-25× | ✗ | ✗ |
| Memory | 5 MB | 5 MB | 20-50 MB | 20-50 MB |
| VCF/BCF | ✗ | **✓** | ✓ | ✓ |
| Annotations | ✗ | **✓** | ✓ | ✓ |
| Streaming | ✓ | ✓ | Partial | Partial |

\* Depends on decompression library investigation outcome

---

### Honest Positioning Statement

**biometal is**:
- ✅ ARM-native with NEON optimizations (4-25× on specific operations)
- ✅ Constant-memory streaming (10-200× lower than alternatives)
- ✅ Evidence-based (catches false assumptions, reports honest results)
- ✅ Competitive BAM parsing (55 MiB/s, comparable to samtools)
- ✅ Complete workflow toolkit (after VCF/annotations)

**biometal is NOT**:
- ❌ 16× faster than samtools at BAM parsing (Rules 3+4 didn't work)
- ❌ Going to achieve massive parallel speedups (conflicts with streaming)
- ❌ Trying to replace HTSlib comprehensively (focused toolkit)

**Unique value**: Evidence-based ARM-native toolkit with streaming architecture and honest performance claims

---

## Lessons Learned: Evidence-Based Optimization

### What Worked

✅ **Multi-scale testing** revealed Rule 3 fails at all scales (not just small files)
✅ **Bottleneck profiling** identified decompression (98.7%) as actual problem
✅ **DAG framework** prevented wasted effort (pruned Rule 3 at <1.5× threshold)
✅ **Amdahl's Law** explained why Rule 4 provides minimal benefit
✅ **Domain analysis** caught context dependency (Entry 029/032 vs biometal)

### What Didn't Work

❌ **Assuming evidence transfers** without validating context (unbounded vs bounded streaming)
❌ **Targeting wrong bottleneck** (optimized I/O when decompression dominates)
❌ **Combined speedup math** (6.5× × 2.5× = 16.3× assumed both apply fully)

### Key Insights

**1. Context Matters More Than Numbers**
- Entry 029: 6.5× (all-at-once) ≠ biometal: 0.77× (bounded streaming)
- Same operation, different architecture → different results

**2. Bottleneck Analysis is Critical**
- Assumed I/O bottleneck → Actually CPU bottleneck
- Rules 3+4 optimize I/O (1.3% of time) → Negligible impact

**3. Architecture Trade-offs are Real**
- Rule 3 (parallelism) vs Rule 5 (streaming) → Can't have both
- Chose streaming (TB-scale capability) over speed

**4. DAG Methodology Saves Time**
- Multi-scale testing (10 hours) prevented 40-60 hours on failed optimization
- Pruning criteria (<1.5×) provides clear go/no-go decisions

**5. Honest Reporting Builds Credibility**
- Reporting 1.3× achieved (not 16.3× claimed) demonstrates scientific rigor
- Negative results are valuable (document what doesn't work)

---

## Timeline Summary (Revised)

### Decision Point (Week 1-2): Investigate Decompression

**Work**: Benchmark zlib-ng, libdeflate
**Effort**: 20-30 hours
**Outcome**: 2-3× speedup OR accept current performance
**Status**: **RECOMMENDED** (low risk, potentially high reward)

---

### Path A: IF Decompression Improves (Weeks 3-4)

**Work**: Integrate alternative library
**Effort**: 20-30 hours
**Outcome**: 55 → 110-165 MiB/s
**Status**: Continue to Phase 2B (horizontal expansion)

---

### Path B: IF Decompression Doesn't Improve (Weeks 3-4)

**Work**: Re-evaluate Rule 2 (14× CPU speedup)
**Effort**: 40-60 hours (if proceed)
**Outcome**: 14× on CPU operations, but requires code duplication
**Status**: Decision point - accept maintenance burden?

---

### Phase 2B: Horizontal Expansion (Weeks 5-14)

**Work** (in priority order):
1. VCF/BCF format (60-80h) - Complete workflow
2. BED/GFF/GTF (10-20h) - Quick wins
3. Canonical k-mers (5-10h) - Correctness fix
4. Alignment analysis (60-80h) - IF demand
5. K-mer improvements (20-40h) - IF demand

**Total Effort**: 155-250 hours over 10 weeks
**Outcome**: Complete genomic workflow toolkit
**Status**: Prioritize based on feedback

---

## Updated Success Metrics

### Phase 2A Success Criteria (IF Pursued)

**Decompression Investigation**:
- ✓ Benchmarked zlib-ng, libdeflate (N=30)
- ✓ Validated 2-3× claims (or documented failure)
- ✓ Clear go/no-go decision made
- ✓ Integration complexity assessed

**IF Successful**:
- ✓ Alternative library integrated
- ✓ 2-3× speedup validated
- ✓ Cross-platform testing complete
- ✓ No memory regression

---

### Phase 2B Success Criteria

**VCF/BCF Format**:
- ✓ VCF streaming parser working
- ✓ BCF binary parser working
- ✓ TBI index support
- ✓ Genotype operations
- ✓ Round-trip conversion validated

**Annotation Formats**:
- ✓ BED parser working
- ✓ GFF/GTF parsers working
- ✓ Interval operations
- ✓ Integration tests

**K-mer Improvements**:
- ✓ Canonical k-mers implemented
- ✓ All existing tests pass
- ✓ Correctness validated

---

### Overall Success Definition

**Technical**:
- ✓ Decompression bottleneck addressed (or explicitly accepted)
- ✓ Complete FASTQ → BAM → VCF workflow
- ✓ Constant memory maintained (5 MB)
- ✓ All tests passing (target: 700+)
- ✓ Evidence-based claims (no false promises)

**Strategic**:
- ✓ Competitive with samtools on features
- ✓ Superior on ARM NEON performance
- ✓ Honest, validated performance claims
- ✓ Community-ready toolkit

---

## Final Recommendations

### 1. IMMEDIATE (Weeks 1-2): Investigate Decompression

**DO THIS FIRST**: 20-30 hours to benchmark alternative libraries
- Low risk, potentially high reward (2-3× speedup)
- Addresses actual bottleneck (98.7% of time)
- Clear decision point after investigation

**Options after investigation**:
- **IF successful**: Integrate library (20-30h more)
- **IF unsuccessful**: Accept current performance, move to horizontal expansion

---

### 2. PRIMARY PATH (Weeks 3-12): Horizontal Expansion

**FOCUS HERE** (regardless of decompression outcome):

**Week 3-6: VCF/BCF Format** (60-80h)
- Highest priority horizontal work
- Completes genomic workflow
- High user demand

**Week 7-8: Annotation Formats** (10-20h)
- Quick wins
- Broad utility
- Low effort

**Week 9-10: Canonical K-mers** (5-10h)
- Correctness fix
- Should have been done originally
- Quick implementation

**Week 11-14: Demand-Driven** (variable)
- Alignment analysis IF requested
- K-mer improvements IF metagenomics demand
- QC operations IF preprocessing demand

---

### 3. DEFERRED: Rule 2 Block Processing

**Re-evaluate** IF:
- Decompression can't be improved (bottleneck remains)
- Community explicitly requests massive CPU speedups
- Willing to accept maintenance burden (~1,200 lines duplication)

**Don't prioritize** IF:
- Focusing on horizontal expansion
- Prefer feature breadth over vertical optimization
- Want to minimize code complexity

---

## Document Status

**Created**: November 11, 2025 (Original)
**Revised**: November 11, 2025 (Post-Rules 3+4 Investigation)
**Version**: 2.0 (Major revision with evidence-based corrections)
**Status**: Strategic planning with validated priorities
**Next Review**: After decompression investigation (Week 2)
**Owner**: biometal core development

---

## References

- `OPTIMIZATION_RULES.md` - Evidence base (corrected for Rules 3+4)
- `RULE2_INVESTIGATION_FINDINGS.md` - Block processing analysis
- `RULE3_BENCHMARK_RESULTS.md` - Multi-scale parallel BGZF testing
- `RULE4_FINDINGS.md` - Bottleneck analysis and Amdahl's Law
- `RULE3_AND_RULE4_SESSION_SUMMARY.md` - Complete investigation narrative
- `CLAUDE.md` - Development guide
- apple-silicon-bio-bench - Evidence base (1,357 experiments)

---

**Remember**: Evidence first, implementation second. Profile before optimizing. Benchmark with N≥10. **Report honest results**, even when they contradict expectations. Negative results are valuable.

**Key Takeaway**: Systematic testing revealed that Phase 2A's original 16.3× target is not achievable with current architecture (actual: ~1.3×). Focus shifted to addressing actual bottleneck (decompression) or accepting current performance and expanding horizontally.
