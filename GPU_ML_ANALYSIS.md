# GPU/ML Work Analysis for biometal

**Date**: November 15, 2025
**Context**: Evaluating whether to continue with GPU/ML or complete format coverage

---

## What Has Already Been Done (November 4-13, 2025)

### ✅ Complete 2-Week Strategic Pivot

**Week 1: Neural Engine**
- Built complete ONNX Runtime + CoreML integration
- Trained quality prediction model (PyTorch → ONNX)
- Production-ready inference pipeline
- **Result**: ❌ **2,940× SLOWDOWN** for streaming use case
- **Finding**: Neural Engine optimized for batch inference, not per-read streaming

**Week 2: GPU Smith-Waterman**
- Metal compute shader (340 lines)
- Rust GPU dispatch (430 lines)
- All tests passing (430 tests)
- **Result**: ⚠️ **1.2-1.4× speedup** (not 10-50× from literature)
- **Finding**: Would need anti-diagonal parallelization for claimed speedups

### 📦 Infrastructure Preserved (Feature-Gated)

**In Codebase Now**:
- `src/ml/` - ONNX Runtime + Neural Engine integration
- `src/alignment/gpu/` - GPU Smith-Waterman (Metal)
- `src/alignment/metal/` - Metal framework code
- Feature flags: `gpu` and `neural-engine`

**Quality**: Production-ready
- ✅ Comprehensive testing (property tests)
- ✅ Full documentation (API docs + guides)
- ✅ Performance benchmarks (N=30 rigor)
- ✅ Error handling (no panics)

**Activation**:
```bash
cargo build --features gpu           # GPU Smith-Waterman
cargo build --features neural-engine # Neural Engine inference
```

---

## Strategic Decision (November 13, 2025)

**DECISION**: Archive GPU/ML, return to core infrastructure

### Why Archived

1. **Performance Mismatch**:
   - Neural Engine: 2,940× slowdown (batch vs streaming)
   - GPU: 1.2-1.4× speedup (vs 10-50× claimed in literature)
   - Not worth platform-specific (Mac-only) complexity

2. **Architectural Conflict**:
   - biometal = streaming-first (constant memory)
   - Apple Silicon = batch-oriented processing
   - Fundamental mismatch

3. **Better Alternatives Available**:
   - **Rule 3** (Parallel BGZF): Tested and **FAILED** (0.77-0.84× slowdown)
     - Conflicts with streaming architecture
     - Chunking overhead compounds with scale
   - **Rule 4** (Smart mmap): **2.5× speedup** for files ≥50 MB
     - APFS prefetching advantage
     - Evidence-based, validated

4. **Opportunity Cost**:
   - 2 weeks invested in GPU/ML
   - Rule 4 (smart mmap) still unimplemented
   - Format coverage incomplete (CRAM, CSI, BCF)

---

## What PROJECT_TODOS.md Proposes (6-Month GPU/ML Plan)

**Original Plan**: 24 weeks (6 months) of GPU/ML work

### Phase 1: GPU/Metal (Weeks 1-8)
- Week 1-2: Smith-Waterman GPU ✅ **ALREADY DONE**
- Week 3: Pileup generation GPU (not started)
- Week 4: GPU results analysis ✅ **ALREADY DONE**
- Week 5-6: Neural Engine quality prediction ✅ **ALREADY DONE**
- Week 7: Alignment primitives (CPU/NEON)
- Week 8: Variant calling primitives

### Phase 2: ML/BERT (Weeks 9-16)
- Week 9: Streaming BERT data loaders
- Week 10: Quality-aware tokenization
- Week 11: GPU-accelerated tokenization
- Week 12: Example models + documentation
- Week 13: Multi-modal data loaders
- Week 14: Neural Engine adapter detection
- Week 15: Training utilities
- Week 16: Format support (BED, GFF/GTF) ✅ **ALREADY DONE**

### Phase 3: Demonstrations (Weeks 17-24)
- Week 17: Variant calling GPU
- Week 18-19: Assembly primitives
- Week 20: AMX RNA-seq matrices
- Week 21-24: Case studies and publications

**Total Effort**: 480-720 hours (12-18 full-time weeks)

---

## Critical Analysis: Should We Continue?

### ❌ Arguments AGAINST Continuing GPU/ML Work

1. **Already Tested and Failed**:
   - Neural Engine: 2,940× slowdown (not suitable)
   - GPU Smith-Waterman: 1.2-1.4× (not compelling)
   - Would need anti-diagonal parallelization for claimed speedups

2. **Architectural Mismatch**:
   - biometal's core value = streaming (constant memory)
   - GPU/ML requires batching (defeats streaming)
   - Fundamental conflict

3. **Rule 3 Already Failed**:
   - Parallel BGZF: 0.77-0.84× slowdown (tested Nov 11)
   - Chunking overhead compounds
   - Won't achieve 6.5× claimed speedup with streaming

4. **Platform Lock-In**:
   - GPU/Neural Engine = Mac-only
   - Breaks cross-platform promise (Graviton, x86_64)
   - Smaller user base

5. **Time Investment vs Return**:
   - 2 weeks already spent on GPU/ML
   - Modest results (1.2-1.4× GPU, 0.0003× Neural Engine)
   - **Better alternatives**: Rule 4 (2.5× proven), CRAM reader (critical gap)

6. **PROJECT_TODOS.md is Outdated**:
   - Written before strategic pivot (Nov 4-13)
   - Assumes GPU/ML would work well
   - Evidence shows otherwise

### ⚠️ Arguments FOR Continuing (Weak)

1. **Infrastructure Already Built**:
   - Code is production-ready
   - Could be used for different operations
   - But: No compelling use cases identified

2. **Potential for Better GPU Algorithms**:
   - Anti-diagonal Smith-Waterman could hit 10-50× claims
   - But: Requires significant additional work (7-10 days)
   - Still Mac-only

3. **ML Use Cases Unexplored**:
   - Adapter detection (batch classification)
   - Read quality binning
   - But: Speculative, no evidence these would work well

---

## What Actually Needs Doing (Evidence-Based)

### Option A: Complete Format Coverage (Recommended) ✅

**Critical Gaps**:
1. **CRAM Reader** (80-120h, 2-3 weeks) - HIGHEST PRIORITY
   - 1000 Genomes uses CRAM exclusively
   - 3-5× smaller than BAM
   - Blocker for modern datasets

2. **CSI Index** (20-30h, 3-5 days) - Quick win
   - Supports >2GB chromosomes (BAI limited to 512MB)
   - Partially implemented, needs finishing

3. **BCF Format** (30-40h, 5-7 days)
   - Binary VCF with compression
   - Common in pipelines

**Total**: 4-5 weeks to complete all format coverage

**Then**: biometal = complete format library (foundational)

### Option B: Optimize Existing Code (Alternative) ✅

**Rule 4: Smart mmap** (Already validated)
- **Effort**: 10-20 hours (1-2 days)
- **Speedup**: 2.5× for files ≥50 MB
- **Evidence**: Entry 032 (scale validated)
- **Platform**: macOS (validated), Linux (future)
- **Complexity**: Low (threshold-based I/O selection)

**Compression Optimization** (Already done!)
- ✅ cloudflare_zlib: 1.67× decompression (v1.7.0)
- Nothing more to do here

### Option C: Continue GPU/ML (NOT Recommended) ❌

**Why not**:
- Already tested (2 weeks)
- Results not compelling (1.2-1.4× GPU, 0.0003× Neural Engine)
- Architectural mismatch (batching vs streaming)
- Rule 3 failed (0.77-0.84× slowdown)
- Better alternatives available (Rule 4, CRAM)

**What would be needed**:
- Anti-diagonal Smith-Waterman GPU (7-10 days)
- Different ML use cases (adapter detection, 15-25 weeks)
- Risk: Still might not work well

---

## Recommendation: Path Forward

### Path A: Complete Core Infrastructure (RECOMMENDED) ✅

**Phase 1**: Format Coverage (4-5 weeks)
1. CSI Index (3-5 days) - Quick win
2. CRAM Reader (2-3 weeks) - Critical gap
3. BCF Format (5-7 days) - Binary VCF

**Phase 2**: Optimization (1-2 days)
4. Rule 4 (smart mmap) - 2.5× speedup (validated)

**Total**: 5-6 weeks to complete foundational work

**Result**:
- ✅ All major formats supported (CRAM unblocks 1000 Genomes)
- ✅ Rule 4 implemented (2.5× speedup)
- ✅ Solid foundation for production use
- ✅ Then ready for GPU/ML IF compelling use cases emerge

### Path B: Continue GPU/ML (NOT RECOMMENDED) ❌

**Why avoid**:
- Already invested 2 weeks
- Results not compelling
- Architectural mismatch
- Better ROI elsewhere

**If insisting on GPU/ML**:
1. Identify compelling batch use case (not streaming)
2. Accept Mac-only limitation
3. Implement anti-diagonal Smith-Waterman (7-10 days)
4. Test with real workloads
5. Risk: Still might fail

---

## Summary

### What We Know:
- ✅ GPU/ML infrastructure built (2 weeks)
- ✅ Tested and archived (results not compelling)
- ❌ Rule 3 (Parallel BGZF) failed (0.77-0.84× slowdown)
- ✅ Rule 4 (smart mmap) validated (2.5× speedup)
- ❌ CRAM, CSI, BCF still missing (critical gaps)

### What We Should Do:
1. **Complete format coverage** (4-5 weeks)
   - CRAM Reader (critical)
   - CSI Index (quick win)
   - BCF Format (standard)
2. **Implement Rule 4** (1-2 days)
   - 2.5× speedup (validated)
   - Low effort, high return
3. **Archive GPU/ML work** (already done)
   - Infrastructure preserved
   - No active development
   - Available if better use cases emerge

### What We Should NOT Do:
- ❌ Continue GPU/ML work from PROJECT_TODOS.md
- ❌ Spend 6 months on speculative GPU/ML features
- ❌ Ignore critical format gaps (CRAM)

---

**Conclusion**: PROJECT_TODOS.md is outdated. Evidence shows format coverage + Rule 4 is the high-ROI path.
