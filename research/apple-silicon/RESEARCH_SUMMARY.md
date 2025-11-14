# Apple Silicon Research Archive

**Period**: November 4-13, 2025 (2 weeks)
**Status**: Archived - Returning to Core Roadmap
**Decision**: Option A - Focus on Rules 3+4 (16× proven speedup)

---

## Executive Summary

Explored Apple Silicon hardware acceleration (Neural Engine + GPU) for biometal use cases. Built production-quality infrastructure but discovered strategic misalignment with biometal's streaming-first architecture and evidence-based optimization principles.

**Outcome**: Archive as valuable research, return to core roadmap (Phase 1 consolidation + Rules 3+4 implementation).

---

## What Was Built

### 1. Neural Engine Integration (Week 1)
**Location**: `neural-engine/`
- ✅ Complete ONNX Runtime + CoreML integration
- ✅ Quality prediction model (PyTorch → ONNX)
- ✅ Production-ready inference pipeline
- ❌ 2,940× **slowdown** for streaming use case

**Finding**: Neural Engine optimized for batch inference, not per-read streaming

### 2. GPU Smith-Waterman (Week 2)
**Location**: `gpu-smith-waterman/`
- ✅ Metal compute shader (340 lines)
- ✅ Rust GPU dispatch (430 lines)
- ✅ All tests passing (430 tests)
- ⚠️ 1.2-1.4× speedup for batches ≥10
- ❌ Slower than CPU for single alignments

**Finding**: Current implementation achieves 1.2-1.4× vs 10-50× from literature. Would need anti-diagonal parallelization for claimed speedups.

---

## Strategic Assessment

### Core Tension

**biometal's strengths**:
- Streaming-first (constant memory, terabyte-scale)
- Evidence-based (ASBB: 1,357 experiments, N=30)
- Cross-platform (Mac ARM → AWS Graviton → x86_64)

**Apple Silicon requirements**:
- Batch-oriented processing
- Platform-specific (Mac-only)
- Speculative optimization (literature vs measured)

### Results vs Effort

| Component | Time | Result | Strategic Fit |
|-----------|------|--------|---------------|
| Neural Engine | 1 week | 0.0003× (slowdown) | ❌ Mismatch |
| GPU Smith-Waterman | 1 week | 1.2-1.4× speedup | ⚠️ Modest |
| **Total** | **2 weeks** | **Not compelling** | **Low** |

**Opportunity cost**: Rules 3+4 (16× proven speedup) still unimplemented

---

## Why Return to Core Roadmap

### Evidence-Based Decision

**Rules 3+4 offer**:
- 16× **proven** speedup (ASBB entry 029, 032)
- Cross-platform (benefits all users)
- Aligns with streaming architecture
- Completes Phase 1 (all 6 rules)

**Current Phase 1 status**:
- 4/6 rules implemented (67%)
- Rules 3+4 remaining (parallel BGZF + smart mmap)
- Expected: 55 MiB/s → 880 MiB/s BAM parsing

### What We Learned

**Apple Silicon strengths** (for different use cases):
- Batch classification (adapter detection)
- ML inference pipelines (variant calling)
- Large-scale read binning
- Image processing (not biometal's domain)

**biometal's sweet spot**:
- Streaming architecture
- Cross-platform portability
- Evidence-based optimization
- Constant-memory processing

**Mismatch**: Apple Silicon hardware (batch-oriented) vs biometal architecture (streaming-first)

---

## Valuable Outcomes

### Infrastructure Built

**Reusable for future**:
1. ✅ ONNX Runtime integration (`src/ml/`)
2. ✅ Metal GPU framework (`src/alignment/gpu/`)
3. ✅ Feature flag patterns (`gpu`, `neural-engine`)
4. ✅ Comprehensive testing methodology

**If better use cases identified**:
- Adapter detection (batch classification)
- Read quality binning (bulk processing)
- Variant calling assistance (ML inference)

### Code Quality

**All implementations production-ready**:
- ✅ Comprehensive testing (property tests)
- ✅ Full documentation (API docs + guides)
- ✅ Performance benchmarks (N=30 rigor)
- ✅ Error handling (no panics)

**Can be activated with feature flags**:
```bash
cargo build --features gpu           # GPU Smith-Waterman
cargo build --features neural-engine # Neural Engine inference
```

---

## Lessons Learned

### Technical

1. **Hardware-software fit matters**
   - Neural Engine: batch inference, not streaming
   - GPU: batch parallelism, not single-alignment

2. **Literature claims require validation**
   - CUDA 10-50× speedup requires anti-diagonal parallelization
   - Our 1.2-1.4× uses simpler batch-parallel approach

3. **Platform-specific code has costs**
   - Breaks portability promise
   - Smaller user base (Mac-only vs cross-platform)
   - Maintenance overhead

### Strategic

1. **Evidence-based methodology is powerful**
   - ASBB experiments guide correct optimizations
   - Proven speedups (16×) > speculative claims (10-50×)

2. **Opportunity cost is real**
   - 2 weeks on Apple Silicon
   - Rules 3+4 still waiting (16× proven)

3. **Core strengths should guide priorities**
   - biometal = streaming + evidence-based + cross-platform
   - Apple Silicon = batch + speculative + platform-specific
   - Mismatch → deprioritize

---

## Future Considerations

### When to Revisit Apple Silicon

**Good use cases** (batch-oriented):
- Adapter/contamination detection (batch classification)
- Read binning (categorize millions of reads)
- Variant calling assistance (batch ML inference)
- Quality control pipelines (batch processing)

**Poor use cases** (streaming):
- Per-read quality prediction ❌ (proven)
- Single-alignment scoring ❌ (GPU overhead)
- Streaming sequence operations ❌ (NEON better)

### Potential Future Work

**If compelling use case emerges**:
1. GPU anti-diagonal parallelization (10-50× potential)
2. Neural Engine adapter detection (batch classification)
3. Vulkan backend (cross-platform GPU)

**Priority**: After Rules 3+4 + Phase 1 completion

---

## Archive Contents

### Documents
- `RESEARCH_SUMMARY.md` (this file)
- `STRATEGIC_DECISION_NOV_2025.md` - Decision rationale
- `STRATEGIC_PIVOT_PLAN.md` - Original exploration plan

### GPU Smith-Waterman
- `gpu-smith-waterman/SMITH_WATERMAN_GPU_COMPLETE.md` - Full implementation report
- `gpu-smith-waterman/GPU_INTEGRATION_PLAN.md` - Integration plan

### Neural Engine
- `neural-engine/NEURAL_ENGINE_STRATEGIC_DISCUSSION.md` - Strategic analysis

### Code (still in codebase, feature-gated)
- `src/alignment/gpu/` - GPU implementation (feature: `gpu`)
- `src/ml/` - Neural Engine (feature: `neural-engine`)
- `benches/smith_waterman.rs` - GPU benchmarks
- `benches/neural_engine_quality.rs` - Neural Engine benchmarks

---

## Recommendation: Return to Core Roadmap

**Next priorities**:
1. ✅ Archive Apple Silicon research (this document)
2. ✅ Update CLAUDE.md (remove Apple Silicon focus)
3. 🔄 Implement Rules 3+4 (parallel BGZF + smart mmap)
4. 🔄 Complete Phase 1 consolidation
5. 🔄 Community launch (blog post, social media)

**Expected impact**:
- Rules 3+4: 55 MiB/s → 880 MiB/s BAM parsing (16× proven)
- Phase 1 complete: All 6 rules implemented (27× combined)
- Community-ready: Strong performance story for v1.8.0

**Timeline**: 3-4 weeks for Rules 3+4 implementation

---

## Conclusion

**Apple Silicon exploration**: Valuable learning experience
**Strategic fit**: Poor alignment with biometal's core strengths
**Decision**: Return to evidence-based roadmap (Rules 3+4)

**Value preserved**:
- ✅ Infrastructure reusable for future use cases
- ✅ Lessons learned guide better decisions
- ✅ Code quality demonstrates capabilities

**Moving forward**: Focus on proven optimizations that deliver transformative performance for all users, not just Mac users.

---

**Status**: Archived November 13, 2025
**Next**: Phase 1 consolidation + Rules 3+4 implementation
