# Strategic Pivot Summary: GPU + Neural Engine Exploration

**Duration**: Weeks 1-2 of 24-week exploration plan
**Date**: November 12-13, 2025
**Objective**: Evaluate GPU and Neural Engine acceleration for bioinformatics operations
**Status**: ✅ COMPLETE (2/2 weeks)

---

## Executive Summary

Explored GPU (Metal) and Neural Engine (CoreML) acceleration for biometal operations:

**Week 1 - GPU Acceleration**:
- ✅ **Success**: 771× speedup for Smith-Waterman alignment
- Target: Sequence alignment, batch processing
- Technology: Metal compute shaders, Apple Silicon GPU
- Result: **Production-ready, highly recommended**

**Week 2 - Neural Engine**:
- ⚠️ **Mixed Results**: 2,940× slowdown for simple quality filtering
- Target: ML-powered read quality prediction
- Technology: ONNX Runtime + CoreML, Apple Neural Engine
- Result: **Infrastructure ready, but simple use case inappropriate**

---

## Week 1: GPU Acceleration (Smith-Waterman)

### Achievements ✅

**Performance**:
- Sequential CPU: 14.9 µs per alignment
- CPU Batch (16 alignments): 14.1 µs per alignment (1.06× speedup)
- **GPU Batch (1,024 alignments)**: 19.3 ns per alignment (**771× speedup**)

**Implementation**:
- Metal compute shaders for parallel alignment
- Batched processing (1,024 alignments/batch)
- Complete Smith-Waterman with traceback
- Production-ready code (320 lines)

**Documentation**:
- `research/smith-waterman-gpu/GPU_BENCHMARK_RESULTS.md`
- `research/smith-waterman-gpu/IMPLEMENTATION_NOTES.md`
- Comprehensive benchmarks (N=30)

### Key Findings

1. **GPU Excels at Batch Processing**: 771× speedup with 1,024-alignment batches
2. **Overhead Matters**: Small batches (<16) slower than CPU due to GPU dispatch
3. **Sweet Spot**: Batches of 256-1,024 alignments ideal
4. **Use Cases**: Multiple sequence alignment, database search, variant calling

### Recommendation: **ADOPT** ✅

GPU Smith-Waterman is production-ready and delivers massive speedups for batch workflows.

---

## Week 2: Neural Engine (Quality Prediction)

### Achievements ✅

**Infrastructure**:
- ONNX Runtime v2.0 integration with CoreML
- Neural Engine dispatch working
- Complete training pipeline (PyTorch → ONNX)
- Benchmarking framework
- 1,640 lines of production code
- 1,900 lines of documentation

**Performance (Quality Filtering)**:
- Traditional (NEON): 4.4 ns per read
- **Neural Engine: 13.2 µs per read (2,940× slower)**

### Key Findings

1. **Overhead-Dominated**: Feature encoding + ONNX dispatch >> inference time
2. **Use Case Inappropriate**: Simple operations don't benefit from ML
3. **Baseline Too Fast**: NEON-optimized quality filtering (25× scalar) hard to beat
4. **Infrastructure Valuable**: Framework supports any ONNX model (future-ready)

### Better Use Cases (Not Tested)

Neural Engine would excel at:
- **Adapter detection**: Complex pattern recognition
- **Read classification**: Multiclass ML (human/bacterial/viral)
- **Quality prediction**: Predict from sequence alone (no quality scores)
- **Variant calling**: ML-assisted variant detection

### Recommendation: **ARCHIVE** ⚠️

Current use case (quality filtering) inappropriate. Infrastructure ready for future ML tasks.

---

## Strategic Assessment

### What We Learned

1. **GPU**: Massive speedups for embarrassingly parallel operations (771×)
2. **Neural Engine**: Overhead dominates simple operations (2,940× slower)
3. **Use Case Selection**: Critical for performance gains
4. **Infrastructure**: Both frameworks production-ready

### Performance Summary

| Technology | Use Case | Result | Recommendation |
|------------|----------|--------|----------------|
| **GPU (Metal)** | Smith-Waterman alignment | **771× speedup** | ✅ **ADOPT** |
| **Neural Engine** | Simple quality filtering | **2,940× slowdown** | ⚠️ **ARCHIVE** |
| **Neural Engine** | Complex ML tasks | **Not tested** | 🔮 **Future potential** |

### Time Investment vs. Value

**GPU Week 1**:
- Time: ~20 hours development + 5 hours benchmarking
- Result: Production-ready 771× speedup
- Value: **Extremely high** (immediate impact)

**Neural Engine Week 2**:
- Time: ~30 hours infrastructure + 10 hours benchmarking
- Result: Infrastructure ready, use case inappropriate
- Value: **Medium** (future-ready, but no immediate use)

---

## Recommendations

### Immediate Actions (Post Week 2)

1. **Integrate GPU Smith-Waterman into biometal**:
   - Move from `research/` to `src/alignment/gpu.rs`
   - Add feature flag: `gpu` (Metal, macOS only)
   - Document batch API
   - Add to benchmarks
   - Update CHANGELOG

2. **Archive Neural Engine Research**:
   - Document findings (COMPLETE)
   - Preserve infrastructure for future ML tasks
   - Do NOT pursue simple quality filtering
   - Consider complex ML use cases in Phase 3

3. **Return to Phase 1 Core Development**:
   - Week 3: Community building
   - Week 4: Quality assurance
   - Weeks 5-9: Rules 3+4 (parallel BGZF, mmap)

### Future Considerations (Phase 3+)

**Potential Neural Engine Use Cases** (NOT core development):
- Adapter/contamination detection
- Read classification (species identification)
- Basecalling (sequence-only quality prediction)
- Variant calling assistance

**Decision Point**: These are **NEW features** beyond current roadmap. Would require:
- Scope expansion discussion
- User demand validation
- Resource allocation assessment
- Strategic priority evaluation

---

## Archive Contents

### Week 1: GPU Acceleration ✅

```
research/smith-waterman-gpu/
├── GPU_BENCHMARK_RESULTS.md        # Performance analysis
├── IMPLEMENTATION_NOTES.md         # Technical details
├── shaders/
│   └── smith_waterman.metal        # Metal compute shader (320 lines)
└── benches/
    └── smith_waterman.rs           # Benchmarks (N=30)
```

**Status**: Ready for production integration

### Week 2: Neural Engine ✅

```
research/neural-engine/
├── RESEARCH_NOTES.md               # Neural Engine capabilities (500+ lines)
├── README.md                       # Training workflow (350+ lines)
├── DAYS_3_4_SUMMARY.md            # Progress report
├── BENCHMARK_RESULTS.md            # Performance analysis (400+ lines)
├── WEEK_2_COMPLETE.md             # Summary report
├── train_quality_model.py          # PyTorch training (450 lines)
├── test_onnx_inference.py          # ONNX testing (150 lines)
├── generate_training_data.py       # Synthetic data (200 lines)
├── create_minimal_model.py         # Minimal model (120 lines)
├── requirements.txt                # Python dependencies
├── quality_model.onnx              # ONNX model (1.4 KB)
└── training_data.fq.gz             # Synthetic FASTQ (8 MB)
```

**Status**: Infrastructure complete, archived for future ML exploration

---

## Strategic Pivot: Final Verdict

### Continue: GPU Smith-Waterman ✅

**Rationale**:
- 771× speedup is transformative
- Production-ready implementation
- Clear use cases (MSA, database search, variant calling)
- Low maintenance burden
- High user value

**Action**: Integrate into biometal v1.7.0+

### Archive: Neural Engine ⚠️

**Rationale**:
- Simple quality filtering inappropriate (2,940× slower)
- Infrastructure valuable for future ML tasks
- Better use cases exist (not tested)
- Not part of current core roadmap
- Requires scope expansion discussion

**Action**: Archive research, preserve infrastructure, defer to Phase 3+ strategic planning

---

## Lessons for Future Exploration

1. **Validate Use Cases First**: Ensure operation benefits from acceleration technology
2. **Benchmark Early**: Test prototype before full implementation
3. **Consider Overhead**: Dispatch/encoding costs can dominate simple operations
4. **Compare to Optimized Baseline**: NEON/SIMD can be extremely fast
5. **Infrastructure Has Value**: Even "failed" experiments create reusable frameworks
6. **Document Everything**: Comprehensive docs enable future decision-making

---

## Return to Core Development

**Current Status**: Phase 1, Week 2 complete (50% of consolidation phase)

**Remaining Phase 1 Work**:
- Week 3: Community building (blog post, social media, maintainer outreach)
- Week 4: Quality assurance (property testing, fuzz testing, cross-platform validation)

**Phase 2 Timeline** (Weeks 5-9):
- Rule 3: Parallel BGZF decompression (6.5× speedup)
- Rule 4: Smart mmap (2.5× additional speedup)
- Combined: 16× BAM parsing improvement (55 MiB/s → 895 MiB/s)

**Strategic Pivot Status**: ✅ COMPLETE
**Next Focus**: Return to Phase 1 consolidation (community + quality)

---

**Archive Date**: November 13, 2025
**Archive Location**: planning_archive/strategic_pivot/
**Status**: Strategic pivot exploration complete, returning to core development roadmap
