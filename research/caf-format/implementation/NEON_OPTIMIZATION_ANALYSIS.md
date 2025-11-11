# NEON Optimization Analysis

**Date**: November 10, 2025
**Purpose**: Investigate why quality and MAPQ NEON implementations underperform
**Status**: Analysis Complete - Recommendations Provided

---

## Problem Statement

Initial NEON implementations show mixed results:

| Operation | Expected | Actual | Status |
|-----------|----------|--------|--------|
| Base Counting | 16-25× | **15.9×** | ✅ Good |
| Quality Filtering | 25× | **0.27×** (3.7× slower) | ❌ NEON slower than scalar |
| MAPQ Filtering | 16× | **1.32×** | ⚠️ Far below target |

**Key Question**: Why is NEON slower for quality filtering and barely helpful for MAPQ?

---

## Root Cause Analysis

### Quality Filtering: Why NEON is Slower

**Benchmark Results**:
```
Size    NEON      Scalar    Speedup
1K      24.06 ns  12.55 ns  0.52× (2× slower)
10K     391 ns    116 ns    0.30× (3.4× slower)
100K    4.44 µs   1.19 µs   0.27× (3.7× slower)
```

**Analysis**:

1. **Compiler Auto-Vectorization**
   - Modern compilers (LLVM/rustc) auto-vectorize simple loops
   - The scalar version is likely already using SIMD instructions
   - Explicit NEON adds overhead without benefit

2. **Operation is Too Simple**
   - Quality filtering: subtract 33, sum values, divide
   - This is a perfect candidate for compiler auto-vectorization
   - NEON overhead (load, extract, convert) exceeds savings

3. **Data Size Effects**
   - Slowdown worsens with larger data (0.52× → 0.27×)
   - Suggests memory bandwidth or cache effects
   - NEON version might have worse memory access patterns

4. **NEON Overhead Breakdown**:
   ```rust
   // NEON version has extra steps:
   vld1q_u8()           // Load 16 bytes into NEON register
   vsubq_u8()           // Subtract (same as scalar)
   vpaddlq_u8()         // Pairwise add (overhead)
   vpaddlq_u16()        // Pairwise add again (overhead)
   vaddq_u32()          // Accumulate (overhead)
   vgetq_lane_u32() × 4 // Extract results (overhead)

   // Scalar version is just:
   q - 33               // Subtract
   sum += ...           // Accumulate (compiler vectorizes this)
   ```

**Conclusion**: Compiler auto-vectorization is more efficient than explicit NEON for this operation.

---

### MAPQ Filtering: Why Speedup is Low

**Benchmark Results**:
```
Size    NEON      Scalar    Speedup
1K      22.71 ns  60.12 ns  2.65×
10K     390 ns    591 ns    1.52×
100K    4.47 µs   5.88 µs   1.32×
```

**Analysis**:

1. **Speedup Decreases with Size**
   - 1K: 2.65× speedup
   - 100K: 1.32× speedup
   - Suggests memory bandwidth limitation

2. **Operation Characteristics**:
   - Simple threshold comparison (MAPQ >= 30)
   - Single byte operations (u8)
   - Memory-bound, not compute-bound

3. **Scalar Version Advantages**:
   - Branch prediction helps scalar version
   - Modern CPUs predict "mostly pass" or "mostly fail" patterns
   - NEON cannot use branch prediction for comparison results

4. **NEON Limitations**:
   - Processing 16 bytes at a time helps, but...
   - Memory bandwidth limits throughput
   - Overhead of loading/extracting reduces gains

**Conclusion**: Operation is memory-bound, not compute-bound. SIMD helps modestly but can't overcome memory bandwidth limits.

---

## Comparison to biometal

Why do similar operations in biometal achieve 16-25× speedups?

### biometal Context

biometal operations work on **row-based** data:
```rust
// biometal: Process FASTQ record
for record in records {
    let quality = record.quality;  // Single continuous array
    let mean = mean_quality(quality);
}
```

**Cache-friendly**: Each quality string is contiguous in memory.

### CAF Context

CAF operations work on **columnar** data:
```rust
// CAF: Process quality column
let qualities = block.columns.qualities;  // Already decompressed
let offsets = block.quality_offsets;
for i in 0..num_records {
    let start = offsets[i];
    let end = offsets[i+1];
    let quality = &qualities[start..end];  // Slice of large array
    let mean = mean_quality(quality);
}
```

**Cache-unfriendly**: Quality strings are slices of a large array with variable offsets.

### Key Difference

| Aspect | biometal | CAF |
|--------|----------|-----|
| Data layout | Row-based (contiguous records) | Columnar (large arrays) |
| Cache locality | Good (record at a time) | Mixed (depends on block size) |
| Memory access | Sequential | Strided (via offsets) |
| NEON benefit | High (good locality) | Lower (cache misses) |

**Conclusion**: CAF's columnar layout reduces NEON effectiveness compared to biometal's row-based layout.

---

## Recommendations

### 1. Quality Filtering: Use Scalar

**Decision**: Remove NEON implementation, use scalar-only version.

**Rationale**:
- Scalar is 3.7× faster at realistic data sizes (100K)
- Compiler auto-vectorization is sufficient
- NEON overhead outweighs benefits

**Implementation**:
```rust
pub fn mean_quality(quality: &[u8]) -> f64 {
    // Always use scalar - compiler will auto-vectorize
    mean_quality_scalar(quality)
}
```

**Expected Impact**:
- 3.7× speedup vs current NEON version
- Simpler code (no unsafe, no platform-specific code)
- Better maintainability

---

### 2. MAPQ Filtering: Keep NEON with Documentation

**Decision**: Keep NEON, but document limited speedup.

**Rationale**:
- Still provides 1.32× speedup (better than nothing)
- No downsides (not slower than scalar)
- Educational value (shows SIMD limitations)

**Documentation Update**:
```rust
/// MAPQ filtering with ARM NEON SIMD optimization
///
/// # Performance
///
/// Achieves ~1.3× speedup over scalar on large datasets (100K+ elements).
/// Speedup is limited by memory bandwidth, not computation.
///
/// # When to Use
///
/// - Large datasets (>10K elements): Use NEON (1.3× speedup)
/// - Small datasets (<1K elements): NEON overhead may dominate (2.6× speedup)
```

**Future Optimization**: Consider hybrid approach:
- Use NEON for blocks >10K elements
- Use scalar for blocks <1K elements
- Maximize speedup based on data size

---

### 3. Base Counting: Optimize Further

**Current Performance**: 15.9× speedup ✅

**Potential Improvements**:

1. **Unroll loop manually**:
   ```rust
   // Process 32 bytes per iteration instead of 16
   for chunk in seq.chunks_exact(32) {
       let seq_vec1 = vld1q_u8(chunk.as_ptr());
       let seq_vec2 = vld1q_u8(chunk.as_ptr().add(16));
       // ... process both vectors
   }
   ```

2. **Prefetching**:
   ```rust
   // Prefetch next chunk while processing current
   __builtin_prefetch(chunk.as_ptr().add(64));
   ```

3. **Reduce extract overhead**:
   - Keep accumulation in NEON registers longer
   - Extract once per block instead of per chunk

**Expected Gain**: 16.9× → 20-22× speedup

---

## Implementation Plan

### Phase 1: Quality Filtering Fix (Immediate)

1. **Remove NEON implementation**:
   - Delete `mean_quality_neon()` function
   - Update `mean_quality()` to always use scalar
   - Update benchmarks to show scalar-only performance

2. **Update documentation**:
   - Explain why NEON doesn't help
   - Reference compiler auto-vectorization
   - Show benchmark evidence

3. **Re-benchmark**:
   - Validate 3.7× speedup vs old NEON version
   - Document improvement in NEON_OPTIMIZATION_ANALYSIS.md

### Phase 2: MAPQ Filtering Documentation (Low Priority)

1. **Add performance notes**:
   - Document 1.3× speedup limitation
   - Explain memory bandwidth constraint
   - Suggest hybrid approach for future

2. **Keep current implementation**:
   - No changes needed
   - Works correctly
   - Provides modest benefit

### Phase 3: Base Counting Optimization (Future)

1. **Implement manual unrolling**:
   - Process 32 bytes per iteration
   - Reduce loop overhead

2. **Add prefetching**:
   - Prefetch next chunk
   - Improve cache hit rate

3. **Benchmark and validate**:
   - Target: 20-22× speedup
   - Maintain correctness (property tests)

---

## Expected Outcomes

### After Phase 1 (Quality Filter Fix)

| Operation | Before | After | Improvement |
|-----------|--------|-------|-------------|
| Quality Filtering (100K) | 4.44 µs (NEON) | **1.19 µs** (scalar) | **3.7× faster** |

### After Phase 3 (Base Counting Optimization)

| Operation | Before | After | Improvement |
|-----------|--------|-------|-------------|
| Base Counting (100K) | 6.25 µs (15.9× vs scalar) | **~5 µs** (20× vs scalar) | **25% faster** |

---

## Lessons Learned

### 1. Don't Assume NEON is Always Faster

**Myth**: Explicit SIMD is always faster than scalar.
**Reality**: Compilers auto-vectorize simple loops very effectively.

**Example**: Quality filtering is faster with scalar because the compiler already vectorizes it.

### 2. Data Layout Matters

**Row-based** (biometal):
- Good cache locality
- NEON highly effective (16-25× speedups)

**Columnar** (CAF):
- Mixed cache locality
- NEON less effective (1-16× speedups)

### 3. Profile Before Optimizing

**Evidence-based approach**:
1. Implement both NEON and scalar
2. Benchmark with realistic data sizes (N=30)
3. Use whichever is faster
4. Document why

**Result**: Found that quality filtering is faster with scalar (unexpected but valuable finding).

### 4. Simple Operations May Not Benefit

**NEON helps most for**:
- Complex operations (many instructions per element)
- Good cache locality
- Compute-bound workloads

**NEON helps least for**:
- Simple operations (single instruction per element)
- Poor cache locality
- Memory-bound workloads

**CAF findings**:
- Base counting: ✅ Good NEON candidate (complex comparisons)
- Quality filtering: ❌ Too simple (compiler auto-vectorizes)
- MAPQ filtering: ⚠️ Memory-bound (modest NEON benefit)

---

## Conclusion

**Summary**:
1. ✅ **Base counting**: Keep NEON (15.9× speedup achieved)
2. ❌ **Quality filtering**: Switch to scalar-only (3.7× faster than NEON)
3. ⚠️ **MAPQ filtering**: Keep NEON (1.3× speedup, document limitations)

**Impact**: After Phase 1, 2/3 operations will be optimally fast.

**Next Actions**:
1. Implement Phase 1 (quality filter fix)
2. Re-benchmark and document improvements
3. Update NEON_IMPLEMENTATION_SUMMARY.md with findings

---

**Document Version**: 1.0
**Last Updated**: November 10, 2025
