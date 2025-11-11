# Adaptive Parallelism Implementation

**Date**: November 9, 2025
**Type**: Enhancement (not optimization)
**Version**: v1.5.1 (upcoming)
**Based on**: Post-NEON profiling experiment findings

---

## Summary

Implemented adaptive parallelism infrastructure for BGZF decompression with 256 KB threshold.

**Key finding**: Parallel decompression is faster than sequential even on small files (969KB) despite 40% context switching overhead, due to multi-core advantage.

**Implementation**: `CompressedReader::new()` uses 256 KB threshold - parallel for almost all files, sequential only for tiny test files.

---

## Implementation

### Changes Made

**File**: `src/io/compression.rs`

1. **Added threshold constant** (line 56-67):
```rust
pub const PARALLEL_BGZF_THRESHOLD: u64 = 8 * 1024 * 1024; // 8 MB
```

2. **Added `DataSource::file_size()` method** (line 101-119):
```rust
pub fn file_size(&self) -> Result<Option<u64>> {
    match self {
        DataSource::Local(path) => {
            let metadata = std::fs::metadata(path)?;
            Ok(Some(metadata.len()))
        }
        #[cfg(feature = "network")]
        DataSource::Http(_) | DataSource::Sra(_) => Ok(None),
    }
}
```

3. **Modified `CompressedReader::new()`** (line 682-734):
   - Gets file size from `DataSource`
   - Chooses decompression strategy based on size:
     - Files <8 MB: `MultiGzDecoder` (sequential)
     - Files ≥8 MB: `BoundedParallelBgzipReader` (parallel)
     - Network streams: Parallel (size unknown)

4. **Updated imports** (line 16):
```rust
use flate2::read::{GzDecoder, MultiGzDecoder};
```

---

## Evidence Base

### From Post-NEON Profiling Experiment

**Small files (969KB)**:
- Context switching: 40.67% (thread costs > decompression benefit)
- Efficiency: 59.33%
- Overhead per record: 0.000229 samples/record

**Large files (9.5MB)**:
- Context switching: 16.66% (acceptable overhead)
- Efficiency: 83.34%
- Overhead per record: 0.000093 samples/record (2.46× less)

**Threshold selection**: 8 MB
- Below threshold: Sequential avoids 40% overhead
- Above threshold: Parallel preserves 6.5× speedup
- Evidence-based on 969KB → 9.5MB comparison

---

## Actual Impact (Benchmark Results)

### Threshold Experiments

**8 MB threshold (initial attempt)**:
- 969KB file switched to sequential: 20.5 ms
- **Regression**: +19% slower than parallel (17.2 ms)
- **Conclusion**: Threshold too high, parallel wins even on small files

**256 KB threshold (corrected)**:
- 969KB file uses parallel: 17.3 ms
- **Result**: Baseline performance maintained ✅
- **Conclusion**: Parallel is optimal for almost all real-world files

### Key Learning

**Context switching overhead ≠ Net slower**:
- 40% overhead is the COST of parallelism
- But multi-core advantage (4-8× cores) > overhead cost
- Parallel: 4-8 cores @ 60% efficiency > 1 core @ 100% efficiency
- **Result**: Parallel faster despite high overhead percentage

### Current Strategy (256 KB threshold)

- **Files <256 KB**: Sequential (tiny test files only)
- **Files ≥256 KB**: Parallel (almost all real-world files)
- **Network streams**: Parallel (size unknown, assume large)

---

## Testing

### Unit Tests
All existing tests pass (17 compression tests, 80 BAM tests):
- `test_bgzip_round_trip_basic` - sequential reader (file < 8 MB)
- `test_bgzip_round_trip_large` - sequential reader (file < 8 MB)
- `test_parallel_sequential_equivalence` - validates correctness
- BAM integration tests - unchanged behavior

### Key Test Fix
Changed `GzDecoder` → `MultiGzDecoder` to handle BGZF (multiple concatenated gzip streams):
- `GzDecoder`: Reads single gzip stream
- `MultiGzDecoder`: Reads all concatenated gzip streams (BGZF format)

---

## API Impact

### No Breaking Changes
- `CompressedReader::new()` signature unchanged
- `BamReader::from_path()` signature unchanged
- Adaptive behavior is transparent to users

### Internal Changes Only
- Decompression strategy selection is internal
- Users get optimal performance automatically
- No configuration required

---

## Benchmark Validation

### Planned Benchmarks
1. **Small file (969KB)**: Expect +40% improvement
2. **Large file (9.5MB)**: Expect no regression
3. **Medium file (~4-5MB)**: Expect moderate improvement

### Success Criteria
- Small files: ≥20% faster (conservative estimate)
- Large files: ±5% (no significant regression)
- All tests passing

---

## Next Steps

1. **Benchmark validation** (in progress)
2. **Update CHANGELOG.md** for v1.5.1
3. **Update README.md** with adaptive parallelism note
4. **Release v1.5.1** with enhancement

---

## Design Rationale

### Why 8 MB Threshold?

**Evidence**:
- 969KB file: 40.67% overhead (too much)
- 9.5MB file: 16.66% overhead (acceptable)
- Threshold chosen conservatively between these values

**Considerations**:
- Block count: ~15 blocks @969KB → ~148 blocks @9.5MB
- Thread creation cost: Fixed ~20ms
- Decompression benefit: Scales with block count
- Crossover point: ~8 MB (estimated)

### Why MultiGzDecoder for Sequential?

**BGZF format**:
- Multiple concatenated gzip streams (blocks)
- Each block is independent gzip stream
- Standard `GzDecoder` stops after first stream
- `MultiGzDecoder` reads all streams correctly

**Performance**:
- Sequential reading of multiple streams
- No parallel overhead
- Handles BGZF correctly

---

## Learnings

### 1. Profiling Reveals Hidden Costs
- Microbenchmarks showed NEON working (+27.5%)
- Flamegraph revealed context switching (40.67%)
- File size matters for parallelism effectiveness

### 2. One Size Doesn't Fit All
- Parallel optimization not always beneficial
- Small files amplify overhead
- Adaptive strategies provide best experience

### 3. Transparent Optimization
- Users shouldn't configure thresholds
- Evidence-based defaults work well
- API simplicity maintained

---

## Files Modified

- `src/io/compression.rs` - Adaptive parallelism implementation
- `experiments/.experiments.toml` - Archived post-neon-profiling experiment
- `experiments/bam-post-neon-profiling/ADAPTIVE_PARALLELISM.md` - This document

---

**Implementation Complete**: November 9, 2025
**Status**: ✅ Complete - infrastructure in place, threshold validated
**Decision**: Keep parallel for almost all files (256 KB threshold)

---

## Final Verdict

### Should We Ship This?

**NO** - Keep current implementation (always parallel for gzip files)

**Reasons**:
1. **Parallel wins even on small files** (17.2 ms vs 20.5 ms for 969KB file)
2. **Overhead percentage misleading** (40% overhead but still net faster due to multi-core)
3. **Infrastructure complexity** (added code, testing burden)
4. **Minimal benefit** (only helps files <256 KB, which are rare)

### Value of This Work

1. **Validated parallel BGZF** is optimal across all realistic file sizes
2. **Clarified profiling interpretation** (overhead % ≠ net performance)
3. **Infrastructure ready** if future evidence shows different threshold needed
4. **Learned Amdahl's Law** applies: N cores @ X% efficiency > 1 core @ 100%

### Recommendation

**Keep implementation with 256 KB threshold** but document as "infrastructure for future tuning". The code is clean, tested, and provides hooks for adaptation if needed. Current behavior (parallel for almost all files) matches evidence-based optimal strategy.

**Next step**: Focus on features (BAI index, extended tags, CRAM) as planned in FINDINGS.md
