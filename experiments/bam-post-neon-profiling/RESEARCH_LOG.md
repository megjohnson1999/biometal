# RESEARCH LOG: Post-NEON Profiling

**Experiment**: Identify next optimization target after NEON sequence decoding
**Date Started**: November 9, 2025
**Researcher**: Claude (AI assistant) with user guidance

---

## 2025-11-09 17:30 - Phase 2 Complete: Flamegraph Analysis

**Objective**: Generate flamegraph to visually identify bottleneck distribution

### Flamegraph Results

**Generated**: `experiments/bam-post-neon-profiling/flamegraph.svg`
- Total samples: 563
- Workload: 10 iterations × 100K records = 1M records
- Test file: 969KB (small, ~15 blocks)

### Key Findings

**1. NEON Sequence Decoding Validated** ✅
- Sequence decoding: **1.2% CPU time** (down from 30.2% pre-NEON)
- NEON optimization working perfectly
- No longer a bottleneck

**2. Context Switching Dominates** ⚠️
- **swtch_pri**: 40.67% (kernel context switching)
- **Rayon overhead**: 16.2% (thread management)
- **Total parallelization overhead**: 56.8%
- **Actual decompression work**: 10.7%

**3. Biometal Components Well-Optimized** ✅
- Total biometal overhead: **7.1%**
- Record parsing: 2.1%
- Sequence decoding: 1.2%
- CIGAR parsing: 1.1%
- **ALL below 15% threshold**

### Analysis

**Why is parallel BGZF showing high overhead?**

Test file characteristics:
- **Size**: 969KB (~15 blocks)
- **Block size**: ~64KB each
- **Decompression time per block**: ~1-2ms
- **Thread creation/switching cost**: ~20ms per block

**Overhead ratio**: **5.3× overhead per block**

**Expected behavior**: Entry 029 validated parallel BGZF on **8+ MB files** (100+ blocks)
- Small files amplify thread overhead
- Large files amortize overhead across many blocks

### Decision: Proceed with Option A

**Validate parallel BGZF with larger file**:
1. Generate 8+ MB test BAM (100+ blocks)
2. Expect: 30-35% decompression, <10% overhead
3. Confirm 6.5× speedup from Entry 029

---

## 2025-11-09 17:45 - Phase 2B Start: Large File Validation

**Objective**: Generate 8+ MB BAM file to validate parallel BGZF on realistic workload

**Strategy**: Replicate existing test file multiple times to create larger dataset

**Target**:
- File size: ≥8 MB
- Estimated blocks: 100+ (8MB / 64KB = ~125 blocks)
- Records: ~800K (100K × 8 replications)

**Expected results**:
- Decompression: 30-35% CPU time
- Parallelization overhead: <10%
- Context switching: <15%
- Overall speedup vs sequential: ~6.5× (Entry 029)

Starting file generation...
