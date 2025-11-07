# Research Log: SIMD Minimizers Analysis

**Experiment**: simd-minimizers-analysis
**Start Date**: November 6, 2025
**Researcher**: Scott Handley + Claude Code

---

## Day 1: Initial Analysis (November 6, 2025)

### Session 1: Repository Setup & First Impressions

**Time**: 20:30 PST

**Actions Taken**:
1. ✅ Created experiment directory following TEMPLATE structure
2. ✅ Written comprehensive PROPOSAL.md with clear go/no-go criteria
3. ✅ Cloned `rust-seq/simd-minimizers` repository
4. ⏳ Starting code analysis

**Initial Repository Scan**:
```bash
# Repository structure
simd-minimizers/
├── src/
│   ├── lib.rs              # Core library
│   ├── minimizer.rs        # Minimizer implementations?
│   └── ...
├── benches/
│   └── ...                 # Performance benchmarks
├── examples/
│   └── ...                 # Usage examples
└── Cargo.toml
```

**First Questions**:
1. What SIMD backends do they support? (AVX2, NEON, both?)
2. Where is the rolling hash implementation?
3. How do they handle the argmin (find minimum) operation?
4. What's the memory footprint?

**Next Steps**:
- [ ] Read `src/lib.rs` - understand API surface
- [ ] Find SIMD-specific code (look for `#[cfg(target_arch)]`)
- [ ] Locate rolling hash implementation
- [ ] Study benchmark code to understand testing methodology

---

### Session 2: Source Code Analysis

**Status**: Complete ✅

**Time**: 20:45 PST

#### File: src/lib.rs (506 lines)

**Key observations**:
```rust
//! ## Minimizers
//!
//! Minimizers are found as follows:
//! 1. Split the input to 8 chunks that are processed in parallel using SIMD.
//! 2. Compute a 32-bit ntHash rolling hash of the k-mers.
//! 3. Use the 'two stacks' sliding window minimum on the top 16 bits of each hash.
//! 4. Break ties towards the leftmost position by storing the position in the bottom 16 bits.
//! 5. Compute 8 consecutive minimizer positions, and dedup them.
//! 6. Collect the deduplicated minimizer positions from all 8 chunks into a single vector.
```

**API Design**:
- Builder pattern: `minimizers(k, w).run(seq, &mut out_vec)`
- Reusable output vector (avoids allocations)
- Supports custom hashers via `hasher(&hasher)` method

#### Key Findings

**Finding 1: SIMD Architecture Support** ✅
- **CONFIRMED**: Supports both **AVX2 AND NEON**!
- Uses `packed-seq` crate for cross-platform SIMD abstraction
- Uses `wide` crate for portable SIMD types (`u32x8`, `i32x8`)
- `.cargo/config.toml`: `rustflags = ["-C", "target-cpu=native"]`
- **Evidence**: Found `bench/results-neon.json` - they've benchmarked on ARM NEON!
- **No `#[cfg(target_arch)]` blocks**: Portable SIMD via abstraction layer

**Finding 2: Algorithm Structure** ✅
Data flow from sequence → minimizers:

```
Input sequence
    ↓
1. Split to 8 SIMD lanes (packed-seq::Seq::iter_bp)
    ↓
2. Rolling hash (ntHash from seq-hash crate)
   - 32-bit hash per k-mer
   - Top 16 bits: hash value
   - Bottom 16 bits: position (for tie-breaking)
    ↓
3. Sliding window minimum (two stacks algorithm)
   - Finds minimum in each window of w k-mers
   - O(1) amortized per window
   - SIMD: process 8 windows in parallel
    ↓
4. Deduplication (collect::append_unique_vals)
   - SIMD-based dedup using intrinsics
    ↓
5. Collect into output vector
   - All 8 SIMD lanes merged
```

**Key files**:
- `src/lib.rs`: API, Builder pattern (506 lines)
- `src/minimizers.rs`: Core algorithm (200 lines)
- `src/sliding_min.rs`: Two stacks sliding minimum (300+ lines)
- `src/canonical.rs`: Canonical minimizers (62 lines)
- `src/collect.rs`: Deduplication logic (370 lines)

**Finding 3: Rolling Hash** ✅
- **ntHash** from `seq-hash` crate (external dependency)
- **NOT in this repository** - abstracted via `KmerHasher` trait
- Key property: `hash_kmers_simd()` method → processes 8 lanes
- Also supports: `MulHasher`, `AntiLexHasher` (for general text)

**ntHash characteristics** (from documentation):
```rust
//! By default, the library uses the `ntHash` hash function, which maps
//! each DNA base `ACTG` to a pseudo-random value using a table lookup.
//! This hash function is specifically designed to be fast for hashing
//! DNA sequences with input type [`packed_seq::PackedSeq`].
```

**Finding 4: Memory Usage** ✅
- **Reuses output vector**: `out_vec: &mut Vec<u32>` passed in, appended to
- **Thread-local cache**: `thread_local! { static CACHE: ... }`
  - Caches ring buffers for sliding window algorithm
  - Avoids per-call allocations
- **Memory pattern**: Similar to our approach (reusable buffers)
- **BUT**: Appears to buffer all minimizers before returning
  - Not constant-memory streaming like biometal
  - Processes entire sequence at once

---

#### Critical Comparison: SimdMinimizers vs biometal Entry 034

**What They Do Differently**:

1. **ntHash vs FNV-1a**:
   ```rust
   // Our Entry 034: FNV-1a (sequential, not vectorizable)
   let mut hash = FNV_OFFSET;
   for &byte in kmer {
       hash ^= byte as u64;
       hash = hash.wrapping_mul(FNV_PRIME);
   }

   // Their approach: ntHash (rolling, vectorizable)
   // Uses seq-hash crate with hash_kmers_simd() method
   // Processes 8 k-mers in parallel via SIMD
   ```

2. **Two Stacks Sliding Minimum**:
   ```rust
   // Our Entry 034: Linear scan per window
   let mut min_hash = u64::MAX;
   for kmer in window {
       let hash = fnv1a(kmer);
       if hash < min_hash {
           min_hash = hash;
       }
   }

   // Their approach: O(1) amortized two stacks algorithm
   // Maintains prefix/suffix minimums
   // Clever trick: packs hash (16 bits) + position (16 bits) into u32
   ```

3. **SIMD Parallelism**:
   ```
   Our Entry 034: Single-threaded, one window at a time

   Their approach: 8 SIMD lanes processing 8 windows simultaneously
   ```

**Why They Get 9.5× and We Got 1.26×**:

1. ✅ **Rolling hash** (ntHash) is vectorizable → FNV-1a is not
2. ✅ **Two stacks algorithm** reduces window scan from O(w) to O(1) amortized
3. ✅ **8-way SIMD parallelism** via packed-seq lane splitting
4. ✅ **Decoupled operations**: Hash all k-mers first, then find minimums
5. ⚠️ **Trade-off**: Not streaming (buffers entire sequence)

**Key Insight**: They've successfully **decoupled** the SIMD-friendly parts (hashing, argmin) from the data-structure-bound parts (deduplication). This is the algorithmic breakthrough!

---

### Preliminary Hypotheses (UPDATED)

**Hypothesis 1: Decoupled Hash Computation** ✅ CONFIRMED
They compute all hashes first (SIMD-friendly batch), then find minimizers (SIMD argmin), separating these from HashMap operations.

```
Our approach (Entry 034):
for window in windows {
    for kmer in window { hash + compare }  // Coupled, FNV-1a
    hashmap.insert()                        // Data-structure-bound
}
Result: 1.26× (HashMap dominates)

Their approach (CONFIRMED from code):
all_hashes = hash_kmers_simd(sequence)    // ntHash, 8 lanes SIMD
for window in all_hashes {
    min = sliding_min_simd(window)        // Two stacks, 8 lanes SIMD
}
dedup(mins)                               // HashMap isolated
Result: 9.5× (SIMD-friendly parts optimized)
```

**Hypothesis 2: Rolling Hash Formula** ✅ CONFIRMED (ntHash)
ntHash uses rolling hash with table lookups:
- Each base (ACTG) → pseudo-random constant (table)
- Rolling update allows incremental computation
- **Vectorizable** because it's arithmetic + table lookup

FNV-1a (our Entry 034):
- Sequential XOR and multiply
- **Not vectorizable** due to data dependency chain

**Hypothesis 3: SIMD Minimum Finding** ✅ CONFIRMED (Two Stacks)
They use **two stacks sliding minimum** algorithm:
- Maintains prefix and suffix minimums
- O(1) amortized per window (vs O(w) naive scan)
- Processes 8 windows in parallel via SIMD (`u32x8`)
- Clever: packs hash value (16 bits) + position (16 bits) into u32

---

### Questions Arising (UPDATED)

1. **Q**: Do they use NEON or just AVX2?
   **A**: ✅ **BOTH!** Via `packed-seq` + `wide` crates (portable SIMD). Found `bench/results-neon.json`.

2. **Q**: What's the k/w parameter range they support?
   **A**: From benchmarks: k=19-31, w=5-19. Appears flexible (k stored in 16 bits of u32).

3. **Q**: How does memory scale with sequence length?
   **A**: ⚠️ **Linear O(n)** - buffers all hashes, not constant-memory streaming.
   **Trade-off**: Speed vs memory.

4. **Q**: Do they maintain streaming property (constant memory)?
   **A**: ❌ **NO** - processes entire sequence, buffers results.
   **Impact**: Won't work for our 5TB dataset streaming use case.

5. **NEW Q**: Can we adapt their technique to streaming architecture?
   **Status**: Key question for Day 2-3 analysis.
   **Idea**: Block-based processing (Rule 2) with ntHash + two stacks?

6. **NEW Q**: Is ntHash hash quality good enough for minimizer indexing?
   **Status**: Check their paper/tests for collision rates vs FNV-1a.

---

### Comparison to Entry 034

**Our Findings (Entry 034)**:
- Minimizers NEON: 1.02-1.26× speedup
- Conclusion: Data-structure-bound (HashMap 50-60% of runtime)
- Implementation: FNV-1a hash per k-mer, immediate HashMap insert

**Their Claims**:
- 9.5× speedup for w=5
- 4.5× speedup for w=19
- Human genome: 4.1 seconds

**Key Difference to Investigate**:
Why 9.5× vs our 1.26×?

Possible reasons:
1. ✅ Different hash function (rolling vs FNV-1a)
2. ✅ Batch processing (all hashes at once)
3. ✅ SIMD argmin (efficient minimum finding)
4. ❓ Different use case (canonical vs forward minimizers?)
5. ❓ Different hardware (AVX2 wider than NEON?)

---

#### Day 1 Session 2 Summary

**Major Discoveries**:

1. ✅ **NEON support confirmed** - Portable SIMD via `packed-seq` + `wide`
2. ✅ **ntHash rolling hash** - Vectorizable (vs our non-vectorizable FNV-1a)
3. ✅ **Two stacks algorithm** - O(1) amortized sliding minimum (vs O(w) scan)
4. ✅ **8-way SIMD parallelism** - Process 8 windows simultaneously
5. ⚠️ **Not streaming** - Buffers entire sequence (trade-off)

**Why They Achieve 9.5× Speedup**:
- Rolling hash (ntHash) is vectorizable
- Two stacks reduces algorithmic complexity
- 8-lane SIMD parallelism
- Decoupled SIMD-friendly operations from data structures

**Implications for biometal**:
- ✅ Technique IS applicable to ARM NEON (not AVX2-only)
- ⚠️ Memory trade-off: O(n) vs our O(1) streaming
- 🔬 Research question: Can we adapt to block-based streaming?

**Files Analyzed**:
- `src/lib.rs` (506 lines) - API, Builder pattern
- `src/minimizers.rs` (200 lines) - Core algorithm
- `src/sliding_min.rs` (300+ lines) - Two stacks implementation
- `src/canonical.rs` (62 lines) - Canonical minimizers
- `.cargo/config.toml` - Compile flags (target-cpu=native)
- `bench/results-neon.json` - NEON benchmark data EXISTS!

**Next Steps** (Day 2):
- [ ] Examine `seq-hash` crate (ntHash implementation)
- [ ] Understand two stacks algorithm in detail
- [ ] Build and run benchmarks on Mac M-series
- [ ] Compare NEON performance to their published results

---

### Next Session Plan

**Day 2 (Tomorrow)**:
1. ✅ Source code analysis (COMPLETE)
2. Study ntHash implementation (seq-hash crate)
3. Understand two stacks algorithm mechanics
4. Build and run benchmarks on M-series Mac
5. Compare performance to Entry 034

**Goals for Day 2**:
- [ ] Deep dive into ntHash rolling hash formula
- [ ] Understand two stacks algorithm step-by-step
- [ ] Run their benchmarks (cargo bench)
- [ ] Compare results to published 9.5× claim
- [ ] Document algorithm with diagrams

---

### Notes & Observations

**Note 1**: This is their first version (Jan 2025) - very recent work
**Note 2**: Published in SEA 2025 (Symposium on Experimental Algorithms) - peer-reviewed
**Note 3**: Authors from ETH Zurich - strong algorithmic background
**Note 4**: Code is MIT licensed - we can study freely
**Note 5**: Uses portable SIMD (`packed-seq`, `wide`) - not architecture-specific intrinsics

**Key Insight**: The fact that they got 9.5× suggests there IS a SIMD-friendly way to do minimizers, we just haven't found it yet. This is exciting - means improvement is possible!

**Critical Discovery**: We now understand WHY they succeeded where we didn't:
1. **Hash function choice matters**: ntHash (vectorizable) vs FNV-1a (sequential)
2. **Algorithm matters**: Two stacks (O(1)) vs naive scan (O(w))
3. **Decoupling matters**: SIMD-friendly ops separated from data structures

**Trade-off Identified**: Their approach sacrifices constant-memory streaming for speed. Question: Can we combine their technique with our streaming architecture?

---

## Day 2: [To be continued...]

**Status**: Not started
**Planned**: November 7, 2025

---

## Decision Log

| Date | Decision | Rationale |
|------|----------|-----------|
| Nov 6 | Start experiment | SimdMinimizers' 9.5× speedup contradicts our Entry 034 (1.26×) |
| TBD | GO/NO-GO | Based on Week 1 findings |

---

## References Used

**Papers**:
- [SimdMinimizers bioRxiv](https://www.biorxiv.org/content/10.1101/2025.01.27.634998v1)
- Entry 034: K-mer operations analysis

**Code**:
- [rust-seq/simd-minimizers](https://github.com/rust-seq/simd-minimizers)
- biometal: `src/operations/kmer.rs`

**Blog Posts**:
- [CuriousCoding: SIMD Minimizers](https://curiouscoding.nl/posts/simd-minimizers/)

---

**Log maintained by**: Claude Code assisting Scott Handley
**Format**: Daily updates with session-level granularity
**Purpose**: Evidence-based decision making for biometal development
