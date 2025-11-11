# PROPOSAL: Post-NEON Profiling and Next Bottleneck Identification

**Experiment**: Identify next optimization target after NEON sequence decoding (v1.5.0)
**Status**: PROPOSAL
**Date**: November 9, 2025
**Expected Duration**: 2-4 hours
**Researcher**: Claude (AI assistant) with user guidance

---

## Context

We just completed v1.5.0 with ARM NEON sequence decoding:
- **Sequence decoding**: 4.62× faster (30.2% → ~8% CPU time)
- **Overall BAM parsing**: +27.5% faster (43.0 → 55.1 MiB/s)

As predicted in `experiments/bam-simd-sequence-decoding/FINDINGS.md`:
> "Expect Cascading Bottlenecks: Each optimization shifts the bottleneck distribution"

We need to re-profile to identify the current bottleneck.

---

## Hypothesis

**Expected bottleneck distribution** (post-NEON, estimated):
- BGZF decompression: ~30-35% (parallel implementation exists, but needs validation)
- Record parsing: ~15-20% (CIGAR, flags, tags)
- Sequence decoding: ~8% (optimized with NEON)
- Quality scores: ~8-10% (memcpy, no optimization needed)
- Other: ~20-30% (I/O, memory allocation, etc.)

**Primary candidates for next optimization** (Rule 1: ≥15% CPU time):
1. BGZF decompression (~30-35%) - IF parallel isn't working effectively
2. Record parsing (~15-20%) - IF it reaches threshold

**Validation needed**:
- Is parallel BGZF actually providing 6.5× speedup?
- Small test file (969KB, ~1 block) may not show parallel benefit
- Need larger test file (≥8 MB, ≥8 blocks) to validate

---

## Objectives

1. **Profile current BAM parsing** (with NEON enabled)
   - Create microbenchmarks for each operation
   - Measure CPU time percentage for each component
   - Use flamegraph for holistic view

2. **Validate parallel BGZF performance**
   - Test with small file (969KB, ~1 block)
   - Test with larger file (≥8 MB, ≥8 blocks)
   - Verify 6.5× speedup from Entry 029

3. **Identify next optimization target**
   - Apply Rule 1 threshold (≥15% CPU time)
   - Estimate expected speedup
   - Calculate overall impact

4. **Create proposal for next optimization**
   - Detailed implementation plan
   - Go/no-go criteria
   - Expected timeline

---

## Methodology

### Phase 1: Microbenchmark Suite (1 hour)

Create isolated benchmarks for each BAM component:

```rust
// Measure ONLY BGZF decompression
bench_bgzf_only(file) -> Duration

// Measure ONLY record parsing (pre-decoded BGZF blocks)
bench_record_parsing(decompressed_blocks) -> Duration

// Measure ONLY sequence decoding (with NEON)
bench_sequence_neon(packed_sequences) -> Duration

// Measure ONLY quality score extraction
bench_quality_extraction(quality_data) -> Duration

// Measure ONLY CIGAR parsing
bench_cigar_parsing(cigar_data) -> Duration

// Measure ONLY tag parsing
bench_tag_parsing(tag_data) -> Duration
```

**Sample size**: N=30 for statistical significance

### Phase 2: Profiling (1 hour)

1. **Run microbenchmarks**:
   - Measure each component in isolation
   - Calculate CPU time percentage
   - Identify components ≥15% threshold

2. **Generate flamegraph**:
   - `cargo flamegraph --bench bam_parsing`
   - Visual validation of CPU time distribution
   - Identify unexpected bottlenecks

3. **Parallel BGZF validation**:
   - Small file: 969KB test file (current)
   - Large file: Generate 8+ MB test BAM (8+ blocks)
   - Compare parallel vs sequential decompression

### Phase 3: Analysis (1 hour)

1. **Calculate bottleneck distribution**:
   - Sort components by CPU time percentage
   - Identify primary target (highest ≥15%)
   - Estimate optimization potential

2. **Evaluate parallel BGZF**:
   - If small file: ~1× speedup (expected, 1 block)
   - If large file: ~6.5× speedup (expected, 8+ blocks)
   - Decision: Is parallel BGZF working correctly?

3. **Determine next optimization**:
   - Primary: Highest CPU time ≥15%
   - Expected speedup: Based on operation type
   - Overall impact: Calculate using Amdahl's law

### Phase 4: Proposal (1 hour)

Create detailed proposal for next optimization:
- Clear objectives and success criteria
- Implementation approach
- Expected timeline
- Go/no-go criteria

---

## Success Criteria

### Profiling Validation
- ✅ Microbenchmarks isolate each component
- ✅ CPU time percentages sum to ~100%
- ✅ Flamegraph confirms microbenchmark results
- ✅ Statistical significance (N=30, 95% CI)

### Parallel BGZF Validation
- ✅ Small file (1 block): ~1× speedup (no parallel benefit)
- ✅ Large file (8+ blocks): ~6.5× speedup (parallel working)
- ✅ Thread utilization: 8 cores active during decompression

### Next Target Identification
- ✅ Identified component with ≥15% CPU time
- ✅ Expected speedup calculated (based on operation type)
- ✅ Overall impact ≥5% improvement

---

## Go/No-Go Criteria

**GO if**:
1. Profiling identifies component with ≥15% CPU time
2. Expected overall improvement ≥5%
3. Implementation complexity reasonable (<2 weeks)
4. Clear optimization strategy exists

**NO-GO if**:
1. No component ≥15% CPU time (diminishing returns)
2. Expected overall improvement <3%
3. Implementation complexity excessive (>1 month)
4. Unclear optimization strategy

---

## Expected Outcomes

### Scenario 1: Parallel BGZF Not Working
- **Finding**: BGZF still ~30-35% CPU time, no parallel speedup
- **Root cause**: Implementation issue or test file too small
- **Next step**: Fix parallel BGZF or generate larger test files

### Scenario 2: Record Parsing is Bottleneck
- **Finding**: Record parsing ~15-20% CPU time
- **Target**: CIGAR parsing, tag parsing, or field extraction
- **Strategy**: Optimize parsing logic, consider SIMD for array operations

### Scenario 3: Multiple Small Bottlenecks
- **Finding**: No single component ≥15% CPU time
- **Implication**: Diminishing returns, overall architecture good
- **Next step**: Focus on other features (BAI index, extended tags)

### Scenario 4: I/O Bound
- **Finding**: BGZF + parsing <50% CPU time, rest is I/O
- **Implication**: Performance limited by disk/memory bandwidth
- **Next step**: Optimize I/O patterns (madvise, prefetching)

---

## Timeline

| Phase | Duration | Activities |
|-------|----------|------------|
| Phase 1 | 1 hour | Create microbenchmark suite |
| Phase 2 | 1 hour | Run profiling and parallel BGZF validation |
| Phase 3 | 1 hour | Analyze results, identify next target |
| Phase 4 | 1 hour | Create detailed proposal for next optimization |
| **Total** | **4 hours** | **Complete profiling and next step proposal** |

---

## Resources Required

1. **Benchmarking**:
   - Criterion benchmarks for each component
   - Flamegraph for visual profiling

2. **Test Data**:
   - Small file: 969KB (existing `tests/data/test.bam`)
   - Large file: 8+ MB BAM (generate with samtools)

3. **Tools**:
   - `cargo bench` for microbenchmarks
   - `cargo flamegraph` for CPU profiling
   - `htop` for thread utilization monitoring

---

## Risk Assessment

### Low Risk
- Profiling is non-invasive (no code changes)
- Microbenchmarks validate existing implementation
- Clear methodology from previous experiments

### Medium Risk
- Generating large test files may take time
- Parallel BGZF validation requires multi-core machine
- Flamegraph interpretation requires expertise

### Mitigation
- Use existing test data where possible
- Document parallel BGZF expectations clearly
- Cross-validate flamegraph with microbenchmarks

---

## References

### Evidence Base
- **Rule 1** (OPTIMIZATION_RULES.md): ≥15% CPU time threshold for SIMD
- **Entry 029**: Parallel BGZF (6.5× expected speedup)
- **Sequence NEON**: Validated 30.2% → ~8% CPU time reduction

### Previous Experiments
- `experiments/bam-simd-sequence-decoding/FINDINGS.md`: Cascading bottlenecks
- `experiments/native-bam-implementation/DECISION_REVISED.md`: Phase 0 profiling

---

## Approval

**Researcher**: Ready to proceed with Phase 1 (microbenchmark suite creation)
**User**: [Awaiting approval]

---

**Status**: PROPOSAL (awaiting user approval to proceed)
**Next Step**: Create microbenchmark suite for each BAM component
