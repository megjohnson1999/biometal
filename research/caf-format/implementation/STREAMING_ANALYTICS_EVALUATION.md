# CAF Streaming Analytics Evaluation Plan

**Date**: November 11, 2025
**Purpose**: Evaluate CAF's real value proposition - streaming analytics at scale
**Context**: CAF was designed for biometal's streaming architecture, not as a BAM replacement

---

## Hypothesis

**CAF should outperform BAM for streaming analytics** where you only need a subset of columns, because:
1. Column-selective reading (decompress only what you need)
2. NEON optimization for columnar operations (16× for base counting)
3. Constant memory footprint (~5 MB, validated)
4. Network efficiency (fewer bytes transferred)

---

## Critical Benchmarks

### 1. Quality-Only Filtering

**Scenario**: Filter 1M reads by mean quality ≥30, discard rest

**BAM approach**:
```rust
for record in bam_reader {
    // Decompresses: sequence, quality, CIGAR, tags, names
    let quality = record.quality();
    if mean_quality(quality) >= 30 {
        // Only 10-20% of records pass
    }
}
```

**CAF approach**:
```rust
for block in caf_reader {
    // Decompress ONLY quality column
    let qualities = block.quality_column();
    let passing = filter_by_quality(qualities, 30);
    // Never touched sequences/CIGAR/tags!
}
```

**Metrics**:
- Time to filter 1M records
- Bytes decompressed
- Memory footprint
- **Expected**: CAF 3-5× faster (only decompress quality)

---

### 2. MAPQ-Only Filtering

**Scenario**: Filter 1M reads by MAPQ ≥30

**BAM approach**:
- Decompress full records
- Extract MAPQ byte
- Filter

**CAF approach**:
- Decompress ONLY MAPQ column (single byte per record)
- Filter with NEON (1.3× speedup)

**Metrics**:
- Time to filter 1M records
- Bytes decompressed
- **Expected**: CAF 5-10× faster (MAPQ is 1 byte, sequences are 100s of bytes)

---

### 3. Base Composition Analysis

**Scenario**: Count A/C/G/T in 1M reads, no other data needed

**BAM approach**:
- Decompress full records
- Decode 4-bit sequences
- Count bases (scalar)

**CAF approach**:
- Decompress ONLY sequence column
- Count bases with NEON (16× speedup)

**Metrics**:
- Time to count bases in 1M records
- Bytes decompressed
- **Expected**: CAF 10-16× faster (column-only + NEON)

---

### 4. Combined Filtering (Quality + MAPQ)

**Scenario**: Filter by quality ≥30 AND MAPQ ≥30

**BAM approach**:
- Decompress full records
- Check both filters

**CAF approach**:
- Decompress quality + MAPQ columns only
- Filter (avoid sequences/CIGAR/tags)

**Metrics**:
- Time to filter 1M records
- Bytes decompressed
- **Expected**: CAF 2-3× faster (two columns vs full record)

---

### 5. Network Streaming Efficiency

**Scenario**: Stream from S3, filter by quality, save passing records

**BAM approach**:
```
Download: 100 MB (full BAM)
Filter: Decompress all → 10% pass
Result: 10 MB output
Total bytes transferred: 100 MB
```

**CAF approach**:
```
Download block headers: 1 MB
Download quality columns: 10 MB (dictionary compressed)
Filter: 10% pass
Download other columns for passing records: 9 MB
Result: 10 MB output
Total bytes transferred: 20 MB
```

**Metrics**:
- Bytes downloaded from network
- Time to complete workflow
- **Expected**: CAF 5× fewer bytes transferred

---

## What Success Looks Like

### Strong Evidence for CAF

1. **Quality filtering**: 3-5× faster than BAM (column-only vs full record)
2. **MAPQ filtering**: 5-10× faster than BAM (tiny column vs full decode)
3. **Base counting**: 10-16× faster than BAM (column-only + NEON)
4. **Network efficiency**: 3-5× fewer bytes transferred
5. **Memory**: Constant ~5 MB regardless of dataset size

### Weak Evidence / CAF Not Better

1. Quality filtering only 1.2× faster (column overhead dominates)
2. Network bytes similar (CAF 1.6× larger negates column savings)
3. Memory not constant (streaming implementation flawed)
4. End-to-end no better than BAM with samtools

---

## Implementation Plan

### Phase 1: Basic Column Benchmarks (2 hours)

**Create**: `benches/streaming_analytics.rs`

Benchmarks:
1. `quality_only_decode` (BAM vs CAF)
2. `mapq_only_decode` (BAM vs CAF)
3. `sequence_only_decode` (BAM vs CAF)
4. `full_record_decode` (baseline)

**Measure**: Time and bytes decompressed

### Phase 2: End-to-End Workflows (3 hours)

**Create**: `examples/streaming_workflows.rs`

Workflows:
1. Quality filter (input → filter → output)
2. Base composition (input → count → stats)
3. MAPQ filter (input → filter → output)
4. Combined filter (quality + MAPQ)

**Measure**: Total time, memory, bytes processed

### Phase 3: Network Simulation (2 hours)

**Create**: Mock S3 streaming benchmark

Simulate:
1. Network latency (10-50 ms per request)
2. Bandwidth limits (10-100 MB/s)
3. Block-level fetching

**Measure**: Total bytes transferred, time to complete

### Phase 4: Analysis & Report (2 hours)

**Create**: `STREAMING_ANALYTICS_RESULTS.md`

Document:
1. Benchmark results with statistical analysis
2. Use cases where CAF wins
3. Use cases where BAM is fine
4. Break-even analysis (when does CAF pay off?)
5. Recommendations for adoption

---

## Expected Findings

### Where CAF Should Win Big

1. **Selective column analytics** (quality-only, MAPQ-only)
   - Hypothesis: 3-10× faster due to avoiding full decode

2. **Network streaming with filtering**
   - Hypothesis: 3-5× fewer bytes transferred

3. **Base composition with NEON**
   - Hypothesis: 10-16× faster (column + SIMD)

### Where CAF May Not Help

1. **Full record processing** (need all columns anyway)
   - No advantage over BAM

2. **Small datasets** (< 1M records)
   - Overhead may dominate

3. **Random access** (specific genomic regions)
   - BAM+BAI is mature, CAF would need indexing

---

## Use Case Matrix

| Use Case | CAF Advantage | Reason |
|----------|---------------|--------|
| **Quality filtering** | High (3-5×) | Column-only decode |
| **MAPQ filtering** | Very High (5-10×) | Tiny column vs full record |
| **Base composition** | Very High (10-16×) | Column-only + NEON |
| **Sequence stats** | High (5-10×) | Column-only + NEON |
| **Network streaming + filter** | High (3-5×) | Fewer bytes transferred |
| **Full variant calling** | None | Need all columns |
| **Random access (IGV)** | None | Would need indexing |
| **Long-term archival** | Negative | Files 1.6× larger |

---

## Decision Criteria

### CAF is Worth It If:

1. Quality filtering is ≥3× faster than BAM
2. Network streaming transfers ≥3× fewer bytes
3. Memory stays constant (~5 MB) at scale
4. NEON provides genuine speedups (10×+) for analytics

### CAF is NOT Worth It If:

1. Column overhead makes it <2× faster than BAM
2. File size (1.6×) negates network savings
3. Memory scales with dataset size (streaming broken)
4. Use cases require full records anyway

---

## Timeline

**Total**: 9-11 hours

- Phase 1 (Basic benchmarks): 2 hours
- Phase 2 (Workflows): 3 hours
- Phase 3 (Network sim): 2 hours
- Phase 4 (Analysis): 2 hours
- Buffer: 2 hours

**Deliverable**: Comprehensive evaluation showing when CAF beats BAM for streaming analytics

---

## Next Steps

1. Start with Phase 1 benchmarks (quick validation)
2. If promising, proceed to Phase 2-4
3. If not, document findings and pivot

**Goal**: Honest evaluation of CAF's streaming analytics value proposition

---

**Document Version**: 1.0
**Last Updated**: November 11, 2025
