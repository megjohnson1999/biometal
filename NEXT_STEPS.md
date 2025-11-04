# biometal: Next Steps Recommendation

**Current Status**: v0.1.0 tagged and pushed ✅
**Date**: November 4, 2025
**Review**: Planning documents, roadmap, evidence base

---

## Current State Analysis

### ✅ Week 1-2 COMPLETE (v0.1.0)

**Deliverables Shipped**:
- FASTQ/FASTA streaming parsers (constant ~5 MB memory)
- Paired-end support (synchronized R1/R2)
- ARM NEON operations: base_counting (16.7×), gc_content (20.3×), quality_filter (25.1×)
- Intelligent I/O: Parallel bgzip (6.5×) + smart mmap (2.5×) = 16.3× combined
- Evidence-based design: All optimizations validated with N=30 statistical rigor
- **65 tests** (41 unit + 12 property + 12 doctests), **0 failures**
- Comprehensive benchmarks with criterion
- Production-quality: Zero unwrap/expect, full error handling

**Architecture Compliance**:
- ✅ Rule 1: ARM NEON SIMD (16-25× speedup with scalar fallback)
- ✅ Rule 2: Block-based processing (10K records preserve NEON gains)
- ✅ Rule 3: Parallel bgzip (bounded 8-block, constant memory)
- ✅ Rule 4: Smart mmap (≥50 MB threshold, 2.5× speedup)
- ✅ Rule 5: Constant-memory streaming (~5 MB regardless of file size)
- ⏳ Rule 6: Network streaming (NOT YET IMPLEMENTED)

**What's Missing for Full v1.0**:
- Network streaming (Rule 6 - CRITICAL bottleneck mitigation)
- Python bindings (Week 5-6)
- Extended operations coverage
- Cross-platform validation (Graviton, x86_64)

---

## Strategic Options

### Option A: Week 3-4 Network Streaming (RECOMMENDED)

**Timeline**: Nov 4-15 (2 weeks)
**Target**: biometal v0.2.0

**Rationale**:
1. **Rule 6 Validation**: I/O bottleneck is **264-352× slower** than compute (Entry 028)
   - Without network streaming: NEON gives only 1.04-1.08× E2E speedup (I/O masks compute gains)
   - With I/O optimization: Projects to ~17× E2E speedup
   - **Conclusion**: Network streaming is CRITICAL, not optional

2. **Mission Alignment**: "Democratize bioinformatics" means analyzing 5TB datasets WITHOUT 5TB downloads
   - Current: Must download entire dataset (days/weeks on slow connections)
   - With streaming: Start analysis immediately, cache only what's needed
   - Impact: Makes large-scale genomics accessible to LMIC/small labs

3. **Natural Progression**: Builds on v0.1.0 foundation
   - DataSource abstraction already exists (Local only)
   - Compression pipeline ready for network integration
   - Streaming architecture designed for this

4. **Timeline Efficiency**: Gets v0.2.0 shipped before Thanksgiving
   - Week 3-4: Nov 4-15 (network streaming)
   - Week 5-6: Dec 2-13 (Python bindings)
   - v1.0 release: Dec 15 (on schedule)

**Implementation Plan**:

```
Week 3 (Nov 4-8): HTTP Streaming Foundation
├─ Day 1-2: DataSource::Http implementation
│  └─ src/io/network.rs (HTTP client with reqwest)
│  └─ Range request support for partial downloads
│  └─ Error handling (timeouts, retries, network failures)
│
├─ Day 3-4: Smart Caching
│  └─ LRU cache for bgzip blocks (configurable size)
│  └─ Cache eviction policies
│  └─ Memory-bounded guarantees
│
└─ Day 5: Integration + Testing
   └─ HTTP + compression pipeline integration
   └─ Unit tests for network layer
   └─ Error scenarios (offline, 404, timeouts)

Week 4 (Nov 11-15): Optimization + SRA
├─ Day 1-2: Background Prefetching
│  └─ Async prefetch next N blocks
│  └─ Hide network latency during processing
│  └─ Bounded concurrency
│
├─ Day 3-4: SRA Toolkit Integration
│  └─ src/io/sra.rs (SRA accession → HTTP URLs)
│  └─ Optional feature flag (sra)
│  └─ Examples with real SRA data
│
└─ Day 5: Polish + Release
   └─ Property tests for network edge cases
   └─ Benchmarks (local vs HTTP streaming)
   └─ Documentation + examples
   └─ Tag v0.2.0
```

**Key Files to Create**:
1. `src/io/network.rs` (~300 lines): HTTP streaming with range requests
2. `src/io/cache.rs` (~150 lines): LRU cache for bgzip blocks
3. `src/io/prefetch.rs` (~100 lines): Background prefetching
4. `src/io/sra.rs` (~200 lines): SRA accession support
5. `examples/network_streaming.rs`: Demo HTTP streaming
6. `tests/network_tests.rs`: Integration tests

**Dependencies to Add**:
```toml
[dependencies]
reqwest = { version = "0.11", features = ["stream", "blocking"] }
lru = "0.12"
tokio = { version = "1", features = ["rt-multi-thread"] }  # For async
```

**Expected Outcomes**:
- ✅ Rule 6 implemented: Network streaming addresses I/O bottleneck
- ✅ 5TB dataset analysis without 5TB download
- ✅ Smart caching balances memory and network
- ✅ SRA integration enables direct SRA accession streaming
- ✅ v0.2.0 ready for Week 5-6 Python bindings

**Risks**:
- Medium: Network code complexity (async, error handling)
- Low: Dependency on reqwest/tokio (well-established crates)
- Low: SRA API changes (can adapt)

---

### Option B: Week 5-6 Python Bindings (SKIP NETWORK)

**Timeline**: Nov 4-15 (2 weeks)
**Target**: biometal v0.2.0 (Python-ready)

**Rationale**:
- Get to Python ecosystem faster
- Enable Jupyter notebook demos
- Appeal to ML/data science users

**CONCERNS**:
- ❌ **Violates Rule 6**: I/O bottleneck remains unaddressed
- ❌ **Limited utility**: Python users can't access large datasets without download
- ❌ **Incomplete story**: "Fast library that can't handle large files" is weak messaging
- ❌ **Technical debt**: Have to retrofit network streaming later (harder)

**Recommendation**: **NOT RECOMMENDED**. Network streaming is architectural - do it now.

---

### Option C: Jump to v1.1 BAM/SAM (AMBITIOUS)

**Timeline**: 4+ weeks
**Target**: biometal v1.1.0 (BAM support + Metal GPU)

**Rationale**:
- BAM/SAM is high-value format
- Metal GPU breakthrough potential (world-first)
- More technically exciting than network plumbing

**CONCERNS**:
- ❌ **Scope creep**: BAM is complex (CIGAR, flags, optional fields)
- ❌ **Missing foundation**: Network streaming needed for large BAM files anyway
- ❌ **Timeline risk**: Metal GPU is R&D, not guaranteed success
- ❌ **Dependencies**: Requires noodles integration (new complexity)

**Recommendation**: **NOT RECOMMENDED**. Finish v1.0 FASTQ/FASTA first, then BAM.

---

### Option D: Consolidate v0.1.0 (REAL-WORLD TESTING)

**Timeline**: 1-2 weeks
**Target**: biometal v0.1.1 (bug fixes, polish)

**Activities**:
- Test with real ASBB datasets (tiny → vlarge)
- Cross-platform testing (Graviton, x86_64)
- Performance profiling and optimization
- Documentation improvements
- Community feedback (if public)

**CONCERNS**:
- ⚠️ **Pauses forward progress**: Delays v1.0 timeline
- ⚠️ **Limited gains**: v0.1.0 is already production-quality
- ⚠️ **Premature**: Better to consolidate after v0.2.0 (network streaming)

**Recommendation**: **DEFER**. Do this after v0.2.0 instead (more complete feature set).

---

## Final Recommendation: Option A (Week 3-4 Network Streaming)

### Why This is the Right Choice

1. **Evidence-Driven**: Rule 6 shows I/O bottleneck is CRITICAL (264-352×)
   - Network streaming isn't a nice-to-have, it's foundational
   - Addresses the primary performance bottleneck

2. **Mission-Critical**: "Democratize bioinformatics" requires network streaming
   - 5TB datasets become accessible without 5TB downloads
   - LMIC/small labs can analyze large-scale genomics
   - Infrastructure for cloud-native genomics

3. **Natural Progression**: Completes the architectural foundation
   - v0.1.0: Local file streaming ✅
   - v0.2.0: Network streaming ✅
   - v1.0.0: Python bindings + polish ✅

4. **Timeline Efficiency**: Ships v1.0.0 by Dec 15 (on schedule)
   - Week 3-4: Network streaming (Nov 4-15)
   - Week 5-6: Python bindings (Dec 2-13)
   - v1.0.0: Dec 15

5. **Technical Soundness**: Builds on solid v0.1.0 foundation
   - DataSource abstraction already exists
   - Compression pipeline ready for integration
   - Streaming architecture designed for this

### Immediate Next Steps (This Week)

**Day 1 (Today)**:
1. ✅ Review and approve this plan
2. Create Week 3-4 todo list
3. Set up network dependencies in Cargo.toml
4. Create `src/io/network.rs` skeleton

**Day 2-3 (Nov 5-6)**:
1. Implement DataSource::Http with reqwest
2. Add range request support
3. Unit tests for HTTP client

**Day 4-5 (Nov 7-8)**:
1. LRU cache implementation
2. Integration with compression pipeline
3. Property tests for caching

**Week 4 (Nov 11-15)**:
1. Background prefetching
2. SRA integration
3. Examples and documentation
4. Tag v0.2.0

---

## Alternative: 2-Phase Approach (if uncertain)

If you want to validate v0.1.0 before committing to network streaming:

**Phase 1: Quick Validation (2-3 days)**
- Run v0.1.0 against all ASBB datasets (tiny → vlarge)
- Verify performance claims hold
- Test on Graviton/x86_64 (if available)
- Fix any critical bugs

**Phase 2: Network Streaming (10-12 days)**
- Proceed with Option A timeline
- Benefit: Higher confidence in foundation
- Cost: Delays v0.2.0 by 2-3 days (minor)

---

## Recommended Decision

**Proceed with Week 3-4 Network Streaming (Option A)**

**Rationale Summary**:
- Evidence-based (Rule 6 validation)
- Mission-aligned (democratize bioinformatics)
- Architecturally sound (completes foundation)
- Timeline-efficient (v1.0.0 by Dec 15)
- Technically ready (v0.1.0 is solid)

**Expected Outcome**: biometal v0.2.0 by Nov 15 with network streaming, positioning us perfectly for Week 5-6 Python bindings and v1.0.0 release by Dec 15.

---

## Long-Term Vision Alignment

This recommendation aligns with the broader roadmap:

**v1.0 (Dec 15)**: Production FASTQ/FASTA library
- ✅ Streaming parsers
- ✅ ARM NEON operations
- ✅ Network streaming
- ✅ Python bindings

**v1.1 (Jan 12)**: BAM/SAM + Metal GPU breakthrough
- BAM/SAM parsing (noodles wrapper)
- NEON CIGAR optimization
- Metal GPU pileup (world-first)

**v1.2+ (Jan+)**: Comprehensive genomics toolkit
- Coverage, depth, variant calling
- Extended format support (BED, VCF, GTF)
- Neural Engine exploration

---

**Recommendation**: Start Week 3-4 Network Streaming tomorrow (Option A)

Let me know if you'd like to proceed with this plan, or if you'd prefer to explore any of the alternative options!
