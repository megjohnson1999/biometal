# CAF Decision Log

**Purpose**: Record major design and implementation decisions with rationale
**Audience**: Future self, reviewers, community

---

## Decision 1: Pursue CAF as Research Project

**Date**: November 10, 2025
**Context**: Completed BAM implementation, considering next steps

**Decision**: Pursue CAF as independent research project with publication target

**Alternatives considered**:
1. Focus on optimizing existing BAM implementation
2. Add CRAM support
3. Focus on other formats (VCF, BCF)
4. CAF research project ← **CHOSEN**

**Rationale**:
- **Innovation**: First ARM-native columnar alignment format
- **Scientific merit**: Publishable novel contribution
- **biometal mission**: Aligns with ARM-native, evidence-based philosophy
- **Timing**: Good opportunity after BAM completion

**Risks**:
- Community may not adopt (acceptable - research contribution)
- 8 weeks is aggressive (mitigate with phased approach)
- Storage overhead may be >2× (monitor, document)

**Success criteria**:
- 5-10× speedup demonstrated (N=30, p < 0.05)
- Lossless conversion validated
- Manuscript accepted

---

## Decision 2: Block Size = 10,000 Records

**Date**: November 10, 2025 (from prior experiments)
**Context**: Choosing granularity for columnar blocks

**Decision**: Use 10,000 records per block

**Alternatives considered**:
1. 1,000 records (finer granularity)
2. 10,000 records ← **CHOSEN** (from OPTIMIZATION_RULES.md Rule 2)
3. 100,000 records (coarser)

**Rationale**:
- **Evidence**: Entry 027 (1,440 measurements) validates 10K optimal
- **NEON**: Large enough for vectorization benefits
- **Memory**: Small enough for in-memory processing (~5-10 MB)
- **Parallelism**: Good granularity for rayon

**Implementation**:
```rust
const DEFAULT_BLOCK_SIZE: usize = 10_000;
```

**Validation plan**: Benchmark 1K, 10K, 100K in Phase 2

---

## Decision 3: Pre-decoded Sequences (ASCII not 4-bit)

**Date**: November 10, 2025
**Context**: How to store sequences in CAF

**Decision**: Store sequences as ASCII ('A', 'C', 'G', 'T'), not 4-bit encoded

**Alternatives considered**:
1. 4-bit encoding like BAM (saves space, burns CPU)
2. ASCII pre-decoded ← **CHOSEN**
3. 2-bit encoding (4 bases in 2 bits)

**Rationale**:
- **CPU vs Storage trade-off**: Storage is cheap (~$10/TB), CPU expensive
- **NEON optimization**: ASCII enables direct NEON operations (no unpacking)
- **Compression**: lz4 compresses ASCII ~2-3× → similar final size to 4-bit
- **Decompression speed**: lz4 extremely fast (>GB/s)

**Trade-offs**:
- 2× raw storage (mitigated by lz4 compression)
- Simpler code (no unpacking logic)
- **Net result**: 1.5-2× larger files, but 10-25× faster operations

**Validation**:
- Benchmark ASCII vs 4-bit performance (Phase 2)
- Measure actual compression ratios (Phase 3)

---

## Decision 4: Modern Compression (zstd/lz4 not gzip)

**Date**: November 10, 2025
**Context**: Which compression algorithms to use

**Decision**: Use zstd for integers, lz4 for sequences, raw for qualities

**Alternatives considered**:
1. gzip (like BAM) - compatible but slow
2. zstd/lz4 ← **CHOSEN** - modern, faster
3. Snappy, brotli, etc.

**Rationale**:
- **zstd**: 2-3× faster decompression than gzip, similar/better ratio
- **lz4**: Extremely fast (>GB/s), good for random access
- **Platform support**: Both widely available (Rust crates)

**Column-specific choices**:
- Positions, flags, CIGAR: **zstd level 3** (good compression, fast)
- Sequences: **lz4** (speed priority)
- Qualities: **raw** (incompressible, high entropy)

**Validation**:
- Benchmark vs gzip (Phase 3)
- Measure actual ratios on real data

---

## Decision 5: Block-Level Indexing (not base-level)

**Date**: November 10, 2025
**Context**: Random access strategy

**Decision**: Block-level index (10,000 records granularity)

**Alternatives considered**:
1. Base-level index (like BAI) - complex, unnecessary
2. Block-level ← **CHOSEN**
3. No index (sequential only)

**Rationale**:
- **Use case**: CAF for analytical queries, not genome browsers
- **Simplicity**: Block metadata sufficient for chr:start-end queries
- **Overhead**: Minimal (one BlockMeta per block)
- **Query speed**: Fast enough for batch processing

**Not suitable for**: Single-base genome browsing (use BAM for that)

---

## Decision 6: Lossless BAM Conversion Mandatory

**Date**: November 10, 2025
**Context**: Correctness requirements

**Decision**: BAM ↔ CAF conversion MUST be 100% lossless

**Rationale**:
- **Scientific integrity**: Can't lose data
- **Adoption**: Users need confidence in correctness
- **Validation**: Enables differential testing against BAM

**Implementation**:
- Comprehensive round-trip testing
- Property-based testing (proptest)
- Byte-level comparison

**Success criterion**: 100% lossless on 1,000+ diverse BAM files

---

## Decision 7: Target Bioinformatics Journal

**Date**: November 10, 2025
**Context**: Publication strategy

**Decision**: Target *Bioinformatics* (Oxford) - Application Note format

**Alternatives considered**:
1. Bioinformatics (Oxford) ← **CHOSEN** (IF=~7, prestigious)
2. BMC Bioinformatics (IF=~3, open access)
3. GigaScience (IF=~7, reproducibility focus)
4. arXiv only (no peer review)

**Rationale**:
- **Prestige**: Bioinformatics is top journal in field
- **Format**: Application Note (2,000 words) fits CAF scope
- **Audience**: Reaches bioinformatics community
- **Timeline**: 4-6 months review → publication Q2 2026

**Backup plan**: If rejected, submit to BMC Bioinformatics or GigaScience

---

## Decision 8: 8-Week Timeline

**Date**: November 10, 2025
**Context**: Project scheduling

**Decision**: 8-week research timeline (Nov 10 - Jan 10)

**Rationale**:
- **Phases**: Logical breakdown (Impl, NEON, Val, Pub)
- **Target**: Submit before Q1 2026
- **Aggressive but achievable**: Prior experiments provide foundation

**Risks**:
- Implementation may take longer → Can extend to 10-12 weeks
- Benchmarking delays → Automate early
- Publication writing → Start early, iterate

**Mitigation**: Weekly decision points, can adjust timeline

---

## Decision Template

*For future decisions, use this template:*

## Decision N: [Title]

**Date**: YYYY-MM-DD
**Context**: [Why was this decision needed?]

**Decision**: [What did you decide?]

**Alternatives considered**:
1. Option A
2. Option B ← **CHOSEN**
3. Option C

**Rationale**:
- Reason 1
- Reason 2
- Reason 3

**Trade-offs**:
- Pro 1
- Con 1

**Validation**:
- How will you validate this decision?

---

**Status**: Active
**Last Updated**: November 10, 2025
**Next Decision**: TBD
