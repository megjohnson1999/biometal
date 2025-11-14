# biometal Planning Documentation Index

**Last Updated**: November 13, 2025
**Status**: Format Library Sprint approved - Phase 2 begins
**Current Phase**: Phase 2 (Format Development, 12-16 weeks)

---

## Current Strategic Plan (Active)

### 1. PHASE2_FORMAT_LIBRARY_SPRINT.md ⭐ **PRIMARY PLAN**
**Purpose**: Comprehensive format library development plan
**Status**: ✅ **APPROVED** - Ready to begin
**Duration**: 12-16 weeks (Nov 2025 - Mar 2026)
**Covers**:
- Strategic context and rationale
- Format tier system (Tier 1: BED, GFA, VCF, GFF3)
- Format development framework (apply proven optimizations)
- GFA format specifics (path extraction, graph analysis)
- Testing strategy per format
- Optimization decision framework
- Success criteria and deliverables

**Use When**: Planning format work, understanding Phase 2 goals

---

## Strategic Context Documents

### 2. STRATEGIC_PIVOT_PLAN.md
**Purpose**: Analysis of post-Phase 1 strategic options
**Status**: ✅ **DECISION MADE** - Format Library Sprint chosen
**Covers**:
- Current situation assessment
- 4 strategic options analyzed (Feature Expansion, Community, Quality, Pause)
- Rationale for Format Library Sprint
- Why Rules 3+4 optimization is not the path forward

**Use When**: Understanding why we pivoted to format development

### 3. RULES_3_4_REALITY_CHECK.md
**Purpose**: Critical correction to original Phase 2 plan
**Status**: ✅ Complete - Important reference
**Covers**:
- Why Rules 3+4 are not viable (Rule 3: 0.77-0.84× slowdown, Rule 4: ~1% benefit)
- Architectural conflict (streaming vs parallelism)
- Honest assessment of planning error
- Recommendations for alternative directions

**Use When**: Understanding why performance optimization is not Phase 2 focus

### 4. PHASE2_TRANSITION.md
**Purpose**: Documents transition from Apple Silicon research back to core roadmap
**Status**: ⚠️ **ARCHIVED** - Contains incorrect Rules 3+4 claims (see correction at top)
**Covers**:
- Apple Silicon research summary (Neural Engine, GPU)
- Strategic decision to return to core roadmap
- ⚠️ Original (invalid) Phase 2 plan for Rules 3+4

**Use When**: Understanding Apple Silicon research history

---

## Development Guidance

### 5. CLAUDE.md
**Purpose**: Development guide for AI-assisted sessions
**Status**: ✅ Updated with strategic planning focus (Nov 13, 2025)
**Covers**:
- Mission and core principles
- Evidence-based design guidelines
- ARM-native with portable fallback patterns
- Production quality requirements
- Recent research (Apple Silicon archived, CAF archived)
- Recent optimization (cloudflare_zlib v1.7.0)
- Current work: Strategic planning phase
- Future work: Pending Phase 2 direction (now decided)

**Use When**: Starting development sessions, architectural decisions

### 6. OPTIMIZATION_RULES.md
**Purpose**: Evidence-based optimization rules from ASBB
**Status**: ✅ Foundational reference
**Covers**:
- 6 optimization rules (1,357 experiments, N=30)
- Rule 1: NEON SIMD (16-25× speedup) ✅ Implemented
- Rule 2: Block processing ✅ Implemented
- Rule 3: Parallel BGZF ❌ DISABLED (0.77-0.84× slowdown, conflicts with streaming)
- Rule 4: Smart mmap ⏳ Optional (~1% benefit for compressed files)
- Rule 5: Streaming (99.5% memory reduction) ✅ Implemented
- Rule 6: Network streaming ✅ Implemented

**Use When**: Justifying optimization decisions, understanding proven techniques

---

## Completed Phase Documents

### 7. PHASE1_PROGRESS_REPORT.md
**Purpose**: Phase 1 (Consolidation) final report
**Status**: ✅ **COMPLETE** - All 4 weeks finished
**Covers**:
- Week 1: Documentation Sprint (40,000+ words)
- Week 2: Performance Benchmarking (vs samtools/pysam)
- Week 3: Community Building (blog post, GitHub templates)
- Week 4: Quality Assurance (11 property tests, 2 audits)
- 23 major deliverables
- Production-ready status

**Use When**: Understanding Phase 1 accomplishments, referencing deliverables

---

## Research Archives (Completed Work)

### Apple Silicon Research (Nov 4-13, 2025)
**Location**: `research/apple-silicon/`
**Key Document**: `RESEARCH_SUMMARY.md`
**Status**: ✅ Complete - Archived
**Outcome**:
- Neural Engine: 2,940× slowdown (batch-oriented vs streaming mismatch)
- GPU Smith-Waterman: 1.2-1.4× speedup (modest gains)
- Infrastructure preserved (feature-gated: `--features gpu`)
**Lessons**: Hardware-software fit matters, evidence-based > speculative

### CAF Columnar Format Research (Nov 4-11, 2025)
**Location**: `research/caf-format/implementation/`
**Key Document**: `CAF_FINAL_REPORT.md`
**Status**: ✅ Complete - Archived
**Outcome**: 1.6× larger files, 1.4× column-selective speedup (vs target 0.5-1.0×, 3-5×)
**Lessons**: Data layout affects performance, compiler auto-vec can beat explicit SIMD

### Compression Backend Investigation (Nov 12-13, 2025)
**Location**: Root directory
**Key Documents**:
- `BACKEND_COMPARISON_FINDINGS.md`
- `COMPRESSION_INVESTIGATION_FINDINGS.md`
**Status**: ✅ Complete - Integrated in v1.7.0
**Outcome**: cloudflare_zlib provides 1.67× decompression, 2.29× compression speedups

---

## Document Hierarchy

```
Current Phase
└── PHASE2_FORMAT_LIBRARY_SPRINT.md ⭐ PRIMARY (12-16 weeks)

Strategic Context
├── STRATEGIC_PIVOT_PLAN.md (Why format library sprint)
├── RULES_3_4_REALITY_CHECK.md (Why not Rules 3+4)
└── PHASE2_TRANSITION.md (Apple Silicon archive, archived)

Development Guidance
├── CLAUDE.md (Development patterns, current status)
├── OPTIMIZATION_RULES.md (Evidence base, 1,357 experiments)
└── README.md (Public messaging, roadmap)

Completed Phases
└── PHASE1_PROGRESS_REPORT.md (4 weeks, 23 deliverables)

Research Archive
├── research/apple-silicon/ (Neural Engine, GPU)
├── research/caf-format/implementation/ (Columnar format)
└── BACKEND_COMPARISON_FINDINGS.md (Compression)
```

---

## Quick Reference

### Starting Phase 2 Work
1. Read **PHASE2_FORMAT_LIBRARY_SPRINT.md** (format development plan)
2. Review **OPTIMIZATION_RULES.md** (proven optimizations to apply)
3. Check **CLAUDE.md** (development patterns)

### Planning Weekly Work
1. Check **PHASE2_FORMAT_LIBRARY_SPRINT.md** (current week's focus)
2. Review format-specific success criteria
3. Plan tasks (primitives → BED → GFA → VCF → GFF3)

### Making Optimization Decisions
1. Check **OPTIMIZATION_RULES.md** (is this a proven technique?)
2. Review **PHASE2_FORMAT_LIBRARY_SPRINT.md** (optimization decision framework)
3. Profile before optimizing (don't speculate)

### Understanding Strategic Decisions
1. Read **STRATEGIC_PIVOT_PLAN.md** (why format library sprint)
2. Review **RULES_3_4_REALITY_CHECK.md** (why not Rules 3+4)
3. Check **PHASE1_PROGRESS_REPORT.md** (what's already built)

---

## Key Metrics & Status

**Evidence Base**: 1,357 experiments, 40,710 measurements (N=30)
**Released**: v1.0.0 - v1.7.0 (Nov 5-13, 2025)
**Tests**: 582 passing (100% pass rate)
**Performance**: 16-25× NEON speedup, 92 MiB/s BAM parsing
**Current Formats**: FASTQ, FASTA, BAM, SAM, BAI
**Phase 2 Target Formats**: BED, GFA, VCF, GFF3, FAI, TBI
**Expected Completion**: v2.0.0 (March 2026)

---

## Revision History

- **Nov 13, 2025**: Updated for Format Library Sprint
  - Added PHASE2_FORMAT_LIBRARY_SPRINT.md as primary plan
  - Added RULES_3_4_REALITY_CHECK.md (critical correction)
  - Updated STRATEGIC_PIVOT_PLAN.md (decision made)
  - Marked PHASE2_TRANSITION.md as archived (invalid)
  - Reorganized for Phase 2 focus

- **Nov 13, 2025** (earlier): Initial creation post-strategic pivot
  - Consolidated planning documents
  - Archived pre-pivot documents
  - Established documentation hierarchy

---

**Next Review**: After completing format primitives (Week 2, late November 2025)
