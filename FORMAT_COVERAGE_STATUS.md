# biometal File Format Coverage (November 15, 2025)

## Current Status: Reading & Writing

### ✅ Fully Implemented (Read + Write)

| Format | Read | Write | Python Bindings | Status |
|--------|------|-------|-----------------|--------|
| **FASTQ** | ✅ v1.0.0 | ✅ v1.9.0 | ✅ Read/Write | Production |
| **FASTA** | ✅ v1.0.0 | ✅ v1.9.0 | ✅ Read/Write | Production |
| **SAM** | ✅ v1.7.0 | ✅ v1.7.0 | ❌ (PyO3 issue) | Rust only |
| **BAM** | ✅ v1.4.0 | ✅ v1.8.0 **NEW** | Read ✅, Write ❌ (PyO3 issue) | Rust production, Python read-only |
| **GFA** | ✅ v1.8.0 | ✅ v1.8.0 | Read ✅, Write ❌ (PyO3 issue) | Rust production, Python read-only |
| **VCF** | ✅ v1.8.0 | ✅ v1.9.0 | ✅ Read/Write | Production |
| **BED** (3/6/12) | ✅ v1.8.0 | ✅ v1.9.0 | ✅ Read/Write | Production |
| **narrowPeak** | ✅ v1.10.0 | ✅ v1.9.0 | ✅ Read/Write | Production |
| **GFF3** | ✅ v1.8.0 | ✅ v1.9.0 | ✅ Read/Write | Production |
| **GTF** | ✅ v1.10.0 | ✅ v1.9.0 | ✅ Read/Write | Production |
| **PAF** | ✅ v1.10.0 | ✅ v1.9.0 | ✅ Read/Write | Production |

### Index Formats

| Format | Read | Write | Status |
|--------|------|-------|--------|
| **BAI** (BAM index) | ✅ v1.6.0 | ❌ | Read-only (sufficient) |
| **FAI** (FASTA index) | ✅ v1.9.0 | ❌ | Read-only (sufficient) |
| **TBI** (Tabix) | ✅ v1.9.0 | ❌ | Read-only (sufficient) |
| **CSI** (Coordinate-sorted index) | ⏳ Partial | ❌ | Not complete |

### ⏳ Partial Implementation

| Format | Status | Next Steps | Effort |
|--------|--------|------------|--------|
| **CRAM** | Phase 1 Complete (v1.11.0+) | Phase 2: Full reference reconstruction + multi-codec | 20-30h |

### ❌ Not Implemented

| Format | Priority | Effort | Reason |
|--------|----------|--------|---------|
| **BCF** (Binary VCF) | MEDIUM | 30-40h | Binary compression, BGZF |
| **CSI** (Complete) | LOW | 20-30h | Less common than BAI |

---

## What We Just Completed (This Session)

### 🚀 CRAM Reader (v1.11.0) - Phase 2 In Progress ⏳
- **Strategic Decision**: Pivoted to **native, zero-dependency** ARM-optimized implementation
  - ❌ Abandoned noodles-cram (API compatibility issues, no ARM optimization)
  - ✅ Native implementation = full control + NEON optimizations
  - ✅ Minimal dependencies (bzip2, xz2 for additional codecs)
  - ✅ Flagship feature demonstrating biometal's ARM-native capabilities

- **Phase 1 Status**: ✅ COMPLETE (November 15, 2025)
  - Complete module structure (2,100+ lines)
  - File definition parsing (magic, version, file ID)
  - Container and slice structure parsing
  - ITF-8/LTF-8 variable-length integer decoding
  - Block decompression (gzip, bzip2, lzma)
  - Basic record iteration with placeholder data
  - **38 tests passing** (100% pass rate, +2 codec tests)

- **Phase 2 Status**: ⏳ IN PROGRESS (November 15, 2025, ~12 hours invested)
  - ✅ Multi-codec support (gzip, bzip2, lzma) - 4 hours
  - ✅ Reference FASTA integration (FAI index loading) - 4 hours
  - ⏳ Basic reference reconstruction (fetch reference subsequences) - 4 hours
  - ❌ Full CRAM feature decoding (substitutions, insertions, deletions) - 8-10 hours remaining
  - ❌ Quality score decoding from external blocks - 2-3 hours remaining
  - ❌ Full tag support (A, i, f, Z, H, B types) - 3-4 hours remaining
  - **618 total library tests** passing

- **Current Capability**:
  - Read CRAM file structure ✅
  - Decompress blocks (gzip, bzip2, lzma) ✅
  - Load reference FASTA with FAI index ✅
  - Fetch reference sequences for alignment regions ✅
  - Return records with reference-based sequences (simplified) ✅
  - **Limitation**: Sequences are reference subsequences, not actual reads with variations

- **Why Native?**:
  1. **First ARM-optimized CRAM reader** - 16-25× NEON speedup potential
  2. **Minimal external dependencies** - only compression codecs
  3. **Perfect streaming** - designed for constant ~5 MB memory
  4. **Full control** - optimize for biometal's architecture

- **Remaining Work** (20-30 hours, 3-5 days):
  - **Phase 2 Full** (20-30 hours): Complete CRAM feature decoding + tags
  - **Phase 3** (20-30 hours): ARM NEON optimizations (target: 2-3× faster than samtools)

### ✅ BAM Writer (v1.8.0) - Previous Session
- **Rust Implementation**: Production-ready
  - 14/14 tests passing (6 unit + 8 integration)
  - Round-trip verification with 100K+ records
  - Full BGZF compression with cloudflare_zlib
  - Documentation and examples

- **Python Bindings**: ❌ Known PyO3 Issue
  - Implemented but not appearing in module
  - Same issue as SAM reader, GFA writer
  - Documented in KNOWN_ISSUES.md
  - Rust users unaffected

---

## PyO3 Registration Mystery (3 Affected Classes)

**Pattern**: Correctly implemented Python bindings don't appear in compiled module

1. ❌ `PySamReader` (SAM reading)
2. ❌ `PyGfaWriter` (GFA writing)
3. ❌ `PyBamWriter` (BAM writing) **← Just discovered**

**Common factors**:
- All use `Option<T>` wrapping
- All contain enum types internally
- All compile without errors
- Symbols missing from `.so` file

**Impact**:
- Rust users: ✅ Full functionality
- Python users: ❌ Missing 3 features (workaround: use Rust)

**See**: `KNOWN_ISSUES.md` for full analysis

---

## What's Left: Next Priorities

### Option 1: Critical Format Gaps

#### 1.1. **CRAM Reader** (HIGH Priority)
- **Effort**: 80-120 hours (2-3 weeks)
- **Why Critical**:
  - 1000 Genomes uses CRAM exclusively
  - 3-5× smaller than BAM
  - Required for most modern datasets
- **Complexity**: High
  - Reference-based compression
  - Multiple compression codecs
  - Spec compliance (v3.0)
- **Deliverable**: Read-only (writing less important)

#### 1.2. **CSI Index Completion** (MEDIUM Priority)
- **Effort**: 20-30 hours (3-5 days)
- **Why Useful**:
  - Supports >2GB chromosomes (BAI limited to 512MB)
  - Required for some large genomes
  - Relatively straightforward extension of BAI
- **Status**: Partially implemented, needs completion
- **Deliverable**: Read-only (sufficient)

#### 1.3. **BCF Format** (MEDIUM Priority)
- **Effort**: 30-40 hours (5-7 days)
- **Why Useful**:
  - Binary VCF with compression
  - Faster parsing than VCF text
  - Common in variant calling pipelines
- **Complexity**: BGZF + binary encoding

---

### Option 2: Fix PyO3 Issues (MEDIUM Priority)

- **Effort**: 10-20 hours (investigation + fix)
- **Impact**: Unblocks Python users for 3 features
- **Approach**:
  1. Test with different PyO3 versions (0.26, 0.28)
  2. Minimal reproduction case
  3. Compare working vs failing implementations
  4. Report to PyO3 project if needed
- **Risk**: May not be solvable (PyO3 limitation)

---

### Option 3: Polish & Performance (LOW Priority)

- **Performance benchmarking** (BAM writer vs samtools)
- **Memory profiling** (ensure constant 5 MB)
- **Cross-platform validation** (Graviton, x86_64)
- **Documentation improvements**

---

## Recommendation

Based on current status, I recommend **Option 1.2 + 1.1** (in that order):

### Phase 1: CSI Index Completion (Quick Win)
- **Timeline**: 3-5 days
- **Effort**: 20-30 hours
- **Value**: Completes index format coverage
- **Risk**: Low (extension of existing BAI code)

### Phase 2: CRAM Reader (High Value)
- **Timeline**: 2-3 weeks
- **Effort**: 80-120 hours
- **Value**: Unblocks 1000 Genomes dataset access
- **Risk**: Medium (complex but well-specified)

### Phase 3: Python Issues (Optional)
- Investigate PyO3 issue after CRAM
- Lower priority than format coverage
- May resolve itself with PyO3 updates

---

## Strategic Consideration

**Two paths forward**:

**Path A: Complete Core Format Coverage** (Recommended)
- ✅ CSI Index (3-5 days)
- ✅ CRAM Reader (2-3 weeks)
- ✅ BCF Format (5-7 days)
- **Result**: biometal supports ALL major formats
- **Timeline**: 4-5 weeks total
- **Then**: Pivot to GPU/ML work (PROJECT_TODOS.md)

**Path B: Skip to GPU/ML Work** (Alternative)
- Start Week 1 of PROJECT_TODOS.md (Smith-Waterman GPU)
- Accept missing CRAM, BCF, CSI
- Users must convert formats (inconvenient)
- **Risk**: Incomplete core library

**My Recommendation**: Path A
- Complete format coverage is foundational
- CRAM is critical for modern datasets
- 4-5 weeks to finish vs years of tech debt
- Then GPU/ML work builds on solid base

---

## Summary

**What's Done**:
- ✅ 11 formats with full read/write support
- ✅ 3 index formats (read-only)
- ✅ BAM writer (Rust, just completed)

**What's Missing** (Critical):
- ❌ CRAM reader (1000 Genomes blocker)
- ❌ CSI index (partial implementation)
- ❌ BCF format (binary VCF)

**Python Bindings**:
- ✅ 9/12 formats working
- ❌ 3/12 formats blocked by PyO3 issue

**Next Action**:
Complete CSI Index (quick win, 3-5 days), then CRAM Reader (2-3 weeks)
