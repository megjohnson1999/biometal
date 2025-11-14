# Strategic Discussion: Neural Engine Use Cases

**Date**: November 13, 2025
**Context**: Week 2 Neural Engine research complete
**Question**: Are ML use cases (adapter detection, read classification, variant calling) part of biometal's core development?

---

## The Question

Week 2 Neural Engine research identified better use cases than simple quality filtering:

**Potential Neural Engine Applications**:
- ✅ Adapter/contamination detection (complex pattern recognition)
- ✅ Read classification (multiclass ML: human/bacterial/viral)
- ✅ Sequence-only quality prediction (predict quality from sequence + position)
- ✅ Variant calling assistance (ML-assisted variant detection)

**Should we pursue these?** Or are they out of scope for biometal's core mission?

---

## Current Roadmap Analysis

### What IS Core Development

**Phase 1: Consolidation** (Weeks 1-4, current):
- Documentation (✅ complete)
- Benchmarking (✅ complete)
- Community building (⏳ Week 3)
- Quality assurance (⏳ Week 4)

**Phase 2: High-ROI Performance** (Weeks 5-9):
- Rule 3: Parallel BGZF decompression (6.5× speedup)
- Rule 4: Smart mmap (2.5× additional speedup)
- Combined: 16× BAM parsing improvement (55 → 895 MiB/s)

**Phase 3: Strategic Expansion** (Weeks 10-14):
- Format expansion (CRAM, VCF, CSI index)
- Horizontal expansion (GFF/GTF, BED)
- Community-driven features

### What is NOT in Current Roadmap

**ML/AI Features**:
- ❌ Adapter detection (not mentioned)
- ❌ Read classification (not mentioned)
- ❌ Quality prediction (not mentioned)
- ❌ Variant calling (not mentioned)

**GPU Features**:
- ❌ Smith-Waterman GPU (not mentioned, but Week 1 showed 771× speedup!)

**Strategic Pivot Exploration** (Week 1-2):
- This was a 2-week exploration outside the core roadmap
- Not planned in original Phase 1-3 timeline

---

## biometal's Core Mission

From CLAUDE.md:

> **Mission**: Democratize bioinformatics by enabling 5TB dataset analysis on consumer hardware through:
> - Streaming architecture: Constant ~5 MB memory
> - ARM-native performance: 16-25× NEON speedup
> - Network streaming: Analyze without downloading
> - Evidence-based optimization: Every rule validated experimentally

**Focus**: Infrastructure for efficient data access and processing
**Not**: Bioinformatics algorithms or ML-powered analysis

---

## Scope Analysis

### In Scope (Core Infrastructure)

1. **I/O Acceleration**:
   - Streaming parsers (FASTQ, FASTA, BAM)
   - Compression (BGZF, parallel decompression)
   - Network streaming (HTTP, SRA)
   - Indexing (BAI, CSI)

2. **ARM Optimization**:
   - NEON SIMD operations
   - Base counting, GC content, quality filtering
   - Sequence operations (reverse complement)

3. **Performance Rules**:
   - Rule 1: ARM NEON (✅ implemented)
   - Rule 2: Block-based processing (✅ implemented)
   - Rule 3: Parallel BGZF (⏳ planned Phase 2)
   - Rule 4: Smart mmap (⏳ planned Phase 2)
   - Rule 5: Streaming (✅ implemented)
   - Rule 6: Network streaming (✅ implemented)

### Out of Scope (Analysis Algorithms)

1. **Variant Calling**:
   - Complex bioinformatics algorithm
   - Requires domain expertise
   - Competitive tools exist (GATK, bcftools, DeepVariant)
   - Not an infrastructure problem

2. **Read Classification**:
   - Species identification (Kraken, Centrifuge)
   - Contamination detection (FastQ Screen, Kraken)
   - Specialized domain knowledge required

3. **Adapter Detection**:
   - Specialized task (Cutadapt, Trimmomatic)
   - Not infrastructure-level

4. **ML/AI Features**:
   - Requires training data curation
   - Model maintenance and updates
   - Specialized ML expertise
   - Different value proposition

---

## Strategic Options

### Option 1: Stay Focused (Recommended) ✅

**What**: Archive Neural Engine research, return to core roadmap

**Rationale**:
- biometal = infrastructure library (parsers, I/O, NEON optimization)
- ML features are analysis algorithms (different scope)
- Limited resources (solo development)
- Core roadmap already ambitious (14 weeks planned)
- Phase 2 (Rules 3+4) offers 16× speedup (high ROI)

**Action**:
1. ✅ Archive Neural Engine research (complete)
2. Integrate GPU Smith-Waterman (Week 1 success, 771× speedup)
3. Continue Phase 1: Community (Week 3) + QA (Week 4)
4. Execute Phase 2: Rules 3+4 (Weeks 5-9)

**Outcome**: Focused library, clear value proposition, achievable timeline

### Option 2: Expand Scope (High Risk) ⚠️

**What**: Add ML features as "Phase 4: ML-Powered Analysis"

**Rationale**:
- Neural Engine infrastructure already built
- Unique differentiation (no other tools use Neural Engine)
- Modern ML techniques for bioinformatics

**Challenges**:
1. **Resource Burden**: Training data, model maintenance, updates
2. **Scope Creep**: Moving from infrastructure to algorithms
3. **Competition**: DeepVariant, Kraken, established tools
4. **Timeline Impact**: Delays core roadmap (Rules 3+4)
5. **Mission Drift**: biometal becomes "ML analysis tool" not "infrastructure library"

**Risks**: Overextension, delayed core features, unclear value proposition

### Option 3: Hybrid Approach (Compromise) 🤔

**What**: Provide ML infrastructure, let users bring models

**Implementation**:
- Keep `src/ml/` module (ONNX Runtime integration)
- Document how to deploy custom ONNX models
- Provide examples (adapter detection, classification)
- Users train their own models, biometal provides inference

**Rationale**:
- Infrastructure-focused (aligned with mission)
- Enables ML without biometal maintaining models
- Users can experiment with Neural Engine
- Low maintenance burden

**Action**:
1. Keep Neural Engine infrastructure in codebase
2. Document ONNX model deployment
3. Provide reference examples
4. Users responsible for training/models

**Outcome**: ML-capable infrastructure, no algorithm maintenance

---

## Recommendation: **Option 1 (Stay Focused)** ✅

### Why Stay Focused?

1. **Clear Mission**: biometal = infrastructure library
2. **Resource Reality**: Solo development, ambitious core roadmap
3. **High-ROI Path**: Rules 3+4 offer 16× speedup (proven approach)
4. **Competitive Position**: Fast I/O + ARM optimization (unique)
5. **User Value**: Enable analysis workflows, don't provide algorithms

### What to Do with Neural Engine?

**Archive but Preserve**:
- ✅ Documentation complete (how to use Neural Engine if needed)
- ✅ Infrastructure code stays in `src/ml/` (feature-gated)
- ✅ Examples available (neural_quality.rs)
- Users can deploy custom ONNX models if desired
- biometal doesn't maintain ML models

**Focus on GPU Smith-Waterman**:
- 771× speedup is transformative
- Production-ready (Week 1 success)
- Clear use case (sequence alignment)
- Integrate into v1.7.0+

### What to Do with ML Use Cases?

**Defer to Community/Ecosystem**:
- Adapter detection → Cutadapt, Trimmomatic
- Read classification → Kraken, Centrifuge
- Variant calling → GATK, DeepVariant
- biometal provides fast I/O for these tools

---

## Discussion Points

### 1. Is biometal an Infrastructure Library or Analysis Tool?

**Infrastructure** (current):
- Parsers, I/O, streaming, ARM optimization
- Enables other tools to run faster
- Clear, narrow scope

**Analysis** (expansion):
- Algorithms, ML models, variant calling
- Competes with GATK, Kraken, DeepVariant
- Broad, complex scope

**Question**: Which path aligns with biometal's mission and resources?

### 2. What is biometal's Unique Value?

**Current Strengths**:
- ARM NEON optimization (16-25× speedup)
- Streaming architecture (5 MB constant memory)
- Evidence-based design (1,357 experiments)
- Network streaming
- Python bindings

**Adding ML**:
- Neural Engine integration (unique)
- But: Performance not competitive for simple ops
- Better ML tools exist (DeepVariant, etc.)

**Question**: Does ML enhance or dilute biometal's value proposition?

### 3. Resource Allocation

**Core Roadmap** (14 weeks):
- Phase 1: 4 weeks (50% complete)
- Phase 2: 5 weeks (Rules 3+4, 16× speedup)
- Phase 3: 5 weeks (format expansion)

**Adding ML Features** (estimated):
- Model training: 2-4 weeks per use case
- Validation: 1-2 weeks per model
- Maintenance: Ongoing
- Total: 15-25 weeks for 4 ML features

**Question**: Can solo developer maintain both infrastructure and ML models?

### 4. User Demand

**Evidence**:
- No user requests for ML features (yet)
- Core roadmap based on common workflows (I/O, parsing, alignment)
- GPU Smith-Waterman has clear demand (MSA, database search)

**Question**: Should we add ML features speculatively, or wait for user demand?

---

## My Recommendation

**Archive Neural Engine research, stay focused on core infrastructure.**

**Rationale**:
1. **Mission Alignment**: biometal = infrastructure, not algorithms
2. **Resource Reality**: Solo dev can't maintain both
3. **High ROI**: Rules 3+4 offer proven 16× speedup
4. **Clear Value**: Fast I/O + ARM optimization is unique
5. **Ecosystem**: Let specialized tools handle ML (DeepVariant, Kraken)

**GPU Smith-Waterman**: Integrate (771× speedup, production-ready)
**Neural Engine**: Archive (infrastructure preserved, no active development)
**ML Features**: Defer to community/ecosystem

**Next Steps**:
1. Integrate GPU Smith-Waterman into v1.7.0
2. Continue Phase 1: Community (Week 3) + QA (Week 4)
3. Execute Phase 2: Rules 3+4 (Weeks 5-9)

---

## Questions for Discussion

1. **Scope**: Is biometal an infrastructure library or analysis tool?
2. **Resources**: Can solo dev maintain infrastructure + ML models?
3. **Value**: Does ML enhance or dilute biometal's unique strengths?
4. **Demand**: Should we add ML speculatively or wait for user requests?
5. **GPU**: Should we integrate Smith-Waterman (771× speedup)?

---

**Status**: Awaiting strategic decision before proceeding
**Options**: Stay focused (recommended) | Expand scope | Hybrid approach
**Next Action**: Decide direction, then update roadmap accordingly
