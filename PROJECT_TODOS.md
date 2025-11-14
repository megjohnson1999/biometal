# biometal: Detailed Project TODOs (6-Month Plan)

**Strategic Pivot Date**: November 13, 2025
**Duration**: 6 months (24 weeks)
**Goal**: Apple Silicon-native bioinformatics + ML/BERT integration

**Overall Progress**: 0/24 weeks complete (0%)

---

## Quick Status Overview

| Phase | Weeks | Status | Progress |
|-------|-------|--------|----------|
| **Phase 1: GPU/Metal** | 1-8 | 🔴 Not Started | 0/8 weeks |
| **Phase 2: ML/BERT** | 9-16 | 🔴 Not Started | 0/8 weeks |
| **Phase 3: Primitives + Demo** | 17-24 | 🔴 Not Started | 0/8 weeks |

**Legend**: 🔴 Not Started | 🟡 In Progress | 🟢 Complete | ⏸️ Blocked | ❌ Cancelled

---

## Phase 1: GPU/Metal Exploration (Months 1-2, Weeks 1-8)

**Objective**: Identify Apple Silicon GPU wins for high-complexity operations

**Key Insight**: GPU only helps when complexity >0.55 AND NEON <2× AND batch ≥10K
- ASBB tested: complexity 0.20-0.61 (too simple for GPU)
- **New targets**: complexity >0.70 (alignment, pileup, variant calling)

**Decision Point**: End of Week 4
- ✅ If GPU achieves 10-50× speedup → Continue GPU work (Weeks 5+)
- ❌ If GPU <2× speedup → Pivot to primitives-focused (skip GPU work)

---

### Month 1: GPU Operations (Weeks 1-4)

#### **Week 1: Smith-Waterman Alignment on GPU (Part 1)** ⭐⭐⭐

**Status**: 🔴 Not Started
**Priority**: HIGH (Strategic pivot core work)
**Estimated Effort**: 20-25 hours

**Why this week**:
- Complexity >0.70 (dynamic programming)
- NEON limited (<2× due to irregular memory access)
- Proven in CUDA literature (10-50× speedup)
- Unified memory advantage (zero-copy CPU↔GPU)

**Tasks**:
- [ ] **Research** (3-4 hours):
  - [ ] Review CUDA Smith-Waterman implementations (literature survey)
  - [ ] Review Metal compute shader documentation
  - [ ] Design parallel Smith-Waterman algorithm
  - [ ] Determine data structures (scoring matrix, traceback)

- [ ] **Implementation** (10-12 hours):
  - [ ] Create `src/alignment/` module structure
  - [ ] Implement naive CPU Smith-Waterman (baseline)
  - [ ] Implement Metal compute shader (parallel version):
    ```metal
    kernel void smith_waterman_parallel(
        device const uint8_t* query_seqs [[buffer(0)]],
        device const uint8_t* reference_seqs [[buffer(1)]],
        device int* alignment_scores [[buffer(2)]],
        device uint8_t* cigar_operations [[buffer(3)]],
        uint2 gid [[thread_position_in_grid]]
    );
    ```
  - [ ] Implement Rust wrapper for Metal shader
  - [ ] Handle unified memory (MTLResourceOptions::StorageModeShared)
  - [ ] Implement scoring matrix (match, mismatch, gap)

- [ ] **Testing** (4-5 hours):
  - [ ] Unit tests: Known sequence pairs (expected scores)
  - [ ] Property test: GPU score == CPU score (proptest)
  - [ ] Edge cases: Empty sequences, identical sequences, no matches
  - [ ] Integration test: 100 random pairs

- [ ] **Documentation** (2-3 hours):
  - [ ] API documentation with examples
  - [ ] Algorithm explanation (dynamic programming)
  - [ ] Performance characteristics
  - [ ] Platform support (Mac vs Graviton/x86)

**Deliverables**:
- Naive CPU Smith-Waterman (baseline)
- Metal GPU Smith-Waterman (parallel)
- Unit tests (10+)
- Property tests (GPU == CPU)

**Blockers/Risks**:
- Metal shader complexity (mitigation: start simple, iterate)
- Unified memory issues (mitigation: follow Apple docs)

**Success Criteria**:
- [ ] GPU implementation compiles and runs
- [ ] All tests passing
- [ ] GPU score matches CPU score (correctness)

---

#### **Week 2: Smith-Waterman Alignment on GPU (Part 2)** ⭐⭐⭐

**Status**: 🔴 Not Started
**Priority**: HIGH
**Estimated Effort**: 20-25 hours
**Dependencies**: Week 1 complete

**Tasks**:
- [ ] **NEON Optimization** (6-8 hours):
  - [ ] Implement NEON-optimized CPU version (striped algorithm)
  - [ ] Benchmark: Naive vs NEON (N=30)
  - [ ] Analyze: Is NEON <2×? (Expect yes due to irregular access)

- [ ] **GPU Optimization** (6-8 hours):
  - [ ] Optimize Metal shader (reduce memory access, coalescing)
  - [ ] Batch processing (align multiple pairs in parallel)
  - [ ] Tune thread group size (optimal for Apple GPU)

- [ ] **Comprehensive Benchmarking** (6-8 hours):
  - [ ] Benchmark suite (N=30 for each):
    - [ ] Query length: 100bp, 500bp, 1000bp, 5000bp
    - [ ] Batch size: 1, 10, 100, 1000, 10000 pairs
    - [ ] Similarity: 90%, 70%, 50% (affects complexity)
  - [ ] Measure throughput (alignments/sec)
  - [ ] Measure latency (ms per alignment)
  - [ ] Statistical analysis: 95% CI, Cohen's d

- [ ] **CIGAR Generation** (2-3 hours):
  - [ ] Implement traceback for CIGAR string
  - [ ] Both CPU and GPU versions
  - [ ] Test: CIGAR correctness

**Deliverables**:
- NEON-optimized Smith-Waterman
- Optimized GPU Smith-Waterman
- Comprehensive benchmarks (N=30)
- CIGAR string generation

**Success Criteria**:
- [ ] GPU achieves 10-50× speedup vs naive CPU (target met)
- [ ] GPU achieves >5× speedup vs NEON CPU (GPU justified)
- [ ] Benchmarks: Mean, std dev, 95% CI, Cohen's d calculated

**Decision Point**:
- ✅ If GPU ≥10× faster: **HUGE SUCCESS** → Continue GPU work
- ⚠️ If GPU 5-10× faster: **MODERATE SUCCESS** → Continue with caution
- ❌ If GPU <5× faster: **FAIL** → Reconsider GPU strategy

---

#### **Week 3: Pileup Generation on GPU** ⭐⭐⭐

**Status**: 🔴 Not Started
**Priority**: HIGH
**Estimated Effort**: 20-25 hours
**Dependencies**: BAM parser (✅ v1.7.0), Smith-Waterman results

**Why this week**:
- Complexity ~0.60 (accumulation + filtering)
- Massive parallelism (millions of genomic positions)
- Large batches (100K+ reads >> 10K GPU threshold)
- Not tested in ASBB (unique to genomics)

**Tasks**:
- [ ] **Design** (3-4 hours):
  - [ ] Review samtools mpileup algorithm
  - [ ] Design parallel pileup (atomic operations)
  - [ ] Memory layout: Position array (atomic counters)
  - [ ] Filtering: MAPQ, quality thresholds

- [ ] **CPU Implementation** (5-6 hours):
  - [ ] Naive CPU pileup (baseline)
  - [ ] Streaming interface (constant memory)
  - [ ] Per-position: count, depth, quality stats
  - [ ] Integration with BAM parser

- [ ] **GPU Implementation** (8-10 hours):
  - [ ] Metal compute shader:
    ```metal
    kernel void pileup_accumulate(
        device const BamRecord* records [[buffer(0)]],
        device atomic_uint* pileup [[buffer(1)]],
        device atomic_uint* quality_sum [[buffer(2)]],
        constant uint& num_records [[buffer(3)]],
        constant uint& start_position [[buffer(4)]],
        uint gid [[thread_position_in_grid]]
    );
    ```
  - [ ] Atomic operations (thread-safe accumulation)
  - [ ] Parallel across reads (each thread processes one read)
  - [ ] Unified memory (zero-copy BAM records)

- [ ] **Testing & Benchmarking** (4-5 hours):
  - [ ] Unit tests: Known BAM regions (validate counts)
  - [ ] Compare: GPU pileup vs samtools mpileup (correctness)
  - [ ] Benchmark (N=30):
    - [ ] Region sizes: 1Kb, 10Kb, 100Kb, 1Mb, 10Mb
    - [ ] Coverage: 10×, 30×, 50×, 100×
    - [ ] Throughput: Positions/sec
  - [ ] Compare vs samtools mpileup (speedup)

**Deliverables**:
- CPU pileup generation (naive + NEON)
- GPU pileup generation (Metal)
- Integration with BAM parser
- Benchmarks vs samtools

**Success Criteria**:
- [ ] GPU pileup matches samtools mpileup (correctness)
- [ ] GPU achieves 20-100× speedup vs CPU (target met)
- [ ] Constant ~5 MB memory (streaming maintained)

---

#### **Week 4: GPU Results Analysis & Blog Post** ⭐⭐⭐

**Status**: 🔴 Not Started
**Priority**: HIGH (Decision point week)
**Estimated Effort**: 15-20 hours
**Dependencies**: Weeks 1-3 complete

**Tasks**:
- [ ] **Statistical Analysis** (4-5 hours):
  - [ ] Compile all benchmark results (Weeks 1-3)
  - [ ] Calculate statistics:
    - [ ] Mean, median, std dev for each operation
    - [ ] 95% confidence intervals
    - [ ] Cohen's d effect sizes
  - [ ] Create visualizations:
    - [ ] Speedup charts (naive vs NEON vs GPU)
    - [ ] Scaling curves (batch size, sequence length)
    - [ ] Energy consumption (if measured)

- [ ] **Comparative Analysis** (3-4 hours):
  - [ ] Compare vs CUDA literature (is Apple Silicon competitive?)
  - [ ] Compare vs samtools/other tools
  - [ ] Identify sweet spots (when GPU helps most)
  - [ ] Identify limitations (when GPU doesn't help)

- [ ] **Documentation** (3-4 hours):
  - [ ] Update CLAUDE.md (Phase 1 results)
  - [ ] Update STRATEGIC_PIVOT_PLAN.md (decision)
  - [ ] Document lessons learned
  - [ ] Document negative results (if any)

- [ ] **Blog Post** (5-7 hours):
  - [ ] Title: "GPU-Accelerated Bioinformatics on Apple Silicon: First Results"
  - [ ] Sections:
    - [ ] Motivation (why GPU for genomics?)
    - [ ] Methods (Smith-Waterman, pileup)
    - [ ] Results (10-50× speedup, visualizations)
    - [ ] Lessons learned (unified memory advantage)
    - [ ] Future work (variant calling, more operations)
  - [ ] Code examples
  - [ ] Publish to blog, Reddit, Twitter

**Deliverables**:
- Statistical analysis report
- Updated planning documents
- Blog post published

**CRITICAL DECISION POINT**:

**If GPU achieves 10-50× speedup** ✅:
- **SUCCESS!** Continue GPU work in Weeks 5-8
- Expand to variant calling GPU (Week 17)
- Write paper: "Apple Silicon GPU for Genomics"

**If GPU achieves 2-10× speedup** ⚠️:
- **MODERATE SUCCESS** - GPU helps but not revolutionary
- Complete GPU work as planned
- Focus more on ML/BERT primitives (higher novelty)

**If GPU achieves <2× speedup** ❌:
- **FAIL** - GPU doesn't justify complexity
- Document negative result (valuable for community)
- Pivot: Skip GPU work in Weeks 5+
- Focus: ML/BERT primitives + core primitives (higher ROI)

---

### Month 2: Neural Engine + Core Primitives (Weeks 5-8)

#### **Week 5: Neural Engine Quality Prediction (Part 1)** ⭐⭐⭐

**Status**: 🔴 Not Started
**Priority**: HIGH (Novel, not tested in ASBB)
**Estimated Effort**: 20-25 hours
**Dependencies**: None (independent of GPU results)

**Why this week**:
- Neural Engine designed for transformer/BERT inference
- 15.8 TOPS, <1ms latency, 0.5W energy
- Novel: Replace/augment Phred scores with ML predictions
- Unique to Apple Silicon (no equivalent on Graviton/x86)

**Tasks**:
- [ ] **Research & Design** (4-5 hours):
  - [ ] Review quality score prediction literature
  - [ ] Review CoreML documentation (Neural Engine)
  - [ ] Design model architecture:
    - [ ] Input: Sequence context (±50 bases)
    - [ ] Output: Quality score prediction (0-60)
    - [ ] Architecture: LSTM or Transformer?
  - [ ] Identify training data sources

- [ ] **Data Preparation** (5-6 hours):
  - [ ] Generate training data:
    - [ ] High-quality reads (aligned to reference)
    - [ ] Ground truth: Consensus quality from alignment
    - [ ] Features: Sequence context, position
  - [ ] Split: Train (80%) / Val (10%) / Test (10%)
  - [ ] Data augmentation: Reverse complement

- [ ] **Model Training** (6-8 hours):
  - [ ] Implement PyTorch model:
    ```python
    class QualityPredictor(nn.Module):
        def __init__(self):
            super().__init__()
            self.embedding = nn.Embedding(5, 32)  # ACGTN
            self.lstm = nn.LSTM(32, 128, num_layers=2, bidirectional=True)
            self.fc = nn.Linear(256, 1)  # Predict quality

        def forward(self, sequence):
            embedded = self.embedding(sequence)
            lstm_out, _ = self.lstm(embedded)
            quality = self.fc(lstm_out).squeeze(-1)
            return quality
    ```
  - [ ] Train model (10-20 epochs)
  - [ ] Validate on held-out data
  - [ ] Measure accuracy: MAE, RMSE, R²

- [ ] **CoreML Conversion** (2-3 hours):
  - [ ] Convert PyTorch → CoreML:
    ```python
    import coremltools as ct

    traced_model = torch.jit.trace(model, example_input)
    coreml_model = ct.convert(
        traced_model,
        inputs=[ct.TensorType(shape=(1, 100), dtype=int32)],
        convert_to="neuralnetwork",  # Neural Engine format
    )
    coreml_model.save("quality_predictor.mlmodel")
    ```
  - [ ] Test CoreML model (predictions match PyTorch)

- [ ] **Documentation** (2-3 hours):
  - [ ] Document training process
  - [ ] Document CoreML conversion
  - [ ] API documentation (Rust side)

**Deliverables**:
- Trained quality prediction model (PyTorch)
- CoreML model (Neural Engine format)
- Training/validation data
- Documentation

**Success Criteria**:
- [ ] Model achieves MAE <5 (vs Phred scores)
- [ ] CoreML model loads and runs
- [ ] Predictions: CoreML == PyTorch (correctness)

---

#### **Week 6: Neural Engine Quality Prediction (Part 2)** ⭐⭐⭐

**Status**: 🔴 Not Started
**Priority**: HIGH
**Estimated Effort**: 20-25 hours
**Dependencies**: Week 5 complete

**Tasks**:
- [ ] **Rust Integration** (8-10 hours):
  - [ ] Create `src/ml/` module
  - [ ] Objective-C bridge for CoreML:
    ```rust
    #[cfg(all(target_arch = "aarch64", target_os = "macos"))]
    mod coreml {
        use objc::*;
        pub struct CoreMLModel { /* ... */ }
        impl CoreMLModel {
            pub fn load(path: &str) -> Result<Self>;
            pub fn predict(&self, sequence: &[u8]) -> Result<Vec<f32>>;
        }
    }
    ```
  - [ ] Streaming integration:
    ```rust
    pub fn predict_quality_streaming(
        fastq_path: &str,
        model_path: &str,
    ) -> impl Iterator<Item = Result<(FastqRecord, Vec<f32>)>>;
    ```
  - [ ] CPU fallback (for non-Mac platforms)

- [ ] **Benchmarking** (6-8 hours):
  - [ ] Latency benchmarks (N=30):
    - [ ] Per-sequence prediction time
    - [ ] Batch prediction time (10, 100, 1000 sequences)
  - [ ] Throughput benchmarks:
    - [ ] Sequences/sec (Neural Engine vs CPU)
  - [ ] Energy measurement:
    - [ ] `powermetrics` on Mac (measure watts)
    - [ ] Compare: Neural Engine vs CPU inference
  - [ ] Memory: Constant 5 MB maintained?

- [ ] **Comparative Analysis** (4-5 hours):
  - [ ] Compare ML predictions vs Phred scores:
    - [ ] Agreement (correlation, MAE)
    - [ ] Disagreement (where ML differs)
    - [ ] Error analysis (systematic biases?)
  - [ ] Filtering comparison:
    - [ ] Quality filter using Phred
    - [ ] Quality filter using ML predictions
    - [ ] Which is more accurate? (validate against alignment)

- [ ] **Documentation** (2-3 hours):
  - [ ] Jupyter notebook: "Neural Engine Quality Prediction"
  - [ ] Usage examples (Rust + Python)
  - [ ] Performance characteristics
  - [ ] Comparison: Phred vs ML

**Deliverables**:
- Rust integration (CoreML bridge)
- Streaming quality prediction API
- Benchmarks (latency, throughput, energy)
- Jupyter notebook tutorial

**Success Criteria**:
- [ ] Neural Engine achieves <1ms latency per sequence
- [ ] Throughput: >1000 sequences/sec
- [ ] Energy: <1W (300× less than GPU cluster)
- [ ] ML predictions: Correlation >0.7 with Phred

**Blog Post** (if time):
- Title: "Neural Engine-Accelerated Genomics: Real-Time Quality Prediction"
- Focus: Energy efficiency, laptop-level ML inference

---

#### **Week 7: Alignment Primitives (CPU/NEON)** ⭐⭐⭐

**Status**: 🔴 Not Started
**Priority**: HIGH (Core primitives - enables tool building)
**Estimated Effort**: 25-30 hours
**Dependencies**: Smith-Waterman GPU (Week 1-2) provides foundation

**Tasks**:
- [ ] **Module Structure** (2-3 hours):
  - [ ] Create `src/alignment/` (if not exists from Week 1)
  - [ ] Module organization:
    ```
    src/alignment/
    ├── mod.rs              # Public API
    ├── smith_waterman.rs   # Local alignment
    ├── needleman_wunsch.rs # Global alignment
    ├── banded.rs           # Banded alignment
    ├── scoring.rs          # Scoring matrices
    ├── cigar.rs            # CIGAR operations
    └── neon.rs             # NEON optimizations
    ```

- [ ] **Smith-Waterman** (6-8 hours):
  - [ ] Refactor Week 1-2 code into clean API
  - [ ] CPU-only interface (no GPU dependency)
  - [ ] NEON-optimized striped algorithm
  - [ ] Configurable scoring (match, mismatch, gap)
  - [ ] Affine gap penalties
  - [ ] Traceback for CIGAR

- [ ] **Needleman-Wunsch** (6-8 hours):
  - [ ] Implement global alignment
  - [ ] Similar structure to Smith-Waterman
  - [ ] NEON optimization
  - [ ] Benchmarks (N=30)

- [ ] **Banded Alignment** (5-6 hours):
  - [ ] Implement constrained diagonal band
  - [ ] For similar sequences (faster, O(kn) vs O(mn))
  - [ ] Configurable band width
  - [ ] NEON optimization

- [ ] **Testing** (4-5 hours):
  - [ ] Unit tests: Known sequence pairs
  - [ ] Property tests: Score consistency
  - [ ] Edge cases: Empty, identical, no match
  - [ ] Integration tests: Various scoring matrices

- [ ] **Documentation** (2-3 hours):
  - [ ] API documentation with examples
  - [ ] Algorithm explanations
  - [ ] Usage patterns (how to build aligner)
  - [ ] Performance characteristics

**Deliverables**:
- Smith-Waterman primitive (local alignment)
- Needleman-Wunsch primitive (global alignment)
- Banded alignment primitive
- Comprehensive tests (50+)
- Documentation

**Success Criteria**:
- [ ] All alignment primitives working
- [ ] Tests passing (100%)
- [ ] Documentation complete
- [ ] Example: "Build aligner in 100 lines using biometal"

---

#### **Week 8: Variant Calling Primitives** ⭐⭐⭐

**Status**: 🔴 Not Started
**Priority**: HIGH (Core primitives - enables tool building)
**Estimated Effort**: 25-30 hours
**Dependencies**: Pileup GPU (Week 3) provides foundation, BAM parser (✅ v1.7.0)

**Tasks**:
- [ ] **Module Structure** (2-3 hours):
  - [ ] Create `src/variant/` module:
    ```
    src/variant/
    ├── mod.rs          # Public API
    ├── pileup.rs       # Pileup generation (CPU/GPU)
    ├── vcf.rs          # VCF parsing
    ├── bcf.rs          # BCF parsing (binary VCF)
    ├── vcf_writer.rs   # VCF writing
    ├── statistics.rs   # Statistical models
    └── filters.rs      # Quality/MAPQ filtering
    ```

- [ ] **Pileup Primitives** (5-6 hours):
  - [ ] Refactor Week 3 pileup into clean API
  - [ ] CPU-only interface (no GPU dependency for portability)
  - [ ] Per-position statistics:
    - [ ] Depth (read count)
    - [ ] Base counts (A, C, G, T, N)
    - [ ] Quality distribution
    - [ ] Strand bias
    - [ ] MAPQ distribution
  - [ ] Streaming interface (constant memory)

- [ ] **VCF/BCF Parsing** (8-10 hours):
  - [ ] VCF format parser (text):
    - [ ] Header parsing (##contig, ##FORMAT, etc.)
    - [ ] Record parsing (CHROM, POS, REF, ALT, QUAL, FILTER, INFO, FORMAT)
    - [ ] Genotype parsing (GT, GQ, DP, etc.)
    - [ ] Streaming parser (constant memory)
  - [ ] BCF format parser (binary):
    - [ ] Binary VCF (compressed)
    - [ ] Lazy parsing (don't decompress unless needed)
  - [ ] VCF writer:
    - [ ] Header generation
    - [ ] Record writing
    - [ ] Bgzip compression

- [ ] **Statistical Models** (6-8 hours):
  - [ ] Binomial test (variant detection):
    ```rust
    pub fn binomial_test(
        ref_count: u32,
        alt_count: u32,
        error_rate: f64,
    ) -> f64; // p-value
    ```
  - [ ] Beta-binomial test (overdispersion):
    ```rust
    pub fn beta_binomial_test(
        ref_count: u32,
        alt_count: u32,
        alpha: f64,
        beta: f64,
    ) -> f64;
    ```
  - [ ] Quality score recalibration
  - [ ] Allele frequency estimation

- [ ] **Testing** (4-5 hours):
  - [ ] Unit tests: Known VCF files
  - [ ] Property tests: Round-trip (read→write→read)
  - [ ] Integration tests: Pileup → variant calling
  - [ ] Compare: biometal VCF vs samtools/bcftools

- [ ] **Documentation** (2-3 hours):
  - [ ] API documentation
  - [ ] Usage examples: "Build variant caller in 200 lines"
  - [ ] Statistical model explanations
  - [ ] Performance characteristics

**Deliverables**:
- Pileup generation primitive (CPU/GPU)
- VCF/BCF parsing and writing
- Statistical models (binomial, beta-binomial)
- Comprehensive tests (50+)
- Documentation + usage examples

**Success Criteria**:
- [ ] VCF parsing: biometal matches bcftools (correctness)
- [ ] Pileup: biometal matches samtools mpileup
- [ ] Statistical models: p-values validated
- [ ] Example: "Build variant caller using biometal primitives"

---

## Phase 2: ML/BERT Integration (Months 3-4, Weeks 9-16)

**Objective**: Build world's first streaming ML data loaders + quality-aware tokenization for genomics

**Key Innovation**: No other library provides genomics-aware ML preprocessing with streaming

---

### Month 3: Core ML Primitives (Weeks 9-12)

#### **Week 9: Streaming BERT Data Loaders** ⭐⭐⭐

**Status**: 🔴 Not Started
**Priority**: HIGH (NOVEL - high publication value)
**Estimated Effort**: 25-30 hours
**Dependencies**: FASTQ parser (✅ v1.7.0), k-mer operations (✅ v1.0.0)

**Why this week**:
- Current ML pipelines load entire datasets (OOM on TB-scale)
- Novel: Streaming data loaders with constant 5 MB memory
- NO OTHER LIBRARY DOES THIS

**Tasks**:
- [ ] **Module Structure** (2-3 hours):
  - [ ] Create `src/ml/` module:
    ```
    src/ml/
    ├── mod.rs              # Public API
    ├── tokenizer.rs        # K-mer tokenization
    ├── vocab.rs            # Vocabulary (k-mer → token ID)
    ├── data_loader.rs      # Streaming data loader
    ├── quality_aware.rs    # Quality-aware tokenization (Week 10)
    ├── multimodal.rs       # Multi-modal inputs (Week 13)
    └── coreml.rs           # CoreML integration (Week 5-6)
    ```

- [ ] **K-mer Tokenizer** (6-8 hours):
  - [ ] Sliding window k-mer extraction:
    ```rust
    pub struct DnaBertTokenizer {
        k: usize,
        stride: usize,
        vocab: Vocabulary,
    }

    impl DnaBertTokenizer {
        pub fn tokenize(&self, sequence: &[u8]) -> Vec<u32> {
            // Extract overlapping k-mers
            // Map to vocabulary IDs
        }
    }
    ```
  - [ ] Vocabulary generation (all k-mers of size k)
  - [ ] Special tokens: [PAD], [CLS], [SEP], [MASK], [UNK]
  - [ ] Reverse complement handling

- [ ] **Streaming Data Loader** (10-12 hours):
  - [ ] Iterator-based streaming:
    ```rust
    pub struct DnaBertDataLoader<R: BufRead> {
        fastq_stream: FastqStream<R>,
        tokenizer: DnaBertTokenizer,
        max_length: usize,
        batch_size: usize,
        buffer: Vec<TokenizedSequence>,
    }

    impl<R: BufRead> Iterator for DnaBertDataLoader<R> {
        type Item = Result<BertBatch>;
        // Yields batches of tokenized sequences
    }
    ```
  - [ ] Batching (configurable batch size)
  - [ ] Padding (sequences → same length)
  - [ ] Attention masks (1 for real tokens, 0 for padding)
  - [ ] Memory: Constant ~5 MB (reuse buffers)

- [ ] **PyTorch/TensorFlow Integration** (4-5 hours):
  - [ ] Python bindings (PyO3):
    ```python
    import biometal

    dataset = biometal.ml.DnaBertDataset(
        fastq_path="5TB_metagenome.fq.gz",
        k=6,
        max_length=512,
        stride=256,
    )

    # PyTorch DataLoader (just works!)
    from torch.utils.data import DataLoader
    dataloader = DataLoader(dataset, batch_size=32, num_workers=4)

    for batch in dataloader:
        loss = model(batch.input_ids, batch.attention_mask, batch.labels)
        loss.backward()
    ```
  - [ ] TensorFlow tf.data integration
  - [ ] Test: Load 1M sequences (constant memory)

- [ ] **Testing** (2-3 hours):
  - [ ] Unit tests: Tokenization correctness
  - [ ] Memory tests: Constant 5 MB (1K, 10K, 100K, 1M sequences)
  - [ ] Integration tests: PyTorch/TensorFlow

- [ ] **Documentation** (2-3 hours):
  - [ ] API documentation
  - [ ] Usage examples (Rust + Python)
  - [ ] Jupyter notebook: "Training DNABERT on TB-scale Data"

**Deliverables**:
- K-mer tokenizer
- Streaming data loader (constant memory)
- PyTorch/TensorFlow integration
- Documentation + Jupyter notebook

**Success Criteria**:
- [ ] Memory constant at 5 MB (regardless of dataset size)
- [ ] Compatible with PyTorch DataLoader
- [ ] Throughput: >10K sequences/sec
- [ ] Example: Train DNABERT on 1TB using laptop

---

#### **Week 10: Quality-Aware Tokenization** ⭐⭐⭐

**Status**: 🔴 Not Started
**Priority**: HIGH (NOVEL - no one has done this!)
**Estimated Effort**: 20-25 hours
**Dependencies**: Week 9 complete

**Why this is HUGE**:
- Current BERT models ignore quality scores (treat all bases equally - WRONG!)
- Quality-aware tokenization: NO ONE HAS DONE THIS
- Research opportunity: "Quality-Aware Genomic Language Models" paper (Nature Methods?)

**Tasks**:
- [ ] **Design** (3-4 hours):
  - [ ] Research: How to incorporate quality into tokenization?
  - [ ] Approach 1: Masking (low-quality → [MASK])
  - [ ] Approach 2: Weighting (embed quality in representation)
  - [ ] Approach 3: Probabilistic (multiple token candidates)
  - [ ] Choose approach (or combine)

- [ ] **Implementation** (10-12 hours):
  - [ ] Quality masking:
    ```rust
    pub fn quality_aware_tokenize(
        sequence: &[u8],
        quality: &[u8],
        k: usize,
        quality_threshold: u8,
    ) -> QualityAwareTokens {
        let mut tokens = Vec::new();
        for kmer in sequence.windows(k) {
            let min_quality = kmer.iter()
                .zip(quality.iter())
                .map(|(_, &q)| q)
                .min()
                .unwrap_or(0);

            if min_quality < quality_threshold {
                tokens.push(MASK_TOKEN);  // Mask low-quality k-mers
            } else {
                tokens.push(vocab.get(kmer));
            }
        }
        tokens
    }
    ```
  - [ ] Quality weighting:
    ```rust
    pub struct WeightedToken {
        token_id: u32,
        quality_weight: f32,  // 0.0-1.0
    }
    ```
  - [ ] Probabilistic tokenization:
    ```rust
    // Generate multiple token candidates for ambiguous bases
    pub fn probabilistic_tokenize(...) -> Vec<Vec<(u32, f32)>> {
        // Returns: token_id → probability
    }
    ```

- [ ] **Integration** (3-4 hours):
  - [ ] Update `DnaBertDataLoader` to support quality-aware mode
  - [ ] Python API:
    ```python
    dataset = biometal.ml.DnaBertDataset(
        fastq_path="data.fq.gz",
        k=6,
        quality_aware=True,  # NEW!
        quality_threshold=20,
    )
    ```

- [ ] **Benchmarking** (3-4 hours):
  - [ ] Compare: Standard vs quality-aware tokenization
  - [ ] Masked bases: How many? (at different thresholds)
  - [ ] Training: Does quality-aware improve model? (need to train model)

- [ ] **Testing** (2-3 hours):
  - [ ] Unit tests: Low-quality → [MASK]
  - [ ] Property tests: High-quality preserved
  - [ ] Integration tests: End-to-end pipeline

- [ ] **Documentation** (2-3 hours):
  - [ ] API documentation
  - [ ] Jupyter notebook: "Quality-Aware BERT Tokenization"
  - [ ] Paper outline: "Quality-Aware Genomic Language Models"

**Deliverables**:
- Quality-aware tokenization (masking + weighting + probabilistic)
- Integration with data loader
- Benchmarks (standard vs quality-aware)
- Paper outline (for publication)

**Success Criteria**:
- [ ] Low-quality bases correctly masked
- [ ] Quality weighting integrated
- [ ] Memory still constant at 5 MB
- [ ] Paper outline: Ready for training experiments

**Research Impact**:
- This is PUBLISHABLE (Nature Methods, Bioinformatics, PLOS Comp Bio)
- NO ONE has done quality-aware genomic language models
- Could improve model accuracy 5-10%

---

#### **Week 11: GPU-Accelerated Tokenization** ⭐⭐

**Status**: 🔴 Not Started
**Priority**: MEDIUM (Apple Silicon advantage, but not essential)
**Estimated Effort**: 20-25 hours
**Dependencies**: Week 9 complete, GPU results (Week 4)

**Tasks**:
- [ ] **Metal Compute Shader** (10-12 hours):
  - [ ] Parallel k-mer extraction:
    ```metal
    kernel void extract_kmers_gpu(
        device const uint8_t* sequences [[buffer(0)]],
        device const uint32_t* vocab [[buffer(1)]],
        device uint32_t* tokens [[buffer(2)]],
        constant uint& k [[buffer(3)]],
        uint gid [[thread_position_in_grid]]
    ) {
        // Each GPU thread tokenizes one sequence
        // Parallel k-mer extraction + vocab lookup
    }
    ```
  - [ ] Vocabulary lookup on GPU (hash table or linear search)
  - [ ] Batch processing (tokenize 1000s of sequences in parallel)

- [ ] **Integration** (4-5 hours):
  - [ ] Rust wrapper for Metal shader
  - [ ] Unified memory (zero-copy sequences)
  - [ ] Fallback to CPU tokenization (non-Mac platforms)

- [ ] **Benchmarking** (4-5 hours):
  - [ ] Compare: CPU vs GPU tokenization (N=30)
  - [ ] Throughput: Sequences/sec
  - [ ] Latency: ms per sequence
  - [ ] Batch size effects (1, 10, 100, 1000, 10000)

- [ ] **Testing** (2-3 hours):
  - [ ] Property test: GPU tokens == CPU tokens
  - [ ] Edge cases

**Deliverables**:
- GPU-accelerated tokenization (Metal)
- Benchmarks (CPU vs GPU)
- Integration with data loader

**Success Criteria**:
- [ ] GPU achieves 10-100× speedup vs CPU tokenization
- [ ] Zero-copy advantage demonstrated (vs CUDA PCIe bottleneck)

---

#### **Week 12: Example Models + Documentation** ⭐⭐⭐

**Status**: 🔴 Not Started
**Priority**: HIGH (Deliverables for community)
**Estimated Effort**: 20-25 hours
**Dependencies**: Weeks 9-11 complete

**Tasks**:
- [ ] **Convert DNABERT to CoreML** (5-6 hours):
  - [ ] Download pretrained DNABERT from Hugging Face
  - [ ] Convert to CoreML (for Neural Engine):
    ```python
    from transformers import AutoModel
    import coremltools as ct

    model = AutoModel.from_pretrained("zhihan1996/DNA_bert_6")
    traced = torch.jit.trace(model, example_input)
    coreml_model = ct.convert(traced, convert_to="neuralnetwork")
    coreml_model.save("dnabert_6mer.mlmodel")
    ```
  - [ ] Test: CoreML predictions match PyTorch

- [ ] **Train Promoter Prediction Model** (5-6 hours):
  - [ ] Dataset: Promoter vs non-promoter sequences
  - [ ] Fine-tune DNABERT for classification
  - [ ] Evaluate: AUROC, AUPRC
  - [ ] Convert to CoreML

- [ ] **Train Quality Prediction Model** (Already done Week 5-6):
  - [ ] Document usage with data loaders

- [ ] **Jupyter Notebooks** (6-8 hours):
  - [ ] Notebook 1: "Streaming BERT Data Loaders" (Week 9)
  - [ ] Notebook 2: "Quality-Aware Tokenization" (Week 10)
  - [ ] Notebook 3: "GPU-Accelerated Tokenization" (Week 11)
  - [ ] Notebook 4: "Training DNABERT on TB-scale Data"
  - [ ] Each: 30-45 minutes hands-on, executable code

- [ ] **Blog Post** (4-5 hours):
  - [ ] Title: "From FASTQ to BERT on Your Laptop"
  - [ ] Sections:
    - [ ] Problem: TB-scale data doesn't fit in memory
    - [ ] Solution: Streaming data loaders (constant 5 MB)
    - [ ] Innovation: Quality-aware tokenization
    - [ ] Results: Train DNABERT on 5TB using MacBook
    - [ ] Code examples
  - [ ] Publish: Blog, Reddit, Twitter, Biostars

**Deliverables**:
- DNABERT CoreML model (Neural Engine)
- Promoter prediction model
- 4 Jupyter notebooks
- Blog post published

**Success Criteria**:
- [ ] Models working on Neural Engine
- [ ] Notebooks executable and clear
- [ ] Blog post published (500+ views)

---

### Month 4: Advanced ML Features (Weeks 13-16)

#### **Week 13: Multi-Modal Data Loaders** ⭐⭐⭐

**Status**: 🔴 Not Started
**Priority**: HIGH (NOVEL - research opportunity)
**Estimated Effort**: 25-30 hours
**Dependencies**: Week 3 (pileup GPU), Week 8 (VCF primitives), Week 9 (data loaders)

**Why this matters**:
- BERT models could benefit from richer context than sequence alone
- Novel combinations: Sequence + coverage + variants + quality
- NO OTHER LIBRARY provides unified multi-modal loading

**Tasks**:
- [ ] **Design** (4-5 hours):
  - [ ] Identify multi-modal use cases:
    - [ ] Sequence + coverage → variant calling
    - [ ] Sequence + ChIP-seq → TF binding prediction
    - [ ] Sequence + RNA-seq → expression prediction
    - [ ] Sequence + variants + quality → variant effect prediction
  - [ ] Design embedding strategies:
    - [ ] Coverage: Normalize, bin, embed
    - [ ] Variants: Binary mask or count
    - [ ] Quality: Distribution embedding

- [ ] **Implementation** (12-15 hours):
  - [ ] `MultiModalBamLoader`:
    ```rust
    pub struct MultiModalBamLoader {
        bam: BamReader,
        bai: BaiIndex,
        vcf: Option<VcfReader>,
        reference: Option<ReferenceGenome>,
    }

    impl MultiModalBamLoader {
        pub fn load_region_with_context(
            &self,
            chrom: &str,
            start: i32,
            end: i32,
        ) -> Result<MultiModalInput> {
            // 1. Sequence from reference
            let sequence = self.reference.fetch(chrom, start, end)?;

            // 2. Coverage from BAM (GPU-accelerated!)
            let coverage = pileup_gpu(&self.bam, &self.bai, chrom, start, end)?;

            // 3. Variants from VCF
            let variants = self.vcf.query(chrom, start, end)?;

            // 4. Quality distribution
            let quality_dist = quality_distribution(&self.bam, chrom, start, end)?;

            Ok(MultiModalInput {
                sequence_tokens: tokenize(&sequence, k=6),
                coverage_embedding: coverage.normalize().embed(),
                variant_mask: variants.to_binary_mask(),
                quality_embedding: quality_dist.embed(),
            })
        }
    }
    ```
  - [ ] Embedding functions:
    - [ ] Coverage normalization (log-scale, quantile)
    - [ ] Variant mask (binary or count)
    - [ ] Quality distribution (histogram)
  - [ ] Integration with data loader (streaming)

- [ ] **Python API** (4-5 hours):
  - [ ] PyO3 bindings:
    ```python
    import biometal

    loader = biometal.ml.MultiModalBamLoader(
        bam_path="alignments.bam",
        vcf_path="variants.vcf.gz",
        reference_path="hg38.fa",
    )

    for region in regions:
        multimodal_input = loader.load_region(region.chrom, region.start, region.end)
        # multimodal_input.sequence_tokens
        # multimodal_input.coverage_embedding
        # multimodal_input.variant_mask
        # multimodal_input.quality_embedding
    ```

- [ ] **Example Use Case** (3-4 hours):
  - [ ] Variant effect prediction:
    - [ ] Load region around variant
    - [ ] Multi-modal input
    - [ ] BERT prediction (benign vs pathogenic)

- [ ] **Testing** (2-3 hours):
  - [ ] Unit tests: Each component
  - [ ] Integration tests: End-to-end
  - [ ] Memory: Constant 5 MB

- [ ] **Documentation** (2-3 hours):
  - [ ] API documentation
  - [ ] Jupyter notebook: "Multi-Modal Genomic Models"

**Deliverables**:
- `MultiModalBamLoader` (sequence + coverage + variants + quality)
- Python API
- Example: Variant effect prediction
- Jupyter notebook

**Success Criteria**:
- [ ] Multi-modal loading working
- [ ] Integration with GPU pileup (fast coverage)
- [ ] Memory constant at 5 MB
- [ ] Example model trained and evaluated

**Research Impact**:
- Novel approach (no existing work on multi-modal genomic BERT)
- Publishable: "Multi-Modal Transformers for Variant Effect Prediction"

---

#### **Week 14: Neural Engine Adapter Detection** ⭐⭐

**Status**: 🔴 Not Started
**Priority**: MEDIUM (Novel, but lower impact than quality-aware)
**Estimated Effort**: 20-25 hours
**Dependencies**: Week 5-6 (Neural Engine integration)

**Tasks**:
- [ ] **Dataset Preparation** (5-6 hours):
  - [ ] Adapter sequences (Illumina, Nextera, etc.)
  - [ ] Generate training data:
    - [ ] Positive: Reads with adapters
    - [ ] Negative: Reads without adapters
    - [ ] Augmentation: Different positions

- [ ] **Model Training** (6-8 hours):
  - [ ] Classification model (adapter yes/no)
  - [ ] PyTorch → CoreML

- [ ] **Integration** (5-6 hours):
  - [ ] Real-time detection during streaming:
    ```rust
    pub fn detect_adapters_streaming(
        fastq_path: &str,
        model_path: &str,
    ) -> impl Iterator<Item = Result<(FastqRecord, bool)>> {
        // Real-time adapter detection using Neural Engine
    }
    ```
  - [ ] Automatic trimming

- [ ] **Benchmarking** (3-4 hours):
  - [ ] Compare vs k-mer matching (Cutadapt)
  - [ ] Accuracy: Precision, recall, F1
  - [ ] Speed: Sequences/sec

- [ ] **Documentation** (2-3 hours):
  - [ ] API documentation
  - [ ] Usage examples

**Deliverables**:
- Adapter detection model (Neural Engine)
- Real-time detection API
- Benchmarks vs Cutadapt

**Success Criteria**:
- [ ] Accuracy: F1 >0.95
- [ ] Real-time: <1ms per sequence (Neural Engine)

---

#### **Week 15: Training Utilities** ⭐⭐

**Status**: 🔴 Not Started
**Priority**: MEDIUM
**Estimated Effort**: 15-20 hours

**Tasks**:
- [ ] **Model Fine-Tuning Helpers** (5-6 hours):
  - [ ] Load pretrained DNABERT
  - [ ] Fine-tune on custom task
  - [ ] Save/load checkpoints

- [ ] **Loss Functions** (4-5 hours):
  - [ ] Binary classification (promoter yes/no)
  - [ ] Regression (quality prediction)
  - [ ] Sequence-to-sequence (e.g., translation)

- [ ] **Data Augmentation** (3-4 hours):
  - [ ] Reverse complement
  - [ ] Random masking (MLM)
  - [ ] Sequence perturbation

- [ ] **Evaluation Metrics** (3-4 hours):
  - [ ] AUROC, AUPRC
  - [ ] Precision, recall, F1
  - [ ] MAE, RMSE, R²

**Deliverables**:
- Training utilities
- Documentation

---

#### **Week 16: Format Support (BED, GFF/GTF)** ⭐⭐⭐

**Status**: 🔴 Not Started
**Priority**: HIGH (Core primitives)
**Estimated Effort**: 20-25 hours

**Tasks**:
- [ ] **BED Format** (8-10 hours):
  - [ ] Streaming parser
  - [ ] Interval operations (overlap, merge, subtract)
  - [ ] Integration with BAI queries:
    ```rust
    // Filter BAM by BED regions
    let regions = BedStream::from_path("exons.bed")?;
    for region in regions {
        for record in bam.query(&index, &region)? {
            // Process reads in exons only
        }
    }
    ```

- [ ] **GFF/GTF Format** (8-10 hours):
  - [ ] Streaming parser
  - [ ] Feature hierarchy (gene → transcript → exon)
  - [ ] Attribute parsing (gene_id, transcript_id)
  - [ ] Integration with RNA-seq workflows

- [ ] **Testing** (3-4 hours):
  - [ ] Unit tests
  - [ ] Integration tests

- [ ] **Documentation** (2-3 hours):
  - [ ] API documentation
  - [ ] Usage examples

**Deliverables**:
- BED parser (streaming)
- GFF/GTF parser (streaming)
- Interval operations
- Documentation

**Success Criteria**:
- [ ] BED/GFF parsing matches bedtools/gffread
- [ ] Interval operations correct
- [ ] Integration with BAM queries working

---

## Phase 3: More Primitives + Demonstration (Months 5-6, Weeks 17-24)

### Month 5: Additional Primitives (Weeks 17-20)

#### **Week 17: Variant Calling GPU (Statistical Models)** ⭐⭐

**Status**: 🔴 Not Started
**Priority**: MEDIUM (Build on Week 3 pileup)
**Estimated Effort**: 20-25 hours
**Dependencies**: Week 3 (pileup GPU), Week 8 (statistical models)

**Tasks**:
- [ ] **GPU Statistical Inference** (12-15 hours):
  - [ ] Parallel binomial tests (one thread per position)
  - [ ] Parallel beta-binomial tests
  - [ ] Metal compute shader
  - [ ] Integration with pileup GPU

- [ ] **Benchmarking** (5-6 hours):
  - [ ] Compare: GPU vs CPU statistical inference
  - [ ] Throughput: Positions/sec
  - [ ] Target: 10-20× speedup

- [ ] **Testing** (3-4 hours):
  - [ ] Validate: GPU p-values == CPU p-values

**Deliverables**:
- GPU-accelerated variant calling
- Benchmarks

**Success Criteria**:
- [ ] GPU achieves 10-20× speedup
- [ ] Correctness validated

---

#### **Week 18-19: Assembly Primitives** ⭐⭐

**Status**: 🔴 Not Started
**Priority**: MEDIUM (Enables assembler building)
**Estimated Effort**: 40-50 hours (2 weeks)

**Tasks**:
- [ ] **De Bruijn Graph** (15-20 hours):
  - [ ] Graph construction from k-mers
  - [ ] Node/edge representation
  - [ ] Memory-efficient storage

- [ ] **Graph Operations** (15-20 hours):
  - [ ] Eulerian path traversal
  - [ ] Tip clipping (remove dead ends)
  - [ ] Bubble popping (merge variants)
  - [ ] Contig extraction

- [ ] **Overlap Detection** (8-10 hours):
  - [ ] Suffix-prefix matching (OLC)
  - [ ] Overlap graph construction

- [ ] **Testing** (4-5 hours):
  - [ ] Unit tests: Known assemblies
  - [ ] Integration tests: Real data

- [ ] **Documentation** (3-4 hours):
  - [ ] API documentation
  - [ ] Usage: "Build assembler using biometal"

**Deliverables**:
- De Bruijn graph construction
- Graph operations (Eulerian path, tip clipping, bubble popping)
- Overlap detection (OLC)
- Documentation

**Success Criteria**:
- [ ] Graph construction working
- [ ] Assembly primitives functional
- [ ] Example: "Build assembler in 300 lines"

---

#### **Week 20: AMX RNA-seq Matrices** ⭐⭐

**Status**: 🔴 Not Started
**Priority**: EXPLORATORY (May fail, document negative result)
**Estimated Effort**: 15-20 hours

**Why try (despite ASBB failure)**:
- ASBB tested small matrices (<100×100)
- RNA-seq has LARGE matrices (20K × 500 = 10M elements)
- May amortize overhead

**Tasks**:
- [ ] **Implementation** (8-10 hours):
  - [ ] Matrix operations using Accelerate framework (AMX)
  - [ ] Count matrix normalization (log, CPM, TPM)
  - [ ] PCA/dimensionality reduction
  - [ ] Distance matrices

- [ ] **Benchmarking** (5-6 hours):
  - [ ] Compare: AMX vs NEON vs scalar
  - [ ] Matrix sizes: 1K×1K, 10K×10K, 20K×500
  - [ ] Operations: Multiplication, normalization, PCA

- [ ] **Documentation** (2-3 hours):
  - [ ] If AMX helps: Document usage
  - [ ] If AMX fails: Document negative result (valuable!)

**Deliverables**:
- AMX RNA-seq operations (or negative result documentation)

**Success Criteria**:
- [ ] If AMX >2× faster: Success, document usage
- [ ] If AMX <2× faster: Negative result, document why (overhead still dominates)

---

### Month 6: Demonstration & Publication (Weeks 21-24)

#### **Week 21: Case Study 1 - "500GB RNA-seq on MacBook Air"** ⭐⭐⭐

**Status**: 🔴 Not Started
**Priority**: HIGH (Proof of democratization vision)
**Estimated Effort**: 20-25 hours

**Tasks**:
- [ ] **Dataset** (3-4 hours):
  - [ ] 50 RNA-seq samples × 10GB each = 500GB
  - [ ] Public data (SRA) or simulated

- [ ] **Analysis Pipeline** (8-10 hours):
  - [ ] Streaming FASTQ parser (constant 5 MB memory)
  - [ ] Read counting (gene expression)
  - [ ] Count matrix generation (20K genes × 50 samples)
  - [ ] Normalization (CPM, TPM)
  - [ ] PCA/clustering
  - [ ] AMX-accelerated (if Week 20 successful)

- [ ] **Comparative Benchmarking** (6-8 hours):
  - [ ] Mac Studio ($5K): Time, energy, cost
  - [ ] HPC cluster ($50K+): Time, energy, cost, queue wait
  - [ ] Traditional tools (Salmon, STAR): Memory usage

- [ ] **Documentation** (3-4 hours):
  - [ ] Blog post: "500GB RNA-seq on a MacBook Air"
  - [ ] Jupyter notebook (reproducible)
  - [ ] Metrics: Time, energy, cost, memory

**Deliverables**:
- Complete RNA-seq pipeline (streaming)
- Comparative analysis (Mac vs HPC)
- Blog post
- Jupyter notebook

**Success Criteria**:
- [ ] Memory constant at 5 MB (not 500GB)
- [ ] Analysis completes on MacBook Air
- [ ] Mac competitive with HPC (time/cost)

---

#### **Week 22: Case Study 2 - "GPU-Accelerated Variant Calling"** ⭐⭐⭐

**Status**: 🔴 Not Started
**Priority**: HIGH (If GPU successful in Week 1-4)
**Estimated Effort**: 20-25 hours
**Dependencies**: Weeks 1-4 (GPU), Week 8 (variant calling)

**Tasks**:
- [ ] **Dataset** (2-3 hours):
  - [ ] 50× coverage human genome (150GB BAM)
  - [ ] Public data (1000 Genomes) or simulated

- [ ] **Variant Calling Pipeline** (10-12 hours):
  - [ ] Streaming BAM parser
  - [ ] GPU-accelerated pileup generation (Week 3)
  - [ ] GPU-accelerated statistical inference (Week 17)
  - [ ] VCF output

- [ ] **Comparative Benchmarking** (6-8 hours):
  - [ ] Mac Studio + GPU: Time, energy
  - [ ] HPC cluster CPU-only: Time, energy
  - [ ] Commercial cloud GPU (CUDA): Time, cost
  - [ ] Compare vs samtools/GATK/DeepVariant

- [ ] **Documentation** (3-4 hours):
  - [ ] Blog post: "GPU-Accelerated Variant Calling on Apple Silicon"
  - [ ] Jupyter notebook
  - [ ] Metrics

**Deliverables**:
- Complete variant calling pipeline (GPU)
- Comparative analysis (Mac vs HPC vs cloud)
- Blog post
- Jupyter notebook

**Success Criteria**:
- [ ] GPU pileup 20-100× faster than CPU
- [ ] Variant calling competitive with GATK
- [ ] Memory constant at 5 MB

---

#### **Week 23: Case Study 3 - "Training DNABERT on 5TB Metagenome"** ⭐⭐⭐

**Status**: 🔴 Not Started
**Priority**: HIGH (Proof of ML streaming)
**Estimated Effort**: 20-25 hours
**Dependencies**: Weeks 9-12 (ML data loaders)

**Tasks**:
- [ ] **Dataset** (3-4 hours):
  - [ ] 1,000 metagenome samples from SRA (5TB if downloaded)
  - [ ] Don't download! Stream directly from NCBI

- [ ] **Training Pipeline** (10-12 hours):
  - [ ] Streaming data loader (constant 5 MB memory)
  - [ ] Quality-aware tokenization
  - [ ] Train DNABERT (or fine-tune pretrained)
  - [ ] Monitor: Memory usage, throughput

- [ ] **Comparative Benchmarking** (5-6 hours):
  - [ ] Traditional: Download 5TB + load all (FAIL on laptop)
  - [ ] biometal: Stream from network (constant 5 MB memory)
  - [ ] Time: Traditional vs biometal
  - [ ] Storage: 5TB vs <100MB

- [ ] **Documentation** (3-4 hours):
  - [ ] Blog post: "Training DNABERT on 5TB Without Downloading"
  - [ ] Jupyter notebook
  - [ ] Metrics

**Deliverables**:
- Training pipeline (streaming from SRA)
- Trained DNABERT model
- Comparative analysis
- Blog post
- Jupyter notebook

**Success Criteria**:
- [ ] Training completes on laptop (constant 5 MB memory)
- [ ] Never downloads 5TB
- [ ] Model performance competitive with traditional training

---

#### **Week 24: Publications & Community** ⭐⭐⭐

**Status**: 🔴 Not Started
**Priority**: HIGH (Dissemination)
**Estimated Effort**: 25-30 hours

**Tasks**:
- [ ] **Blog Posts** (8-10 hours):
  - [ ] Post 1: "GPU-Accelerated Bioinformatics on Apple Silicon" (Week 4)
  - [ ] Post 2: "Neural Engine for Genomics" (Week 6)
  - [ ] Post 3: "From FASTQ to BERT on Your Laptop" (Week 12)
  - [ ] Post 4: Case studies (Weeks 21-23)
  - [ ] Publish on:
    - [ ] Blog
    - [ ] Reddit (r/bioinformatics, r/MachineLearning, r/apple)
    - [ ] Twitter/X
    - [ ] Biostars
    - [ ] LinkedIn

- [ ] **Paper 1: biometal Architecture** (6-8 hours):
  - [ ] Title: "biometal: Apple Silicon-Native Bioinformatics with ML Integration"
  - [ ] Venue: Bioinformatics (Application Note) or JOSS
  - [ ] Sections:
    - [ ] Abstract (250 words)
    - [ ] Introduction (motivation, gap)
    - [ ] Methods (architecture, streaming, GPU, Neural Engine)
    - [ ] Results (benchmarks, case studies)
    - [ ] Discussion (democratization)
    - [ ] Availability (GitHub, PyPI)
  - [ ] Status: Draft complete

- [ ] **Paper 2: Quality-Aware Models** (6-8 hours):
  - [ ] Title: "Quality-Aware Genomic Language Models"
  - [ ] Venue: Nature Methods, Bioinformatics, or PLOS Comp Bio
  - [ ] Sections:
    - [ ] Abstract
    - [ ] Introduction (BERT models ignore quality - wrong!)
    - [ ] Methods (quality-aware tokenization, training)
    - [ ] Results (accuracy improvement, benchmarks)
    - [ ] Discussion (implications)
  - [ ] Status: Draft complete

- [ ] **Paper 3: Unified Memory** (4-5 hours):
  - [ ] Title: "Unified Memory Architectures for Bioinformatics ML: Apple Silicon Case Study"
  - [ ] Venue: ISMB, NeurIPS (Computational Biology track)
  - [ ] Status: Outline complete

- [ ] **Conference Submissions** (3-4 hours):
  - [ ] ISMB 2026 (abstract deadline: Feb 2026)
  - [ ] NeurIPS 2026 Comp Bio track (deadline: May 2026)
  - [ ] BOSC 2026 (abstract deadline: April 2026)

- [ ] **Community Outreach** (2-3 hours):
  - [ ] GitHub Discussions (announce papers, gather feedback)
  - [ ] Engage with researchers (reply to issues, questions)
  - [ ] Social media (promote case studies)

**Deliverables**:
- 4 blog posts published
- 3 paper drafts complete
- 2-3 conference submissions
- Community engagement

**Success Criteria**:
- [ ] Blog posts: 1000+ views total
- [ ] Papers: Submitted to journals
- [ ] Conferences: Abstracts accepted
- [ ] GitHub stars: 500+
- [ ] Community: 10+ users building with biometal

---

## Risk Management & Contingency Plans

### GPU Experiments (Weeks 1-4)

**Risk**: GPU doesn't achieve 10-50× speedup
**Probability**: 30-40%
**Impact**: High (affects Weeks 5-8 GPU work, Week 17 variant calling GPU)

**Contingency Plan**:
- **If GPU <2× speedup** ❌:
  - Document negative result (valuable for community)
  - Skip GPU work in Weeks 5-8, 17
  - Reallocate time:
    - Week 5-8: More ML primitives, more core primitives
    - Week 17: Skip GPU variant calling
  - Pivot: Focus on ML/BERT (higher novelty)

**Mitigation**:
- Target operations with complexity >0.70 (not 0.20-0.61 like ASBB)
- Smith-Waterman proven in CUDA literature
- Start simple, iterate

### Neural Engine (Weeks 5-6, 14)

**Risk**: Neural Engine too slow or complex to integrate
**Probability**: 20-30%
**Impact**: Medium (affects Neural Engine work)

**Contingency Plan**:
- **If Neural Engine fails** ❌:
  - Document negative result
  - Use CPU inference instead (still works, just slower)
  - Focus on streaming data loaders (Week 9-12) - higher impact

**Mitigation**:
- Use CoreML (Apple's official framework)
- Start with quality prediction (simpler than BERT)
- CPU fallback always available

### Quality-Aware Tokenization (Week 10)

**Risk**: Quality-aware models don't improve accuracy
**Probability**: 40-50%
**Impact**: Medium (affects paper novelty)

**Contingency Plan**:
- **If no accuracy improvement** ❌:
  - Still document approach (novel contribution)
  - Publish negative result (valuable for community)
  - Focus on other ML primitives

**Success**: Even if doesn't improve accuracy, still novel (no one has tried this!)

### AMX (Week 20)

**Risk**: AMX fails again (like ASBB)
**Probability**: 60-70%
**Impact**: Low (exploratory, not essential)

**Contingency Plan**:
- **If AMX fails** ❌:
  - Document negative result (confirms ASBB finding)
  - Use NEON/scalar instead
  - 1 week wasted, but valuable negative result

---

## Weekly Time Commitment

**Estimated hours per week**: 20-30 hours

**Breakdown**:
- Implementation: 10-15 hours
- Testing: 4-5 hours
- Benchmarking: 3-4 hours
- Documentation: 2-3 hours
- Research/planning: 2-3 hours

**Total for 24 weeks**: 480-720 hours (12-18 full-time weeks)

**Recommended pace**: 25 hours/week (sustainable, allows for other work)

---

## Dependencies Graph

```
Week 1-2: Smith-Waterman GPU
├── Week 3: Pileup GPU (independent)
│   └── Week 8: Variant Calling Primitives
│       └── Week 17: Variant Calling GPU
│       └── Week 13: Multi-Modal Loaders
├── Week 4: GPU Results (decision point)
└── Week 7: Alignment Primitives (CPU/NEON)

Week 5-6: Neural Engine Quality
├── Week 14: Neural Engine Adapter (builds on Week 5-6)
└── Week 12: Example Models (converts to CoreML)

Week 9: Streaming Data Loaders
├── Week 10: Quality-Aware Tokenization
├── Week 11: GPU Tokenization (depends on Week 4 GPU results)
├── Week 12: Example Models
└── Week 13: Multi-Modal Loaders
    └── Week 23: Case Study 3 (BERT training)

Week 18-19: Assembly Primitives (independent)

Week 20: AMX (independent, exploratory)

Week 21-24: Case Studies (depend on Weeks 1-20)
```

---

## Tracking Progress

**Update this document weekly**:
- [ ] Mark tasks complete (✅)
- [ ] Update status (🔴 → 🟡 → 🟢)
- [ ] Document blockers (⏸️)
- [ ] Update contingency plans (if needed)

**Weekly checklist**:
- [ ] Review last week's progress
- [ ] Update PROJECT_TODOS.md
- [ ] Update CHANGELOG.md
- [ ] Document decisions/lessons learned
- [ ] Plan next week's tasks

---

**Document Created**: November 13, 2025
**Last Updated**: November 13, 2025
**Status**: Complete - Ready for execution
**Next Action**: Begin Week 1 (Smith-Waterman GPU implementation)
