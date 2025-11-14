# Neural Engine Research Notes

**Date**: November 13, 2025
**Goal**: Explore Apple Neural Engine for bioinformatics acceleration
**Status**: Week 2 of Strategic Pivot - Initial Research

---

## Executive Summary

**Key Finding**: Neural Engine is completely untapped in bioinformatics - **zero existing work found**

**Opportunity**: First-mover advantage for ML-accelerated sequence analysis on Apple Silicon

**Technical Path**: Core ML + Rust bindings → DNA sequence classification/embedding

---

## Apple Neural Engine Architecture

### Performance Specifications

| Platform | Neural Engine Cores | TOPS | Process Node | Year |
|----------|---------------------|------|--------------|------|
| M1 | 16 cores | 11 TOPS | 5nm | 2020 |
| M2 | 16 cores | 15.8 TOPS | 5nm | 2022 |
| M3 | 16 cores | 18 TOPS | 3nm | 2023 |
| M4 | 16 cores | 38 TOPS | 3nm (2nd gen) | 2024 |
| A17 Pro | 16 cores | 35 TOPS | 3nm | 2023 |

**Mac Studio (M3 Ultra)** (2025): Can run 600B+ parameter LLMs entirely in memory

### Architecture Characteristics

1. **Dedicated ML Accelerator**
   - Separate from CPU and GPU
   - Optimized for matrix operations and neural network inference
   - Extremely power-efficient (vs GPU)

2. **Unified Memory Architecture**
   - Zero-copy access to system memory
   - No explicit data transfers needed (unlike CUDA)
   - Shared memory pool with CPU/GPU

3. **Core ML Integration**
   - Automatic hardware dispatch (CPU/GPU/Neural Engine)
   - Framework handles optimization decisions
   - Developer specifies model, not execution target

4. **Transformer Optimization**
   - Apple published specific optimizations for Transformer models
   - Efficient attention mechanisms
   - Quantization support (INT8, FP16)

---

## Core ML Framework

### Overview

Core ML is Apple's unified ML framework for on-device inference:
- Supports TensorFlow, PyTorch, ONNX model conversion
- Automatic hardware selection (CPU/GPU/Neural Engine)
- Optimization: quantization, pruning, model compilation
- Privacy-preserving (all inference on-device)

### Rust Integration Options

| Crate | Status | Approach | Maturity |
|-------|--------|----------|----------|
| `coreml-rs` | Experimental | Direct CoreML bindings | Early (WIP) |
| `objc2-core-ml` | Maintained | Objective-C bridge | Mature |
| `ort` (ONNX Runtime) | Production | CoreML execution provider | Production-ready ✅ |

**Recommendation**: Start with `ort` crate (ONNX Runtime with CoreML backend)
- Most mature and production-ready
- Standard ONNX model format
- Proven CoreML integration
- Active maintenance

### Model Deployment Flow

```
Train Model (PyTorch/TensorFlow)
    ↓
Export to ONNX format
    ↓
Load with ort crate in Rust
    ↓
ort automatically uses CoreML backend
    ↓
CoreML dispatches to Neural Engine
```

---

## DNA Sequence ML: State of the Art

### Foundation Models (2024)

1. **Nucleotide Transformer (NT)** - November 2024
   - 500M to 2.5B parameters
   - Pre-trained on human genome + 3,202 diverse genomes + 850 species
   - Uses k=6 tokenization (6-mers as "words")
   - Masked language modeling (MLM) objective
   - Applications: gene prediction, variant calling, regulatory element identification

2. **DNABERT** - 2020 (still widely used)
   - BERT architecture adapted for DNA
   - K-mer tokenization (k=3 to k=6)
   - Pre-trained on human genome
   - Fine-tunable for specific tasks

3. **Sentence Transformers for DNA** - January 2025 (very recent!)
   - Fine-tuned on DNA k-mers (k=3, n=3000 samples)
   - Compared against DNABERT and NT
   - Evaluated across 8 benchmark tasks

### K-mer Representation

**Standard approach**: Tokenize DNA as k-mers
- **Overlapping k-mers**: Better context but information leakage
- **Non-overlapping k-mers**: Efficient but may miss patterns
- **Common k values**: 3, 4, 5, 6 (tradeoff: vocabulary size vs context)

**Example** (k=3, overlapping):
```
DNA:  ACGTACGT
K-mers: ACG, CGT, GTA, TAC, ACG, CGT
```

### Applications in Bioinformatics

| Task | Model Type | Input | Output | Difficulty |
|------|------------|-------|--------|------------|
| Taxonomic classification | Transformer | Sequence | Species ID | Medium |
| Gene prediction | Transformer | Sequence | Gene regions | High |
| Quality prediction | CNN/MLP | Sequence + QC | Pass/Fail | **Low** ✅ |
| Basecalling | RNN/Transformer | Signal | Sequence | High |
| K-mer embedding | Embedding layer | K-mer | Vector | **Low** ✅ |
| Variant calling | Transformer | Alignment | Variants | High |

---

## Proposed Neural Engine Tasks (Prioritized)

### Option 1: Read Quality Prediction (RECOMMENDED)

**Goal**: Predict if a FASTQ read passes quality thresholds

**Input**:
- DNA sequence (encoded as integers: A=0, C=1, G=2, T=3)
- Quality scores (Phred scores)
- Sequence length

**Output**:
- Binary classification: PASS/FAIL
- Probability score

**Model Architecture**:
- Simple CNN or MLP
- Input: Concatenated sequence + quality (e.g., 300 features for 150bp read)
- Hidden layers: 2-3 layers, 64-128 neurons
- Output: Sigmoid activation → probability

**Why This is Perfect**:
- ✅ **Simple model**: Can train in hours, not days
- ✅ **Clear labels**: Ground truth from existing QC tools
- ✅ **Practical value**: Replace CPU-based quality filtering with Neural Engine
- ✅ **Measurable**: Compare accuracy vs current quality_filter_neon()
- ✅ **Small size**: Model <1MB, perfect for on-device inference

**Expected Performance**:
- Training: 10,000-100,000 FASTQ reads with quality labels
- Inference: Neural Engine should be 10-100× faster than CPU
- Accuracy target: >95% (match or beat rule-based filtering)

**Implementation Complexity**: LOW (achievable in Week 2)

---

### Option 2: K-mer Embedding

**Goal**: Learn dense vector representations of k-mers

**Input**: K-mer sequence (e.g., "ACG" → 0, "CGT" → 1, ...)

**Output**: Dense vector embedding (e.g., 64 or 128 dimensions)

**Model**: Simple embedding layer + optional dense layers

**Why Interesting**:
- ✅ Foundation for downstream tasks
- ✅ Can compare to DNABERT embeddings
- ✅ Useful for k-mer similarity, clustering, search

**Challenges**:
- ⚠️ Requires labeled data or self-supervised learning
- ⚠️ Less immediately practical than quality prediction
- ⚠️ Need to define learning objective

**Implementation Complexity**: MEDIUM

---

### Option 3: GC Content Prediction (PROOF OF CONCEPT)

**Goal**: Predict GC content from sequence

**Input**: DNA sequence (encoded)

**Output**: GC content percentage (regression)

**Why Interesting**:
- ✅ **Trivial ground truth**: Can generate infinite training data
- ✅ **Sanity check**: Verifies Neural Engine integration works
- ✅ **Baseline**: Compare Neural Engine vs current gc_content_neon()

**Why NOT Recommended**:
- ❌ Too simple - rule-based algorithm is already optimal
- ❌ No real value - ML is overkill for this task
- ❌ Doesn't demonstrate ML capabilities

**Use Case**: Initial integration test only

**Implementation Complexity**: VERY LOW

---

## Recommended Path Forward (Week 2)

### Day 1-2: Core ML + ONNX Integration

1. **Set up ONNX Runtime with CoreML**
   ```toml
   [dependencies]
   ort = { version = "2.0", features = ["coreml"] }
   ```

2. **Create minimal example**
   - Load pre-trained ONNX model
   - Run inference on dummy data
   - Verify Neural Engine execution

3. **Benchmark Neural Engine vs CPU**
   - Simple matrix multiplication model
   - Measure latency and throughput
   - Confirm Neural Engine engagement

### Day 3-4: Read Quality Prediction Model

1. **Data preparation**
   - Extract 10,000 FASTQ reads from existing test data
   - Label as PASS/FAIL using current quality_filter()
   - Encode sequences and quality scores as arrays

2. **Model training (Python/PyTorch)**
   ```python
   import torch
   import torch.nn as nn

   class ReadQualityClassifier(nn.Module):
       def __init__(self, seq_length=150):
           super().__init__()
           self.fc1 = nn.Linear(seq_length * 2, 128)  # seq + quality
           self.fc2 = nn.Linear(128, 64)
           self.fc3 = nn.Linear(64, 1)
           self.sigmoid = nn.Sigmoid()

       def forward(self, x):
           x = torch.relu(self.fc1(x))
           x = torch.relu(self.fc2(x))
           x = self.sigmoid(self.fc3(x))
           return x
   ```

3. **Export to ONNX**
   ```python
   torch.onnx.export(model, dummy_input, "read_quality.onnx")
   ```

4. **Integrate in Rust**
   - Load ONNX model with `ort`
   - Create inference function
   - Add to biometal's quality filtering pipeline

### Day 5: Benchmarking and Analysis

1. **Performance comparison**
   - Neural Engine inference time
   - CPU quality_filter_neon() time
   - GPU metal time (if applicable)

2. **Accuracy evaluation**
   - Confusion matrix vs rule-based filter
   - Precision, recall, F1 score
   - ROC curve analysis

3. **Documentation**
   - Neural Engine integration guide
   - Model training pipeline
   - Performance analysis report

---

## Open Questions

1. **Model Size Limits**: What's the maximum model size for Neural Engine?
   - Answer: Depends on available memory, but M4 can handle 600B+ parameters

2. **Batch Processing**: How does Neural Engine perform on batch inference?
   - Need to benchmark: single vs batch (10, 50, 100 reads)

3. **Quantization**: Does INT8 quantization improve Neural Engine performance?
   - Core ML supports FP16 and INT8
   - May provide 2-4× speedup with minimal accuracy loss

4. **Real-time Constraints**: Can Neural Engine handle streaming reads?
   - Need to measure latency: target <1ms per read

5. **Model Update Strategy**: How to distribute and update ML models?
   - Ship ONNX model with crate
   - Allow user-provided custom models
   - Model versioning strategy

---

## Success Criteria (Week 2)

### Minimum Viable Success ✅
- [ ] ONNX Runtime + CoreML integration working
- [ ] Simple model running on Neural Engine (verified)
- [ ] Benchmark showing Neural Engine speedup vs CPU

### Target Success 🎯
- [ ] Read quality prediction model trained (>90% accuracy)
- [ ] Integrated into biometal's quality filtering
- [ ] 10-100× speedup vs CPU demonstrated
- [ ] Documentation and examples complete

### Stretch Goals 🚀
- [ ] K-mer embedding model
- [ ] Comparison to DNABERT embeddings
- [ ] Multi-task model (quality + GC content + complexity)
- [ ] Model quantization analysis (FP32 vs FP16 vs INT8)

---

## Strategic Value

**Why Neural Engine Matters**:

1. **Power Efficiency**: 10-100× more efficient than GPU for ML inference
2. **Latency**: Optimized for real-time inference (<1ms)
3. **Untapped Resource**: Every Apple Silicon Mac has Neural Engine, almost never used
4. **Differentiation**: Zero existing bioinformatics tools use Neural Engine
5. **Scalability**: Enable ML-powered analysis on laptops, not just servers

**Impact on biometal**:
- First bioinformatics library to leverage Neural Engine
- Enables sophisticated ML models in production pipelines
- Foundation for Week 3-4: BERT integration for taxonomic classification
- Proves viability of on-device ML for genomics

---

## Next Steps

1. ✅ Complete this research document
2. ⏳ Choose specific task (RECOMMENDED: Read quality prediction)
3. ⏳ Set up ONNX Runtime + CoreML integration
4. ⏳ Create minimal working example
5. ⏳ Begin model training pipeline

**Decision Point**: Confirm read quality prediction as Week 2 focus before proceeding to implementation.
