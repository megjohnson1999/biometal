# Neural Engine Integration for biometal

This directory contains research, training scripts, and documentation for integrating Apple Neural Engine ML inference into biometal.

## Overview

**Goal**: Enable ML-powered bioinformatics operations using Apple's Neural Engine for sub-millisecond inference latency.

**Status**: Week 2 - Days 1-4 (ONNX Runtime integration + model training)

## Architecture

```
PyTorch Model Training → ONNX Export → ort crate → CoreML → Neural Engine
```

- **Training**: Python (PyTorch) - `train_quality_model.py`
- **Inference**: Rust (ort + CoreML) - `src/ml/neural_engine.rs`
- **Format**: ONNX (opset 13) for CoreML compatibility

## Use Case: Read Quality Prediction

**Problem**: FASTQ quality filtering is compute-intensive for large datasets

**Traditional Approach**: Iterate over quality scores, compute average, compare threshold
- Latency: ~10-50 µs per read (CPU scalar)
- Throughput: ~20K-100K reads/sec

**ML Approach**: Encode sequence + quality → Neural Engine inference → PASS/FAIL
- Latency: <1 ms per read (includes encoding)
- Throughput: >10K reads/sec (expected)
- Advantage: Can learn complex quality patterns beyond simple averages

### Model Architecture

```
Input Layer:     300 features (150bp read × 2 for sequence + quality)
Hidden Layer 1:  128 neurons + ReLU + Dropout(0.2)
Hidden Layer 2:  64 neurons + ReLU + Dropout(0.2)
Output Layer:    1 neuron + Sigmoid → probability (PASS/FAIL)
```

**Total Parameters**: ~51,000 (lightweight for Neural Engine)

### Feature Encoding

**Sequence Encoding**:
- A → 0.0, C → 0.25, G → 0.5, T → 0.75, N → 1.0
- Normalized to [0, 1] for neural network input

**Quality Encoding**:
- Phred scores (ASCII - 33) → [0, 93]
- Normalized to [0, 1]: `phred / 93.0`

**Final Features**:
- First 150 values: Encoded sequence
- Next 150 values: Encoded quality scores
- Total: 300 features

## Quick Start

### 1. Install Dependencies

```bash
# Python dependencies (training)
pip install torch numpy scikit-learn onnx onnxruntime-coreml

# Rust dependencies (inference)
cargo build --features neural-engine
```

### 2. Train Model

```bash
# Train on existing FASTQ data
python research/neural-engine/train_quality_model.py \
    --fastq tests/data/sample.fastq.gz \
    --max-records 50000 \
    --epochs 20 \
    --output quality_model.onnx
```

**Expected Output**:
```
Epoch 20/20 | Train Loss: 0.1234 | Val Loss: 0.1456
Test Set Evaluation
Accuracy:  0.9500
Precision: 0.9400
Recall:    0.9600
F1 Score:  0.9500
Model exported to ONNX: quality_model.onnx
```

### 3. Test ONNX Model

```bash
# Verify ONNX model works correctly
python research/neural-engine/test_onnx_inference.py \
    --model quality_model.onnx
```

**Expected Output**:
```
Execution provider: CoreMLExecutionProvider
✅ Neural Engine available!

Test 1: High quality read (should PASS)
Probability: 0.9842
Prediction:  ✅ PASS
Result: ✅ CORRECT
```

### 4. Test Rust Integration

```bash
# Test Neural Engine integration in Rust
cargo test --features neural-engine test_quality_predictor

# Run example
cargo run --example neural_quality --features neural-engine
```

### 5. Benchmark Performance

```bash
# Benchmark Neural Engine vs CPU
cargo bench --features neural-engine neural_engine
```

## Files

### Training Scripts

- **`train_quality_model.py`**: Train PyTorch model and export to ONNX
  - Parses FASTQ files
  - Generates training data with quality_filter labels
  - Trains simple MLP classifier
  - Exports to ONNX format

- **`test_onnx_inference.py`**: Test ONNX model inference
  - Verifies CoreML execution provider works
  - Tests sample inputs
  - Validates predictions

### Research Documentation

- **`RESEARCH_NOTES.md`**: Comprehensive research on Neural Engine capabilities
  - Neural Engine architecture (M4: 38 TOPS, 16 cores)
  - Core ML framework overview
  - DNA sequence ML state-of-the-art
  - Task prioritization and feasibility analysis

### Rust Implementation

- **`src/ml/neural_engine.rs`**: ONNX Runtime wrapper with CoreML backend
  - Session creation with Neural Engine dispatch
  - Tensor input/output handling
  - Error handling and metadata extraction

- **`src/ml/quality.rs`**: Read quality predictor using Neural Engine
  - Feature encoding (sequence + quality)
  - Inference with <1ms latency
  - Probability and binary prediction outputs

## Training Data

The training script uses biometal's existing `quality_filter()` function as ground truth:

```python
def quality_filter_python(sequence: str, quality: str, min_quality: int = 20) -> bool:
    """Returns True if average quality >= min_quality"""
    phred_scores = [ord(q) - 33 for q in quality]
    avg_quality = sum(phred_scores) / len(phred_scores)
    return avg_quality >= min_quality
```

**Labeling**:
- PASS (1.0): Average quality >= 20
- FAIL (0.0): Average quality < 20

**Dataset Size**: Recommend 10,000-100,000 reads for training
- Train: 70% (stratified split)
- Val: 15% (stratified split)
- Test: 15% (stratified split)

## Performance Expectations

### Neural Engine Specifications (M4)

- **Compute**: 38 TOPS (Tera Operations Per Second)
- **Cores**: 16 Neural Engine cores
- **Memory**: Unified memory (zero-copy with CPU)
- **Power**: 10-100× more efficient than GPU

### Expected Performance

| Metric | Traditional (CPU) | Neural Engine | Improvement |
|--------|------------------|---------------|-------------|
| Latency (per read) | ~10-50 µs | <1 ms | ~1-5× slower |
| Throughput (batch) | ~20K-100K/s | >10K/s | Competitive |
| Power efficiency | Baseline | 10-100× better | Significant |
| Complex patterns | Limited | Excellent | New capability |

**Key Insight**: Neural Engine is best for:
- Batch inference (amortize overhead)
- Complex pattern recognition (beyond simple thresholds)
- Power-constrained environments
- Learning from labeled data

## ONNX Runtime v2.0 API

### Session Creation

```rust
use ort::{
    execution_providers::CoreMLExecutionProvider,
    session::{Session, builder::GraphOptimizationLevel},
    value::Tensor,
};

let session = Session::builder()?
    .with_execution_providers([
        CoreMLExecutionProvider::default()
            .build()
            .error_on_failure()
    ])?
    .with_optimization_level(GraphOptimizationLevel::Level3)?
    .commit_from_file(model_path)?;
```

### Inference

```rust
// Create input tensor
let input_tensor = Tensor::from_array((
    vec![1, 300],  // shape: [batch_size, features]
    features.to_vec()  // data: Vec<f32>
))?;

// Run inference
let outputs = session.run(ort::inputs![input_tensor])?;

// Extract output
let output_array: ndarray::ArrayViewD<f32> = outputs[0].try_extract_array()?;
let probability = output_array[[0, 0]];
```

## Next Steps

### Week 2 Completion (Days 3-5)

- [x] Create training script
- [x] Create testing script
- [ ] Train model on real FASTQ data
- [ ] Test ONNX model inference
- [ ] Integrate with Rust `QualityPredictor`
- [ ] Create Rust example
- [ ] Benchmark Neural Engine vs CPU
- [ ] Document findings

### Future Work (Week 3+)

- **Read embedding**: Encode sequences for downstream ML tasks
- **Variant calling**: ML-based variant detection
- **Quality score prediction**: Predict quality from sequence alone
- **Adapter detection**: ML-based adapter identification
- **K-mer classification**: Neural embeddings for k-mer analysis

## References

- **Apple Neural Engine**: https://github.com/hollance/neural-engine
- **Core ML**: https://developer.apple.com/documentation/coreml
- **ONNX Runtime**: https://onnxruntime.ai/
- **ort crate**: https://ort.pyke.io/
- **PyTorch**: https://pytorch.org/

## Platform Support

- **macOS ARM64 (M1/M2/M3/M4)**: Full Neural Engine acceleration
- **macOS x86_64**: CPU fallback via ONNX Runtime
- **Other platforms**: Feature disabled (compile-time)

## License

MIT OR Apache-2.0 (same as biometal)
