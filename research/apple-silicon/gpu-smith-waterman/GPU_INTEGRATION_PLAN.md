# GPU Smith-Waterman Integration Plan

**Target**: biometal v1.7.0
**Timeline**: 1-2 weeks
**Status**: Ready to begin
**Rationale**: 771× speedup proven, production-ready code from Week 1 research

---

## Overview

Integrate GPU-accelerated Smith-Waterman alignment from strategic pivot Week 1 research into biometal core library.

**Performance**: 771× speedup for batch alignment (1,024 sequences)
**Technology**: Metal compute shaders (Apple Silicon GPU)
**Source**: `planning_archive/strategic_pivot/smith-waterman-gpu/`

---

## Integration Steps

### 1. Code Organization

**Move GPU Code**:
```bash
# Create GPU module structure
mkdir -p src/alignment/gpu/shaders

# Move Metal shader
cp planning_archive/strategic_pivot/smith-waterman-gpu/shaders/smith_waterman.metal \
   src/alignment/gpu/shaders/

# Create Rust module files (new)
touch src/alignment/gpu/mod.rs
touch src/alignment/gpu/batch.rs
```

**Directory Structure**:
```
src/alignment/
├── mod.rs                          # Update: Export GPU module
├── smith_waterman.rs               # Existing CPU implementation
└── gpu/
    ├── mod.rs                      # GPU module organization
    ├── batch.rs                    # Batch processing API
    └── shaders/
        └── smith_waterman.metal    # Metal compute shader (320 lines)
```

### 2. Cargo.toml Updates

**Add GPU Feature Flag**:
```toml
[features]
default = ["network"]
network = ["dep:reqwest", "dep:lru", "dep:tokio", "dep:bytes"]
python = ["dep:pyo3"]
simd = ["dep:simd-minimizers", "dep:packed-seq"]
gpu = ["dep:metal"]  # NEW: GPU acceleration (macOS only)
neural-engine = ["dep:ort", "dep:ndarray"]

[target.'cfg(target_os = "macos")'.dependencies]
metal = { version = "0.29", optional = true }  # Already present
```

**Add GPU Benchmark**:
```toml
[[bench]]
name = "smith_waterman_gpu"
harness = false
required-features = ["gpu"]
```

### 3. API Design

**Public API** (`src/alignment/gpu/mod.rs`):
```rust
//! GPU-accelerated sequence alignment using Metal
//!
//! Provides batch Smith-Waterman alignment with 771× speedup
//! on Apple Silicon GPUs.

#[cfg(feature = "gpu")]
pub mod batch;

#[cfg(feature = "gpu")]
pub use batch::{smith_waterman_batch_gpu, GpuAlignmentBatch};
```

**Batch API** (`src/alignment/gpu/batch.rs`):
```rust
use crate::alignment::{Alignment, ScoringMatrix};
use metal::*;

/// GPU-accelerated batch Smith-Waterman alignment
///
/// # Performance
///
/// - Small batches (<16): ~1× (GPU overhead dominates)
/// - Medium batches (256): ~400× speedup
/// - Large batches (1024): ~771× speedup
///
/// # Optimal Use Cases
///
/// - Multiple sequence alignment (MSA)
/// - Database searches (BLAST-like)
/// - Variant calling (many reads vs reference)
///
/// # Example
///
/// ```no_run
/// use biometal::alignment::gpu::smith_waterman_batch_gpu;
/// use biometal::alignment::ScoringMatrix;
///
/// let queries = vec![b"ACGT", b"TGCA", b"GGCC"];
/// let targets = vec![b"ACGG", b"TGCA", b"GGCT"];
/// let scoring = ScoringMatrix::default();
///
/// let alignments = smith_waterman_batch_gpu(&queries, &targets, &scoring)?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn smith_waterman_batch_gpu(
    queries: &[&[u8]],
    targets: &[&[u8]],
    scoring: &ScoringMatrix,
) -> Result<Vec<Alignment>, String> {
    // Implementation from Week 1 research
    // ...
}

/// Reusable GPU batch processor
///
/// Avoids recreating Metal resources for multiple batches.
pub struct GpuAlignmentBatch {
    device: Device,
    pipeline: ComputePipelineState,
    command_queue: CommandQueue,
}

impl GpuAlignmentBatch {
    /// Create a new GPU batch processor
    pub fn new() -> Result<Self, String> {
        // Initialize Metal resources
        // ...
    }

    /// Process a batch of alignments
    pub fn align_batch(
        &self,
        queries: &[&[u8]],
        targets: &[&[u8]],
        scoring: &ScoringMatrix,
    ) -> Result<Vec<Alignment>, String> {
        // Batch processing implementation
        // ...
    }
}
```

### 4. Module Integration

**Update `src/alignment/mod.rs`**:
```rust
//! Sequence alignment algorithms
//!
//! Provides CPU and GPU-accelerated Smith-Waterman alignment.

pub mod smith_waterman;

#[cfg(feature = "gpu")]
pub mod gpu;

pub use smith_waterman::{smith_waterman, Alignment, CigarOp, ScoringMatrix};

#[cfg(feature = "gpu")]
pub use gpu::{smith_waterman_batch_gpu, GpuAlignmentBatch};
```

**Update `src/lib.rs`**:
```rust
pub mod alignment;

// Re-export commonly used types
pub use alignment::{smith_waterman, Alignment, CigarOp, ScoringMatrix};

#[cfg(feature = "gpu")]
pub use alignment::{smith_waterman_batch_gpu, GpuAlignmentBatch};
```

### 5. Testing

**Property Tests** (`src/alignment/gpu/batch.rs`):
```rust
#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn test_gpu_availability() {
        // Check if GPU is available on macOS
        #[cfg(target_os = "macos")]
        assert!(Device::system_default().is_some());
    }

    proptest! {
        #[test]
        fn test_gpu_matches_cpu(
            query in "[ACGT]{10,50}",
            target in "[ACGT]{10,50}"
        ) {
            let scoring = ScoringMatrix::default();

            // CPU result
            let cpu_result = smith_waterman(
                query.as_bytes(),
                target.as_bytes(),
                &scoring
            );

            // GPU result (batch of 1)
            let gpu_result = smith_waterman_batch_gpu(
                &[query.as_bytes()],
                &[target.as_bytes()],
                &scoring
            ).unwrap();

            prop_assert_eq!(cpu_result.score, gpu_result[0].score);
        }
    }
}
```

**Integration Tests** (`tests/gpu_alignment.rs`):
```rust
#![cfg(feature = "gpu")]

use biometal::alignment::gpu::smith_waterman_batch_gpu;
use biometal::alignment::ScoringMatrix;

#[test]
fn test_batch_alignment() {
    let queries = vec![
        b"ACGTACGT".as_slice(),
        b"TGCATGCA".as_slice(),
        b"GGCCGGCC".as_slice(),
    ];
    let targets = vec![
        b"ACGGACGG".as_slice(),
        b"TGCATGCA".as_slice(),
        b"GGCTGGCT".as_slice(),
    ];
    let scoring = ScoringMatrix::default();

    let results = smith_waterman_batch_gpu(&queries, &targets, &scoring).unwrap();

    assert_eq!(results.len(), 3);
    assert!(results[1].score > results[0].score); // Perfect match
}

#[test]
fn test_large_batch() {
    // Test 1,024 alignments (optimal batch size)
    let queries: Vec<&[u8]> = (0..1024)
        .map(|_| b"ACGTACGTACGTACGT".as_slice())
        .collect();
    let targets: Vec<&[u8]> = (0..1024)
        .map(|_| b"ACGTACGTACGTACGT".as_slice())
        .collect();
    let scoring = ScoringMatrix::default();

    let results = smith_waterman_batch_gpu(&queries, &targets, &scoring).unwrap();

    assert_eq!(results.len(), 1024);
}
```

### 6. Benchmarking

**Move Benchmark** (`benches/smith_waterman_gpu.rs`):
```rust
use biometal::alignment::{smith_waterman, smith_waterman_batch_gpu, ScoringMatrix};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

fn generate_sequences(count: usize, length: usize) -> Vec<Vec<u8>> {
    // Generate random DNA sequences
    // ...
}

fn bench_cpu_vs_gpu(c: &mut Criterion) {
    let mut group = c.benchmark_group("smith_waterman_cpu_vs_gpu");

    for &batch_size in &[1, 16, 256, 1024] {
        group.throughput(Throughput::Elements(batch_size));

        let queries = generate_sequences(batch_size as usize, 150);
        let targets = generate_sequences(batch_size as usize, 150);
        let scoring = ScoringMatrix::default();

        // CPU (sequential)
        group.bench_with_input(
            BenchmarkId::new("cpu_sequential", batch_size),
            &batch_size,
            |b, _| {
                b.iter(|| {
                    for (query, target) in queries.iter().zip(targets.iter()) {
                        black_box(smith_waterman(query, target, &scoring));
                    }
                });
            },
        );

        // GPU (batch)
        group.bench_with_input(
            BenchmarkId::new("gpu_batch", batch_size),
            &batch_size,
            |b, _| {
                let query_refs: Vec<&[u8]> = queries.iter().map(|q| q.as_slice()).collect();
                let target_refs: Vec<&[u8]> = targets.iter().map(|t| t.as_slice()).collect();

                b.iter(|| {
                    black_box(smith_waterman_batch_gpu(&query_refs, &target_refs, &scoring).unwrap());
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_cpu_vs_gpu);
criterion_main!(benches);
```

### 7. Documentation

**User Guide Section** (`docs/USER_GUIDE.md`):
```markdown
### GPU-Accelerated Alignment

biometal provides GPU-accelerated Smith-Waterman alignment for batch processing.

**When to Use GPU**:
- Batch size ≥256 alignments (400-771× speedup)
- Multiple sequence alignment (MSA)
- Database searches
- Variant calling (many reads vs reference)

**When to Use CPU**:
- Single alignments or small batches (<16)
- GPU overhead dominates for small workloads

**Example**:
```python
import biometal

# GPU batch alignment (requires gpu feature)
queries = [b"ACGT" * 50 for _ in range(1024)]
targets = [b"ACGG" * 50 for _ in range(1024)]

alignments = biometal.smith_waterman_batch_gpu(queries, targets)
# 771× faster than sequential CPU!
```

**Performance**:
- 1 alignment: ~1× (GPU overhead)
- 256 alignments: ~400× speedup
- 1,024 alignments: ~771× speedup

**Requirements**:
- macOS with Apple Silicon (M1/M2/M3/M4)
- `gpu` feature enabled: `pip install biometal-rs[gpu]`
```

**API Documentation** (`src/alignment/gpu/batch.rs`):
```rust
//! GPU-accelerated batch Smith-Waterman alignment
//!
//! # Performance Characteristics
//!
//! Based on benchmarks (N=30) on Apple M4 Max:
//!
//! | Batch Size | CPU Time | GPU Time | Speedup |
//! |------------|----------|----------|---------|
//! | 1          | 14.9 µs  | 13.5 µs  | 1.1×    |
//! | 16         | 236 µs   | 164 µs   | 1.4×    |
//! | 256        | 3.81 ms  | 9.45 µs  | 403×    |
//! | 1,024      | 15.3 ms  | 19.8 µs  | 771×    |
//!
//! # Optimal Use Cases
//!
//! - **Multiple Sequence Alignment (MSA)**: Align many sequences to consensus
//! - **Database Search**: BLAST-like searches (query vs many targets)
//! - **Variant Calling**: Align many reads to reference genome
//!
//! # Example
//!
//! [example code]
```

### 8. Examples

**Create GPU Example** (`examples/gpu_alignment.rs`):
```rust
//! GPU-accelerated Smith-Waterman alignment example
//!
//! Demonstrates batch alignment with 771× speedup.

#[cfg(not(feature = "gpu"))]
fn main() {
    eprintln!("Error: This example requires the 'gpu' feature");
    eprintln!("Run with: cargo run --example gpu_alignment --features gpu");
    std::process::exit(1);
}

#[cfg(feature = "gpu")]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    use biometal::alignment::{smith_waterman_batch_gpu, ScoringMatrix};
    use std::time::Instant;

    println!("GPU Smith-Waterman Batch Alignment Example");
    println!("===========================================\n");

    // Generate test data
    let batch_size = 1024;
    let queries: Vec<&[u8]> = (0..batch_size)
        .map(|_| b"ACGTACGTACGTACGT".as_slice())
        .collect();
    let targets: Vec<&[u8]> = (0..batch_size)
        .map(|_| b"ACGTACGTACGTACGT".as_slice())
        .collect();

    println!("Batch size: {} alignments", batch_size);
    println!("Sequence length: 16 bp\n");

    // GPU batch alignment
    let scoring = ScoringMatrix::default();
    let start = Instant::now();
    let results = smith_waterman_batch_gpu(&queries, &targets, &scoring)?;
    let gpu_time = start.elapsed();

    println!("GPU Results:");
    println!("  Time: {:?}", gpu_time);
    println!("  Throughput: {:.0} alignments/sec", batch_size as f64 / gpu_time.as_secs_f64());
    println!("  Per-alignment: {:?}", gpu_time / batch_size);
    println!("  First alignment score: {}", results[0].score);

    Ok(())
}
```

### 9. Python Bindings

**Update `src/python.rs`**:
```rust
#[cfg(feature = "gpu")]
#[pyfunction]
fn smith_waterman_batch_gpu(
    queries: Vec<Vec<u8>>,
    targets: Vec<Vec<u8>>,
) -> PyResult<Vec<PyAlignment>> {
    use crate::alignment::gpu::smith_waterman_batch_gpu as gpu_batch;

    let query_refs: Vec<&[u8]> = queries.iter().map(|q| q.as_slice()).collect();
    let target_refs: Vec<&[u8]> = targets.iter().map(|t| t.as_slice()).collect();

    let scoring = ScoringMatrix::default();
    let results = gpu_batch(&query_refs, &target_refs, &scoring)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e))?;

    Ok(results.into_iter().map(PyAlignment::from).collect())
}

#[pymodule]
fn biometal(_py: Python, m: &PyModule) -> PyResult<()> {
    // ... existing functions ...

    #[cfg(feature = "gpu")]
    m.add_function(wrap_pyfunction!(smith_waterman_batch_gpu, m)?)?;

    Ok(())
}
```

### 10. CHANGELOG Update

**Add to CHANGELOG.md** (v1.7.0):
```markdown
## [1.7.0] - 2025-11-XX

### Added

- **GPU-Accelerated Smith-Waterman Alignment** (`gpu` feature):
  - 771× speedup for batch alignment (1,024 sequences)
  - Metal compute shaders for Apple Silicon GPU
  - Batch API: `smith_waterman_batch_gpu(queries, targets, scoring)`
  - Optimal for MSA, database search, variant calling
  - Python bindings: `biometal.smith_waterman_batch_gpu()`
  - Comprehensive benchmarks (N=30)
  - Documentation and examples

- **Neural Engine Infrastructure** (`neural-engine` feature, experimental):
  - ONNX Runtime integration with CoreML backend
  - Custom ML model deployment capability
  - Example: Read quality prediction
  - Preserved for future ML use cases

### Performance

- GPU Smith-Waterman: 771× speedup (1,024-alignment batches)
- Compression (cloudflare_zlib): 1.67× decompression, 2.29× compression (from v1.6.0)
- BAM parsing: 92 MiB/s (cloudflare_zlib backend)

### Documentation

- GPU alignment user guide and API docs
- Strategic pivot research archived
- Integration examples (Rust + Python)
```

---

## Implementation Checklist

### Phase 1: Code Migration (Day 1)

- [ ] Create `src/alignment/gpu/` directory structure
- [ ] Copy Metal shader to `src/alignment/gpu/shaders/smith_waterman.metal`
- [ ] Create `src/alignment/gpu/mod.rs` (module organization)
- [ ] Create `src/alignment/gpu/batch.rs` (batch API implementation)
- [ ] Update `src/alignment/mod.rs` (export GPU module)
- [ ] Update `src/lib.rs` (re-export GPU API)

### Phase 2: Build System (Day 1-2)

- [ ] Add `gpu` feature flag to `Cargo.toml`
- [ ] Verify `metal` dependency (already present for alignment)
- [ ] Add GPU benchmark entry to `Cargo.toml`
- [ ] Test compilation: `cargo build --features gpu`
- [ ] Test without feature: `cargo build` (should still work)

### Phase 3: Testing (Day 2-3)

- [ ] Add property tests (GPU matches CPU)
- [ ] Add integration tests (`tests/gpu_alignment.rs`)
- [ ] Test batch sizes: 1, 16, 256, 1024
- [ ] Test edge cases: empty sequences, long sequences
- [ ] Run benchmarks: `cargo bench --features gpu smith_waterman`
- [ ] Verify 771× speedup maintained

### Phase 4: Documentation (Day 3-4)

- [ ] Update user guide (`docs/USER_GUIDE.md`)
- [ ] Add API documentation to `src/alignment/gpu/batch.rs`
- [ ] Create GPU example (`examples/gpu_alignment.rs`)
- [ ] Update CHANGELOG.md for v1.7.0
- [ ] Document when to use GPU vs CPU

### Phase 5: Python Bindings (Day 4-5)

- [ ] Add `smith_waterman_batch_gpu()` to `src/python.rs`
- [ ] Test Python bindings: `import biometal; biometal.smith_waterman_batch_gpu(...)`
- [ ] Update Python documentation
- [ ] Create Python example notebook

### Phase 6: Quality Assurance (Day 5-6)

- [ ] Run full test suite: `cargo test --all-features`
- [ ] Run benchmarks: `cargo bench --features gpu`
- [ ] Property testing: `cargo test --features gpu proptest`
- [ ] Cross-platform test (x86_64 fallback)
- [ ] Memory leak check (instruments on macOS)

### Phase 7: Release Prep (Day 7)

- [ ] Final CHANGELOG.md review
- [ ] Version bump: v1.6.0 → v1.7.0
- [ ] Create release notes
- [ ] Tag release: `git tag v1.7.0`
- [ ] Publish to crates.io: `cargo publish`
- [ ] Publish Python bindings: `maturin publish`

---

## Success Criteria

- ✅ GPU code integrated into `src/alignment/gpu/`
- ✅ `gpu` feature flag working
- ✅ 771× speedup maintained in benchmarks
- ✅ All tests passing (property tests, integration tests)
- ✅ Documentation complete (user guide, API docs, examples)
- ✅ Python bindings functional
- ✅ CHANGELOG.md updated for v1.7.0
- ✅ Ready for crates.io + PyPI release

---

## Timeline

**Week 1** (Days 1-5):
- Days 1-2: Code migration + build system
- Days 2-3: Testing
- Days 3-4: Documentation
- Days 4-5: Python bindings
- Day 5: QA

**Week 2** (Days 6-7):
- Day 6: Final testing and polish
- Day 7: Release prep and publish

**Total**: 7-10 days (1-2 weeks)

---

## Next Steps

1. Begin Phase 1: Code migration
2. Set up GPU feature flag
3. Implement batch API
4. Run property tests (GPU matches CPU)
5. Benchmark to verify 771× speedup

**Status**: Ready to begin integration
**Target**: v1.7.0 release in 1-2 weeks
