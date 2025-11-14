//! GPU-accelerated sequence alignment using Metal
//!
//! This module provides GPU-accelerated Smith-Waterman alignment for Apple Silicon
//! using the Metal framework.
//!
//! # Performance
//!
//! GPU acceleration provides modest but consistent speedups for batch alignment
//! (Measured on Apple M-series, 500bp sequences, N=30):
//!
//! | Batch Size | CPU Time | GPU Time | Speedup |
//! |------------|----------|----------|---------|
//! | 1 (100bp)  | 32 µs    | 2.5 ms   | 0.01× (overhead) |
//! | 1 (500bp)  | 430 µs   | 964 µs   | 0.45× (overhead) |
//! | 1 (1000bp) | 1.72 ms  | 2.14 ms  | 0.80× |
//! | 10         | 4.33 ms  | 3.17 ms  | **1.36×** ✓ |
//! | 50         | 21.9 ms  | 18.4 ms  | **1.19×** ✓ |
//! | 100        | 44.2 ms  | 33.7 ms  | **1.31×** ✓ |
//! | 500        | 218 ms   | 174 ms   | **1.26×** ✓ |
//! | 1,000      | 435 ms   | 345 ms   | **1.26×** ✓ |
//!
//! **Dispatch overhead**: ~598 µs (batch creation, one-time cost)
//!
//! **Note**: Current implementation uses one GPU thread per alignment (serial DP within
//! each alignment). Speedups are modest (1.2-1.4×) but consistent. Future optimization
//! with anti-diagonal parallelization could achieve 10-50× from literature.
//!
//! # When to Use GPU
//!
//! GPU acceleration is most effective for:
//! - **Batch processing**: ≥10 alignments (amortizes dispatch overhead)
//! - **Medium-long sequences**: ≥500bp (more computation per alignment)
//! - **Repeated batches**: Create `GpuAlignmentBatch` once, reuse many times
//!
//! Use CPU for:
//! - Single alignments: GPU overhead dominates (0.4-0.8× slower)
//! - Very short sequences: <100bp (not enough work to justify GPU)
//! - Small batches: <10 alignments (marginal benefit)
//!
//! # Platform Support
//!
//! - **macOS with Apple Silicon** (M1/M2/M3/M4): Full GPU acceleration
//! - **macOS with Intel**: Falls back to CPU (Metal available but slower)
//! - **Other platforms**: Feature disabled at compile time
//!
//! # Examples
//!
//! ## Quick batch alignment
//!
//! ```no_run
//! use biometal::alignment::gpu::smith_waterman_batch_gpu;
//! use biometal::alignment::ScoringMatrix;
//!
//! let queries = vec![b"ACGTACGT", b"TGCATGCA"];
//! let targets = vec![b"ACGGACGG", b"TGCATGCA"];
//! let query_refs: Vec<&[u8]> = queries.iter().map(|q| q.as_slice()).collect();
//! let target_refs: Vec<&[u8]> = targets.iter().map(|t| t.as_slice()).collect();
//!
//! let scoring = ScoringMatrix::default();
//! let alignments = smith_waterman_batch_gpu(&query_refs, &target_refs, &scoring)?;
//!
//! for (i, alignment) in alignments.iter().enumerate() {
//!     println!("Alignment {}: score = {}", i, alignment.score);
//! }
//! # Ok::<(), String>(())
//! ```
//!
//! ## Reusable batch processor (faster for multiple batches)
//!
//! ```no_run
//! use biometal::alignment::gpu::GpuAlignmentBatch;
//! use biometal::alignment::ScoringMatrix;
//!
//! // Create GPU batch processor once
//! let mut gpu_batch = GpuAlignmentBatch::new()?;
//! let scoring = ScoringMatrix::default();
//!
//! // Process multiple batches without recreating Metal resources
//! for batch_id in 0..10 {
//!     let queries = vec![b"ACGTACGT", b"TGCATGCA"];
//!     let targets = vec![b"ACGGACGG", b"TGCATGCA"];
//!     let query_refs: Vec<&[u8]> = queries.iter().map(|q| q.as_slice()).collect();
//!     let target_refs: Vec<&[u8]> = targets.iter().map(|t| t.as_slice()).collect();
//!
//!     let alignments = gpu_batch.align_batch(&query_refs, &target_refs, &scoring)?;
//!     println!("Batch {}: {} alignments", batch_id, alignments.len());
//! }
//! # Ok::<(), String>(())
//! ```
//!
//! # Implementation Details
//!
//! The GPU implementation uses Metal compute shaders to parallelize Smith-Waterman:
//!
//! 1. **Sequence packing**: Concatenate all sequences into contiguous buffers
//! 2. **GPU dispatch**: Launch one thread per alignment
//! 3. **DP computation**: Each thread computes its own DP matrix
//! 4. **Traceback**: Find max score and traceback to start position
//! 5. **CIGAR generation**: Convert traceback path to CIGAR string
//!
//! **Memory layout**: Each alignment gets its own DP matrix to avoid synchronization.
//! This trades memory for parallelism (essential for GPU performance).
//!
//! **Shader**: `shaders/smith_waterman.metal` implements the DP algorithm in Metal Shading Language.

pub mod batch;

pub use batch::{smith_waterman_batch_gpu, GpuAlignmentBatch};
