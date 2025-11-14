# Smith-Waterman GPU Implementation - Research Notes

**Date**: November 13, 2025
**Goal**: Implement Smith-Waterman alignment on Apple Silicon GPU using Metal
**Target**: 10-50× speedup vs CPU (based on CUDA literature)

---

## Research Findings

### 1. CUDA Smith-Waterman Literature

**Key Papers**:
- CUDASW++ (2009): First major GPU implementation
- CUDASW++ 2.0 (2010): SIMT abstraction-based optimization
- CUDASW++ 3.0 (2013): Coupling CPU and GPU SIMD instructions

**Reported Speedups**: 10-50× over CPU implementations (varies by sequence length, GPU hardware)

**Key Parallelization Strategies**:

1. **Divide and Conquer**:
   - Split alignment matrices into sub-matrices
   - Each GPU thread processes one sub-matrix
   - Run sub-matrices in parallel on GPU hardware

2. **Thread Warps**:
   - Group threads into warps for efficient execution
   - Use padding techniques to reduce synchronization costs
   - Minimize conditional branches (loop unrolling)

3. **SIMT (Single Instruction, Multiple Thread)**:
   - All threads execute same instruction on different data
   - Ideal for Smith-Waterman (same algorithm, different sequence pairs)

4. **Memory Optimization**:
   - Coalesced memory access patterns
   - Shared memory for scoring matrices
   - Minimize global memory transfers

**Why GPU Works for Smith-Waterman**:
- Complexity >0.70 (dynamic programming with nested loops)
- Highly parallel (align 1000s of sequence pairs independently)
- Regular memory access (scoring matrix lookups)
- Limited branching (mostly arithmetic operations)

---

### 2. Metal Compute Shader Capabilities

**Metal vs CUDA**:
- Similar capabilities (massively parallel compute)
- Metal: `kernel` functions (analogous to CUDA `__global__`)
- Metal: Threadgroup (analogous to CUDA block)
- Metal: Thread position in grid (analogous to threadIdx, blockIdx)

**Rust + Metal Integration**:
- `metal-rs` crate: Rust wrapper for Metal API
- Example: LambdaClass FFT implementation (2024)
- Pattern:
  1. Write compute shader in Metal Shading Language (MSL)
  2. Compile shader at runtime or build-time
  3. Create Metal buffers for input/output
  4. Dispatch threads to GPU
  5. Read results back from GPU

**Apple Silicon Advantage**:
- **Unified Memory**: Zero-copy CPU↔GPU transfers
- **Low Latency**: ~1-3ms dispatch overhead (vs 3-5ms on discrete GPUs)
- **Integrated**: No PCIe bottleneck
- **Power Efficient**: ~10-15W GPU power vs 150-300W discrete GPUs

---

### 3. Smith-Waterman Algorithm Overview

**Purpose**: Local sequence alignment (find best matching region)

**Algorithm** (Classic Dynamic Programming):
```
H(i,j) = max(
    H(i-1, j-1) + S(a_i, b_j),  // Match/mismatch
    H(i-1, j) + gap_penalty,     // Deletion
    H(i, j-1) + gap_penalty,     // Insertion
    0                             // Start new alignment
)

where:
- H(i,j) = alignment score at position (i,j)
- S(a,b) = scoring matrix (match=2, mismatch=-1)
- gap_penalty = -1 (linear gap penalty)
```

**Traceback**: Follow max scores backwards to reconstruct alignment (generates CIGAR string)

**Complexity**: O(m × n) where m, n are sequence lengths
- For 1000bp sequences: 1,000,000 operations per alignment
- **Key insight**: Each alignment is independent → perfect for GPU parallelization

---

## Implementation Strategy

### Phase 1: CPU Baseline (This Week)

**Goal**: Establish correct implementation and baseline performance

**Implementations**:
1. **Naive Scalar** (Rust):
   - Classic DP algorithm (nested loops)
   - No optimizations
   - Reference for correctness

2. **NEON-Optimized** (Rust + ARM SIMD):
   - Striped algorithm (process multiple cells in parallel)
   - Expected: 2-4× speedup vs naive (irregular memory access limits NEON)
   - Portable to Linux ARM (Graviton)

**Output**: Alignment score + CIGAR string

---

### Phase 2: GPU Implementation (This Week)

**Goal**: Massively parallel alignment using Metal

**Design**:
```metal
kernel void smith_waterman_parallel(
    device const uint8_t* query_seqs [[buffer(0)]],      // N query sequences (packed)
    device const uint32_t* query_lengths [[buffer(1)]],  // N query lengths
    device const uint8_t* ref_seqs [[buffer(2)]],        // N reference sequences
    device const uint32_t* ref_lengths [[buffer(3)]],    // N reference lengths
    device int* scores [[buffer(4)]],                    // N output scores
    device uint8_t* cigars [[buffer(5)]],                // N CIGAR strings
    constant ScoringMatrix& scoring [[buffer(6)]],       // Match/mismatch scores
    uint gid [[thread_position_in_grid]]                 // Thread ID (0 to N-1)
) {
    // Each GPU thread aligns one sequence pair independently
    uint query_offset = ...;  // Calculate from gid
    uint ref_offset = ...;

    // Allocate DP matrix in threadgroup memory (shared)
    threadgroup int H[MAX_LEN * MAX_LEN];

    // Smith-Waterman DP
    for (uint i = 0; i < query_len; i++) {
        for (uint j = 0; j < ref_len; j++) {
            // Compute H(i,j) using scoring matrix
            // ...
        }
    }

    // Traceback to generate CIGAR
    // ...

    // Write results
    scores[gid] = max_score;
    cigars[gid * MAX_CIGAR_LEN ... ] = cigar_string;
}
```

**Parallelization Strategy**:
- **Thread Grid**: N threads (one per sequence pair)
- **Thread Assignment**: Thread `gid` aligns query[gid] to ref[gid]
- **Memory**: Each thread has private DP matrix (threadgroup memory)
- **No Synchronization**: Threads are independent (embarrassingly parallel)

**Expected Performance**:
- **Naive CPU**: ~1,000 alignments/sec (1ms per alignment)
- **NEON CPU**: ~2,000-4,000 alignments/sec (2-4× speedup)
- **GPU**: ~20,000-50,000 alignments/sec (10-50× speedup)

**Why This Should Work**:
1. **High complexity** (>0.70): DP dominates, overhead is minimal
2. **Large batch size**: 1000+ alignments >> 10K GPU threshold (amortizes dispatch)
3. **Independent operations**: No thread communication needed
4. **Regular memory access**: Scoring matrix is small, cacheable
5. **Apple Silicon advantage**: Unified memory eliminates PCIe bottleneck

---

### Phase 3: Benchmarking (This Week)

**Methodology** (N=30 statistical rigor):
```rust
// Benchmark configurations
let configs = [
    (100, 100),   // Short sequences (100bp × 100bp)
    (500, 500),   // Medium (500bp × 500bp)
    (1000, 1000), // Long (1000bp × 1000bp)
    (5000, 1000), // Asymmetric (query >> ref)
];

let batch_sizes = [1, 10, 100, 1000, 10000];

for (query_len, ref_len) in configs {
    for batch_size in batch_sizes {
        // Benchmark naive, NEON, GPU
        // N=30 runs each
        // Calculate: mean, std dev, 95% CI, Cohen's d
    }
}
```

**Metrics**:
- Throughput: Alignments per second
- Latency: Milliseconds per alignment
- Speedup: GPU / CPU ratio
- Scaling: Throughput vs batch size
- Energy: Watts per 1000 alignments (if measurable)

**Success Criteria**:
- ✅ **GPU ≥10× faster than naive CPU**: Major success
- ⚠️ **GPU 5-10× faster**: Moderate success
- ❌ **GPU <5× faster**: Investigate or pivot

---

## Key Differences: CUDA vs Metal

| Feature | CUDA | Metal |
|---------|------|-------|
| **API Language** | C/C++ | Metal Shading Language (MSL) |
| **Kernel Declaration** | `__global__` | `kernel` |
| **Thread Organization** | Block + Grid | Threadgroup + Grid |
| **Thread ID** | `threadIdx`, `blockIdx` | `thread_position_in_threadgroup`, `thread_position_in_grid` |
| **Shared Memory** | `__shared__` | `threadgroup` |
| **Memory Transfer** | `cudaMemcpy` (explicit) | Unified memory (implicit) |
| **Dispatch Overhead** | 3-5ms (PCIe) | 1-3ms (unified memory) |
| **Platform** | NVIDIA GPUs | Apple Silicon |

**Apple Silicon Advantage**:
- No PCIe bottleneck (unified memory)
- Lower dispatch overhead (1-3ms vs 3-5ms)
- Better for smaller batch sizes (break-even at 100-1000 alignments vs 10K)

---

## Implementation Risks

### Risk 1: GPU Dispatch Overhead
**Problem**: Metal dispatch costs 1-3ms per kernel launch
**Mitigation**: Batch 1000+ alignments per dispatch (amortize overhead)
**Fallback**: If batch too small, use CPU

### Risk 2: Memory Limits
**Problem**: DP matrix size = query_len × ref_len
**Constraint**: Threadgroup memory ~32KB on Apple Silicon
**Solution**: Limit max sequence length to ~1000bp (1MB per thread)
**Fallback**: For longer sequences, use CPU or tiled GPU algorithm

### Risk 3: Complexity Too Low
**Problem**: If DP overhead dominates compute, GPU won't help
**Evidence**: ASBB showed complexity >0.55 needed for GPU wins
**Smith-Waterman**: Complexity ~0.70-0.80 (DP dominates) ✅
**Mitigation**: Start with 1000bp sequences (1M operations each)

### Risk 4: Metal Debugging Difficulty
**Problem**: GPU debugging harder than CPU
**Mitigation**:
- Start with simple test cases (short sequences)
- Validate GPU output matches CPU output (property tests)
- Use Metal debugger/profiler if available

---

## Next Steps

### This Week (Week 1):
1. ✅ Research complete (this document)
2. ⏳ Design parallel algorithm (detailed pseudocode)
3. ⏳ Create `src/alignment/` module structure
4. ⏳ Implement naive CPU Smith-Waterman
5. ⏳ Implement NEON-optimized version
6. ⏳ Implement Metal GPU version
7. ⏳ Property tests (GPU == CPU correctness)
8. ⏳ Benchmarks (N=30)

### Next Week (Week 2):
- Optimize Metal shader (memory access, thread configuration)
- CIGAR generation (traceback algorithm)
- Comprehensive benchmarking (all configs)
- Blog post: "GPU-Accelerated Smith-Waterman on Apple Silicon"

### Decision Point (End of Week 4):
- **If GPU ≥10× faster**: Continue GPU work (pileup, variant calling)
- **If GPU <5× faster**: Document negative result, pivot to primitives

---

## References

**CUDA Smith-Waterman**:
- CUDASW++: Optimizing Smith-Waterman sequence database searches for CUDA-enabled GPUs (2009)
- CUDASW++ 2.0: Enhanced Smith-Waterman protein database search on CUDA-enabled GPUs (2010)
- CUDASW++ 3.0: Accelerating Smith-Waterman by coupling CPU and GPU SIMD instructions (2013)

**Metal Compute**:
- Metal Shading Language Specification (Apple Developer)
- Introduction to Compute Programming in Metal (Metal by Example)
- Using Metal and Rust to make FFT even faster (LambdaClass, 2024)

**Apple Silicon**:
- Optimize Metal Performance for Apple Silicon Macs (WWDC 2020)
- Apple GPU microarchitecture (GitHub: philipturner/metal-benchmarks)

---

**Research Status**: ✅ Complete
**Next**: Design detailed parallel algorithm
**Date**: November 13, 2025
