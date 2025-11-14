# Smith-Waterman Algorithm Design - Detailed Pseudocode

**Date**: November 13, 2025
**Goal**: Detailed algorithm design for naive CPU, NEON CPU, and Metal GPU implementations
**Status**: Ready for implementation

---

## Table of Contents

1. [Data Structures](#data-structures)
2. [Algorithm Overview](#algorithm-overview)
3. [CPU Naive Implementation](#cpu-naive-implementation)
4. [NEON Optimized Implementation](#neon-optimized-implementation)
5. [Metal GPU Implementation](#metal-gpu-implementation)
6. [Memory Layout](#memory-layout)
7. [Thread Organization](#thread-organization)
8. [Implementation Notes](#implementation-notes)

---

## Data Structures

### Scoring Configuration
```rust
pub struct ScoringMatrix {
    pub match_score: i32,      // +2 (bases match)
    pub mismatch_score: i32,   // -1 (bases don't match)
    pub gap_open: i32,         // -2 (start a gap)
    pub gap_extend: i32,       // -1 (extend existing gap)
}

impl Default for ScoringMatrix {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_score: -1,
            gap_open: -2,
            gap_extend: -1,
        }
    }
}
```

### Alignment Result
```rust
pub struct Alignment {
    pub score: i32,                    // Maximum alignment score
    pub query_start: usize,            // Start position in query
    pub query_end: usize,              // End position in query
    pub ref_start: usize,              // Start position in reference
    pub ref_end: usize,                // End position in reference
    pub cigar: Vec<CigarOp>,          // CIGAR string operations
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CigarOp {
    Match(usize),      // M: alignment match (could be match or mismatch)
    Insertion(usize),  // I: insertion to reference
    Deletion(usize),   // D: deletion from reference
}
```

### DP Matrix Cell
```rust
struct Cell {
    score: i32,           // H(i,j) - alignment score
    direction: Direction, // Traceback direction
}

enum Direction {
    Diagonal,  // Match/mismatch (from H(i-1,j-1))
    Up,        // Deletion (from H(i-1,j))
    Left,      // Insertion (from H(i,j-1))
    None,      // Start of alignment (score = 0)
}
```

---

## Algorithm Overview

### Smith-Waterman Dynamic Programming

The Smith-Waterman algorithm finds the optimal **local** alignment between two sequences using dynamic programming.

**Recurrence Relation**:
```
H(i,j) = max(
    H(i-1, j-1) + score(query[i], ref[j]),  // Diagonal: match/mismatch
    H(i-1, j) + gap_penalty,                 // Up: deletion
    H(i, j-1) + gap_penalty,                 // Left: insertion
    0                                         // None: start new alignment
)

where:
- H(i,j) = alignment score at position (i,j)
- score(a,b) = match_score if a==b else mismatch_score
- gap_penalty = gap_open (first gap) or gap_extend (continuation)
- 0 = local alignment can start anywhere
```

**Two Phases**:
1. **Forward Pass**: Fill DP matrix, track maximum score
2. **Traceback**: Follow directions from max score to reconstruct alignment

**Complexity**: O(m × n) time, O(m × n) space
- m = query length
- n = reference length

---

## CPU Naive Implementation

### Algorithm: Classical Dynamic Programming

```rust
pub fn smith_waterman_naive(
    query: &[u8],
    reference: &[u8],
    scoring: &ScoringMatrix,
) -> Alignment {
    let m = query.len();
    let n = reference.len();

    // Step 1: Initialize DP matrix (m+1 × n+1)
    let mut H = vec![vec![Cell { score: 0, direction: Direction::None }; n + 1]; m + 1];

    // First row and column are all zeros (local alignment)
    // No initialization needed (already zero)

    // Step 2: Fill DP matrix (forward pass)
    let mut max_score = 0;
    let mut max_i = 0;
    let mut max_j = 0;

    for i in 1..=m {
        for j in 1..=n {
            // Calculate scores from three directions
            let match_score = if query[i-1] == reference[j-1] {
                scoring.match_score
            } else {
                scoring.mismatch_score
            };

            let diagonal = H[i-1][j-1].score + match_score;
            let up = H[i-1][j].score + scoring.gap_open;  // Simplified: always gap_open
            let left = H[i][j-1].score + scoring.gap_open;

            // Take maximum (including 0 for local alignment)
            let (score, direction) = max4(
                (diagonal, Direction::Diagonal),
                (up, Direction::Up),
                (left, Direction::Left),
                (0, Direction::None)
            );

            H[i][j] = Cell { score, direction };

            // Track maximum score for traceback
            if score > max_score {
                max_score = score;
                max_i = i;
                max_j = j;
            }
        }
    }

    // Step 3: Traceback to reconstruct alignment
    let cigar = traceback(&H, max_i, max_j);

    Alignment {
        score: max_score,
        query_start: max_i - cigar_length(&cigar),
        query_end: max_i,
        ref_start: max_j - cigar_length(&cigar),
        ref_end: max_j,
        cigar,
    }
}

fn max4(a: (i32, Direction), b: (i32, Direction),
        c: (i32, Direction), d: (i32, Direction)) -> (i32, Direction) {
    let max_ab = if a.0 >= b.0 { a } else { b };
    let max_cd = if c.0 >= d.0 { c } else { d };
    if max_ab.0 >= max_cd.0 { max_ab } else { max_cd }
}

fn traceback(H: &[Vec<Cell>], start_i: usize, start_j: usize) -> Vec<CigarOp> {
    let mut cigar = Vec::new();
    let mut i = start_i;
    let mut j = start_j;

    while i > 0 && j > 0 && H[i][j].direction != Direction::None {
        match H[i][j].direction {
            Direction::Diagonal => {
                cigar.push(CigarOp::Match(1));
                i -= 1;
                j -= 1;
            }
            Direction::Up => {
                cigar.push(CigarOp::Deletion(1));
                i -= 1;
            }
            Direction::Left => {
                cigar.push(CigarOp::Insertion(1));
                j -= 1;
            }
            Direction::None => break,
        }
    }

    // Reverse (traceback is backwards)
    cigar.reverse();

    // Compress consecutive operations (e.g., M M M → 3M)
    compress_cigar(cigar)
}

fn compress_cigar(cigar: Vec<CigarOp>) -> Vec<CigarOp> {
    if cigar.is_empty() {
        return cigar;
    }

    let mut compressed = Vec::new();
    let mut current = cigar[0];
    let mut count = 1;

    for &op in &cigar[1..] {
        if std::mem::discriminant(&op) == std::mem::discriminant(&current) {
            count += 1;
        } else {
            compressed.push(match current {
                CigarOp::Match(_) => CigarOp::Match(count),
                CigarOp::Insertion(_) => CigarOp::Insertion(count),
                CigarOp::Deletion(_) => CigarOp::Deletion(count),
            });
            current = op;
            count = 1;
        }
    }

    // Push final operation
    compressed.push(match current {
        CigarOp::Match(_) => CigarOp::Match(count),
        CigarOp::Insertion(_) => CigarOp::Insertion(count),
        CigarOp::Deletion(_) => CigarOp::Deletion(count),
    });

    compressed
}
```

**Performance**: ~1,000 alignments/sec (1ms per 1000bp × 1000bp alignment)

**Memory**: O(m × n) - full DP matrix

---

## NEON Optimized Implementation

### Algorithm: Striped SIMD

The NEON implementation uses a "striped" approach where multiple DP cells are computed in parallel using SIMD.

**Key Idea**: Process diagonal stripes of the DP matrix in parallel
- Each SIMD lane processes one cell in the diagonal
- 16 cells processed simultaneously (ARM NEON processes 16 bytes at once)

**Challenges**:
- **Data dependencies**: H(i,j) depends on H(i-1,j-1), H(i-1,j), H(i,j-1)
- **Irregular memory access**: Diagonal access pattern is not cache-friendly
- **Expected speedup**: 2-4× (less than 16× due to dependencies)

```rust
#[cfg(target_arch = "aarch64")]
pub unsafe fn smith_waterman_neon(
    query: &[u8],
    reference: &[u8],
    scoring: &ScoringMatrix,
) -> Alignment {
    use std::arch::aarch64::*;

    let m = query.len();
    let n = reference.len();

    // Allocate DP matrix (only store scores, not directions for now)
    let mut H = vec![vec![0i32; n + 1]; m + 1];

    // SIMD registers for scoring
    let match_vec = vdupq_n_s32(scoring.match_score);
    let mismatch_vec = vdupq_n_s32(scoring.mismatch_score);
    let gap_vec = vdupq_n_s32(scoring.gap_open);
    let zero_vec = vdupq_n_s32(0);

    let mut max_score = 0;
    let mut max_i = 0;
    let mut max_j = 0;

    // Process in 4-element chunks (NEON processes 4 x i32 per operation)
    let chunks = n / 4;

    for i in 1..=m {
        for chunk in 0..chunks {
            let j_start = chunk * 4 + 1;

            // Load query base (broadcast to all lanes)
            let query_base = query[i - 1];

            // Load 4 reference bases
            let ref_bases = [
                reference[j_start - 1],
                reference[j_start],
                reference[j_start + 1],
                reference[j_start + 2],
            ];

            // Calculate match/mismatch scores
            let mut scores = [0i32; 4];
            for k in 0..4 {
                scores[k] = if query_base == ref_bases[k] {
                    scoring.match_score
                } else {
                    scoring.mismatch_score
                };
            }
            let score_vec = vld1q_s32(scores.as_ptr());

            // Load previous values
            let diag = vld1q_s32(&H[i-1][j_start-1] as *const i32);
            let up = vld1q_s32(&H[i-1][j_start] as *const i32);
            let left = vld1q_s32(&H[i][j_start-1] as *const i32);

            // Calculate three options
            let diagonal = vaddq_s32(diag, score_vec);
            let deletion = vaddq_s32(up, gap_vec);
            let insertion = vaddq_s32(left, gap_vec);

            // Take maximum of all four options (diagonal, up, left, 0)
            let max_du = vmaxq_s32(diagonal, deletion);
            let max_il = vmaxq_s32(insertion, zero_vec);
            let result = vmaxq_s32(max_du, max_il);

            // Store results
            vst1q_s32(&mut H[i][j_start] as *mut i32, result);

            // Track maximum (scalar fallback for simplicity)
            for k in 0..4 {
                let score = H[i][j_start + k];
                if score > max_score {
                    max_score = score;
                    max_i = i;
                    max_j = j_start + k;
                }
            }
        }

        // Handle remaining elements (scalar)
        for j in (chunks * 4 + 1)..=n {
            let match_score = if query[i-1] == reference[j-1] {
                scoring.match_score
            } else {
                scoring.mismatch_score
            };

            let diagonal = H[i-1][j-1] + match_score;
            let up = H[i-1][j] + scoring.gap_open;
            let left = H[i][j-1] + scoring.gap_open;

            let score = diagonal.max(up).max(left).max(0);
            H[i][j] = score;

            if score > max_score {
                max_score = score;
                max_i = i;
                max_j = j;
            }
        }
    }

    // Traceback (scalar - same as naive)
    // Note: Would need to store directions for proper traceback
    // For now, return score only

    Alignment {
        score: max_score,
        query_start: 0,  // Would calculate from traceback
        query_end: max_i,
        ref_start: 0,
        ref_end: max_j,
        cigar: vec![],  // Would generate from traceback
    }
}
```

**Performance**: ~2,000-4,000 alignments/sec (2-4× speedup)

**Limitations**:
- Data dependencies limit parallelism
- Irregular memory access reduces cache efficiency
- Traceback still needs scalar implementation

---

## Metal GPU Implementation

### Algorithm: Massively Parallel (One Thread Per Alignment)

Each GPU thread independently aligns one sequence pair. No synchronization needed.

### Metal Shading Language (MSL) Kernel

```metal
// smith_waterman.metal

#include <metal_stdlib>
using namespace metal;

struct ScoringMatrix {
    int match_score;
    int mismatch_score;
    int gap_open;
    int gap_extend;
};

struct AlignmentResult {
    int score;
    uint query_start;
    uint query_end;
    uint ref_start;
    uint ref_end;
};

// Main kernel: Each thread aligns one sequence pair
kernel void smith_waterman_kernel(
    // Input sequences (packed: all queries, then all references)
    device const uint8_t* query_seqs [[buffer(0)]],
    device const uint32_t* query_offsets [[buffer(1)]],
    device const uint32_t* query_lengths [[buffer(2)]],

    device const uint8_t* ref_seqs [[buffer(3)]],
    device const uint32_t* ref_offsets [[buffer(4)]],
    device const uint32_t* ref_lengths [[buffer(5)]],

    // Scoring configuration
    constant ScoringMatrix& scoring [[buffer(6)]],

    // Output results
    device AlignmentResult* results [[buffer(7)]],

    // Thread position
    uint gid [[thread_position_in_grid]]
) {
    // Get this thread's sequence pair
    uint query_offset = query_offsets[gid];
    uint query_len = query_lengths[gid];

    uint ref_offset = ref_offsets[gid];
    uint ref_len = ref_lengths[gid];

    // Pointers to this thread's sequences
    device const uint8_t* query = query_seqs + query_offset;
    device const uint8_t* ref = ref_seqs + ref_offset;

    // Allocate DP matrix in thread-private memory
    // Max supported: 512×512 = 256KB (adjust based on GPU limits)
    const uint MAX_LEN = 512;

    if (query_len > MAX_LEN || ref_len > MAX_LEN) {
        // Sequences too long for GPU, return error score
        results[gid].score = -1;
        return;
    }

    // DP matrix (threadgroup shared memory for potential optimization)
    // For now, use thread-private arrays
    int H[(MAX_LEN + 1) * (MAX_LEN + 1)];

    // Initialize first row and column to 0
    for (uint i = 0; i <= query_len; i++) {
        H[i * (ref_len + 1) + 0] = 0;
    }
    for (uint j = 0; j <= ref_len; j++) {
        H[0 * (ref_len + 1) + j] = 0;
    }

    // Fill DP matrix
    int max_score = 0;
    uint max_i = 0;
    uint max_j = 0;

    for (uint i = 1; i <= query_len; i++) {
        for (uint j = 1; j <= ref_len; j++) {
            // Calculate match/mismatch score
            int match_score = (query[i-1] == ref[j-1]) ?
                scoring.match_score : scoring.mismatch_score;

            // Three options
            int diagonal = H[(i-1) * (ref_len + 1) + (j-1)] + match_score;
            int up = H[(i-1) * (ref_len + 1) + j] + scoring.gap_open;
            int left = H[i * (ref_len + 1) + (j-1)] + scoring.gap_open;

            // Take maximum (including 0 for local alignment)
            int score = max(max(diagonal, up), max(left, 0));

            H[i * (ref_len + 1) + j] = score;

            // Track maximum
            if (score > max_score) {
                max_score = score;
                max_i = i;
                max_j = j;
            }
        }
    }

    // Write results
    results[gid].score = max_score;
    results[gid].query_start = 0;  // Would calculate from traceback
    results[gid].query_end = max_i;
    results[gid].ref_start = 0;
    results[gid].ref_end = max_j;
}
```

### Rust Host Code (metal-rs wrapper)

```rust
use metal::*;

pub struct MetalSmithWaterman {
    device: Device,
    command_queue: CommandQueue,
    pipeline: ComputePipelineState,
}

impl MetalSmithWaterman {
    pub fn new() -> Result<Self> {
        // Get GPU device
        let device = Device::system_default()
            .ok_or("No Metal-capable device found")?;

        // Create command queue
        let command_queue = device.new_command_queue();

        // Compile Metal shader
        let shader_source = include_str!("smith_waterman.metal");
        let library = device.new_library_with_source(shader_source, &CompileOptions::new())?;

        let kernel = library.get_function("smith_waterman_kernel", None)?;
        let pipeline = device.new_compute_pipeline_state_with_function(&kernel)?;

        Ok(Self {
            device,
            command_queue,
            pipeline,
        })
    }

    pub fn align_batch(
        &self,
        query_seqs: &[&[u8]],
        ref_seqs: &[&[u8]],
        scoring: &ScoringMatrix,
    ) -> Result<Vec<Alignment>> {
        let batch_size = query_seqs.len();

        // Step 1: Pack sequences into contiguous buffers
        let (query_buffer, query_offsets, query_lengths) = pack_sequences(&self.device, query_seqs);
        let (ref_buffer, ref_offsets, ref_lengths) = pack_sequences(&self.device, ref_seqs);

        // Step 2: Create scoring matrix buffer
        let scoring_buffer = self.device.new_buffer_with_data(
            scoring as *const ScoringMatrix as *const c_void,
            mem::size_of::<ScoringMatrix>() as u64,
            MTLResourceOptions::CPUCacheModeDefaultCache,
        );

        // Step 3: Create output buffer
        let results_size = batch_size * mem::size_of::<AlignmentResult>();
        let results_buffer = self.device.new_buffer(
            results_size as u64,
            MTLResourceOptions::StorageModeShared,  // CPU-readable
        );

        // Step 4: Create command buffer
        let command_buffer = self.command_queue.new_command_buffer();
        let encoder = command_buffer.new_compute_command_encoder();

        encoder.set_compute_pipeline_state(&self.pipeline);

        // Bind buffers
        encoder.set_buffer(0, Some(&query_buffer), 0);
        encoder.set_buffer(1, Some(&query_offsets), 0);
        encoder.set_buffer(2, Some(&query_lengths), 0);
        encoder.set_buffer(3, Some(&ref_buffer), 0);
        encoder.set_buffer(4, Some(&ref_offsets), 0);
        encoder.set_buffer(5, Some(&ref_lengths), 0);
        encoder.set_buffer(6, Some(&scoring_buffer), 0);
        encoder.set_buffer(7, Some(&results_buffer), 0);

        // Step 5: Dispatch threads (one per alignment)
        let grid_size = MTLSize {
            width: batch_size as u64,
            height: 1,
            depth: 1,
        };

        let threadgroup_size = MTLSize {
            width: 32,  // Warp size on Apple GPU
            height: 1,
            depth: 1,
        };

        encoder.dispatch_threads(grid_size, threadgroup_size);
        encoder.end_encoding();

        // Step 6: Execute and wait
        command_buffer.commit();
        command_buffer.wait_until_completed();

        // Step 7: Read results
        let results_ptr = results_buffer.contents() as *const AlignmentResult;
        let results = unsafe {
            std::slice::from_raw_parts(results_ptr, batch_size)
        };

        // Convert to Rust Alignment structs
        Ok(results.iter().map(|r| Alignment {
            score: r.score,
            query_start: r.query_start as usize,
            query_end: r.query_end as usize,
            ref_start: r.ref_start as usize,
            ref_end: r.ref_end as usize,
            cigar: vec![],  // Would generate in separate pass
        }).collect())
    }
}

fn pack_sequences(device: &Device, seqs: &[&[u8]]) -> (Buffer, Buffer, Buffer) {
    // Pack all sequences into one contiguous buffer
    let total_len: usize = seqs.iter().map(|s| s.len()).sum();
    let mut packed = Vec::with_capacity(total_len);
    let mut offsets = Vec::with_capacity(seqs.len());
    let mut lengths = Vec::with_capacity(seqs.len());

    let mut offset = 0;
    for seq in seqs {
        offsets.push(offset as u32);
        lengths.push(seq.len() as u32);
        packed.extend_from_slice(seq);
        offset += seq.len();
    }

    // Create Metal buffers (unified memory - zero copy!)
    let data_buffer = device.new_buffer_with_data(
        packed.as_ptr() as *const c_void,
        packed.len() as u64,
        MTLResourceOptions::StorageModeShared,
    );

    let offsets_buffer = device.new_buffer_with_data(
        offsets.as_ptr() as *const c_void,
        (offsets.len() * 4) as u64,
        MTLResourceOptions::StorageModeShared,
    );

    let lengths_buffer = device.new_buffer_with_data(
        lengths.as_ptr() as *const c_void,
        (lengths.len() * 4) as u64,
        MTLResourceOptions::StorageModeShared,
    );

    (data_buffer, offsets_buffer, lengths_buffer)
}
```

**Performance**: ~20,000-50,000 alignments/sec (10-50× speedup)

**Key Advantages**:
- **Unified memory**: Zero-copy transfers (Apple Silicon exclusive)
- **Low dispatch overhead**: 1-3ms (vs 3-5ms on discrete GPUs)
- **Embarrassingly parallel**: No thread synchronization needed
- **Memory efficient**: Each thread uses private DP matrix

---

## Memory Layout

### CPU (Naive & NEON)
```
DP Matrix H:
[                ] ← Row 0 (all zeros)
[  H(1,1) ... ] ← Row 1
[  H(2,1) ... ] ← Row 2
...
[  H(m,1) ... ] ← Row m

Total: (m+1) × (n+1) × sizeof(Cell) bytes
Example: 1000×1000 = 1M cells × 8 bytes = 8 MB per alignment
```

### GPU (Metal)
```
Global Memory:
- query_seqs[]: Concatenated query sequences
- query_offsets[]: Start offset for each query
- query_lengths[]: Length of each query
- ref_seqs[]: Concatenated reference sequences
- ref_offsets[]: Start offset for each reference
- ref_lengths[]: Length of each reference
- results[]: Output alignment results

Thread-Private Memory (each thread):
- H[]: DP matrix (512×512 max = 256KB)
- Local variables (query_base, scores, etc.)

Total GPU memory:
- Input: ~2 MB (1000 sequences × 1000bp × 2)
- Output: ~20 KB (1000 results × 20 bytes)
- Per-thread: ~256 KB × num_active_threads
```

---

## Thread Organization

### GPU Thread Grid

```
Grid: 1D array of threads
├── Thread 0: Aligns query[0] with ref[0]
├── Thread 1: Aligns query[1] with ref[1]
├── Thread 2: Aligns query[2] with ref[2]
...
└── Thread N-1: Aligns query[N-1] with ref[N-1]

Threadgroup Size: 32 (Apple GPU warp size)
- Threadgroup 0: Threads 0-31
- Threadgroup 1: Threads 32-63
- ...

Total Threads: Batch size (e.g., 1000 for 1000 alignments)
```

**Why This Works**:
- Each thread is independent (no communication)
- No synchronization barriers needed
- Maximum GPU occupancy (thousands of threads in flight)

---

## Implementation Notes

### 1. Sequence Length Limits

**CPU**: No practical limit (limited by RAM)
- 1000bp × 1000bp = 8 MB (manageable)
- 10,000bp × 10,000bp = 800 MB (still ok)

**GPU**: Limited by thread-private memory
- **Max sequence length**: 512bp (512×512 = 256KB per thread)
- For longer sequences: Use tiled algorithm or fallback to CPU

### 2. Batch Size Considerations

**Small Batch** (<100 alignments):
- GPU dispatch overhead dominates
- **Use CPU** (faster for small batches)

**Medium Batch** (100-1000 alignments):
- GPU starts to win
- **Breakeven**: ~200-500 alignments

**Large Batch** (>1000 alignments):
- GPU significantly faster (10-50×)
- **Sweet spot**: 1000-10,000 alignments per dispatch

### 3. Memory Access Patterns

**CPU**:
- Sequential row access (cache-friendly)
- Prefetching helps

**NEON**:
- Diagonal stripes (less cache-friendly)
- Limited speedup due to irregular access

**GPU**:
- Each thread accesses private memory (no contention)
- Coalesced global memory reads (packed sequences)

### 4. Traceback Optimization

**Current**: Store only scores (no traceback)
- Faster forward pass
- Must recompute for CIGAR

**Future**: Store directions
- Slower forward pass
- Instant CIGAR generation

**Hybrid**: Store scores, recompute traceback on CPU
- GPU does expensive forward pass
- CPU does cheap traceback

---

## Testing Strategy

### Unit Tests
```rust
#[test]
fn test_known_alignment() {
    let query = b"ACGT";
    let reference = b"ACGT";
    let scoring = ScoringMatrix::default();

    let result = smith_waterman_naive(query, reference, &scoring);

    assert_eq!(result.score, 8);  // 4 matches × 2 = 8
    assert_eq!(result.cigar, vec![CigarOp::Match(4)]);
}

#[test]
fn test_mismatch() {
    let query = b"AAAA";
    let reference = b"TTTT";
    let scoring = ScoringMatrix::default();

    let result = smith_waterman_naive(query, reference, &scoring);

    assert_eq!(result.score, 0);  // No good alignment
}
```

### Property Tests
```rust
proptest! {
    #[test]
    fn test_gpu_matches_cpu(
        query in prop::collection::vec(prop::sample::select(b"ACGT"), 10..100),
        reference in prop::collection::vec(prop::sample::select(b"ACGT"), 10..100)
    ) {
        let scoring = ScoringMatrix::default();

        let cpu_result = smith_waterman_naive(&query, &reference, &scoring);
        let gpu_result = smith_waterman_gpu_single(&query, &reference, &scoring);

        prop_assert_eq!(cpu_result.score, gpu_result.score);
    }
}
```

---

## Next Steps

1. **Implement CPU Naive** (baseline correctness)
2. **Implement NEON** (ARM portable optimization)
3. **Implement Metal GPU** (Apple Silicon exclusive)
4. **Add Tests** (unit + property tests)
5. **Benchmark** (N=30 statistical rigor)
6. **Analyze Results** (achieve 10-50× GPU speedup?)

---

**Status**: ✅ Algorithm design complete
**Next**: Begin implementation (start with CPU naive)
**Date**: November 13, 2025
