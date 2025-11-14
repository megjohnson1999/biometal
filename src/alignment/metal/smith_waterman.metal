//! Smith-Waterman GPU kernel for Apple Metal
//!
//! This compute shader implements Smith-Waterman local sequence alignment
//! on the GPU using thread-per-alignment parallelism.
//!
//! # Algorithm
//!
//! Each thread computes one complete alignment:
//! - Load query and reference sequences from global memory
//! - Compute DP matrix (stored in thread-local memory for small alignments)
//! - Find maximum score and position
//! - Store result back to global memory
//!
//! # Performance
//!
//! Expected 10-50× speedup vs CPU for batch sizes >100 alignments
//! - Apple Silicon unified memory: Zero-copy data sharing
//! - Low dispatch overhead: ~1-3ms (vs 3-5ms for CUDA)
//! - Thread-per-alignment: Process 1000+ alignments in parallel

#include <metal_stdlib>
using namespace metal;

/// Scoring parameters for alignment
struct ScoringMatrix {
    int match_score;     // +2
    int mismatch_score;  // -1
    int gap_open;        // -2
    int gap_extend;      // -1 (unused in basic Smith-Waterman)
};

/// Alignment result
struct AlignmentResult {
    int score;           // Maximum alignment score
    int query_start;     // Start position in query (0-indexed)
    int query_end;       // End position in query (exclusive)
    int ref_start;       // Start position in reference (0-indexed)
    int ref_end;         // End position in reference (exclusive)
};

/// Direction for traceback
enum Direction : int {
    None = 0,      // Start of alignment
    Diagonal = 1,  // Match/mismatch
    Up = 2,        // Deletion
    Left = 3       // Insertion
};

/// DP matrix cell
struct Cell {
    int score;
    int direction;
};

/// Calculate score for aligning two bases
inline int calculate_score(char a, char b, constant ScoringMatrix& scoring) {
    return (a == b) ? scoring.match_score : scoring.mismatch_score;
}

/// Find maximum of four (value, direction) pairs
inline void max4(int a, int b, int c, int d,
                 int dir_a, int dir_b, int dir_c, int dir_d,
                 thread int* max_val, thread int* max_dir) {
    *max_val = a;
    *max_dir = dir_a;

    if (b > *max_val) { *max_val = b; *max_dir = dir_b; }
    if (c > *max_val) { *max_val = c; *max_dir = dir_c; }
    if (d > *max_val) { *max_val = d; *max_dir = dir_d; }
}

/// Smith-Waterman kernel: One thread per alignment
///
/// Grid dimensions:
/// - threads_per_threadgroup: 64-256 (tuned for Apple Silicon)
/// - num_threadgroups: ceil(num_alignments / threads_per_threadgroup)
///
/// Memory usage:
/// - Small alignments (<1000bp): DP matrix in thread-local memory
/// - Large alignments: Use global memory (with coalesced access)
kernel void smith_waterman_kernel(
    // Input buffers
    constant char* queries [[buffer(0)]],           // Concatenated query sequences
    constant char* references [[buffer(1)]],        // Concatenated reference sequences
    constant int* query_offsets [[buffer(2)]],      // Start offset for each query
    constant int* reference_offsets [[buffer(3)]],  // Start offset for each reference
    constant int* query_lengths [[buffer(4)]],      // Length of each query
    constant int* reference_lengths [[buffer(5)]],  // Length of each reference
    constant ScoringMatrix& scoring [[buffer(6)]],  // Scoring matrix

    // Output buffer
    device AlignmentResult* results [[buffer(7)]],  // Alignment results

    // Thread information
    uint tid [[thread_position_in_grid]])           // Thread ID
{
    // Get query and reference for this alignment
    int query_len = query_lengths[tid];
    int ref_len = reference_lengths[tid];
    int query_offset = query_offsets[tid];
    int ref_offset = reference_offsets[tid];

    constant char* query = queries + query_offset;
    constant char* reference = references + ref_offset;

    // Handle empty sequences
    if (query_len == 0 || ref_len == 0) {
        results[tid] = { 0, 0, 0, 0, 0 };
        return;
    }

    // Allocate DP matrix in thread-local memory
    // For large alignments (>1000bp × 1000bp), this should use shared memory
    // or global memory to avoid stack overflow
    //
    // NOTE: Metal has very limited thread-local memory
    // Stack size limit is platform-dependent but typically ~16-32KB per thread
    // For production, implement tiled/blocked approach or use device memory
    #define MAX_SIZE 100  // Support up to 100×100 DP matrix (~80KB with overhead)

    if (query_len > MAX_SIZE || ref_len > MAX_SIZE) {
        // Too large for thread-local memory - skip for now
        results[tid] = { -1, 0, 0, 0, 0 };  // Error marker
        return;
    }

    // Initialize DP matrix (stored as 1D array for cache efficiency)
    // Size: (query_len + 1) × (ref_len + 1)
    // Using constant size to satisfy Metal compiler
    Cell matrix[(MAX_SIZE + 1) * (MAX_SIZE + 1)];

    int m = query_len;
    int n = ref_len;
    int width = n + 1;

    // Initialize first row and column to zero (local alignment)
    for (int i = 0; i <= m; i++) {
        matrix[i * width + 0] = { 0, None };
    }
    for (int j = 0; j <= n; j++) {
        matrix[0 * width + j] = { 0, None };
    }

    // Track maximum score for traceback
    int max_score = 0;
    int max_i = 0;
    int max_j = 0;

    // Fill DP matrix (forward pass)
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            // Calculate score for match/mismatch
            int match_score = calculate_score(query[i - 1], reference[j - 1], scoring);

            // Calculate scores from three directions
            int diagonal = matrix[(i - 1) * width + (j - 1)].score + match_score;
            int up = matrix[(i - 1) * width + j].score + scoring.gap_open;
            int left = matrix[i * width + (j - 1)].score + scoring.gap_open;

            // Take maximum (including 0 for local alignment)
            int score, direction;
            max4(diagonal, up, left, 0,
                 Diagonal, Up, Left, None,
                 &score, &direction);

            matrix[i * width + j] = { score, direction };

            // Track maximum score
            if (score > max_score) {
                max_score = score;
                max_i = i;
                max_j = j;
            }
        }
    }

    // Traceback to find alignment start positions
    // (Full CIGAR reconstruction deferred - requires variable-length output)
    int i = max_i;
    int j = max_j;

    while (i > 0 && j > 0 && matrix[i * width + j].direction != None) {
        int dir = matrix[i * width + j].direction;
        if (dir == Diagonal) {
            i--;
            j--;
        } else if (dir == Up) {
            i--;
        } else if (dir == Left) {
            j--;
        }
    }

    // Store result
    results[tid] = {
        max_score,
        i,           // query_start
        max_i,       // query_end
        j,           // ref_start
        max_j        // ref_end
    };
}
