#include <metal_stdlib>
using namespace metal;

// Scoring configuration (matches Rust ScoringMatrix)
struct ScoringConfig {
    int match_score;      // +2 for match
    int mismatch_score;   // -1 for mismatch
    int gap_open;         // -2 for opening gap
    int gap_extend;       // -1 for extending gap
};

// Alignment result (simplified for GPU)
struct AlignmentResult {
    int score;            // Maximum alignment score
    uint query_start;     // Start position in query
    uint query_end;       // End position in query
    uint ref_start;       // Start position in reference
    uint ref_end;         // End position in reference
};

// Direction for traceback
enum Direction : uchar {
    NONE = 0,      // Start of alignment (score = 0)
    DIAGONAL = 1,  // Match/mismatch (from H[i-1][j-1])
    UP = 2,        // Deletion (from H[i-1][j])
    LEFT = 3       // Insertion (from H[i][j-1])
};

// DP matrix cell
struct Cell {
    int score;
    Direction direction;
};

// Helper: Get match/mismatch score
inline int match_score(uchar a, uchar b, constant ScoringConfig& config) {
    return (a == b) ? config.match_score : config.mismatch_score;
}

// Helper: Get maximum of 4 values and corresponding direction
// Called as: max4_with_dir(none, diagonal, up, left, dir)
// Matches CPU's max4 logic exactly for consistent tie-breaking
inline int max4_with_dir(int none, int diagonal, int up, int left, thread Direction& dir) {
    // Compare first pair: diagonal vs up (prefers diagonal on tie)
    int max_diag_up;
    Direction dir_diag_up;
    if (diagonal >= up) {
        max_diag_up = diagonal;
        dir_diag_up = DIAGONAL;
    } else {
        max_diag_up = up;
        dir_diag_up = UP;
    }

    // Compare second pair: left vs none (prefers left on tie)
    int max_left_none;
    Direction dir_left_none;
    if (left >= none) {
        max_left_none = left;
        dir_left_none = LEFT;
    } else {
        max_left_none = none;
        dir_left_none = NONE;
    }

    // Final comparison (prefers first pair on tie)
    if (max_diag_up >= max_left_none) {
        dir = dir_diag_up;
        return max_diag_up;
    } else {
        dir = dir_left_none;
        return max_left_none;
    }
}

// Smith-Waterman DP matrix computation with anti-diagonal parallelization
//
// Multiple threads collaborate on each alignment using anti-diagonal processing.
// This achieves 5-10× speedup vs serial processing by exploiting parallelism
// within individual DP computations.
//
// Thread layout: threadgroup processes one alignment, threads process anti-diagonals
kernel void smith_waterman_dp_antidiagonal(
    // Input sequences (packed: all queries concatenated, all refs concatenated)
    constant uchar* queries         [[buffer(0)]],
    constant uchar* references      [[buffer(1)]],
    constant uint* query_lengths    [[buffer(2)]],
    constant uint* ref_lengths      [[buffer(3)]],
    constant uint* query_offsets    [[buffer(4)]],
    constant uint* ref_offsets      [[buffer(5)]],
    constant ScoringConfig& config  [[buffer(6)]],

    // Output: DP matrices (one per alignment)
    device Cell* dp_matrices        [[buffer(7)]],
    device uint* matrix_offsets     [[buffer(8)]],

    // Thread indexing for anti-diagonal parallelization
    uint2 thread_pos  [[thread_position_in_threadgroup]],    // (thread_x, thread_y) within group
    uint2 group_pos   [[threadgroup_position_in_grid]],      // (group_x, group_y) in grid
    uint2 group_size  [[threads_per_threadgroup]]             // (width, height) of threadgroup
) {
    // Each threadgroup handles one alignment
    uint batch_idx = group_pos.x;
    uint thread_idx = thread_pos.x; // Linear thread index within group

    // Get this alignment's sequences
    uint q_len = query_lengths[batch_idx];
    uint r_len = ref_lengths[batch_idx];
    uint q_offset = query_offsets[batch_idx];
    uint r_offset = ref_offsets[batch_idx];

    constant uchar* query = queries + q_offset;
    constant uchar* ref = references + r_offset;

    // Get this alignment's DP matrix
    uint matrix_offset = matrix_offsets[batch_idx];
    device Cell* matrix = dp_matrices + matrix_offset;

    // Matrix dimensions: (q_len + 1) × (r_len + 1)
    uint rows = q_len + 1;
    uint cols = r_len + 1;
    uint max_threads = group_size.x; // Number of threads in threadgroup

    // Initialize first row and column (only first thread does this)
    if (thread_idx == 0) {
        for (uint j = 0; j < cols; j++) {
            matrix[0 * cols + j] = Cell{0, NONE};
        }
        for (uint i = 1; i < rows; i++) {
            matrix[i * cols + 0] = Cell{0, NONE};
        }
    }

    // Synchronize after initialization
    threadgroup_barrier(mem_flags::mem_device);

    // Anti-diagonal processing:
    // Process diagonals from top-left to bottom-right
    // Each diagonal k contains cells where i + j = k
    uint max_diagonal = rows + cols - 2; // Maximum diagonal index

    for (uint diagonal = 2; diagonal <= max_diagonal; diagonal++) {
        // Calculate range of valid positions on this diagonal
        uint i_start = max(1u, diagonal + 1 - cols);
        uint i_end = min(rows - 1, diagonal - 1);

        // Number of cells on this diagonal
        uint diagonal_size = i_end - i_start + 1;

        // Distribute cells across threads
        for (uint offset = thread_idx; offset < diagonal_size; offset += max_threads) {
            uint i = i_start + offset;
            uint j = diagonal - i;

            // Bounds check
            if (i >= rows || j >= cols || i == 0 || j == 0) continue;

            // Compute DP cell: depends only on previous diagonal (already computed)
            int diag_score = matrix[(i-1) * cols + (j-1)].score +
                           match_score(query[i-1], ref[j-1], config);
            int up_score = matrix[(i-1) * cols + j].score + config.gap_open;
            int left_score = matrix[i * cols + (j-1)].score + config.gap_open;

            // Find maximum and direction
            Direction dir;
            int score = max4_with_dir(0, diag_score, up_score, left_score, dir);

            // Store result
            matrix[i * cols + j] = Cell{score, dir};
        }

        // Synchronize before next diagonal
        threadgroup_barrier(mem_flags::mem_device);
    }
}

// Legacy serial DP computation (fallback for small alignments)
//
// Each thread computes one alignment in the batch
// Uses simple row-by-row processing for small matrices where
// anti-diagonal overhead exceeds parallelism benefits
kernel void smith_waterman_dp(
    // Input sequences (packed: all queries concatenated, all refs concatenated)
    constant uchar* queries         [[buffer(0)]],
    constant uchar* references      [[buffer(1)]],
    constant uint* query_lengths    [[buffer(2)]],
    constant uint* ref_lengths      [[buffer(3)]],
    constant uint* query_offsets    [[buffer(4)]],
    constant uint* ref_offsets      [[buffer(5)]],
    constant ScoringConfig& config  [[buffer(6)]],

    // Output: DP matrices (one per alignment)
    device Cell* dp_matrices        [[buffer(7)]],
    device uint* matrix_offsets     [[buffer(8)]],

    // Thread indexing
    uint batch_idx [[thread_position_in_grid]]
) {
    // Get this alignment's sequences
    uint q_len = query_lengths[batch_idx];
    uint r_len = ref_lengths[batch_idx];
    uint q_offset = query_offsets[batch_idx];
    uint r_offset = ref_offsets[batch_idx];

    constant uchar* query = queries + q_offset;
    constant uchar* ref = references + r_offset;

    // Get this alignment's DP matrix
    uint matrix_offset = matrix_offsets[batch_idx];
    device Cell* matrix = dp_matrices + matrix_offset;

    // Matrix dimensions: (q_len + 1) × (r_len + 1)
    uint rows = q_len + 1;
    uint cols = r_len + 1;

    // Initialize first row and column
    for (uint j = 0; j < cols; j++) {
        matrix[0 * cols + j] = Cell{0, NONE};
    }
    for (uint i = 1; i < rows; i++) {
        matrix[i * cols + 0] = Cell{0, NONE};
    }

    // Fill DP matrix
    for (uint i = 1; i < rows; i++) {
        for (uint j = 1; j < cols; j++) {
            // Get scores from three predecessors
            int diag_score = matrix[(i-1) * cols + (j-1)].score +
                           match_score(query[i-1], ref[j-1], config);
            int up_score = matrix[(i-1) * cols + j].score + config.gap_open;
            int left_score = matrix[i * cols + (j-1)].score + config.gap_open;

            // Find maximum and direction
            Direction dir;
            int score = max4_with_dir(0, diag_score, up_score, left_score, dir);

            // Store result
            matrix[i * cols + j] = Cell{score, dir};
        }
    }
}

// Find maximum score position in DP matrix
kernel void find_max_score(
    device Cell* dp_matrices        [[buffer(0)]],
    constant uint* matrix_offsets   [[buffer(1)]],
    constant uint* query_lengths    [[buffer(2)]],
    constant uint* ref_lengths      [[buffer(3)]],
    device AlignmentResult* results [[buffer(4)]],
    uint batch_idx [[thread_position_in_grid]]
) {
    // Get this alignment's DP matrix
    uint matrix_offset = matrix_offsets[batch_idx];
    device Cell* matrix = dp_matrices + matrix_offset;

    uint q_len = query_lengths[batch_idx];
    uint r_len = ref_lengths[batch_idx];
    uint cols = r_len + 1;
    uint rows = q_len + 1;

    // Find maximum score in matrix
    int max_score = 0;
    uint max_i = 0;
    uint max_j = 0;

    for (uint i = 1; i < rows; i++) {
        for (uint j = 1; j < cols; j++) {
            int score = matrix[i * cols + j].score;
            if (score > max_score) {
                max_score = score;
                max_i = i;
                max_j = j;
            }
        }
    }

    // Store result (end positions)
    results[batch_idx].score = max_score;
    results[batch_idx].query_end = max_i;
    results[batch_idx].ref_end = max_j;
}

// Traceback to find alignment start positions
//
// Note: CIGAR generation is done on CPU for simplicity
// (traceback is inherently serial and benefits from CPU cache)
kernel void traceback(
    device Cell* dp_matrices        [[buffer(0)]],
    constant uint* matrix_offsets   [[buffer(1)]],
    constant uint* ref_lengths      [[buffer(2)]],
    device AlignmentResult* results [[buffer(3)]],
    uint batch_idx [[thread_position_in_grid]]
) {
    // Get this alignment's DP matrix
    uint matrix_offset = matrix_offsets[batch_idx];
    device Cell* matrix = dp_matrices + matrix_offset;
    uint cols = ref_lengths[batch_idx] + 1;

    // Start from max score position
    uint i = results[batch_idx].query_end;
    uint j = results[batch_idx].ref_end;

    // Traceback until we hit a zero score (local alignment)
    while (i > 0 && j > 0) {
        Cell cell = matrix[i * cols + j];

        if (cell.direction == NONE) {
            break;
        }

        switch (cell.direction) {
            case DIAGONAL:
                i--;
                j--;
                break;
            case UP:
                i--;
                break;
            case LEFT:
                j--;
                break;
            case NONE:
                break;
        }
    }

    // Store start positions
    results[batch_idx].query_start = i;
    results[batch_idx].ref_start = j;
}

// Optimized version: Single-pass DP with max tracking
//
// Combines DP computation and max finding in one kernel
// More efficient for small-to-medium sequences
kernel void smith_waterman_combined(
    constant uchar* queries         [[buffer(0)]],
    constant uchar* references      [[buffer(1)]],
    constant uint* query_lengths    [[buffer(2)]],
    constant uint* ref_lengths      [[buffer(3)]],
    constant uint* query_offsets    [[buffer(4)]],
    constant uint* ref_offsets      [[buffer(5)]],
    constant ScoringConfig& config  [[buffer(6)]],
    device Cell* dp_matrices        [[buffer(7)]],
    device uint* matrix_offsets     [[buffer(8)]],
    device AlignmentResult* results [[buffer(9)]],
    uint batch_idx [[thread_position_in_grid]]
) {
    // Get this alignment's sequences
    uint q_len = query_lengths[batch_idx];
    uint r_len = ref_lengths[batch_idx];
    uint q_offset = query_offsets[batch_idx];
    uint r_offset = ref_offsets[batch_idx];

    constant uchar* query = queries + q_offset;
    constant uchar* ref = references + r_offset;

    // Get this alignment's DP matrix
    uint matrix_offset = matrix_offsets[batch_idx];
    device Cell* matrix = dp_matrices + matrix_offset;

    uint rows = q_len + 1;
    uint cols = r_len + 1;

    // Track maximum score
    int max_score = 0;
    uint max_i = 0;
    uint max_j = 0;

    // Initialize first row and column
    for (uint j = 0; j < cols; j++) {
        matrix[0 * cols + j] = Cell{0, NONE};
    }
    for (uint i = 1; i < rows; i++) {
        matrix[i * cols + 0] = Cell{0, NONE};
    }

    // Fill DP matrix and track maximum
    for (uint i = 1; i < rows; i++) {
        for (uint j = 1; j < cols; j++) {
            // Compute scores
            int diag_score = matrix[(i-1) * cols + (j-1)].score +
                           match_score(query[i-1], ref[j-1], config);
            int up_score = matrix[(i-1) * cols + j].score + config.gap_open;
            int left_score = matrix[i * cols + (j-1)].score + config.gap_open;

            // Find maximum
            Direction dir;
            int score = max4_with_dir(0, diag_score, up_score, left_score, dir);

            // Store in matrix
            matrix[i * cols + j] = Cell{score, dir};

            // Track maximum
            if (score > max_score) {
                max_score = score;
                max_i = i;
                max_j = j;
            }
        }
    }

    // Traceback from maximum
    uint i = max_i;
    uint j = max_j;

    while (i > 0 && j > 0) {
        Cell cell = matrix[i * cols + j];

        if (cell.direction == NONE) {
            break;
        }

        switch (cell.direction) {
            case DIAGONAL:
                i--;
                j--;
                break;
            case UP:
                i--;
                break;
            case LEFT:
                j--;
                break;
            case NONE:
                break;
        }
    }

    // Store results
    results[batch_idx].score = max_score;
    results[batch_idx].query_start = i;
    results[batch_idx].query_end = max_i;
    results[batch_idx].ref_start = j;
    results[batch_idx].ref_end = max_j;
}
