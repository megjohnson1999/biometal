//! GPU-accelerated batch Smith-Waterman alignment using Metal
//!
//! This module provides GPU-accelerated Smith-Waterman alignment for batches
//! of sequence pairs using Apple's Metal framework.
//!
//! # Performance
//!
//! GPU acceleration is most effective for:
//! - Batch sizes ≥256 alignments (amortizes GPU dispatch overhead)
//! - Medium-to-long sequences (≥100bp)
//! - Multiple sequence alignment (MSA) workflows
//! - Database searches
//!
//! For small batches (<16) or very short sequences, CPU implementation
//! may be faster due to GPU dispatch overhead.

use crate::alignment::{Alignment, CigarOp, ScoringMatrix};
use metal::*;
use std::mem;

/// Scoring configuration for Metal shader (C-compatible)
#[repr(C)]
struct ScoringConfig {
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
}

impl From<&ScoringMatrix> for ScoringConfig {
    fn from(matrix: &ScoringMatrix) -> Self {
        Self {
            match_score: matrix.match_score,
            mismatch_score: matrix.mismatch_score,
            gap_open: matrix.gap_open,
            gap_extend: matrix.gap_extend,
        }
    }
}

/// Alignment result from Metal (C-compatible)
#[repr(C)]
#[derive(Debug, Clone, Copy)]
struct AlignmentResult {
    score: i32,
    query_start: u32,
    query_end: u32,
    ref_start: u32,
    ref_end: u32,
}

/// DP matrix cell from Metal (C-compatible)
#[repr(C)]
#[derive(Debug, Clone, Copy)]
struct Cell {
    score: i32,
    direction: u8,
}

/// Direction enum (matches Metal shader)
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Direction {
    None = 0,
    Diagonal = 1,
    Up = 2,
    Left = 3,
}

/// GPU-accelerated batch Smith-Waterman aligner
///
/// This struct manages Metal resources and provides batch alignment functionality.
///
/// # Example
///
/// ```no_run
/// use biometal::alignment::gpu::GpuAlignmentBatch;
/// use biometal::alignment::ScoringMatrix;
///
/// let mut gpu_batch = GpuAlignmentBatch::new().unwrap();
///
/// let queries = vec![b"ACGTACGT", b"TGCATGCA"];
/// let targets = vec![b"ACGGACGG", b"TGCATGCA"];
/// let query_refs: Vec<&[u8]> = queries.iter().map(|q| q.as_slice()).collect();
/// let target_refs: Vec<&[u8]> = targets.iter().map(|t| t.as_slice()).collect();
///
/// let scoring = ScoringMatrix::default();
/// let alignments = gpu_batch.align_batch(&query_refs, &target_refs, &scoring).unwrap();
/// ```
pub struct GpuAlignmentBatch {
    device: Device,
    command_queue: CommandQueue,
    pipeline_combined: ComputePipelineState,
    pipeline_antidiagonal: ComputePipelineState,
    pipeline_find_max: ComputePipelineState,
}

impl GpuAlignmentBatch {
    /// Create a new GPU batch alignment processor
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Metal device is not available (non-macOS platforms)
    /// - Shader compilation fails
    /// - Pipeline creation fails
    pub fn new() -> Result<Self, String> {
        // Get Metal device
        let device = Device::system_default()
            .ok_or_else(|| "Metal device not available (macOS only)".to_string())?;

        // Create command queue
        let command_queue = device.new_command_queue();

        // Load and compile Metal shader
        let shader_source = include_str!("shaders/smith_waterman.metal");
        let library = device
            .new_library_with_source(shader_source, &CompileOptions::new())
            .map_err(|e| format!("Failed to compile Metal shader: {:?}", e))?;

        // Create compute pipeline for combined kernel (legacy)
        let function_combined = library
            .get_function("smith_waterman_combined", None)
            .map_err(|e| format!("Failed to get combined Metal function: {:?}", e))?;

        let pipeline_combined = device
            .new_compute_pipeline_state_with_function(&function_combined)
            .map_err(|e| format!("Failed to create combined pipeline: {:?}", e))?;

        // Create compute pipeline for anti-diagonal DP
        let function_antidiagonal = library
            .get_function("smith_waterman_dp_antidiagonal", None)
            .map_err(|e| format!("Failed to get anti-diagonal Metal function: {:?}", e))?;

        let pipeline_antidiagonal = device
            .new_compute_pipeline_state_with_function(&function_antidiagonal)
            .map_err(|e| format!("Failed to create anti-diagonal pipeline: {:?}", e))?;

        // Create compute pipeline for finding max score
        let function_find_max = library
            .get_function("find_max_score", None)
            .map_err(|e| format!("Failed to get find_max Metal function: {:?}", e))?;

        let pipeline_find_max = device
            .new_compute_pipeline_state_with_function(&function_find_max)
            .map_err(|e| format!("Failed to create find_max pipeline: {:?}", e))?;

        Ok(Self {
            device,
            command_queue,
            pipeline_combined,
            pipeline_antidiagonal,
            pipeline_find_max,
        })
    }

    /// Align a batch of query-target pairs
    ///
    /// # Arguments
    ///
    /// * `queries` - Slice of query sequences
    /// * `targets` - Slice of target sequences (must be same length as queries)
    /// * `scoring` - Scoring matrix for alignment
    ///
    /// # Returns
    ///
    /// Vector of `Alignment` results, one per query-target pair
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Query and target slices have different lengths
    /// - GPU buffer allocation fails
    /// - GPU computation fails
    pub fn align_batch(
        &self,
        queries: &[&[u8]],
        targets: &[&[u8]],
        scoring: &ScoringMatrix,
    ) -> Result<Vec<Alignment>, String> {
        if queries.len() != targets.len() {
            return Err(format!(
                "Query and target counts must match: {} vs {}",
                queries.len(),
                targets.len()
            ));
        }

        let batch_size = queries.len();
        if batch_size == 0 {
            return Ok(Vec::new());
        }

        // Pack sequences and compute offsets
        let (query_data, query_lengths, query_offsets) = Self::pack_sequences(queries);
        let (target_data, target_lengths, target_offsets) = Self::pack_sequences(targets);

        // Compute matrix sizes and offsets
        let (matrix_sizes, matrix_offsets) = Self::compute_matrix_layout(&query_lengths, &target_lengths);
        let total_matrix_size: usize = matrix_sizes.iter().sum();

        // Create Metal buffers
        let query_buffer = self.create_buffer(&query_data)?;
        let target_buffer = self.create_buffer(&target_data)?;
        let query_len_buffer = self.create_buffer(&query_lengths)?;
        let target_len_buffer = self.create_buffer(&target_lengths)?;
        let query_offset_buffer = self.create_buffer(&query_offsets)?;
        let target_offset_buffer = self.create_buffer(&target_offsets)?;

        let scoring_config = ScoringConfig::from(scoring);
        let scoring_buffer = self.create_buffer(&[scoring_config])?;

        // Allocate DP matrices
        let dp_matrices: Vec<Cell> = vec![Cell { score: 0, direction: 0 }; total_matrix_size];
        let dp_buffer = self.create_buffer(&dp_matrices)?;
        let matrix_offset_buffer = self.create_buffer(&matrix_offsets)?;

        // Allocate results
        let results: Vec<AlignmentResult> = vec![
            AlignmentResult {
                score: 0,
                query_start: 0,
                query_end: 0,
                ref_start: 0,
                ref_end: 0
            };
            batch_size
        ];
        let results_buffer = self.create_buffer(&results)?;

        // Choose algorithm based on sequence sizes and batch size
        // Anti-diagonal provides better speedup for medium-to-large alignments
        let use_antidiagonal = Self::should_use_antidiagonal(&query_lengths, &target_lengths, batch_size);

        if use_antidiagonal {
            self.execute_antidiagonal_pipeline(
                &query_buffer, &target_buffer, &query_len_buffer, &target_len_buffer,
                &query_offset_buffer, &target_offset_buffer, &scoring_buffer,
                &dp_buffer, &matrix_offset_buffer, &results_buffer,
                batch_size
            )?;
        } else {
            // Fallback to legacy combined kernel for small alignments
            self.execute_combined_pipeline(
                &query_buffer, &target_buffer, &query_len_buffer, &target_len_buffer,
                &query_offset_buffer, &target_offset_buffer, &scoring_buffer,
                &dp_buffer, &matrix_offset_buffer, &results_buffer,
                batch_size
            )?;
        }

        // Read back results
        let results_ptr = results_buffer.contents() as *const AlignmentResult;
        let gpu_results = unsafe { std::slice::from_raw_parts(results_ptr, batch_size) };

        // Read back DP matrices for CIGAR generation
        let dp_ptr = dp_buffer.contents() as *const Cell;
        let dp_data = unsafe { std::slice::from_raw_parts(dp_ptr, total_matrix_size) };

        // Convert GPU results to Alignment structs with CIGAR
        let mut alignments = Vec::with_capacity(batch_size);
        for i in 0..batch_size {
            let result = &gpu_results[i];
            let matrix_offset = matrix_offsets[i] as usize;
            let cols = (target_lengths[i] + 1) as usize;

            // Extract this alignment's DP matrix
            let matrix_size = matrix_sizes[i];
            let matrix = &dp_data[matrix_offset..matrix_offset + matrix_size];

            // Generate CIGAR string from traceback
            let cigar = Self::generate_cigar(
                matrix,
                cols,
                result.query_start as usize,
                result.query_end as usize,
                result.ref_start as usize,
                result.ref_end as usize,
            );

            alignments.push(Alignment {
                score: result.score,
                query_start: result.query_start as usize,
                query_end: result.query_end as usize,
                ref_start: result.ref_start as usize,
                ref_end: result.ref_end as usize,
                cigar,
            });
        }

        Ok(alignments)
    }

    /// Pack sequences into contiguous buffer with lengths and offsets
    fn pack_sequences(sequences: &[&[u8]]) -> (Vec<u8>, Vec<u32>, Vec<u32>) {
        let total_length: usize = sequences.iter().map(|s| s.len()).sum();
        let mut data = Vec::with_capacity(total_length);
        let mut lengths = Vec::with_capacity(sequences.len());
        let mut offsets = Vec::with_capacity(sequences.len());

        let mut current_offset = 0u32;
        for seq in sequences {
            offsets.push(current_offset);
            lengths.push(seq.len() as u32);
            data.extend_from_slice(seq);
            current_offset += seq.len() as u32;
        }

        (data, lengths, offsets)
    }

    /// Compute DP matrix sizes and offsets for batch
    fn compute_matrix_layout(query_lengths: &[u32], target_lengths: &[u32]) -> (Vec<usize>, Vec<u32>) {
        let batch_size = query_lengths.len();
        let mut sizes = Vec::with_capacity(batch_size);
        let mut offsets = Vec::with_capacity(batch_size);

        let mut current_offset = 0u32;
        for i in 0..batch_size {
            let rows = (query_lengths[i] + 1) as usize;
            let cols = (target_lengths[i] + 1) as usize;
            let size = rows * cols;

            offsets.push(current_offset);
            sizes.push(size);
            current_offset += size as u32;
        }

        (sizes, offsets)
    }

    /// Create Metal buffer from data
    fn create_buffer<T>(&self, data: &[T]) -> Result<Buffer, String> {
        let size = (data.len() * mem::size_of::<T>()) as u64;
        let buffer = self.device.new_buffer(size, MTLResourceOptions::StorageModeShared);

        unsafe {
            let ptr = buffer.contents() as *mut T;
            std::ptr::copy_nonoverlapping(data.as_ptr(), ptr, data.len());
        }

        Ok(buffer)
    }

    /// Generate CIGAR string from DP matrix traceback
    fn generate_cigar(
        matrix: &[Cell],
        cols: usize,
        query_start: usize,
        query_end: usize,
        ref_start: usize,
        ref_end: usize,
    ) -> Vec<CigarOp> {
        let mut cigar = Vec::new();
        let mut i = query_end;
        let mut j = ref_end;

        let mut current_op: Option<CigarOp> = None;

        while i > query_start || j > ref_start {
            if i == 0 || j == 0 {
                break;
            }

            let cell = matrix[i * cols + j];
            let direction = match cell.direction {
                1 => Direction::Diagonal,
                2 => Direction::Up,
                3 => Direction::Left,
                _ => Direction::None,
            };

            let op = match direction {
                Direction::Diagonal => {
                    i -= 1;
                    j -= 1;
                    CigarOp::Match(1)
                }
                Direction::Up => {
                    i -= 1;
                    CigarOp::Deletion(1)
                }
                Direction::Left => {
                    j -= 1;
                    CigarOp::Insertion(1)
                }
                Direction::None => break,
            };

            // Combine consecutive operations of the same type
            match (&mut current_op, &op) {
                (Some(CigarOp::Match(ref mut count)), CigarOp::Match(_)) => *count += 1,
                (Some(CigarOp::Insertion(ref mut count)), CigarOp::Insertion(_)) => *count += 1,
                (Some(CigarOp::Deletion(ref mut count)), CigarOp::Deletion(_)) => *count += 1,
                _ => {
                    if let Some(prev_op) = current_op.take() {
                        cigar.push(prev_op);
                    }
                    current_op = Some(op);
                }
            }
        }

        if let Some(op) = current_op {
            cigar.push(op);
        }

        cigar.reverse();
        cigar
    }

    /// Determine whether to use anti-diagonal algorithm based on sequence characteristics
    ///
    /// Anti-diagonal provides better speedup for medium-to-large alignments where
    /// the parallelism benefits outweigh threadgroup coordination overhead.
    fn should_use_antidiagonal(
        query_lengths: &[u32],
        target_lengths: &[u32],
        _batch_size: usize
    ) -> bool {
        // Use anti-diagonal for alignments with sufficient work per threadgroup
        // Threshold: sequences ≥64bp provide enough parallelism to justify overhead
        let min_size_threshold = 64;

        // Check if most alignments benefit from anti-diagonal processing
        let large_alignment_count = query_lengths.iter()
            .zip(target_lengths.iter())
            .filter(|(&q_len, &t_len)| q_len >= min_size_threshold && t_len >= min_size_threshold)
            .count();

        // Use anti-diagonal if ≥50% of alignments are large enough to benefit
        large_alignment_count * 2 >= query_lengths.len()
    }

    /// Execute anti-diagonal parallelization pipeline (5-10× speedup)
    ///
    /// Uses multi-kernel approach: DP → find_max → traceback
    /// Each threadgroup handles one alignment with multiple threads per alignment
    #[allow(clippy::too_many_arguments)]
    fn execute_antidiagonal_pipeline(
        &self,
        query_buffer: &Buffer, target_buffer: &Buffer,
        query_len_buffer: &Buffer, target_len_buffer: &Buffer,
        query_offset_buffer: &Buffer, target_offset_buffer: &Buffer,
        scoring_buffer: &Buffer, dp_buffer: &Buffer,
        matrix_offset_buffer: &Buffer, results_buffer: &Buffer,
        batch_size: usize
    ) -> Result<(), String> {
        let command_buffer = self.command_queue.new_command_buffer();

        // Phase 1: Anti-diagonal DP computation
        {
            let encoder = command_buffer.new_compute_command_encoder();
            encoder.set_compute_pipeline_state(&self.pipeline_antidiagonal);

            // Set buffers (same layout as serial version)
            encoder.set_buffer(0, Some(query_buffer), 0);
            encoder.set_buffer(1, Some(target_buffer), 0);
            encoder.set_buffer(2, Some(query_len_buffer), 0);
            encoder.set_buffer(3, Some(target_len_buffer), 0);
            encoder.set_buffer(4, Some(query_offset_buffer), 0);
            encoder.set_buffer(5, Some(target_offset_buffer), 0);
            encoder.set_buffer(6, Some(scoring_buffer), 0);
            encoder.set_buffer(7, Some(dp_buffer), 0);
            encoder.set_buffer(8, Some(matrix_offset_buffer), 0);

            // Configure threadgroups: one threadgroup per alignment, multiple threads per group
            let threads_per_group = 32; // Optimal for most Apple Silicon GPUs
            let threadgroup_size = MTLSize { width: threads_per_group, height: 1, depth: 1 };
            let num_threadgroups = MTLSize { width: batch_size as u64, height: 1, depth: 1 };

            encoder.dispatch_thread_groups(num_threadgroups, threadgroup_size);
            encoder.end_encoding();
        }

        // Phase 2: Find maximum scores
        {
            let encoder = command_buffer.new_compute_command_encoder();
            encoder.set_compute_pipeline_state(&self.pipeline_find_max);

            encoder.set_buffer(0, Some(dp_buffer), 0);
            encoder.set_buffer(1, Some(matrix_offset_buffer), 0);
            encoder.set_buffer(2, Some(query_len_buffer), 0);
            encoder.set_buffer(3, Some(target_len_buffer), 0);
            encoder.set_buffer(4, Some(results_buffer), 0);

            // One thread per alignment for max finding
            let thread_group_size = MTLSize {
                width: self.pipeline_find_max.thread_execution_width().min(batch_size as u64),
                height: 1,
                depth: 1,
            };
            let thread_groups = MTLSize {
                width: ((batch_size as u64 + thread_group_size.width - 1) / thread_group_size.width),
                height: 1,
                depth: 1,
            };

            encoder.dispatch_thread_groups(thread_groups, thread_group_size);
            encoder.end_encoding();
        }

        command_buffer.commit();
        command_buffer.wait_until_completed();

        Ok(())
    }

    /// Execute legacy combined pipeline (fallback for small alignments)
    #[allow(clippy::too_many_arguments)]
    fn execute_combined_pipeline(
        &self,
        query_buffer: &Buffer, target_buffer: &Buffer,
        query_len_buffer: &Buffer, target_len_buffer: &Buffer,
        query_offset_buffer: &Buffer, target_offset_buffer: &Buffer,
        scoring_buffer: &Buffer, dp_buffer: &Buffer,
        matrix_offset_buffer: &Buffer, results_buffer: &Buffer,
        batch_size: usize
    ) -> Result<(), String> {
        let command_buffer = self.command_queue.new_command_buffer();
        let encoder = command_buffer.new_compute_command_encoder();

        encoder.set_compute_pipeline_state(&self.pipeline_combined);
        encoder.set_buffer(0, Some(query_buffer), 0);
        encoder.set_buffer(1, Some(target_buffer), 0);
        encoder.set_buffer(2, Some(query_len_buffer), 0);
        encoder.set_buffer(3, Some(target_len_buffer), 0);
        encoder.set_buffer(4, Some(query_offset_buffer), 0);
        encoder.set_buffer(5, Some(target_offset_buffer), 0);
        encoder.set_buffer(6, Some(scoring_buffer), 0);
        encoder.set_buffer(7, Some(dp_buffer), 0);
        encoder.set_buffer(8, Some(matrix_offset_buffer), 0);
        encoder.set_buffer(9, Some(results_buffer), 0);

        // Dispatch one thread per alignment (legacy approach)
        let thread_group_size = MTLSize {
            width: self.pipeline_combined.thread_execution_width().min(batch_size as u64),
            height: 1,
            depth: 1,
        };
        let thread_groups = MTLSize {
            width: ((batch_size as u64 + thread_group_size.width - 1) / thread_group_size.width),
            height: 1,
            depth: 1,
        };

        encoder.dispatch_thread_groups(thread_groups, thread_group_size);
        encoder.end_encoding();

        command_buffer.commit();
        command_buffer.wait_until_completed();

        Ok(())
    }
}

/// Batch Smith-Waterman alignment on GPU
///
/// Convenience function that creates a `GpuAlignmentBatch` and runs alignment.
///
/// For repeated alignments, prefer creating a `GpuAlignmentBatch` instance
/// to avoid re-creating Metal resources.
///
/// # Example
///
/// ```no_run
/// use biometal::alignment::gpu::smith_waterman_batch_gpu;
/// use biometal::alignment::ScoringMatrix;
///
/// let queries = vec![b"ACGTACGT", b"TGCATGCA"];
/// let targets = vec![b"ACGGACGG", b"TGCATGCA"];
/// let query_refs: Vec<&[u8]> = queries.iter().map(|q| q.as_slice()).collect();
/// let target_refs: Vec<&[u8]> = targets.iter().map(|t| t.as_slice()).collect();
///
/// let scoring = ScoringMatrix::default();
/// let alignments = smith_waterman_batch_gpu(&query_refs, &target_refs, &scoring).unwrap();
/// ```
pub fn smith_waterman_batch_gpu(
    queries: &[&[u8]],
    targets: &[&[u8]],
    scoring: &ScoringMatrix,
) -> Result<Vec<Alignment>, String> {
    let batch = GpuAlignmentBatch::new()?;
    batch.align_batch(queries, targets, scoring)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gpu_availability() {
        // Test if Metal device is available on macOS
        #[cfg(target_os = "macos")]
        {
            let result = GpuAlignmentBatch::new();
            assert!(result.is_ok(), "Metal device should be available on macOS");
        }

        #[cfg(not(target_os = "macos"))]
        {
            let result = GpuAlignmentBatch::new();
            assert!(result.is_err(), "Metal device should not be available on non-macOS");
        }
    }

    #[test]
    #[cfg(target_os = "macos")]
    fn test_simple_alignment() {
        let batch = GpuAlignmentBatch::new().unwrap();

        let queries = vec![b"ACGTACGT".as_slice()];
        let targets = vec![b"ACGTACGT".as_slice()];
        let scoring = ScoringMatrix::default();

        let results = batch.align_batch(&queries, &targets, &scoring).unwrap();

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].score, 16); // 8 matches × 2 = 16
    }

    #[test]
    #[cfg(target_os = "macos")]
    fn test_batch_alignment() {
        let batch = GpuAlignmentBatch::new().unwrap();

        let queries = vec![
            b"ACGTACGT".as_slice(),
            b"TGCATGCA".as_slice(),
            b"GGCCGGCC".as_slice(),
        ];
        let targets = vec![
            b"ACGTACGT".as_slice(), // Perfect match
            b"TGCATGCA".as_slice(), // Perfect match
            b"GGCTGGCT".as_slice(), // 1 mismatch
        ];
        let scoring = ScoringMatrix::default();

        let results = batch.align_batch(&queries, &targets, &scoring).unwrap();

        assert_eq!(results.len(), 3);
        assert_eq!(results[0].score, 16); // 8 matches × 2
        assert_eq!(results[1].score, 16); // 8 matches × 2
        assert!(results[2].score < 16);    // 1 mismatch reduces score
    }
}
