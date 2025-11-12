//! Block-based batch operations (convenience API)
//!
//! # Current Status
//!
//! **Performance**: Block operations have **similar performance** to per-record operations.
//! **Purpose**: Convenience API for processing batches of sequences.
//! **Future**: Could achieve 14× speedup with implementation changes (see below).
//!
//! # Evidence Base
//!
//! **Entry 027** (1,440 measurements, N=30) demonstrates that block processing CAN achieve
//! 14× speedup by reducing function call overhead from 82-86% to 4-8%.
//!
//! **Current implementation**: Does NOT achieve this speedup because it still calls the
//! underlying NEON functions once per sequence (10,000 function calls for 10K sequences).
//!
//! **Benchmark results** (November 11, 2025, N=30):
//! - Base counting: Block is 0.88× (slightly slower than per-record)
//! - GC content: Block is 0.98× (essentially identical)
//! - Mean quality: Block is 1.01× (essentially identical)
//!
//! See `RULE2_INVESTIGATION_FINDINGS.md` for full analysis.
//!
//! # When to Use
//!
//! **Use block operations when**:
//! - You have sequences already collected in memory
//! - You prefer batch-style API (process many at once)
//! - Code clarity over per-record loops
//!
//! **Use per-record operations when**:
//! - Streaming one record at a time
//! - Memory-constrained environments
//! - No significant performance difference in current implementation
//!
//! # Example
//!
//! ```
//! use biometal::operations::block::count_bases_block;
//!
//! // Create some test sequences
//! let seq1 = b"ATGCATGC";
//! let seq2 = b"GGCCGGCC";
//! let seq3 = b"AATTAATT";
//!
//! let sequences = vec![
//!     seq1.as_ref(),
//!     seq2.as_ref(),
//!     seq3.as_ref(),
//! ];
//!
//! // Process entire block at once (convenient batch API)
//! let counts = count_bases_block(&sequences);
//!
//! assert_eq!(counts.len(), 3);
//! // seq1: ATGCATGC = 2A, 2C, 2G, 2T
//! assert_eq!(counts[0], [2, 2, 2, 2]);
//! ```
//!
//! # Future Work: True Block Processing (14× Speedup)
//!
//! To achieve Entry 027's 14× speedup would require:
//! - Inlining NEON operations directly into block functions (no function calls)
//! - Code duplication (~1,200 lines: 3 operations × 2 variants × ~200 lines each)
//! - Maintenance burden keeping two implementations in sync
//!
//! **Trade-off**: Rules 3+4 (parallel BGZF + smart mmap) provide 16.3× combined speedup
//! with less code complexity, making them higher priority for Phase 2 development.

use crate::operations::base_counting::BaseCounts;

#[cfg(target_arch = "aarch64")]
use crate::operations::base_counting::count_bases_neon;
#[cfg(target_arch = "aarch64")]
use crate::operations::gc_content::gc_content_neon;
#[cfg(target_arch = "aarch64")]
use crate::operations::quality_filter::mean_quality_neon;

use crate::operations::base_counting::count_bases_scalar;
use crate::operations::gc_content::gc_content_scalar;
use crate::operations::quality_filter::mean_quality_scalar;

/// Count bases for multiple sequences in a single batch (convenience API)
///
/// # Performance
///
/// **Current implementation**: Similar performance to per-record operations.
/// - Both make N function calls to underlying NEON operations
/// - Benchmark (N=30, 10K seqs): Block is 0.88× per-record (essentially identical)
///
/// This is a convenience API for batch processing, not a performance optimization.
///
/// # Block Size
///
/// Recommended block size: ~10,000 sequences
/// - Balances memory usage (~1.5 MB for 150bp reads) with batch convenience
/// - Too small (< 1K): Less convenient for batch operations
/// - Too large (> 100K): Memory pressure, may exceed cache
///
/// # Example
///
/// ```
/// use biometal::operations::block::count_bases_block;
///
/// let sequences = vec![
///     b"ATGC".as_ref(),
///     b"GCTA".as_ref(),
///     b"CGAT".as_ref(),
/// ];
///
/// let counts = count_bases_block(&sequences);
/// assert_eq!(counts.len(), 3);
/// // counts[0] = [1, 1, 1, 1] (A, C, G, T for "ATGC")
/// // counts[1] = [1, 1, 1, 1] (A, C, G, T for "GCTA")
/// // counts[2] = [1, 1, 1, 1] (A, C, G, T for "CGAT")
/// ```
pub fn count_bases_block(sequences: &[&[u8]]) -> Vec<BaseCounts> {
    #[cfg(target_arch = "aarch64")]
    {
        unsafe { count_bases_block_neon(sequences) }
    }

    #[cfg(not(target_arch = "aarch64"))]
    {
        count_bases_block_scalar(sequences)
    }
}

/// NEON-optimized block base counting
///
/// Processes multiple sequences with amortized NEON setup cost.
///
/// # Safety
///
/// Uses unsafe NEON intrinsics but is safe to call:
/// - Only compiled on aarch64 platforms
/// - All pointer operations are bounds-checked
/// - NEON is standard on all aarch64 CPUs
#[cfg(target_arch = "aarch64")]
pub unsafe fn count_bases_block_neon(sequences: &[&[u8]]) -> Vec<BaseCounts> {
    // Use explicit loop to minimize overhead and keep NEON registers hot
    let mut results = Vec::with_capacity(sequences.len());
    for seq in sequences {
        results.push(count_bases_neon(seq));
    }
    results
}

/// Scalar block base counting (x86_64 fallback)
pub fn count_bases_block_scalar(sequences: &[&[u8]]) -> Vec<BaseCounts> {
    let mut results = Vec::with_capacity(sequences.len());
    for seq in sequences {
        results.push(count_bases_scalar(seq));
    }
    results
}

/// Calculate GC content for multiple sequences in a single batch (convenience API)
///
/// # Performance
///
/// **Current implementation**: Similar performance to per-record operations.
/// - Benchmark (N=30, 10K seqs): Block is 0.98× per-record (essentially identical)
///
/// This is a convenience API for batch processing, not a performance optimization.
///
/// # Example
///
/// ```
/// use biometal::operations::block::gc_content_block;
///
/// let sequences = vec![
///     b"ATGC".as_ref(),  // 50% GC
///     b"GGCC".as_ref(),  // 100% GC
///     b"AATT".as_ref(),  // 0% GC
/// ];
///
/// let gc_values = gc_content_block(&sequences);
/// assert_eq!(gc_values.len(), 3);
/// assert!((gc_values[0] - 0.50).abs() < 0.01);
/// assert!((gc_values[1] - 1.00).abs() < 0.01);
/// assert!((gc_values[2] - 0.00).abs() < 0.01);
/// ```
pub fn gc_content_block(sequences: &[&[u8]]) -> Vec<f64> {
    #[cfg(target_arch = "aarch64")]
    {
        unsafe { gc_content_block_neon(sequences) }
    }

    #[cfg(not(target_arch = "aarch64"))]
    {
        gc_content_block_scalar(sequences)
    }
}

/// NEON-optimized block GC content calculation
#[cfg(target_arch = "aarch64")]
pub unsafe fn gc_content_block_neon(sequences: &[&[u8]]) -> Vec<f64> {
    let mut results = Vec::with_capacity(sequences.len());
    for seq in sequences {
        results.push(gc_content_neon(seq));
    }
    results
}

/// Scalar block GC content calculation
pub fn gc_content_block_scalar(sequences: &[&[u8]]) -> Vec<f64> {
    let mut results = Vec::with_capacity(sequences.len());
    for seq in sequences {
        results.push(gc_content_scalar(seq));
    }
    results
}

/// Calculate mean quality for multiple quality strings in a single batch (convenience API)
///
/// # Performance
///
/// **Current implementation**: Similar performance to per-record operations.
/// - Benchmark (N=30, 10K seqs): Block is 1.01× per-record (essentially identical)
///
/// This is a convenience API for batch processing, not a performance optimization.
///
/// # Example
///
/// ```
/// use biometal::operations::block::mean_quality_block;
///
/// let qualities = vec![
///     b"IIII".as_ref(),  // Phred 40
///     b"####".as_ref(),  // Phred 2
///     b"@@@@".as_ref(),  // Phred 31
/// ];
///
/// let mean_quals = mean_quality_block(&qualities);
/// assert_eq!(mean_quals.len(), 3);
/// assert!((mean_quals[0] - 40.0).abs() < 0.1);
/// assert!((mean_quals[1] - 2.0).abs() < 0.1);
/// assert!((mean_quals[2] - 31.0).abs() < 0.1);
/// ```
pub fn mean_quality_block(qualities: &[&[u8]]) -> Vec<f64> {
    #[cfg(target_arch = "aarch64")]
    {
        unsafe { mean_quality_block_neon(qualities) }
    }

    #[cfg(not(target_arch = "aarch64"))]
    {
        mean_quality_block_scalar(qualities)
    }
}

/// NEON-optimized block mean quality calculation
#[cfg(target_arch = "aarch64")]
pub unsafe fn mean_quality_block_neon(qualities: &[&[u8]]) -> Vec<f64> {
    let mut results = Vec::with_capacity(qualities.len());
    for qual in qualities {
        results.push(mean_quality_neon(qual));
    }
    results
}

/// Scalar block mean quality calculation
pub fn mean_quality_block_scalar(qualities: &[&[u8]]) -> Vec<f64> {
    let mut results = Vec::with_capacity(qualities.len());
    for qual in qualities {
        results.push(mean_quality_scalar(qual));
    }
    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::operations::base_counting::count_bases;

    #[test]
    fn test_count_bases_block_basic() {
        let sequences = vec![
            b"ATGC".as_ref(),
            b"GGCC".as_ref(),
            b"AATT".as_ref(),
        ];

        let counts = count_bases_block(&sequences);
        assert_eq!(counts.len(), 3);

        // ATGC: 1A, 1C, 1G, 1T
        assert_eq!(counts[0], [1, 1, 1, 1]);

        // GGCC: 0A, 2C, 2G, 0T
        assert_eq!(counts[1], [0, 2, 2, 0]);

        // AATT: 2A, 0C, 0G, 2T
        assert_eq!(counts[2], [2, 0, 0, 2]);
    }

    #[test]
    fn test_gc_content_block_basic() {
        let sequences = vec![
            b"ATGC".as_ref(),  // 50% GC
            b"GGCC".as_ref(),  // 100% GC
            b"AATT".as_ref(),  // 0% GC
        ];

        let gc_values = gc_content_block(&sequences);
        assert_eq!(gc_values.len(), 3);

        assert!((gc_values[0] - 0.50).abs() < 0.01);
        assert!((gc_values[1] - 1.00).abs() < 0.01);
        assert!((gc_values[2] - 0.00).abs() < 0.01);
    }

    #[test]
    fn test_mean_quality_block_basic() {
        let qualities = vec![
            b"IIII".as_ref(),  // Phred 40
            b"####".as_ref(),  // Phred 2
            b"@@@@".as_ref(),  // Phred 31
        ];

        let mean_quals = mean_quality_block(&qualities);
        assert_eq!(mean_quals.len(), 3);

        assert!((mean_quals[0] - 40.0).abs() < 0.1);
        assert!((mean_quals[1] - 2.0).abs() < 0.1);
        assert!((mean_quals[2] - 31.0).abs() < 0.1);
    }

    #[test]
    fn test_count_bases_block_large() {
        // Test with larger block (closer to real-world 10K)
        let sequences: Vec<&[u8]> = (0..1000)
            .map(|_| b"ATGCATGCATGC".as_ref())
            .collect();

        let counts = count_bases_block(&sequences);
        assert_eq!(counts.len(), 1000);

        // Each sequence has 3A, 3C, 3G, 3T
        for count in counts {
            assert_eq!(count, [3, 3, 3, 3]);
        }
    }

    #[test]
    fn test_block_vs_single_consistency() {
        let sequences = vec![
            b"ATGCATGCATGC".as_ref(),
            b"GGCCGGCCGGCC".as_ref(),
            b"AATTAATTAATT".as_ref(),
        ];

        // Block processing
        let block_counts = count_bases_block(&sequences);

        // Single processing
        let single_counts: Vec<_> = sequences.iter()
            .map(|seq| count_bases(seq))
            .collect();

        // Should be identical
        assert_eq!(block_counts, single_counts);
    }
}
