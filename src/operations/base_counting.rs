//! Base counting with ARM NEON SIMD optimization (Rule 1)
//!
//! # Evidence
//!
//! Entry 020 (Lab Notebook):
//! - **Speedup**: 16.7× faster than scalar
//! - **Statistical rigor**: Cohen's d = 4.82 (very large effect)
//! - **Cross-platform**: Mac M4 Max, AWS Graviton 3
//!
//! # Architecture
//!
//! This module provides both NEON (ARM) and scalar (portable) implementations:
//! - NEON: Processes 16 bytes at a time with SIMD instructions
//! - Scalar: Sequential processing (x86_64 fallback)
//!
//! The public API automatically selects the best implementation for the platform.

/// Base count result [A, C, G, T]
pub type BaseCounts = [u32; 4];

/// Count occurrences of each DNA base (A, C, G, T)
///
/// # Platform-Specific Optimization
///
/// - **ARM (aarch64)**: Uses NEON SIMD (16.7× speedup)
/// - **x86_64**: Uses scalar fallback (portable)
///
/// # Evidence
///
/// Entry 020: NEON base counting achieves 16.7× speedup with Cohen's d = 4.82
///
/// # Example
///
/// ```
/// use biometal::operations::base_counting::count_bases;
///
/// let sequence = b"GATTACAGATTACA";
/// let counts = count_bases(sequence);
/// // counts = [A_count, C_count, G_count, T_count]
/// ```
pub fn count_bases(seq: &[u8]) -> BaseCounts {
    #[cfg(target_arch = "aarch64")]
    {
        unsafe { count_bases_neon(seq) }
    }

    #[cfg(not(target_arch = "aarch64"))]
    {
        count_bases_scalar(seq)
    }
}

/// NEON-optimized base counting (16.7× faster than scalar)
///
/// # Evidence
///
/// Entry 020, Cohen's d = 4.82 (very large effect)
///
/// # Safety
///
/// This function uses unsafe NEON intrinsics but is safe to call:
/// - Only called on aarch64 platforms (compile-time check)
/// - NEON is standard on all aarch64 CPUs
/// - Pointer operations are bounds-checked via chunks_exact
#[cfg(target_arch = "aarch64")]
pub unsafe fn count_bases_neon(seq: &[u8]) -> BaseCounts {
    use std::arch::aarch64::*;

    let mut counts = [0u32; 4];

    // NEON registers for ACGT counts
    let mut vcounts = [vdupq_n_u32(0); 4];

    // Process 16 bytes at a time
    let chunks = seq.chunks_exact(16);
    let remainder = chunks.remainder();

    for chunk in chunks {
        let seq_vec = vld1q_u8(chunk.as_ptr());

        // Compare against A, C, G, T (case-insensitive: 8 NEON comparisons in parallel)
        // Result: 0xFF for match, 0x00 for non-match
        let a_upper_mask = vceqq_u8(seq_vec, vdupq_n_u8(b'A'));
        let a_lower_mask = vceqq_u8(seq_vec, vdupq_n_u8(b'a'));
        let a_mask = vorrq_u8(a_upper_mask, a_lower_mask);

        let c_upper_mask = vceqq_u8(seq_vec, vdupq_n_u8(b'C'));
        let c_lower_mask = vceqq_u8(seq_vec, vdupq_n_u8(b'c'));
        let c_mask = vorrq_u8(c_upper_mask, c_lower_mask);

        let g_upper_mask = vceqq_u8(seq_vec, vdupq_n_u8(b'G'));
        let g_lower_mask = vceqq_u8(seq_vec, vdupq_n_u8(b'g'));
        let g_mask = vorrq_u8(g_upper_mask, g_lower_mask);

        let t_upper_mask = vceqq_u8(seq_vec, vdupq_n_u8(b'T'));
        let t_lower_mask = vceqq_u8(seq_vec, vdupq_n_u8(b't'));
        let t_mask = vorrq_u8(t_upper_mask, t_lower_mask);

        // Convert 0xFF to 0x01 by shifting right 7 bits
        // This converts mask (0xFF = match, 0x00 = no match) to count (0x01 or 0x00)
        let a_count = vshrq_n_u8(a_mask, 7);
        let c_count = vshrq_n_u8(c_mask, 7);
        let g_count = vshrq_n_u8(g_mask, 7);
        let t_count = vshrq_n_u8(t_mask, 7);

        // Accumulate counts
        // vpaddlq_u8: pairwise add and widen u8 -> u16
        // vpaddlq_u16: pairwise add and widen u16 -> u32
        vcounts[0] = vaddq_u32(vcounts[0], vpaddlq_u16(vpaddlq_u8(a_count)));
        vcounts[1] = vaddq_u32(vcounts[1], vpaddlq_u16(vpaddlq_u8(c_count)));
        vcounts[2] = vaddq_u32(vcounts[2], vpaddlq_u16(vpaddlq_u8(g_count)));
        vcounts[3] = vaddq_u32(vcounts[3], vpaddlq_u16(vpaddlq_u8(t_count)));
    }

    // Extract counts from NEON registers (each register holds 4x u32)
    for i in 0..4 {
        counts[i] = vgetq_lane_u32(vcounts[i], 0)
            + vgetq_lane_u32(vcounts[i], 1)
            + vgetq_lane_u32(vcounts[i], 2)
            + vgetq_lane_u32(vcounts[i], 3);
    }

    // Handle remainder with scalar code (case-insensitive)
    for &base in remainder {
        match base {
            b'A' | b'a' => counts[0] += 1,
            b'C' | b'c' => counts[1] += 1,
            b'G' | b'g' => counts[2] += 1,
            b'T' | b't' => counts[3] += 1,
            _ => {} // Ignore non-ACGT characters
        }
    }

    counts
}

/// Scalar fallback for non-ARM platforms
///
/// This provides a portable implementation for x86_64 and other architectures.
pub fn count_bases_scalar(seq: &[u8]) -> BaseCounts {
    let mut counts = [0u32; 4];

    for &base in seq {
        match base {
            b'A' | b'a' => counts[0] += 1,
            b'C' | b'c' => counts[1] += 1,
            b'G' | b'g' => counts[2] += 1,
            b'T' | b't' => counts[3] += 1,
            _ => {} // Ignore non-ACGT characters
        }
    }

    counts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_count_bases_basic() {
        let seq = b"GATTACA";
        let counts = count_bases(seq);
        assert_eq!(counts, [3, 1, 1, 2]); // A=3, C=1, G=1, T=2
    }

    #[test]
    fn test_count_bases_all_same() {
        let seq = b"AAAAAAAAAA";
        let counts = count_bases(seq);
        assert_eq!(counts, [10, 0, 0, 0]);
    }

    #[test]
    fn test_count_bases_with_n() {
        let seq = b"GATNACAN"; // N should be ignored
        let counts = count_bases(seq);
        // G=1, A=3, T=1, N=2(ignored), C=1, A=1(counted), N=1(ignored)
        assert_eq!(counts, [3, 1, 1, 1]); // A=3, C=1, G=1, T=1
    }

    #[test]
    fn test_count_bases_mixed_case() {
        // Test case-insensitive counting (the bug fix)
        let seq = b"AAAAaaaaTTTTtttt"; // 4A + 4a + 4T + 4t = 8A, 8T
        let counts = count_bases(seq);
        assert_eq!(counts, [8, 0, 0, 8]); // A=8, C=0, G=0, T=8
    }

    #[test]
    fn test_count_bases_soft_masked_genome() {
        // Realistic soft-masked genome sequence
        let seq = b"ATCGNNNNatcgtatcgATCGNNNNatcgtatcg";
        let counts = count_bases(seq);
        // Expected: A=6, T=8, C=6, G=6 (case-insensitive, N ignored)
        assert_eq!(counts, [6, 6, 6, 8]); // A=6, C=6, G=6, T=8
    }

    #[test]
    fn test_count_bases_all_lowercase() {
        let seq = b"gattaca"; // All lowercase
        let counts = count_bases(seq);
        assert_eq!(counts, [3, 1, 1, 2]); // a=3, c=1, g=1, t=2
    }

    #[test]
    fn test_count_bases_empty() {
        let seq = b"";
        let counts = count_bases(seq);
        assert_eq!(counts, [0, 0, 0, 0]);
    }

    #[test]
    fn test_count_bases_large() {
        // Test with sequence larger than 16 bytes (tests NEON chunking)
        let seq = b"GATTACAGATTACAGATTACAGATTACA"; // 28 bytes
        let counts = count_bases(seq);
        // G=4, A=12, T=8, C=4
        assert_eq!(counts, [12, 4, 4, 8]);
    }

    #[cfg(target_arch = "aarch64")]
    #[test]
    fn test_neon_matches_scalar() {
        // Verify NEON and scalar implementations produce identical results
        let sequences = vec![
            b"GATTACA".as_slice(),
            b"AAAA".as_slice(),
            b"CCCC".as_slice(),
            b"GGGG".as_slice(),
            b"TTTT".as_slice(),
            b"ACGTACGTACGTACGT".as_slice(), // Exactly 16 bytes
            b"ACGTACGTACGTACGTACGT".as_slice(), // >16 bytes
            b"NNNACGTNNNN".as_slice(),       // With non-ACGT
        ];

        for seq in sequences {
            let neon_result = unsafe { count_bases_neon(seq) };
            let scalar_result = count_bases_scalar(seq);
            assert_eq!(
                neon_result, scalar_result,
                "NEON and scalar results differ for sequence: {:?}",
                std::str::from_utf8(seq).unwrap()
            );
        }
    }

    #[cfg(target_arch = "aarch64")]
    mod proptests {
        use super::*;
        use proptest::prelude::*;

        proptest! {
            #[test]
            fn test_neon_matches_scalar_proptest(seq in "[ACGTNacgtn]{0,1000}") {
                // Property: NEON and scalar implementations must produce identical results
                // NOW INCLUDES LOWERCASE for case-insensitive testing
                let neon_result = unsafe { count_bases_neon(seq.as_bytes()) };
                let scalar_result = count_bases_scalar(seq.as_bytes());
                prop_assert_eq!(
                    neon_result,
                    scalar_result,
                    "NEON and scalar results differ for sequence: {}",
                    seq
                );
            }

            #[test]
            fn test_count_bases_sum_equals_length(seq in "[ACGTacgt]{0,1000}") {
                // Property: Total count should equal sequence length (case-insensitive, excluding N)
                // NOW INCLUDES LOWERCASE for case-insensitive testing
                let counts = count_bases(seq.as_bytes());
                let total: u32 = counts.iter().sum();
                prop_assert_eq!(total, seq.len() as u32, "Total count should equal sequence length");
            }

            #[test]
            fn test_count_bases_monotonicity(seq in "[ACGTacgt]{0,100}") {
                // Property: Adding more bases should increase counts
                // NOW INCLUDES LOWERCASE for case-insensitive testing
                let counts1 = count_bases(seq.as_bytes());
                let extended = format!("{}AAAA", seq);
                let counts2 = count_bases(extended.as_bytes());

                // A count should increase by 4
                prop_assert_eq!(counts2[0], counts1[0] + 4, "A count should increase");
                // Other counts unchanged
                prop_assert_eq!(counts2[1], counts1[1], "C count should be unchanged");
                prop_assert_eq!(counts2[2], counts1[2], "G count should be unchanged");
                prop_assert_eq!(counts2[3], counts1[3], "T count should be unchanged");
            }
        }
    }
}
