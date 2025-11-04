//! Quality filtering with ARM NEON SIMD optimization (Rule 1)
//!
//! # Evidence
//!
//! Entry 020-025 (Lab Notebook):
//! - **Speedup**: 25.1× faster than scalar
//! - **Statistical rigor**: Cohen's d = 5.87 (very large effect)
//! - **Cross-platform**: Mac M4 Max, AWS Graviton 3
//!
//! # Architecture
//!
//! Computes mean Phred quality scores using NEON SIMD for rapid filtering.

/// Calculate mean Phred quality score
///
/// Quality scores are Phred+33 encoded. This function:
/// 1. Converts ASCII to numeric quality (Q = ASCII - 33)
/// 2. Computes mean quality
///
/// # Platform-Specific Optimization
///
/// - **ARM (aarch64)**: Uses NEON SIMD (25.1× speedup)
/// - **x86_64**: Uses scalar fallback (portable)
///
/// # Example
///
/// ```
/// use biometal::operations::mean_quality;
///
/// let quality = b"IIIIIIIII"; // Phred+33, Q=40 for each
/// let mean_q = mean_quality(quality);
/// assert!((mean_q - 40.0).abs() < 0.1);
/// ```
pub fn mean_quality(quality: &[u8]) -> f64 {
    #[cfg(target_arch = "aarch64")]
    {
        unsafe { mean_quality_neon(quality) }
    }

    #[cfg(not(target_arch = "aarch64"))]
    {
        mean_quality_scalar(quality)
    }
}

/// NEON-optimized mean quality calculation (25.1× faster than scalar)
///
/// # Evidence
///
/// Entry 020-025, Cohen's d = 5.87 (very large effect)
///
/// # Safety
///
/// This function uses unsafe NEON intrinsics but is safe to call:
/// - Only called on aarch64 platforms (compile-time check)
/// - NEON is standard on all aarch64 CPUs
/// - Pointer operations are bounds-checked via chunks_exact
#[cfg(target_arch = "aarch64")]
pub unsafe fn mean_quality_neon(quality: &[u8]) -> f64 {
    use std::arch::aarch64::*;

    if quality.is_empty() {
        return 0.0;
    }

    let offset_vec = vdupq_n_u8(33); // Phred+33 offset
    let mut sum_vcount = vdupq_n_u32(0);

    let chunks = quality.chunks_exact(16);
    let remainder = chunks.remainder();

    for chunk in chunks {
        let qual_vec = vld1q_u8(chunk.as_ptr());

        // Subtract 33 to get quality scores
        let q_vec = vsubq_u8(qual_vec, offset_vec);

        // Widen to u32 and accumulate
        sum_vcount = vaddq_u32(sum_vcount, vpaddlq_u16(vpaddlq_u8(q_vec)));
    }

    // Extract sum
    let mut sum = 0u32;
    sum += vgetq_lane_u32(sum_vcount, 0);
    sum += vgetq_lane_u32(sum_vcount, 1);
    sum += vgetq_lane_u32(sum_vcount, 2);
    sum += vgetq_lane_u32(sum_vcount, 3);

    // Handle remainder
    for &q in remainder {
        sum += (q - 33) as u32;
    }

    sum as f64 / quality.len() as f64
}

/// Scalar fallback for non-ARM platforms
///
/// This provides a portable implementation for x86_64 and other architectures.
pub fn mean_quality_scalar(quality: &[u8]) -> f64 {
    if quality.is_empty() {
        return 0.0;
    }

    let sum: u32 = quality.iter().map(|&q| (q - 33) as u32).sum();
    sum as f64 / quality.len() as f64
}

/// Filter FASTQ record by minimum mean quality
pub fn passes_quality_filter(quality: &[u8], min_quality: f64) -> bool {
    mean_quality(quality) >= min_quality
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mean_quality_high() {
        let quality = b"IIIIII"; // Q=40
        let mean_q = mean_quality(quality);
        assert!((mean_q - 40.0).abs() < 0.1);
    }

    #[test]
    fn test_mean_quality_low() {
        let quality = b"!!!!!!"; // Q=0
        let mean_q = mean_quality(quality);
        assert!((mean_q - 0.0).abs() < 0.1);
    }

    #[test]
    fn test_mean_quality_mixed() {
        let quality = b"!I"; // Q=0 and Q=40, mean=20
        let mean_q = mean_quality(quality);
        assert!((mean_q - 20.0).abs() < 0.1);
    }

    #[test]
    fn test_passes_filter() {
        let high_quality = b"IIIIII"; // Q=40
        assert!(passes_quality_filter(high_quality, 30.0));
        assert!(!passes_quality_filter(high_quality, 50.0));
    }

    #[cfg(target_arch = "aarch64")]
    #[test]
    fn test_neon_matches_scalar() {
        let qualities = vec![
            b"IIIIII".as_slice(),
            b"!!!!!!".as_slice(),
            b"IIIIIIIIIIIIIIIIIII".as_slice(), // >16 bytes
        ];

        for qual in qualities {
            let neon = unsafe { mean_quality_neon(qual) };
            let scalar = mean_quality_scalar(qual);
            assert!((neon - scalar).abs() < 0.1);
        }
    }
}
