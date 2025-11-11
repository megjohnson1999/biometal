//! Compression strategies (zstd, lz4, RLE).
//!
//! This module implements columnar compression optimized for different data types:
//! - **Raw**: No compression (for incompressible data like quality scores)
//! - **Zstd**: High compression for integers/positions (level 3, ~10× ratio)
//! - **Lz4**: Fast compression for sequences (>1 GB/s decompression)
//! - **RLE**: Run-length encoding for repetitive values (e.g., MAPQ)
//!
//! # Evidence Base
//!
//! - Zstd level 3: Optimal compression/speed trade-off (SPECIFICATION.md)
//! - Lz4 fast: >1 GB/s decompression for sequences (RFC 8878, SPECIFICATION.md)
//! - Raw for qualities: High entropy, incompressible (SPECIFICATION.md)
//!
//! # Compression Parameters
//!
//! | Column | Compression | Level | Rationale |
//! |--------|-------------|-------|-----------|
//! | Positions/Flags | zstd | 3 | Sorted integers, delta encoding |
//! | Sequences | lz4 | fast | >1 GB/s decompression |
//! | Qualities | Raw | none | High entropy, incompressible |
//! | MAPQ | Adaptive RLE | - | Often repetitive values |

use crate::{types::CompressionType, CafError, Result};

/// Default zstd compression level (from SPECIFICATION.md).
pub const ZSTD_LEVEL: i32 = 3;

/// Default dictionary size for quality score compression (110 KB).
/// Based on zstd documentation: 100-110 KB is optimal for most use cases.
pub const DICTIONARY_SIZE: usize = 110_000;

/// Maximum decompressed size for RLE (100 MB) to prevent decompression bombs.
const MAX_RLE_DECOMPRESSED_SIZE: usize = 100_000_000;

/// Compress data using the specified compression type.
///
/// # Example
///
/// ```
/// use caf::compression::compress;
/// use caf::types::CompressionType;
///
/// let data = b"ACGTACGTACGTACGT".repeat(100);
/// let compressed = compress(&data, CompressionType::Zstd).unwrap();
/// assert!(compressed.len() < data.len());
/// ```
pub fn compress(data: &[u8], compression_type: CompressionType) -> Result<Vec<u8>> {
    match compression_type {
        CompressionType::Raw => compress_raw(data),
        CompressionType::Zstd => compress_zstd(data, ZSTD_LEVEL),
        CompressionType::Lz4 => compress_lz4(data),
        CompressionType::Rle => compress_rle(data),
    }
}

/// Decompress data using the specified compression type.
///
/// # Example
///
/// ```
/// use caf::compression::{compress, decompress};
/// use caf::types::CompressionType;
///
/// let data = b"ACGTACGTACGTACGT".repeat(100);
/// let compressed = compress(&data, CompressionType::Zstd).unwrap();
/// let decompressed = decompress(&compressed, CompressionType::Zstd, data.len()).unwrap();
/// assert_eq!(decompressed, data);
/// ```
pub fn decompress(
    data: &[u8],
    compression_type: CompressionType,
    uncompressed_size: usize,
) -> Result<Vec<u8>> {
    match compression_type {
        CompressionType::Raw => decompress_raw(data),
        CompressionType::Zstd => decompress_zstd(data, uncompressed_size),
        CompressionType::Lz4 => decompress_lz4(data, uncompressed_size),
        CompressionType::Rle => decompress_rle(data),
    }
}

/// Raw compression (no-op, just copy data).
///
/// Used for incompressible data like quality scores.
fn compress_raw(data: &[u8]) -> Result<Vec<u8>> {
    Ok(data.to_vec())
}

/// Raw decompression (no-op, just copy data).
fn decompress_raw(data: &[u8]) -> Result<Vec<u8>> {
    Ok(data.to_vec())
}

/// Compress data using zstd at specified level.
///
/// Default level 3 provides good compression ratio (~10×) with fast decompression.
fn compress_zstd(data: &[u8], level: i32) -> Result<Vec<u8>> {
    zstd::encode_all(data, level).map_err(|e| CafError::CompressionError {
        column: "unknown".to_string(),
        block_id: 0,
        source: Box::new(e),
    })
}

/// Decompress zstd-compressed data.
fn decompress_zstd(data: &[u8], _uncompressed_size: usize) -> Result<Vec<u8>> {
    zstd::decode_all(data).map_err(|e| CafError::DecompressionError {
        column: "unknown".to_string(),
        block_id: 0,
        source: Box::new(e),
    })
}

/// Compress data using lz4 (fast mode).
///
/// Optimized for sequences: >1 GB/s decompression speed.
fn compress_lz4(data: &[u8]) -> Result<Vec<u8>> {
    lz4::block::compress(data, None, false).map_err(|e| CafError::CompressionError {
        column: "unknown".to_string(),
        block_id: 0,
        source: Box::new(std::io::Error::new(std::io::ErrorKind::Other, e)),
    })
}

/// Decompress lz4-compressed data.
fn decompress_lz4(data: &[u8], uncompressed_size: usize) -> Result<Vec<u8>> {
    lz4::block::decompress(data, Some(uncompressed_size as i32)).map_err(|e| {
        CafError::DecompressionError {
            column: "unknown".to_string(),
            block_id: 0,
            source: Box::new(std::io::Error::new(std::io::ErrorKind::Other, e)),
        }
    })
}

/// Compress data using run-length encoding.
///
/// Format: [value, count, value, count, ...]
/// Efficient for repetitive data like MAPQ scores.
///
/// # Example
///
/// Input:  [60, 60, 60, 60, 30, 30, 60]
/// Output: [60, 4, 30, 2, 60, 1]
fn compress_rle(data: &[u8]) -> Result<Vec<u8>> {
    if data.is_empty() {
        return Ok(Vec::new());
    }

    let mut compressed = Vec::with_capacity(data.len() / 2); // Estimate
    let mut current = data[0];
    let mut count: u32 = 1;

    for &byte in &data[1..] {
        if byte == current && count < u32::MAX {
            count += 1;
        } else {
            // Emit run
            compressed.push(current);
            // Encode count as varint (simplified: use 4 bytes for now)
            compressed.extend_from_slice(&count.to_le_bytes());
            current = byte;
            count = 1;
        }
    }

    // Emit final run
    compressed.push(current);
    compressed.extend_from_slice(&count.to_le_bytes());

    Ok(compressed)
}

/// Decompress RLE-encoded data.
fn decompress_rle(data: &[u8]) -> Result<Vec<u8>> {
    if data.is_empty() {
        return Ok(Vec::new());
    }

    if data.len() % 5 != 0 {
        return Err(CafError::DecompressionError {
            column: "unknown".to_string(),
            block_id: 0,
            source: Box::new(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Invalid RLE data length",
            )),
        });
    }

    let mut decompressed = Vec::new();

    for chunk in data.chunks_exact(5) {
        let value = chunk[0];
        let count = u32::from_le_bytes([chunk[1], chunk[2], chunk[3], chunk[4]]);

        // Check for decompression bomb
        if decompressed.len() + count as usize > MAX_RLE_DECOMPRESSED_SIZE {
            return Err(CafError::DecompressionError {
                column: "unknown".to_string(),
                block_id: 0,
                source: Box::new(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "RLE decompression would exceed {} bytes (decompression bomb?)",
                        MAX_RLE_DECOMPRESSED_SIZE
                    ),
                )),
            });
        }

        for _ in 0..count {
            decompressed.push(value);
        }
    }

    Ok(decompressed)
}

/// Train a Zstandard dictionary from sample data.
///
/// The dictionary learns common patterns in the data to improve compression.
/// Recommended for quality scores, read names, and other structured data.
///
/// # Arguments
///
/// * `samples` - Multiple samples of data to train on (e.g., quality scores from multiple blocks)
/// * `dict_size` - Target dictionary size (recommended: 100-110 KB)
///
/// # Example
///
/// ```no_run
/// use caf::compression::{train_dictionary, DICTIONARY_SIZE};
///
/// let quality_samples = vec![
///     vec![30, 30, 35, 37, 40, /* ... */],
///     vec![25, 28, 32, 35, 38, /* ... */],
///     // More samples...
/// ];
/// let dict = train_dictionary(&quality_samples, DICTIONARY_SIZE).unwrap();
/// ```
pub fn train_dictionary(samples: &[Vec<u8>], dict_size: usize) -> Result<Vec<u8>> {
    if samples.is_empty() {
        return Err(CafError::CompressionError {
            column: "dictionary".to_string(),
            block_id: 0,
            source: Box::new(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "Cannot train dictionary on empty samples",
            )),
        });
    }

    // Flatten samples into continuous buffer and track sizes
    let sample_sizes: Vec<usize> = samples.iter().map(|s| s.len()).collect();
    let total_size: usize = sample_sizes.iter().sum();
    let mut continuous = Vec::with_capacity(total_size);
    for sample in samples {
        continuous.extend_from_slice(sample);
    }

    // Train dictionary using from_continuous with sample sizes
    zstd::dict::from_continuous(&continuous, &sample_sizes, dict_size).map_err(|e| {
        CafError::CompressionError {
            column: "dictionary".to_string(),
            block_id: 0,
            source: Box::new(e),
        }
    })
}

/// Compress data using zstd with a trained dictionary.
///
/// Dictionary compression finds common patterns learned during training,
/// improving compression ratio by 1.5-3× for structured data.
pub fn compress_zstd_dict(data: &[u8], level: i32, dictionary: &[u8]) -> Result<Vec<u8>> {
    // Create encoder with dictionary (level is set during creation)
    let mut encoder = zstd::Encoder::with_dictionary(Vec::new(), level, dictionary)
        .map_err(|e| CafError::CompressionError {
            column: "unknown".to_string(),
            block_id: 0,
            source: Box::new(e),
        })?;

    std::io::copy(&mut std::io::Cursor::new(data), &mut encoder)
        .map_err(|e| CafError::CompressionError {
            column: "unknown".to_string(),
            block_id: 0,
            source: Box::new(e),
        })?;

    encoder.finish().map_err(|e| CafError::CompressionError {
        column: "unknown".to_string(),
        block_id: 0,
        source: Box::new(e),
    })
}

/// Decompress zstd-compressed data using a trained dictionary.
pub fn decompress_zstd_dict(data: &[u8], dictionary: &[u8]) -> Result<Vec<u8>> {
    let mut decoder = zstd::Decoder::with_dictionary(data, dictionary).map_err(|e| {
        CafError::DecompressionError {
            column: "unknown".to_string(),
            block_id: 0,
            source: Box::new(e),
        }
    })?;

    let mut decompressed = Vec::new();
    std::io::copy(&mut decoder, &mut decompressed).map_err(|e| {
        CafError::DecompressionError {
            column: "unknown".to_string(),
            block_id: 0,
            source: Box::new(e),
        }
    })?;

    Ok(decompressed)
}

/// Estimate compression ratio for given data and compression type.
///
/// This is used for adaptive compression selection.
pub fn estimate_compression_ratio(data: &[u8], compression_type: CompressionType) -> f64 {
    if data.is_empty() {
        return 1.0;
    }

    match compress(data, compression_type) {
        Ok(compressed) => data.len() as f64 / compressed.len() as f64,
        Err(_) => 1.0, // No compression on error
    }
}

/// Select best compression type for data.
///
/// Tries multiple compression types and returns the one with best ratio,
/// unless the improvement is marginal (< 10% better than raw).
pub fn select_compression(data: &[u8]) -> CompressionType {
    if data.is_empty() {
        return CompressionType::Raw;
    }

    let mut best_type = CompressionType::Raw;
    let mut best_ratio = 1.0;

    // Try zstd
    if let Ok(compressed) = compress_zstd(data, ZSTD_LEVEL) {
        let ratio = data.len() as f64 / compressed.len() as f64;
        if ratio > best_ratio {
            best_ratio = ratio;
            best_type = CompressionType::Zstd;
        }
    }

    // Try lz4
    if let Ok(compressed) = compress_lz4(data) {
        let ratio = data.len() as f64 / compressed.len() as f64;
        if ratio > best_ratio {
            best_ratio = ratio;
            best_type = CompressionType::Lz4;
        }
    }

    // Try RLE
    if let Ok(compressed) = compress_rle(data) {
        let ratio = data.len() as f64 / compressed.len() as f64;
        if ratio > best_ratio {
            best_ratio = ratio;
            best_type = CompressionType::Rle;
        }
    }

    // Only use compression if ratio > 1.1 (10% improvement)
    if best_ratio > 1.1 {
        best_type
    } else {
        CompressionType::Raw
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn test_raw_roundtrip() {
        let data = b"Hello, World!";
        let compressed = compress_raw(data).unwrap();
        assert_eq!(compressed, data);

        let decompressed = decompress_raw(&compressed).unwrap();
        assert_eq!(decompressed, data);
    }

    #[test]
    fn test_zstd_roundtrip() {
        let data = b"ACGTACGTACGTACGT".repeat(100);
        let compressed = compress_zstd(&data, ZSTD_LEVEL).unwrap();
        assert!(compressed.len() < data.len(), "Should compress");

        let decompressed = decompress_zstd(&compressed, data.len()).unwrap();
        assert_eq!(decompressed, data);
    }

    #[test]
    fn test_lz4_roundtrip() {
        let data = b"ACGTACGTACGTACGT".repeat(100);
        let compressed = compress_lz4(&data).unwrap();
        assert!(compressed.len() < data.len(), "Should compress");

        let decompressed = decompress_lz4(&compressed, data.len()).unwrap();
        assert_eq!(decompressed, data);
    }

    #[test]
    fn test_rle_roundtrip() {
        let data = vec![60, 60, 60, 60, 30, 30, 60];
        let compressed = compress_rle(&data).unwrap();
        // Compressed should be: [60, 4, 30, 2, 60, 1] = [60, 4,0,0,0, 30, 2,0,0,0, 60, 1,0,0,0]
        assert_eq!(compressed.len(), 15); // 3 runs × 5 bytes each

        let decompressed = decompress_rle(&compressed).unwrap();
        assert_eq!(decompressed, data);
    }

    #[test]
    fn test_rle_uniform() {
        let data = vec![42; 1000];
        let compressed = compress_rle(&data).unwrap();
        assert_eq!(compressed.len(), 5); // 1 value + 4 bytes count

        let decompressed = decompress_rle(&compressed).unwrap();
        assert_eq!(decompressed, data);
    }

    #[test]
    fn test_rle_empty() {
        let data: Vec<u8> = vec![];
        let compressed = compress_rle(&data).unwrap();
        assert!(compressed.is_empty());

        let decompressed = decompress_rle(&compressed).unwrap();
        assert!(decompressed.is_empty());
    }

    #[test]
    fn test_compress_decompress_api() {
        let data = b"ACGTACGTACGTACGT".repeat(100);

        // Test zstd
        let compressed = compress(&data, CompressionType::Zstd).unwrap();
        let decompressed = decompress(&compressed, CompressionType::Zstd, data.len()).unwrap();
        assert_eq!(decompressed, data);

        // Test lz4
        let compressed = compress(&data, CompressionType::Lz4).unwrap();
        let decompressed = decompress(&compressed, CompressionType::Lz4, data.len()).unwrap();
        assert_eq!(decompressed, data);

        // Test raw
        let compressed = compress(&data, CompressionType::Raw).unwrap();
        let decompressed = decompress(&compressed, CompressionType::Raw, data.len()).unwrap();
        assert_eq!(decompressed, data);
    }

    #[test]
    fn test_compression_ratios() {
        // Highly compressible data
        let data = b"AAAAAAAAAA".repeat(1000);
        let ratio = estimate_compression_ratio(&data, CompressionType::Zstd);
        assert!(ratio > 10.0, "Should compress very well");

        // Sequential data (has patterns)
        let sequential: Vec<u8> = (0..1000).map(|i| (i % 256) as u8).collect();
        let ratio = estimate_compression_ratio(&sequential, CompressionType::Zstd);
        // Sequential data compresses moderately well due to patterns
        assert!(ratio >= 1.0, "Ratio should be valid");
    }

    #[test]
    fn test_select_compression() {
        // Highly repetitive: should select zstd or rle
        let data = vec![60; 1000];
        let selected = select_compression(&data);
        assert!(
            selected == CompressionType::Zstd || selected == CompressionType::Rle,
            "Should select compression for repetitive data"
        );

        // Sequential data: may or may not compress well
        let sequential: Vec<u8> = (0..100).map(|i| i as u8).collect();
        let _selected = select_compression(&sequential);
        // May select raw or compression depending on data characteristics
    }

    #[test]
    fn test_empty_data() {
        let data: Vec<u8> = vec![];

        for compression_type in [
            CompressionType::Raw,
            CompressionType::Zstd,
            CompressionType::Lz4,
            CompressionType::Rle,
        ] {
            let compressed = compress(&data, compression_type).unwrap();
            let decompressed = decompress(&compressed, compression_type, 0).unwrap();
            assert_eq!(decompressed, data);
        }
    }

    #[test]
    fn test_large_data() {
        // Test with larger data (10 MB)
        let data = b"ACGTACGTACGTACGT".repeat(625_000); // ~10 MB
        let compressed = compress(&data, CompressionType::Zstd).unwrap();
        let ratio = data.len() as f64 / compressed.len() as f64;
        assert!(ratio > 5.0, "Should achieve good compression on repetitive data");

        let decompressed = decompress(&compressed, CompressionType::Zstd, data.len()).unwrap();
        assert_eq!(decompressed.len(), data.len());
        assert_eq!(&decompressed[..100], &data[..100]);
    }

    // Property-based tests
    proptest! {
        #[test]
        fn prop_raw_roundtrip(data in prop::collection::vec(any::<u8>(), 0..10000)) {
            let compressed = compress_raw(&data)?;
            let decompressed = decompress_raw(&compressed)?;
            prop_assert_eq!(decompressed, data);
        }

        #[test]
        fn prop_zstd_roundtrip(data in prop::collection::vec(any::<u8>(), 0..10000)) {
            let compressed = compress_zstd(&data, ZSTD_LEVEL)?;
            let decompressed = decompress_zstd(&compressed, data.len())?;
            prop_assert_eq!(decompressed, data);
        }

        #[test]
        fn prop_lz4_roundtrip(data in prop::collection::vec(any::<u8>(), 0..10000)) {
            let compressed = compress_lz4(&data)?;
            let decompressed = decompress_lz4(&compressed, data.len())?;
            prop_assert_eq!(decompressed, data);
        }

        #[test]
        fn prop_rle_roundtrip(data in prop::collection::vec(any::<u8>(), 0..1000)) {
            let compressed = compress_rle(&data)?;
            let decompressed = decompress_rle(&compressed)?;
            prop_assert_eq!(decompressed, data);
        }
    }
}
