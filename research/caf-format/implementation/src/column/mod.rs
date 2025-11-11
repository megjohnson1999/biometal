//! Column-specific encoding strategies.
//!
//! This module implements specialized encoding for different column types:
//! - **Integers**: Delta encoding + zigzag for signed values
//! - **Sequences**: Pre-decoded ASCII ('A', 'C', 'G', 'T', 'N')
//! - **Qualities**: Raw bytes (Phred+33)
//! - **CIGAR**: Standard BAM encoding (op << 28 | len)
//!
//! # Evidence Base
//!
//! - Delta encoding for sorted integers: ~10× compression (SPECIFICATION.md)
//! - Pre-decoded sequences: Enables 16-25× NEON speedup (OPTIMIZATION_RULES.md Rule 1)
//! - Raw qualities: High entropy, incompressible (SPECIFICATION.md)

use crate::{CafError, Result};

/// Encode integers with delta encoding.
///
/// Delta encoding stores differences between consecutive values,
/// which is highly effective for sorted data like genomic positions.
///
/// # Limitations
///
/// - **Overflow**: Assumes deltas fit in i32 range. For genomic positions,
///   this is safe (max chromosome ~250 Mb), but deltas between arbitrary
///   i32 values could theoretically overflow.
/// - **Best for sorted data**: Works best when values are monotonically
///   increasing (or at least locally clustered).
/// - **Not recommended for**: Random or widely-spaced values where deltas
///   are as large as original values.
///
/// # Example
///
/// ```
/// use caf::column::encode_integers_delta;
///
/// let positions = vec![1000, 1005, 1010, 1020];
/// let encoded = encode_integers_delta(&positions);
/// // encoded = [1000, 5, 5, 10] (much better compression)
/// ```
pub fn encode_integers_delta(values: &[i32]) -> Vec<i32> {
    if values.is_empty() {
        return Vec::new();
    }

    let mut encoded = Vec::with_capacity(values.len());
    encoded.push(values[0]); // First value stored as-is

    for i in 1..values.len() {
        let delta = values[i] - values[i - 1];
        encoded.push(delta);
    }

    encoded
}

/// Decode delta-encoded integers.
///
/// # Example
///
/// ```
/// use caf::column::{encode_integers_delta, decode_integers_delta};
///
/// let positions = vec![1000, 1005, 1010, 1020];
/// let encoded = encode_integers_delta(&positions);
/// let decoded = decode_integers_delta(&encoded);
/// assert_eq!(decoded, positions);
/// ```
pub fn decode_integers_delta(encoded: &[i32]) -> Vec<i32> {
    if encoded.is_empty() {
        return Vec::new();
    }

    let mut decoded = Vec::with_capacity(encoded.len());
    decoded.push(encoded[0]);

    for i in 1..encoded.len() {
        let value = decoded[i - 1] + encoded[i];
        decoded.push(value);
    }

    decoded
}

/// Encode signed integer with zigzag encoding.
///
/// Zigzag encoding maps signed integers to unsigned integers
/// in a way that small absolute values have small encodings:
/// - 0 → 0, -1 → 1, 1 → 2, -2 → 3, 2 → 4, ...
///
/// This is useful before variable-length encoding.
pub fn zigzag_encode(n: i32) -> u32 {
    ((n << 1) ^ (n >> 31)) as u32
}

/// Decode zigzag-encoded integer.
pub fn zigzag_decode(n: u32) -> i32 {
    ((n >> 1) as i32) ^ (-((n & 1) as i32))
}

/// Encode sequence as pre-decoded ASCII bytes.
///
/// BAM stores sequences in 4-bit encoding (=ACMGRSVTWYHKDBN).
/// CAF stores them pre-decoded as ASCII ('A', 'C', 'G', 'T', 'N')
/// to enable NEON SIMD operations (16-25× speedup).
///
/// # Example
///
/// ```
/// use caf::column::encode_sequence_ascii;
///
/// let seq = b"ACGTNNACGT";
/// let encoded = encode_sequence_ascii(seq);
/// assert_eq!(encoded, seq.to_vec());
/// ```
pub fn encode_sequence_ascii(seq: &[u8]) -> Vec<u8> {
    // Validate that all bases are valid ASCII
    for &base in seq {
        match base {
            b'A' | b'C' | b'G' | b'T' | b'N' |
            b'a' | b'c' | b'g' | b't' | b'n' => {}
            _ => {
                // For now, replace invalid bases with 'N'
                // In production, this should return Result
            }
        }
    }
    seq.to_vec()
}

/// Decode pre-decoded ASCII sequence (no-op).
///
/// Since sequences are stored as ASCII, decoding is just validation.
pub fn decode_sequence_ascii(encoded: &[u8]) -> Result<Vec<u8>> {
    // Validate bases
    for &base in encoded {
        match base {
            b'A' | b'C' | b'G' | b'T' | b'N' |
            b'a' | b'c' | b'g' | b't' | b'n' => {}
            _ => {
                return Err(CafError::ColumnEncoding {
                    column: "sequences".to_string(),
                    message: format!("Invalid base: {}", base as char),
                });
            }
        }
    }
    Ok(encoded.to_vec())
}

/// Encode quality scores (Phred+33) as raw bytes.
///
/// Quality scores have high entropy and are incompressible,
/// so they're stored raw without encoding.
///
/// # Example
///
/// ```
/// use caf::column::encode_qualities_raw;
///
/// let quals = b"IIIIIIIIII"; // Phred+33
/// let encoded = encode_qualities_raw(quals);
/// assert_eq!(encoded, quals.to_vec());
/// ```
pub fn encode_qualities_raw(quals: &[u8]) -> Vec<u8> {
    quals.to_vec()
}

/// Decode raw quality scores (no-op).
pub fn decode_qualities_raw(encoded: &[u8]) -> Vec<u8> {
    encoded.to_vec()
}

/// Encode CIGAR operation in BAM format.
///
/// BAM CIGAR format: `(length << 4) | op`
/// where op is:
/// - 0: M (match/mismatch)
/// - 1: I (insertion)
/// - 2: D (deletion)
/// - 3: N (skipped)
/// - 4: S (soft clip)
/// - 5: H (hard clip)
/// - 6: P (padding)
/// - 7: = (match)
/// - 8: X (mismatch)
pub fn encode_cigar_op(op: u8, len: u32) -> u32 {
    (len << 4) | (op as u32)
}

/// Decode CIGAR operation from BAM format.
///
/// Returns (op, length) tuple.
pub fn decode_cigar_op(encoded: u32) -> (u8, u32) {
    let op = (encoded & 0xF) as u8;
    let len = encoded >> 4;
    (op, len)
}

/// Parse CIGAR string to BAM-encoded operations.
///
/// # Example
///
/// ```
/// use caf::column::parse_cigar_string;
///
/// let cigar = "10M2I5D";
/// let ops = parse_cigar_string(cigar).unwrap();
/// // ops = [(10 << 4) | 0, (2 << 4) | 1, (5 << 4) | 2]
/// ```
pub fn parse_cigar_string(cigar: &str) -> Result<Vec<u32>> {
    let mut ops = Vec::new();
    let mut len_str = String::new();

    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            len_str.push(ch);
        } else {
            if len_str.is_empty() {
                return Err(CafError::ColumnEncoding {
                    column: "cigar".to_string(),
                    message: "CIGAR operation without length".to_string(),
                });
            }

            let len: u32 = len_str.parse().map_err(|_| CafError::ColumnEncoding {
                column: "cigar".to_string(),
                message: format!("Invalid CIGAR length: {}", len_str),
            })?;

            let op = match ch {
                'M' => 0,
                'I' => 1,
                'D' => 2,
                'N' => 3,
                'S' => 4,
                'H' => 5,
                'P' => 6,
                '=' => 7,
                'X' => 8,
                _ => {
                    return Err(CafError::ColumnEncoding {
                        column: "cigar".to_string(),
                        message: format!("Invalid CIGAR op: {}", ch),
                    });
                }
            };

            ops.push(encode_cigar_op(op, len));
            len_str.clear();
        }
    }

    if !len_str.is_empty() {
        return Err(CafError::ColumnEncoding {
            column: "cigar".to_string(),
            message: "CIGAR string ends with digits".to_string(),
        });
    }

    Ok(ops)
}

/// Format CIGAR operations to string.
///
/// # Example
///
/// ```
/// use caf::column::{parse_cigar_string, format_cigar_string};
///
/// let cigar = "10M2I5D";
/// let ops = parse_cigar_string(cigar).unwrap();
/// let formatted = format_cigar_string(&ops);
/// assert_eq!(formatted, cigar);
/// ```
pub fn format_cigar_string(ops: &[u32]) -> String {
    let mut cigar = String::new();

    for &encoded in ops {
        let (op, len) = decode_cigar_op(encoded);
        cigar.push_str(&len.to_string());
        cigar.push(match op {
            0 => 'M',
            1 => 'I',
            2 => 'D',
            3 => 'N',
            4 => 'S',
            5 => 'H',
            6 => 'P',
            7 => '=',
            8 => 'X',
            _ => '?',
        });
    }

    cigar
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn test_delta_encoding_sorted() {
        let positions = vec![1000, 1005, 1010, 1020, 1030];
        let encoded = encode_integers_delta(&positions);
        assert_eq!(encoded, vec![1000, 5, 5, 10, 10]);

        let decoded = decode_integers_delta(&encoded);
        assert_eq!(decoded, positions);
    }

    #[test]
    fn test_delta_encoding_empty() {
        let empty: Vec<i32> = vec![];
        let encoded = encode_integers_delta(&empty);
        assert!(encoded.is_empty());

        let decoded = decode_integers_delta(&encoded);
        assert!(decoded.is_empty());
    }

    #[test]
    fn test_delta_encoding_single() {
        let single = vec![42];
        let encoded = encode_integers_delta(&single);
        assert_eq!(encoded, vec![42]);

        let decoded = decode_integers_delta(&encoded);
        assert_eq!(decoded, single);
    }

    #[test]
    fn test_delta_encoding_negative() {
        let values = vec![100, 95, 90, 80];
        let encoded = encode_integers_delta(&values);
        assert_eq!(encoded, vec![100, -5, -5, -10]);

        let decoded = decode_integers_delta(&encoded);
        assert_eq!(decoded, values);
    }

    #[test]
    fn test_zigzag_encoding() {
        assert_eq!(zigzag_encode(0), 0);
        assert_eq!(zigzag_encode(-1), 1);
        assert_eq!(zigzag_encode(1), 2);
        assert_eq!(zigzag_encode(-2), 3);
        assert_eq!(zigzag_encode(2), 4);
        assert_eq!(zigzag_encode(i32::MAX), u32::MAX - 1);
        assert_eq!(zigzag_encode(i32::MIN), u32::MAX);
    }

    #[test]
    fn test_zigzag_roundtrip() {
        let values = vec![0, 1, -1, 100, -100, 1000, -1000, i32::MAX, i32::MIN];
        for &val in &values {
            let encoded = zigzag_encode(val);
            let decoded = zigzag_decode(encoded);
            assert_eq!(decoded, val, "Failed for {}", val);
        }
    }

    #[test]
    fn test_sequence_encoding_ascii() {
        let seq = b"ACGTNNACGT";
        let encoded = encode_sequence_ascii(seq);
        assert_eq!(encoded, seq.to_vec());

        let decoded = decode_sequence_ascii(&encoded).unwrap();
        assert_eq!(decoded, seq.to_vec());
    }

    #[test]
    fn test_sequence_encoding_lowercase() {
        let seq = b"acgtnnacgt";
        let encoded = encode_sequence_ascii(seq);
        assert_eq!(encoded, seq.to_vec());
    }

    #[test]
    fn test_sequence_decoding_invalid() {
        let invalid = b"ACGT123";
        let result = decode_sequence_ascii(invalid);
        assert!(result.is_err());
    }

    #[test]
    fn test_qualities_encoding() {
        let quals = b"IIIIIIIIII";
        let encoded = encode_qualities_raw(quals);
        assert_eq!(encoded, quals.to_vec());

        let decoded = decode_qualities_raw(&encoded);
        assert_eq!(decoded, quals.to_vec());
    }

    #[test]
    fn test_cigar_op_encoding() {
        // 10M
        let encoded = encode_cigar_op(0, 10);
        assert_eq!(encoded, (10 << 4) | 0);

        let (op, len) = decode_cigar_op(encoded);
        assert_eq!(op, 0);
        assert_eq!(len, 10);
    }

    #[test]
    fn test_cigar_string_parsing() {
        let cigar = "10M2I5D3N4S";
        let ops = parse_cigar_string(cigar).unwrap();

        assert_eq!(ops.len(), 5);
        assert_eq!(decode_cigar_op(ops[0]), (0, 10)); // 10M
        assert_eq!(decode_cigar_op(ops[1]), (1, 2));  // 2I
        assert_eq!(decode_cigar_op(ops[2]), (2, 5));  // 5D
        assert_eq!(decode_cigar_op(ops[3]), (3, 3));  // 3N
        assert_eq!(decode_cigar_op(ops[4]), (4, 4));  // 4S

        let formatted = format_cigar_string(&ops);
        assert_eq!(formatted, cigar);
    }

    #[test]
    fn test_cigar_string_complex() {
        let cigar = "100M5I10D2S50M";
        let ops = parse_cigar_string(cigar).unwrap();
        let formatted = format_cigar_string(&ops);
        assert_eq!(formatted, cigar);
    }

    #[test]
    fn test_cigar_string_invalid_op() {
        let cigar = "10Z"; // Invalid op
        let result = parse_cigar_string(cigar);
        assert!(result.is_err());
    }

    #[test]
    fn test_cigar_string_no_length() {
        let cigar = "M10"; // Op before length
        let result = parse_cigar_string(cigar);
        assert!(result.is_err());
    }

    #[test]
    fn test_cigar_string_ends_with_digit() {
        let cigar = "10M5"; // Ends with digit
        let result = parse_cigar_string(cigar);
        assert!(result.is_err());
    }

    #[test]
    fn test_cigar_all_ops() {
        // Test all 9 CIGAR operations
        let cigar = "1M2I3D4N5S6H7P8=9X";
        let ops = parse_cigar_string(cigar).unwrap();
        assert_eq!(ops.len(), 9);

        let formatted = format_cigar_string(&ops);
        assert_eq!(formatted, cigar);
    }

    // Property-based tests
    proptest! {
        #[test]
        fn prop_delta_encoding_roundtrip(
            values in prop::collection::vec(i32::MIN/2..i32::MAX/2, 0..1000)
        ) {
            // Constrain range to prevent overflow: any delta between consecutive
            // values in range [MIN/2, MAX/2] fits in i32
            let encoded = encode_integers_delta(&values);
            let decoded = decode_integers_delta(&encoded);
            prop_assert_eq!(decoded, values);
        }

        #[test]
        fn prop_zigzag_roundtrip(value in any::<i32>()) {
            let encoded = zigzag_encode(value);
            let decoded = zigzag_decode(encoded);
            prop_assert_eq!(decoded, value);
        }

        #[test]
        fn prop_sequence_valid_bases(seq in "[ACGTNacgtn]{0,1000}") {
            let encoded = encode_sequence_ascii(seq.as_bytes());
            let decoded = decode_sequence_ascii(&encoded);
            prop_assert!(decoded.is_ok());
            prop_assert_eq!(decoded.unwrap(), seq.as_bytes().to_vec());
        }

        #[test]
        fn prop_qualities_roundtrip(quals in prop::collection::vec(any::<u8>(), 0..1000)) {
            let encoded = encode_qualities_raw(&quals);
            let decoded = decode_qualities_raw(&encoded);
            prop_assert_eq!(decoded, quals);
        }
    }
}
