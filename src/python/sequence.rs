//! Python bindings for sequence operations
//!
//! Provides Python access to biometal's sequence manipulation primitives
//! (Phase 4: reverse complement, complement, reverse, validation).

use pyo3::prelude::*;
use crate::operations::sequence;

/// Reverse complement a DNA/RNA sequence
///
/// Returns the reverse complement of the input sequence. Supports both DNA
/// (A, C, G, T) and RNA (A, C, G, U) with full IUPAC ambiguity codes.
///
/// Args:
///     sequence (bytes): DNA or RNA sequence
///
/// Returns:
///     bytes: Reverse complemented sequence
///
/// Example:
///     >>> import biometal
///     >>> biometal.reverse_complement(b"ATGC")
///     b'GCAT'
///     >>> biometal.reverse_complement(b"AUGC")  # RNA
///     b'GCAU'
///
/// Note:
///     This is a scalar implementation (Phase 4). NEON deferred based on
///     evidence (estimated <2× speedup, fails ≥5× threshold).
#[pyfunction(name = "reverse_complement")]
pub fn py_reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence::reverse_complement(sequence)
}

/// Complement a DNA/RNA sequence
///
/// Returns the complement (without reversing) of the input sequence.
/// A ↔ T/U, C ↔ G, with full IUPAC support.
///
/// Args:
///     sequence (bytes): DNA or RNA sequence
///
/// Returns:
///     bytes: Complemented sequence
///
/// Example:
///     >>> import biometal
///     >>> biometal.complement(b"ATGC")
///     b'TACG'
///
/// Note:
///     For reverse complement, use `reverse_complement()` which is more
///     commonly needed in bioinformatics.
#[pyfunction(name = "complement")]
pub fn py_complement(sequence: &[u8]) -> Vec<u8> {
    sequence::complement(sequence)
}

/// Reverse a sequence
///
/// Returns the reversed sequence (without complementing).
///
/// Args:
///     sequence (bytes): Any sequence
///
/// Returns:
///     bytes: Reversed sequence
///
/// Example:
///     >>> import biometal
///     >>> biometal.reverse(b"ATGC")
///     b'CGTA'
///
/// Note:
///     Benchmark shows this is 3× faster than reverse_complement
///     (11.1 GiB/s vs 3.7 GiB/s), proving table lookup is the bottleneck.
#[pyfunction(name = "reverse")]
pub fn py_reverse(sequence: &[u8]) -> Vec<u8> {
    sequence::reverse(sequence)
}

/// Check if sequence is valid DNA
///
/// Returns True if sequence contains only valid DNA bases: A, C, G, T, N
/// and IUPAC ambiguity codes (case-insensitive).
///
/// Args:
///     sequence (bytes): Sequence to validate
///
/// Returns:
///     bool: True if valid DNA, False otherwise
///
/// Example:
///     >>> import biometal
///     >>> biometal.is_valid_dna(b"ATGCN")
///     True
///     >>> biometal.is_valid_dna(b"AUGC")  # U is RNA, not DNA
///     False
///
/// Note:
///     For RNA validation, use `is_valid_rna()`.
#[pyfunction(name = "is_valid_dna")]
pub fn py_is_valid_dna(sequence: &[u8]) -> bool {
    sequence::is_valid_dna(sequence)
}

/// Check if sequence is valid RNA
///
/// Returns True if sequence contains only valid RNA bases: A, C, G, U, N
/// and IUPAC ambiguity codes (case-insensitive).
///
/// Args:
///     sequence (bytes): Sequence to validate
///
/// Returns:
///     bool: True if valid RNA, False otherwise
///
/// Example:
///     >>> import biometal
///     >>> biometal.is_valid_rna(b"AUGCN")
///     True
///     >>> biometal.is_valid_rna(b"ATGC")  # T is DNA, not RNA
///     False
///
/// Warning:
///     RNA support has limitations. Reverse complement works but may not
///     handle all RNA-specific modifications. For production RNA analysis,
///     validate results carefully.
#[pyfunction(name = "is_valid_rna")]
pub fn py_is_valid_rna(sequence: &[u8]) -> bool {
    sequence::is_valid_rna(sequence)
}

/// Count invalid bases in sequence
///
/// Returns the number of bases that are not valid IUPAC codes.
/// Useful for quality checking sequences.
///
/// Args:
///     sequence (bytes): Sequence to check
///
/// Returns:
///     int: Number of invalid bases
///
/// Example:
///     >>> import biometal
///     >>> biometal.count_invalid_bases(b"ATGC")
///     0
///     >>> biometal.count_invalid_bases(b"ATGC123")
///     3
///
/// Note:
///     Valid IUPAC codes: A, C, G, T, U, N, R, Y, S, W, K, M, B, D, H, V
///     (case-insensitive)
#[pyfunction(name = "count_invalid_bases")]
pub fn py_count_invalid_bases(sequence: &[u8]) -> usize {
    sequence::count_invalid_bases(sequence)
}
