//! Python bindings for record-level operations
//!
//! Provides Python access to FastqRecord manipulation functions
//! (Phase 4: extract_region, reverse_complement_record, length filtering).

use pyo3::prelude::*;
use crate::operations::record_ops;
use crate::python::records::PyFastqRecord;
use crate::python::records::PyFastaRecord;

/// Extract a region from a FASTQ record
///
/// Returns a new record containing only the specified region [start, end).
/// Both sequence and quality scores are extracted.
///
/// Args:
///     record (FastqRecord): Input FASTQ record
///     start (int): Start position (0-based, inclusive)
///     end (int): End position (0-based, exclusive)
///
/// Returns:
///     FastqRecord: New record with extracted region
///
/// Raises:
///     ValueError: If start >= end or range is out of bounds
///
/// Example:
///     >>> import biometal
///     >>> stream = biometal.FastqStream.from_path("data.fq")
///     >>> for record in stream:
///     ...     # Extract bases 10-50
///     ...     region = biometal.extract_region(record, 10, 50)
///     ...     break
///
/// Note:
///     Quality scores are preserved and extracted with the sequence.
#[pyfunction(name = "extract_region")]
pub fn py_extract_region(
    record: &PyFastqRecord,
    start: usize,
    end: usize,
) -> PyResult<PyFastqRecord> {
    record_ops::extract_region(&record.to_fastq_record(), start, end)
        .map(PyFastqRecord::from)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
}

/// Reverse complement a FASTQ record
///
/// Returns a new record with reverse complemented sequence and reversed
/// quality scores (to maintain alignment).
///
/// Args:
///     record (FastqRecord): Input FASTQ record
///
/// Returns:
///     FastqRecord: New record with reverse complemented sequence
///
/// Example:
///     >>> import biometal
///     >>> stream = biometal.FastqStream.from_path("data.fq")
///     >>> for record in stream:
///     ...     rc = biometal.reverse_complement_record(record)
///     ...     print(rc.sequence)
///     ...     break
///
/// Note:
///     Quality scores are reversed to maintain position alignment.
///     For sequence-only RC, use `reverse_complement(record.sequence)`.
#[pyfunction(name = "reverse_complement_record")]
pub fn py_reverse_complement_record(record: &PyFastqRecord) -> PyFastqRecord {
    PyFastqRecord::from(record_ops::reverse_complement_record(&record.to_fastq_record()))
}

/// Get sequence length from a FASTQ record
///
/// Returns the length of the sequence in the record.
///
/// Args:
///     record (FastqRecord): Input FASTQ record
///
/// Returns:
///     int: Length of sequence
///
/// Example:
///     >>> import biometal
///     >>> stream = biometal.FastqStream.from_path("data.fq")
///     >>> for record in stream:
///     ...     length = biometal.sequence_length(record)
///     ...     print(f"Read length: {length}")
///     ...     break
///
/// Note:
///     This is equivalent to `len(record.sequence)` but provided for
///     API completeness.
#[pyfunction(name = "sequence_length")]
pub fn py_sequence_length(record: &PyFastqRecord) -> usize {
    record_ops::sequence_length(&record.to_fastq_record())
}

/// Check if record meets length requirements
///
/// Returns True if record length is within [min_len, max_len] inclusive.
///
/// Args:
///     record (FastqRecord): Input FASTQ record
///     min_len (int): Minimum length (inclusive)
///     max_len (int): Maximum length (inclusive)
///
/// Returns:
///     bool: True if length requirement met, False otherwise
///
/// Example:
///     >>> import biometal
///     >>> stream = biometal.FastqStream.from_path("data.fq")
///     >>> for record in stream:
///     ...     # Keep only reads 50-150 bp
///     ...     if biometal.meets_length_requirement(record, 50, 150):
///     ...         print("Keep this read")
///     ...     break
///
/// Note:
///     Common in quality control pipelines to filter reads by length.
#[pyfunction(name = "meets_length_requirement")]
pub fn py_meets_length_requirement(
    record: &PyFastqRecord,
    min_len: usize,
    max_len: usize,
) -> bool {
    record_ops::meets_length_requirement(&record.to_fastq_record(), min_len, max_len)
}

/// Convert FASTQ record to FASTA record
///
/// Drops quality scores and returns a FASTA record with just ID and sequence.
///
/// Args:
///     record (FastqRecord): Input FASTQ record
///
/// Returns:
///     FastaRecord: FASTA record (no quality scores)
///
/// Example:
///     >>> import biometal
///     >>> fastq_stream = biometal.FastqStream.from_path("data.fq")
///     >>> for record in fastq_stream:
///     ...     fasta = biometal.to_fasta_record(record)
///     ...     print(f">{fasta.id}")
///     ...     print(fasta.sequence)
///     ...     break
///
/// Use Case:
///     Converting FASTQ to FASTA after quality filtering.
#[pyfunction(name = "to_fasta_record")]
pub fn py_to_fasta_record(record: &PyFastqRecord) -> PyFastaRecord {
    PyFastaRecord::from(record_ops::to_fasta_record(&record.to_fastq_record()))
}
