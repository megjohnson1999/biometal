//! Python bindings for masking operations
//!
//! Provides Python access to biometal's quality-based masking functions
//! (Phase 4: replace low-quality bases with 'N').

use pyo3::prelude::*;
use crate::operations::masking;
use crate::python::records::PyFastqRecord;

/// Mask low-quality bases in a record (in-place for Python, returns new record)
///
/// Replaces bases with quality below threshold with 'N'. Unlike trimming,
/// this preserves read length.
///
/// Args:
///     record (FastqRecord): Input FASTQ record
///     min_quality (int): Minimum quality threshold (Phred score)
///
/// Returns:
///     FastqRecord: Record with low-quality bases masked as 'N'
///
/// Raises:
///     ValueError: If sequence and quality length mismatch
///
/// Example:
///     >>> import biometal
///     >>> stream = biometal.FastqStream.from_path("data.fq")
///     >>> for record in stream:
///     ...     # Mask bases with Q < 20
///     ...     masked = biometal.mask_low_quality(record, 20)
///     ...     print(masked.sequence)  # Will contain 'N' for low-quality bases
///     ...     break
///
/// Use Case:
///     Variant calling pipelines often prefer masking over trimming to
///     preserve read structure and alignment positions.
///
/// Note:
///     In Python, this returns a new record (Rust in-place not possible due
///     to Python bytes immutability). For explicit copy semantics, this is
///     equivalent to the Rust `mask_low_quality_copy()`.
#[pyfunction(name = "mask_low_quality")]
pub fn py_mask_low_quality(record: &PyFastqRecord, min_quality: u8) -> PyResult<PyFastqRecord> {
    masking::mask_low_quality_copy(&record.to_fastq_record(), min_quality)
        .map(PyFastqRecord::from)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
}

/// Count masked bases (N's) in a record
///
/// Returns the number of 'N' or 'n' bases in the sequence. Useful for
/// QC metrics after masking.
///
/// Args:
///     record (FastqRecord): Input FASTQ record
///
/// Returns:
///     int: Number of masked (N) bases
///
/// Example:
///     >>> import biometal
///     >>> stream = biometal.FastqStream.from_path("data.fq")
///     >>> for record in stream:
///     ...     masked = biometal.mask_low_quality(record, 20)
///     ...     n_count = biometal.count_masked_bases(masked)
///     ...     if n_count / len(masked.sequence) < 0.1:  # <10% masked
///     ...         print("Pass QC")
///     ...     break
///
/// Note:
///     Counts both 'N' and 'n' (case-insensitive).
#[pyfunction(name = "count_masked_bases")]
pub fn py_count_masked_bases(record: &PyFastqRecord) -> usize {
    masking::count_masked_bases(&record.to_fastq_record())
}
