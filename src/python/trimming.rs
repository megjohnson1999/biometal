//! Python bindings for trimming operations
//!
//! Provides Python access to biometal's read trimming functions
//! (Phase 4: fixed position and quality-based trimming).

use pyo3::prelude::*;
use crate::operations::trimming;
use crate::python::records::PyFastqRecord;

/// Trim bases from the start of a record
///
/// Removes the specified number of bases from the start (5' end).
///
/// Args:
///     record (FastqRecord): Input FASTQ record
///     bases (int): Number of bases to trim from start
///
/// Returns:
///     FastqRecord: Trimmed record
///
/// Raises:
///     ValueError: If trimming would result in empty record
///
/// Example:
///     >>> import biometal
///     >>> stream = biometal.FastqStream.from_path("data.fq")
///     >>> for record in stream:
///     ...     # Remove first 10 bases (e.g., adapter)
///     ...     trimmed = biometal.trim_start(record, 10)
///     ...     break
#[pyfunction(name = "trim_start")]
pub fn py_trim_start(record: &PyFastqRecord, bases: usize) -> PyResult<PyFastqRecord> {
    trimming::trim_start(&record.to_fastq_record(), bases)
        .map(PyFastqRecord::from)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
}

/// Trim bases from the end of a record
///
/// Removes the specified number of bases from the end (3' end).
///
/// Args:
///     record (FastqRecord): Input FASTQ record
///     bases (int): Number of bases to trim from end
///
/// Returns:
///     FastqRecord: Trimmed record
///
/// Raises:
///     ValueError: If trimming would result in empty record
///
/// Example:
///     >>> import biometal
///     >>> stream = biometal.FastqStream.from_path("data.fq")
///     >>> for record in stream:
///     ...     # Remove last 5 bases
///     ...     trimmed = biometal.trim_end(record, 5)
///     ...     break
#[pyfunction(name = "trim_end")]
pub fn py_trim_end(record: &PyFastqRecord, bases: usize) -> PyResult<PyFastqRecord> {
    trimming::trim_end(&record.to_fastq_record(), bases)
        .map(PyFastqRecord::from)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
}

/// Trim bases from both ends of a record
///
/// Removes specified bases from both start (5') and end (3').
///
/// Args:
///     record (FastqRecord): Input FASTQ record
///     start_bases (int): Number of bases to trim from start
///     end_bases (int): Number of bases to trim from end
///
/// Returns:
///     FastqRecord: Trimmed record
///
/// Raises:
///     ValueError: If trimming would result in empty record
///
/// Example:
///     >>> import biometal
///     >>> stream = biometal.FastqStream.from_path("data.fq")
///     >>> for record in stream:
///     ...     # Trim 10 from start, 5 from end
///     ...     trimmed = biometal.trim_both(record, 10, 5)
///     ...     break
#[pyfunction(name = "trim_both")]
pub fn py_trim_both(
    record: &PyFastqRecord,
    start_bases: usize,
    end_bases: usize,
) -> PyResult<PyFastqRecord> {
    trimming::trim_both(&record.to_fastq_record(), start_bases, end_bases)
        .map(PyFastqRecord::from)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
}

/// Trim low-quality bases from the end
///
/// Removes bases from the 3' end while quality is below threshold.
/// Uses Phred+33 encoding (Illumina 1.8+).
///
/// Args:
///     record (FastqRecord): Input FASTQ record
///     min_quality (int): Minimum quality threshold (Phred score)
///
/// Returns:
///     FastqRecord: Trimmed record (may be empty if all bases low quality)
///
/// Raises:
///     ValueError: If sequence and quality length mismatch
///
/// Example:
///     >>> import biometal
///     >>> stream = biometal.FastqStream.from_path("data.fq")
///     >>> for record in stream:
///     ...     # Trim bases with Q < 20 from end
///     ...     trimmed = biometal.trim_quality_end(record, 20)
///     ...     if not trimmed.is_empty():
///     ...         print("Kept trimmed read")
///     ...     break
///
/// Note:
///     Common quality threshold values:
///     - Q20: 99% base call accuracy
///     - Q30: 99.9% base call accuracy
#[pyfunction(name = "trim_quality_end")]
pub fn py_trim_quality_end(record: &PyFastqRecord, min_quality: u8) -> PyResult<PyFastqRecord> {
    trimming::trim_quality_end(&record.to_fastq_record(), min_quality)
        .map(PyFastqRecord::from)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
}

/// Trim low-quality bases from the start
///
/// Removes bases from the 5' end while quality is below threshold.
/// Uses Phred+33 encoding (Illumina 1.8+).
///
/// Args:
///     record (FastqRecord): Input FASTQ record
///     min_quality (int): Minimum quality threshold (Phred score)
///
/// Returns:
///     FastqRecord: Trimmed record (may be empty if all bases low quality)
///
/// Raises:
///     ValueError: If sequence and quality length mismatch
///
/// Example:
///     >>> import biometal
///     >>> stream = biometal.FastqStream.from_path("data.fq")
///     >>> for record in stream:
///     ...     # Trim bases with Q < 20 from start
///     ...     trimmed = biometal.trim_quality_start(record, 20)
///     ...     break
#[pyfunction(name = "trim_quality_start")]
pub fn py_trim_quality_start(record: &PyFastqRecord, min_quality: u8) -> PyResult<PyFastqRecord> {
    trimming::trim_quality_start(&record.to_fastq_record(), min_quality)
        .map(PyFastqRecord::from)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
}

/// Trim low-quality bases from both ends
///
/// Removes bases from both 5' and 3' ends while quality is below threshold.
/// Uses Phred+33 encoding (Illumina 1.8+).
///
/// Args:
///     record (FastqRecord): Input FASTQ record
///     min_quality (int): Minimum quality threshold (Phred score)
///
/// Returns:
///     FastqRecord: Trimmed record (may be empty if all bases low quality)
///
/// Raises:
///     ValueError: If sequence and quality length mismatch
///
/// Example:
///     >>> import biometal
///     >>> stream = biometal.FastqStream.from_path("data.fq")
///     >>> for record in stream:
///     ...     # Trim Q < 20 from both ends
///     ...     trimmed = biometal.trim_quality_both(record, 20)
///     ...     if len(trimmed.sequence) >= 50:  # Keep if ≥50bp remain
///     ...         print("Pass QC")
///     ...     break
///
/// Note:
///     Single-pass optimized implementation (Phase 4 improvement).
#[pyfunction(name = "trim_quality_both")]
pub fn py_trim_quality_both(record: &PyFastqRecord, min_quality: u8) -> PyResult<PyFastqRecord> {
    trimming::trim_quality_both(&record.to_fastq_record(), min_quality)
        .map(PyFastqRecord::from)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
}

/// Trim using sliding window quality (Trimmomatic-style)
///
/// Scans with a sliding window and trims when average quality drops below
/// threshold. More aggressive than simple end trimming.
///
/// Args:
///     record (FastqRecord): Input FASTQ record
///     min_quality (int): Minimum average quality threshold
///     window_size (int): Window size for averaging (e.g., 4)
///
/// Returns:
///     FastqRecord: Trimmed record (may be empty if all windows fail)
///
/// Raises:
///     ValueError: If sequence and quality length mismatch or window_size = 0
///
/// Example:
///     >>> import biometal
///     >>> stream = biometal.FastqStream.from_path("data.fq")
///     >>> for record in stream:
///     ...     # Trimmomatic-style: 4bp window, Q20 average
///     ...     trimmed = biometal.trim_quality_window(record, 20, 4)
///     ...     break
///
/// Note:
///     This is the Trimmomatic SLIDINGWINDOW algorithm. Window slides from
///     5' to 3', trimming when average quality drops below threshold.
#[pyfunction(name = "trim_quality_window")]
pub fn py_trim_quality_window(
    record: &PyFastqRecord,
    min_quality: u8,
    window_size: usize,
) -> PyResult<PyFastqRecord> {
    trimming::trim_quality_window(&record.to_fastq_record(), min_quality, window_size)
        .map(PyFastqRecord::from)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
}
