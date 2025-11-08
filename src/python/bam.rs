//! Python bindings for BAM/SAM parser
//!
//! Provides Python access to biometal's high-performance BAM parser with:
//! - Parallel BGZF decompression (4× speedup, Rule 3)
//! - Streaming architecture (constant ~5 MB memory, Rule 5)
//! - 4.54 million records/sec throughput (43.0 MiB/s)
//!
//! # Evidence Base
//!
//! - Rule 3 (Parallel BGZF): Entry 029, 6.5× validated speedup
//! - Rule 5 (Streaming): Entry 026, constant memory architecture
//! - Performance: 4× overall speedup (experiments/native-bam-implementation/PHASE_3_BENCHMARKS.md)

use pyo3::prelude::*;
use pyo3::exceptions::PyStopIteration;
use std::path::PathBuf;
use crate::io::bam::{BamReader, Header};
use crate::io::compression::CompressedReader;

/// BAM record with alignment information
///
/// Attributes:
///     name (str): Read name/identifier
///     reference_id (int | None): Reference sequence ID (-1 for unmapped)
///     position (int | None): 0-based leftmost position (-1 for unmapped)
///     mapq (int | None): Mapping quality (255 for unavailable)
///     flags (int): SAM flags (bitwise)
///     mate_reference_id (int | None): Mate reference ID
///     mate_position (int | None): Mate position
///     template_length (int): Template length
///     sequence (bytes): Read sequence
///     quality (bytes): Phred quality scores
///
/// Example:
///     >>> bam = biometal.BamReader.from_path("alignments.bam")
///     >>> for record in bam:
///     ...     if record.mapq and record.mapq >= 30:
///     ...         print(f"{record.name}: chr{record.reference_id}:{record.position}")
///
/// Note:
///     CIGAR and tags are available as raw data. Full parsing coming in future version.
#[pyclass(name = "BamRecord")]
#[derive(Clone)]
pub struct PyBamRecord {
    /// Read name/identifier
    #[pyo3(get)]
    pub name: String,

    /// Reference sequence ID (None for unmapped)
    #[pyo3(get)]
    pub reference_id: Option<usize>,

    /// 0-based leftmost position (None for unmapped)
    #[pyo3(get)]
    pub position: Option<i32>,

    /// Mapping quality (None if unavailable)
    #[pyo3(get)]
    pub mapq: Option<u8>,

    /// SAM flags (bitwise)
    #[pyo3(get)]
    pub flags: u16,

    /// Mate reference ID (None for unmapped mate)
    #[pyo3(get)]
    pub mate_reference_id: Option<usize>,

    /// Mate position (None for unmapped mate)
    #[pyo3(get)]
    pub mate_position: Option<i32>,

    /// Template length
    #[pyo3(get)]
    pub template_length: i32,

    /// Read sequence
    #[pyo3(get)]
    pub sequence: Vec<u8>,

    /// Phred quality scores
    #[pyo3(get)]
    pub quality: Vec<u8>,
}

#[pymethods]
impl PyBamRecord {
    /// Check if read is mapped
    #[getter]
    fn is_mapped(&self) -> bool {
        (self.flags & 0x4) == 0
    }

    /// Check if read is reverse complement
    #[getter]
    fn is_reverse(&self) -> bool {
        (self.flags & 0x10) != 0
    }

    /// Check if read is first in pair
    #[getter]
    fn is_first(&self) -> bool {
        (self.flags & 0x40) != 0
    }

    /// Check if read is second in pair
    #[getter]
    fn is_second(&self) -> bool {
        (self.flags & 0x80) != 0
    }

    /// Check if read is paired
    #[getter]
    fn is_paired(&self) -> bool {
        (self.flags & 0x1) != 0
    }

    /// Check if alignment is primary
    #[getter]
    fn is_primary(&self) -> bool {
        (self.flags & 0x900) == 0  // Not secondary, not supplementary
    }

    /// Get sequence as string
    #[getter]
    fn sequence_str(&self) -> String {
        String::from_utf8_lossy(&self.sequence).to_string()
    }

    /// Get quality as string
    #[getter]
    fn quality_str(&self) -> String {
        String::from_utf8_lossy(&self.quality).to_string()
    }

    fn __repr__(&self) -> String {
        let ref_str = self.reference_id
            .map(|r| r.to_string())
            .unwrap_or_else(|| "*".to_string());
        let pos_str = self.position
            .map(|p| p.to_string())
            .unwrap_or_else(|| "*".to_string());
        let mapq_str = self.mapq
            .map(|m| m.to_string())
            .unwrap_or_else(|| "*".to_string());

        format!(
            "BamRecord(name='{}', ref={}, pos={}, mapq={}, flags={})",
            self.name, ref_str, pos_str, mapq_str, self.flags
        )
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

impl From<crate::io::bam::Record> for PyBamRecord {
    fn from(record: crate::io::bam::Record) -> Self {
        PyBamRecord {
            name: record.name,
            reference_id: record.reference_id,
            position: record.position,
            mapq: record.mapq,
            flags: record.flags,
            mate_reference_id: record.mate_reference_id,
            mate_position: record.mate_position,
            template_length: record.template_length,
            sequence: record.sequence,
            quality: record.quality,
        }
    }
}

/// Stream BAM records with constant memory and parallel BGZF decompression
///
/// High-performance streaming BAM parser with automatic parallel BGZF decompression.
/// Maintains constant ~5 MB memory regardless of file size (terabyte-scale capable).
///
/// Args:
///     path (str): Path to BAM file (.bam, .bam.gz, or uncompressed)
///
/// Performance:
///     - 4.54 million records/sec throughput
///     - 43.0 MiB/s compressed file processing
///     - 4× speedup via parallel BGZF decompression
///     - Constant ~5 MB memory (streams terabyte-scale files)
///
/// Example:
///     >>> import biometal
///     >>> bam = biometal.BamReader.from_path("alignments.bam")
///     >>>
///     >>> # Access header info
///     >>> print(f"References: {bam.reference_count}")
///     >>>
///     >>> # Stream records with constant memory
///     >>> mapped_count = 0
///     >>> for record in bam:
///     ...     if record.is_mapped and record.mapq and record.mapq >= 30:
///     ...         mapped_count += 1
///     ...         # Process high-quality alignments
///     >>>
///     >>> print(f"High-quality alignments: {mapped_count}")
///
/// Note:
///     Automatically detects and uses parallel BGZF decompression for compressed BAM files.
///     Memory footprint remains constant at ~5 MB even for TB-scale alignments.
#[pyclass(name = "BamReader", unsendable)]
pub struct PyBamReader {
    inner: Option<BamReader<CompressedReader>>,
    header: PyBamHeader,
}

/// BAM header information
///
/// Attributes:
///     text (str): SAM header text
///     reference_count (int): Number of reference sequences
///
/// Note:
///     Reference sequences can be accessed via the header text.
///     Full reference sequence API coming in future version.
#[pyclass(name = "BamHeader")]
#[derive(Clone)]
pub struct PyBamHeader {
    /// SAM header text
    #[pyo3(get)]
    pub text: String,

    /// Number of reference sequences
    #[pyo3(get)]
    pub reference_count: usize,
}

#[pymethods]
impl PyBamHeader {
    fn __repr__(&self) -> String {
        format!(
            "BamHeader(references={}, text_len={})",
            self.reference_count,
            self.text.len()
        )
    }
}

impl From<&Header> for PyBamHeader {
    fn from(header: &Header) -> Self {
        PyBamHeader {
            text: header.text.clone(),
            reference_count: header.reference_count(),
        }
    }
}

#[pymethods]
impl PyBamReader {
    /// Open BAM file with automatic parallel BGZF decompression
    ///
    /// Args:
    ///     path (str): Path to BAM file
    ///
    /// Returns:
    ///     BamReader: Streaming iterator with ~5 MB constant memory
    ///
    /// Raises:
    ///     IOError: If file cannot be opened
    ///     ValueError: If file format is invalid
    ///
    /// Performance:
    ///     Automatically uses parallel BGZF decompression (4× speedup).
    ///     Processes 4.54 million records/sec with constant memory.
    ///
    /// Example:
    ///     >>> bam = biometal.BamReader.from_path("alignments.bam")
    ///     >>> print(f"Opened BAM with {bam.reference_count} references")
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let reader = BamReader::from_path(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        // Clone header info for Python access
        let header = PyBamHeader::from(reader.header());

        Ok(PyBamReader {
            inner: Some(reader),
            header,
        })
    }

    /// Get BAM header
    ///
    /// Returns:
    ///     BamHeader: Header with reference information
    #[getter]
    fn header(&self) -> PyBamHeader {
        self.header.clone()
    }

    /// Get number of reference sequences
    ///
    /// Returns:
    ///     int: Number of reference sequences in BAM header
    #[getter]
    fn reference_count(&self) -> usize {
        self.header.reference_count
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyBamRecord> {
        if let Some(ref mut reader) = slf.inner {
            match reader.read_record() {
                Ok(Some(record)) => Ok(record.into()),
                Ok(None) => Err(PyStopIteration::new_err("no more records")),
                Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())),
            }
        } else {
            Err(PyStopIteration::new_err("reader exhausted"))
        }
    }

    fn __repr__(&self) -> String {
        format!("BamReader(references={})", self.header.reference_count)
    }
}
