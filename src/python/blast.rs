//! Python bindings for BLAST tabular format (outfmt 6/7)

use pyo3::prelude::*;
use pyo3::exceptions::PyStopIteration;
use crate::formats::blast::{BlastRecord, BlastTabularParser};
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

/// BLAST tabular alignment record (outfmt 6/7)
///
/// Represents a single alignment between a query and subject sequence from BLAST output.
///
/// Attributes:
///     qseqid (str): Query sequence ID
///     sseqid (str): Subject (reference) sequence ID
///     pident (float): Percentage of identical matches (0-100)
///     length (int): Alignment length
///     mismatch (int): Number of mismatches
///     gapopen (int): Number of gap openings
///     qstart (int): Query start position (1-based)
///     qend (int): Query end position (1-based, inclusive)
///     sstart (int): Subject start position (1-based)
///     send (int): Subject end position (1-based, inclusive)
///     evalue (float): Expect value (statistical significance, lower is better)
///     bitscore (float): Bit score (alignment quality, higher is better)
///
/// Example:
///     >>> for record in parser:
///     ...     if record.evalue < 1e-10 and record.pident > 95.0:
///     ...         print(f"{record.qseqid} -> {record.sseqid}: {record.pident:.1f}% identity")
#[pyclass(name = "BlastRecord")]
#[derive(Clone)]
pub struct PyBlastRecord {
    #[pyo3(get)]
    pub qseqid: String,
    #[pyo3(get)]
    pub sseqid: String,
    #[pyo3(get)]
    pub pident: f64,
    #[pyo3(get)]
    pub length: u32,
    #[pyo3(get)]
    pub mismatch: u32,
    #[pyo3(get)]
    pub gapopen: u32,
    #[pyo3(get)]
    pub qstart: u32,
    #[pyo3(get)]
    pub qend: u32,
    #[pyo3(get)]
    pub sstart: u32,
    #[pyo3(get)]
    pub send: u32,
    #[pyo3(get)]
    pub evalue: f64,
    #[pyo3(get)]
    pub bitscore: f64,
}

impl From<BlastRecord> for PyBlastRecord {
    fn from(r: BlastRecord) -> Self {
        PyBlastRecord {
            qseqid: r.qseqid,
            sseqid: r.sseqid,
            pident: r.pident,
            length: r.length,
            mismatch: r.mismatch,
            gapopen: r.gapopen,
            qstart: r.qstart,
            qend: r.qend,
            sstart: r.sstart,
            send: r.send,
            evalue: r.evalue,
            bitscore: r.bitscore,
        }
    }
}

#[pymethods]
impl PyBlastRecord {
    fn __repr__(&self) -> String {
        format!(
            "BlastRecord(qseqid='{}', sseqid='{}', pident={:.2}, evalue={:.2e})",
            self.qseqid, self.sseqid, self.pident, self.evalue
        )
    }

    /// Calculate percent identity as a fraction (0.0-1.0)
    ///
    /// Returns:
    ///     float: Percent identity divided by 100
    ///
    /// Example:
    ///     >>> identity = record.identity()  # 0.98 for 98% identity
    fn identity(&self) -> f64 {
        self.pident / 100.0
    }

    /// Check if alignment is high quality (low evalue, high identity)
    ///
    /// Args:
    ///     min_pident (float): Minimum percent identity (0-100)
    ///     max_evalue (float): Maximum E-value
    ///
    /// Returns:
    ///     bool: True if both criteria are met
    ///
    /// Example:
    ///     >>> if record.is_high_quality(min_pident=95.0, max_evalue=1e-10):
    ///     ...     print("High quality alignment")
    fn is_high_quality(&self, min_pident: f64, max_evalue: f64) -> bool {
        self.pident >= min_pident && self.evalue <= max_evalue
    }

    /// Calculate query alignment coverage
    ///
    /// Returns:
    ///     int: Number of query bases aligned
    ///
    /// Example:
    ///     >>> coverage = record.query_coverage()
    fn query_coverage(&self) -> u32 {
        self.qend.saturating_sub(self.qstart) + 1
    }

    /// Calculate subject alignment coverage
    ///
    /// Returns:
    ///     int: Number of subject bases aligned
    ///
    /// Example:
    ///     >>> coverage = record.subject_coverage()
    fn subject_coverage(&self) -> u32 {
        self.send.saturating_sub(self.sstart) + 1
    }
}

/// Read BLAST tabular files with streaming architecture
///
/// Streaming parser for BLAST tabular output formats (outfmt 6 and 7).
/// Provides constant memory usage through iterator-based architecture.
///
/// Args:
///     path (str): Path to BLAST tabular file (.blast, .outfmt6, .outfmt7, .txt)
///
/// Features:
///     - Streaming architecture (constant memory)
///     - Supports both outfmt 6 (no headers) and outfmt 7 (with comments)
///     - Automatic comment line skipping
///     - 12-column standard BLAST format
///
/// Example:
///     >>> import biometal
///     >>> parser = biometal.BlastTabularParser.from_path("alignments.blast")
///     >>>
///     >>> # Filter high-quality alignments
///     >>> for record in parser:
///     ...     if record.evalue < 1e-10 and record.pident > 95.0:
///     ...         print(f"{record.qseqid} -> {record.sseqid}")
///     ...         print(f"  Identity: {record.pident:.1f}%")
///     ...         print(f"  E-value: {record.evalue:.2e}")
///     ...         print(f"  Bit score: {record.bitscore}")
///     >>>
///     >>> # Group by query
///     >>> hits_by_query = {}
///     >>> for record in parser:
///     ...     if record.qseqid not in hits_by_query:
///     ...         hits_by_query[record.qseqid] = []
///     ...     hits_by_query[record.qseqid].append(record)
///     >>>
///     >>> # Find best hit per query
///     >>> for query, hits in hits_by_query.items():
///     ...     best = min(hits, key=lambda h: h.evalue)
///     ...     print(f"{query}: best hit = {best.sseqid} (E={best.evalue:.2e})")
///
/// Note:
///     - Memory footprint remains constant at ~5 MB
///     - Comment lines (starting with #) are automatically skipped
///     - Both outfmt 6 and outfmt 7 formats are supported
#[pyclass(name = "BlastTabularParser", unsendable)]
pub struct PyBlastTabularParser {
    inner: Option<BlastTabularParser<File>>,
}

#[pymethods]
impl PyBlastTabularParser {
    /// Open BLAST tabular file for streaming
    ///
    /// Args:
    ///     path (str): Path to BLAST tabular file
    ///
    /// Returns:
    ///     BlastTabularParser: Streaming iterator with constant memory
    ///
    /// Raises:
    ///     IOError: If file cannot be opened
    ///     ValueError: If file format is invalid
    ///
    /// Example:
    ///     >>> parser = biometal.BlastTabularParser.from_path("alignments.blast")
    ///     >>> for record in parser:
    ///     ...     print(f"{record.qseqid}: {record.pident:.1f}% identity")
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let parser = BlastTabularParser::from_path(&PathBuf::from(path))
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        Ok(PyBlastTabularParser {
            inner: Some(parser),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyBlastRecord> {
        if let Some(ref mut parser) = slf.inner {
            match parser.next() {
                Some(Ok(record)) => Ok(record.into()),
                Some(Err(e)) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())),
                None => Err(PyStopIteration::new_err("no more records")),
            }
        } else {
            Err(PyStopIteration::new_err("iterator exhausted"))
        }
    }

    fn __repr__(&self) -> String {
        "BlastTabularParser()".to_string()
    }
}
