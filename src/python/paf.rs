//! Python wrappers for PAF format (minimap2 alignments)

use pyo3::prelude::*;
use crate::formats::paf::{PafRecord, PafParser};
use crate::formats::TabDelimitedRecord;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::io::{BufReader, Read};
use flate2::read::MultiGzDecoder;

/// Helper function to open a file, detecting and handling gzip compression
fn open_file(path: &Path) -> std::io::Result<Box<dyn Read>> {
    let file = File::open(path)?;

    // Check if file is gzipped by extension
    if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        Ok(Box::new(MultiGzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

/// PAF alignment record
///
/// Pairwise mapping record from minimap2 output.
///
/// Attributes:
///     query_name (str): Query sequence name
///     query_length (int): Query sequence length
///     query_start (int): Query start (0-based)
///     query_end (int): Query end (0-based, exclusive)
///     strand (str): Relative strand ('+' or '-')
///     target_name (str): Target sequence name
///     target_length (int): Target sequence length
///     target_start (int): Target start (0-based)
///     target_end (int): Target end (0-based, exclusive)
///     num_matches (int): Number of matching bases
///     alignment_length (int): Total alignment length (matches + mismatches + gaps)
///     mapq (int): Mapping quality (0-255, 255 = missing)
///
/// Example:
///     >>> line = "read1\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t60"
///     >>> record = PafRecord.from_line(line)
///     >>> print(f"{record.query_name}: {record.identity():.2%} identity")
///     read1: 96.94% identity
#[pyclass(name = "PafRecord")]
#[derive(Clone)]
pub struct PyPafRecord {
    inner: PafRecord,
}

#[pymethods]
impl PyPafRecord {
    #[staticmethod]
    fn from_line(line: &str) -> PyResult<Self> {
        let record = PafRecord::from_line(line)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;
        Ok(record.into())
    }

    // Property getters
    #[getter]
    fn query_name(&self) -> &str {
        &self.inner.query_name
    }

    #[getter]
    fn query_length(&self) -> u64 {
        self.inner.query_length
    }

    #[getter]
    fn query_start(&self) -> u64 {
        self.inner.query_start
    }

    #[getter]
    fn query_end(&self) -> u64 {
        self.inner.query_end
    }

    #[getter]
    fn strand(&self) -> char {
        self.inner.strand
    }

    #[getter]
    fn target_name(&self) -> &str {
        &self.inner.target_name
    }

    #[getter]
    fn target_length(&self) -> u64 {
        self.inner.target_length
    }

    #[getter]
    fn target_start(&self) -> u64 {
        self.inner.target_start
    }

    #[getter]
    fn target_end(&self) -> u64 {
        self.inner.target_end
    }

    #[getter]
    fn num_matches(&self) -> u64 {
        self.inner.num_matches
    }

    #[getter]
    fn alignment_length(&self) -> u64 {
        self.inner.alignment_length
    }

    #[getter]
    fn mapq(&self) -> u8 {
        self.inner.mapq
    }

    fn to_line(&self) -> String {
        self.inner.to_line()
    }

    fn identity(&self) -> f64 {
        self.inner.identity()
    }

    fn query_coverage(&self) -> f64 {
        self.inner.query_coverage()
    }

    fn target_coverage(&self) -> f64 {
        self.inner.target_coverage()
    }

    fn is_high_quality(&self, min_mapq: u8) -> bool {
        self.inner.is_high_quality(min_mapq)
    }

    fn is_forward(&self) -> bool {
        self.inner.is_forward()
    }

    fn query_aligned_length(&self) -> u64 {
        self.inner.query_aligned_length()
    }

    fn target_aligned_length(&self) -> u64 {
        self.inner.target_aligned_length()
    }

    fn __repr__(&self) -> String {
        format!(
            "PafRecord(query='{}', target='{}', identity={:.2}%)",
            self.inner.query_name, self.inner.target_name, self.identity() * 100.0
        )
    }

    fn __str__(&self) -> String {
        format!(
            "{} -> {} ({}): {:.2}% identity, mapq={}",
            self.inner.query_name,
            self.inner.target_name,
            self.inner.strand,
            self.identity() * 100.0,
            self.inner.mapq
        )
    }
}

impl From<PafRecord> for PyPafRecord {
    fn from(record: PafRecord) -> Self {
        PyPafRecord {
            inner: record,
        }
    }
}

/// Stream PAF records with constant memory
///
/// Streaming iterator for minimap2 PAF alignment files.
///
/// Args:
///     path (str): Path to PAF file (.paf or .paf.gz)
///
/// Example:
///     >>> stream = PafStream.from_path("alignments.paf")
///     >>> high_quality = []
///     >>> for record in stream:
///     ...     if record.is_high_quality(20) and record.identity() >= 0.95:
///     ...         high_quality.append(record)
///     >>> print(f"Found {len(high_quality)} high-quality alignments")
#[pyclass(name = "PafStream", unsendable)]
pub struct PyPafStream {
    inner: Option<PafParser<BufReader<Box<dyn Read>>>>,
}

#[pymethods]
impl PyPafStream {
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let reader = open_file(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        let buf_reader = BufReader::new(reader);
        let parser = PafParser::new(buf_reader);

        Ok(PyPafStream {
            inner: Some(parser),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyPafRecord> {
        if let Some(ref mut parser) = slf.inner {
            match parser.next() {
                Some(Ok(record)) => Ok(record.into()),
                Some(Err(e)) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())),
                None => Err(pyo3::exceptions::PyStopIteration::new_err("no more records")),
            }
        } else {
            Err(pyo3::exceptions::PyStopIteration::new_err("stream exhausted"))
        }
    }

    fn __repr__(&self) -> String {
        "PafStream(...)".to_string()
    }
}
