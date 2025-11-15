//! Python wrappers for GTF format (gene annotations, RNA-seq)

use pyo3::prelude::*;
use crate::formats::gtf::{GtfRecord, GtfParser};
use crate::formats::primitives::Strand;
use crate::formats::TabDelimitedRecord;
use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::io::Read;
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

/// GTF feature record
///
/// Represents a genomic feature annotation (gene, transcript, exon, CDS, etc.)
/// Coordinates are 1-based, inclusive [start, end].
///
/// Attributes:
///     seqname (str): Sequence/chromosome name
///     source (str): Annotation source (e.g., "HAVANA", "Ensembl")
///     feature (str): Feature type (e.g., "gene", "transcript", "exon", "CDS")
///     start (int): Start position (1-based, inclusive)
///     end (int): End position (1-based, inclusive)
///     score (float | None): Feature score
///     strand (str): Strand ('+', '-', or '.')
///     frame (int | None): CDS reading frame (0, 1, or 2)
///     attributes (dict[str, str]): Feature attributes (gene_id, transcript_id, etc.)
///
/// Example:
///     >>> line = 'chr1\tHAVANA\texon\t11869\t12227\t.\t+\t.\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328";'
///     >>> record = GtfRecord.from_line(line)
///     >>> print(f"{record.gene_id()}: {record.feature} on {record.seqname}")
///     ENSG00000223972: exon on chr1
///     >>> print(f"Location: {record.start}-{record.end} ({record.strand})")
///     Location: 11869-12227 (+)
#[pyclass(name = "GtfRecord")]
#[derive(Clone)]
pub struct PyGtfRecord {
    inner: GtfRecord,
}

#[pymethods]
impl PyGtfRecord {
    #[staticmethod]
    fn from_line(line: &str) -> PyResult<Self> {
        let record = GtfRecord::from_line(line)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;
        Ok(record.into())
    }

    // Property getters
    #[getter]
    fn seqname(&self) -> &str {
        &self.inner.seqname
    }

    #[getter]
    fn source(&self) -> &str {
        &self.inner.source
    }

    #[getter]
    fn feature(&self) -> &str {
        &self.inner.feature
    }

    #[getter]
    fn start(&self) -> u64 {
        self.inner.start
    }

    #[getter]
    fn end(&self) -> u64 {
        self.inner.end
    }

    #[getter]
    fn score(&self) -> Option<f64> {
        self.inner.score
    }

    #[getter]
    fn strand(&self) -> String {
        match self.inner.strand {
            Strand::Forward => "+".to_string(),
            Strand::Reverse => "-".to_string(),
            Strand::Unknown => ".".to_string(),
        }
    }

    #[getter]
    fn frame(&self) -> Option<u8> {
        self.inner.frame
    }

    #[getter]
    fn attributes(&self) -> HashMap<String, String> {
        self.inner.attributes.clone()
    }

    fn to_line(&self) -> String {
        self.inner.to_line()
    }

    fn length(&self) -> u64 {
        self.inner.length()
    }

    fn gene_id(&self) -> String {
        self.inner.gene_id().to_string()
    }

    fn transcript_id(&self) -> String {
        self.inner.transcript_id().to_string()
    }

    fn gene_name(&self) -> Option<String> {
        self.inner.gene_name().map(|s| s.to_string())
    }

    fn gene_biotype(&self) -> Option<String> {
        self.inner.gene_biotype().map(|s| s.to_string())
    }

    fn transcript_name(&self) -> Option<String> {
        self.inner.attributes.get("transcript_name").cloned()
    }

    fn transcript_biotype(&self) -> Option<String> {
        self.inner.attributes.get("transcript_biotype").cloned()
    }

    /// Convert to 0-based half-open interval [start, end)
    ///
    /// Returns:
    ///     tuple[int, int]: (start, end) in 0-based coordinates
    ///
    /// Example:
    ///     >>> record = GtfRecord.from_line('chr1\t.\texon\t1000\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";')
    ///     >>> start, end = record.to_0based()
    ///     >>> print(f"0-based: [{start}, {end})")
    ///     0-based: [999, 2000)
    fn to_0based(&self) -> (u64, u64) {
        (self.inner.start - 1, self.inner.end)
    }

    fn __repr__(&self) -> String {
        format!(
            "GtfRecord(seqname='{}', feature='{}', start={}, end={})",
            self.inner.seqname, self.inner.feature, self.inner.start, self.inner.end
        )
    }

    fn __str__(&self) -> String {
        let gene_name = self.gene_name().unwrap_or_else(|| self.gene_id());
        format!(
            "{} ({}): {}:{}-{} [{}]",
            gene_name, self.inner.feature, self.inner.seqname, self.inner.start, self.inner.end, self.strand()
        )
    }
}

impl From<GtfRecord> for PyGtfRecord {
    fn from(record: GtfRecord) -> Self {
        PyGtfRecord {
            inner: record,
        }
    }
}

/// Stream GTF records with constant memory
///
/// Streaming iterator that processes GTF annotation files one record at a time.
/// Automatically skips comment lines (# prefix).
///
/// Args:
///     path (str): Path to GTF file (.gtf or .gtf.gz)
///
/// Example:
///     >>> stream = GtfStream.from_path("annotations.gtf")
///     >>> genes = []
///     >>> exons = []
///     >>> for record in stream:
///     ...     if record.feature == "gene":
///     ...         genes.append(record)
///     ...     elif record.feature == "exon":
///     ...         exons.append(record)
///     >>> print(f"Found {len(genes)} genes, {len(exons)} exons")
///
///     # Group features by gene
///     >>> from collections import defaultdict
///     >>> stream = GtfStream.from_path("annotations.gtf")
///     >>> genes = defaultdict(list)
///     >>> for record in stream:
///     ...     gene_id = record.gene_id()
///     ...     if gene_id:
///     ...         genes[gene_id].append(record)
///     >>> for gene_id, features in genes.items():
///     ...     print(f"Gene {gene_id}: {len(features)} features")
#[pyclass(name = "GtfStream", unsendable)]
pub struct PyGtfStream {
    inner: Option<GtfParser<Box<dyn Read>>>,
}

#[pymethods]
impl PyGtfStream {
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let reader = open_file(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        let parser = GtfParser::new(reader);

        Ok(PyGtfStream {
            inner: Some(parser),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyGtfRecord> {
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
        "GtfStream(...)".to_string()
    }
}

// ============================================================================
// Writer
// ============================================================================

use crate::formats::primitives::TabDelimitedWriter;

/// Write GTF records to a file
///
/// Streaming writer for GTF format with automatic compression support.
///
/// Methods:
///     create(path: str) -> GtfWriter: Create writer for a file path
///     stdout() -> GtfWriter: Create writer for stdout
///     write_record(record: GtfRecord): Write a single record
///     records_written() -> int: Get count of written records
///     finish(): Flush and close the writer
///
/// Example:
///     >>> writer = GtfWriter.create("output.gtf.gz")
///     >>> for record in GtfStream.from_path("input.gtf"):
///     ...     if record.feature == "exon":
///     ...         writer.write_record(record)
///     >>> writer.finish()
#[pyclass(name = "GtfWriter", unsendable)]
pub struct PyGtfWriter {
    inner: Option<TabDelimitedWriter<GtfRecord>>,
}

#[pymethods]
impl PyGtfWriter {
    #[staticmethod]
    fn create(path: String) -> PyResult<Self> {
        let writer = TabDelimitedWriter::create(&PathBuf::from(path))
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        Ok(PyGtfWriter { inner: Some(writer) })
    }

    #[staticmethod]
    fn stdout() -> PyResult<Self> {
        let writer = TabDelimitedWriter::stdout()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        Ok(PyGtfWriter { inner: Some(writer) })
    }

    fn write_record(&mut self, record: &PyGtfRecord) -> PyResult<()> {
        if let Some(ref mut writer) = self.inner {
            writer.write_record(&record.inner)
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyIOError, _>("writer already finished"))
        }
    }

    fn records_written(&self) -> PyResult<usize> {
        if let Some(ref writer) = self.inner {
            Ok(writer.records_written())
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyIOError, _>("writer already finished"))
        }
    }

    fn finish(&mut self) -> PyResult<()> {
        if let Some(writer) = self.inner.take() {
            writer.finish()
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyIOError, _>("writer already finished"))
        }
    }

    fn __repr__(&self) -> String {
        format!("GtfWriter(records_written={})",
            self.inner.as_ref().map(|w| w.records_written()).unwrap_or(0))
    }
}
