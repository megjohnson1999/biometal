//! Python wrappers for VCF format (variant calling)

use pyo3::prelude::*;
use crate::formats::vcf::{VcfHeader, VcfRecord, VcfParser};
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

/// VCF file header
///
/// Contains metadata from VCF header lines (##fileformat, ##INFO, etc.)
///
/// Attributes:
///     fileformat (str): VCF version (e.g., "VCFv4.2")
///     info_fields (dict[str, str]): INFO field definitions
///     format_fields (dict[str, str]): FORMAT field definitions
///     filters (dict[str, str]): FILTER definitions
///     contigs (dict[str, int | None]): Contig names and lengths
///     samples (list[str]): Sample IDs
///
/// Example:
///     >>> parser = VcfStream.from_path("variants.vcf")
///     >>> header = parser.header()
///     >>> print(f"VCF version: {header.fileformat}")
///     >>> print(f"Samples: {', '.join(header.samples)}")
#[pyclass(name = "VcfHeader")]
#[derive(Clone)]
pub struct PyVcfHeader {
    #[pyo3(get, set)]
    pub fileformat: String,
    #[pyo3(get, set)]
    pub info_fields: HashMap<String, String>,
    #[pyo3(get, set)]
    pub format_fields: HashMap<String, String>,
    #[pyo3(get, set)]
    pub filters: HashMap<String, String>,
    #[pyo3(get, set)]
    pub contigs: HashMap<String, Option<u64>>,
    #[pyo3(get, set)]
    pub samples: Vec<String>,
}

#[pymethods]
impl PyVcfHeader {
    #[new]
    fn new(fileformat: String) -> Self {
        PyVcfHeader {
            fileformat,
            info_fields: HashMap::new(),
            format_fields: HashMap::new(),
            filters: HashMap::new(),
            contigs: HashMap::new(),
            samples: Vec::new(),
        }
    }

    /// Add an INFO field definition
    fn add_info(&mut self, id: String, description: String) {
        self.info_fields.insert(id, description);
    }

    /// Add a FORMAT field definition
    fn add_format(&mut self, id: String, description: String) {
        self.format_fields.insert(id, description);
    }

    /// Add a FILTER definition
    fn add_filter(&mut self, id: String, description: String) {
        self.filters.insert(id, description);
    }

    /// Add a contig
    fn add_contig(&mut self, id: String, length: Option<u64>) {
        self.contigs.insert(id, length);
    }

    fn __repr__(&self) -> String {
        format!(
            "VcfHeader(fileformat='{}', samples={})",
            self.fileformat,
            self.samples.len()
        )
    }

    fn __str__(&self) -> String {
        format!(
            "VCF {}: {} samples, {} contigs",
            self.fileformat,
            self.samples.len(),
            self.contigs.len()
        )
    }
}

impl From<VcfHeader> for PyVcfHeader {
    fn from(header: VcfHeader) -> Self {
        PyVcfHeader {
            fileformat: header.fileformat,
            info_fields: header.info_fields,
            format_fields: header.format_fields,
            filters: header.filters,
            contigs: header.contigs,
            samples: header.samples,
        }
    }
}

/// VCF variant record
///
/// Represents a single genetic variant.
///
/// Attributes:
///     chrom (str): Chromosome name
///     pos (int): Position (1-based)
///     id (str | None): Variant ID (e.g., rs12345)
///     reference (str): Reference allele
///     alternate (list[str]): Alternate alleles
///     quality (float | None): Variant quality score
///     filter (str | None): Filter status (e.g., "PASS")
///     info (dict[str, str]): INFO field key-value pairs
///     format (str | None): FORMAT field specification
///     samples (list[str]): Sample genotype data
///
/// Example:
///     >>> record = VcfRecord.from_line("chr1\t12345\trs123\tA\tT\t30\tPASS\tDP=100")
///     >>> print(f"{record.chrom}:{record.pos} {record.reference}>{record.alternate[0]}")
///     chr1:12345 A>T
///     >>> if record.filter == "PASS":
///     ...     depth = record.info.get("DP")
///     ...     print(f"Depth: {depth}")
///     Depth: 100
#[pyclass(name = "VcfRecord")]
#[derive(Clone)]
pub struct PyVcfRecord {
    #[pyo3(get)]
    pub chrom: String,
    #[pyo3(get)]
    pub pos: u64,
    #[pyo3(get)]
    pub id: Option<String>,
    #[pyo3(get)]
    pub reference: String,
    #[pyo3(get)]
    pub alternate: Vec<String>,
    #[pyo3(get)]
    pub quality: Option<f64>,
    #[pyo3(get)]
    pub filter: Option<String>,
    #[pyo3(get)]
    pub info: HashMap<String, String>,
    #[pyo3(get)]
    pub format: Option<String>,
    #[pyo3(get)]
    pub samples: Vec<String>,
}

#[pymethods]
impl PyVcfRecord {
    #[staticmethod]
    fn from_line(line: &str) -> PyResult<Self> {
        let record = VcfRecord::from_line(line)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;
        Ok(record.into())
    }

    fn to_line(&self) -> String {
        VcfRecord {
            chrom: self.chrom.clone(),
            pos: self.pos,
            id: self.id.clone(),
            reference: self.reference.clone(),
            alternate: self.alternate.clone(),
            quality: self.quality,
            filter: self.filter.clone(),
            info: self.info.clone(),
            format: self.format.clone(),
            samples: self.samples.clone(),
        }
        .to_line()
    }

    fn is_snp(&self) -> bool {
        self.reference.len() == 1 && self.alternate.iter().all(|a| a.len() == 1)
    }

    fn is_insertion(&self) -> bool {
        self.alternate.iter().any(|a| a.len() > self.reference.len())
    }

    fn is_deletion(&self) -> bool {
        self.alternate.iter().any(|a| a.len() < self.reference.len())
    }

    fn is_indel(&self) -> bool {
        self.is_insertion() || self.is_deletion()
    }

    fn __repr__(&self) -> String {
        format!(
            "VcfRecord(chrom='{}', pos={}, ref='{}', alt={:?})",
            self.chrom, self.pos, self.reference, self.alternate
        )
    }

    fn __str__(&self) -> String {
        format!(
            "{}:{} {}->{}",
            self.chrom,
            self.pos,
            self.reference,
            self.alternate.join(",")
        )
    }
}

impl From<VcfRecord> for PyVcfRecord {
    fn from(record: VcfRecord) -> Self {
        PyVcfRecord {
            chrom: record.chrom,
            pos: record.pos,
            id: record.id,
            reference: record.reference,
            alternate: record.alternate,
            quality: record.quality,
            filter: record.filter,
            info: record.info,
            format: record.format,
            samples: record.samples,
        }
    }
}

/// Stream VCF records with constant memory
///
/// Streaming iterator that processes VCF variant files one record at a time.
/// Header must be parsed before iterating records.
///
/// Args:
///     path (str): Path to VCF file
///
/// Example:
///     >>> stream = VcfStream.from_path("variants.vcf")
///     >>> header = stream.header()
///     >>> print(f"VCF version: {header.fileformat}")
///     >>>
///     >>> snps = 0
///     >>> indels = 0
///     >>> for record in stream:
///     ...     if record.is_snp():
///     ...         snps += 1
///     ...     elif record.is_indel():
///     ...         indels += 1
///     >>> print(f"SNPs: {snps}, Indels: {indels}")
#[pyclass(name = "VcfStream", unsendable)]
pub struct PyVcfStream {
    inner: Option<VcfParser<Box<dyn Read>>>,
    header: Option<PyVcfHeader>,
}

#[pymethods]
impl PyVcfStream {
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let reader = open_file(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        let mut parser = VcfParser::new(reader);

        // Parse header immediately
        let header = parser
            .parse_header()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;

        Ok(PyVcfStream {
            inner: Some(parser),
            header: Some(header.into()),
        })
    }

    fn header(&self) -> PyResult<PyVcfHeader> {
        self.header
            .clone()
            .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyValueError, _>("Header not parsed"))
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyVcfRecord> {
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
        "VcfStream(...)".to_string()
    }
}

// ============================================================================
// Writer
// ============================================================================

use crate::formats::vcf::VcfWriter;
use crate::io::DataSink;

impl From<PyVcfHeader> for VcfHeader {
    fn from(py_header: PyVcfHeader) -> Self {
        VcfHeader {
            fileformat: py_header.fileformat,
            info_fields: py_header.info_fields,
            format_fields: py_header.format_fields,
            filters: py_header.filters,
            contigs: py_header.contigs,
            metadata: Vec::new(), // Python interface doesn't expose raw metadata yet
            samples: py_header.samples,
        }
    }
}

/// Write VCF records to a file
///
/// Streaming writer for VCF (Variant Call Format) with header support.
/// Automatically handles compression based on file extension (.gz, .bgz).
///
/// Methods:
///     create(path: str, header: VcfHeader) -> VcfWriter: Create writer for a file path
///     stdout(header: VcfHeader) -> VcfWriter: Create writer for stdout
///     write_record(record: VcfRecord): Write a single variant record
///     records_written() -> int: Get count of written records
///     finish(): Flush and close the writer
///
/// Example:
///     >>> # Create header
///     >>> header = VcfHeader("VCFv4.2")
///     >>> header.info_fields["DP"] = "Total Depth"
///     >>> header.samples = ["sample1", "sample2"]
///     >>>
///     >>> # Write VCF
///     >>> writer = VcfWriter.create("output.vcf.gz", header)
///     >>> for record in VcfStream.from_path("input.vcf"):
///     ...     if record.filter == "PASS":
///     ...         writer.write_record(record)
///     >>> writer.finish()
#[pyclass(name = "VcfWriter", unsendable)]
pub struct PyVcfWriter {
    inner: Option<VcfWriter>,
}

#[pymethods]
impl PyVcfWriter {
    #[staticmethod]
    fn create(path: String, header: PyVcfHeader) -> PyResult<Self> {
        let rust_header: VcfHeader = header.into();
        let writer = VcfWriter::create(&PathBuf::from(path), rust_header)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        Ok(PyVcfWriter { inner: Some(writer) })
    }

    #[staticmethod]
    fn stdout(header: PyVcfHeader) -> PyResult<Self> {
        let rust_header: VcfHeader = header.into();
        let writer = VcfWriter::stdout(rust_header)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        Ok(PyVcfWriter { inner: Some(writer) })
    }

    fn write_record(&mut self, record: &PyVcfRecord) -> PyResult<()> {
        if let Some(ref mut writer) = self.inner {
            let rust_record = VcfRecord {
                chrom: record.chrom.clone(),
                pos: record.pos,
                id: record.id.clone(),
                reference: record.reference.clone(),
                alternate: record.alternate.clone(),
                quality: record.quality,
                filter: record.filter.clone(),
                info: record.info.clone(),
                format: record.format.clone(),
                samples: record.samples.clone(),
            };

            writer.write_record(&rust_record)
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
        format!("VcfWriter(records_written={})",
            self.inner.as_ref().map(|w| w.records_written()).unwrap_or(0))
    }
}
