//! Python wrappers for BED format records

use pyo3::prelude::*;
use crate::formats::bed::{Bed3Record, Bed6Record, Bed12Record, NarrowPeakRecord};
use crate::formats::primitives::Strand;
use crate::formats::{TabDelimitedRecord, TabDelimitedParser};
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

/// BED3 record (chrom, start, end)
///
/// Minimal BED format with genomic coordinates.
/// Coordinates are 0-based, half-open intervals [start, end).
///
/// Attributes:
///     chrom (str): Chromosome name
///     start (int): Start position (0-based, inclusive)
///     end (int): End position (0-based, exclusive)
///
/// Example:
///     >>> record = Bed3Record.from_line("chr1\t1000\t2000")
///     >>> print(f"{record.chrom}:{record.start}-{record.end}")
///     chr1:1000-2000
#[pyclass(name = "Bed3Record")]
#[derive(Clone)]
pub struct PyBed3Record {
    #[pyo3(get)]
    pub chrom: String,
    #[pyo3(get)]
    pub start: u64,
    #[pyo3(get)]
    pub end: u64,
}

#[pymethods]
impl PyBed3Record {
    #[staticmethod]
    fn from_line(line: &str) -> PyResult<Self> {
        let record = Bed3Record::from_line(line)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;
        Ok(record.into())
    }

    fn to_line(&self) -> String {
        Bed3Record {
            interval: crate::formats::primitives::GenomicInterval {
                chrom: self.chrom.clone(),
                start: self.start,
                end: self.end,
            },
        }.to_line()
    }

    fn length(&self) -> u64 {
        self.end - self.start
    }

    fn __repr__(&self) -> String {
        format!("Bed3Record(chrom='{}', start={}, end={})", self.chrom, self.start, self.end)
    }

    fn __str__(&self) -> String {
        format!("{}:{}-{}", self.chrom, self.start, self.end)
    }
}

impl From<Bed3Record> for PyBed3Record {
    fn from(record: Bed3Record) -> Self {
        PyBed3Record {
            chrom: record.interval.chrom,
            start: record.interval.start,
            end: record.interval.end,
        }
    }
}

/// BED6 record (chrom, start, end, name, score, strand)
///
/// Extended BED format with name, score, and strand information.
///
/// Attributes:
///     chrom (str): Chromosome name
///     start (int): Start position (0-based, inclusive)
///     end (int): End position (0-based, exclusive)
///     name (str | None): Feature name
///     score (int | None): Score (0-1000)
///     strand (str | None): Strand ('+', '-', or '.')
///
/// Example:
///     >>> record = Bed6Record.from_line("chr1\t1000\t2000\tgene1\t500\t+")
///     >>> print(f"{record.name}: {record.chrom}:{record.start}-{record.end} ({record.strand})")
///     gene1: chr1:1000-2000 (+)
#[pyclass(name = "Bed6Record")]
#[derive(Clone)]
pub struct PyBed6Record {
    #[pyo3(get)]
    pub chrom: String,
    #[pyo3(get)]
    pub start: u64,
    #[pyo3(get)]
    pub end: u64,
    #[pyo3(get)]
    pub name: Option<String>,
    #[pyo3(get)]
    pub score: Option<u32>,
    #[pyo3(get)]
    pub strand: Option<String>,
}

#[pymethods]
impl PyBed6Record {
    #[staticmethod]
    fn from_line(line: &str) -> PyResult<Self> {
        let record = Bed6Record::from_line(line)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;
        Ok(record.into())
    }

    fn to_line(&self) -> String {
        let strand = self.strand.as_ref().and_then(|s| match s.as_str() {
            "+" => Some(Strand::Forward),
            "-" => Some(Strand::Reverse),
            "." => Some(Strand::Unknown),
            _ => None,
        });

        Bed6Record {
            bed3: Bed3Record {
                interval: crate::formats::primitives::GenomicInterval {
                    chrom: self.chrom.clone(),
                    start: self.start,
                    end: self.end,
                },
            },
            name: self.name.clone(),
            score: self.score,
            strand,
        }.to_line()
    }

    fn length(&self) -> u64 {
        self.end - self.start
    }

    fn __repr__(&self) -> String {
        format!(
            "Bed6Record(chrom='{}', start={}, end={}, name={:?}, score={:?}, strand={:?})",
            self.chrom, self.start, self.end, self.name, self.score, self.strand
        )
    }

    fn __str__(&self) -> String {
        let mut s = format!("{}:{}-{}", self.chrom, self.start, self.end);
        if let Some(name) = &self.name {
            s.push_str(&format!(" ({})", name));
        }
        if let Some(strand) = &self.strand {
            s.push_str(&format!(" [{}]", strand));
        }
        s
    }
}

impl From<Bed6Record> for PyBed6Record {
    fn from(record: Bed6Record) -> Self {
        PyBed6Record {
            chrom: record.bed3.interval.chrom,
            start: record.bed3.interval.start,
            end: record.bed3.interval.end,
            name: record.name,
            score: record.score,
            strand: record.strand.map(|s| s.to_string()),
        }
    }
}

/// BED12 record (full BED format with blocks)
///
/// Complete BED format with thick start/end, RGB color, and block structure.
///
/// Attributes:
///     chrom (str): Chromosome name
///     start (int): Start position (0-based, inclusive)
///     end (int): End position (0-based, exclusive)
///     name (str | None): Feature name
///     score (int | None): Score (0-1000)
///     strand (str | None): Strand ('+', '-', or '.')
///     thick_start (int | None): Thick region start
///     thick_end (int | None): Thick region end
///     item_rgb (str | None): RGB color string
///     block_count (int | None): Number of blocks
///     block_sizes (list[int] | None): Block sizes
///     block_starts (list[int] | None): Block starts (relative to start)
///
/// Example:
///     >>> line = "chr1\t1000\t5000\tgene1\t500\t+\t1200\t4800\t255,0,0\t3\t500,500,500\t0,2000,4000"
///     >>> record = Bed12Record.from_line(line)
///     >>> print(f"{record.name}: {record.block_count} exons")
///     gene1: 3 exons
#[pyclass(name = "Bed12Record")]
#[derive(Clone)]
pub struct PyBed12Record {
    #[pyo3(get)]
    pub chrom: String,
    #[pyo3(get)]
    pub start: u64,
    #[pyo3(get)]
    pub end: u64,
    #[pyo3(get)]
    pub name: Option<String>,
    #[pyo3(get)]
    pub score: Option<u32>,
    #[pyo3(get)]
    pub strand: Option<String>,
    #[pyo3(get)]
    pub thick_start: Option<u64>,
    #[pyo3(get)]
    pub thick_end: Option<u64>,
    #[pyo3(get)]
    pub item_rgb: Option<String>,
    #[pyo3(get)]
    pub block_count: Option<u32>,
    #[pyo3(get)]
    pub block_sizes: Option<Vec<u32>>,
    #[pyo3(get)]
    pub block_starts: Option<Vec<u32>>,
}

#[pymethods]
impl PyBed12Record {
    #[staticmethod]
    fn from_line(line: &str) -> PyResult<Self> {
        let record = Bed12Record::from_line(line)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;
        Ok(record.into())
    }

    fn to_line(&self) -> String {
        let strand = self.strand.as_ref().and_then(|s| match s.as_str() {
            "+" => Some(Strand::Forward),
            "-" => Some(Strand::Reverse),
            "." => Some(Strand::Unknown),
            _ => None,
        });

        Bed12Record {
            bed6: Bed6Record {
                bed3: Bed3Record {
                    interval: crate::formats::primitives::GenomicInterval {
                        chrom: self.chrom.clone(),
                        start: self.start,
                        end: self.end,
                    },
                },
                name: self.name.clone(),
                score: self.score,
                strand,
            },
            thick_start: self.thick_start,
            thick_end: self.thick_end,
            item_rgb: self.item_rgb.clone(),
            block_count: self.block_count,
            block_sizes: self.block_sizes.clone(),
            block_starts: self.block_starts.clone(),
        }.to_line()
    }

    fn length(&self) -> u64 {
        self.end - self.start
    }

    fn __repr__(&self) -> String {
        format!(
            "Bed12Record(chrom='{}', start={}, end={}, name={:?}, blocks={:?})",
            self.chrom, self.start, self.end, self.name, self.block_count
        )
    }

    fn __str__(&self) -> String {
        let mut s = format!("{}:{}-{}", self.chrom, self.start, self.end);
        if let Some(name) = &self.name {
            s.push_str(&format!(" ({})", name));
        }
        if let Some(blocks) = self.block_count {
            s.push_str(&format!(" [{}b]", blocks));
        }
        s
    }
}

impl From<Bed12Record> for PyBed12Record {
    fn from(record: Bed12Record) -> Self {
        PyBed12Record {
            chrom: record.bed6.bed3.interval.chrom,
            start: record.bed6.bed3.interval.start,
            end: record.bed6.bed3.interval.end,
            name: record.bed6.name,
            score: record.bed6.score,
            strand: record.bed6.strand.map(|s| s.to_string()),
            thick_start: record.thick_start,
            thick_end: record.thick_end,
            item_rgb: record.item_rgb,
            block_count: record.block_count,
            block_sizes: record.block_sizes,
            block_starts: record.block_starts,
        }
    }
}

/// Stream BED3 records with constant memory
///
/// Streaming iterator that processes BED3 files one record at a time.
///
/// Args:
///     path (str): Path to BED file
///
/// Example:
///     >>> stream = Bed3Stream.from_path("peaks.bed")
///     >>> for record in stream:
///     ...     print(f"{record.chrom}:{record.start}-{record.end}")
#[pyclass(name = "Bed3Stream", unsendable)]
pub struct PyBed3Stream {
    inner: Option<TabDelimitedParser<Box<dyn Read>, Bed3Record>>,
}

#[pymethods]
impl PyBed3Stream {
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let reader = open_file(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        let parser = TabDelimitedParser::<Box<dyn Read>, Bed3Record>::new(reader);

        Ok(PyBed3Stream {
            inner: Some(parser),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyBed3Record> {
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
        "Bed3Stream(...)".to_string()
    }
}

/// Stream BED6 records with constant memory
///
/// Streaming iterator that processes BED6 files one record at a time.
///
/// Args:
///     path (str): Path to BED file
///
/// Example:
///     >>> stream = Bed6Stream.from_path("genes.bed")
///     >>> for record in stream:
///     ...     print(f"{record.name}: {record.chrom}:{record.start}-{record.end}")
#[pyclass(name = "Bed6Stream", unsendable)]
pub struct PyBed6Stream {
    inner: Option<TabDelimitedParser<Box<dyn Read>, Bed6Record>>,
}

#[pymethods]
impl PyBed6Stream {
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let reader = open_file(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        let parser = TabDelimitedParser::<Box<dyn Read>, Bed6Record>::new(reader);

        Ok(PyBed6Stream {
            inner: Some(parser),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyBed6Record> {
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
        "Bed6Stream(...)".to_string()
    }
}

/// Stream BED12 records with constant memory
///
/// Streaming iterator that processes BED12 files one record at a time.
///
/// Args:
///     path (str): Path to BED file
///
/// Example:
///     >>> stream = Bed12Stream.from_path("transcripts.bed")
///     >>> for record in stream:
///     ...     print(f"{record.name}: {record.block_count} exons")
#[pyclass(name = "Bed12Stream", unsendable)]
pub struct PyBed12Stream {
    inner: Option<TabDelimitedParser<Box<dyn Read>, Bed12Record>>,
}

#[pymethods]
impl PyBed12Stream {
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let reader = open_file(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        let parser = TabDelimitedParser::<Box<dyn Read>, Bed12Record>::new(reader);

        Ok(PyBed12Stream {
            inner: Some(parser),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyBed12Record> {
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
        "Bed12Stream(...)".to_string()
    }
}

/// narrowPeak record (ENCODE ChIP-seq peaks)
///
/// BED6+4 format for called peaks of signal enrichment from ChIP-seq data.
/// Extends BED6 with peak-calling statistics.
///
/// Attributes:
///     chrom (str): Chromosome name
///     start (int): Start position (0-based, inclusive)
///     end (int): End position (0-based, exclusive)
///     name (str | None): Peak identifier
///     score (int | None): Score (0-1000)
///     strand (str | None): Strand (+/-)
///     signal_value (float): Overall enrichment measurement
///     p_value (float): Statistical significance as -log10 (-1 if unavailable)
///     q_value (float): FDR significance as -log10 (-1 if unavailable)
///     peak (int | None): Offset from start to peak summit (None if -1)
///
/// Example:
///     >>> record = NarrowPeakRecord.from_line("chr1\t1000\t2000\tpeak1\t100\t.\t12.5\t8.3\t5.2\t450")
///     >>> print(f"{record.name}: signal={record.signal_value}, summit={record.peak_position()}")
///     peak1: signal=12.5, summit=1450
#[pyclass(name = "NarrowPeakRecord")]
#[derive(Clone)]
pub struct PyNarrowPeakRecord {
    inner: NarrowPeakRecord,
}

#[pymethods]
impl PyNarrowPeakRecord {
    #[staticmethod]
    fn from_line(line: &str) -> PyResult<Self> {
        let record = NarrowPeakRecord::from_line(line)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;
        Ok(record.into())
    }

    // Property getters
    #[getter]
    fn chrom(&self) -> &str {
        &self.inner.bed6.bed3.interval.chrom
    }

    #[getter]
    fn start(&self) -> u64 {
        self.inner.bed6.bed3.interval.start
    }

    #[getter]
    fn end(&self) -> u64 {
        self.inner.bed6.bed3.interval.end
    }

    #[getter]
    fn name(&self) -> Option<String> {
        self.inner.bed6.name.clone()
    }

    #[getter]
    fn score(&self) -> Option<u32> {
        self.inner.bed6.score
    }

    #[getter]
    fn strand(&self) -> Option<String> {
        self.inner.bed6.strand.as_ref().map(|s| s.to_string())
    }

    #[getter]
    fn signal_value(&self) -> f64 {
        self.inner.signal_value
    }

    #[getter]
    fn p_value(&self) -> f64 {
        self.inner.p_value
    }

    #[getter]
    fn q_value(&self) -> f64 {
        self.inner.q_value
    }

    #[getter]
    fn peak(&self) -> Option<i32> {
        self.inner.peak
    }

    fn to_line(&self) -> String {
        self.inner.to_line()
    }

    fn length(&self) -> u64 {
        self.inner.bed6.bed3.interval.end - self.inner.bed6.bed3.interval.start
    }

    fn peak_position(&self) -> Option<u64> {
        self.inner.peak_position()
    }

    fn is_significant(&self, p_threshold: f64, q_threshold: f64) -> bool {
        self.inner.is_significant(p_threshold, q_threshold)
    }

    fn __repr__(&self) -> String {
        format!(
            "NarrowPeakRecord(chrom='{}', start={}, end={}, signal={})",
            self.inner.bed6.bed3.interval.chrom,
            self.inner.bed6.bed3.interval.start,
            self.inner.bed6.bed3.interval.end,
            self.inner.signal_value
        )
    }

    fn __str__(&self) -> String {
        let name = self.inner.bed6.name.as_ref().map(|s| s.as_str()).unwrap_or("?");
        format!(
            "{} ({}:{}-{}): signal={:.2}, p={:.2}, q={:.2}",
            name,
            self.inner.bed6.bed3.interval.chrom,
            self.inner.bed6.bed3.interval.start,
            self.inner.bed6.bed3.interval.end,
            self.inner.signal_value,
            self.inner.p_value,
            self.inner.q_value
        )
    }
}

impl From<NarrowPeakRecord> for PyNarrowPeakRecord {
    fn from(record: NarrowPeakRecord) -> Self {
        PyNarrowPeakRecord {
            inner: record,
        }
    }
}

/// Stream narrowPeak records with constant memory
///
/// Streaming iterator for ENCODE narrowPeak (ChIP-seq peak) files.
///
/// Args:
///     path (str): Path to narrowPeak file (.narrowPeak or .narrowPeak.gz)
///
/// Example:
///     >>> stream = NarrowPeakStream.from_path("peaks.narrowPeak")
///     >>> significant_peaks = []
///     >>> for record in stream:
///     ...     if record.is_significant(2.0, 1.3):  # p<0.01, q<0.05
///     ...         significant_peaks.append(record)
///     >>> print(f"Found {len(significant_peaks)} significant peaks")
#[pyclass(name = "NarrowPeakStream", unsendable)]
pub struct PyNarrowPeakStream {
    inner: Option<TabDelimitedParser<Box<dyn Read>, NarrowPeakRecord>>,
}

#[pymethods]
impl PyNarrowPeakStream {
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let reader = open_file(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        let parser = TabDelimitedParser::new(reader);

        Ok(PyNarrowPeakStream {
            inner: Some(parser),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyNarrowPeakRecord> {
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
        "NarrowPeakStream(...)".to_string()
    }
}
