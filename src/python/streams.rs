//! Python wrappers for streaming classes

use pyo3::prelude::*;
use pyo3::exceptions::PyStopIteration;
use std::path::PathBuf;
use crate::io::{FastqStream, FastaStream, FastqWriter, FastaWriter, CompressedReader};
use crate::python::records::{PyFastqRecord, PyFastaRecord};

/// Stream FASTQ records with constant memory
///
/// Streaming iterator that processes FASTQ files one record at a time,
/// maintaining constant ~5 MB memory regardless of file size.
///
/// Args:
///     path (str): Path to FASTQ file (.fq, .fastq, .fq.gz, .fastq.gz)
///
/// Example:
///     >>> stream = biometal.FastqStream.from_path("data.fq.gz")
///     >>> for record in stream:
///     ...     gc = biometal.gc_content(record.sequence)
///     ...     print(f"{record.id}: {gc:.2%}")
///
/// Note:
///     Memory footprint remains constant at ~5 MB even for TB-scale files.
///     This enables analysis on consumer hardware without downloading.
#[pyclass(name = "FastqStream", unsendable)]
pub struct PyFastqStream {
    inner: Option<FastqStream<CompressedReader>>,
}

#[pymethods]
impl PyFastqStream {
    /// Create FASTQ stream from file path
    ///
    /// Args:
///         path (str): Path to FASTQ file
    ///
    /// Returns:
    ///     FastqStream: Streaming iterator
    ///
    /// Raises:
    ///     IOError: If file cannot be opened
    ///     ValueError: If file format is invalid
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let stream = FastqStream::from_path(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        Ok(PyFastqStream {
            inner: Some(stream),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyFastqRecord> {
        if let Some(ref mut stream) = slf.inner {
            match stream.next() {
                Some(Ok(record)) => Ok(record.into()),
                Some(Err(e)) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())),
                None => Err(PyStopIteration::new_err("no more records")),
            }
        } else {
            Err(PyStopIteration::new_err("stream exhausted"))
        }
    }

    fn __repr__(&self) -> String {
        "FastqStream(...)".to_string()
    }
}

/// Stream FASTA records with constant memory
///
/// Streaming iterator that processes FASTA files one record at a time,
/// maintaining constant ~5 MB memory regardless of file size.
///
/// Args:
///     path (str): Path to FASTA file (.fa, .fasta, .fa.gz, .fasta.gz)
///
/// Example:
///     >>> stream = biometal.FastaStream.from_path("genome.fa.gz")
///     >>> for record in stream:
///     ...     length = len(record.sequence)
///     ...     print(f"{record.id}: {length} bp")
///
/// Note:
///     Memory footprint remains constant at ~5 MB even for entire genomes.
#[pyclass(name = "FastaStream", unsendable)]
pub struct PyFastaStream {
    inner: Option<FastaStream<CompressedReader>>,
}

#[pymethods]
impl PyFastaStream {
    /// Create FASTA stream from file path
    ///
    /// Args:
    ///     path (str): Path to FASTA file
    ///
    /// Returns:
    ///     FastaStream: Streaming iterator
    ///
    /// Raises:
    ///     IOError: If file cannot be opened
    ///     ValueError: If file format is invalid
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let stream = FastaStream::from_path(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        Ok(PyFastaStream {
            inner: Some(stream),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyFastaRecord> {
        if let Some(ref mut stream) = slf.inner {
            match stream.next() {
                Some(Ok(record)) => Ok(record.into()),
                Some(Err(e)) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())),
                None => Err(PyStopIteration::new_err("no more records")),
            }
        } else {
            Err(PyStopIteration::new_err("stream exhausted"))
        }
    }

    fn __repr__(&self) -> String {
        "FastaStream(...)".to_string()
    }
}

/// Write FASTQ records with compression support
///
/// Writer that supports automatic compression based on file extension.
/// Supports .fq, .fastq, .fq.gz, .fastq.gz, .fq.bgz, .fastq.bgz
///
/// Args:
///     path (str): Path to output FASTQ file
///
/// Example:
///     >>> writer = biometal.FastqWriter.create("output.fq.gz")
///     >>> record = biometal.FastqRecord(...)
///     >>> writer.write_record(record)
///     >>> writer.finish()  # IMPORTANT: Flush and close
///
/// Note:
///     The finish() method MUST be called to ensure all data is written.
///     Alternatively, use context manager (not yet implemented).
#[pyclass(name = "FastqWriter", unsendable)]
pub struct PyFastqWriter {
    inner: Option<FastqWriter>,
}

#[pymethods]
impl PyFastqWriter {
    /// Create FASTQ writer from file path
    ///
    /// Args:
    ///     path (str): Path to output file
    ///
    /// Returns:
    ///     FastqWriter: Writer instance
    ///
    /// Raises:
    ///     IOError: If file cannot be created
    ///
    /// Example:
    ///     >>> writer = biometal.FastqWriter.create("output.fq.gz")
    #[staticmethod]
    fn create(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let writer = FastqWriter::create(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        Ok(PyFastqWriter {
            inner: Some(writer),
        })
    }

    /// Create FASTQ writer to stdout
    ///
    /// Returns:
    ///     FastqWriter: Writer instance
    ///
    /// Raises:
    ///     IOError: If stdout cannot be accessed
    ///
    /// Example:
    ///     >>> writer = biometal.FastqWriter.stdout()
    #[staticmethod]
    fn stdout() -> PyResult<Self> {
        let writer = FastqWriter::stdout()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        Ok(PyFastqWriter {
            inner: Some(writer),
        })
    }

    /// Write a single FASTQ record
    ///
    /// Args:
    ///     record (FastqRecord): Record to write
    ///
    /// Raises:
    ///     ValueError: If record is invalid
    ///     IOError: If write fails
    ///
    /// Example:
    ///     >>> writer.write_record(record)
    fn write_record(&mut self, record: &PyFastqRecord) -> PyResult<()> {
        if let Some(ref mut writer) = self.inner {
            let rust_record = record.to_fastq_record();
            writer
                .write_record(&rust_record)
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                "writer already finished",
            ))
        }
    }

    /// Get the number of records written so far
    ///
    /// Returns:
    ///     int: Number of records written
    ///
    /// Example:
    ///     >>> count = writer.records_written()
    fn records_written(&self) -> PyResult<usize> {
        if let Some(ref writer) = self.inner {
            Ok(writer.records_written())
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                "writer already finished",
            ))
        }
    }

    /// Flush buffered data to disk
    ///
    /// Raises:
    ///     IOError: If flush fails
    ///
    /// Example:
    ///     >>> writer.flush()
    fn flush(&mut self) -> PyResult<()> {
        if let Some(ref mut writer) = self.inner {
            writer
                .flush()
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                "writer already finished",
            ))
        }
    }

    /// Finish writing and flush all data
    ///
    /// This method MUST be called to ensure all data is written to disk.
    ///
    /// Raises:
    ///     IOError: If finish fails
    ///
    /// Example:
    ///     >>> writer.finish()
    fn finish(&mut self) -> PyResult<()> {
        if let Some(writer) = self.inner.take() {
            writer
                .finish()
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                "writer already finished",
            ))
        }
    }

    fn __repr__(&self) -> String {
        "FastqWriter(...)".to_string()
    }
}

/// Write FASTA records with compression support
///
/// Writer that supports automatic compression based on file extension.
/// Supports .fa, .fasta, .fa.gz, .fasta.gz, .fa.bgz, .fasta.bgz
///
/// Args:
///     path (str): Path to output FASTA file
///
/// Example:
///     >>> writer = biometal.FastaWriter.create("output.fa.gz")
///     >>> record = biometal.FastaRecord(...)
///     >>> writer.write_record(record)
///     >>> writer.finish()  # IMPORTANT: Flush and close
///
/// Note:
///     The finish() method MUST be called to ensure all data is written.
///     By default, sequences are wrapped at 80 characters per line.
#[pyclass(name = "FastaWriter", unsendable)]
pub struct PyFastaWriter {
    inner: Option<FastaWriter>,
}

#[pymethods]
impl PyFastaWriter {
    /// Create FASTA writer from file path
    ///
    /// Args:
    ///     path (str): Path to output file
    ///
    /// Returns:
    ///     FastaWriter: Writer instance
    ///
    /// Raises:
    ///     IOError: If file cannot be created
    ///
    /// Example:
    ///     >>> writer = biometal.FastaWriter.create("output.fa.gz")
    #[staticmethod]
    fn create(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let writer = FastaWriter::create(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        Ok(PyFastaWriter {
            inner: Some(writer),
        })
    }

    /// Create FASTA writer to stdout
    ///
    /// Returns:
    ///     FastaWriter: Writer instance
    ///
    /// Raises:
    ///     IOError: If stdout cannot be accessed
    ///
    /// Example:
    ///     >>> writer = biometal.FastaWriter.stdout()
    #[staticmethod]
    fn stdout() -> PyResult<Self> {
        let writer = FastaWriter::stdout()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        Ok(PyFastaWriter {
            inner: Some(writer),
        })
    }

    /// Set the line width for sequence wrapping
    ///
    /// Args:
    ///     width (int): Number of characters per line (use sys.maxsize to disable wrapping)
    ///
    /// Returns:
    ///     FastaWriter: Self (for method chaining)
    ///
    /// Example:
    ///     >>> writer = biometal.FastaWriter.create("output.fa").with_line_width(60)
    fn with_line_width(mut slf: PyRefMut<'_, Self>, width: usize) -> PyRefMut<'_, Self> {
        if let Some(writer) = slf.inner.take() {
            slf.inner = Some(writer.with_line_width(width));
        }
        slf
    }

    /// Write a single FASTA record
    ///
    /// Args:
    ///     record (FastaRecord): Record to write
    ///
    /// Raises:
    ///     ValueError: If record is invalid
    ///     IOError: If write fails
    ///
    /// Example:
    ///     >>> writer.write_record(record)
    fn write_record(&mut self, record: &PyFastaRecord) -> PyResult<()> {
        if let Some(ref mut writer) = self.inner {
            let rust_record = record.to_fasta_record();
            writer
                .write_record(&rust_record)
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                "writer already finished",
            ))
        }
    }

    /// Get the number of records written so far
    ///
    /// Returns:
    ///     int: Number of records written
    ///
    /// Example:
    ///     >>> count = writer.records_written()
    fn records_written(&self) -> PyResult<usize> {
        if let Some(ref writer) = self.inner {
            Ok(writer.records_written())
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                "writer already finished",
            ))
        }
    }

    /// Flush buffered data to disk
    ///
    /// Raises:
    ///     IOError: If flush fails
    ///
    /// Example:
    ///     >>> writer.flush()
    fn flush(&mut self) -> PyResult<()> {
        if let Some(ref mut writer) = self.inner {
            writer
                .flush()
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                "writer already finished",
            ))
        }
    }

    /// Finish writing and flush all data
    ///
    /// This method MUST be called to ensure all data is written to disk.
    ///
    /// Raises:
    ///     IOError: If finish fails
    ///
    /// Example:
    ///     >>> writer.finish()
    fn finish(&mut self) -> PyResult<()> {
        if let Some(writer) = self.inner.take() {
            writer
                .finish()
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                "writer already finished",
            ))
        }
    }

    fn __repr__(&self) -> String {
        "FastaWriter(...)".to_string()
    }
}
