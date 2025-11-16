//! Python bindings for CSI (Coordinate-Sorted Index) format

use crate::formats::index::CsiIndex;
use pyo3::exceptions::{PyFileNotFoundError, PyIOError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyList;
use std::path::PathBuf;

/// CSI chunk representing a contiguous range in the file
///
/// Each chunk has a start and end virtual file offset.
#[pyclass(name = "CsiChunk")]
#[derive(Clone)]
pub struct PyCsiChunk {
    /// Virtual file offset where chunk starts
    #[pyo3(get)]
    pub start: u64,
    /// Virtual file offset where chunk ends
    #[pyo3(get)]
    pub end: u64,
}

#[pymethods]
impl PyCsiChunk {
    fn __repr__(&self) -> String {
        format!("CsiChunk(start={:#x}, end={:#x})", self.start, self.end)
    }
}

/// CSI (Coordinate-Sorted Index) reader
///
/// Provides fast random access to sorted genomic files with configurable
/// binning parameters. CSI is the successor to BAI/TBI and supports larger
/// reference sequences.
///
/// # Example
///
/// ```python
/// from biometal import CsiIndex
///
/// # Load CSI index
/// index = CsiIndex.from_path("alignments.bam.csi")
///
/// print(f"Min shift: {index.min_shift}")
/// print(f"Depth: {index.depth}")
/// print(f"References: {len(index.references())}")
///
/// # Query region by index (if names not available)
/// chunks = index.query_by_index(0, 1000000, 2000000)
/// if chunks is not None:
///     print(f"Found {len(chunks)} chunks")
///
/// # Query by name (if names available in aux data)
/// chunks = index.query("chr1", 1000000, 2000000)
/// if chunks is not None:
///     for chunk in chunks:
///         print(f"Chunk: {chunk.start:#x} - {chunk.end:#x}")
/// ```
#[pyclass(name = "CsiIndex", unsendable)]
pub struct PyCsiIndex {
    inner: CsiIndex,
}

#[pymethods]
impl PyCsiIndex {
    /// Load CSI index from a file
    ///
    /// # Arguments
    ///
    /// * `path` - Path to .csi file
    ///
    /// # Returns
    ///
    /// CsiIndex object
    ///
    /// # Example
    ///
    /// ```python
    /// index = CsiIndex.from_path("alignments.bam.csi")
    /// ```
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let index = CsiIndex::from_path(&PathBuf::from(&path)).map_err(|e| match e {
            crate::error::BiometalError::Io(io_err)
                if io_err.kind() == std::io::ErrorKind::NotFound =>
            {
                PyFileNotFoundError::new_err(format!("CSI file not found: {}", path))
            }
            _ => PyIOError::new_err(e.to_string()),
        })?;

        Ok(PyCsiIndex { inner: index })
    }

    /// Get binning depth parameter
    ///
    /// Typically 14 (= 16kb base bin size)
    #[getter]
    fn min_shift(&self) -> i32 {
        self.inner.min_shift()
    }

    /// Get number of binning levels
    ///
    /// Typically 5
    #[getter]
    fn depth(&self) -> i32 {
        self.inner.depth()
    }

    /// Get auxiliary data
    ///
    /// May contain reference names or other metadata
    #[getter]
    fn aux_data(&self) -> Vec<u8> {
        self.inner.aux_data().to_vec()
    }

    /// Get number of references in the index
    ///
    /// # Example
    ///
    /// ```python
    /// index = CsiIndex.from_path("data.bam.csi")
    /// print(f"Number of references: {len(index.references())}")
    /// ```
    fn references(&self) -> usize {
        self.inner.references().len()
    }

    /// Get reference names (if available in auxiliary data)
    ///
    /// Returns list of reference names or None if not available
    ///
    /// # Example
    ///
    /// ```python
    /// index = CsiIndex.from_path("data.bam.csi")
    /// names = index.reference_names()
    /// if names:
    ///     for name in names:
    ///         print(name)
    /// ```
    fn reference_names(&self) -> Option<Vec<String>> {
        let refs = self.inner.references();
        let names: Vec<_> = refs
            .iter()
            .filter_map(|r| r.name.clone())
            .collect();

        if names.is_empty() {
            None
        } else {
            Some(names)
        }
    }

    /// Query region by reference name
    ///
    /// Returns None if reference name is not found in the index.
    ///
    /// # Arguments
    ///
    /// * `ref_name` - Reference sequence name
    /// * `start` - Start position (0-based, inclusive)
    /// * `end` - End position (0-based, exclusive)
    ///
    /// # Returns
    ///
    /// List of chunks or None if reference not found
    ///
    /// # Example
    ///
    /// ```python
    /// index = CsiIndex.from_path("data.bam.csi")
    /// chunks = index.query("chr1", 1000000, 2000000)
    /// if chunks is not None:
    ///     print(f"Found {len(chunks)} chunks")
    /// ```
    fn query(&self, ref_name: String, start: u32, end: u32) -> PyResult<Option<Vec<PyCsiChunk>>> {
        let result = self
            .inner
            .query(&ref_name, start, end)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(result.map(|chunks| {
            chunks
                .into_iter()
                .map(|chunk| PyCsiChunk {
                    start: chunk.start.as_raw(),
                    end: chunk.end.as_raw(),
                })
                .collect()
        }))
    }

    /// Query region by reference index
    ///
    /// Use this when reference names are not available in the index.
    ///
    /// # Arguments
    ///
    /// * `ref_idx` - Reference sequence index (0-based)
    /// * `start` - Start position (0-based, inclusive)
    /// * `end` - End position (0-based, exclusive)
    ///
    /// # Returns
    ///
    /// List of chunks or None if reference index out of bounds
    ///
    /// # Example
    ///
    /// ```python
    /// index = CsiIndex.from_path("data.bam.csi")
    /// chunks = index.query_by_index(0, 1000000, 2000000)
    /// if chunks is not None:
    ///     print(f"Found {len(chunks)} chunks")
    /// ```
    fn query_by_index(&self, ref_idx: usize, start: u32, end: u32) -> PyResult<Option<Vec<PyCsiChunk>>> {
        let result = self
            .inner
            .query_by_index(ref_idx, start, end)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(result.map(|chunks| {
            chunks
                .into_iter()
                .map(|chunk| PyCsiChunk {
                    start: chunk.start.as_raw(),
                    end: chunk.end.as_raw(),
                })
                .collect()
        }))
    }

    fn __repr__(&self) -> String {
        format!(
            "CsiIndex(min_shift={}, depth={}, n_ref={})",
            self.inner.min_shift(),
            self.inner.depth(),
            self.inner.references().len()
        )
    }
}
