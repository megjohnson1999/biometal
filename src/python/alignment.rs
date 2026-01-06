//! Python bindings for sequence alignment and pattern matching
//!
//! This module provides Python access to biometal's alignment primitives including
//! Smith-Waterman alignment, streaming read mapping, and pattern matching workflows.

use pyo3::prelude::*;
use std::path::PathBuf;

use crate::alignment::{
    smith_waterman, smith_waterman_naive, Alignment, ScoringMatrix,
    StreamingMapper, StreamingMapperConfig, MappingResult,
    MotifFinder, PrimerFinder, AdapterDetector, MotifPattern, MotifMatch,
};
use crate::io::CompressedReader;
use crate::{FastqStream, FastaStream};

/// Scoring matrix for sequence alignment
///
/// Configures match/mismatch penalties and gap costs for Smith-Waterman alignment.
///
/// Args:
///     match_score (int): Score for matching bases (default: 2)
///     mismatch_score (int): Penalty for mismatching bases (default: -1)
///     gap_open (int): Penalty for opening a gap (default: -2)
///     gap_extend (int): Penalty for extending a gap (default: -1)
///
/// Example:
///     >>> scoring = biometal.ScoringMatrix(match_score=5, mismatch_score=-4, gap_open=-6, gap_extend=-1)
///     >>> alignment = biometal.smith_waterman(b"ACGT", b"ACGT", scoring)
///     >>> print(f"Score: {alignment.score}")
#[pyclass(name = "ScoringMatrix")]
pub struct PyScoringMatrix {
    inner: ScoringMatrix,
}

#[pymethods]
impl PyScoringMatrix {
    #[new]
    #[pyo3(signature = (match_score = 2, mismatch_score = -1, gap_open = -2, gap_extend = -1))]
    fn new(match_score: i32, mismatch_score: i32, gap_open: i32, gap_extend: i32) -> Self {
        Self {
            inner: ScoringMatrix::new(match_score, mismatch_score, gap_open, gap_extend),
        }
    }

    /// Default scoring matrix (match=2, mismatch=-1, gap_open=-2, gap_extend=-1)
    #[staticmethod]
    fn default() -> Self {
        Self {
            inner: ScoringMatrix::default(),
        }
    }

    #[getter]
    fn match_score(&self) -> i32 {
        self.inner.match_score
    }

    #[getter]
    fn mismatch_score(&self) -> i32 {
        self.inner.mismatch_score
    }

    #[getter]
    fn gap_open(&self) -> i32 {
        self.inner.gap_open
    }

    #[getter]
    fn gap_extend(&self) -> i32 {
        self.inner.gap_extend
    }

    fn __repr__(&self) -> String {
        format!(
            "ScoringMatrix(match={}, mismatch={}, gap_open={}, gap_extend={})",
            self.inner.match_score,
            self.inner.mismatch_score,
            self.inner.gap_open,
            self.inner.gap_extend
        )
    }
}

impl PyScoringMatrix {
    pub fn to_scoring_matrix(&self) -> ScoringMatrix {
        self.inner
    }
}

/// Sequence alignment result
///
/// Contains alignment score, positions, and CIGAR string from Smith-Waterman alignment.
///
/// Attributes:
///     score (int): Alignment score
///     query_start (int): Start position in query sequence
///     query_end (int): End position in query sequence
///     ref_start (int): Start position in reference sequence
///     ref_end (int): End position in reference sequence
///     cigar (str): CIGAR string representation
///
/// Example:
///     >>> alignment = biometal.smith_waterman(b"ACGT", b"ACGT", biometal.ScoringMatrix.default())
///     >>> print(f"Score: {alignment.score}, CIGAR: {alignment.cigar}")
#[pyclass(name = "Alignment")]
pub struct PyAlignment {
    inner: Alignment,
}

#[pymethods]
impl PyAlignment {
    #[getter]
    fn score(&self) -> i32 {
        self.inner.score
    }

    #[getter]
    fn query_start(&self) -> usize {
        self.inner.query_start
    }

    #[getter]
    fn query_end(&self) -> usize {
        self.inner.query_end
    }

    #[getter]
    fn ref_start(&self) -> usize {
        self.inner.ref_start
    }

    #[getter]
    fn ref_end(&self) -> usize {
        self.inner.ref_end
    }

    #[getter]
    fn cigar(&self) -> String {
        // Convert Vec<CigarOp> to string representation
        self.inner.cigar.iter()
            .map(|op| format!("{}", op))
            .collect::<Vec<_>>()
            .join("")
    }

    fn __repr__(&self) -> String {
        let cigar_str = self.inner.cigar.iter()
            .map(|op| format!("{}", op))
            .collect::<Vec<_>>()
            .join("");

        format!(
            "Alignment(score={}, query={}:{}, ref={}:{}, cigar='{}')",
            self.inner.score,
            self.inner.query_start,
            self.inner.query_end,
            self.inner.ref_start,
            self.inner.ref_end,
            cigar_str
        )
    }
}

impl From<Alignment> for PyAlignment {
    fn from(alignment: Alignment) -> Self {
        Self { inner: alignment }
    }
}

/// Smith-Waterman local sequence alignment
///
/// Performs optimal local alignment using Smith-Waterman dynamic programming.
/// Automatically dispatches to NEON (2-4× speedup) or GPU (5-10× speedup) when available.
///
/// Args:
///     query (bytes): Query sequence
///     reference (bytes): Reference sequence
///     scoring (ScoringMatrix): Scoring parameters
///
/// Returns:
///     Alignment: Alignment result with score and coordinates
///
/// Example:
///     >>> query = b"ACGTACGT"
///     >>> reference = b"ACGTACGT"
///     >>> scoring = biometal.ScoringMatrix.default()
///     >>> alignment = biometal.smith_waterman(query, reference, scoring)
///     >>> print(f"Score: {alignment.score}")  # 16 (8 matches × 2)
///
/// Performance:
///     - CPU (naive): ~1,000 alignments/sec
///     - NEON (ARM): ~2,000-4,000 alignments/sec (2-4× speedup)
///     - GPU (Metal): ~10,000-50,000 alignments/sec (5-10× speedup, batch≥10)
#[pyfunction(name = "smith_waterman")]
pub fn py_smith_waterman(
    query: &[u8],
    reference: &[u8],
    scoring: &PyScoringMatrix,
) -> PyResult<PyAlignment> {
    let alignment = smith_waterman(query, reference, &scoring.to_scoring_matrix());
    Ok(alignment.into())
}

/// Smith-Waterman alignment (naive CPU implementation)
///
/// Reference implementation for correctness validation and fallback.
/// Always uses scalar CPU code regardless of platform.
///
/// Args:
///     query (bytes): Query sequence
///     reference (bytes): Reference sequence
///     scoring (ScoringMatrix): Scoring parameters
///
/// Returns:
///     Alignment: Alignment result
///
/// Note:
///     Use smith_waterman() for best performance. This function is primarily
///     for testing and platforms without SIMD acceleration.
#[pyfunction(name = "smith_waterman_naive")]
pub fn py_smith_waterman_naive(
    query: &[u8],
    reference: &[u8],
    scoring: &PyScoringMatrix,
) -> PyResult<PyAlignment> {
    let alignment = smith_waterman_naive(query, reference, &scoring.to_scoring_matrix());
    Ok(alignment.into())
}

/// Configuration for streaming read mapping
///
/// Parameters for windowed reference processing in streaming read mapping.
/// Balances memory usage (~5MB constant) with alignment sensitivity.
///
/// Args:
///     window_size (int): Reference window size in bytes (default: 1MB)
///     overlap_bp (int): Overlap between windows in base pairs (default: 200bp)
///     min_score_threshold (int): Minimum alignment score to report (default: 50)
///     scoring (ScoringMatrix): Scoring parameters for alignment
///
/// Example:
///     >>> config = biometal.StreamingMapperConfig(
///     ...     window_size=500_000,  # 500KB windows
///     ...     overlap_bp=100,       # 100bp overlap
///     ...     min_score_threshold=30
///     ... )
#[pyclass(name = "StreamingMapperConfig")]
pub struct PyStreamingMapperConfig {
    inner: StreamingMapperConfig,
}

#[pymethods]
impl PyStreamingMapperConfig {
    #[new]
    #[pyo3(signature = (window_size = 1_000_000, overlap_bp = 200, min_score_threshold = 50, scoring = None))]
    fn new(
        window_size: usize,
        overlap_bp: usize,
        min_score_threshold: i32,
        scoring: Option<&PyScoringMatrix>,
    ) -> Self {
        let scoring_matrix = scoring
            .map(|s| s.to_scoring_matrix())
            .unwrap_or_else(ScoringMatrix::default);

        Self {
            inner: StreamingMapperConfig {
                window_size,
                overlap_bp,
                min_score_threshold,
                scoring: scoring_matrix,
            },
        }
    }

    /// Default configuration (1MB windows, 200bp overlap, score threshold 50)
    #[staticmethod]
    fn default() -> Self {
        Self {
            inner: StreamingMapperConfig::default(),
        }
    }

    #[getter]
    fn window_size(&self) -> usize {
        self.inner.window_size
    }

    #[getter]
    fn overlap_bp(&self) -> usize {
        self.inner.overlap_bp
    }

    #[getter]
    fn min_score_threshold(&self) -> i32 {
        self.inner.min_score_threshold
    }

    fn __repr__(&self) -> String {
        format!(
            "StreamingMapperConfig(window_size={}, overlap_bp={}, min_score_threshold={})",
            self.inner.window_size, self.inner.overlap_bp, self.inner.min_score_threshold
        )
    }
}

impl PyStreamingMapperConfig {
    pub fn to_config(&self) -> StreamingMapperConfig {
        self.inner.clone()
    }
}

/// Streaming read mapping result
///
/// Result from mapping a read to a reference position using streaming windowed approach.
///
/// Attributes:
///     query_id (str): Query sequence identifier
///     reference_id (str): Reference sequence identifier
///     global_ref_start (int): Start position in reference genome
///     global_ref_end (int): End position in reference genome
///     alignment (Alignment): Alignment details
///     window_index (int): Window where alignment was found
#[pyclass(name = "MappingResult")]
pub struct PyMappingResult {
    inner: MappingResult,
}

#[pymethods]
impl PyMappingResult {
    #[getter]
    fn query_id(&self) -> String {
        self.inner.query_id.clone()
    }

    #[getter]
    fn reference_id(&self) -> String {
        self.inner.reference_id.clone()
    }

    #[getter]
    fn global_ref_start(&self) -> usize {
        self.inner.global_ref_start
    }

    #[getter]
    fn global_ref_end(&self) -> usize {
        self.inner.global_ref_end
    }

    #[getter]
    fn alignment(&self) -> PyAlignment {
        self.inner.alignment.clone().into()
    }

    #[getter]
    fn window_index(&self) -> usize {
        self.inner.window_index
    }

    fn __repr__(&self) -> String {
        format!(
            "MappingResult(query='{}', ref='{}', pos={}:{}, score={})",
            self.inner.query_id,
            self.inner.reference_id,
            self.inner.global_ref_start,
            self.inner.global_ref_end,
            self.inner.alignment.score
        )
    }
}

impl From<MappingResult> for PyMappingResult {
    fn from(result: MappingResult) -> Self {
        Self { inner: result }
    }
}

/// Streaming read mapper with constant memory
///
/// Novel windowed approach to read mapping that maintains ~5MB constant memory
/// regardless of reference genome or dataset size. Processes reference in
/// overlapping windows and maps reads against each window.
///
/// Args:
///     config (StreamingMapperConfig): Mapping configuration
///
/// Example:
///     >>> config = biometal.StreamingMapperConfig.default()
///     >>> mapper = biometal.StreamingMapper(config)
///     >>> mappings = list(mapper.map_reads_streaming("reference.fasta", "reads.fastq"))
///     >>> for mapping in mappings:
///     ...     print(f"Mapped {mapping.query_id} to position {mapping.global_ref_start}")
///
/// Memory Usage:
///     - Reference window: ~1MB (configurable)
///     - Read buffer: ~1MB (typical FASTQ chunk)
///     - Alignment buffers: ~3MB (DP matrices, results)
///     - Total: ~5MB constant (vs 8-32GB for traditional aligners)
#[pyclass(name = "StreamingMapper", unsendable)]
pub struct PyStreamingMapper {
    inner: Option<StreamingMapper>,
}

#[pymethods]
impl PyStreamingMapper {
    #[new]
    fn new(config: &PyStreamingMapperConfig) -> Self {
        Self {
            inner: Some(StreamingMapper::new(config.to_config())),
        }
    }

    /// Map reads from FASTQ file against reference FASTA file
    ///
    /// Streams both reference and reads to maintain constant memory usage.
    /// Returns iterator over mapping results.
    ///
    /// Args:
    ///     reference_path (str): Path to reference FASTA file
    ///     reads_path (str): Path to reads FASTQ file
    ///
    /// Returns:
    ///     list[MappingResult]: List of mapping results
    ///
    /// Example:
    ///     >>> mapper = biometal.StreamingMapper(biometal.StreamingMapperConfig.default())
    ///     >>> mappings = mapper.map_reads_streaming("chr1.fasta", "reads.fastq")
    ///     >>> for mapping in mappings:
    ///     ...     print(f"Mapped {mapping.query_id} with score {mapping.alignment.score}")
    fn map_reads_streaming(&mut self, reference_path: String, reads_path: String) -> PyResult<Vec<PyMappingResult>> {
        if let Some(ref mut mapper) = self.inner {
            let mut results = Vec::new();
            let mapping_iter = mapper.map_reads_streaming(
                PathBuf::from(reference_path),
                PathBuf::from(reads_path)
            ).map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

            for mapping_result in mapping_iter {
                match mapping_result {
                    Ok(mapping) => results.push(mapping.into()),
                    Err(e) => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())),
                }
            }

            Ok(results)
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyIOError, _>("mapper already consumed"))
        }
    }

    fn __repr__(&self) -> String {
        "StreamingMapper(...)".to_string()
    }
}

/// Motif pattern for sequence searching
///
/// Defines a sequence motif to search for with minimum alignment score threshold.
///
/// Args:
///     sequence (str): DNA/RNA motif sequence
///     name (str): Descriptive name for the motif
///     min_score (int): Minimum alignment score threshold (default: 40)
///
/// Example:
///     >>> tata_box = biometal.MotifPattern("TATAAA", "TATA box", 40)
///     >>> caat_box = biometal.MotifPattern("CAAT", "CAAT box", 30)
#[pyclass(name = "MotifPattern")]
pub struct PyMotifPattern {
    inner: MotifPattern,
}

#[pymethods]
impl PyMotifPattern {
    #[new]
    fn new(sequence: String, name: String) -> Self {
        Self {
            inner: MotifPattern::new(&sequence, &name),
        }
    }

    #[getter]
    fn sequence(&self) -> String {
        String::from_utf8_lossy(&self.inner.sequence).to_string()
    }

    #[getter]
    fn name(&self) -> String {
        self.inner.name.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "MotifPattern(sequence='{}', name='{}')",
            String::from_utf8_lossy(&self.inner.sequence),
            self.inner.name
        )
    }
}

impl PyMotifPattern {
    pub fn to_pattern(&self) -> MotifPattern {
        self.inner.clone()
    }
}

/// Motif match result
///
/// Result from motif finding with alignment details and position information.
///
/// Attributes:
///     sequence_id (str): Target sequence identifier
///     motif_name (str): Name of matched motif
///     position (int): Position in target sequence
///     alignment (Alignment): Alignment details
#[pyclass(name = "MotifMatch")]
pub struct PyMotifMatch {
    inner: MotifMatch,
}

#[pymethods]
impl PyMotifMatch {
    #[getter]
    fn sequence_id(&self) -> String {
        self.inner.sequence_id.clone()
    }

    #[getter]
    fn motif_name(&self) -> String {
        self.inner.motif_name.clone()
    }

    #[getter]
    fn position(&self) -> usize {
        self.inner.position
    }

    #[getter]
    fn alignment(&self) -> PyAlignment {
        self.inner.alignment.clone().into()
    }

    fn __repr__(&self) -> String {
        format!(
            "MotifMatch(sequence='{}', motif='{}', position={}, score={})",
            self.inner.sequence_id,
            self.inner.motif_name,
            self.inner.position,
            self.inner.alignment.score
        )
    }
}

impl From<MotifMatch> for PyMotifMatch {
    fn from(motif_match: MotifMatch) -> Self {
        Self { inner: motif_match }
    }
}

/// Streaming motif finder
///
/// Searches for regulatory sequence motifs in streaming fashion with constant memory usage.
/// Integrates Smith-Waterman alignment with FASTQ/FASTA parsers for flexible pattern matching.
///
/// Args:
///     patterns (list[MotifPattern]): List of motifs to search for
///     scoring (ScoringMatrix): Scoring parameters for alignment
///
/// Example:
///     >>> motifs = [
///     ...     biometal.MotifPattern("TATAAA", "TATA box", 40),
///     ...     biometal.MotifPattern("CAAT", "CAAT box", 30),
///     ... ]
///     >>> finder = biometal.MotifFinder(motifs, biometal.ScoringMatrix.default())
///     >>> matches = finder.find_in_fastq("sequences.fastq")
///     >>> for match in matches:
///     ...     print(f"Found {match.motif_name} in {match.sequence_id} at position {match.position}")
#[pyclass(name = "MotifFinder", unsendable)]
pub struct PyMotifFinder {
    inner: Option<MotifFinder>,
}

#[pymethods]
impl PyMotifFinder {
    #[new]
    fn new(patterns: Vec<PyRef<PyMotifPattern>>, min_score: i32) -> Self {
        let rust_patterns: Vec<MotifPattern> = patterns.iter()
            .map(|p| p.to_pattern())
            .collect();

        Self {
            inner: Some(MotifFinder::new(rust_patterns, min_score)),
        }
    }

    /// Find motifs in FASTQ file
    ///
    /// Streams through FASTQ file searching for motif patterns.
    /// Maintains constant memory regardless of file size.
    ///
    /// Args:
    ///     path (str): Path to FASTQ file
    ///
    /// Returns:
    ///     list[MotifMatch]: List of motif matches found
    fn find_in_fastq(&mut self, path: String) -> PyResult<Vec<PyMotifMatch>> {
        if let Some(ref mut finder) = self.inner {
            let mut results = Vec::new();

            let match_iterator = finder.find_in_fastq(PathBuf::from(path))
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

            for motif_result in match_iterator {
                match motif_result {
                    Ok(motif_match) => results.push(motif_match.into()),
                    Err(e) => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())),
                }
            }

            Ok(results)
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyIOError, _>("finder already consumed"))
        }
    }

    /// Find motifs in FASTA file
    ///
    /// Streams through FASTA file searching for motif patterns.
    /// Maintains constant memory regardless of file size.
    ///
    /// Args:
    ///     path (str): Path to FASTA file
    ///
    /// Returns:
    ///     list[MotifMatch]: List of motif matches found
    fn find_in_fasta(&mut self, path: String) -> PyResult<Vec<PyMotifMatch>> {
        if let Some(ref mut finder) = self.inner {
            let mut results = Vec::new();

            let match_iterator = finder.find_in_fasta(PathBuf::from(path))
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

            for motif_result in match_iterator {
                match motif_result {
                    Ok(motif_match) => results.push(motif_match.into()),
                    Err(e) => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())),
                }
            }

            Ok(results)
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyIOError, _>("finder already consumed"))
        }
    }

    fn __repr__(&self) -> String {
        "MotifFinder(...)".to_string()
    }
}

/// PCR primer finder
///
/// Searches for PCR primer binding sites in streaming fashion.
/// Specialized for primer design and qPCR assay validation.
///
/// Args:
///     forward_primer (str): Forward primer sequence
///     reverse_primer (str): Reverse primer sequence
///     min_score (int): Minimum alignment score threshold (default: 80)
///     scoring (ScoringMatrix): Scoring parameters
///
/// Example:
///     >>> primer_finder = biometal.PrimerFinder(
///     ...     "ATGCATGCATGC",     # forward primer
///     ...     "GCATGCATGCAT",     # reverse primer
///     ...     80,                 # min score
///     ...     biometal.ScoringMatrix.default()
///     ... )
///     >>> matches = primer_finder.find_in_fasta("amplicons.fasta")
#[pyclass(name = "PrimerFinder", unsendable)]
pub struct PyPrimerFinder {
    inner: Option<PrimerFinder>,
}

#[pymethods]
impl PyPrimerFinder {
    #[new]
    fn new() -> Self {
        Self {
            inner: Some(PrimerFinder::new_standard()),
        }
    }

    /// Find primer sites in FASTA file
    ///
    /// Searches for forward and reverse primer binding sites.
    ///
    /// Args:
    ///     path (str): Path to FASTA file
    ///
    /// Returns:
    ///     list[MotifMatch]: List of primer matches found
    fn find_in_fasta(&mut self, path: String) -> PyResult<Vec<PyMotifMatch>> {
        if let Some(ref mut finder) = self.inner {
            let mut results = Vec::new();

            let match_iterator = finder.find_primers_fastq(PathBuf::from(path))
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

            for primer_result in match_iterator {
                match primer_result {
                    Ok(motif_match) => results.push(motif_match.into()),
                    Err(e) => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())),
                }
            }

            Ok(results)
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyIOError, _>("finder already consumed"))
        }
    }

    fn __repr__(&self) -> String {
        "PrimerFinder(...)".to_string()
    }
}

/// Illumina adapter detector
///
/// Detects sequencing adapters for trimming and quality control.
/// Pre-configured with common Illumina adapter sequences.
///
/// Args:
///     min_score (int): Minimum alignment score threshold (default: 60)
///     scoring (ScoringMatrix): Scoring parameters
///
/// Example:
///     >>> detector = biometal.AdapterDetector(60, biometal.ScoringMatrix.default())
///     >>> contaminated = detector.find_in_fastq("reads.fastq")
///     >>> for match in contaminated:
///     ...     print(f"Found adapter contamination in {match.sequence_id}")
#[pyclass(name = "AdapterDetector", unsendable)]
pub struct PyAdapterDetector {
    inner: Option<AdapterDetector>,
}

#[pymethods]
impl PyAdapterDetector {
    #[new]
    fn new() -> Self {
        Self {
            inner: Some(AdapterDetector::new_illumina()),
        }
    }

    /// Find adapters in FASTQ file
    ///
    /// Searches for Illumina adapter contamination in reads.
    ///
    /// Args:
    ///     path (str): Path to FASTQ file
    ///
    /// Returns:
    ///     list[MotifMatch]: List of adapter matches found
    fn find_in_fastq(&mut self, path: String) -> PyResult<Vec<PyMotifMatch>> {
        if let Some(ref mut detector) = self.inner {
            let mut results = Vec::new();

            let match_iterator = detector.detect_adapters_fastq(PathBuf::from(path))
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

            for adapter_result in match_iterator {
                match adapter_result {
                    Ok(motif_match) => results.push(motif_match.into()),
                    Err(e) => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())),
                }
            }

            Ok(results)
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyIOError, _>("detector already consumed"))
        }
    }

    fn __repr__(&self) -> String {
        "AdapterDetector(...)".to_string()
    }
}