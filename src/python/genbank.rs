//! Python bindings for GenBank format (NCBI genetic sequence database)

use pyo3::prelude::*;
use pyo3::exceptions::PyStopIteration;
use crate::formats::genbank::{GenBankParser, GenBankRecord, Reference, Feature};
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

/// GenBank feature annotation (gene, CDS, etc.)
///
/// Represents biological features like genes, coding sequences, regulatory elements.
///
/// Attributes:
///     key (str): Feature type (source, gene, CDS, exon, etc.)
///     location (str): Genomic location (e.g., "1..100", "join(1..50,60..100)")
///     qualifiers (list[tuple[str, str]]): Feature attributes (/gene="ABC", /product="XYZ")
///
/// Example:
///     >>> for record in parser:
///     ...     for feature in record.features:
///     ...         if feature.key == "CDS":
///     ...             gene = next((v for k, v in feature.qualifiers if k == "gene"), None)
///     ...             print(f"Gene: {gene}, Location: {feature.location}")
#[pyclass(name = "GenBankFeature")]
#[derive(Clone)]
pub struct PyGenBankFeature {
    #[pyo3(get)]
    pub key: String,
    #[pyo3(get)]
    pub location: String,
    #[pyo3(get)]
    pub qualifiers: Vec<(String, String)>,
}

impl From<Feature> for PyGenBankFeature {
    fn from(f: Feature) -> Self {
        PyGenBankFeature {
            key: f.key,
            location: f.location,
            qualifiers: f.qualifiers,
        }
    }
}

#[pymethods]
impl PyGenBankFeature {
    fn __repr__(&self) -> String {
        format!("GenBankFeature(key='{}', location='{}')", self.key, self.location)
    }

    /// Get qualifier value by key
    ///
    /// Args:
    ///     key (str): Qualifier key to search for
    ///
    /// Returns:
    ///     str | None: Qualifier value if found, None otherwise
    ///
    /// Example:
    ///     >>> gene_name = feature.get_qualifier("gene")
    ///     >>> product = feature.get_qualifier("product")
    fn get_qualifier(&self, key: &str) -> Option<String> {
        self.qualifiers
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.clone())
    }
}

/// GenBank reference citation
///
/// Represents a literature reference or direct submission.
///
/// Attributes:
///     number (int): Reference number
///     location (str): Base range this reference applies to
///     authors (str): Author list
///     title (str): Publication title
///     journal (str): Journal or submission info
///     pubmed (str | None): PubMed ID if available
///
/// Example:
///     >>> for reference in record.references:
///     ...     print(f"Reference {reference.number}: {reference.title}")
///     ...     if reference.pubmed:
///     ...         print(f"  PubMed: {reference.pubmed}")
#[pyclass(name = "GenBankReference")]
#[derive(Clone)]
pub struct PyGenBankReference {
    #[pyo3(get)]
    pub number: usize,
    #[pyo3(get)]
    pub location: String,
    #[pyo3(get)]
    pub authors: String,
    #[pyo3(get)]
    pub title: String,
    #[pyo3(get)]
    pub journal: String,
    #[pyo3(get)]
    pub pubmed: Option<String>,
}

impl From<Reference> for PyGenBankReference {
    fn from(r: Reference) -> Self {
        PyGenBankReference {
            number: r.number,
            location: r.location,
            authors: r.authors,
            title: r.title,
            journal: r.journal,
            pubmed: r.pubmed,
        }
    }
}

#[pymethods]
impl PyGenBankReference {
    fn __repr__(&self) -> String {
        format!("GenBankReference(number={}, title='{}')", self.number, self.title)
    }
}

/// GenBank record representing a complete sequence entry
///
/// Contains all information from a GenBank record: metadata, annotations, and sequence.
///
/// Attributes:
///     locus (str): LOCUS identifier
///     length (int): Sequence length in base pairs
///     molecule_type (str): Molecule type (DNA, RNA, etc.)
///     topology (str): Topology (linear or circular)
///     division (str): GenBank division (PLN, BCT, VRT, etc.)
///     date (str): Modification date
///     definition (str): Brief sequence description
///     accession (str): Primary accession number
///     version (str): Version identifier (accession.version)
///     keywords (list[str]): Keywords for indexing
///     organism (str): Organism name
///     taxonomy (list[str]): Taxonomic classification
///     references (list[GenBankReference]): Literature references
///     features (list[GenBankFeature]): Biological features (genes, CDS, etc.)
///     sequence (str): Nucleotide sequence (lowercase)
///
/// Example:
///     >>> for record in parser:
///     ...     print(f"Locus: {record.locus}")
///     ...     print(f"Organism: {record.organism}")
///     ...     print(f"Length: {record.length} bp")
///     ...     print(f"Features: {len(record.features)}")
///     ...     print(f"GC content: {(record.sequence.count('g') + record.sequence.count('c')) / len(record.sequence) * 100:.1f}%")
#[pyclass(name = "GenBankRecord")]
#[derive(Clone)]
pub struct PyGenBankRecord {
    #[pyo3(get)]
    pub locus: String,
    #[pyo3(get)]
    pub length: usize,
    #[pyo3(get)]
    pub molecule_type: String,
    #[pyo3(get)]
    pub topology: String,
    #[pyo3(get)]
    pub division: String,
    #[pyo3(get)]
    pub date: String,
    #[pyo3(get)]
    pub definition: String,
    #[pyo3(get)]
    pub accession: String,
    #[pyo3(get)]
    pub version: String,
    #[pyo3(get)]
    pub keywords: Vec<String>,
    #[pyo3(get)]
    pub organism: String,
    #[pyo3(get)]
    pub taxonomy: Vec<String>,
    #[pyo3(get)]
    pub references: Vec<PyGenBankReference>,
    #[pyo3(get)]
    pub features: Vec<PyGenBankFeature>,
    #[pyo3(get)]
    pub sequence: String,
}

impl From<GenBankRecord> for PyGenBankRecord {
    fn from(r: GenBankRecord) -> Self {
        PyGenBankRecord {
            locus: r.locus,
            length: r.length,
            molecule_type: r.molecule_type,
            topology: r.topology,
            division: r.division,
            date: r.date,
            definition: r.definition,
            accession: r.accession,
            version: r.version,
            keywords: r.keywords,
            organism: r.organism,
            taxonomy: r.taxonomy,
            references: r.references.into_iter().map(Into::into).collect(),
            features: r.features.into_iter().map(Into::into).collect(),
            sequence: r.sequence,
        }
    }
}

#[pymethods]
impl PyGenBankRecord {
    fn __repr__(&self) -> String {
        format!(
            "GenBankRecord(locus='{}', length={}, organism='{}')",
            self.locus, self.length, self.organism
        )
    }

    /// Get features by type
    ///
    /// Args:
    ///     key (str): Feature key to filter by (e.g., "gene", "CDS", "exon")
    ///
    /// Returns:
    ///     list[GenBankFeature]: Features matching the specified type
    ///
    /// Example:
    ///     >>> genes = record.get_features_by_type("gene")
    ///     >>> cds = record.get_features_by_type("CDS")
    fn get_features_by_type(&self, key: &str) -> Vec<PyGenBankFeature> {
        self.features
            .iter()
            .filter(|f| f.key == key)
            .cloned()
            .collect()
    }
}

/// Read GenBank files with streaming architecture
///
/// Streaming GenBank parser for NCBI sequence database format. Provides constant
/// memory usage regardless of file size through iterator-based architecture.
///
/// Args:
///     path (str): Path to GenBank file (.gb, .gbk, .genbank)
///
/// Features:
///     - Streaming architecture (constant memory)
///     - Complete GenBank format support
///     - Multi-record file support
///     - Feature table parsing (genes, CDS, regulatory elements)
///     - Reference parsing (literature citations)
///     - Taxonomic classification
///
/// Example:
///     >>> import biometal
///     >>> parser = biometal.GenBankParser.from_path("sequence.gb")
///     >>>
///     >>> # Stream records with constant memory
///     >>> for record in parser:
///     ...     print(f"Locus: {record.locus}")
///     ...     print(f"Organism: {record.organism}")
///     ...     print(f"Sequence length: {record.length} bp")
///     ...
///     ...     # Extract genes
///     ...     genes = record.get_features_by_type("gene")
///     ...     print(f"Genes: {len(genes)}")
///     ...
///     ...     # Analyze sequence
///     ...     gc_content = (record.sequence.count('g') + record.sequence.count('c')) / len(record.sequence)
///     ...     print(f"GC content: {gc_content * 100:.1f}%")
///     >>>
///     >>> # Find specific features
///     >>> for record in parser:
///     ...     for feature in record.features:
///     ...         if feature.key == "CDS":
///     ...             gene = feature.get_qualifier("gene")
///     ...             product = feature.get_qualifier("product")
///     ...             print(f"Gene: {gene}, Product: {product}")
///
/// Note:
///     - Memory footprint remains constant at ~5 MB
///     - Supports multi-record GenBank files
///     - Parses all GenBank sections (LOCUS, FEATURES, ORIGIN, etc.)
#[pyclass(name = "GenBankParser", unsendable)]
pub struct PyGenBankParser {
    inner: Option<GenBankParser<BufReader<File>>>,
}

#[pymethods]
impl PyGenBankParser {
    /// Open GenBank file for streaming
    ///
    /// Args:
    ///     path (str): Path to GenBank file
    ///
    /// Returns:
    ///     GenBankParser: Streaming iterator with constant memory
    ///
    /// Raises:
    ///     IOError: If file cannot be opened
    ///     ValueError: If file format is invalid
    ///
    /// Example:
    ///     >>> parser = biometal.GenBankParser.from_path("sequence.gb")
    ///     >>> for record in parser:
    ///     ...     print(f"{record.locus}: {record.length} bp")
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let parser = GenBankParser::from_path(&PathBuf::from(path))
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        Ok(PyGenBankParser {
            inner: Some(parser),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyGenBankRecord> {
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
        "GenBankParser()".to_string()
    }
}
