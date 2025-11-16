// GenBank Format Parser
//!
//! GenBank is the NIH genetic sequence database, containing an annotated collection
//! of all publicly available DNA sequences.
//!
//! # Format Structure
//!
//! A GenBank record consists of several sections:
//! - **LOCUS**: Identifier, length, molecule type, date
//! - **DEFINITION**: Brief description
//! - **ACCESSION**: Primary accession number
//! - **VERSION**: Accession.version identifier
//! - **KEYWORDS**: Keywords for indexing
//! - **SOURCE**: Organism information
//! - **REFERENCE**: Citation information
//! - **FEATURES**: Biological annotations (genes, CDS, etc.)
//! - **ORIGIN**: Nucleotide sequence
//! - **//**: Record terminator
//!
//! # Example
//!
//! ```no_run
//! use biometal::formats::genbank::GenBankParser;
//!
//! # fn main() -> biometal::Result<()> {
//! let parser = GenBankParser::from_path("sequence.gb")?;
//!
//! for record in parser {
//!     let record = record?;
//!     println!("Locus: {}", record.locus);
//!     println!("Length: {} bp", record.sequence.len());
//!     println!("Features: {}", record.features.len());
//! }
//! # Ok(())
//! # }
//! ```

use crate::{BiometalError, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// GenBank record representing a complete entry
#[derive(Debug, Clone, PartialEq)]
pub struct GenBankRecord {
    /// LOCUS identifier
    pub locus: String,
    /// Sequence length in base pairs
    pub length: usize,
    /// Molecule type (DNA, RNA, etc.)
    pub molecule_type: String,
    /// Topology (linear or circular)
    pub topology: String,
    /// Division code (PLN, BCT, VRT, etc.)
    pub division: String,
    /// Modification date
    pub date: String,
    /// Brief description
    pub definition: String,
    /// Primary accession number
    pub accession: String,
    /// Version identifier
    pub version: String,
    /// Keywords for indexing
    pub keywords: Vec<String>,
    /// Organism name
    pub organism: String,
    /// Taxonomic classification
    pub taxonomy: Vec<String>,
    /// References (citations)
    pub references: Vec<Reference>,
    /// Biological features (genes, CDS, etc.)
    pub features: Vec<Feature>,
    /// Nucleotide sequence (lowercase)
    pub sequence: String,
}

/// Reference citation
#[derive(Debug, Clone, PartialEq)]
pub struct Reference {
    /// Reference number
    pub number: usize,
    /// Base range
    pub location: String,
    /// Authors
    pub authors: String,
    /// Title
    pub title: String,
    /// Journal
    pub journal: String,
    /// PubMed ID (optional)
    pub pubmed: Option<String>,
}

/// Biological feature annotation
#[derive(Debug, Clone, PartialEq)]
pub struct Feature {
    /// Feature key (source, gene, CDS, etc.)
    pub key: String,
    /// Location (e.g., "1..100", "join(1..50,60..100)")
    pub location: String,
    /// Feature qualifiers (/gene="ABC", /product="XYZ")
    pub qualifiers: Vec<(String, String)>,
}

/// Streaming GenBank parser
///
/// Reads GenBank files record by record with constant memory usage.
pub struct GenBankParser<R: BufRead> {
    reader: R,
    current_line: String,
    done: bool,
}

impl GenBankParser<BufReader<File>> {
    /// Create parser from file path
    ///
    /// # Arguments
    ///
    /// * `path` - Path to GenBank file (.gb, .gbk, .genbank)
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::genbank::GenBankParser;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let parser = GenBankParser::from_path("sequence.gb")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path.as_ref())?;
        Ok(Self::new(BufReader::new(file)))
    }
}

impl<R: BufRead> GenBankParser<R> {
    /// Create parser from buffered reader
    pub fn new(reader: R) -> Self {
        GenBankParser {
            reader,
            current_line: String::new(),
            done: false,
        }
    }

    /// Read next line into buffer
    fn read_line(&mut self) -> Result<bool> {
        self.current_line.clear();
        let bytes_read = self
            .reader
            .read_line(&mut self.current_line)?;
        Ok(bytes_read > 0)
    }

    /// Parse LOCUS line
    fn parse_locus(&self) -> Result<(String, usize, String, String, String, String)> {
        // LOCUS format: LOCUS       NAME      LENGTH bp  MOLTYPE  TOPOLOGY DIV  DATE
        let parts: Vec<&str> = self.current_line.split_whitespace().collect();

        if parts.len() < 7 {
            return Err(BiometalError::InvalidGenBankFormat {
                msg: format!("Invalid LOCUS line: {}", self.current_line),
            });
        }

        let locus = parts[1].to_string();
        let length = parts[2]
            .parse::<usize>()
            .map_err(|_| BiometalError::InvalidGenBankFormat {
                msg: format!("Invalid length in LOCUS: {}", parts[2]),
            })?;
        let molecule_type = parts[4].to_string();
        let topology = parts[5].to_string();
        let division = parts[6].to_string();
        let date = parts.get(7).unwrap_or(&"").to_string();

        Ok((locus, length, molecule_type, topology, division, date))
    }

    /// Read multi-line field value (inline, no peeking)
    /// Returns the value and whether we should skip the next read
    fn read_continuing_lines(&mut self, initial_value: String) -> Result<(String, bool)> {
        let mut value = initial_value;
        let mut has_next_line = false;

        loop {
            if !self.read_line()? {
                break;
            }

            // Check if continuation line (starts with whitespace)
            if self.current_line.starts_with("  ") {
                value.push(' ');
                value.push_str(self.current_line.trim());
            } else {
                // Not a continuation - we've read the next section's line
                has_next_line = true;
                break;
            }
        }

        Ok((value, has_next_line))
    }

    /// Parse FEATURES section
    fn parse_features(&mut self) -> Result<Vec<Feature>> {
        let mut features = Vec::new();
        let mut current_feature: Option<Feature> = None;

        while self.read_line()? {
            let line = &self.current_line;

            // End of FEATURES section
            if line.starts_with("ORIGIN") || line.starts_with("//") {
                if let Some(feature) = current_feature.take() {
                    features.push(feature);
                }
                break;
            }

            // Skip header line
            if line.starts_with("FEATURES") {
                continue;
            }

            // Feature key and location (starts at column 5, no leading whitespace in key)
            if line.len() > 5 && !line.chars().nth(5).unwrap_or(' ').is_whitespace() {
                // Save previous feature
                if let Some(feature) = current_feature.take() {
                    features.push(feature);
                }

                // Parse new feature
                let parts: Vec<&str> = line[5..].split_whitespace().collect();
                if parts.len() >= 2 {
                    current_feature = Some(Feature {
                        key: parts[0].to_string(),
                        location: parts[1..].join(" "),
                        qualifiers: Vec::new(),
                    });
                }
            }
            // Qualifier (starts with /)
            else if line.trim().starts_with('/') {
                if let Some(ref mut feature) = current_feature {
                    let qualifier = line.trim();
                    if let Some((key, value)) = qualifier.split_once('=') {
                        let key = key[1..].to_string(); // Remove leading /
                        let value = value.trim_matches('"').to_string();
                        feature.qualifiers.push((key, value));
                    } else {
                        // Boolean qualifier (no value)
                        let key = qualifier[1..].to_string();
                        feature.qualifiers.push((key, String::new()));
                    }
                }
            }
            // Continuation of location or qualifier value
            else if let Some(ref mut feature) = current_feature {
                let trimmed = line.trim();
                if !trimmed.starts_with('/') && !trimmed.is_empty() {
                    // Continue location
                    feature.location.push(' ');
                    feature.location.push_str(trimmed);
                }
            }
        }

        if let Some(feature) = current_feature {
            features.push(feature);
        }

        Ok(features)
    }

    /// Parse ORIGIN section (sequence data)
    fn parse_origin(&mut self) -> Result<String> {
        let mut sequence = String::new();

        while self.read_line()? {
            let line = &self.current_line;

            // End of sequence
            if line.starts_with("//") {
                break;
            }

            // Skip ORIGIN header
            if line.starts_with("ORIGIN") {
                continue;
            }

            // Parse sequence line: "   1 acgtacgtac acgtacgtac"
            // Remove line numbers and whitespace, keep only bases
            for word in line.split_whitespace() {
                // Skip line numbers (digits)
                if word.chars().all(|c| c.is_ascii_digit()) {
                    continue;
                }
                // Append sequence (lowercase)
                sequence.push_str(&word.to_lowercase());
            }
        }

        Ok(sequence)
    }
}

impl<R: BufRead> Iterator for GenBankParser<R> {
    type Item = Result<GenBankRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        // Initialize record fields
        let mut locus = String::new();
        let mut length = 0;
        let mut molecule_type = String::new();
        let mut topology = String::new();
        let mut division = String::new();
        let mut date = String::new();
        let mut definition = String::new();
        let mut accession = String::new();
        let mut version = String::new();
        let mut keywords = Vec::new();
        let mut organism = String::new();
        let mut taxonomy = Vec::new();
        let mut references = Vec::new();
        let mut features = Vec::new();
        let mut sequence = String::new();

        // Read record line by line
        let mut skip_read = false;
        loop {
            // Read next line unless we've already read ahead
            if !skip_read {
                match self.read_line() {
                    Ok(false) => {
                        // EOF
                        self.done = true;
                        return None;
                    }
                    Err(e) => return Some(Err(e)),
                    Ok(true) => {}
                }
            }
            skip_read = false; // Reset flag

            let line = self.current_line.trim_end();

            // Record terminator
            if line.starts_with("//") {
                break;
            }

            // Parse each section
            if line.starts_with("LOCUS") {
                match self.parse_locus() {
                    Ok((l, len, mol, topo, div, d)) => {
                        locus = l;
                        length = len;
                        molecule_type = mol;
                        topology = topo;
                        division = div;
                        date = d;
                    }
                    Err(e) => return Some(Err(e)),
                }
            } else if line.starts_with("DEFINITION") {
                let value = line["DEFINITION".len()..].trim().to_string();
                match self.read_continuing_lines(value) {
                    Ok((v, has_next)) => {
                        definition = v;
                        skip_read = has_next;
                    }
                    Err(e) => return Some(Err(e)),
                }
            } else if line.starts_with("ACCESSION") {
                accession = line["ACCESSION".len()..].trim().to_string();
            } else if line.starts_with("VERSION") {
                version = line["VERSION".len()..].trim().to_string();
            } else if line.starts_with("KEYWORDS") {
                let kw = line["KEYWORDS".len()..].trim();
                if kw != "." {
                    keywords = kw.split(';').map(|s| s.trim().to_string()).collect();
                }
            } else if line.contains("ORGANISM") {
                organism = line.split_once("ORGANISM")
                    .map(|(_, o)| o.trim().to_string())
                    .unwrap_or_default();

                // Read taxonomy lines
                loop {
                    match self.read_line() {
                        Ok(false) => break,
                        Ok(true) => {},
                        Err(e) => return Some(Err(e)),
                    }

                    // Check if next section starts (uppercase at start)
                    if !self.current_line.is_empty()
                        && self.current_line.chars().next().unwrap().is_uppercase() {
                        skip_read = true; // We've read the next section
                        break;
                    }
                    if self.current_line.trim().ends_with('.') {
                        let tax_line = self.current_line.trim().trim_end_matches('.');
                        taxonomy.extend(tax_line.split(';').map(|s| s.trim().to_string()));
                        break;
                    }
                }
            } else if line.starts_with("FEATURES") {
                match self.parse_features() {
                    Ok(f) => features = f,
                    Err(e) => return Some(Err(e)),
                }
                skip_read = true; // parse_features reads ahead
            } else if line.starts_with("ORIGIN") {
                match self.parse_origin() {
                    Ok(seq) => sequence = seq,
                    Err(e) => return Some(Err(e)),
                }
                // ORIGIN is last section before //
                break;
            }
        }

        // Return complete record
        Some(Ok(GenBankRecord {
            locus,
            length,
            molecule_type,
            topology,
            division,
            date,
            definition,
            accession,
            version,
            keywords,
            organism,
            taxonomy,
            references,
            features,
            sequence,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_simple_genbank() {
        let parser = GenBankParser::from_path("tests/data/genbank/test_simple.gb")
            .expect("Failed to open test file");

        let records: Vec<_> = parser.collect();
        assert_eq!(records.len(), 1);

        let record = records[0].as_ref().expect("Failed to parse record");
        assert_eq!(record.locus, "TEST01");
        assert_eq!(record.length, 100);
        assert_eq!(record.accession, "TEST01");
        assert_eq!(record.organism, "Synthetic construct");
        assert!(!record.sequence.is_empty());
        assert_eq!(record.sequence.len(), 100);
    }

    #[test]
    fn test_parse_features() {
        let parser = GenBankParser::from_path("tests/data/genbank/test_simple.gb")
            .expect("Failed to open test file");

        let record = parser.into_iter().next()
            .expect("No record found")
            .expect("Failed to parse record");

        assert_eq!(record.features.len(), 3); // source, gene, CDS

        // Check feature types
        assert_eq!(record.features[0].key, "source");
        assert_eq!(record.features[1].key, "gene");
        assert_eq!(record.features[2].key, "CDS");

        // Check qualifiers
        let gene_feature = &record.features[1];
        let gene_qual = gene_feature.qualifiers.iter()
            .find(|(k, _)| k == "gene")
            .expect("No gene qualifier");
        assert_eq!(gene_qual.1, "test_gene");
    }
}
