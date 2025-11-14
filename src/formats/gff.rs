//! GFF3 (General Feature Format) parser.
//!
//! GFF3 is a standard format for gene annotations and genomic features:
//! - **Genes**: Gene models with exons, introns, UTRs
//! - **Transcripts**: mRNA, tRNA, ncRNA with isoforms
//! - **Hierarchical features**: Parent-child relationships
//!
//! # Format Specification
//!
//! GFF3 uses 9 tab-delimited columns:
//! 1. **seqid**: Chromosome/contig name
//! 2. **source**: Annotation source (e.g., NCBI, Ensembl)
//! 3. **type**: Feature type (gene, mRNA, exon, CDS)
//! 4. **start**: Start position (1-based, inclusive)
//! 5. **end**: End position (1-based, inclusive)
//! 6. **score**: Confidence score (or `.`)
//! 7. **strand**: `+`, `-`, or `.`
//! 8. **phase**: CDS phase (0, 1, 2, or `.`)
//! 9. **attributes**: Semicolon-separated key=value pairs
//!
//! # Examples
//!
//! ## Basic feature record
//!
//! ```
//! use biometal::formats::gff::Gff3Record;
//! use biometal::formats::{TabDelimitedRecord, Strand};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let line = "chr1\tEnsembl\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=ABC1";
//! let record = Gff3Record::from_line(line)?;
//!
//! assert_eq!(record.seqid, "chr1");
//! assert_eq!(record.source, "Ensembl");
//! assert_eq!(record.feature_type, "gene");
//! assert_eq!(record.start, 1000);
//! assert_eq!(record.end, 2000);
//! assert_eq!(record.strand, Strand::Forward);
//! assert_eq!(record.get_id(), Some("gene1"));
//! assert_eq!(record.get_name(), Some("ABC1"));
//! # Ok(())
//! # }
//! ```
//!
//! ## Hierarchical features
//!
//! ```
//! use biometal::formats::gff::Gff3Record;
//! use biometal::formats::TabDelimitedRecord;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Gene
//! let gene = Gff3Record::from_line("chr1\t.\tgene\t1000\t5000\t.\t+\t.\tID=gene1")?;
//! assert_eq!(gene.get_id(), Some("gene1"));
//!
//! // mRNA (child of gene)
//! let mrna = Gff3Record::from_line("chr1\t.\tmRNA\t1000\t5000\t.\t+\t.\tID=mRNA1;Parent=gene1")?;
//! assert_eq!(mrna.get_parent(), Some("gene1"));
//!
//! // Exon (child of mRNA)
//! let exon = Gff3Record::from_line("chr1\t.\texon\t1000\t1500\t.\t+\t.\tID=exon1;Parent=mRNA1")?;
//! assert_eq!(exon.get_parent(), Some("mRNA1"));
//! # Ok(())
//! # }
//! ```
//!
//! ## CDS with phase
//!
//! ```
//! use biometal::formats::gff::Gff3Record;
//! use biometal::formats::TabDelimitedRecord;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let cds = Gff3Record::from_line("chr1\t.\tCDS\t1200\t1800\t.\t+\t0\tID=cds1;Parent=mRNA1")?;
//! assert_eq!(cds.feature_type, "CDS");
//! assert_eq!(cds.phase, Some(0));
//! # Ok(())
//! # }
//! ```
//!
//! ## Streaming parser
//!
//! ```no_run
//! use biometal::formats::gff::Gff3Parser;
//! use std::fs::File;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let file = File::open("genes.gff3")?;
//! let parser = Gff3Parser::new(file);
//!
//! for result in parser {
//!     let record = result?;
//!     if record.feature_type == "gene" {
//!         println!("Gene: {} at {}:{}-{}",
//!                  record.get_id().unwrap_or("unknown"),
//!                  record.seqid, record.start, record.end);
//!     }
//! }
//! # Ok(())
//! # }
//! ```

use crate::formats::primitives::{
    fields::{parse_attributes, parse_optional, parse_required, split_fields},
    GenomicInterval, Result, Strand, TabDelimitedRecord,
};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};
use std::str::FromStr;

/// GFF3 feature record.
///
/// Represents a single genomic feature (gene, exon, CDS, etc.).
///
/// # Examples
///
/// ```
/// use biometal::formats::gff::Gff3Record;
/// use biometal::formats::TabDelimitedRecord;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let line = "chr1\tEnsembl\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=ABC1";
/// let record = Gff3Record::from_line(line)?;
/// assert_eq!(record.feature_type, "gene");
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Gff3Record {
    /// Chromosome/contig name
    pub seqid: String,
    /// Annotation source (e.g., NCBI, Ensembl)
    pub source: String,
    /// Feature type (gene, mRNA, exon, CDS, etc.)
    pub feature_type: String,
    /// Start position (1-based, inclusive)
    pub start: u64,
    /// End position (1-based, inclusive)
    pub end: u64,
    /// Confidence score (None if `.`)
    pub score: Option<f64>,
    /// Strand (+, -, .)
    pub strand: Strand,
    /// CDS phase (0, 1, 2, or None if `.`)
    pub phase: Option<u8>,
    /// Attributes (key=value pairs)
    pub attributes: HashMap<String, String>,
}

impl Gff3Record {
    /// Gets the ID attribute (unique identifier).
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::formats::gff::Gff3Record;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let record = Gff3Record::from_line("chr1\t.\tgene\t100\t200\t.\t+\t.\tID=gene1;Name=ABC")?;
    /// assert_eq!(record.get_id(), Some("gene1"));
    /// # Ok(())
    /// # }
    /// ```
    pub fn get_id(&self) -> Option<&str> {
        self.attributes.get("ID").map(|s| s.as_str())
    }

    /// Gets the Parent attribute (hierarchical relationship).
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::formats::gff::Gff3Record;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let record = Gff3Record::from_line("chr1\t.\texon\t100\t200\t.\t+\t.\tID=exon1;Parent=mRNA1")?;
    /// assert_eq!(record.get_parent(), Some("mRNA1"));
    /// # Ok(())
    /// # }
    /// ```
    pub fn get_parent(&self) -> Option<&str> {
        self.attributes.get("Parent").map(|s| s.as_str())
    }

    /// Gets the Name attribute (human-readable name).
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::formats::gff::Gff3Record;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let record = Gff3Record::from_line("chr1\t.\tgene\t100\t200\t.\t+\t.\tID=gene1;Name=ABC1")?;
    /// assert_eq!(record.get_name(), Some("ABC1"));
    /// # Ok(())
    /// # }
    /// ```
    pub fn get_name(&self) -> Option<&str> {
        self.attributes.get("Name").map(|s| s.as_str())
    }

    /// Returns the genomic interval (0-based, half-open).
    ///
    /// Converts GFF3's 1-based inclusive coordinates to 0-based half-open.
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::formats::gff::Gff3Record;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let record = Gff3Record::from_line("chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1")?;
    /// let interval = record.interval()?;
    /// assert_eq!(interval.start, 999);  // 0-based
    /// assert_eq!(interval.end, 2000);   // half-open
    /// assert_eq!(interval.length(), 1001);
    /// # Ok(())
    /// # }
    /// ```
    pub fn interval(&self) -> Result<GenomicInterval> {
        // GFF3 uses 1-based inclusive coordinates
        // Convert to 0-based half-open: [start-1, end)
        GenomicInterval::new(self.seqid.clone(), self.start - 1, self.end)
    }

    /// Returns the feature length in base pairs.
    pub fn length(&self) -> u64 {
        self.end - self.start + 1
    }
}

impl TabDelimitedRecord for Gff3Record {
    fn from_line(line: &str) -> Result<Self> {
        let fields = split_fields(line, Some(9), 0)?;

        let seqid = fields[0].to_string();
        let source = fields[1].to_string();
        let feature_type = fields[2].to_string();
        let start: u64 = parse_required(fields[3], "start", 0)?;
        let end: u64 = parse_required(fields[4], "end", 0)?;
        let score: Option<f64> = parse_optional(fields[5], "score", 0)?;
        let strand = Strand::from_str(fields[6])?;

        // Parse phase (0, 1, 2, or .)
        let phase: Option<u8> = if fields[7] == "." {
            None
        } else {
            Some(parse_required(fields[7], "phase", 0)?)
        };

        // Parse attributes
        let attributes = parse_attributes(fields[8], 0)?;

        Ok(Gff3Record {
            seqid,
            source,
            feature_type,
            start,
            end,
            score,
            strand,
            phase,
            attributes,
        })
    }

    fn to_line(&self) -> String {
        let mut fields = vec![
            self.seqid.clone(),
            self.source.clone(),
            self.feature_type.clone(),
            self.start.to_string(),
            self.end.to_string(),
            self.score
                .map(|s| s.to_string())
                .unwrap_or_else(|| ".".to_string()),
            self.strand.to_string(),
            self.phase
                .map(|p| p.to_string())
                .unwrap_or_else(|| ".".to_string()),
        ];

        // Serialize attributes
        let attr_str = if self.attributes.is_empty() {
            ".".to_string()
        } else {
            let mut attr_parts: Vec<String> = self
                .attributes
                .iter()
                .map(|(k, v)| format!("{}={}", k, v))
                .collect();
            attr_parts.sort();
            attr_parts.join(";")
        };
        fields.push(attr_str);

        fields.join("\t")
    }

    fn expected_fields() -> Option<usize> {
        Some(9)
    }
}

/// GFF3 parser with header support.
///
/// Parses GFF3 files including ##gff-version and other metadata.
pub struct Gff3Parser<R: Read> {
    reader: BufReader<R>,
    line_buf: String,
    line_number: usize,
}

impl<R: Read> Gff3Parser<R> {
    /// Creates a new GFF3 parser.
    pub fn new(reader: R) -> Self {
        Gff3Parser {
            reader: BufReader::new(reader),
            line_buf: String::new(),
            line_number: 0,
        }
    }

    /// Creates a GFF3 parser from a file path.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_path(path: impl AsRef<std::path::Path>) -> Result<Gff3Parser<std::fs::File>> {
        let file = std::fs::File::open(path)?;
        Ok(Gff3Parser::new(file))
    }
}

impl<R: Read> Iterator for Gff3Parser<R> {
    type Item = Result<Gff3Record>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.line_buf.clear();
            match self.reader.read_line(&mut self.line_buf) {
                Ok(0) => return None, // EOF
                Ok(_) => {
                    self.line_number += 1;
                    let line = self.line_buf.trim_end();

                    // Skip empty lines
                    if line.is_empty() {
                        continue;
                    }

                    // Skip comment/header lines (##gff-version, ##sequence-region, etc.)
                    if line.starts_with('#') {
                        continue;
                    }

                    return Some(Gff3Record::from_line(line));
                }
                Err(e) => return Some(Err(e.into())),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gff3_basic() {
        let line = "chr1\tEnsembl\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=ABC1";
        let record = Gff3Record::from_line(line).unwrap();

        assert_eq!(record.seqid, "chr1");
        assert_eq!(record.source, "Ensembl");
        assert_eq!(record.feature_type, "gene");
        assert_eq!(record.start, 1000);
        assert_eq!(record.end, 2000);
        assert_eq!(record.strand, Strand::Forward);
        assert_eq!(record.get_id(), Some("gene1"));
    }

    #[test]
    fn test_gff3_missing_values() {
        let line = "chr1\t.\tgene\t1000\t2000\t.\t.\t.\t.";
        let record = Gff3Record::from_line(line).unwrap();

        assert_eq!(record.source, ".");
        assert_eq!(record.score, None);
        assert_eq!(record.strand, Strand::Unknown);
        assert_eq!(record.phase, None);
    }

    #[test]
    fn test_gff3_cds_phase() {
        let line = "chr1\t.\tCDS\t1000\t2000\t.\t+\t0\tID=cds1";
        let record = Gff3Record::from_line(line).unwrap();

        assert_eq!(record.feature_type, "CDS");
        assert_eq!(record.phase, Some(0));
    }

    #[test]
    fn test_gff3_attributes() {
        let line = "chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=ABC1;biotype=protein_coding";
        let record = Gff3Record::from_line(line).unwrap();

        assert_eq!(record.get_id(), Some("gene1"));
        assert_eq!(record.attributes.get("biotype"), Some(&"protein_coding".to_string()));
    }

    #[test]
    fn test_gff3_parent_relationship() {
        let mrna = Gff3Record::from_line("chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=mRNA1;Parent=gene1").unwrap();
        assert_eq!(mrna.get_parent(), Some("gene1"));
    }

    #[test]
    fn test_gff3_interval_conversion() {
        let record = Gff3Record::from_line("chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1").unwrap();
        let interval = record.interval().unwrap();

        // GFF3: 1-based inclusive [1000, 2000]
        // Converted to 0-based half-open [999, 2000)
        assert_eq!(interval.start, 999);
        assert_eq!(interval.end, 2000);
        assert_eq!(interval.length(), 1001);
    }

    #[test]
    fn test_gff3_length() {
        let record = Gff3Record::from_line("chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1").unwrap();
        assert_eq!(record.length(), 1001); // 1000-2000 inclusive
    }

    #[test]
    fn test_gff3_round_trip() {
        let original = "chr1\tEnsembl\tgene\t1000\t2000\t50\t+\t.\tID=gene1;Name=ABC1";
        let record = Gff3Record::from_line(original).unwrap();
        let output = record.to_line();

        let record2 = Gff3Record::from_line(&output).unwrap();
        assert_eq!(record, record2);
    }
}
