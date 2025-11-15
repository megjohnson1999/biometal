//! VCF (Variant Call Format) parser.
//!
//! VCF is a text format for representing genetic variants:
//! - **SNPs**: Single nucleotide polymorphisms
//! - **Indels**: Insertions and deletions
//! - **Structural variants**: Large-scale genomic alterations
//!
//! # Format Specification
//!
//! VCF files contain:
//! - **Header lines** (`##`): Metadata (fileformat, INFO, FORMAT, FILTER, contig)
//! - **Column header** (`#`): Field names and sample IDs
//! - **Data records**: Tab-delimited variant calls
//!
//! # Examples
//!
//! ## Basic variant record
//!
//! ```
//! use biometal::formats::vcf::VcfRecord;
//! use biometal::formats::TabDelimitedRecord;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let line = "chr1\t12345\trs123\tA\tT\t30.0\tPASS\tDP=100;AF=0.5";
//! let record = VcfRecord::from_line(line)?;
//!
//! assert_eq!(record.chrom, "chr1");
//! assert_eq!(record.pos, 12345);
//! assert_eq!(record.id, Some("rs123".to_string()));
//! assert_eq!(record.reference, "A");
//! assert_eq!(record.alternate.len(), 1);
//! assert_eq!(record.alternate[0], "T");
//! assert_eq!(record.quality, Some(30.0));
//! assert_eq!(record.filter, Some("PASS".to_string()));
//! # Ok(())
//! # }
//! ```
//!
//! ## INFO field parsing
//!
//! ```
//! use biometal::formats::vcf::VcfRecord;
//! use biometal::formats::TabDelimitedRecord;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let line = "chr1\t100\t.\tA\tT\t.\t.\tDP=50;AF=0.25;DB";
//! let record = VcfRecord::from_line(line)?;
//!
//! assert_eq!(record.info.get("DP"), Some(&"50".to_string()));
//! assert_eq!(record.info.get("AF"), Some(&"0.25".to_string()));
//! assert_eq!(record.info.get("DB"), Some(&"".to_string())); // Flag
//! # Ok(())
//! # }
//! ```
//!
//! ## Multiple alternate alleles
//!
//! ```
//! use biometal::formats::vcf::VcfRecord;
//! use biometal::formats::TabDelimitedRecord;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let line = "chr1\t100\t.\tA\tT,G,C\t.\t.\t.";
//! let record = VcfRecord::from_line(line)?;
//!
//! assert_eq!(record.alternate.len(), 3);
//! assert_eq!(record.alternate, vec!["T", "G", "C"]);
//! # Ok(())
//! # }
//! ```
//!
//! ## Streaming parser
//!
//! ```no_run
//! use biometal::formats::vcf::VcfParser;
//! use std::fs::File;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let file = File::open("variants.vcf")?;
//! let mut parser = VcfParser::new(file);
//!
//! // Parse header
//! let header = parser.parse_header()?;
//! println!("VCF version: {}", header.fileformat);
//! println!("Samples: {:?}", header.samples);
//!
//! // Parse variants
//! for result in parser {
//!     let record = result?;
//!     println!("{}: {} {} -> {}",
//!              record.chrom, record.pos, record.reference,
//!              record.alternate.join(","));
//! }
//! # Ok(())
//! # }
//! ```

use crate::formats::primitives::{
    fields::{parse_comma_list, parse_optional, parse_required, split_fields},
    Result, TabDelimitedRecord,
};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};

/// VCF file header with metadata.
///
/// Contains metadata lines (##) and the column header (#CHROM...).
///
/// # Examples
///
/// ```
/// use biometal::formats::vcf::VcfHeader;
///
/// let header = VcfHeader::new("VCFv4.2".to_string());
/// assert_eq!(header.fileformat, "VCFv4.2");
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct VcfHeader {
    /// VCF version (e.g., "VCFv4.2")
    pub fileformat: String,
    /// INFO field definitions (key -> description)
    pub info_fields: HashMap<String, String>,
    /// FORMAT field definitions (key -> description)
    pub format_fields: HashMap<String, String>,
    /// FILTER definitions (ID -> description)
    pub filters: HashMap<String, String>,
    /// Contig definitions (ID -> length)
    pub contigs: HashMap<String, Option<u64>>,
    /// Other metadata lines
    pub metadata: Vec<String>,
    /// Sample names (from #CHROM header)
    pub samples: Vec<String>,
}

impl VcfHeader {
    /// Creates a new VCF header with the given format version.
    pub fn new(fileformat: String) -> Self {
        VcfHeader {
            fileformat,
            info_fields: HashMap::new(),
            format_fields: HashMap::new(),
            filters: HashMap::new(),
            contigs: HashMap::new(),
            metadata: Vec::new(),
            samples: Vec::new(),
        }
    }

    /// Converts the header to VCF header lines.
    ///
    /// Returns all header lines including the column header.
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::formats::vcf::VcfHeader;
    ///
    /// let mut header = VcfHeader::new("VCFv4.2".to_string());
    /// header.samples = vec!["sample1".to_string()];
    /// let lines = header.to_header_lines();
    /// assert!(lines[0].starts_with("##fileformat="));
    /// ```
    pub fn to_header_lines(&self) -> Vec<String> {
        let mut lines = Vec::new();

        // ##fileformat must be first
        lines.push(format!("##fileformat={}", self.fileformat));

        // Write INFO fields (sorted for consistency)
        let mut info_keys: Vec<&String> = self.info_fields.keys().collect();
        info_keys.sort();
        for key in info_keys {
            let desc = &self.info_fields[key];
            lines.push(format!("##INFO=<ID={},Description=\"{}\">", key, desc));
        }

        // Write FORMAT fields (sorted)
        let mut format_keys: Vec<&String> = self.format_fields.keys().collect();
        format_keys.sort();
        for key in format_keys {
            let desc = &self.format_fields[key];
            lines.push(format!("##FORMAT=<ID={},Description=\"{}\">", key, desc));
        }

        // Write FILTER fields (sorted)
        let mut filter_keys: Vec<&String> = self.filters.keys().collect();
        filter_keys.sort();
        for key in filter_keys {
            let desc = &self.filters[key];
            lines.push(format!("##FILTER=<ID={},Description=\"{}\">", key, desc));
        }

        // Write contig fields (sorted)
        let mut contig_keys: Vec<&String> = self.contigs.keys().collect();
        contig_keys.sort();
        for key in contig_keys {
            if let Some(length) = self.contigs[key] {
                lines.push(format!("##contig=<ID={},length={}>", key, length));
            } else {
                lines.push(format!("##contig=<ID={}>", key));
            }
        }

        // Write other metadata
        for meta in &self.metadata {
            lines.push(format!("##{}", meta));
        }

        // Write column header
        let mut header_line = String::from("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        if !self.samples.is_empty() {
            header_line.push_str("\tFORMAT");
            for sample in &self.samples {
                header_line.push('\t');
                header_line.push_str(sample);
            }
        }
        lines.push(header_line);

        lines
    }

    /// Parses a VCF header metadata line.
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::formats::vcf::VcfHeader;
    ///
    /// let mut header = VcfHeader::new("VCFv4.2".to_string());
    /// header.parse_metadata_line("##fileformat=VCFv4.2");
    /// assert_eq!(header.fileformat, "VCFv4.2");
    /// ```
    pub fn parse_metadata_line(&mut self, line: &str) {
        if !line.starts_with("##") {
            return;
        }

        let line = &line[2..]; // Remove ##

        if let Some(eq_pos) = line.find('=') {
            let key = &line[..eq_pos];
            let value = &line[eq_pos + 1..];

            match key {
                "fileformat" => {
                    self.fileformat = value.to_string();
                }
                "INFO" => {
                    if let Some((id, desc)) = parse_vcf_structured_field(value) {
                        self.info_fields.insert(id, desc);
                    }
                }
                "FORMAT" => {
                    if let Some((id, desc)) = parse_vcf_structured_field(value) {
                        self.format_fields.insert(id, desc);
                    }
                }
                "FILTER" => {
                    if let Some((id, desc)) = parse_vcf_structured_field(value) {
                        self.filters.insert(id, desc);
                    }
                }
                "contig" => {
                    if let Some((id, length_str)) = parse_vcf_structured_field(value) {
                        let length = length_str.parse::<u64>().ok();
                        self.contigs.insert(id, length);
                    }
                }
                _ => {
                    self.metadata.push(line.to_string());
                }
            }
        }
    }

    /// Parses the column header line (#CHROM POS ...).
    ///
    /// Extracts sample names from the header.
    pub fn parse_column_header(&mut self, line: &str) {
        if !line.starts_with("#CHROM") {
            return;
        }

        let fields: Vec<&str> = line.split('\t').collect();

        // Standard VCF columns: #CHROM POS ID REF ALT QUAL FILTER INFO [FORMAT sample1 sample2 ...]
        // Samples start at index 9 if FORMAT is present
        if fields.len() > 9 {
            self.samples = fields[9..].iter().map(|s| s.to_string()).collect();
        }
    }
}

/// VCF variant record.
///
/// Represents a single variant with all associated information.
///
/// # Examples
///
/// ```
/// use biometal::formats::vcf::VcfRecord;
/// use biometal::formats::TabDelimitedRecord;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let line = "chr1\t12345\trs123\tA\tT\t30.0\tPASS\tDP=100";
/// let record = VcfRecord::from_line(line)?;
///
/// assert_eq!(record.chrom, "chr1");
/// assert_eq!(record.pos, 12345);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct VcfRecord {
    /// Chromosome/contig name
    pub chrom: String,
    /// Position (1-based)
    pub pos: u64,
    /// Variant ID (e.g., dbSNP rsID), None if "."
    pub id: Option<String>,
    /// Reference allele
    pub reference: String,
    /// Alternate alleles
    pub alternate: Vec<String>,
    /// Phred-scaled quality score, None if "."
    pub quality: Option<f64>,
    /// Filter status (PASS or filter IDs), None if "."
    pub filter: Option<String>,
    /// INFO field key-value pairs
    pub info: HashMap<String, String>,
    /// FORMAT field (optional)
    pub format: Option<String>,
    /// Sample genotypes (optional, one per sample)
    pub samples: Vec<String>,
}

impl TabDelimitedRecord for VcfRecord {
    fn from_line(line: &str) -> Result<Self> {
        let fields = split_fields(line, Some(8), 0)?;

        let chrom = fields[0].to_string();
        let pos: u64 = parse_required(fields[1], "POS", 0)?;
        let id: Option<String> = parse_optional(fields[2], "ID", 0)?;
        let reference = fields[3].to_string();

        // Parse alternate alleles (comma-separated)
        let alternate: Vec<String> = parse_comma_list(fields[4])
            .iter()
            .map(|s| s.to_string())
            .collect();

        let quality: Option<f64> = parse_optional(fields[5], "QUAL", 0)?;
        let filter: Option<String> = parse_optional(fields[6], "FILTER", 0)?;

        // Parse INFO field
        let info = parse_info_field(fields[7]);

        // Parse optional FORMAT and sample columns
        let (format, samples) = if fields.len() > 8 {
            let format = if fields[8] == "." {
                None
            } else {
                Some(fields[8].to_string())
            };
            let samples = fields[9..]
                .iter()
                .map(|s| s.to_string())
                .collect();
            (format, samples)
        } else {
            (None, Vec::new())
        };

        Ok(VcfRecord {
            chrom,
            pos,
            id,
            reference,
            alternate,
            quality,
            filter,
            info,
            format,
            samples,
        })
    }

    fn to_line(&self) -> String {
        let mut fields = vec![
            self.chrom.clone(),
            self.pos.to_string(),
            self.id.as_deref().unwrap_or(".").to_string(),
            self.reference.clone(),
            if self.alternate.is_empty() {
                ".".to_string()
            } else {
                self.alternate.join(",")
            },
            self.quality
                .map(|q| q.to_string())
                .unwrap_or_else(|| ".".to_string()),
            self.filter.as_deref().unwrap_or(".").to_string(),
        ];

        // Serialize INFO field
        let info_str = if self.info.is_empty() {
            ".".to_string()
        } else {
            let mut info_parts: Vec<String> = self
                .info
                .iter()
                .map(|(k, v)| {
                    if v.is_empty() {
                        k.clone()
                    } else {
                        format!("{}={}", k, v)
                    }
                })
                .collect();
            info_parts.sort();
            info_parts.join(";")
        };
        fields.push(info_str);

        // Add FORMAT and samples if present
        if let Some(format) = &self.format {
            fields.push(format.clone());
            fields.extend(self.samples.iter().cloned());
        }

        fields.join("\t")
    }

    fn expected_fields() -> Option<usize> {
        Some(8) // Minimum required fields
    }
}

/// VCF parser with header support.
///
/// Parses VCF files including header metadata and variant records.
pub struct VcfParser<R: Read> {
    reader: BufReader<R>,
    line_buf: String,
    line_number: usize,
    header_parsed: bool,
}

impl<R: Read> VcfParser<R> {
    /// Creates a new VCF parser.
    pub fn new(reader: R) -> Self {
        VcfParser {
            reader: BufReader::new(reader),
            line_buf: String::new(),
            line_number: 0,
            header_parsed: false,
        }
    }

    /// Creates a VCF parser from a file path.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_path(path: impl AsRef<std::path::Path>) -> Result<VcfParser<std::fs::File>> {
        let file = std::fs::File::open(path)?;
        Ok(VcfParser::new(file))
    }

    /// Parses the VCF header.
    ///
    /// This must be called before iterating over records.
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::formats::vcf::VcfParser;
    /// use std::io::Cursor;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let data = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    /// let mut parser = VcfParser::new(Cursor::new(data.as_bytes()));
    /// let header = parser.parse_header()?;
    /// assert_eq!(header.fileformat, "VCFv4.2");
    /// # Ok(())
    /// # }
    /// ```
    pub fn parse_header(&mut self) -> Result<VcfHeader> {
        let mut header = VcfHeader::new("VCFv4.2".to_string());

        loop {
            self.line_buf.clear();
            let bytes_read = self.reader.read_line(&mut self.line_buf)?;

            if bytes_read == 0 {
                break; // EOF
            }

            self.line_number += 1;
            let line = self.line_buf.trim_end();

            // Skip empty lines
            if line.is_empty() {
                continue;
            }

            if line.starts_with("##") {
                // Metadata line
                header.parse_metadata_line(line);
            } else if line.starts_with("#CHROM") {
                // Column header
                header.parse_column_header(line);
                break; // End of header
            } else {
                // No column header found, this is a data line
                // We need to "put it back" but can't with BufRead
                // So we'll mark header as parsed and expect the iterator to handle it
                break;
            }
        }

        self.header_parsed = true;
        Ok(header)
    }
}

impl<R: Read> Iterator for VcfParser<R> {
    type Item = Result<VcfRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.header_parsed {
            // Header must be parsed first
            return Some(Err(crate::formats::primitives::FormatError::InvalidField {
                field: "header".to_string(),
                line: 0,
                reason: "Header must be parsed before reading records".to_string(),
            }));
        }

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

                    // Skip header/comment lines (in case header wasn't fully consumed)
                    if line.starts_with('#') {
                        continue;
                    }

                    return Some(VcfRecord::from_line(line));
                }
                Err(e) => return Some(Err(e.into())),
            }
        }
    }
}

/// VCF writer with header support.
///
/// Writes VCF files with metadata header and variant records.
///
/// # Examples
///
/// ```
/// use biometal::formats::vcf::{VcfHeader, VcfRecord, VcfWriter};
/// use biometal::formats::TabDelimitedRecord;
/// use std::path::PathBuf;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// // Create header
/// let mut header = VcfHeader::new("VCFv4.2".to_string());
/// header.info_fields.insert("DP".to_string(), "Total Depth".to_string());
///
/// // Create writer
/// # let temp_dir = tempfile::tempdir()?;
/// # let temp_path = temp_dir.path().join("output.vcf");
/// let mut writer = VcfWriter::create(&temp_path, header)?;
///
/// // Write records
/// let line = "chr1\t12345\trs123\tA\tT\t30.0\tPASS\tDP=100";
/// let record = VcfRecord::from_line(line)?;
/// writer.write_record(&record)?;
///
/// writer.finish()?;
/// # Ok(())
/// # }
/// ```
pub struct VcfWriter {
    writer: crate::io::compression::CompressedWriter,
    records_written: usize,
}

impl VcfWriter {
    /// Creates a VCF writer from a data sink with header.
    ///
    /// Writes the header immediately on creation.
    pub fn new(sink: crate::io::DataSink, header: VcfHeader) -> Result<Self> {
        use std::io::Write;

        let mut writer = crate::io::compression::CompressedWriter::new(sink)
            .map_err(|e| crate::formats::primitives::FormatError::from(e))?;

        // Write header lines
        for line in header.to_header_lines() {
            writeln!(writer, "{}", line)
                .map_err(|e| crate::formats::primitives::FormatError::from(e))?;
        }

        Ok(VcfWriter {
            writer,
            records_written: 0,
        })
    }

    /// Creates a VCF writer for the given path with header.
    ///
    /// Automatically detects compression from file extension (.gz, .bgz).
    /// Writes the header immediately on creation.
    pub fn create(path: impl AsRef<std::path::Path>, header: VcfHeader) -> Result<Self> {
        Self::new(crate::io::DataSink::from_path(path), header)
    }

    /// Creates a VCF writer for stdout with header.
    ///
    /// Writes the header immediately on creation.
    pub fn stdout(header: VcfHeader) -> Result<Self> {
        Self::new(crate::io::DataSink::stdout(), header)
    }

    /// Writes a VCF record.
    pub fn write_record(&mut self, record: &VcfRecord) -> Result<()> {
        use std::io::Write;

        let line = record.to_line();
        writeln!(self.writer, "{}", line)
            .map_err(|e| crate::formats::primitives::FormatError::from(e))?;
        self.records_written += 1;
        Ok(())
    }

    /// Returns the number of records written.
    pub fn records_written(&self) -> usize {
        self.records_written
    }

    /// Flushes and closes the writer.
    pub fn finish(mut self) -> Result<()> {
        use std::io::Write;

        self.writer.flush()
            .map_err(|e| crate::formats::primitives::FormatError::from(e))
    }
}

/// Parses VCF INFO field (semicolon-separated key=value pairs).
///
/// Flags (keys without values) are stored with empty string values.
fn parse_info_field(info_str: &str) -> HashMap<String, String> {
    let mut info = HashMap::new();

    if info_str == "." || info_str.is_empty() {
        return info;
    }

    for pair in info_str.split(';') {
        if let Some(eq_pos) = pair.find('=') {
            let key = pair[..eq_pos].to_string();
            let value = pair[eq_pos + 1..].to_string();
            info.insert(key, value);
        } else {
            // Flag field (no value)
            info.insert(pair.to_string(), String::new());
        }
    }

    info
}

/// Parses VCF structured field (e.g., ##INFO=<ID=DP,Number=1,Type=Integer,Description="...">).
///
/// Returns (ID, Description) tuple.
fn parse_vcf_structured_field(field_str: &str) -> Option<(String, String)> {
    // Extract content between < and >
    let start = field_str.find('<')?;
    let end = field_str.rfind('>')?;
    let content = &field_str[start + 1..end];

    // Parse key=value pairs
    let mut id = None;
    let mut description = String::new();

    for part in content.split(',') {
        if let Some(eq_pos) = part.find('=') {
            let key = &part[..eq_pos];
            let value = &part[eq_pos + 1..];

            match key {
                "ID" => id = Some(value.to_string()),
                "Description" => {
                    // Remove quotes
                    description = value.trim_matches('"').to_string();
                }
                "length" => {
                    // For contig
                    description = value.to_string();
                }
                _ => {}
            }
        }
    }

    id.map(|id| (id, description))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vcf_record_basic() {
        let line = "chr1\t12345\trs123\tA\tT\t30.0\tPASS\tDP=100";
        let record = VcfRecord::from_line(line).unwrap();

        assert_eq!(record.chrom, "chr1");
        assert_eq!(record.pos, 12345);
        assert_eq!(record.id, Some("rs123".to_string()));
        assert_eq!(record.reference, "A");
        assert_eq!(record.alternate, vec!["T"]);
        assert_eq!(record.quality, Some(30.0));
        assert_eq!(record.filter, Some("PASS".to_string()));
    }

    #[test]
    fn test_vcf_record_missing_values() {
        let line = "chr1\t100\t.\tA\tT\t.\t.\t.";
        let record = VcfRecord::from_line(line).unwrap();

        assert_eq!(record.id, None);
        assert_eq!(record.quality, None);
        assert_eq!(record.filter, None);
        assert!(record.info.is_empty());
    }

    #[test]
    fn test_vcf_multiple_alts() {
        let line = "chr1\t100\t.\tA\tT,G,C\t.\t.\t.";
        let record = VcfRecord::from_line(line).unwrap();

        assert_eq!(record.alternate.len(), 3);
        assert_eq!(record.alternate, vec!["T", "G", "C"]);
    }

    #[test]
    fn test_vcf_info_parsing() {
        let line = "chr1\t100\t.\tA\tT\t.\t.\tDP=50;AF=0.25;DB";
        let record = VcfRecord::from_line(line).unwrap();

        assert_eq!(record.info.get("DP"), Some(&"50".to_string()));
        assert_eq!(record.info.get("AF"), Some(&"0.25".to_string()));
        assert_eq!(record.info.get("DB"), Some(&"".to_string()));
    }

    #[test]
    fn test_vcf_round_trip() {
        let original = "chr1\t12345\trs123\tA\tT\t30\tPASS\tAF=0.5;DP=100";
        let record = VcfRecord::from_line(original).unwrap();
        let output = record.to_line();

        // Parse again to verify
        let record2 = VcfRecord::from_line(&output).unwrap();
        assert_eq!(record, record2);
    }

    #[test]
    fn test_vcf_header_fileformat() {
        let mut header = VcfHeader::new("VCFv4.2".to_string());
        header.parse_metadata_line("##fileformat=VCFv4.3");
        assert_eq!(header.fileformat, "VCFv4.3");
    }

    #[test]
    fn test_vcf_header_column() {
        let mut header = VcfHeader::new("VCFv4.2".to_string());
        header.parse_column_header("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2");

        assert_eq!(header.samples.len(), 2);
        assert_eq!(header.samples[0], "sample1");
        assert_eq!(header.samples[1], "sample2");
    }

    #[test]
    fn test_info_field_flags() {
        let info = parse_info_field("DP=100;DB;AF=0.5");
        assert_eq!(info.get("DP"), Some(&"100".to_string()));
        assert_eq!(info.get("DB"), Some(&"".to_string()));
        assert_eq!(info.get("AF"), Some(&"0.5".to_string()));
    }

    #[test]
    fn test_vcf_header_serialization() {
        let mut header = VcfHeader::new("VCFv4.2".to_string());
        header.info_fields.insert("DP".to_string(), "Total Depth".to_string());
        header.info_fields.insert("AF".to_string(), "Allele Frequency".to_string());
        header.filters.insert("PASS".to_string(), "All filters passed".to_string());
        header.samples = vec!["sample1".to_string(), "sample2".to_string()];

        let lines = header.to_header_lines();

        // Check first line is fileformat
        assert_eq!(lines[0], "##fileformat=VCFv4.2");

        // Check that we have INFO lines
        assert!(lines.iter().any(|l| l.starts_with("##INFO=<ID=AF,")));
        assert!(lines.iter().any(|l| l.starts_with("##INFO=<ID=DP,")));

        // Check that we have FILTER line
        assert!(lines.iter().any(|l| l.starts_with("##FILTER=<ID=PASS,")));

        // Check column header includes samples
        let col_header = lines.last().unwrap();
        assert!(col_header.starts_with("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"));
        assert!(col_header.contains("sample1"));
        assert!(col_header.contains("sample2"));
    }

    #[test]
    fn test_vcf_writer() -> Result<()> {
        use std::io::Cursor;

        // Create header
        let mut header = VcfHeader::new("VCFv4.2".to_string());
        header.info_fields.insert("DP".to_string(), "Total Depth".to_string());

        // Create temp file
        let temp_dir = tempfile::tempdir()?;
        let temp_path = temp_dir.path().join("test.vcf");

        // Write VCF
        let mut writer = VcfWriter::create(&temp_path, header)?;

        let line1 = "chr1\t12345\trs123\tA\tT\t30.0\tPASS\tDP=100";
        let record1 = VcfRecord::from_line(line1)?;
        writer.write_record(&record1)?;

        let line2 = "chr1\t12346\trs124\tG\tC\t40.0\tPASS\tDP=150";
        let record2 = VcfRecord::from_line(line2)?;
        writer.write_record(&record2)?;

        assert_eq!(writer.records_written(), 2);
        writer.finish()?;

        // Read back and verify
        let content = std::fs::read_to_string(&temp_path)?;
        let lines: Vec<&str> = content.lines().collect();

        // Should have header lines + 2 records
        assert!(lines.len() >= 4);

        // First line should be fileformat
        assert!(lines[0].starts_with("##fileformat="));

        // Should have INFO line
        assert!(lines.iter().any(|l| l.starts_with("##INFO=<ID=DP,")));

        // Should have column header
        assert!(lines.iter().any(|l| l.starts_with("#CHROM\tPOS")));

        // Should have both records
        assert!(content.contains("chr1\t12345\trs123"));
        assert!(content.contains("chr1\t12346\trs124"));

        Ok(())
    }

    #[test]
    fn test_vcf_writer_round_trip() -> Result<()> {
        // Create header with samples
        let mut header = VcfHeader::new("VCFv4.2".to_string());
        header.info_fields.insert("DP".to_string(), "Total Depth".to_string());
        header.info_fields.insert("AF".to_string(), "Allele Frequency".to_string());
        header.samples = vec!["sample1".to_string()];

        // Create temp file
        let temp_dir = tempfile::tempdir()?;
        let temp_path = temp_dir.path().join("test_roundtrip.vcf");

        // Write VCF
        let mut writer = VcfWriter::create(&temp_path, header.clone())?;

        let line1 = "chr1\t100\t.\tA\tT\t.\t.\tDP=50;AF=0.5";
        let record1 = VcfRecord::from_line(line1)?;
        writer.write_record(&record1)?;

        writer.finish()?;

        // Read back with parser
        let file = std::fs::File::open(&temp_path)?;
        let mut parser = VcfParser::new(file);
        let read_header = parser.parse_header()?;

        // Verify header
        assert_eq!(read_header.fileformat, "VCFv4.2");
        assert_eq!(read_header.info_fields.get("DP"), Some(&"Total Depth".to_string()));
        assert_eq!(read_header.samples, vec!["sample1"]);

        // Verify record
        let records: Vec<_> = parser.collect::<Result<Vec<_>>>()?;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[0].pos, 100);

        Ok(())
    }
}
