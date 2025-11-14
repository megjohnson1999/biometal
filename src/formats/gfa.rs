//! GFA (Graphical Fragment Assembly) format parser.
//!
//! GFA is a text format for representing assembly graphs and pangenomes:
//! - **Segments (S)**: Graph nodes (sequences)
//! - **Links (L)**: Edges between segments
//! - **Paths (P)**: Ordered traversals through the graph
//! - **Walks (W)**: GFA v2 path representation
//!
//! # Format Specification
//!
//! GFA uses tab-delimited records with a single-character record type:
//! - **H**: Header
//! - **S**: Segment (node with sequence)
//! - **L**: Link (edge between segments)
//! - **P**: Path (ordered walk through graph)
//! - **W**: Walk (GFA v2 alternative to Path)
//! - **C**: Containment (one segment contained in another)
//!
//! # Examples
//!
//! ## Segment Records
//!
//! ```
//! use biometal::formats::gfa::{GfaSegment, GfaRecord};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let line = "S\tctg1\tACGTACGT\tLN:i:8\tRC:i:42";
//! let record = GfaRecord::from_line(line)?;
//!
//! if let GfaRecord::Segment(seg) = record {
//!     assert_eq!(seg.name, "ctg1");
//!     assert_eq!(seg.sequence, "ACGTACGT");
//!     assert_eq!(seg.length(), 8);
//!     assert_eq!(seg.tags.get("LN"), Some(&"i:8".to_string()));
//! }
//! # Ok(())
//! # }
//! ```
//!
//! ## Link Records
//!
//! ```
//! use biometal::formats::gfa::{GfaLink, GfaRecord, Orientation};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let line = "L\tctg1\t+\tctg2\t-\t4M";
//! let record = GfaRecord::from_line(line)?;
//!
//! if let GfaRecord::Link(link) = record {
//!     assert_eq!(link.from_segment, "ctg1");
//!     assert_eq!(link.from_orient, Orientation::Forward);
//!     assert_eq!(link.to_segment, "ctg2");
//!     assert_eq!(link.to_orient, Orientation::Reverse);
//!     assert_eq!(link.overlap, "4M");
//! }
//! # Ok(())
//! # }
//! ```
//!
//! ## Path Records
//!
//! ```
//! use biometal::formats::gfa::{GfaPath, GfaRecord};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let line = "P\tpath1\tctg1+,ctg2-,ctg3+\t4M,5M";
//! let record = GfaRecord::from_line(line)?;
//!
//! if let GfaRecord::Path(path) = record {
//!     assert_eq!(path.name, "path1");
//!     assert_eq!(path.segments.len(), 3);
//!     assert_eq!(path.segments[0], "ctg1+");
//!     assert_eq!(path.overlaps.len(), 2);
//! }
//! # Ok(())
//! # }
//! ```
//!
//! ## Streaming parser
//!
//! ```no_run
//! use biometal::formats::gfa::{GfaRecord, GfaParser};
//! use std::fs::File;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let file = File::open("assembly.gfa")?;
//! let parser = GfaParser::new(file);
//!
//! for result in parser {
//!     let record = result?;
//!     match record {
//!         GfaRecord::Segment(seg) => println!("Segment: {} (len={})", seg.name, seg.length()),
//!         GfaRecord::Link(link) => println!("Link: {} -> {}", link.from_segment, link.to_segment),
//!         GfaRecord::Path(path) => println!("Path: {} ({} segments)", path.name, path.segments.len()),
//!         _ => {}
//!     }
//! }
//! # Ok(())
//! # }
//! ```

use crate::formats::primitives::{fields::split_fields, Result, TabDelimitedRecord};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};
use std::str::FromStr;

/// Segment/edge orientation in GFA.
///
/// - `+`: Forward orientation
/// - `-`: Reverse orientation (reverse complement)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Orientation {
    /// Forward orientation (+)
    Forward,
    /// Reverse orientation (-)
    Reverse,
}

impl FromStr for Orientation {
    type Err = crate::formats::primitives::FormatError;

    fn from_str(s: &str) -> Result<Self> {
        match s {
            "+" => Ok(Orientation::Forward),
            "-" => Ok(Orientation::Reverse),
            _ => Err(crate::formats::primitives::FormatError::InvalidField {
                field: "orientation".to_string(),
                line: 0,
                reason: format!("Invalid orientation: {} (expected '+' or '-')", s),
            }),
        }
    }
}

impl std::fmt::Display for Orientation {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Orientation::Forward => write!(f, "+"),
            Orientation::Reverse => write!(f, "-"),
        }
    }
}

/// GFA Segment record (S line).
///
/// Represents a node in the assembly graph with an associated sequence.
///
/// # Format
///
/// ```text
/// S  name    sequence   [tags]
/// S  ctg1    ACGTACGT   LN:i:8  RC:i:42
/// ```
///
/// # Examples
///
/// ```
/// use biometal::formats::gfa::GfaSegment;
/// use biometal::formats::TabDelimitedRecord;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let line = "S\tctg1\tACGTACGT\tLN:i:8";
/// let seg = GfaSegment::from_line(line)?;
///
/// assert_eq!(seg.name, "ctg1");
/// assert_eq!(seg.sequence, "ACGTACGT");
/// assert_eq!(seg.length(), 8);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GfaSegment {
    /// Segment name/ID
    pub name: String,
    /// Segment sequence (can be * for absent)
    pub sequence: String,
    /// Optional tags (e.g., LN:i:8, RC:i:42)
    pub tags: HashMap<String, String>,
}

impl GfaSegment {
    /// Returns the length of the segment sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::formats::gfa::GfaSegment;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let seg = GfaSegment::from_line("S\tctg1\tACGTACGT")?;
    /// assert_eq!(seg.length(), 8);
    /// # Ok(())
    /// # }
    /// ```
    pub fn length(&self) -> usize {
        if self.sequence == "*" {
            // Try to get length from LN tag (format: "i:100")
            self.tags
                .get("LN")
                .and_then(|s| s.split(':').last())
                .and_then(|s| s.parse::<usize>().ok())
                .unwrap_or(0)
        } else {
            self.sequence.len()
        }
    }
}

impl TabDelimitedRecord for GfaSegment {
    fn from_line(line: &str) -> Result<Self> {
        let fields = split_fields(line, Some(3), 0)?;

        if fields[0] != "S" {
            return Err(crate::formats::primitives::FormatError::InvalidField {
                field: "record_type".to_string(),
                line: 0,
                reason: format!("Expected 'S', got '{}'", fields[0]),
            });
        }

        let name = fields[1].to_string();
        let sequence = fields[2].to_string();

        // Parse optional tags
        let tags = parse_gfa_tags(&fields[3..]);

        Ok(GfaSegment {
            name,
            sequence,
            tags,
        })
    }

    fn to_line(&self) -> String {
        let mut line = format!("S\t{}\t{}", self.name, self.sequence);

        // Append tags in sorted order for consistency
        let mut tag_keys: Vec<_> = self.tags.keys().collect();
        tag_keys.sort();
        for key in tag_keys {
            line.push('\t');
            line.push_str(key);
            line.push(':');
            line.push_str(&self.tags[key]);
        }

        line
    }

    fn expected_fields() -> Option<usize> {
        Some(3) // At least: type, name, sequence
    }
}

/// GFA Link record (L line).
///
/// Represents an edge between two segments in the assembly graph.
///
/// # Format
///
/// ```text
/// L  from_seg  from_orient  to_seg  to_orient  overlap  [tags]
/// L  ctg1      +            ctg2    -          4M
/// ```
///
/// # Examples
///
/// ```
/// use biometal::formats::gfa::{GfaLink, Orientation};
/// use biometal::formats::TabDelimitedRecord;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let line = "L\tctg1\t+\tctg2\t-\t4M";
/// let link = GfaLink::from_line(line)?;
///
/// assert_eq!(link.from_segment, "ctg1");
/// assert_eq!(link.from_orient, Orientation::Forward);
/// assert_eq!(link.to_segment, "ctg2");
/// assert_eq!(link.to_orient, Orientation::Reverse);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GfaLink {
    /// Source segment name
    pub from_segment: String,
    /// Source segment orientation
    pub from_orient: Orientation,
    /// Target segment name
    pub to_segment: String,
    /// Target segment orientation
    pub to_orient: Orientation,
    /// Overlap between segments (CIGAR format or *)
    pub overlap: String,
    /// Optional tags
    pub tags: HashMap<String, String>,
}

impl TabDelimitedRecord for GfaLink {
    fn from_line(line: &str) -> Result<Self> {
        let fields = split_fields(line, Some(6), 0)?;

        if fields[0] != "L" {
            return Err(crate::formats::primitives::FormatError::InvalidField {
                field: "record_type".to_string(),
                line: 0,
                reason: format!("Expected 'L', got '{}'", fields[0]),
            });
        }

        let from_segment = fields[1].to_string();
        let from_orient = Orientation::from_str(fields[2])?;
        let to_segment = fields[3].to_string();
        let to_orient = Orientation::from_str(fields[4])?;
        let overlap = fields[5].to_string();

        let tags = parse_gfa_tags(&fields[6..]);

        Ok(GfaLink {
            from_segment,
            from_orient,
            to_segment,
            to_orient,
            overlap,
            tags,
        })
    }

    fn to_line(&self) -> String {
        let mut line = format!(
            "L\t{}\t{}\t{}\t{}\t{}",
            self.from_segment,
            self.from_orient,
            self.to_segment,
            self.to_orient,
            self.overlap
        );

        let mut tag_keys: Vec<_> = self.tags.keys().collect();
        tag_keys.sort();
        for key in tag_keys {
            line.push('\t');
            line.push_str(key);
            line.push(':');
            line.push_str(&self.tags[key]);
        }

        line
    }

    fn expected_fields() -> Option<usize> {
        Some(6)
    }
}

/// GFA Path record (P line).
///
/// Represents an ordered traversal through the assembly graph.
///
/// # Format
///
/// ```text
/// P  path_name  segment_names  overlaps  [tags]
/// P  path1      ctg1+,ctg2-    4M,5M
/// ```
///
/// # Examples
///
/// ```
/// use biometal::formats::gfa::GfaPath;
/// use biometal::formats::TabDelimitedRecord;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let line = "P\tpath1\tctg1+,ctg2-,ctg3+\t4M,5M";
/// let path = GfaPath::from_line(line)?;
///
/// assert_eq!(path.name, "path1");
/// assert_eq!(path.segments.len(), 3);
/// assert_eq!(path.overlaps.len(), 2);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GfaPath {
    /// Path name/ID
    pub name: String,
    /// Ordered list of segment names with orientations (e.g., "ctg1+")
    pub segments: Vec<String>,
    /// Overlaps between consecutive segments (CIGAR format)
    pub overlaps: Vec<String>,
    /// Optional tags
    pub tags: HashMap<String, String>,
}

impl TabDelimitedRecord for GfaPath {
    fn from_line(line: &str) -> Result<Self> {
        let fields = split_fields(line, Some(4), 0)?;

        if fields[0] != "P" {
            return Err(crate::formats::primitives::FormatError::InvalidField {
                field: "record_type".to_string(),
                line: 0,
                reason: format!("Expected 'P', got '{}'", fields[0]),
            });
        }

        let name = fields[1].to_string();

        // Parse comma-separated segment names
        let segments: Vec<String> = fields[2]
            .split(',')
            .map(|s| s.trim().to_string())
            .collect();

        // Parse comma-separated overlaps
        let overlaps: Vec<String> = if fields[3] == "*" || fields[3].is_empty() {
            vec![]
        } else {
            fields[3]
                .split(',')
                .map(|s| s.trim().to_string())
                .collect()
        };

        let tags = parse_gfa_tags(&fields[4..]);

        Ok(GfaPath {
            name,
            segments,
            overlaps,
            tags,
        })
    }

    fn to_line(&self) -> String {
        let segments_str = self.segments.join(",");
        let overlaps_str = if self.overlaps.is_empty() {
            "*".to_string()
        } else {
            self.overlaps.join(",")
        };

        let mut line = format!("P\t{}\t{}\t{}", self.name, segments_str, overlaps_str);

        let mut tag_keys: Vec<_> = self.tags.keys().collect();
        tag_keys.sort();
        for key in tag_keys {
            line.push('\t');
            line.push_str(key);
            line.push(':');
            line.push_str(&self.tags[key]);
        }

        line
    }

    fn expected_fields() -> Option<usize> {
        Some(4)
    }
}

/// GFA Header record (H line).
///
/// Contains metadata about the GFA file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GfaHeader {
    /// Header tags (e.g., VN:Z:1.0)
    pub tags: HashMap<String, String>,
}

impl TabDelimitedRecord for GfaHeader {
    fn from_line(line: &str) -> Result<Self> {
        let fields = split_fields(line, Some(1), 0)?;

        if fields[0] != "H" {
            return Err(crate::formats::primitives::FormatError::InvalidField {
                field: "record_type".to_string(),
                line: 0,
                reason: format!("Expected 'H', got '{}'", fields[0]),
            });
        }

        let tags = parse_gfa_tags(&fields[1..]);

        Ok(GfaHeader { tags })
    }

    fn to_line(&self) -> String {
        let mut line = "H".to_string();

        let mut tag_keys: Vec<_> = self.tags.keys().collect();
        tag_keys.sort();
        for key in tag_keys {
            line.push('\t');
            line.push_str(key);
            line.push(':');
            line.push_str(&self.tags[key]);
        }

        line
    }

    fn expected_fields() -> Option<usize> {
        Some(1)
    }
}

/// Unified GFA record type.
///
/// Represents any type of GFA record (Header, Segment, Link, Path, etc.).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GfaRecord {
    /// Header record (H)
    Header(GfaHeader),
    /// Segment record (S)
    Segment(GfaSegment),
    /// Link record (L)
    Link(GfaLink),
    /// Path record (P)
    Path(GfaPath),
    /// Comment or unrecognized record
    Comment(String),
}

impl GfaRecord {
    /// Parses a GFA record from a line.
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::formats::gfa::GfaRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let record = GfaRecord::from_line("S\tctg1\tACGT")?;
    /// assert!(matches!(record, GfaRecord::Segment(_)));
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_line(line: &str) -> Result<Self> {
        let line = line.trim();

        // Skip empty lines and comments
        if line.is_empty() || line.starts_with('#') {
            return Ok(GfaRecord::Comment(line.to_string()));
        }

        // Determine record type from first character
        let record_type = line.chars().next().unwrap();

        match record_type {
            'H' => Ok(GfaRecord::Header(GfaHeader::from_line(line)?)),
            'S' => Ok(GfaRecord::Segment(GfaSegment::from_line(line)?)),
            'L' => Ok(GfaRecord::Link(GfaLink::from_line(line)?)),
            'P' => Ok(GfaRecord::Path(GfaPath::from_line(line)?)),
            _ => Ok(GfaRecord::Comment(line.to_string())),
        }
    }

    /// Converts the record back to a GFA line.
    pub fn to_line(&self) -> String {
        match self {
            GfaRecord::Header(h) => h.to_line(),
            GfaRecord::Segment(s) => s.to_line(),
            GfaRecord::Link(l) => l.to_line(),
            GfaRecord::Path(p) => p.to_line(),
            GfaRecord::Comment(c) => c.clone(),
        }
    }
}

/// GFA parser with streaming support.
///
/// Parses GFA files record by record with constant memory usage.
pub struct GfaParser<R: Read> {
    reader: BufReader<R>,
    line_buf: String,
    line_number: usize,
}

impl<R: Read> GfaParser<R> {
    /// Creates a new GFA parser from a reader.
    pub fn new(reader: R) -> Self {
        GfaParser {
            reader: BufReader::new(reader),
            line_buf: String::new(),
            line_number: 0,
        }
    }

    /// Creates a GFA parser from a file path.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_path(path: impl AsRef<std::path::Path>) -> Result<GfaParser<std::fs::File>> {
        let file = std::fs::File::open(path)?;
        Ok(GfaParser::new(file))
    }
}

impl<R: Read> Iterator for GfaParser<R> {
    type Item = Result<GfaRecord>;

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

                    // Skip comment lines (but don't skip record types we don't handle)
                    if line.starts_with('#') {
                        continue;
                    }

                    return Some(GfaRecord::from_line(line));
                }
                Err(e) => return Some(Err(e.into())),
            }
        }
    }
}

/// Parses GFA optional tags.
///
/// GFA tags have format: `TAG:TYPE:VALUE`
/// - TAG: Two-character tag name
/// - TYPE: Single-character type (A, i, f, Z, J, H, B)
/// - VALUE: The value
///
/// This simplified version stores `TYPE:VALUE` for each tag.
fn parse_gfa_tags(fields: &[&str]) -> HashMap<String, String> {
    let mut tags = HashMap::new();

    for field in fields {
        // GFA tags: TAG:TYPE:VALUE
        let parts: Vec<&str> = field.splitn(3, ':').collect();
        if parts.len() >= 2 {
            let tag_name = parts[0].to_string();
            let tag_value = if parts.len() == 3 {
                format!("{}:{}", parts[1], parts[2])
            } else {
                parts[1].to_string()
            };
            tags.insert(tag_name, tag_value);
        }
    }

    tags
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_orientation_parse() {
        assert_eq!(Orientation::from_str("+").unwrap(), Orientation::Forward);
        assert_eq!(Orientation::from_str("-").unwrap(), Orientation::Reverse);
        assert!(Orientation::from_str("x").is_err());
    }

    #[test]
    fn test_segment_basic() {
        let line = "S\tctg1\tACGTACGT";
        let seg = GfaSegment::from_line(line).unwrap();

        assert_eq!(seg.name, "ctg1");
        assert_eq!(seg.sequence, "ACGTACGT");
        assert_eq!(seg.length(), 8);
    }

    #[test]
    fn test_segment_with_tags() {
        let line = "S\tctg1\tACGT\tLN:i:4\tRC:i:10";
        let seg = GfaSegment::from_line(line).unwrap();

        assert_eq!(seg.tags.get("LN"), Some(&"i:4".to_string()));
        assert_eq!(seg.tags.get("RC"), Some(&"i:10".to_string()));
    }

    #[test]
    fn test_segment_absent_sequence() {
        let line = "S\tctg1\t*\tLN:i:100";
        let seg = GfaSegment::from_line(line).unwrap();

        assert_eq!(seg.sequence, "*");
        assert_eq!(seg.length(), 100); // From LN tag
    }

    #[test]
    fn test_link_basic() {
        let line = "L\tctg1\t+\tctg2\t-\t4M";
        let link = GfaLink::from_line(line).unwrap();

        assert_eq!(link.from_segment, "ctg1");
        assert_eq!(link.from_orient, Orientation::Forward);
        assert_eq!(link.to_segment, "ctg2");
        assert_eq!(link.to_orient, Orientation::Reverse);
        assert_eq!(link.overlap, "4M");
    }

    #[test]
    fn test_path_basic() {
        let line = "P\tpath1\tctg1+,ctg2-,ctg3+\t4M,5M";
        let path = GfaPath::from_line(line).unwrap();

        assert_eq!(path.name, "path1");
        assert_eq!(path.segments, vec!["ctg1+", "ctg2-", "ctg3+"]);
        assert_eq!(path.overlaps, vec!["4M", "5M"]);
    }

    #[test]
    fn test_header() {
        let line = "H\tVN:Z:1.0";
        let header = GfaHeader::from_line(line).unwrap();

        assert_eq!(header.tags.get("VN"), Some(&"Z:1.0".to_string()));
    }

    #[test]
    fn test_gfa_record_dispatch() {
        let lines = vec![
            "H\tVN:Z:1.0",
            "S\tctg1\tACGT",
            "L\tctg1\t+\tctg2\t-\t4M",
            "P\tpath1\tctg1+,ctg2-\t4M",
        ];

        for line in lines {
            let record = GfaRecord::from_line(line).unwrap();
            match record {
                GfaRecord::Header(_) => assert!(line.starts_with('H')),
                GfaRecord::Segment(_) => assert!(line.starts_with('S')),
                GfaRecord::Link(_) => assert!(line.starts_with('L')),
                GfaRecord::Path(_) => assert!(line.starts_with('P')),
                _ => panic!("Unexpected record type"),
            }
        }
    }

    #[test]
    fn test_round_trip_segment() {
        let original = "S\tctg1\tACGT\tLN:i:4";
        let seg = GfaSegment::from_line(original).unwrap();
        let output = seg.to_line();
        assert_eq!(output, original);
    }

    #[test]
    fn test_round_trip_link() {
        let original = "L\tctg1\t+\tctg2\t-\t4M";
        let link = GfaLink::from_line(original).unwrap();
        let output = link.to_line();
        assert_eq!(output, original);
    }

    #[test]
    fn test_round_trip_path() {
        let original = "P\tpath1\tctg1+,ctg2-\t4M";
        let path = GfaPath::from_line(original).unwrap();
        let output = path.to_line();
        assert_eq!(output, original);
    }
}
