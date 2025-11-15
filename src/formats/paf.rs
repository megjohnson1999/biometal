//! PAF (Pairwise mApping Format) parser.
//!
//! PAF is the default output format of minimap2, used for describing approximate
//! mapping positions between two sets of sequences (e.g., reads to reference,
//! assembly to assembly).
//!
//! # Format Specification
//!
//! - **Columns**: 12 required tab-delimited fields + optional SAM-like tags
//! - **Coordinates**: 0-based, half-open `[start, end)` (BED-like)
//! - **Strand**: `+` (same strand) or `-` (opposite strand)
//! - **Use cases**: Long-read alignment, assembly comparison, read mapping
//!
//! # Examples
//!
//! ## Basic PAF record
//!
//! ```
//! use biometal::formats::paf::PafRecord;
//! use biometal::formats::TabDelimitedRecord;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let line = "read1\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t60";
//! let record = PafRecord::from_line(line)?;
//!
//! assert_eq!(record.query_name, "read1");
//! assert_eq!(record.query_length, 10000);
//! assert_eq!(record.target_name, "chr1");
//! assert_eq!(record.strand, '+');
//! assert_eq!(record.mapq, 60);
//!
//! // Calculate identity
//! let identity = record.identity();
//! assert!((identity - 0.9693877).abs() < 0.0001); // 9500/9800
//! # Ok(())
//! # }
//! ```
//!
//! ## Streaming parser
//!
//! ```no_run
//! use biometal::formats::paf::PafParser;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let parser = PafParser::from_path("alignments.paf")?;
//!
//! for result in parser {
//!     let record = result?;
//!     if record.mapq >= 20 && record.identity() >= 0.95 {
//!         println!("{}: {}% identity", record.query_name, record.identity() * 100.0);
//!     }
//! }
//! # Ok(())
//! # }
//! ```

use crate::formats::primitives::{
    fields::{parse_required, split_fields},
    Result, TabDelimitedRecord,
};
use std::io::{BufRead, BufReader, Read};

/// PAF alignment record.
///
/// Represents a pairwise mapping between query and target sequences.
///
/// # Format
///
/// ```text
/// query_name  query_len  q_start  q_end  strand  target_name  target_len  t_start  t_end  matches  aln_len  mapq
/// read1       10000      100      9900   +       chr1         50000000    1000     10900  9500     9800     60
/// ```
///
/// # Coordinates
///
/// - **Query**: 0-based, half-open `[q_start, q_end)`
/// - **Target**: 0-based, half-open `[t_start, t_end)`
/// - **Strand**: `+` or `-` (query/target orientation)
///
/// # Examples
///
/// ```
/// use biometal::formats::paf::PafRecord;
/// use biometal::formats::TabDelimitedRecord;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let line = "read1\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t60";
/// let record = PafRecord::from_line(line)?;
///
/// // Alignment coverage
/// let query_coverage = record.query_coverage();
/// assert_eq!(query_coverage, 0.98); // 9800/10000
///
/// // Alignment identity
/// let identity = record.identity();
/// assert!((identity - 0.9693877).abs() < 0.0001);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PafRecord {
    /// Query sequence name
    pub query_name: String,
    /// Query sequence length
    pub query_length: u64,
    /// Query start coordinate (0-based)
    pub query_start: u64,
    /// Query end coordinate (0-based, exclusive)
    pub query_end: u64,
    /// Relative strand: '+' (same) or '-' (opposite)
    pub strand: char,
    /// Target sequence name
    pub target_name: String,
    /// Target sequence length
    pub target_length: u64,
    /// Target start on original strand (0-based)
    pub target_start: u64,
    /// Target end on original strand (0-based, exclusive)
    pub target_end: u64,
    /// Number of matching bases
    pub num_matches: u64,
    /// Alignment block length (matches + mismatches + gaps)
    pub alignment_length: u64,
    /// Mapping quality (0-255, 255 = missing)
    pub mapq: u8,
}

impl PafRecord {
    /// Calculate alignment identity (matches / alignment_length).
    ///
    /// Returns a value between 0.0 and 1.0.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::formats::paf::PafRecord;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let record = PafRecord::from_line("read1\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t60")?;
    /// let identity = record.identity();
    /// assert!((identity - 0.9693877).abs() < 0.0001); // 9500/9800
    /// # Ok(())
    /// # }
    /// ```
    pub fn identity(&self) -> f64 {
        if self.alignment_length == 0 {
            0.0
        } else {
            self.num_matches as f64 / self.alignment_length as f64
        }
    }

    /// Calculate query coverage (aligned bases / query length).
    ///
    /// Returns a value between 0.0 and 1.0.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::formats::paf::PafRecord;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let record = PafRecord::from_line("read1\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t60")?;
    /// assert_eq!(record.query_coverage(), 0.98); // 9800/10000
    /// # Ok(())
    /// # }
    /// ```
    pub fn query_coverage(&self) -> f64 {
        if self.query_length == 0 {
            0.0
        } else {
            (self.query_end - self.query_start) as f64 / self.query_length as f64
        }
    }

    /// Calculate target coverage (aligned bases / target length).
    ///
    /// Returns a value between 0.0 and 1.0.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::formats::paf::PafRecord;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let record = PafRecord::from_line("read1\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t60")?;
    /// let target_cov = record.target_coverage();
    /// assert!((target_cov - 0.000198).abs() < 0.000001); // 9900/50000000
    /// # Ok(())
    /// # }
    /// ```
    pub fn target_coverage(&self) -> f64 {
        if self.target_length == 0 {
            0.0
        } else {
            (self.target_end - self.target_start) as f64 / self.target_length as f64
        }
    }

    /// Check if the mapping quality passes a threshold.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::formats::paf::PafRecord;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let record = PafRecord::from_line("read1\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t60")?;
    /// assert!(record.is_high_quality(20)); // mapq=60 > 20
    /// assert!(!record.is_high_quality(255)); // 255 means missing
    /// # Ok(())
    /// # }
    /// ```
    pub fn is_high_quality(&self, min_mapq: u8) -> bool {
        self.mapq >= min_mapq && self.mapq < 255
    }

    /// Check if alignment is on the forward strand.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::formats::paf::PafRecord;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let fwd = PafRecord::from_line("read1\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t60")?;
    /// let rev = PafRecord::from_line("read2\t10000\t100\t9900\t-\tchr1\t50000000\t1000\t10900\t9500\t9800\t60")?;
    ///
    /// assert!(fwd.is_forward());
    /// assert!(!rev.is_forward());
    /// # Ok(())
    /// # }
    /// ```
    pub fn is_forward(&self) -> bool {
        self.strand == '+'
    }

    /// Get the alignment length on the query.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::formats::paf::PafRecord;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let record = PafRecord::from_line("read1\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t60")?;
    /// assert_eq!(record.query_aligned_length(), 9800);
    /// # Ok(())
    /// # }
    /// ```
    pub fn query_aligned_length(&self) -> u64 {
        self.query_end - self.query_start
    }

    /// Get the alignment length on the target.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::formats::paf::PafRecord;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let record = PafRecord::from_line("read1\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t60")?;
    /// assert_eq!(record.target_aligned_length(), 9900);
    /// # Ok(())
    /// # }
    /// ```
    pub fn target_aligned_length(&self) -> u64 {
        self.target_end - self.target_start
    }
}

impl TabDelimitedRecord for PafRecord {
    fn from_line(line: &str) -> Result<Self> {
        // PAF has 12 required fields + optional SAM-like tags
        // We only parse the 12 required fields
        let fields = split_fields(line, Some(12), 0)?;

        let query_name = fields[0].to_string();
        let query_length: u64 = parse_required(fields[1], "query_length", 0)?;
        let query_start: u64 = parse_required(fields[2], "query_start", 0)?;
        let query_end: u64 = parse_required(fields[3], "query_end", 0)?;

        let strand_str = fields[4];
        let strand = if strand_str == "+" {
            '+'
        } else if strand_str == "-" {
            '-'
        } else {
            return Err(crate::formats::primitives::FormatError::InvalidField {
                field: "strand".to_string(),
                line: 0,
                reason: format!("Invalid strand '{}', expected '+' or '-'", strand_str),
            });
        };

        let target_name = fields[5].to_string();
        let target_length: u64 = parse_required(fields[6], "target_length", 0)?;
        let target_start: u64 = parse_required(fields[7], "target_start", 0)?;
        let target_end: u64 = parse_required(fields[8], "target_end", 0)?;
        let num_matches: u64 = parse_required(fields[9], "num_matches", 0)?;
        let alignment_length: u64 = parse_required(fields[10], "alignment_length", 0)?;
        let mapq: u8 = parse_required(fields[11], "mapq", 0)?;

        Ok(PafRecord {
            query_name,
            query_length,
            query_start,
            query_end,
            strand,
            target_name,
            target_length,
            target_start,
            target_end,
            num_matches,
            alignment_length,
            mapq,
        })
    }

    fn to_line(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_name,
            self.query_length,
            self.query_start,
            self.query_end,
            self.strand,
            self.target_name,
            self.target_length,
            self.target_start,
            self.target_end,
            self.num_matches,
            self.alignment_length,
            self.mapq
        )
    }

    fn expected_fields() -> Option<usize> {
        Some(12)
    }
}

/// PAF parser with streaming support.
///
/// Iterates through PAF alignment records one at a time for constant memory usage.
///
/// # Examples
///
/// ```no_run
/// use biometal::formats::paf::PafParser;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let parser = PafParser::from_path("alignments.paf")?;
///
/// let mut high_quality = 0;
/// for result in parser {
///     let record = result?;
///     if record.mapq >= 20 && record.identity() >= 0.95 {
///         high_quality += 1;
///     }
/// }
/// println!("High-quality alignments: {}", high_quality);
/// # Ok(())
/// # }
/// ```
pub struct PafParser<R: BufRead> {
    reader: R,
}

impl PafParser<BufReader<std::fs::File>> {
    /// Create a parser from a file path.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::paf::PafParser;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let parser = PafParser::from_path("alignments.paf")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_path(path: impl AsRef<std::path::Path>) -> std::io::Result<Self> {
        let file = std::fs::File::open(path)?;
        let reader = BufReader::new(file);
        Ok(PafParser { reader })
    }
}

impl<R: BufRead> PafParser<R> {
    /// Create a parser from a buffered reader.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::formats::paf::PafParser;
    /// use std::io::Cursor;
    ///
    /// let data = "read1\t100\t0\t100\t+\tchr1\t1000\t0\t100\t95\t100\t60\n";
    /// let cursor = Cursor::new(data);
    /// let parser = PafParser::new(cursor);
    /// ```
    pub fn new(reader: R) -> Self {
        PafParser { reader }
    }
}

impl<R: BufRead> Iterator for PafParser<R> {
    type Item = Result<PafRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line = String::new();

        loop {
            line.clear();
            match self.reader.read_line(&mut line) {
                Ok(0) => return None, // EOF
                Ok(_) => {
                    let trimmed = line.trim();
                    // Skip empty lines
                    if trimmed.is_empty() {
                        continue;
                    }
                    return Some(PafRecord::from_line(trimmed));
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
    fn test_paf_basic() {
        let line = "read1\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t60";
        let record = PafRecord::from_line(line).unwrap();

        assert_eq!(record.query_name, "read1");
        assert_eq!(record.query_length, 10000);
        assert_eq!(record.query_start, 100);
        assert_eq!(record.query_end, 9900);
        assert_eq!(record.strand, '+');
        assert_eq!(record.target_name, "chr1");
        assert_eq!(record.target_length, 50000000);
        assert_eq!(record.target_start, 1000);
        assert_eq!(record.target_end, 10900);
        assert_eq!(record.num_matches, 9500);
        assert_eq!(record.alignment_length, 9800);
        assert_eq!(record.mapq, 60);
    }

    #[test]
    fn test_paf_identity() {
        let line = "read1\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t60";
        let record = PafRecord::from_line(line).unwrap();

        let identity = record.identity();
        assert!((identity - 0.9693877).abs() < 0.0001);
    }

    #[test]
    fn test_paf_coverage() {
        let line = "read1\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t60";
        let record = PafRecord::from_line(line).unwrap();

        assert_eq!(record.query_coverage(), 0.98);
        let target_cov = record.target_coverage();
        assert!((target_cov - 0.000198).abs() < 0.000001);
    }

    #[test]
    fn test_paf_strand() {
        let fwd = "read1\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t60";
        let rev = "read2\t10000\t100\t9900\t-\tchr1\t50000000\t1000\t10900\t9500\t9800\t60";

        let record1 = PafRecord::from_line(fwd).unwrap();
        let record2 = PafRecord::from_line(rev).unwrap();

        assert!(record1.is_forward());
        assert!(!record2.is_forward());
        assert_eq!(record1.strand, '+');
        assert_eq!(record2.strand, '-');
    }

    #[test]
    fn test_paf_quality() {
        let high = "read1\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t60";
        let low = "read2\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t10";
        let missing = "read3\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t255";

        let record1 = PafRecord::from_line(high).unwrap();
        let record2 = PafRecord::from_line(low).unwrap();
        let record3 = PafRecord::from_line(missing).unwrap();

        assert!(record1.is_high_quality(20));
        assert!(!record2.is_high_quality(20));
        assert!(!record3.is_high_quality(0)); // 255 is treated as missing
    }

    #[test]
    fn test_paf_aligned_lengths() {
        let line = "read1\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t60";
        let record = PafRecord::from_line(line).unwrap();

        assert_eq!(record.query_aligned_length(), 9800);
        assert_eq!(record.target_aligned_length(), 9900);
    }

    #[test]
    fn test_paf_round_trip() {
        let original = "read1\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t60";
        let record = PafRecord::from_line(original).unwrap();
        let output = record.to_line();

        let record2 = PafRecord::from_line(&output).unwrap();
        assert_eq!(record, record2);
    }

    #[test]
    fn test_paf_invalid_strand() {
        let line = "read1\t10000\t100\t9900\tX\tchr1\t50000000\t1000\t10900\t9500\t9800\t60";
        let result = PafRecord::from_line(line);
        assert!(result.is_err());
    }

    #[test]
    fn test_paf_parser_iteration() {
        let paf_data = "read1\t10000\t100\t9900\t+\tchr1\t50000000\t1000\t10900\t9500\t9800\t60
read2\t8000\t50\t7950\t-\tchr2\t40000000\t2000\t9900\t7700\t7900\t55
read3\t12000\t200\t11800\t+\tchr1\t50000000\t5000\t16600\t11300\t11600\t58
";

        let cursor = std::io::Cursor::new(paf_data);
        let parser = PafParser::new(cursor);

        let records: Vec<_> = parser.collect::<Result<Vec<_>>>().unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].query_name, "read1");
        assert_eq!(records[1].query_name, "read2");
        assert_eq!(records[2].query_name, "read3");
    }

    #[test]
    fn test_paf_zero_length_alignment() {
        let line = "read1\t10000\t100\t100\t+\tchr1\t50000000\t1000\t1000\t0\t0\t0";
        let record = PafRecord::from_line(line).unwrap();

        assert_eq!(record.query_aligned_length(), 0);
        assert_eq!(record.target_aligned_length(), 0);
        assert_eq!(record.identity(), 0.0);
    }
}
