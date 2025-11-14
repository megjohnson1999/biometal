//! BED (Browser Extensible Data) format parser.
//!
//! BED is a tab-delimited format for representing genomic features:
//! - **BED3**: Minimal format (chrom, start, end)
//! - **BED6**: Standard format (+ name, score, strand)
//! - **BED12**: Full format (+ exon/block information)
//!
//! # Format Specification
//!
//! - **Coordinates**: 0-based, half-open `[start, end)`
//! - **Columns**: Tab-delimited (3, 6, or 12 columns)
//! - **Comments**: Lines starting with `#` or `track`/`browser` directives
//!
//! # Examples
//!
//! ## BED3 - Minimal genomic intervals
//!
//! ```
//! use biometal::formats::bed::Bed3Record;
//! use biometal::formats::TabDelimitedRecord;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let line = "chr1\t1000\t2000";
//! let record = Bed3Record::from_line(line)?;
//!
//! assert_eq!(record.interval.chrom, "chr1");
//! assert_eq!(record.interval.start, 1000);
//! assert_eq!(record.interval.end, 2000);
//! assert_eq!(record.interval.length(), 1000);
//! # Ok(())
//! # }
//! ```
//!
//! ## BED6 - Standard annotations
//!
//! ```
//! use biometal::formats::bed::Bed6Record;
//! use biometal::formats::{TabDelimitedRecord, Strand};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let line = "chr1\t1000\t2000\tgene1\t100\t+";
//! let record = Bed6Record::from_line(line)?;
//!
//! assert_eq!(record.bed3.interval.chrom, "chr1");
//! assert_eq!(record.name, Some("gene1".to_string()));
//! assert_eq!(record.score, Some(100));
//! assert_eq!(record.strand, Some(Strand::Forward));
//! # Ok(())
//! # }
//! ```
//!
//! ## BED12 - Full gene models
//!
//! ```
//! use biometal::formats::bed::Bed12Record;
//! use biometal::formats::TabDelimitedRecord;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let line = "chr1\t1000\t2000\tgene1\t100\t+\t1200\t1800\t0,0,255\t2\t400,400\t0,600";
//! let record = Bed12Record::from_line(line)?;
//!
//! assert_eq!(record.bed6.bed3.interval.chrom, "chr1");
//! assert_eq!(record.thick_start, Some(1200));
//! assert_eq!(record.thick_end, Some(1800));
//! assert_eq!(record.block_count, Some(2));
//! # Ok(())
//! # }
//! ```
//!
//! ## Streaming parser
//!
//! ```no_run
//! use biometal::formats::bed::Bed6Record;
//! use biometal::formats::TabDelimitedParser;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let parser = TabDelimitedParser::<_, Bed6Record>::from_path("data.bed")?;
//!
//! for result in parser {
//!     let record = result?;
//!     println!("{}: {}-{}", record.bed3.interval.chrom,
//!              record.bed3.interval.start, record.bed3.interval.end);
//! }
//! # Ok(())
//! # }
//! ```

use crate::formats::primitives::{
    fields::{parse_comma_list, parse_optional, parse_required, split_fields},
    GenomicInterval, Result, Strand, TabDelimitedRecord,
};
use std::str::FromStr;

/// BED3 record: minimal genomic interval (chrom, start, end).
///
/// This is the simplest BED format with only 3 required fields.
///
/// # Format
///
/// ```text
/// chrom  start  end
/// chr1   1000   2000
/// ```
///
/// # Coordinates
///
/// - **0-based, half-open**: `[start, end)`
/// - `start` is inclusive, `end` is exclusive
/// - Length = `end - start`
///
/// # Examples
///
/// ```
/// use biometal::formats::bed::Bed3Record;
/// use biometal::formats::TabDelimitedRecord;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let record = Bed3Record::from_line("chr1\t1000\t2000")?;
/// assert_eq!(record.interval.length(), 1000);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bed3Record {
    /// Genomic interval (chrom, start, end)
    pub interval: GenomicInterval,
}

impl TabDelimitedRecord for Bed3Record {
    fn from_line(line: &str) -> Result<Self> {
        let fields = split_fields(line, Some(3), 0)?;

        let chrom = fields[0].to_string();
        let start: u64 = parse_required(fields[1], "start", 0)?;
        let end: u64 = parse_required(fields[2], "end", 0)?;

        let interval = GenomicInterval::new(chrom, start, end)?;

        Ok(Bed3Record { interval })
    }

    fn to_line(&self) -> String {
        format!(
            "{}\t{}\t{}",
            self.interval.chrom, self.interval.start, self.interval.end
        )
    }

    fn expected_fields() -> Option<usize> {
        Some(3)
    }
}

/// BED6 record: standard genomic annotation (BED3 + name, score, strand).
///
/// This is the most commonly used BED format.
///
/// # Format
///
/// ```text
/// chrom  start  end    name    score  strand
/// chr1   1000   2000   gene1   100    +
/// ```
///
/// # Fields
///
/// - **name**: Feature name (optional, `.` for missing)
/// - **score**: 0-1000 quality score (optional, `.` for missing)
/// - **strand**: `+`, `-`, or `.` (unknown/not applicable)
///
/// # Examples
///
/// ```
/// use biometal::formats::bed::Bed6Record;
/// use biometal::formats::{TabDelimitedRecord, Strand};
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// // With all fields
/// let record = Bed6Record::from_line("chr1\t1000\t2000\tgene1\t500\t+")?;
/// assert_eq!(record.name, Some("gene1".to_string()));
/// assert_eq!(record.score, Some(500));
/// assert_eq!(record.strand, Some(Strand::Forward));
///
/// // With missing values
/// let record = Bed6Record::from_line("chr1\t1000\t2000\t.\t.\t.")?;
/// assert_eq!(record.name, None);
/// assert_eq!(record.score, None);
/// assert_eq!(record.strand, Some(Strand::Unknown));
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bed6Record {
    /// BED3 fields (chrom, start, end)
    pub bed3: Bed3Record,
    /// Feature name (optional)
    pub name: Option<String>,
    /// Score 0-1000 (optional)
    pub score: Option<u32>,
    /// Strand (optional, Unknown if `.`)
    pub strand: Option<Strand>,
}

impl TabDelimitedRecord for Bed6Record {
    fn from_line(line: &str) -> Result<Self> {
        let fields = split_fields(line, Some(6), 0)?;

        // Parse BED3 fields
        let chrom = fields[0].to_string();
        let start: u64 = parse_required(fields[1], "start", 0)?;
        let end: u64 = parse_required(fields[2], "end", 0)?;
        let interval = GenomicInterval::new(chrom, start, end)?;
        let bed3 = Bed3Record { interval };

        // Parse BED6 additional fields
        let name: Option<String> = parse_optional(fields[3], "name", 0)?;
        let score: Option<u32> = parse_optional(fields[4], "score", 0)?;
        let strand = Some(Strand::from_str(fields[5])?);

        Ok(Bed6Record {
            bed3,
            name,
            score,
            strand,
        })
    }

    fn to_line(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.bed3.interval.chrom,
            self.bed3.interval.start,
            self.bed3.interval.end,
            self.name
                .as_ref()
                .map(|s| s.as_str())
                .unwrap_or("."),
            self.score
                .map(|s| s.to_string())
                .unwrap_or_else(|| ".".to_string()),
            self.strand
                .map(|s| s.to_string())
                .unwrap_or_else(|| ".".to_string())
        )
    }

    fn expected_fields() -> Option<usize> {
        Some(6)
    }
}

/// BED12 record: full gene model with exon/block information.
///
/// This format is used for representing genes with multiple exons,
/// coding sequences, and visual display information.
///
/// # Format
///
/// ```text
/// chrom start end name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts
/// chr1  1000  5000 gene1 100  +      1200       4800     0,0,255 3          400,600,400 0,1000,4600
/// ```
///
/// # Fields
///
/// - **thickStart/thickEnd**: Coding sequence (CDS) region
/// - **itemRgb**: Display color (R,G,B)
/// - **blockCount**: Number of exons/blocks
/// - **blockSizes**: Comma-separated block lengths
/// - **blockStarts**: Comma-separated block start positions (relative to chrom start)
///
/// # Examples
///
/// ```
/// use biometal::formats::bed::Bed12Record;
/// use biometal::formats::TabDelimitedRecord;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let line = "chr1\t1000\t5000\tgene1\t100\t+\t1200\t4800\t0,0,255\t3\t400,600,400\t0,1000,4600";
/// let record = Bed12Record::from_line(line)?;
///
/// assert_eq!(record.block_count, Some(3));
/// assert_eq!(record.block_sizes, Some(vec![400, 600, 400]));
/// assert_eq!(record.block_starts, Some(vec![0, 1000, 4600]));
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bed12Record {
    /// BED6 fields
    pub bed6: Bed6Record,
    /// Start of coding sequence (optional)
    pub thick_start: Option<u64>,
    /// End of coding sequence (optional)
    pub thick_end: Option<u64>,
    /// RGB color for display (e.g., "255,0,0") (optional)
    pub item_rgb: Option<String>,
    /// Number of blocks/exons (optional)
    pub block_count: Option<u32>,
    /// Comma-separated block sizes (optional)
    pub block_sizes: Option<Vec<u32>>,
    /// Comma-separated block starts (relative to chromStart) (optional)
    pub block_starts: Option<Vec<u32>>,
}

impl TabDelimitedRecord for Bed12Record {
    fn from_line(line: &str) -> Result<Self> {
        let fields = split_fields(line, Some(12), 0)?;

        // Parse BED6 fields
        let bed6 = Bed6Record::from_line(&fields[..6].join("\t"))?;

        // Parse BED12 additional fields
        let thick_start: Option<u64> = parse_optional(fields[6], "thickStart", 0)?;
        let thick_end: Option<u64> = parse_optional(fields[7], "thickEnd", 0)?;
        let item_rgb: Option<String> = parse_optional(fields[8], "itemRgb", 0)?;
        let block_count: Option<u32> = parse_optional(fields[9], "blockCount", 0)?;

        // Parse comma-separated lists
        let block_sizes = if fields[10] == "." || fields[10].is_empty() {
            None
        } else {
            Some(
                parse_comma_list(fields[10])
                    .iter()
                    .map(|s| s.parse::<u32>())
                    .collect::<std::result::Result<Vec<_>, _>>()
                    .map_err(|e| crate::formats::primitives::FormatError::InvalidField {
                        field: "blockSizes".to_string(),
                        line: 0,
                        reason: e.to_string(),
                    })?,
            )
        };

        let block_starts = if fields[11] == "." || fields[11].is_empty() {
            None
        } else {
            Some(
                parse_comma_list(fields[11])
                    .iter()
                    .map(|s| s.parse::<u32>())
                    .collect::<std::result::Result<Vec<_>, _>>()
                    .map_err(|e| crate::formats::primitives::FormatError::InvalidField {
                        field: "blockStarts".to_string(),
                        line: 0,
                        reason: e.to_string(),
                    })?,
            )
        };

        Ok(Bed12Record {
            bed6,
            thick_start,
            thick_end,
            item_rgb,
            block_count,
            block_sizes,
            block_starts,
        })
    }

    fn to_line(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.bed6.to_line(),
            self.thick_start
                .map(|v| v.to_string())
                .unwrap_or_else(|| ".".to_string()),
            self.thick_end
                .map(|v| v.to_string())
                .unwrap_or_else(|| ".".to_string()),
            self.item_rgb.as_ref().map(|s| s.as_str()).unwrap_or("."),
            self.block_count
                .map(|v| v.to_string())
                .unwrap_or_else(|| ".".to_string()),
            format!(
                "{}\t{}",
                self.block_sizes
                    .as_ref()
                    .map(|v| v
                        .iter()
                        .map(|n| n.to_string())
                        .collect::<Vec<_>>()
                        .join(","))
                    .unwrap_or_else(|| ".".to_string()),
                self.block_starts
                    .as_ref()
                    .map(|v| v
                        .iter()
                        .map(|n| n.to_string())
                        .collect::<Vec<_>>()
                        .join(","))
                    .unwrap_or_else(|| ".".to_string())
            )
        )
    }

    fn expected_fields() -> Option<usize> {
        Some(12)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // BED3 tests
    #[test]
    fn test_bed3_basic() {
        let line = "chr1\t1000\t2000";
        let record = Bed3Record::from_line(line).unwrap();

        assert_eq!(record.interval.chrom, "chr1");
        assert_eq!(record.interval.start, 1000);
        assert_eq!(record.interval.end, 2000);
        assert_eq!(record.interval.length(), 1000);
    }

    #[test]
    fn test_bed3_round_trip() {
        let original = "chr1\t1000\t2000";
        let record = Bed3Record::from_line(original).unwrap();
        let output = record.to_line();
        assert_eq!(output, original);
    }

    #[test]
    fn test_bed3_invalid_interval() {
        let line = "chr1\t2000\t1000"; // start > end
        let result = Bed3Record::from_line(line);
        assert!(result.is_err());
    }

    // BED6 tests
    #[test]
    fn test_bed6_all_fields() {
        let line = "chr1\t1000\t2000\tgene1\t500\t+";
        let record = Bed6Record::from_line(line).unwrap();

        assert_eq!(record.bed3.interval.chrom, "chr1");
        assert_eq!(record.bed3.interval.start, 1000);
        assert_eq!(record.bed3.interval.end, 2000);
        assert_eq!(record.name, Some("gene1".to_string()));
        assert_eq!(record.score, Some(500));
        assert_eq!(record.strand, Some(Strand::Forward));
    }

    #[test]
    fn test_bed6_missing_values() {
        let line = "chr1\t1000\t2000\t.\t.\t.";
        let record = Bed6Record::from_line(line).unwrap();

        assert_eq!(record.name, None);
        assert_eq!(record.score, None);
        assert_eq!(record.strand, Some(Strand::Unknown));
    }

    #[test]
    fn test_bed6_reverse_strand() {
        let line = "chr2\t3000\t4000\tgene2\t800\t-";
        let record = Bed6Record::from_line(line).unwrap();

        assert_eq!(record.bed3.interval.chrom, "chr2");
        assert_eq!(record.strand, Some(Strand::Reverse));
    }

    #[test]
    fn test_bed6_round_trip() {
        let original = "chr1\t1000\t2000\tgene1\t500\t+";
        let record = Bed6Record::from_line(original).unwrap();
        let output = record.to_line();
        assert_eq!(output, original);
    }

    // BED12 tests
    #[test]
    fn test_bed12_full() {
        let line = "chr1\t1000\t5000\tgene1\t100\t+\t1200\t4800\t0,0,255\t3\t400,600,400\t0,1000,4600";
        let record = Bed12Record::from_line(line).unwrap();

        assert_eq!(record.bed6.bed3.interval.chrom, "chr1");
        assert_eq!(record.bed6.bed3.interval.start, 1000);
        assert_eq!(record.bed6.bed3.interval.end, 5000);
        assert_eq!(record.bed6.name, Some("gene1".to_string()));
        assert_eq!(record.thick_start, Some(1200));
        assert_eq!(record.thick_end, Some(4800));
        assert_eq!(record.item_rgb, Some("0,0,255".to_string()));
        assert_eq!(record.block_count, Some(3));
        assert_eq!(record.block_sizes, Some(vec![400, 600, 400]));
        assert_eq!(record.block_starts, Some(vec![0, 1000, 4600]));
    }

    #[test]
    fn test_bed12_minimal() {
        let line = "chr1\t1000\t2000\tgene1\t100\t+\t.\t.\t.\t.\t.\t.";
        let record = Bed12Record::from_line(line).unwrap();

        assert_eq!(record.thick_start, None);
        assert_eq!(record.thick_end, None);
        assert_eq!(record.item_rgb, None);
        assert_eq!(record.block_count, None);
        assert_eq!(record.block_sizes, None);
        assert_eq!(record.block_starts, None);
    }

    #[test]
    fn test_bed12_single_exon() {
        let line = "chr1\t1000\t2000\tgene1\t100\t+\t1100\t1900\t255,0,0\t1\t1000\t0";
        let record = Bed12Record::from_line(line).unwrap();

        assert_eq!(record.block_count, Some(1));
        assert_eq!(record.block_sizes, Some(vec![1000]));
        assert_eq!(record.block_starts, Some(vec![0]));
    }

    #[test]
    fn test_bed12_round_trip() {
        let original =
            "chr1\t1000\t5000\tgene1\t100\t+\t1200\t4800\t0,0,255\t3\t400,600,400\t0,1000,4600";
        let record = Bed12Record::from_line(original).unwrap();
        let output = record.to_line();
        assert_eq!(output, original);
    }
}
