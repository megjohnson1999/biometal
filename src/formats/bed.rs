//! BED (Browser Extensible Data) format parser.
//!
//! BED is a tab-delimited format for representing genomic features:
//! - **BED3**: Minimal format (chrom, start, end)
//! - **BED6**: Standard format (+ name, score, strand)
//! - **BED12**: Full format (+ exon/block information)
//! - **narrowPeak**: ENCODE ChIP-seq peaks (BED6 + signal, pValue, qValue, peak offset)
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

/// narrowPeak record: ENCODE ChIP-seq peaks (BED6 + signalValue, pValue, qValue, peak).
///
/// This format is used by ENCODE for called peaks of signal enrichment based on
/// pooled, normalized ChIP-seq data. It extends BED6 with peak-calling statistics.
///
/// # Format
///
/// ```text
/// chrom  start  end    name   score  strand  signalValue  pValue  qValue  peak
/// chr1   1000   2000   peak1  100    .       12.5         8.3     5.2     450
/// ```
///
/// # Fields (BED6+4)
///
/// **BED6 fields:**
/// - **chrom**: Chromosome name
/// - **start**: 0-based start position
/// - **end**: Exclusive end position
/// - **name**: Peak identifier (use `.` if unassigned)
/// - **score**: Display darkness 0-1000 (optional, `.` for missing)
/// - **strand**: `+`, `-`, or `.` (not applicable for most ChIP-seq)
///
/// **narrowPeak extensions:**
/// - **signalValue**: Overall enrichment measurement (averaged)
/// - **pValue**: Statistical significance as -log10 (use -1 if unavailable)
/// - **qValue**: FDR significance as -log10 (use -1 if unavailable)
/// - **peak**: 0-based offset from start to peak summit (use -1 if not called)
///
/// # Examples
///
/// ```
/// use biometal::formats::bed::NarrowPeakRecord;
/// use biometal::formats::TabDelimitedRecord;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// // Full narrowPeak with all fields
/// let line = "chr1\t1000\t2000\tpeak1\t100\t.\t12.5\t8.3\t5.2\t450";
/// let record = NarrowPeakRecord::from_line(line)?;
///
/// assert_eq!(record.bed6.bed3.interval.chrom, "chr1");
/// assert_eq!(record.signal_value, 12.5);
/// assert_eq!(record.p_value, 8.3);
/// assert_eq!(record.q_value, 5.2);
/// assert_eq!(record.peak, Some(450));
///
/// // Peak summit absolute position
/// let summit = record.peak_position();
/// assert_eq!(summit, Some(1450)); // start(1000) + peak(450)
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct NarrowPeakRecord {
    /// BED6 fields (chrom, start, end, name, score, strand)
    pub bed6: Bed6Record,
    /// Overall enrichment measurement (averaged)
    pub signal_value: f64,
    /// Statistical significance as -log10 (-1 if unavailable)
    pub p_value: f64,
    /// FDR significance as -log10 (-1 if unavailable)
    pub q_value: f64,
    /// 0-based offset from start to peak summit (None if -1)
    pub peak: Option<i32>,
}

impl NarrowPeakRecord {
    /// Get the absolute genomic position of the peak summit.
    ///
    /// Returns `None` if peak offset is not available (-1 in file).
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::formats::bed::NarrowPeakRecord;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let record = NarrowPeakRecord::from_line("chr1\t1000\t2000\tpeak1\t100\t.\t12.5\t8.3\t5.2\t450")?;
    /// assert_eq!(record.peak_position(), Some(1450)); // 1000 + 450
    /// # Ok(())
    /// # }
    /// ```
    pub fn peak_position(&self) -> Option<u64> {
        self.peak.map(|offset| {
            self.bed6.bed3.interval.start + (offset as u64)
        })
    }

    /// Check if this peak passes significance thresholds.
    ///
    /// # Arguments
    ///
    /// * `p_threshold` - Minimum -log10(p-value) (e.g., 2.0 for p < 0.01)
    /// * `q_threshold` - Minimum -log10(q-value) (e.g., 1.3 for q < 0.05)
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::formats::bed::NarrowPeakRecord;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let record = NarrowPeakRecord::from_line("chr1\t1000\t2000\tpeak1\t100\t.\t12.5\t8.3\t5.2\t450")?;
    ///
    /// // Significant at p < 0.01 (10^-8.3) and q < 0.05 (10^-5.2)
    /// assert!(record.is_significant(2.0, 1.3));
    /// # Ok(())
    /// # }
    /// ```
    pub fn is_significant(&self, p_threshold: f64, q_threshold: f64) -> bool {
        self.p_value >= p_threshold && self.q_value >= q_threshold
    }
}

impl TabDelimitedRecord for NarrowPeakRecord {
    fn from_line(line: &str) -> Result<Self> {
        let fields = split_fields(line, Some(10), 0)?;

        // Parse BED6 fields
        let chrom = fields[0].to_string();
        let start: u64 = parse_required(fields[1], "start", 0)?;
        let end: u64 = parse_required(fields[2], "end", 0)?;
        let interval = GenomicInterval::new(chrom, start, end)?;
        let bed3 = Bed3Record { interval };

        let name: Option<String> = parse_optional(fields[3], "name", 0)?;
        let score: Option<u32> = parse_optional(fields[4], "score", 0)?;
        let strand = Some(Strand::from_str(fields[5])?);
        let bed6 = Bed6Record {
            bed3,
            name,
            score,
            strand,
        };

        // Parse narrowPeak-specific fields
        let signal_value: f64 = parse_required(fields[6], "signalValue", 0)?;
        let p_value: f64 = parse_required(fields[7], "pValue", 0)?;
        let q_value: f64 = parse_required(fields[8], "qValue", 0)?;
        let peak_offset: i32 = parse_required(fields[9], "peak", 0)?;
        let peak = if peak_offset == -1 {
            None
        } else {
            Some(peak_offset)
        };

        Ok(NarrowPeakRecord {
            bed6,
            signal_value,
            p_value,
            q_value,
            peak,
        })
    }

    fn to_line(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.bed6.bed3.interval.chrom,
            self.bed6.bed3.interval.start,
            self.bed6.bed3.interval.end,
            self.bed6.name.as_ref().map(|s| s.as_str()).unwrap_or("."),
            self.bed6.score.map(|s| s.to_string()).unwrap_or_else(|| ".".to_string()),
            self.bed6.strand.map(|s| s.to_string()).unwrap_or_else(|| ".".to_string()),
            self.signal_value,
            self.p_value,
            self.q_value,
            self.peak.unwrap_or(-1)
        )
    }

    fn expected_fields() -> Option<usize> {
        Some(10)
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

    // narrowPeak tests
    #[test]
    fn test_narrowpeak_all_fields() {
        let line = "chr1\t1000\t2000\tpeak1\t100\t.\t12.5\t8.3\t5.2\t450";
        let record = NarrowPeakRecord::from_line(line).unwrap();

        assert_eq!(record.bed6.bed3.interval.chrom, "chr1");
        assert_eq!(record.bed6.bed3.interval.start, 1000);
        assert_eq!(record.bed6.bed3.interval.end, 2000);
        assert_eq!(record.bed6.name, Some("peak1".to_string()));
        assert_eq!(record.bed6.score, Some(100));
        assert_eq!(record.signal_value, 12.5);
        assert_eq!(record.p_value, 8.3);
        assert_eq!(record.q_value, 5.2);
        assert_eq!(record.peak, Some(450));
    }

    #[test]
    fn test_narrowpeak_peak_position() {
        let line = "chr1\t1000\t2000\tpeak1\t100\t.\t12.5\t8.3\t5.2\t450";
        let record = NarrowPeakRecord::from_line(line).unwrap();

        // Peak summit = start + offset = 1000 + 450 = 1450
        assert_eq!(record.peak_position(), Some(1450));
    }

    #[test]
    fn test_narrowpeak_no_peak() {
        let line = "chr1\t1000\t2000\tpeak1\t100\t.\t12.5\t8.3\t5.2\t-1";
        let record = NarrowPeakRecord::from_line(line).unwrap();

        assert_eq!(record.peak, None);
        assert_eq!(record.peak_position(), None);
    }

    #[test]
    fn test_narrowpeak_significance() {
        let line = "chr1\t1000\t2000\tpeak1\t100\t.\t12.5\t8.3\t5.2\t450";
        let record = NarrowPeakRecord::from_line(line).unwrap();

        // p=8.3 means 10^-8.3, q=5.2 means 10^-5.2
        assert!(record.is_significant(2.0, 1.3)); // p<0.01, q<0.05
        assert!(record.is_significant(8.0, 5.0)); // Strict thresholds
        assert!(!record.is_significant(9.0, 5.0)); // p threshold too high
        assert!(!record.is_significant(8.0, 6.0)); // q threshold too high
    }

    #[test]
    fn test_narrowpeak_unavailable_values() {
        let line = "chr1\t1000\t2000\tpeak1\t100\t.\t12.5\t-1\t-1\t-1";
        let record = NarrowPeakRecord::from_line(line).unwrap();

        assert_eq!(record.p_value, -1.0);
        assert_eq!(record.q_value, -1.0);
        assert_eq!(record.peak, None);
    }

    #[test]
    fn test_narrowpeak_round_trip() {
        let original = "chr1\t1000\t2000\tpeak1\t100\t.\t12.5\t8.3\t5.2\t450";
        let record = NarrowPeakRecord::from_line(original).unwrap();
        let output = record.to_line();

        let record2 = NarrowPeakRecord::from_line(&output).unwrap();
        assert_eq!(record, record2);
    }

    #[test]
    fn test_narrowpeak_missing_bed6_values() {
        let line = "chr1\t1000\t2000\t.\t.\t.\t12.5\t8.3\t5.2\t450";
        let record = NarrowPeakRecord::from_line(line).unwrap();

        assert_eq!(record.bed6.name, None);
        assert_eq!(record.bed6.score, None);
        assert_eq!(record.bed6.strand, Some(Strand::Unknown));
    }

    #[test]
    fn test_narrowpeak_different_strands() {
        let forward = "chr1\t1000\t2000\tpeak1\t100\t+\t12.5\t8.3\t5.2\t450";
        let reverse = "chr1\t1000\t2000\tpeak2\t100\t-\t12.5\t8.3\t5.2\t450";
        let unknown = "chr1\t1000\t2000\tpeak3\t100\t.\t12.5\t8.3\t5.2\t450";

        let record1 = NarrowPeakRecord::from_line(forward).unwrap();
        let record2 = NarrowPeakRecord::from_line(reverse).unwrap();
        let record3 = NarrowPeakRecord::from_line(unknown).unwrap();

        assert_eq!(record1.bed6.strand, Some(Strand::Forward));
        assert_eq!(record2.bed6.strand, Some(Strand::Reverse));
        assert_eq!(record3.bed6.strand, Some(Strand::Unknown));
    }

    #[test]
    fn test_narrowpeak_score_ranges() {
        let low = "chr1\t1000\t2000\tpeak1\t100\t.\t1.0\t2.0\t1.5\t450";
        let high = "chr1\t1000\t2000\tpeak2\t900\t.\t100.0\t50.0\t40.0\t450";

        let record1 = NarrowPeakRecord::from_line(low).unwrap();
        let record2 = NarrowPeakRecord::from_line(high).unwrap();

        assert_eq!(record1.bed6.score, Some(100));
        assert_eq!(record2.bed6.score, Some(900));
        assert_eq!(record1.signal_value, 1.0);
        assert_eq!(record2.signal_value, 100.0);
    }

    #[test]
    fn test_narrowpeak_zero_offset_peak() {
        let line = "chr1\t1000\t2000\tpeak1\t100\t.\t12.5\t8.3\t5.2\t0";
        let record = NarrowPeakRecord::from_line(line).unwrap();

        // Peak at start of region
        assert_eq!(record.peak, Some(0));
        assert_eq!(record.peak_position(), Some(1000));
    }
}
