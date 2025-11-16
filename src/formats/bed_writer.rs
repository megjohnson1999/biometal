//! BED format writer with support for BED3, BED6, BED12, and narrowPeak formats
//!
//! # Format Specification
//!
//! BED (Browser Extensible Data) is a tab-delimited format for genomic features:
//! - **BED3**: chrom, start, end (minimal)
//! - **BED6**: BED3 + name, score, strand
//! - **BED12**: BED6 + thick_start, thick_end, item_rgb, block_count, block_sizes, block_starts
//! - **narrowPeak**: BED6 + signal_value, p_value, q_value, peak
//!
//! # Architecture
//!
//! This writer follows biometal patterns:
//! - Automatic compression based on file extension (.gz, .bgz)
//! - cloudflare_zlib backend (1.67× faster decompression)
//! - Validation (coordinates, required fields)
//! - Streaming write (constant memory)
//!
//! # Example
//!
//! ```no_run
//! use biometal::formats::bed::{Bed3Record, BedWriter};
//! use biometal::types::GenomicInterval;
//!
//! # fn main() -> biometal::Result<()> {
//! let mut writer = BedWriter::create("output.bed.gz")?;
//!
//! let record = Bed3Record {
//!     interval: GenomicInterval::new("chr1".to_string(), 1000, 2000),
//! };
//!
//! writer.write_bed3(&record)?;
//! writer.finish()?;
//! # Ok(())
//! # }
//! ```

use crate::error::{BiometalError, Result};
use crate::formats::bed::{Bed12Record, Bed3Record, Bed6Record, NarrowPeakRecord};
use crate::formats::primitives::Strand;
use crate::io::compression::CompressedWriter;
use crate::io::sink::DataSink;
use std::io::Write;
use std::path::Path;

/// BED format writer with compression support
///
/// # Features
///
/// - Automatic compression (gzip, bgzip) based on file extension
/// - Supports all BED variants (BED3, BED6, BED12, narrowPeak)
/// - Validation (coordinates, required fields)
/// - Streaming write (constant memory)
///
/// # Example
///
/// ```no_run
/// use biometal::formats::bed::{Bed6Record, Bed3Record, BedWriter};
/// use biometal::types::{GenomicInterval, Strand};
///
/// # fn main() -> biometal::Result<()> {
/// let mut writer = BedWriter::create("output.bed.gz")?;
///
/// let record = Bed6Record {
///     bed3: Bed3Record {
///         interval: GenomicInterval::new("chr1".to_string(), 1000, 2000),
///     },
///     name: Some("feature1".to_string()),
///     score: Some(500),
///     strand: Some(Strand::Forward),
/// };
///
/// writer.write_bed6(&record)?;
/// writer.finish()?;
/// # Ok(())
/// # }
/// ```
pub struct BedWriter {
    writer: CompressedWriter,
    records_written: usize,
}

impl BedWriter {
    /// Create a new BED writer from a data sink
    ///
    /// Automatically detects compression from file extension:
    /// - `.gz` → gzip compression
    /// - `.bgz` → bgzip compression
    /// - other → uncompressed
    pub fn new(sink: DataSink) -> Result<Self> {
        let writer = CompressedWriter::new(sink)
            .map_err(|e| BiometalError::Io(e))?;
        Ok(Self {
            writer,
            records_written: 0,
        })
    }

    /// Create a BED writer from a file path
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::bed::BedWriter;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut writer = BedWriter::create("output.bed.gz")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::new(DataSink::from_path(path))
    }

    /// Create a BED writer to stdout
    ///
    /// Useful for streaming pipelines:
    /// ```bash
    /// biometal filter input.bed.gz | bedtools intersect -a stdin -b regions.bed
    /// ```
    pub fn stdout() -> Result<Self> {
        Self::new(DataSink::stdout())
    }

    /// Write a BED3 record (chrom, start, end)
    ///
    /// # Arguments
    ///
    /// * `record` - BED3 record to write
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Chromosome name is empty
    /// - Start >= end
    /// - An I/O error occurs
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::bed::{Bed3Record, BedWriter};
    /// use biometal::types::GenomicInterval;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut writer = BedWriter::create("output.bed")?;
    ///
    /// let record = Bed3Record {
    ///     interval: GenomicInterval::new("chr1".to_string(), 1000, 2000),
    /// };
    ///
    /// writer.write_bed3(&record)?;
    /// writer.finish()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_bed3(&mut self, record: &Bed3Record) -> Result<()> {
        // Validate record
        if record.interval.chrom.is_empty() {
            return Err(BiometalError::InvalidInput {
                msg: "BED: Chromosome name cannot be empty".to_string(),
            });
        }

        if record.interval.start >= record.interval.end {
            return Err(BiometalError::InvalidInput {
                msg: format!(
                    "BED: Invalid interval: start ({}) >= end ({})",
                    record.interval.start, record.interval.end
                ),
            });
        }

        // Write BED3 line: chrom\tstart\tend
        writeln!(
            self.writer,
            "{}\t{}\t{}",
            record.interval.chrom,
            record.interval.start,
            record.interval.end
        )
        .map_err(|e| BiometalError::Io(e))?;

        self.records_written += 1;
        Ok(())
    }

    /// Write a BED6 record (BED3 + name, score, strand)
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::bed::{Bed3Record, Bed6Record, BedWriter};
    /// use biometal::types::{GenomicInterval, Strand};
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut writer = BedWriter::create("output.bed")?;
    ///
    /// let record = Bed6Record {
    ///     bed3: Bed3Record {
    ///         interval: GenomicInterval::new("chr1".to_string(), 1000, 2000),
    ///     },
    ///     name: Some("feature1".to_string()),
    ///     score: Some(500),
    ///     strand: Some(Strand::Forward),
    /// };
    ///
    /// writer.write_bed6(&record)?;
    /// writer.finish()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_bed6(&mut self, record: &Bed6Record) -> Result<()> {
        // Validate via BED3
        if record.bed3.interval.chrom.is_empty() {
            return Err(BiometalError::InvalidInput {
                msg: "Chromosome name cannot be empty".to_string(),
            });
        }

        if record.bed3.interval.start >= record.bed3.interval.end {
            return Err(BiometalError::InvalidInput {
                msg: format!(
                    "Invalid interval: start ({}) >= end ({})",
                    record.bed3.interval.start, record.bed3.interval.end
                ),
            });
        }

        // Format optional fields
        let name = record.name.as_ref().map_or(".", |s| s.as_str());
        let score = record.score.map_or(".".to_string(), |s| s.to_string());
        let strand = match record.strand {
            Some(Strand::Forward) => "+",
            Some(Strand::Reverse) => "-",
            Some(Strand::Unknown) | None => ".",
        };

        // Write BED6 line
        writeln!(
            self.writer,
            "{}\t{}\t{}\t{}\t{}\t{}",
            record.bed3.interval.chrom,
            record.bed3.interval.start,
            record.bed3.interval.end,
            name,
            score,
            strand
        )
        .map_err(|e| BiometalError::Io(e))?;

        self.records_written += 1;
        Ok(())
    }

    /// Write a BED12 record (full gene model with exons)
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::bed::{Bed3Record, Bed6Record, Bed12Record, BedWriter};
    /// use biometal::types::{GenomicInterval, Strand};
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut writer = BedWriter::create("output.bed")?;
    ///
    /// let record = Bed12Record {
    ///     bed6: Bed6Record {
    ///         bed3: Bed3Record {
    ///             interval: GenomicInterval::new("chr1".to_string(), 1000, 2000),
    ///         },
    ///         name: Some("transcript1".to_string()),
    ///         score: Some(1000),
    ///         strand: Some(Strand::Forward),
    ///     },
    ///     thick_start: Some(1100),
    ///     thick_end: Some(1900),
    ///     item_rgb: Some("255,0,0".to_string()),
    ///     block_count: Some(2),
    ///     block_sizes: Some(vec![100, 100]),
    ///     block_starts: Some(vec![0, 900]),
    /// };
    ///
    /// writer.write_bed12(&record)?;
    /// writer.finish()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_bed12(&mut self, record: &Bed12Record) -> Result<()> {
        // Validate via BED6
        if record.bed6.bed3.interval.chrom.is_empty() {
            return Err(BiometalError::InvalidInput {
                msg: "Chromosome name cannot be empty".to_string(),
            });
        }

        if record.bed6.bed3.interval.start >= record.bed6.bed3.interval.end {
            return Err(BiometalError::InvalidInput {
                msg: format!(
                    "Invalid interval: start ({}) >= end ({})",
                    record.bed6.bed3.interval.start, record.bed6.bed3.interval.end
                ),
            });
        }

        // Validate block consistency
        if let Some(count) = record.block_count {
            if let Some(ref sizes) = record.block_sizes {
                if sizes.len() != count as usize {
                    return Err(BiometalError::InvalidInput {
                        msg: format!(
                            "BED: block_count ({}) doesn't match block_sizes length ({})",
                            count,
                            sizes.len()
                        ),
                    });
                }
            }
            if let Some(ref starts) = record.block_starts {
                if starts.len() != count as usize {
                    return Err(BiometalError::InvalidInput {
                        msg: format!(
                            "BED: block_count ({}) doesn't match block_starts length ({})",
                            count,
                            starts.len()
                        ),
                    });
                }
            }
        }

        // Format BED6 fields
        let name = record.bed6.name.as_ref().map_or(".", |s| s.as_str());
        let score = record.bed6.score.map_or(".".to_string(), |s| s.to_string());
        let strand = match record.bed6.strand {
            Some(Strand::Forward) => "+",
            Some(Strand::Reverse) => "-",
            Some(Strand::Unknown) | None => ".",
        };

        // Format BED12 fields
        let thick_start = record.thick_start.map_or(".".to_string(), |v| v.to_string());
        let thick_end = record.thick_end.map_or(".".to_string(), |v| v.to_string());
        let item_rgb = record.item_rgb.as_ref().map_or(".", |s| s.as_str());
        let block_count = record.block_count.map_or(".".to_string(), |v| v.to_string());

        let block_sizes = record.block_sizes.as_ref().map_or(".".to_string(), |sizes| {
            sizes.iter().map(|s| s.to_string()).collect::<Vec<_>>().join(",")
        });

        let block_starts = record.block_starts.as_ref().map_or(".".to_string(), |starts| {
            starts.iter().map(|s| s.to_string()).collect::<Vec<_>>().join(",")
        });

        // Write BED12 line
        writeln!(
            self.writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.bed6.bed3.interval.chrom,
            record.bed6.bed3.interval.start,
            record.bed6.bed3.interval.end,
            name,
            score,
            strand,
            thick_start,
            thick_end,
            item_rgb,
            block_count,
            block_sizes,
            block_starts
        )
        .map_err(|e| BiometalError::Io(e))?;

        self.records_written += 1;
        Ok(())
    }

    /// Write a narrowPeak record (ENCODE ChIP-seq peak format)
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::bed::{Bed3Record, Bed6Record, NarrowPeakRecord, BedWriter};
    /// use biometal::types::{GenomicInterval, Strand};
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut writer = BedWriter::create("peaks.narrowPeak.gz")?;
    ///
    /// let record = NarrowPeakRecord {
    ///     bed6: Bed6Record {
    ///         bed3: Bed3Record {
    ///             interval: GenomicInterval::new("chr1".to_string(), 1000, 2000),
    ///         },
    ///         name: Some("peak1".to_string()),
    ///         score: Some(800),
    ///         strand: Some(Strand::Forward),
    ///     },
    ///     signal_value: 15.3,
    ///     p_value: 10.5,
    ///     q_value: 8.2,
    ///     peak: Some(500),
    /// };
    ///
    /// writer.write_narrowpeak(&record)?;
    /// writer.finish()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_narrowpeak(&mut self, record: &NarrowPeakRecord) -> Result<()> {
        // Validate via BED6
        if record.bed6.bed3.interval.chrom.is_empty() {
            return Err(BiometalError::InvalidInput {
                msg: "Chromosome name cannot be empty".to_string(),
            });
        }

        if record.bed6.bed3.interval.start >= record.bed6.bed3.interval.end {
            return Err(BiometalError::InvalidInput {
                msg: format!(
                    "Invalid interval: start ({}) >= end ({})",
                    record.bed6.bed3.interval.start, record.bed6.bed3.interval.end
                ),
            });
        }

        // Format BED6 fields
        let name = record.bed6.name.as_ref().map_or(".", |s| s.as_str());
        let score = record.bed6.score.map_or(".".to_string(), |s| s.to_string());
        let strand = match record.bed6.strand {
            Some(Strand::Forward) => "+",
            Some(Strand::Reverse) => "-",
            Some(Strand::Unknown) | None => ".",
        };

        // Format narrowPeak fields
        let signal = if record.signal_value < 0.0 {
            "-1".to_string()
        } else {
            format!("{:.2}", record.signal_value)
        };

        let p_val = if record.p_value < 0.0 {
            "-1".to_string()
        } else {
            format!("{:.2}", record.p_value)
        };

        let q_val = if record.q_value < 0.0 {
            "-1".to_string()
        } else {
            format!("{:.2}", record.q_value)
        };

        let peak = record.peak.map_or("-1".to_string(), |p| p.to_string());

        // Write narrowPeak line
        writeln!(
            self.writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.bed6.bed3.interval.chrom,
            record.bed6.bed3.interval.start,
            record.bed6.bed3.interval.end,
            name,
            score,
            strand,
            signal,
            p_val,
            q_val,
            peak
        )
        .map_err(|e| BiometalError::Io(e))?;

        self.records_written += 1;
        Ok(())
    }

    /// Get the number of records written so far
    pub fn records_written(&self) -> usize {
        self.records_written
    }

    /// Flush buffered data to disk
    ///
    /// You typically don't need to call this explicitly as `finish()`
    /// will flush automatically. However, it can be useful for
    /// long-running processes to ensure data is persisted.
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()
            .map_err(|e| BiometalError::Io(e))
    }

    /// Finish writing and flush all data
    ///
    /// This method MUST be called to ensure all data is written to disk.
    /// It flushes the internal buffers and closes the compression stream.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::bed::{Bed3Record, BedWriter};
    /// use biometal::types::GenomicInterval;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut writer = BedWriter::create("output.bed.gz")?;
    ///
    /// let record = Bed3Record {
    ///     interval: GenomicInterval::new("chr1".to_string(), 1000, 2000),
    /// };
    ///
    /// writer.write_bed3(&record)?;
    /// writer.finish()?;  // IMPORTANT: Flush and close
    /// # Ok(())
    /// # }
    /// ```
    pub fn finish(mut self) -> Result<()> {
        self.writer.flush()
            .map_err(|e| BiometalError::Io(e))?;
        Ok(())
    }
}
