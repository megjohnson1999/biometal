//! GFF3 format writer with compression support
//!
//! # Format Specification
//!
//! GFF3 (General Feature Format version 3) is a 9-column tab-delimited format:
//! 1. seqid - Chromosome/contig name
//! 2. source - Annotation source (e.g., NCBI, Ensembl)
//! 3. type - Feature type (gene, mRNA, exon, CDS, etc.)
//! 4. start - Start position (1-based, inclusive)
//! 5. end - End position (1-based, inclusive)
//! 6. score - Confidence score or "." if missing
//! 7. strand - +, -, or .
//! 8. phase - CDS phase (0, 1, 2) or "." if missing
//! 9. attributes - Semicolon-separated key=value pairs (ID=gene1;Name=ABC)
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
//! use biometal::formats::gff::Gff3Record;
//! use biometal::formats::gff_writer::Gff3Writer;
//! use biometal::formats::primitives::Strand;
//! use std::collections::HashMap;
//!
//! # fn main() -> biometal::Result<()> {
//! let mut writer = Gff3Writer::create("output.gff3.gz")?;
//!
//! let mut attributes = HashMap::new();
//! attributes.insert("ID".to_string(), "gene1".to_string());
//! attributes.insert("Name".to_string(), "ABC".to_string());
//!
//! let record = Gff3Record {
//!     seqid: "chr1".to_string(),
//!     source: "Ensembl".to_string(),
//!     feature_type: "gene".to_string(),
//!     start: 1000,
//!     end: 2000,
//!     score: Some(100.0),
//!     strand: Strand::Forward,
//!     phase: None,
//!     attributes,
//! };
//!
//! writer.write_record(&record)?;
//! writer.finish()?;
//! # Ok(())
//! # }
//! ```

use crate::error::{BiometalError, Result};
use crate::formats::gff::Gff3Record;
use crate::formats::primitives::Strand;
use crate::io::compression::CompressedWriter;
use crate::io::sink::DataSink;
use std::io::Write;
use std::path::Path;

/// GFF3 format writer with compression support
///
/// # Features
///
/// - Automatic compression (gzip, bgzip) based on file extension
/// - Validation (coordinates, required fields)
/// - Streaming write (constant memory)
/// - Proper attribute formatting (key=value;key2=value2)
///
/// # Example
///
/// ```no_run
/// use biometal::formats::gff::Gff3Record;
/// use biometal::formats::gff_writer::Gff3Writer;
/// use biometal::formats::primitives::Strand;
/// use std::collections::HashMap;
///
/// # fn main() -> biometal::Result<()> {
/// let mut writer = Gff3Writer::create("annotations.gff3.gz")?;
///
/// let mut attributes = HashMap::new();
/// attributes.insert("ID".to_string(), "gene1".to_string());
///
/// let record = Gff3Record {
///     seqid: "chr1".to_string(),
///     source: "NCBI".to_string(),
///     feature_type: "gene".to_string(),
///     start: 1000,
///     end: 2000,
///     score: None,
///     strand: Strand::Forward,
///     phase: None,
///     attributes,
/// };
///
/// writer.write_record(&record)?;
/// writer.finish()?;
/// # Ok(())
/// # }
/// ```
pub struct Gff3Writer {
    writer: CompressedWriter,
    records_written: usize,
    header_written: bool,
}

impl Gff3Writer {
    /// Create a new GFF3 writer from a data sink
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
            header_written: false,
        })
    }

    /// Create a GFF3 writer from a file path
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::gff_writer::Gff3Writer;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut writer = Gff3Writer::create("output.gff3.gz")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::new(DataSink::from_path(path))
    }

    /// Create a GFF3 writer to stdout
    ///
    /// Useful for streaming pipelines:
    /// ```bash
    /// biometal filter input.gff3.gz | bedtools intersect -a stdin -b regions.bed
    /// ```
    pub fn stdout() -> Result<Self> {
        Self::new(DataSink::stdout())
    }

    /// Write the GFF3 header directive
    ///
    /// This writes "##gff-version 3" to the file. Called automatically
    /// on first record write if not already written.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::gff_writer::Gff3Writer;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut writer = Gff3Writer::create("output.gff3")?;
    /// writer.write_header()?;  // Explicit header
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_header(&mut self) -> Result<()> {
        if !self.header_written {
            writeln!(self.writer, "##gff-version 3")
                .map_err(|e| BiometalError::Io(e))?;
            self.header_written = true;
        }
        Ok(())
    }

    /// Write a single GFF3 record
    ///
    /// # Arguments
    ///
    /// * `record` - GFF3 record to write
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - seqid is empty
    /// - source is empty
    /// - feature_type is empty
    /// - start >= end
    /// - phase is invalid (not 0, 1, or 2)
    /// - An I/O error occurs
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::gff::Gff3Record;
    /// use biometal::formats::gff_writer::Gff3Writer;
    /// use biometal::formats::primitives::Strand;
    /// use std::collections::HashMap;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut writer = Gff3Writer::create("output.gff3")?;
    ///
    /// let mut attributes = HashMap::new();
    /// attributes.insert("ID".to_string(), "gene1".to_string());
    ///
    /// let record = Gff3Record {
    ///     seqid: "chr1".to_string(),
    ///     source: "Ensembl".to_string(),
    ///     feature_type: "gene".to_string(),
    ///     start: 1000,
    ///     end: 2000,
    ///     score: Some(50.0),
    ///     strand: Strand::Forward,
    ///     phase: None,
    ///     attributes,
    /// };
    ///
    /// writer.write_record(&record)?;
    /// writer.finish()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_record(&mut self, record: &Gff3Record) -> Result<()> {
        // Write header if not already written
        if !self.header_written {
            self.write_header()?;
        }

        // Validate record
        if record.seqid.is_empty() {
            return Err(BiometalError::InvalidInput {
                msg: "GFF3: seqid cannot be empty".to_string(),
            });
        }

        if record.source.is_empty() {
            return Err(BiometalError::InvalidInput {
                msg: "GFF3: source cannot be empty".to_string(),
            });
        }

        if record.feature_type.is_empty() {
            return Err(BiometalError::InvalidInput {
                msg: "GFF3: feature_type cannot be empty".to_string(),
            });
        }

        if record.start >= record.end {
            return Err(BiometalError::InvalidInput {
                msg: format!(
                    "GFF3: Invalid interval: start ({}) >= end ({})",
                    record.start, record.end
                ),
            });
        }

        // Validate phase if present
        if let Some(phase) = record.phase {
            if phase > 2 {
                return Err(BiometalError::InvalidInput {
                    msg: format!("GFF3: Invalid phase: {} (must be 0, 1, or 2)", phase),
                });
            }
        }

        // Format fields
        let score_str = record.score.map_or(".".to_string(), |s| format!("{:.2}", s));

        let strand_str = match record.strand {
            Strand::Forward => "+",
            Strand::Reverse => "-",
            Strand::Unknown => ".",
        };

        let phase_str = record.phase.map_or(".".to_string(), |p| p.to_string());

        // Format attributes as semicolon-separated key=value pairs
        let attributes_str = if record.attributes.is_empty() {
            ".".to_string()
        } else {
            let mut attrs: Vec<String> = record.attributes
                .iter()
                .map(|(k, v)| format!("{}={}", k, v))
                .collect();
            // Sort for consistent output
            attrs.sort();
            attrs.join(";")
        };

        // Write GFF3 line
        writeln!(
            self.writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.seqid,
            record.source,
            record.feature_type,
            record.start,
            record.end,
            score_str,
            strand_str,
            phase_str,
            attributes_str
        )
        .map_err(|e| BiometalError::Io(e))?;

        self.records_written += 1;
        Ok(())
    }

    /// Write multiple GFF3 records from an iterator
    ///
    /// Convenience method for writing many records. The iterator can be
    /// any type that yields `Result<Gff3Record>`.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::gff_writer::Gff3Writer;
    /// use biometal::formats::primitives::TabDelimitedParser;
    /// use biometal::formats::gff::Gff3Record;
    /// use std::fs::File;
    /// use std::io::BufReader;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let file = File::open("input.gff3")?;
    /// let reader = BufReader::new(file);
    /// let parser = TabDelimitedParser::<_, Gff3Record>::new(reader);
    ///
    /// let mut writer = Gff3Writer::create("output.gff3.gz")?;
    /// writer.write_all(parser)?;
    /// writer.finish()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_all<I>(&mut self, records: I) -> Result<()>
    where
        I: IntoIterator<Item = Result<Gff3Record>>,
    {
        for record in records {
            self.write_record(&record?)?;
        }
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
    /// use biometal::formats::gff::Gff3Record;
    /// use biometal::formats::gff_writer::Gff3Writer;
    /// use biometal::formats::primitives::Strand;
    /// use std::collections::HashMap;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut writer = Gff3Writer::create("output.gff3.gz")?;
    ///
    /// let mut attributes = HashMap::new();
    /// attributes.insert("ID".to_string(), "gene1".to_string());
    ///
    /// let record = Gff3Record {
    ///     seqid: "chr1".to_string(),
    ///     source: "Ensembl".to_string(),
    ///     feature_type: "gene".to_string(),
    ///     start: 1000,
    ///     end: 2000,
    ///     score: None,
    ///     strand: Strand::Forward,
    ///     phase: None,
    ///     attributes,
    /// };
    ///
    /// writer.write_record(&record)?;
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
