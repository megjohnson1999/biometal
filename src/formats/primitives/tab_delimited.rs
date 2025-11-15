//! Generic tab-delimited file parsing.
//!
//! This module provides generic infrastructure for parsing tab-delimited formats
//! like BED, GFF, VCF, and GFA. All these formats share common patterns:
//! - Tab-delimited fields
//! - Comment lines (starting with `#`)
//! - Line-based records
//!
//! # Design
//!
//! The [`TabDelimitedRecord`] trait defines the interface for parsing records.
//! The [`TabDelimitedParser`] provides a generic streaming parser that works
//! with any type implementing this trait.
//!
//! # Examples
//!
//! ```
//! use biometal::formats::primitives::{TabDelimitedRecord, TabDelimitedParser, Result};
//!
//! // Define a custom record type
//! #[derive(Debug, PartialEq)]
//! struct SimpleRecord {
//!     chrom: String,
//!     start: u64,
//!     end: u64,
//! }
//!
//! impl TabDelimitedRecord for SimpleRecord {
//!     fn from_line(line: &str) -> Result<Self> {
//!         let fields: Vec<_> = line.split('\t').collect();
//!         if fields.len() < 3 {
//!             return Err(biometal::formats::primitives::FormatError::FieldCount {
//!                 expected: 3,
//!                 actual: fields.len(),
//!                 line: 0,
//!             });
//!         }
//!
//!         Ok(SimpleRecord {
//!             chrom: fields[0].to_string(),
//!             start: fields[1].parse().unwrap(),
//!             end: fields[2].parse().unwrap(),
//!         })
//!     }
//!
//!     fn to_line(&self) -> String {
//!         format!("{}\t{}\t{}", self.chrom, self.start, self.end)
//!     }
//! }
//!
//! // Parse from string
//! let data = "chr1\t100\t200\nchr2\t300\t400\n";
//! let parser = TabDelimitedParser::<_, SimpleRecord>::new(data.as_bytes());
//!
//! let records: Vec<_> = parser.collect::<Result<_>>().unwrap();
//! assert_eq!(records.len(), 2);
//! assert_eq!(records[0].chrom, "chr1");
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::formats::primitives::Result;
#[cfg(test)]
use crate::formats::primitives::FormatError;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::marker::PhantomData;
use std::path::Path;

/// Trait for types that can be parsed from tab-delimited lines.
///
/// Implement this trait to create custom parsers for tab-delimited formats.
///
/// # Examples
///
/// ```
/// use biometal::formats::primitives::{TabDelimitedRecord, Result};
///
/// #[derive(Debug, PartialEq)]
/// struct BedRecord {
///     chrom: String,
///     start: u64,
///     end: u64,
/// }
///
/// impl TabDelimitedRecord for BedRecord {
///     fn from_line(line: &str) -> Result<Self> {
///         let fields: Vec<_> = line.split('\t').collect();
///         Ok(BedRecord {
///             chrom: fields[0].to_string(),
///             start: fields[1].parse().unwrap(),
///             end: fields[2].parse().unwrap(),
///         })
///     }
///
///     fn to_line(&self) -> String {
///         format!("{}\t{}\t{}", self.chrom, self.start, self.end)
///     }
/// }
/// ```
pub trait TabDelimitedRecord: Sized {
    /// Parse a record from a tab-delimited line.
    ///
    /// The line should not include the trailing newline.
    ///
    /// # Errors
    ///
    /// Returns an error if the line is malformed or contains invalid data.
    fn from_line(line: &str) -> Result<Self>;

    /// Serialize this record to a tab-delimited line.
    ///
    /// The returned string should not include a trailing newline.
    fn to_line(&self) -> String;

    /// Expected number of tab-delimited fields.
    ///
    /// Returns `None` if the number of fields is variable.
    /// Default implementation returns `None`.
    fn expected_fields() -> Option<usize> {
        None
    }
}

/// Generic streaming parser for tab-delimited formats.
///
/// Parses records one at a time with constant memory usage.
/// Automatically skips:
/// - Empty lines
/// - Comment lines (starting with `#`)
///
/// # Type Parameters
///
/// - `R`: The underlying reader (anything implementing `Read`)
/// - `T`: The record type (must implement `TabDelimitedRecord`)
///
/// # Examples
///
/// ## Parse from file
///
/// ```no_run
/// use biometal::formats::primitives::{TabDelimitedParser, TabDelimitedRecord, Result};
///
/// # #[derive(Debug)]
/// # struct BedRecord { chrom: String, start: u64, end: u64 }
/// # impl TabDelimitedRecord for BedRecord {
/// #     fn from_line(line: &str) -> Result<Self> { todo!() }
/// #     fn to_line(&self) -> String { todo!() }
/// # }
/// # fn main() -> Result<()> {
/// let parser = TabDelimitedParser::<_, BedRecord>::from_path("data.bed")?;
///
/// for record in parser {
///     let record = record?;
///     // Process record
/// }
/// # Ok(())
/// # }
/// ```
///
/// ## Parse from compressed file
///
/// ```no_run
/// # use biometal::formats::primitives::{TabDelimitedParser, TabDelimitedRecord, Result};
/// # #[derive(Debug)]
/// # struct BedRecord { chrom: String, start: u64, end: u64 }
/// # impl TabDelimitedRecord for BedRecord {
/// #     fn from_line(line: &str) -> Result<Self> { todo!() }
/// #     fn to_line(&self) -> String { todo!() }
/// # }
/// # fn main() -> Result<()> {
/// // Automatically decompresses with cloudflare_zlib
/// let parser = TabDelimitedParser::<_, BedRecord>::from_bgzip_path("data.bed.gz")?;
///
/// for record in parser {
///     let record = record?;
///     // Process record
/// }
/// # Ok(())
/// # }
/// ```
///
/// ## Parse from HTTP
///
/// ```no_run
/// # use biometal::formats::primitives::{TabDelimitedParser, TabDelimitedRecord, Result};
/// # #[derive(Debug)]
/// # struct BedRecord { chrom: String, start: u64, end: u64 }
/// # impl TabDelimitedRecord for BedRecord {
/// #     fn from_line(line: &str) -> Result<Self> { todo!() }
/// #     fn to_line(&self) -> String { todo!() }
/// # }
/// # fn main() -> Result<()> {
/// let parser = TabDelimitedParser::<_, BedRecord>::from_url(
///     "https://example.com/data.bed"
/// )?;
///
/// for record in parser {
///     let record = record?;
///     // Process record
/// }
/// # Ok(())
/// # }
/// ```
pub struct TabDelimitedParser<R: Read, T: TabDelimitedRecord> {
    reader: BufReader<R>,
    line_buf: String,
    line_number: usize,
    _phantom: PhantomData<T>,
}

impl<R: Read, T: TabDelimitedRecord> TabDelimitedParser<R, T> {
    /// Creates a new parser from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::formats::primitives::{TabDelimitedParser, TabDelimitedRecord, Result};
    ///
    /// # #[derive(Debug)]
    /// # struct SimpleRecord { data: String }
    /// # impl TabDelimitedRecord for SimpleRecord {
    /// #     fn from_line(line: &str) -> Result<Self> { Ok(SimpleRecord { data: line.to_string() }) }
    /// #     fn to_line(&self) -> String { self.data.clone() }
    /// # }
    /// let data = "field1\tfield2\tfield3\n";
    /// let parser = TabDelimitedParser::<_, SimpleRecord>::new(data.as_bytes());
    /// ```
    pub fn new(reader: R) -> Self {
        TabDelimitedParser {
            reader: BufReader::new(reader),
            line_buf: String::with_capacity(1024),
            line_number: 0,
            _phantom: PhantomData,
        }
    }

    /// Returns the current line number (1-based).
    ///
    /// Useful for error reporting.
    pub fn line_number(&self) -> usize {
        self.line_number
    }
}

impl<T: TabDelimitedRecord> TabDelimitedParser<File, T> {
    /// Creates a parser from a file path.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use biometal::formats::primitives::{TabDelimitedParser, TabDelimitedRecord, Result};
    ///
    /// # #[derive(Debug)]
    /// # struct BedRecord { chrom: String, start: u64, end: u64 }
    /// # impl TabDelimitedRecord for BedRecord {
    /// #     fn from_line(line: &str) -> Result<Self> { todo!() }
    /// #     fn to_line(&self) -> String { todo!() }
    /// # }
    /// # fn main() -> Result<()> {
    /// let parser = TabDelimitedParser::<_, BedRecord>::from_path("data.bed")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self> {
        let file = File::open(path)?;
        Ok(Self::new(file))
    }
}

impl<T: TabDelimitedRecord> TabDelimitedParser<MultiGzDecoder<File>, T> {
    /// Creates a parser from a gzip/bgzip-compressed file.
    ///
    /// Uses cloudflare_zlib backend for optimal decompression performance.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened or is not valid gzip.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use biometal::formats::primitives::{TabDelimitedParser, TabDelimitedRecord, Result};
    ///
    /// # #[derive(Debug)]
    /// # struct BedRecord { chrom: String, start: u64, end: u64 }
    /// # impl TabDelimitedRecord for BedRecord {
    /// #     fn from_line(line: &str) -> Result<Self> { todo!() }
    /// #     fn to_line(&self) -> String { todo!() }
    /// # }
    /// # fn main() -> Result<()> {
    /// let parser = TabDelimitedParser::<_, BedRecord>::from_bgzip_path("data.bed.gz")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_bgzip_path(path: impl AsRef<Path>) -> Result<Self> {
        let file = File::open(path)?;
        let decoder = MultiGzDecoder::new(file);
        Ok(Self::new(decoder))
    }
}

#[cfg(feature = "network")]
impl<T: TabDelimitedRecord> TabDelimitedParser<Box<dyn Read>, T> {
    /// Creates a parser from an HTTP/HTTPS URL.
    ///
    /// Supports streaming from remote files without downloading.
    ///
    /// # Errors
    ///
    /// Returns an error if the URL is invalid or cannot be fetched.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use biometal::formats::primitives::{TabDelimitedParser, TabDelimitedRecord, Result};
    ///
    /// # #[derive(Debug)]
    /// # struct BedRecord { chrom: String, start: u64, end: u64 }
    /// # impl TabDelimitedRecord for BedRecord {
    /// #     fn from_line(line: &str) -> Result<Self> { todo!() }
    /// #     fn to_line(&self) -> String { todo!() }
    /// # }
    /// # fn main() -> Result<()> {
    /// let parser = TabDelimitedParser::<_, BedRecord>::from_url(
    ///     "https://example.com/data.bed"
    /// )?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_url(url: &str) -> Result<Self> {
        let response = reqwest::blocking::get(url)?;
        let reader: Box<dyn Read> = Box::new(response);
        Ok(Self::new(reader))
    }
}

impl<R: Read, T: TabDelimitedRecord> Iterator for TabDelimitedParser<R, T> {
    type Item = Result<T>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.line_buf.clear();

            match self.reader.read_line(&mut self.line_buf) {
                Ok(0) => return None, // EOF
                Ok(_) => {
                    self.line_number += 1;

                    // Trim trailing newline
                    let line = self.line_buf.trim_end();

                    // Skip empty lines
                    if line.is_empty() {
                        continue;
                    }

                    // Skip comments (lines starting with #)
                    if line.starts_with('#') {
                        continue;
                    }

                    // Parse record
                    return Some(T::from_line(line));
                }
                Err(e) => return Some(Err(e.into())),
            }
        }
    }
}

/// Generic writer for tab-delimited formats.
///
/// Writes records one at a time with automatic compression support.
/// Works with any type implementing [`TabDelimitedRecord`].
///
/// # Type Parameters
///
/// - `T`: The record type (must implement `TabDelimitedRecord`)
///
/// # Features
///
/// - Automatic compression (gzip, bgzip) based on file extension
/// - Streaming write (constant memory)
/// - Validation via record's `to_line()` method
///
/// # Examples
///
/// ## Write BED records
///
/// ```no_run
/// use biometal::formats::bed::Bed3Record;
/// use biometal::formats::primitives::{TabDelimitedWriter, TabDelimitedRecord};
/// use biometal::formats::primitives::GenomicInterval;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let mut writer = TabDelimitedWriter::create("output.bed.gz")?;
///
/// let interval = GenomicInterval::new("chr1".to_string(), 1000, 2000)?;
/// let record = Bed3Record { interval };
///
/// writer.write_record(&record)?;
/// writer.finish()?;
/// # Ok(())
/// # }
/// ```
///
/// ## Write from stream
///
/// ```no_run
/// use biometal::formats::bed::Bed6Record;
/// use biometal::formats::primitives::{TabDelimitedParser, TabDelimitedWriter};
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let parser = TabDelimitedParser::<_, Bed6Record>::from_path("input.bed")?;
/// let mut writer = TabDelimitedWriter::create("output.bed.gz")?;
///
/// for record in parser {
///     writer.write_record(&record?)?;
/// }
///
/// writer.finish()?;
/// # Ok(())
/// # }
/// ```
pub struct TabDelimitedWriter<T: TabDelimitedRecord> {
    writer: crate::io::compression::CompressedWriter,
    records_written: usize,
    _phantom: PhantomData<T>,
}

impl<T: TabDelimitedRecord> TabDelimitedWriter<T> {
    /// Create a new writer from a data sink
    ///
    /// Automatically detects compression from file extension:
    /// - `.gz` → gzip compression
    /// - `.bgz` → bgzip compression
    /// - other → uncompressed
    pub fn new(sink: crate::io::sink::DataSink) -> Result<Self> {
        let writer = crate::io::compression::CompressedWriter::new(sink)
            .map_err(|e| crate::formats::primitives::FormatError::Io(e))?;
        Ok(Self {
            writer,
            records_written: 0,
            _phantom: PhantomData,
        })
    }

    /// Create a writer from a file path
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use biometal::formats::bed::Bed3Record;
    /// use biometal::formats::primitives::TabDelimitedWriter;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let writer = TabDelimitedWriter::<Bed3Record>::create("output.bed.gz")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::new(crate::io::sink::DataSink::from_path(path))
    }

    /// Create a writer to stdout
    ///
    /// Useful for streaming pipelines:
    /// ```bash
    /// biometal filter input.bed | biometal stats
    /// ```
    pub fn stdout() -> Result<Self> {
        Self::new(crate::io::sink::DataSink::stdout())
    }

    /// Write a single record
    ///
    /// # Arguments
    ///
    /// * `record` - Record to write
    ///
    /// # Errors
    ///
    /// Returns an error if an I/O error occurs.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use biometal::formats::bed::Bed3Record;
    /// use biometal::formats::primitives::{TabDelimitedWriter, GenomicInterval};
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let mut writer = TabDelimitedWriter::create("output.bed")?;
    ///
    /// let interval = GenomicInterval::new("chr1".to_string(), 1000, 2000)?;
    /// let record = Bed3Record { interval };
    ///
    /// writer.write_record(&record)?;
    /// writer.finish()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_record(&mut self, record: &T) -> Result<()> {
        use std::io::Write;

        let line = record.to_line();
        writeln!(self.writer, "{}", line)
            .map_err(|e| crate::formats::primitives::FormatError::Io(e))?;

        self.records_written += 1;
        Ok(())
    }

    /// Write multiple records from an iterator
    ///
    /// Convenience method for writing many records. The iterator can be
    /// any type that yields `Result<T>`, such as `TabDelimitedParser<_, T>`.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use biometal::formats::bed::Bed6Record;
    /// use biometal::formats::primitives::{TabDelimitedParser, TabDelimitedWriter};
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let input = TabDelimitedParser::<_, Bed6Record>::from_path("input.bed")?;
    /// let mut writer = TabDelimitedWriter::create("output.bed.gz")?;
    ///
    /// writer.write_all(input)?;
    /// writer.finish()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_all<I>(&mut self, records: I) -> Result<()>
    where
        I: IntoIterator<Item = Result<T>>,
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
        use std::io::Write;
        self.writer.flush()
            .map_err(|e| crate::formats::primitives::FormatError::Io(e))
    }

    /// Finish writing and flush all data
    ///
    /// This method MUST be called to ensure all data is written to disk.
    /// It flushes the internal buffers and closes the compression stream.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use biometal::formats::bed::Bed3Record;
    /// use biometal::formats::primitives::{TabDelimitedWriter, GenomicInterval};
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let mut writer = TabDelimitedWriter::create("output.bed.gz")?;
    ///
    /// let interval = GenomicInterval::new("chr1".to_string(), 1000, 2000)?;
    /// let record = Bed3Record { interval };
    ///
    /// writer.write_record(&record)?;
    /// writer.finish()?;  // IMPORTANT: Flush and close
    /// # Ok(())
    /// # }
    /// ```
    pub fn finish(mut self) -> Result<()> {
        use std::io::Write;
        self.writer.flush()
            .map_err(|e| crate::formats::primitives::FormatError::Io(e))?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Test record type
    #[derive(Debug, PartialEq)]
    struct TestRecord {
        chrom: String,
        start: u64,
        end: u64,
    }

    impl TabDelimitedRecord for TestRecord {
        fn from_line(line: &str) -> Result<Self> {
            let fields: Vec<_> = line.split('\t').collect();
            if fields.len() < 3 {
                return Err(FormatError::FieldCount {
                    expected: 3,
                    actual: fields.len(),
                    line: 0,
                });
            }

            Ok(TestRecord {
                chrom: fields[0].to_string(),
                start: fields[1].parse().map_err(|e| FormatError::InvalidField {
                    field: "start".to_string(),
                    line: 0,
                    reason: format!("{}", e),
                })?,
                end: fields[2].parse().map_err(|e| FormatError::InvalidField {
                    field: "end".to_string(),
                    line: 0,
                    reason: format!("{}", e),
                })?,
            })
        }

        fn to_line(&self) -> String {
            format!("{}\t{}\t{}", self.chrom, self.start, self.end)
        }

        fn expected_fields() -> Option<usize> {
            Some(3)
        }
    }

    #[test]
    fn test_parse_basic() {
        let data = "chr1\t100\t200\nchr2\t300\t400\n";
        let parser = TabDelimitedParser::<_, TestRecord>::new(data.as_bytes());

        let records: Vec<_> = parser.collect::<Result<_>>().unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[0].start, 100);
        assert_eq!(records[0].end, 200);
        assert_eq!(records[1].chrom, "chr2");
        assert_eq!(records[1].start, 300);
        assert_eq!(records[1].end, 400);
    }

    #[test]
    fn test_parse_skip_comments() {
        let data = "# This is a comment\nchr1\t100\t200\n# Another comment\nchr2\t300\t400\n";
        let parser = TabDelimitedParser::<_, TestRecord>::new(data.as_bytes());

        let records: Vec<_> = parser.collect::<Result<_>>().unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[1].chrom, "chr2");
    }

    #[test]
    fn test_parse_skip_empty_lines() {
        let data = "chr1\t100\t200\n\n\nchr2\t300\t400\n";
        let parser = TabDelimitedParser::<_, TestRecord>::new(data.as_bytes());

        let records: Vec<_> = parser.collect::<Result<_>>().unwrap();

        assert_eq!(records.len(), 2);
    }

    #[test]
    fn test_parse_mixed() {
        let data = "# Header\n\nchr1\t100\t200\n# Comment\n\nchr2\t300\t400\n\n# Trailer\n";
        let parser = TabDelimitedParser::<_, TestRecord>::new(data.as_bytes());

        let records: Vec<_> = parser.collect::<Result<_>>().unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[1].chrom, "chr2");
    }

    #[test]
    fn test_to_line_round_trip() {
        let original = TestRecord {
            chrom: "chr1".to_string(),
            start: 100,
            end: 200,
        };

        let line = original.to_line();
        let parsed = TestRecord::from_line(&line).unwrap();

        assert_eq!(parsed, original);
    }

    #[test]
    fn test_line_number_tracking() {
        let data = "# Comment\nchr1\t100\t200\nchr2\t300\t400\n";
        let mut parser = TabDelimitedParser::<_, TestRecord>::new(data.as_bytes());

        // Before first record
        assert_eq!(parser.line_number(), 0);

        // After first record (line 2, after comment)
        let _ = parser.next();
        assert_eq!(parser.line_number(), 2);

        // After second record
        let _ = parser.next();
        assert_eq!(parser.line_number(), 3);
    }

    // TabDelimitedWriter tests
    #[test]
    fn test_writer_basic() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        {
            let mut writer = TabDelimitedWriter::<TestRecord>::create(path).unwrap();

            let record = TestRecord {
                chrom: "chr1".to_string(),
                start: 1000,
                end: 2000,
            };

            writer.write_record(&record).unwrap();
            writer.finish().unwrap();
        }

        // Read back and verify
        let parser = TabDelimitedParser::<_, TestRecord>::from_path(path).unwrap();
        let records: Vec<_> = parser.collect::<Result<_>>().unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[0].start, 1000);
        assert_eq!(records[0].end, 2000);
    }

    #[test]
    fn test_writer_multiple_records() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        {
            let mut writer = TabDelimitedWriter::<TestRecord>::create(path).unwrap();

            writer.write_record(&TestRecord {
                chrom: "chr1".to_string(),
                start: 1000,
                end: 2000,
            }).unwrap();

            writer.write_record(&TestRecord {
                chrom: "chr2".to_string(),
                start: 3000,
                end: 4000,
            }).unwrap();

            assert_eq!(writer.records_written(), 2);
            writer.finish().unwrap();
        }

        // Read back and verify
        let parser = TabDelimitedParser::<_, TestRecord>::from_path(path).unwrap();
        let records: Vec<_> = parser.collect::<Result<_>>().unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[1].chrom, "chr2");
    }

    #[test]
    fn test_writer_write_all() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        {
            let mut writer = TabDelimitedWriter::<TestRecord>::create(path).unwrap();

            let records = vec![
                Ok(TestRecord { chrom: "chr1".to_string(), start: 100, end: 200 }),
                Ok(TestRecord { chrom: "chr2".to_string(), start: 300, end: 400 }),
                Ok(TestRecord { chrom: "chr3".to_string(), start: 500, end: 600 }),
            ];

            writer.write_all(records).unwrap();
            assert_eq!(writer.records_written(), 3);
            writer.finish().unwrap();
        }

        // Read back and verify
        let parser = TabDelimitedParser::<_, TestRecord>::from_path(path).unwrap();
        let records: Vec<_> = parser.collect::<Result<_>>().unwrap();
        assert_eq!(records.len(), 3);
    }

    #[test]
    fn test_writer_round_trip() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        // Original records
        let original_records = vec![
            TestRecord { chrom: "chr1".to_string(), start: 1000, end: 2000 },
            TestRecord { chrom: "chr2".to_string(), start: 3000, end: 4000 },
            TestRecord { chrom: "chr3".to_string(), start: 5000, end: 6000 },
        ];

        // Write records
        {
            let mut writer = TabDelimitedWriter::<TestRecord>::create(path).unwrap();

            for record in &original_records {
                writer.write_record(record).unwrap();
            }

            writer.finish().unwrap();
        }

        // Read back
        let parser = TabDelimitedParser::<_, TestRecord>::from_path(path).unwrap();
        let parsed_records: Vec<_> = parser.collect::<Result<_>>().unwrap();

        // Verify
        assert_eq!(parsed_records.len(), original_records.len());
        for (parsed, original) in parsed_records.iter().zip(original_records.iter()) {
            assert_eq!(parsed, original);
        }
    }

    #[test]
    fn test_writer_with_bed3() {
        use tempfile::NamedTempFile;
        use crate::formats::bed::Bed3Record;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        {
            let mut writer = TabDelimitedWriter::<Bed3Record>::create(path).unwrap();

            let interval = crate::formats::primitives::GenomicInterval::new(
                "chr1".to_string(), 1000, 2000
            ).unwrap();
            let record = Bed3Record { interval };

            writer.write_record(&record).unwrap();
            writer.finish().unwrap();
        }

        // Read back and verify
        let parser = TabDelimitedParser::<_, Bed3Record>::from_path(path).unwrap();
        let records: Vec<_> = parser.collect::<Result<_>>().unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].interval.chrom, "chr1");
        assert_eq!(records[0].interval.start, 1000);
        assert_eq!(records[0].interval.end, 2000);
    }

    #[test]
    fn test_writer_with_bed6() {
        use tempfile::NamedTempFile;
        use crate::formats::bed::Bed6Record;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        {
            let mut writer = TabDelimitedWriter::<Bed6Record>::create(path).unwrap();

            let interval = crate::formats::primitives::GenomicInterval::new(
                "chr1".to_string(), 1000, 2000
            ).unwrap();

            let bed3 = crate::formats::bed::Bed3Record { interval };
            let record = Bed6Record {
                bed3,
                name: Some("gene1".to_string()),
                score: Some(100),
                strand: Some(crate::formats::primitives::Strand::Forward),
            };

            writer.write_record(&record).unwrap();
            writer.finish().unwrap();
        }

        // Read back and verify
        let parser = TabDelimitedParser::<_, Bed6Record>::from_path(path).unwrap();
        let records: Vec<_> = parser.collect::<Result<_>>().unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].bed3.interval.chrom, "chr1");
        assert_eq!(records[0].name, Some("gene1".to_string()));
        assert_eq!(records[0].score, Some(100));
    }

    #[test]
    fn test_writer_bed3_round_trip() {
        use tempfile::NamedTempFile;
        use crate::formats::bed::Bed3Record;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let original = Bed3Record {
            interval: crate::formats::primitives::GenomicInterval::new(
                "chr1".to_string(), 1000, 2000
            ).unwrap(),
        };

        // Write
        {
            let mut writer = TabDelimitedWriter::<Bed3Record>::create(path).unwrap();
            writer.write_record(&original).unwrap();
            writer.finish().unwrap();
        }

        // Read back
        let parser = TabDelimitedParser::<_, Bed3Record>::from_path(path).unwrap();
        let parsed: Vec<_> = parser.collect::<Result<_>>().unwrap();

        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].interval.chrom, original.interval.chrom);
        assert_eq!(parsed[0].interval.start, original.interval.start);
        assert_eq!(parsed[0].interval.end, original.interval.end);
    }
}
