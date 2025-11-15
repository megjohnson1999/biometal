//! FASTA format writer with compression support
//!
//! # Format
//!
//! FASTA format consists of:
//! - Header line starting with '>' followed by sequence identifier
//! - Sequence lines (wrapped at 80 characters by default)
//!
//! Example:
//! ```text
//! >sequence1 description
//! GATTACAGATTACATGCATGCAGATTACAGATTACATGCATGCAGATTACAGATTACATGCATGCA
//! GATTACAGATTACA
//! >sequence2
//! ACGTACGT
//! ```
//!
//! # Architecture
//!
//! This writer follows the same patterns as FastqWriter:
//! - Automatic compression based on file extension (.gz, .bgz)
//! - cloudflare_zlib backend (1.67× faster decompression)
//! - Configurable line width (default 80 characters)
//! - Validation (non-empty sequences, valid IDs)

use crate::error::{BiometalError, Result};
use crate::io::compression::CompressedWriter;
use crate::io::sink::DataSink;
use crate::types::FastaRecord;
use std::io::Write;
use std::path::Path;

/// Default line width for sequence wrapping (80 characters)
///
/// Standard FASTA convention wraps sequences at 80 characters.
/// Some tools (samtools, bwa) work better with wrapped sequences.
pub const DEFAULT_LINE_WIDTH: usize = 80;

/// FASTA format writer with compression support
///
/// # Features
///
/// - Automatic compression (gzip, bgzip) based on file extension
/// - Configurable sequence line width (default 80 characters)
/// - Validation (non-empty sequences, valid IDs)
/// - Streaming write (constant memory)
///
/// # Example
///
/// ```no_run
/// use biometal::io::fasta::FastaWriter;
/// use biometal::FastaRecord;
///
/// # fn main() -> biometal::Result<()> {
/// let mut writer = FastaWriter::create("output.fa.gz")?;
///
/// let record = FastaRecord::new(
///     "chr1".to_string(),
///     b"GATTACA".to_vec(),
/// );
///
/// writer.write_record(&record)?;
/// writer.finish()?;
/// # Ok(())
/// # }
/// ```
///
/// # Compression
///
/// Compression is automatically detected from file extension:
/// - `.gz` → gzip compression
/// - `.bgz` → bgzip compression (better for random access)
/// - other → uncompressed
///
/// # Sequence Wrapping
///
/// By default, sequences are wrapped at 80 characters. You can customize this:
///
/// ```no_run
/// use biometal::io::fasta::FastaWriter;
///
/// # fn main() -> biometal::Result<()> {
/// let mut writer = FastaWriter::create("output.fa")?
///     .with_line_width(60);
/// # Ok(())
/// # }
/// ```
///
/// Or disable wrapping entirely:
///
/// ```no_run
/// use biometal::io::fasta::FastaWriter;
///
/// # fn main() -> biometal::Result<()> {
/// let mut writer = FastaWriter::create("output.fa")?
///     .with_line_width(usize::MAX); // No wrapping
/// # Ok(())
/// # }
/// ```
pub struct FastaWriter {
    writer: CompressedWriter,
    line_width: usize,
    records_written: usize,
}

impl FastaWriter {
    /// Create a new FASTA writer from a data sink
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
            line_width: DEFAULT_LINE_WIDTH,
            records_written: 0,
        })
    }

    /// Create a FASTA writer from a file path
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::fasta::FastaWriter;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut writer = FastaWriter::create("output.fa.gz")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::new(DataSink::from_path(path))
    }

    /// Create a FASTA writer to stdout
    ///
    /// Useful for streaming pipelines:
    /// ```bash
    /// biometal filter input.fa.gz | biometal stats
    /// ```
    pub fn stdout() -> Result<Self> {
        Self::new(DataSink::stdout())
    }

    /// Set the line width for sequence wrapping
    ///
    /// # Arguments
    ///
    /// * `width` - Number of characters per line (use `usize::MAX` to disable wrapping)
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::fasta::FastaWriter;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut writer = FastaWriter::create("output.fa")?
    ///     .with_line_width(60);  // Wrap at 60 characters
    /// # Ok(())
    /// # }
    /// ```
    pub fn with_line_width(mut self, width: usize) -> Self {
        self.line_width = width;
        self
    }

    /// Write a single FASTA record
    ///
    /// # Arguments
    ///
    /// * `record` - FASTA record to write
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The sequence is empty
    /// - The ID is empty
    /// - An I/O error occurs
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::fasta::FastaWriter;
    /// use biometal::FastaRecord;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut writer = FastaWriter::create("output.fa")?;
    ///
    /// let record = FastaRecord::new(
    ///     "chr1".to_string(),
    ///     b"GATTACA".to_vec(),
    /// );
    ///
    /// writer.write_record(&record)?;
    /// writer.finish()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_record(&mut self, record: &FastaRecord) -> Result<()> {
        // Validate record
        if record.id.is_empty() {
            return Err(BiometalError::InvalidFastaFormat {
                line: self.records_written + 1,
                msg: "Sequence ID cannot be empty".to_string(),
            });
        }

        if record.sequence.is_empty() {
            return Err(BiometalError::InvalidFastaFormat {
                line: self.records_written + 1,
                msg: "Sequence cannot be empty".to_string(),
            });
        }

        // Write header line
        write!(self.writer, ">{}\n", record.id)
            .map_err(|e| BiometalError::Io(e))?;

        // Write sequence with wrapping
        if self.line_width == usize::MAX {
            // No wrapping - write as single line
            self.writer.write_all(&record.sequence)
                .map_err(|e| BiometalError::Io(e))?;
            self.writer.write_all(b"\n")
                .map_err(|e| BiometalError::Io(e))?;
        } else {
            // Write in chunks of line_width
            for chunk in record.sequence.chunks(self.line_width) {
                self.writer.write_all(chunk)
                    .map_err(|e| BiometalError::Io(e))?;
                self.writer.write_all(b"\n")
                    .map_err(|e| BiometalError::Io(e))?;
            }
        }

        self.records_written += 1;
        Ok(())
    }

    /// Write multiple FASTA records from an iterator
    ///
    /// Convenience method for writing many records. The iterator can be
    /// any type that yields `Result<FastaRecord>`, such as `FastaStream`.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::fasta::{FastaStream, FastaWriter};
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let input = FastaStream::from_path("input.fa.gz")?;
    /// let mut writer = FastaWriter::create("output.fa.gz")?;
    ///
    /// writer.write_all(input)?;
    /// writer.finish()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_all<I>(&mut self, records: I) -> Result<()>
    where
        I: IntoIterator<Item = Result<FastaRecord>>,
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
    /// use biometal::io::fasta::FastaWriter;
    /// use biometal::FastaRecord;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut writer = FastaWriter::create("output.fa.gz")?;
    ///
    /// let record = FastaRecord::new(
    ///     "chr1".to_string(),
    ///     b"GATTACA".to_vec(),
    /// );
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::fasta::FastaStream;
    use proptest::prelude::*;

    #[test]
    fn test_fasta_writer_basic() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        // Write a record
        {
            let mut writer = FastaWriter::create(path).unwrap();

            let record = FastaRecord::new(
                "chr1".to_string(),
                b"GATTACA".to_vec(),
            );

            writer.write_record(&record).unwrap();
            writer.finish().unwrap();
        }

        // Read it back
        let stream = FastaStream::from_path(path).unwrap();
        let records: Vec<_> = stream.collect::<Result<Vec<_>>>().unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id, "chr1");
        assert_eq!(records[0].sequence, b"GATTACA");
    }

    #[test]
    fn test_fasta_writer_wrapping() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        // Long sequence to test wrapping
        let long_seq = b"GATTACA".repeat(20); // 140 bases

        {
            let mut writer = FastaWriter::create(path).unwrap()
                .with_line_width(80);

            let record = FastaRecord::new(
                "chr1".to_string(),
                long_seq.clone(),
            );

            writer.write_record(&record).unwrap();
            writer.finish().unwrap();
        }

        // Read it back (should unwrap automatically)
        let stream = FastaStream::from_path(path).unwrap();
        let records: Vec<_> = stream.collect::<Result<Vec<_>>>().unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].sequence, long_seq);
    }

    #[test]
    fn test_fasta_writer_no_wrapping() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let long_seq = b"GATTACA".repeat(50); // 350 bases

        {
            let mut writer = FastaWriter::create(path).unwrap()
                .with_line_width(usize::MAX); // No wrapping

            let record = FastaRecord::new(
                "chr1".to_string(),
                long_seq.clone(),
            );

            writer.write_record(&record).unwrap();
            writer.finish().unwrap();
        }

        // Read it back
        let stream = FastaStream::from_path(path).unwrap();
        let records: Vec<_> = stream.collect::<Result<Vec<_>>>().unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].sequence, long_seq);
    }

    #[test]
    fn test_fasta_writer_write_all() {
        use tempfile::NamedTempFile;

        let temp_file_in = NamedTempFile::new().unwrap();
        let path_in = temp_file_in.path();
        let temp_file_out = NamedTempFile::new().unwrap();
        let path_out = temp_file_out.path();

        // Create test input
        {
            let mut writer = FastaWriter::create(path_in).unwrap();
            for i in 0..100 {
                let record = FastaRecord::new(
                    format!("chr{}", i),
                    b"GATTACA".to_vec(),
                );
                writer.write_record(&record).unwrap();
            }
            writer.finish().unwrap();
        }

        // Use write_all() to copy
        {
            let input = FastaStream::from_path(path_in).unwrap();
            let mut writer = FastaWriter::create(path_out).unwrap();

            writer.write_all(input).unwrap();

            assert_eq!(writer.records_written(), 100);
            writer.finish().unwrap();
        }

        // Verify output
        let output = FastaStream::from_path(path_out).unwrap();
        let count = output.count();
        assert_eq!(count, 100);
    }

    #[test]
    fn test_fasta_writer_empty_id() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let mut writer = FastaWriter::create(path).unwrap();

        let record = FastaRecord::new(
            "".to_string(), // Empty ID
            b"GATTACA".to_vec(),
        );

        let result = writer.write_record(&record);
        assert!(result.is_err());
    }

    #[test]
    fn test_fasta_writer_empty_sequence() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let mut writer = FastaWriter::create(path).unwrap();

        let record = FastaRecord::new(
            "chr1".to_string(),
            Vec::new(), // Empty sequence
        );

        let result = writer.write_record(&record);
        assert!(result.is_err());
    }

    // Property-based test for round-trip correctness
    proptest! {
        #[test]
        fn test_fasta_writer_roundtrip_property(
            id in "[a-zA-Z0-9_]{1,50}",
            seq in "[ACGT]{10,200}",
        ) {
            use tempfile::NamedTempFile;

            let temp_file = NamedTempFile::new().unwrap();
            let path = temp_file.path();

            let original = FastaRecord::new(
                id.clone(),
                seq.as_bytes().to_vec(),
            );

            // Write
            {
                let sink = DataSink::from_path(path);
                let mut writer = FastaWriter::new(sink).unwrap();
                writer.write_record(&original).unwrap();
                writer.finish().unwrap();
            }

            // Read back
            let stream = FastaStream::from_path(path).unwrap();
            let mut records = stream.collect::<Result<Vec<_>>>().unwrap();
            prop_assert_eq!(records.len(), 1);

            let read_back = records.remove(0);
            prop_assert_eq!(read_back.id, original.id);
            prop_assert_eq!(read_back.sequence, original.sequence);
        }
    }
}
