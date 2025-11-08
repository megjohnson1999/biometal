//! SAM format writer (BAM → SAM conversion).
//!
//! Provides primitives for writing BAM records to SAM format (text).
//! This is the inverse of BAM reading - a fundamental I/O primitive.
//!
//! # Design
//!
//! - Streaming writer (constant memory, Rule 5)
//! - Format conversion (BAM binary → SAM text)
//! - Header and record serialization
//! - No buffering (caller controls buffering)
//!
//! # Example
//!
//! ```no_run
//! use biometal::io::bam::{BamReader, SamWriter};
//! use std::fs::File;
//! use std::io::BufWriter;
//!
//! # fn main() -> biometal::Result<()> {
//! // Read BAM file
//! let mut bam = BamReader::from_path("input.bam")?;
//!
//! // Write to SAM file
//! let output = File::create("output.sam")?;
//! let mut writer = BufWriter::new(output);
//! let mut sam = SamWriter::new(&mut writer);
//!
//! // Write header
//! sam.write_header(bam.header())?;
//!
//! // Stream records (constant memory)
//! for record in bam.records() {
//!     let record = record?;
//!     sam.write_record(&record)?;
//! }
//! # Ok(())
//! # }
//! ```

use super::{Header, Record};
use std::io::{self, Write};

/// SAM format writer.
///
/// Writes BAM records to SAM text format. This is a streaming primitive
/// that maintains constant memory (Rule 5).
pub struct SamWriter<W: Write> {
    /// Underlying writer (typically BufWriter<File>)
    writer: W,
}

impl<W: Write> SamWriter<W> {
    /// Create a new SAM writer.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::io::bam::SamWriter;
    /// use std::io::BufWriter;
    ///
    /// let mut output = Vec::new();
    /// let writer = BufWriter::new(&mut output);
    /// let sam = SamWriter::new(writer);
    /// ```
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    /// Write SAM header.
    ///
    /// Writes the @HD, @SQ, and other header lines to the output.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use biometal::io::bam::{BamReader, SamWriter};
    /// # use std::io::BufWriter;
    /// # fn main() -> biometal::Result<()> {
    /// let mut bam = BamReader::from_path("input.bam")?;
    /// let mut output = Vec::new();
    /// let writer = BufWriter::new(&mut output);
    /// let mut sam = SamWriter::new(writer);
    ///
    /// sam.write_header(bam.header())?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_header(&mut self, header: &Header) -> io::Result<()> {
        // Write SAM header text (if present)
        if !header.text.is_empty() {
            // Header text already contains @HD, @SQ, etc. lines
            write!(self.writer, "{}", header.text)?;
            if !header.text.ends_with('\n') {
                writeln!(self.writer)?;
            }
        } else {
            // Generate minimal header
            writeln!(self.writer, "@HD\tVN:1.6\tSO:unknown")?;

            // Write reference sequences
            for reference in &header.references {
                writeln!(
                    self.writer,
                    "@SQ\tSN:{}\tLN:{}",
                    reference.name, reference.length
                )?;
            }
        }

        Ok(())
    }

    /// Write a single BAM record in SAM format.
    ///
    /// Converts the binary BAM record to SAM text format and writes it.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use biometal::io::bam::{BamReader, SamWriter};
    /// # use std::io::BufWriter;
    /// # fn main() -> biometal::Result<()> {
    /// let mut bam = BamReader::from_path("input.bam")?;
    /// let mut output = Vec::new();
    /// let writer = BufWriter::new(&mut output);
    /// let mut sam = SamWriter::new(writer);
    ///
    /// sam.write_header(bam.header())?;
    ///
    /// for record in bam.records() {
    ///     let record = record?;
    ///     sam.write_record(&record)?;
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        // SAM format: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAGS]

        // 1. QNAME (read name)
        write!(self.writer, "{}\t", record.name)?;

        // 2. FLAG
        write!(self.writer, "{}\t", record.flags)?;

        // 3. RNAME (reference name, "*" if unmapped)
        if let Some(ref_id) = record.reference_id {
            write!(self.writer, "{}\t", ref_id)?; // Note: Should lookup name from header
        } else {
            write!(self.writer, "*\t")?;
        }

        // 4. POS (1-based position, 0 if unmapped)
        if let Some(pos) = record.position {
            write!(self.writer, "{}\t", pos + 1)?; // BAM is 0-based, SAM is 1-based
        } else {
            write!(self.writer, "0\t")?;
        }

        // 5. MAPQ (255 if unavailable)
        write!(self.writer, "{}\t", record.mapq.unwrap_or(255))?;

        // 6. CIGAR
        if record.cigar.is_empty() {
            write!(self.writer, "*\t")?;
        } else {
            for op in &record.cigar {
                write!(self.writer, "{}", op)?;
            }
            write!(self.writer, "\t")?;
        }

        // 7. RNEXT (mate reference, "=" if same, "*" if unmapped)
        if let Some(mate_ref_id) = record.mate_reference_id {
            if Some(mate_ref_id) == record.reference_id {
                write!(self.writer, "=\t")?;
            } else {
                write!(self.writer, "{}\t", mate_ref_id)?; // Note: Should lookup name
            }
        } else {
            write!(self.writer, "*\t")?;
        }

        // 8. PNEXT (mate position, 0 if unmapped)
        if let Some(mate_pos) = record.mate_position {
            write!(self.writer, "{}\t", mate_pos + 1)?; // BAM is 0-based, SAM is 1-based
        } else {
            write!(self.writer, "0\t")?;
        }

        // 9. TLEN (template length)
        write!(self.writer, "{}\t", record.template_length)?;

        // 10. SEQ (sequence, "*" if not stored)
        if record.sequence.is_empty() {
            write!(self.writer, "*\t")?;
        } else {
            write!(self.writer, "{}\t", String::from_utf8_lossy(&record.sequence))?;
        }

        // 11. QUAL (quality scores, "*" if not stored)
        if record.quality.is_empty() {
            write!(self.writer, "*")?;
        } else {
            // Convert Phred scores to ASCII (add 33)
            for &q in &record.quality {
                write!(self.writer, "{}", (q + 33) as char)?;
            }
        }

        // 12. TAGS (optional fields)
        if !record.tags.is_empty() {
            // Parse and write each tag
            if let Ok(tags) = record.tags.iter() {
                for tag in tags {
                    write!(self.writer, "\t{}", tag)?;
                }
            }
        }

        writeln!(self.writer)?;

        Ok(())
    }

    /// Get a mutable reference to the underlying writer.
    ///
    /// Useful for flushing or accessing the writer directly.
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.writer
    }

    /// Consume the writer and return the underlying writer.
    pub fn into_inner(self) -> W {
        self.writer
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::bam::Tags;

    #[test]
    fn test_write_minimal_record() {
        let mut output = Vec::new();
        let mut sam = SamWriter::new(&mut output);

        let record = Record {
            name: "read1".to_string(),
            reference_id: None,
            position: None,
            mapq: Some(0),
            flags: 4, // Unmapped
            mate_reference_id: None,
            mate_position: None,
            template_length: 0,
            sequence: b"ACGT".to_vec(),
            quality: vec![30, 30, 30, 30],
            cigar: vec![],
            tags: Tags::new(),
        };

        sam.write_record(&record).unwrap();

        let sam_text = String::from_utf8(output).unwrap();
        assert!(sam_text.contains("read1"));
        assert!(sam_text.contains("ACGT"));
    }

    #[test]
    fn test_write_header() {
        let mut output = Vec::new();
        let mut sam = SamWriter::new(&mut output);

        let header = Header {
            text: "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n".to_string(),
            references: vec![],
        };

        sam.write_header(&header).unwrap();

        let sam_text = String::from_utf8(output).unwrap();
        assert!(sam_text.contains("@HD"));
        assert!(sam_text.contains("chr1"));
    }
}
