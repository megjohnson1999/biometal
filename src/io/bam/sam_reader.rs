//! SAM format reader (text alignment files).
//!
//! Provides streaming reading of SAM (Sequence Alignment/Map) text files.
//! SAM is the text representation of BAM - same data, different encoding.
//!
//! # Design
//!
//! - Streaming reader (constant memory, Rule 5)
//! - Text parsing → BAM record structures
//! - Header and alignment line parsing
//! - Reuses existing BamRecord infrastructure
//!
//! # Example
//!
//! ```no_run
//! use biometal::io::bam::SamReader;
//!
//! # fn main() -> biometal::Result<()> {
//! // Read SAM file
//! let mut sam = SamReader::from_path("input.sam")?;
//!
//! // Access header
//! let header = sam.header();
//! println!("References: {}", header.references.len());
//!
//! // Stream records (constant memory)
//! for record in sam.records() {
//!     let record = record?;
//!     println!("{}: {} bp", record.name, record.sequence.len());
//! }
//! # Ok(())
//! # }
//! ```

use super::{CigarOp, Header, Record, Reference, Tags};
use crate::error::{BiometalError, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// SAM format reader.
///
/// Streams SAM alignment records with constant memory usage.
pub struct SamReader<R: Read> {
    reader: BufReader<R>,
    header: Header,
    line_buffer: String,
    pending_line: Option<String>,
}

impl SamReader<File> {
    /// Open a SAM file from a path.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::bam::SamReader;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut sam = SamReader::from_path("alignments.sam")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        Self::new(file)
    }
}

impl<R: Read> SamReader<R> {
    /// Create a new SAM reader from any Read source.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::io::bam::SamReader;
    /// use std::io::Cursor;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let data = b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n";
    /// let cursor = Cursor::new(data);
    /// let sam = SamReader::new(cursor)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(reader: R) -> Result<Self> {
        let mut sam_reader = SamReader {
            reader: BufReader::new(reader),
            header: Header {
                text: String::new(),
                references: Vec::new(),
            },
            line_buffer: String::new(),
            pending_line: None,
        };

        // Parse header
        sam_reader.parse_header()?;

        Ok(sam_reader)
    }

    /// Get a reference to the SAM header.
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Read the next record from the SAM file.
    ///
    /// Returns `Ok(Some(record))` if a record was read,
    /// `Ok(None)` if EOF was reached,
    /// or `Err` if an error occurred.
    pub fn read_record(&mut self) -> Result<Option<Record>> {
        // Handle the pending line from header parsing
        if let Some(line) = self.pending_line.take() {
            if !line.trim().is_empty() && !line.starts_with('@') {
                return Ok(Some(parse_sam_record(&line, &self.header)?));
            }
        }

        // Read subsequent lines
        loop {
            self.line_buffer.clear();
            let bytes_read = self.reader.read_line(&mut self.line_buffer)?;

            if bytes_read == 0 {
                return Ok(None); // EOF
            }

            let line = self.line_buffer.trim_end();

            // Skip empty lines and header lines
            if line.is_empty() || line.starts_with('@') {
                continue;
            }

            return Ok(Some(parse_sam_record(line, &self.header)?));
        }
    }

    /// Parse SAM header lines (@HD, @SQ, etc.)
    fn parse_header(&mut self) -> Result<()> {
        let mut header_text = String::new();
        let mut references = Vec::new();

        loop {
            self.line_buffer.clear();
            let bytes_read = self.reader.read_line(&mut self.line_buffer)?;

            if bytes_read == 0 {
                break; // EOF
            }

            let line = self.line_buffer.trim_end();

            if line.starts_with('@') {
                // Header line
                header_text.push_str(line);
                header_text.push('\n');

                // Parse @SQ lines for reference sequences
                if line.starts_with("@SQ") {
                    if let Some(reference) = parse_sq_line(line) {
                        references.push(reference);
                    }
                }
            } else {
                // First non-header line - save it for the iterator
                self.pending_line = Some(line.to_string());
                break;
            }
        }

        self.header = Header {
            text: header_text,
            references,
        };

        Ok(())
    }

    /// Create an iterator over SAM records.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::bam::SamReader;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut sam = SamReader::from_path("input.sam")?;
    ///
    /// for record in sam.records() {
    ///     let record = record?;
    ///     println!("{}", record.name);
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn records(&mut self) -> SamRecordIterator<'_, R> {
        SamRecordIterator { reader: self }
    }
}

/// Iterator over SAM records.
pub struct SamRecordIterator<'a, R: Read> {
    reader: &'a mut SamReader<R>,
}

impl<'a, R: Read> Iterator for SamRecordIterator<'a, R> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        // Handle the pending line from header parsing
        if let Some(line) = self.reader.pending_line.take() {
            if !line.trim().is_empty() && !line.starts_with('@') {
                return Some(parse_sam_record(&line, &self.reader.header));
            }
        }

        // Read subsequent lines
        loop {
            self.reader.line_buffer.clear();
            match self.reader.reader.read_line(&mut self.reader.line_buffer) {
                Ok(0) => return None, // EOF
                Ok(_) => {
                    let line = self.reader.line_buffer.trim_end();

                    // Skip empty lines and header lines
                    if line.is_empty() || line.starts_with('@') {
                        continue;
                    }

                    return Some(parse_sam_record(line, &self.reader.header));
                }
                Err(e) => return Some(Err(e.into())),
            }
        }
    }
}

/// Parse a @SQ header line to extract reference information.
///
/// Format: @SQ\tSN:chr1\tLN:248956422
fn parse_sq_line(line: &str) -> Option<Reference> {
    let mut name = None;
    let mut length = None;

    for field in line.split('\t').skip(1) {
        if let Some((key, value)) = field.split_once(':') {
            match key {
                "SN" => name = Some(value.to_string()),
                "LN" => length = value.parse::<u32>().ok(),
                _ => {}
            }
        }
    }

    if let (Some(name), Some(length)) = (name, length) {
        Some(Reference { name, length })
    } else {
        None
    }
}

/// Parse a SAM alignment line into a BamRecord.
///
/// Format: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAGS]
fn parse_sam_record(line: &str, header: &Header) -> Result<Record> {
    let fields: Vec<&str> = line.split('\t').collect();

    if fields.len() < 11 {
        return Err(BiometalError::InvalidSamFormat {
            msg: format!("Expected at least 11 fields, found {}", fields.len()),
        });
    }

    // 1. QNAME
    let name = fields[0].to_string();

    // 2. FLAG
    let flags: u16 = fields[1].parse().map_err(|_| BiometalError::InvalidSamFormat {
        msg: format!("Invalid FLAG: {}", fields[1]),
    })?;

    // 3. RNAME (reference name)
    let reference_id = if fields[2] == "*" {
        None
    } else {
        // Look up reference ID by name
        header
            .references
            .iter()
            .position(|r| r.name == fields[2])
    };

    // 4. POS (1-based in SAM, 0-based in BAM)
    let position: i32 = fields[3].parse().map_err(|_| BiometalError::InvalidSamFormat {
        msg: format!("Invalid POS: {}", fields[3]),
    })?;
    let position = if position == 0 {
        None
    } else {
        Some(position - 1) // Convert to 0-based
    };

    // 5. MAPQ
    let mapq: u8 = fields[4].parse().map_err(|_| BiometalError::InvalidSamFormat {
        msg: format!("Invalid MAPQ: {}", fields[4]),
    })?;
    let mapq = if mapq == 255 { None } else { Some(mapq) };

    // 6. CIGAR
    let cigar = if fields[5] == "*" {
        Vec::new()
    } else {
        parse_cigar(fields[5])?
    };

    // 7. RNEXT (mate reference)
    let mate_reference_id = if fields[6] == "*" {
        None
    } else if fields[6] == "=" {
        reference_id // Same as main reference
    } else {
        header
            .references
            .iter()
            .position(|r| r.name == fields[6])
    };

    // 8. PNEXT (mate position, 1-based in SAM, 0-based in BAM)
    let mate_pos: i32 = fields[7].parse().map_err(|_| BiometalError::InvalidSamFormat {
        msg: format!("Invalid PNEXT: {}", fields[7]),
    })?;
    let mate_position = if mate_pos == 0 {
        None
    } else {
        Some(mate_pos - 1) // Convert to 0-based
    };

    // 9. TLEN
    let template_length: i32 = fields[8].parse().map_err(|_| BiometalError::InvalidSamFormat {
        msg: format!("Invalid TLEN: {}", fields[8]),
    })?;

    // 10. SEQ
    let sequence = if fields[9] == "*" {
        Vec::new()
    } else {
        fields[9].as_bytes().to_vec()
    };

    // 11. QUAL (ASCII Phred+33, convert to Phred)
    let quality = if fields[10] == "*" {
        Vec::new()
    } else {
        fields[10].bytes().map(|b| b.saturating_sub(33)).collect()
    };

    // 12. Optional tags
    let tags = if fields.len() > 11 {
        parse_tags(&fields[11..])?
    } else {
        Tags::new()
    };

    Ok(Record {
        name,
        reference_id,
        position,
        mapq,
        flags,
        mate_reference_id,
        mate_position,
        template_length,
        sequence,
        quality,
        cigar,
        tags,
    })
}

/// Parse CIGAR string (e.g., "10M2I5M").
fn parse_cigar(cigar_str: &str) -> Result<Vec<CigarOp>> {
    let mut cigar = Vec::new();
    let mut num_str = String::new();

    for ch in cigar_str.chars() {
        if ch.is_ascii_digit() {
            num_str.push(ch);
        } else {
            let length: u32 = num_str.parse().map_err(|_| BiometalError::InvalidSamFormat {
                msg: format!("Invalid CIGAR length: {}", num_str),
            })?;

            let op = match ch {
                'M' => CigarOp::Match(length),
                'I' => CigarOp::Insertion(length),
                'D' => CigarOp::Deletion(length),
                'N' => CigarOp::RefSkip(length),
                'S' => CigarOp::SoftClip(length),
                'H' => CigarOp::HardClip(length),
                'P' => CigarOp::Padding(length),
                '=' => CigarOp::SeqMatch(length),
                'X' => CigarOp::SeqMismatch(length),
                _ => {
                    return Err(BiometalError::InvalidSamFormat {
                        msg: format!("Unknown CIGAR op: {}", ch),
                    })
                }
            };

            cigar.push(op);
            num_str.clear();
        }
    }

    Ok(cigar)
}

/// Parse SAM optional tags.
fn parse_tags(_tag_fields: &[&str]) -> Result<Tags> {
    // For now, store tags as raw data
    // Full tag parsing would require implementing the tag parsing logic
    // which is complex (type:value format with various types)

    // Placeholder: create empty tags for now
    // TODO: Implement full tag parsing when needed
    Ok(Tags::new())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_sq_line() {
        let line = "@SQ\tSN:chr1\tLN:248956422";
        let reference = parse_sq_line(line).unwrap();
        assert_eq!(reference.name, "chr1");
        assert_eq!(reference.length, 248956422);
    }

    #[test]
    fn test_parse_cigar() {
        let cigar = parse_cigar("10M2I5M").unwrap();
        assert_eq!(cigar.len(), 3);
        assert_eq!(cigar[0], CigarOp::Match(10));
        assert_eq!(cigar[1], CigarOp::Insertion(2));
        assert_eq!(cigar[2], CigarOp::Match(5));
    }

    #[test]
    fn test_sam_reader_basic() {
        let sam_data = b"@HD\tVN:1.6\tSO:unsorted\n\
                         @SQ\tSN:chr1\tLN:1000\n\
                         read1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\t**********\n";

        let cursor = Cursor::new(sam_data);
        let mut reader = SamReader::new(cursor).unwrap();

        // Check header
        assert_eq!(reader.header().references.len(), 1);
        assert_eq!(reader.header().references[0].name, "chr1");

        // Check record
        let records: Vec<_> = reader.records().collect();
        assert_eq!(records.len(), 1);

        let record = records[0].as_ref().unwrap();
        assert_eq!(record.name, "read1");
        assert_eq!(record.position, Some(99)); // 0-based
        assert_eq!(record.cigar.len(), 1);
        assert_eq!(record.sequence, b"ACGTACGTAC");
    }

    #[test]
    fn test_sam_reader_unmapped() {
        let sam_data = b"@HD\tVN:1.6\n\
                         read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t****\n";

        let cursor = Cursor::new(sam_data);
        let mut reader = SamReader::new(cursor).unwrap();

        let records: Vec<_> = reader.records().collect();
        assert_eq!(records.len(), 1);

        let record = records[0].as_ref().unwrap();
        assert_eq!(record.name, "read1");
        assert_eq!(record.flags, 4); // Unmapped
        assert_eq!(record.reference_id, None);
        assert_eq!(record.position, None);
    }
}
