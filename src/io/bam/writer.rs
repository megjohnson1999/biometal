//! BAM file writer with BGZF compression.
//!
//! This module provides BAM writing capabilities with:
//! - BGZF (blocked gzip) compression using cloudflare_zlib
//! - Streaming architecture (constant memory)
//! - Full BAM specification compliance
//! - Round-trip compatibility with BamReader
//!
//! # Example
//!
//! ```no_run
//! use biometal::io::bam::{BamWriter, Record, Header, Reference};
//! use std::path::Path;
//!
//! # fn main() -> biometal::Result<()> {
//! // Create header
//! let references = vec![
//!     Reference::new("chr1".to_string(), 248956422),
//!     Reference::new("chr2".to_string(), 242193529),
//! ];
//! let header = Header::new(
//!     "@HD\tVN:1.6\tSO:coordinate\n".to_string(),
//!     references
//! );
//!
//! // Create writer
//! let mut writer = BamWriter::create("output.bam", header)?;
//!
//! // Write records
//! let mut record = Record::new();
//! record.name = "read1".to_string();
//! record.reference_id = Some(0);
//! record.position = Some(1000);
//! record.sequence = b"ACGT".to_vec();
//! record.quality = b"####".to_vec();
//!
//! writer.write_record(&record)?;
//!
//! // Finish writing (flushes buffers, writes EOF marker)
//! writer.finish()?;
//! # Ok(())
//! # }
//! ```

use std::io::{self, Write};
use std::path::Path;
use super::{Header, Record, CigarOp};
use crate::io::compression::CompressedWriter;
use crate::io::sink::DataSink;

/// BAM file writer.
///
/// Writes BAM files with BGZF compression using cloudflare_zlib backend.
/// Maintains constant memory usage regardless of file size.
pub struct BamWriter {
    /// Compressed writer (BGZF format)
    writer: CompressedWriter,
    /// Header (needed for reference lookups)
    header: Header,
    /// Number of records written
    records_written: usize,
}

impl BamWriter {
    /// Create a new BAM writer.
    ///
    /// Writes BAM magic bytes and header before returning.
    ///
    /// # Arguments
    ///
    /// * `path` - Output file path
    /// * `header` - BAM header with references
    ///
    /// # Returns
    ///
    /// Writer ready to accept records.
    ///
    /// # Errors
    ///
    /// Returns error if file cannot be created or header writing fails.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use biometal::io::bam::{BamWriter, Header, Reference};
    /// # fn main() -> biometal::Result<()> {
    /// let header = Header::new(
    ///     "@HD\tVN:1.6\n".to_string(),
    ///     vec![Reference::new("chr1".to_string(), 248956422)]
    /// );
    /// let writer = BamWriter::create("output.bam", header)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn create<P: AsRef<Path>>(path: P, header: Header) -> io::Result<Self> {
        let sink = DataSink::from_path(path);
        let writer = CompressedWriter::new(sink)?;

        let mut bam_writer = Self {
            writer,
            header: header.clone(),
            records_written: 0,
        };

        // Write BAM header
        bam_writer.write_header(&header)?;

        Ok(bam_writer)
    }

    /// Write BAM file header.
    ///
    /// Format:
    /// - Magic: "BAM\1" (4 bytes)
    /// - l_text: SAM header text length (4 bytes, little-endian)
    /// - text: SAM header text (l_text bytes)
    /// - n_ref: Number of references (4 bytes, little-endian)
    /// - For each reference:
    ///   - l_name: Reference name length including null (4 bytes)
    ///   - name: Reference name (null-terminated)
    ///   - l_ref: Reference length (4 bytes)
    fn write_header(&mut self, header: &Header) -> io::Result<()> {
        // Magic bytes
        self.writer.write_all(b"BAM\x01")?;

        // SAM header text
        let text_bytes = header.text.as_bytes();
        self.writer.write_all(&(text_bytes.len() as u32).to_le_bytes())?;
        self.writer.write_all(text_bytes)?;

        // Number of references
        self.writer.write_all(&(header.references.len() as u32).to_le_bytes())?;

        // Write each reference
        for reference in &header.references {
            // Reference name (null-terminated)
            let name_bytes = reference.name.as_bytes();
            let l_name = (name_bytes.len() + 1) as u32; // +1 for null terminator
            self.writer.write_all(&l_name.to_le_bytes())?;
            self.writer.write_all(name_bytes)?;
            self.writer.write_all(&[0])?; // null terminator

            // Reference length
            self.writer.write_all(&reference.length.to_le_bytes())?;
        }

        Ok(())
    }

    /// Write a BAM record.
    ///
    /// Encodes the record in BAM binary format and writes to the compressed stream.
    ///
    /// # Arguments
    ///
    /// * `record` - Record to write
    ///
    /// # Returns
    ///
    /// Ok if record written successfully.
    ///
    /// # Errors
    ///
    /// Returns error if record encoding or writing fails.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use biometal::io::bam::{BamWriter, Record, Header, Reference};
    /// # fn main() -> biometal::Result<()> {
    /// # let header = Header::new("@HD\tVN:1.6\n".to_string(), vec![]);
    /// # let mut writer = BamWriter::create("output.bam", header)?;
    /// let mut record = Record::new();
    /// record.name = "read1".to_string();
    /// record.sequence = b"ACGT".to_vec();
    /// record.quality = b"####".to_vec();
    /// writer.write_record(&record)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        // Encode record to bytes
        let encoded = encode_record(record)?;

        // Write to BGZF stream
        self.writer.write_all(&encoded)?;

        self.records_written += 1;
        Ok(())
    }

    /// Get number of records written.
    pub fn records_written(&self) -> usize {
        self.records_written
    }

    /// Finish writing and close the file.
    ///
    /// This flushes all buffers and writes the BGZF EOF marker.
    /// Must be called to ensure all data is written correctly.
    ///
    /// # Errors
    ///
    /// Returns error if flushing or EOF marker writing fails.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use biometal::io::bam::{BamWriter, Header};
    /// # fn main() -> biometal::Result<()> {
    /// # let header = Header::new("@HD\tVN:1.6\n".to_string(), vec![]);
    /// # let mut writer = BamWriter::create("output.bam", header)?;
    /// // ... write records ...
    /// writer.finish()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn finish(mut self) -> io::Result<()> {
        self.writer.finish()
    }
}

/// Encode a BAM record to bytes.
///
/// BAM record format:
/// - block_size: Size of the record data (excluding this field) (4 bytes)
/// - refID: Reference sequence ID, -1 if unmapped (4 bytes)
/// - pos: 0-based leftmost position, -1 if unmapped (4 bytes)
/// - bin_mq_nl: bin (16 bits) | MAPQ (8 bits) | l_read_name (8 bits) (4 bytes)
/// - flag_nc: FLAG (16 bits) | n_cigar_op (16 bits) (4 bytes)
/// - l_seq: Sequence length (4 bytes)
/// - next_refID: Mate reference ID (4 bytes)
/// - next_pos: Mate position (4 bytes)
/// - tlen: Template length (4 bytes)
/// - read_name: Read name (null-terminated, l_read_name bytes)
/// - cigar: CIGAR operations (n_cigar_op × 4 bytes)
/// - seq: Sequence (4-bit encoding, (l_seq+1)/2 bytes)
/// - qual: Quality scores (l_seq bytes)
/// - tags: Optional tags (variable length)
fn encode_record(record: &Record) -> io::Result<Vec<u8>> {
    let mut buffer = Vec::with_capacity(512); // Pre-allocate reasonable size

    // Calculate sizes
    let read_name_bytes = record.name.as_bytes();
    let l_read_name = (read_name_bytes.len() + 1) as u8; // +1 for null terminator
    let n_cigar_op = record.cigar.len() as u16;
    let l_seq = record.sequence.len() as u32;

    // Encode refID (-1 if unmapped)
    let ref_id: i32 = record.reference_id.map(|id| id as i32).unwrap_or(-1);

    // Encode position (-1 if unmapped)
    let pos: i32 = record.position.unwrap_or(-1);

    // Calculate BAM bin (for now, use bin 0 - proper calculation can be added later)
    let bin: u16 = calculate_bin(record.position.unwrap_or(0),
                                  record.cigar.iter().map(|op| op.reference_length()).sum());

    // Encode bin_mq_nl: bin (16 bits) << 16 | mapq (8 bits) << 8 | l_read_name (8 bits)
    let mapq = record.mapq.unwrap_or(255);
    let bin_mq_nl: u32 = ((bin as u32) << 16) | ((mapq as u32) << 8) | (l_read_name as u32);

    // Encode flag_nc: flag (16 bits) << 16 | n_cigar_op (16 bits)
    let flag_nc: u32 = ((record.flags as u32) << 16) | (n_cigar_op as u32);

    // Encode next_refID and next_pos
    let next_ref_id: i32 = record.mate_reference_id.map(|id| id as i32).unwrap_or(-1);
    let next_pos: i32 = record.mate_position.unwrap_or(-1);

    // Build record data (everything except block_size)
    let mut record_data = Vec::with_capacity(512);

    // Fixed-size fields
    record_data.extend_from_slice(&ref_id.to_le_bytes());
    record_data.extend_from_slice(&pos.to_le_bytes());
    record_data.extend_from_slice(&bin_mq_nl.to_le_bytes());
    record_data.extend_from_slice(&flag_nc.to_le_bytes());
    record_data.extend_from_slice(&l_seq.to_le_bytes());
    record_data.extend_from_slice(&next_ref_id.to_le_bytes());
    record_data.extend_from_slice(&next_pos.to_le_bytes());
    record_data.extend_from_slice(&record.template_length.to_le_bytes());

    // Read name (null-terminated)
    record_data.extend_from_slice(read_name_bytes);
    record_data.push(0); // null terminator

    // CIGAR operations
    for cigar_op in &record.cigar {
        record_data.extend_from_slice(&encode_cigar_op(cigar_op).to_le_bytes());
    }

    // Sequence (4-bit encoding)
    record_data.extend_from_slice(&encode_sequence(&record.sequence));

    // Quality scores
    record_data.extend_from_slice(&record.quality);

    // Tags (already in raw format)
    record_data.extend_from_slice(record.tags.as_raw());

    // Write block_size (length of record_data)
    buffer.extend_from_slice(&(record_data.len() as u32).to_le_bytes());

    // Write record data
    buffer.extend_from_slice(&record_data);

    Ok(buffer)
}

/// Encode a CIGAR operation to BAM format.
///
/// BAM CIGAR format: 32-bit integer with:
/// - Low 4 bits: operation type (0-8)
/// - High 28 bits: operation length
fn encode_cigar_op(op: &CigarOp) -> u32 {
    let (op_code, length) = match op {
        CigarOp::Match(len) => (0, *len),
        CigarOp::Insertion(len) => (1, *len),
        CigarOp::Deletion(len) => (2, *len),
        CigarOp::RefSkip(len) => (3, *len),
        CigarOp::SoftClip(len) => (4, *len),
        CigarOp::HardClip(len) => (5, *len),
        CigarOp::Padding(len) => (6, *len),
        CigarOp::SeqMatch(len) => (7, *len),
        CigarOp::SeqMismatch(len) => (8, *len),
    };

    (length << 4) | op_code
}

/// Encode DNA sequence to BAM 4-bit format.
///
/// BAM stores 2 bases per byte:
/// - High 4 bits: first base
/// - Low 4 bits: second base
///
/// Base encoding:
/// - =: 0, A: 1, C: 2, M: 3, G: 4, R: 5, S: 6, V: 7
/// - T: 8, W: 9, Y: 10, H: 11, K: 12, D: 13, B: 14, N: 15
fn encode_sequence(seq: &[u8]) -> Vec<u8> {
    let mut encoded = Vec::with_capacity((seq.len() + 1) / 2);

    for chunk in seq.chunks(2) {
        let byte = if chunk.len() == 2 {
            (base_to_4bit(chunk[0]) << 4) | base_to_4bit(chunk[1])
        } else {
            (base_to_4bit(chunk[0]) << 4) // Last byte, pad with 0
        };
        encoded.push(byte);
    }

    encoded
}

/// Convert ASCII base to 4-bit BAM encoding.
fn base_to_4bit(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'=' => 0,
        b'A' => 1,
        b'C' => 2,
        b'M' => 3,
        b'G' => 4,
        b'R' => 5,
        b'S' => 6,
        b'V' => 7,
        b'T' => 8,
        b'W' => 9,
        b'Y' => 10,
        b'H' => 11,
        b'K' => 12,
        b'D' => 13,
        b'B' => 14,
        b'N' | _ => 15, // N or unknown
    }
}

/// Calculate BAM bin for a given alignment region.
///
/// Uses hierarchical binning system from BAM specification.
/// For now, returns bin 0 (to be implemented fully later).
fn calculate_bin(start: i32, length: i32) -> u16 {
    let end = start + length;

    // Implement BAM binning according to spec
    // See SAM spec section 4.1.1 for details
    if end > start {
        // Calculate bin using the hierarchical binning formula
        let beg = start as u32;
        let end = (end - 1) as u32;

        if (beg >> 14) == (end >> 14) { return ((beg >> 14) + 4681) as u16; }
        if (beg >> 17) == (end >> 17) { return ((beg >> 17) + 585) as u16; }
        if (beg >> 20) == (end >> 20) { return ((beg >> 20) + 73) as u16; }
        if (beg >> 23) == (end >> 23) { return ((beg >> 23) + 9) as u16; }
        if (beg >> 26) == (end >> 26) { return ((beg >> 26) + 1) as u16; }
    }

    0 // Bin 0 for unmapped or problematic regions
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::bam::{BamReader, Reference, Header};
    use tempfile::NamedTempFile;

    #[test]
    fn test_base_to_4bit() {
        assert_eq!(base_to_4bit(b'A'), 1);
        assert_eq!(base_to_4bit(b'C'), 2);
        assert_eq!(base_to_4bit(b'G'), 4);
        assert_eq!(base_to_4bit(b'T'), 8);
        assert_eq!(base_to_4bit(b'N'), 15);

        // Test lowercase
        assert_eq!(base_to_4bit(b'a'), 1);
        assert_eq!(base_to_4bit(b'c'), 2);
    }

    #[test]
    fn test_encode_sequence() {
        // ACGT = 0x12, 0x48
        let seq = b"ACGT";
        let encoded = encode_sequence(seq);
        assert_eq!(encoded, vec![0x12, 0x48]);

        // ACG (odd length) = 0x12, 0x40
        let seq = b"ACG";
        let encoded = encode_sequence(seq);
        assert_eq!(encoded, vec![0x12, 0x40]);
    }

    #[test]
    fn test_encode_cigar_op() {
        // 100M = (100 << 4) | 0 = 1600
        assert_eq!(encode_cigar_op(&CigarOp::Match(100)), 1600);

        // 5I = (5 << 4) | 1 = 81
        assert_eq!(encode_cigar_op(&CigarOp::Insertion(5)), 81);

        // 10D = (10 << 4) | 2 = 162
        assert_eq!(encode_cigar_op(&CigarOp::Deletion(10)), 162);
    }

    #[test]
    fn test_bam_writer_basic() {
        // Create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.path();

        // Create header
        let references = vec![
            Reference::new("chr1".to_string(), 248956422),
            Reference::new("chr2".to_string(), 242193529),
        ];
        let header = Header::new(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n@SQ\tSN:chr2\tLN:242193529\n".to_string(),
            references
        );

        // Write BAM file
        {
            let mut writer = BamWriter::create(temp_path, header.clone()).unwrap();

            // Create a simple record
            let mut record = Record::new();
            record.name = "read1".to_string();
            record.reference_id = Some(0);
            record.position = Some(1000);
            record.mapq = Some(60);
            record.flags = 0;
            record.sequence = b"ACGTACGTACGT".to_vec();
            record.quality = b"############".to_vec();
            record.cigar = vec![CigarOp::Match(12)];

            writer.write_record(&record).unwrap();

            assert_eq!(writer.records_written(), 1);

            writer.finish().unwrap();
        }

        // Read back and verify
        {
            let mut reader = BamReader::from_path(temp_path).unwrap();

            // Verify header
            assert_eq!(reader.header().references.len(), 2);
            assert_eq!(reader.header().references[0].name, "chr1");
            assert_eq!(reader.header().references[0].length, 248956422);

            // Read first record
            let mut records = reader.records();
            let record = records.next().unwrap().unwrap();

            assert_eq!(record.name, "read1");
            assert_eq!(record.reference_id, Some(0));
            assert_eq!(record.position, Some(1000));
            assert_eq!(record.mapq, Some(60));
            assert_eq!(record.flags, 0);
            assert_eq!(record.sequence, b"ACGTACGTACGT");
            assert_eq!(record.quality, b"############");
            assert_eq!(record.cigar.len(), 1);
            assert_eq!(record.cigar[0], CigarOp::Match(12));

            // Should be no more records
            assert!(records.next().is_none());
        }
    }

    #[test]
    fn test_bam_writer_multiple_records() {
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.path();

        let references = vec![Reference::new("chr1".to_string(), 1000000)];
        let header = Header::new("@HD\tVN:1.6\n".to_string(), references);

        {
            let mut writer = BamWriter::create(temp_path, header).unwrap();

            // Write 3 records
            for i in 0..3 {
                let mut record = Record::new();
                record.name = format!("read{}", i + 1);
                record.reference_id = Some(0);
                record.position = Some(i * 100);
                record.sequence = b"ACGT".to_vec();
                record.quality = b"####".to_vec();

                writer.write_record(&record).unwrap();
            }

            assert_eq!(writer.records_written(), 3);
            writer.finish().unwrap();
        }

        {
            let mut reader = BamReader::from_path(temp_path).unwrap();
            let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();

            assert_eq!(records.len(), 3);
            assert_eq!(records[0].name, "read1");
            assert_eq!(records[1].name, "read2");
            assert_eq!(records[2].name, "read3");
        }
    }

    #[test]
    fn test_bam_writer_unmapped_record() {
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.path();

        let header = Header::new("@HD\tVN:1.6\n".to_string(), vec![]);

        {
            let mut writer = BamWriter::create(temp_path, header).unwrap();

            // Create unmapped record
            let mut record = Record::new();
            record.name = "unmapped".to_string();
            record.reference_id = None;
            record.position = None;
            record.mapq = None;
            record.flags = 4; // unmapped flag
            record.sequence = b"ACGT".to_vec();
            record.quality = b"####".to_vec();

            writer.write_record(&record).unwrap();
            writer.finish().unwrap();
        }

        {
            let mut reader = BamReader::from_path(temp_path).unwrap();
            let record = reader.records().next().unwrap().unwrap();

            assert_eq!(record.name, "unmapped");
            assert_eq!(record.reference_id, None);
            assert_eq!(record.position, None);
            assert_eq!(record.flags, 4);
            assert!(record.is_unmapped());
        }
    }
}
