//! BAM ↔ CAF conversion.
//!
//! This module provides functions to convert between BAM and CAF formats,
//! preserving all alignment information.
//!
//! # Example
//!
//! ```no_run
//! use caf::conversion::bam_to_caf;
//! use std::path::Path;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! bam_to_caf(
//!     Path::new("input.bam"),
//!     Path::new("output.caf"),
//! )?;
//! # Ok(())
//! # }
//! ```

use crate::{
    block::AlignmentRecord,
    io::{CafFileWriter, CafFileReader},
    CafError, Result,
};
use biometal::io::bam::{BamReader, Record as BamRecord, CigarOp};
use std::fs::File;
use std::path::Path;

/// Convert a CAF file to SAM/BAM format.
///
/// Reads a CAF file and converts all records to SAM format,
/// preserving all alignment information.
///
/// # Arguments
///
/// * `caf_path` - Path to input CAF file
/// * `sam_path` - Path to output SAM file (use .bam for BAM output)
///
/// # Errors
///
/// Returns error if:
/// - Input CAF file cannot be read
/// - Output SAM file cannot be created
/// - CAF parsing fails
/// - SAM writing fails
///
/// # Example
///
/// ```no_run
/// use caf::conversion::caf_to_sam;
/// use std::path::Path;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// caf_to_sam(
///     Path::new("input.caf"),
///     Path::new("output.sam"),
/// )?;
/// # Ok(())
/// # }
/// ```
pub fn caf_to_sam<P: AsRef<Path>>(caf_path: P, sam_path: P) -> Result<()> {
    use biometal::io::bam::{SamWriter, Header as BamHeader, Reference};
    use std::io::BufWriter;

    // Open CAF file
    let caf_reader = CafFileReader::open(caf_path.as_ref())?;
    let caf_header = caf_reader.header();

    // Create SAM writer
    let output = File::create(sam_path.as_ref())?;
    let buf_writer = BufWriter::new(output);
    let mut sam_writer = SamWriter::new(buf_writer);

    // Convert header
    let sam_text = String::from_utf8_lossy(&caf_header.sam_header).to_string();

    let references: Vec<Reference> = caf_header.ref_names.iter()
        .zip(caf_header.ref_lengths.iter())
        .map(|(name, length)| Reference {
            name: name.clone(),
            length: *length as u32,
        })
        .collect();

    let bam_header = BamHeader {
        text: sam_text,
        references,
    };

    // Write header
    sam_writer.write_header(&bam_header)
        .map_err(|e| CafError::Other(format!("Failed to write SAM header: {}", e)))?;

    // Convert records
    let mut records_converted = 0u64;

    for result in caf_reader.records()? {
        let caf_record = result?;
        let bam_record = alignment_record_to_bam_record(&caf_record)?;
        sam_writer.write_record(&bam_record)
            .map_err(|e| CafError::Other(format!("Failed to write SAM record: {}", e)))?;
        records_converted += 1;

        if records_converted % 100_000 == 0 {
            log::info!("Converted {} records", records_converted);
        }
    }

    log::info!(
        "Conversion complete: {} records written to {:?}",
        records_converted,
        sam_path.as_ref()
    );

    Ok(())
}

/// Convert a BAM file to CAF format.
///
/// Reads a BAM file and converts all records to CAF format,
/// preserving all alignment information including:
/// - Alignment positions and mapping quality
/// - Sequences and quality scores
/// - CIGAR strings
/// - Read names
/// - Mate pair information
/// - Template lengths
///
/// # Arguments
///
/// * `bam_path` - Path to input BAM file
/// * `caf_path` - Path to output CAF file
///
/// # Errors
///
/// Returns error if:
/// - Input BAM file cannot be read
/// - Output CAF file cannot be created
/// - BAM parsing fails
/// - CAF writing fails
///
/// # Example
///
/// ```no_run
/// use caf::conversion::bam_to_caf;
/// use std::path::Path;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// bam_to_caf(
///     Path::new("input.bam"),
///     Path::new("output.caf"),
/// )?;
/// # Ok(())
/// # }
/// ```
pub fn bam_to_caf<P: AsRef<Path>>(bam_path: P, caf_path: P) -> Result<()> {
    // Phase 1: Collect quality score samples for dictionary training
    log::info!("Phase 1: Collecting quality score samples for dictionary training");
    let mut quality_samples: Vec<Vec<u8>> = Vec::new();
    const SAMPLE_SIZE: usize = 30_000; // Collect 30K records for training

    {
        let mut sample_reader = BamReader::from_path(bam_path.as_ref())
            .map_err(|e| CafError::Other(format!("Failed to open BAM file for sampling: {}", e)))?;

        let mut sample_count = 0;
        while let Some(bam_record) = sample_reader.read_record()
            .map_err(|e| CafError::Other(format!("Failed to read BAM record: {}", e)))?
        {
            if !bam_record.quality.is_empty() {
                quality_samples.push(bam_record.quality.clone());
            }
            sample_count += 1;

            if sample_count >= SAMPLE_SIZE {
                break;
            }
        }
    }

    // Train dictionary if we have samples
    let quality_dict = if !quality_samples.is_empty() {
        log::info!("Training dictionary on {} quality score samples", quality_samples.len());
        let dict = crate::compression::train_dictionary(
            &quality_samples,
            crate::compression::DICTIONARY_SIZE,
        )?;
        log::info!("Dictionary trained: {} bytes", dict.len());
        Some(dict)
    } else {
        None
    };

    // Phase 2: Full conversion with trained dictionary
    log::info!("Phase 2: Converting BAM → CAF");

    // Open BAM file (automatically handles BGZF decompression)
    let mut bam_reader = BamReader::from_path(bam_path.as_ref())
        .map_err(|e| CafError::Other(format!("Failed to open BAM file: {}", e)))?;

    // Create CAF writer
    let mut caf_writer = CafFileWriter::create(caf_path.as_ref())?;

    // Set quality dictionary if trained
    if let Some(ref dict) = quality_dict {
        caf_writer.set_quality_dict(dict.clone())?;
    }

    // Extract header information
    let header = bam_reader.header();

    // Set SAM header in CAF
    caf_writer.set_sam_header(header.text.clone().into_bytes())?;

    // Extract reference sequences
    let ref_names: Vec<String> = header
        .references
        .iter()
        .map(|r| r.name.clone())
        .collect();

    let ref_lengths: Vec<i32> = header
        .references
        .iter()
        .map(|r| r.length as i32)
        .collect();

    if !ref_names.is_empty() {
        caf_writer.set_references(ref_names, ref_lengths)?;
    }

    // Convert records
    let mut records_converted = 0u64;

    while let Some(bam_record) = bam_reader.read_record()
        .map_err(|e| CafError::Other(format!("Failed to read BAM record: {}", e)))?
    {
        let caf_record = bam_record_to_alignment_record(&bam_record)?;
        caf_writer.add_record(caf_record)?;
        records_converted += 1;

        if records_converted % 100_000 == 0 {
            log::info!("Converted {} records", records_converted);
        }
    }

    // Finalize CAF file
    caf_writer.finalize()?;

    log::info!(
        "Conversion complete: {} records written to {:?}",
        records_converted,
        caf_path.as_ref()
    );

    Ok(())
}

/// Convert a single BAM record to CAF AlignmentRecord.
///
/// This function handles all field conversions, including:
/// - Converting Option<usize> to i32 (-1 for None)
/// - Converting Option<i32> to i32 (-1 for None)
/// - Converting Option<u8> to u8 (255 for None per SAM spec)
/// - Converting CigarOp enum to (op_code, length) tuples
/// - Converting String names to Vec<u8>
pub fn bam_record_to_alignment_record(bam: &BamRecord) -> Result<AlignmentRecord> {
    // Convert reference IDs (Option<usize> -> i32, -1 for None)
    let ref_id = match bam.reference_id {
        Some(id) => id as i32,
        None => -1,
    };

    let mate_ref_id = match bam.mate_reference_id {
        Some(id) => id as i32,
        None => -1,
    };

    // Convert positions (Option<i32> -> i32, -1 for None)
    let position = bam.position.unwrap_or(-1);
    let mate_position = bam.mate_position.unwrap_or(-1);

    // Convert mapping quality (Option<u8> -> u8, 255 for None per SAM spec)
    let mapq = bam.mapq.unwrap_or(255);

    // Convert CIGAR operations
    let cigar = bam.cigar.iter().map(cigar_op_to_tuple).collect();

    // Convert read name to bytes (with null terminator for BAM compatibility)
    let mut read_name = bam.name.as_bytes().to_vec();
    if !read_name.is_empty() && read_name[read_name.len() - 1] != 0 {
        read_name.push(0); // Add null terminator
    }

    Ok(AlignmentRecord {
        ref_id,
        position,
        mapq,
        flags: bam.flags,
        sequence: bam.sequence.clone(),
        qualities: bam.quality.clone(),
        cigar,
        read_name,
        mate_ref_id,
        mate_position,
        template_length: bam.template_length,
    })
}

/// Convert CAF AlignmentRecord to BAM Record.
///
/// This is the inverse of `bam_record_to_alignment_record`,
/// converting from CAF's columnar format back to BAM's row format.
fn alignment_record_to_bam_record(caf: &AlignmentRecord) -> Result<BamRecord> {
    use biometal::io::bam::Tags;

    // Convert reference IDs (i32 -> Option<usize>, -1 for unmapped)
    let reference_id = if caf.ref_id >= 0 {
        Some(caf.ref_id as usize)
    } else {
        None
    };

    let mate_reference_id = if caf.mate_ref_id >= 0 {
        Some(caf.mate_ref_id as usize)
    } else {
        None
    };

    // Convert positions (i32 -> Option<i32>, -1 for unmapped)
    let position = if caf.position >= 0 {
        Some(caf.position)
    } else {
        None
    };

    let mate_position = if caf.mate_position >= 0 {
        Some(caf.mate_position)
    } else {
        None
    };

    // Convert mapping quality (u8 -> Option<u8>, 255 for unavailable)
    let mapq = if caf.mapq == 255 {
        None
    } else {
        Some(caf.mapq)
    };

    // Convert CIGAR tuples back to CigarOp enum
    let cigar = caf.cigar.iter().map(|(op, len)| tuple_to_cigar_op(*op, *len)).collect::<Result<Vec<_>>>()?;

    // Convert read name from bytes (strip null terminator)
    let name = {
        let bytes = if caf.read_name.ends_with(&[0]) {
            &caf.read_name[..caf.read_name.len() - 1]
        } else {
            &caf.read_name[..]
        };
        String::from_utf8(bytes.to_vec())
            .map_err(|e| CafError::Other(format!("Invalid read name UTF-8: {}", e)))?
    };

    Ok(BamRecord {
        name,
        reference_id,
        position,
        mapq,
        flags: caf.flags,
        mate_reference_id,
        mate_position,
        template_length: caf.template_length,
        sequence: caf.sequence.clone(),
        quality: caf.qualities.clone(),
        cigar,
        tags: Tags::new(), // TODO: Preserve auxiliary tags
    })
}

/// Convert (op_code, length) tuple to CigarOp.
///
/// CIGAR operation codes (per SAM spec):
/// - 0: Match/Mismatch (M)
/// - 1: Insertion (I)
/// - 2: Deletion (D)
/// - 3: RefSkip (N)
/// - 4: SoftClip (S)
/// - 5: HardClip (H)
/// - 6: Padding (P)
/// - 7: SeqMatch (=)
/// - 8: SeqMismatch (X)
fn tuple_to_cigar_op(op_code: u8, length: u32) -> Result<CigarOp> {
    match op_code {
        0 => Ok(CigarOp::Match(length)),
        1 => Ok(CigarOp::Insertion(length)),
        2 => Ok(CigarOp::Deletion(length)),
        3 => Ok(CigarOp::RefSkip(length)),
        4 => Ok(CigarOp::SoftClip(length)),
        5 => Ok(CigarOp::HardClip(length)),
        6 => Ok(CigarOp::Padding(length)),
        7 => Ok(CigarOp::SeqMatch(length)),
        8 => Ok(CigarOp::SeqMismatch(length)),
        _ => Err(CafError::Other(format!("Invalid CIGAR op code: {}", op_code))),
    }
}

/// Convert a CigarOp to (op_code, length) tuple.
///
/// CIGAR operation codes (per SAM spec):
/// - 0: Match/Mismatch (M)
/// - 1: Insertion (I)
/// - 2: Deletion (D)
/// - 3: RefSkip (N)
/// - 4: SoftClip (S)
/// - 5: HardClip (H)
/// - 6: Padding (P)
/// - 7: SeqMatch (=)
/// - 8: SeqMismatch (X)
fn cigar_op_to_tuple(op: &CigarOp) -> (u8, u32) {
    match op {
        CigarOp::Match(len) => (0, *len),
        CigarOp::Insertion(len) => (1, *len),
        CigarOp::Deletion(len) => (2, *len),
        CigarOp::RefSkip(len) => (3, *len),
        CigarOp::SoftClip(len) => (4, *len),
        CigarOp::HardClip(len) => (5, *len),
        CigarOp::Padding(len) => (6, *len),
        CigarOp::SeqMatch(len) => (7, *len),
        CigarOp::SeqMismatch(len) => (8, *len),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use biometal::io::bam::CigarOp;
    use tempfile::NamedTempFile;

    #[test]
    fn test_cigar_op_to_tuple() {
        assert_eq!(cigar_op_to_tuple(&CigarOp::Match(10)), (0, 10));
        assert_eq!(cigar_op_to_tuple(&CigarOp::Insertion(5)), (1, 5));
        assert_eq!(cigar_op_to_tuple(&CigarOp::Deletion(3)), (2, 3));
        assert_eq!(cigar_op_to_tuple(&CigarOp::RefSkip(100)), (3, 100));
        assert_eq!(cigar_op_to_tuple(&CigarOp::SoftClip(2)), (4, 2));
        assert_eq!(cigar_op_to_tuple(&CigarOp::HardClip(1)), (5, 1));
        assert_eq!(cigar_op_to_tuple(&CigarOp::Padding(4)), (6, 4));
        assert_eq!(cigar_op_to_tuple(&CigarOp::SeqMatch(20)), (7, 20));
        assert_eq!(cigar_op_to_tuple(&CigarOp::SeqMismatch(8)), (8, 8));
    }

    #[test]
    fn test_bam_record_conversion_mapped() {
        use biometal::io::bam::Record;

        let bam_record = Record {
            name: "read1".to_string(),
            reference_id: Some(0),
            position: Some(1000),
            mapq: Some(60),
            flags: 99,
            mate_reference_id: Some(0),
            mate_position: Some(1100),
            template_length: 200,
            sequence: b"ACGT".to_vec(),
            quality: b"IIII".to_vec(),
            cigar: vec![CigarOp::Match(4)],
            tags: Default::default(),
        };

        let caf_record = bam_record_to_alignment_record(&bam_record).unwrap();

        assert_eq!(caf_record.ref_id, 0);
        assert_eq!(caf_record.position, 1000);
        assert_eq!(caf_record.mapq, 60);
        assert_eq!(caf_record.flags, 99);
        assert_eq!(caf_record.mate_ref_id, 0);
        assert_eq!(caf_record.mate_position, 1100);
        assert_eq!(caf_record.template_length, 200);
        assert_eq!(caf_record.sequence, b"ACGT");
        assert_eq!(caf_record.qualities, b"IIII");
        assert_eq!(caf_record.cigar, vec![(0, 4)]);
        assert_eq!(&caf_record.read_name[..5], b"read1");
    }

    #[test]
    fn test_bam_record_conversion_unmapped() {
        use biometal::io::bam::Record;

        let bam_record = Record {
            name: "unmapped".to_string(),
            reference_id: None,
            position: None,
            mapq: None,
            flags: 4, // Unmapped flag
            mate_reference_id: None,
            mate_position: None,
            template_length: 0,
            sequence: b"ACGT".to_vec(),
            quality: b"IIII".to_vec(),
            cigar: vec![],
            tags: Default::default(),
        };

        let caf_record = bam_record_to_alignment_record(&bam_record).unwrap();

        assert_eq!(caf_record.ref_id, -1);
        assert_eq!(caf_record.position, -1);
        assert_eq!(caf_record.mapq, 255); // Unavailable per SAM spec
        assert_eq!(caf_record.flags, 4);
        assert_eq!(caf_record.mate_ref_id, -1);
        assert_eq!(caf_record.mate_position, -1);
    }

    #[test]
    fn test_bam_record_conversion_complex_cigar() {
        use biometal::io::bam::Record;

        let bam_record = Record {
            name: "complex".to_string(),
            reference_id: Some(0),
            position: Some(1000),
            mapq: Some(30),
            flags: 99,
            mate_reference_id: Some(0),
            mate_position: Some(1200),
            template_length: 300,
            sequence: b"ACGTACGTAC".to_vec(),
            quality: b"IIIIIIIIII".to_vec(),
            cigar: vec![
                CigarOp::Match(4),
                CigarOp::Insertion(2),
                CigarOp::Match(4),
            ],
            tags: Default::default(),
        };

        let caf_record = bam_record_to_alignment_record(&bam_record).unwrap();

        assert_eq!(caf_record.cigar, vec![(0, 4), (1, 2), (0, 4)]);
    }
}
