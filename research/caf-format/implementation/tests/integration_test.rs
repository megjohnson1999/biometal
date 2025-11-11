//! Integration tests for CAF format write → read round-trip.
//!
//! These tests validate the complete lifecycle:
//! 1. Create CAF file with CafFileWriter
//! 2. Read back with CafFileReader
//! 3. Verify all records match exactly

use caf::{
    io::{CafFileWriter, CafFileReader},
    block::AlignmentRecord,
};
use tempfile::NamedTempFile;

fn create_test_record(ref_id: i32, position: i32, seq_num: u32) -> AlignmentRecord {
    let read_name = format!("read{}", seq_num);
    AlignmentRecord {
        ref_id,
        position,
        mapq: 60,
        flags: 99,
        sequence: b"ACGTACGTACGT".to_vec(),
        qualities: b"IIIIIIIIIIII".to_vec(),
        cigar: vec![(0, 12)], // 12M
        read_name: read_name.into_bytes(),
        mate_ref_id: ref_id,
        mate_position: position + 100,
        template_length: 200,
    }
}

#[test]
fn test_write_read_single_block() {
    let temp = NamedTempFile::new().unwrap();

    // Write
    let mut writer = CafFileWriter::create(temp.path()).unwrap();
    writer.set_references(
        vec!["chr1".to_string()],
        vec![248956422],
    ).unwrap();

    for i in 0..100 {
        writer.add_record(create_test_record(0, 1000 + i, i as u32)).unwrap();
    }

    writer.finalize().unwrap();

    // Read back
    let reader = CafFileReader::open(temp.path()).unwrap();
    assert_eq!(reader.num_blocks(), 1);
    assert_eq!(reader.total_records(), 100);

    // Verify all records
    let records: Vec<_> = reader.records().unwrap()
        .map(|r| r.unwrap())
        .collect();

    assert_eq!(records.len(), 100);

    for (i, record) in records.iter().enumerate() {
        assert_eq!(record.ref_id, 0);
        assert_eq!(record.position, 1000 + i as i32);
        assert_eq!(record.mapq, 60);
        assert_eq!(record.sequence, b"ACGTACGTACGT");
        assert_eq!(record.read_name, format!("read{}", i).as_bytes());
    }
}

#[test]
fn test_write_read_multiple_blocks() {
    let temp = NamedTempFile::new().unwrap();

    // Write 25,000 records (2.5 blocks)
    let mut writer = CafFileWriter::create(temp.path()).unwrap();
    writer.set_references(
        vec!["chr1".to_string(), "chr2".to_string()],
        vec![248956422, 242193529],
    ).unwrap();

    for i in 0..25_000 {
        let ref_id = (i / 10_000) % 2; // Alternate chromosomes
        writer.add_record(create_test_record(ref_id, 1000 + (i % 10000), i as u32)).unwrap();
    }

    writer.finalize().unwrap();

    // Read back
    let reader = CafFileReader::open(temp.path()).unwrap();
    assert_eq!(reader.num_blocks(), 3); // 10K + 10K + 5K
    assert_eq!(reader.total_records(), 25_000);

    // Verify record count
    let records: Vec<_> = reader.records().unwrap()
        .map(|r| r.unwrap())
        .collect();

    assert_eq!(records.len(), 25_000);

    // Spot check first, middle, and last records
    assert_eq!(records[0].ref_id, 0);
    assert_eq!(records[0].position, 1000);

    assert_eq!(records[12_500].position, 1000 + (12_500 % 10000));

    assert_eq!(records[24_999].position, 1000 + (24_999 % 10000));
}

#[test]
fn test_random_block_access() {
    let temp = NamedTempFile::new().unwrap();

    // Write 15,000 records (1.5 blocks)
    let mut writer = CafFileWriter::create(temp.path()).unwrap();

    for i in 0..15_000 {
        writer.add_record(create_test_record(0, 1000 + i, i as u32)).unwrap();
    }

    writer.finalize().unwrap();

    // Read back with random access
    let mut reader = CafFileReader::open(temp.path()).unwrap();
    assert_eq!(reader.num_blocks(), 2);

    // Access second block directly
    let block_reader = reader.read_block(1).unwrap();
    assert_eq!(block_reader.num_records(), 5_000);

    let first_record = block_reader.get_record(0).unwrap();
    assert_eq!(first_record.position, 1000 + 10_000); // First record of second block

    // Access first block
    let block_reader = reader.read_block(0).unwrap();
    assert_eq!(block_reader.num_records(), 10_000);

    let first_record = block_reader.get_record(0).unwrap();
    assert_eq!(first_record.position, 1000);
}

#[test]
fn test_region_query() {
    let temp = NamedTempFile::new().unwrap();

    // Write records at different positions
    let mut writer = CafFileWriter::create(temp.path()).unwrap();
    writer.set_references(
        vec!["chr1".to_string()],
        vec![248956422],
    ).unwrap();

    // Create records: 1000, 1100, 1200, ..., 10900
    for i in 0..100 {
        writer.add_record(create_test_record(0, 1000 + i * 100, i as u32)).unwrap();
    }

    writer.finalize().unwrap();

    // Query region 2000-5000
    let reader = CafFileReader::open(temp.path()).unwrap();
    let records: Vec<_> = reader
        .query_region(0, 2000, 5000)
        .unwrap()
        .map(|r| r.unwrap())
        .collect();

    // Verify all records in range
    for record in &records {
        assert_eq!(record.ref_id, 0);
        assert!(record.position >= 2000);
        assert!(record.position < 5000);
    }

    // Should get positions: 2000, 2100, 2200, ..., 4900
    assert_eq!(records.len(), 30); // (5000 - 2000) / 100
}

#[test]
fn test_header_metadata() {
    let temp = NamedTempFile::new().unwrap();

    // Write with references
    let mut writer = CafFileWriter::create(temp.path()).unwrap();
    writer.set_references(
        vec!["chr1".to_string(), "chr2".to_string(), "chrX".to_string()],
        vec![248956422, 242193529, 156040895],
    ).unwrap();

    writer.add_record(create_test_record(0, 1000, 0)).unwrap();
    writer.finalize().unwrap();

    // Read back and verify header
    let reader = CafFileReader::open(temp.path()).unwrap();
    let header = reader.header();

    assert_eq!(header.num_refs, 3);
    assert_eq!(header.ref_names, vec!["chr1", "chr2", "chrX"]);
    assert_eq!(header.ref_lengths, vec![248956422, 242193529, 156040895]);
}

#[test]
fn test_sam_header() {
    let temp = NamedTempFile::new().unwrap();

    // Write with SAM header
    let mut writer = CafFileWriter::create(temp.path()).unwrap();
    let sam_header = b"@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:248956422\n".to_vec();
    writer.set_sam_header(sam_header.clone()).unwrap();

    writer.add_record(create_test_record(0, 1000, 0)).unwrap();
    writer.finalize().unwrap();

    // Read back and verify
    let reader = CafFileReader::open(temp.path()).unwrap();
    let header = reader.header();

    assert_eq!(header.sam_header, sam_header);
}

#[test]
fn test_empty_file() {
    let temp = NamedTempFile::new().unwrap();

    // Write empty file (no records)
    let mut writer = CafFileWriter::create(temp.path()).unwrap();
    writer.finalize().unwrap();

    // Read back
    let reader = CafFileReader::open(temp.path()).unwrap();
    assert_eq!(reader.num_blocks(), 0);
    assert_eq!(reader.total_records(), 0);

    let records: Vec<_> = reader.records().unwrap().collect();
    assert_eq!(records.len(), 0);
}

#[test]
fn test_diverse_records() {
    let temp = NamedTempFile::new().unwrap();

    // Write diverse records with different field values
    let mut writer = CafFileWriter::create(temp.path()).unwrap();

    // Record 1: Simple
    writer.add_record(AlignmentRecord {
        ref_id: 0,
        position: 1000,
        mapq: 60,
        flags: 99,
        sequence: b"ACGT".to_vec(),
        qualities: b"IIII".to_vec(),
        cigar: vec![(0, 4)],
        read_name: b"read1".to_vec(),
        mate_ref_id: 0,
        mate_position: 1100,
        template_length: 200,
    }).unwrap();

    // Record 2: Complex CIGAR
    writer.add_record(AlignmentRecord {
        ref_id: 1,
        position: 5000,
        mapq: 30,
        flags: 147,
        sequence: b"ACGTACGTACGT".to_vec(),
        qualities: b"HHHHHHHHHHHH".to_vec(),
        cigar: vec![(0, 4), (1, 2), (0, 6)], // 4M2I6M
        read_name: b"read2".to_vec(),
        mate_ref_id: 1,
        mate_position: 5200,
        template_length: 300,
    }).unwrap();

    // Record 3: Long sequence
    writer.add_record(AlignmentRecord {
        ref_id: 0,
        position: 2000,
        mapq: 0,
        flags: 4,
        sequence: b"ACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        qualities: b"IIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec(),
        cigar: vec![(0, 28)],
        read_name: b"long_read_name_test".to_vec(),
        mate_ref_id: -1,
        mate_position: 0,
        template_length: 0,
    }).unwrap();

    writer.finalize().unwrap();

    // Read back and verify
    let reader = CafFileReader::open(temp.path()).unwrap();
    let records: Vec<_> = reader.records().unwrap()
        .map(|r| r.unwrap())
        .collect();

    assert_eq!(records.len(), 3);

    // Verify record 1
    assert_eq!(records[0].ref_id, 0);
    assert_eq!(records[0].position, 1000);
    assert_eq!(records[0].sequence, b"ACGT");
    assert_eq!(records[0].cigar, vec![(0, 4)]);

    // Verify record 2
    assert_eq!(records[1].ref_id, 1);
    assert_eq!(records[1].position, 5000);
    assert_eq!(records[1].cigar, vec![(0, 4), (1, 2), (0, 6)]);

    // Verify record 3
    assert_eq!(records[2].ref_id, 0);
    assert_eq!(records[2].sequence.len(), 28);
    assert_eq!(records[2].read_name, b"long_read_name_test");
}
