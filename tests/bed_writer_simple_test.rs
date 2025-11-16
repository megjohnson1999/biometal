//! Simple BED writer tests to verify basic functionality

use biometal::formats::bed::{Bed3Record, Bed6Record};
use biometal::formats::bed_writer::BedWriter;
use biometal::formats::primitives::{GenomicInterval, Strand, TabDelimitedParser};
use std::io::BufReader;
use std::fs::File;
use tempfile::NamedTempFile;

#[test]
fn test_bed3_basic_roundtrip() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    // Write
    {
        let mut writer = BedWriter::create(path).unwrap();
        let interval = GenomicInterval { chrom: "chr1".to_string(), start: 1000, end: 2000 };
        let record = Bed3Record { interval };
        writer.write_bed3(&record).unwrap();
        writer.finish().unwrap();
    }

    // Read back
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let parser = TabDelimitedParser::<_, Bed3Record>::new(reader);
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 1);
    assert_eq!(records[0].interval.chrom, "chr1");
    assert_eq!(records[0].interval.start, 1000);
    assert_eq!(records[0].interval.end, 2000);
}

#[test]
fn test_bed6_basic_roundtrip() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    // Write
    {
        let mut writer = BedWriter::create(path).unwrap();
        let interval = GenomicInterval { chrom: "chr1".to_string(), start: 1000, end: 2000 };
        let record = Bed6Record {
            bed3: Bed3Record { interval },
            name: Some("feature1".to_string()),
            score: Some(500),
            strand: Some(Strand::Forward),
        };
        writer.write_bed6(&record).unwrap();
        writer.finish().unwrap();
    }

    // Read back
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let parser = TabDelimitedParser::<_, Bed6Record>::new(reader);
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 1);
    assert_eq!(records[0].bed3.interval.chrom, "chr1");
    assert_eq!(records[0].name, Some("feature1".to_string()));
    assert_eq!(records[0].strand, Some(Strand::Forward));
}

#[test]
fn test_bed_writer_validation() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();
    let mut writer = BedWriter::create(path).unwrap();

    // Invalid interval: start >= end
    let interval = GenomicInterval { chrom: "chr1".to_string(), start: 2000, end: 1000 };
    let record = Bed3Record { interval };
    let result = writer.write_bed3(&record);
    assert!(result.is_err());
}
