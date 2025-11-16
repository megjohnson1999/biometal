//! Integration tests for GFF3 writer

use biometal::formats::gff::Gff3Record;
use biometal::formats::gff_writer::Gff3Writer;
use biometal::formats::primitives::{Strand, TabDelimitedParser};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use tempfile::NamedTempFile;

#[test]
fn test_gff3_writer_basic_roundtrip() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    // Write
    {
        let mut writer = Gff3Writer::create(path).unwrap();

        let mut attributes = HashMap::new();
        attributes.insert("ID".to_string(), "gene1".to_string());
        attributes.insert("Name".to_string(), "ABC".to_string());

        let record = Gff3Record {
            seqid: "chr1".to_string(),
            source: "Ensembl".to_string(),
            feature_type: "gene".to_string(),
            start: 1000,
            end: 2000,
            score: Some(100.0),
            strand: Strand::Forward,
            phase: None,
            attributes,
        };

        writer.write_record(&record).unwrap();
        writer.finish().unwrap();
    }

    // Read back
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let parser = TabDelimitedParser::<_, Gff3Record>::new(reader);
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 1);
    assert_eq!(records[0].seqid, "chr1");
    assert_eq!(records[0].source, "Ensembl");
    assert_eq!(records[0].feature_type, "gene");
    assert_eq!(records[0].start, 1000);
    assert_eq!(records[0].end, 2000);
    assert!((records[0].score.unwrap() - 100.0).abs() < 0.01);
    assert_eq!(records[0].strand, Strand::Forward);
    assert_eq!(records[0].phase, None);
    assert_eq!(records[0].attributes.get("ID"), Some(&"gene1".to_string()));
    assert_eq!(records[0].attributes.get("Name"), Some(&"ABC".to_string()));
}

#[test]
fn test_gff3_writer_header() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    // Write with automatic header
    {
        let mut writer = Gff3Writer::create(path).unwrap();

        let record = Gff3Record {
            seqid: "chr1".to_string(),
            source: "test".to_string(),
            feature_type: "gene".to_string(),
            start: 100,
            end: 200,
            score: None,
            strand: Strand::Forward,
            phase: None,
            attributes: HashMap::new(),
        };

        writer.write_record(&record).unwrap();
        writer.finish().unwrap();
    }

    // Verify header was written
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let lines: Vec<_> = reader.lines().collect::<Result<_, _>>().unwrap();

    assert_eq!(lines[0], "##gff-version 3");
}

#[test]
fn test_gff3_writer_optional_fields() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    // Write with all optional fields missing
    {
        let mut writer = Gff3Writer::create(path).unwrap();

        let record = Gff3Record {
            seqid: "chr1".to_string(),
            source: "test".to_string(),
            feature_type: "exon".to_string(),
            start: 100,
            end: 200,
            score: None,
            strand: Strand::Unknown,
            phase: None,
            attributes: HashMap::new(),
        };

        writer.write_record(&record).unwrap();
        writer.finish().unwrap();
    }

    // Read back
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let parser = TabDelimitedParser::<_, Gff3Record>::new(reader);
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 1);
    assert_eq!(records[0].score, None);
    assert_eq!(records[0].strand, Strand::Unknown);
    assert_eq!(records[0].phase, None);
    assert!(records[0].attributes.is_empty());
}

#[test]
fn test_gff3_writer_with_phase() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    // Write CDS with phase
    {
        let mut writer = Gff3Writer::create(path).unwrap();

        for phase in 0..=2 {
            let record = Gff3Record {
                seqid: "chr1".to_string(),
                source: "test".to_string(),
                feature_type: "CDS".to_string(),
                start: 100 + phase as u64 * 100,
                end: 200 + phase as u64 * 100,
                score: None,
                strand: Strand::Forward,
                phase: Some(phase),
                attributes: HashMap::new(),
            };
            writer.write_record(&record).unwrap();
        }

        writer.finish().unwrap();
    }

    // Read back
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let parser = TabDelimitedParser::<_, Gff3Record>::new(reader);
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 3);
    assert_eq!(records[0].phase, Some(0));
    assert_eq!(records[1].phase, Some(1));
    assert_eq!(records[2].phase, Some(2));
}

#[test]
fn test_gff3_writer_validation_empty_seqid() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();
    let mut writer = Gff3Writer::create(path).unwrap();

    let record = Gff3Record {
        seqid: "".to_string(),  // Empty seqid
        source: "test".to_string(),
        feature_type: "gene".to_string(),
        start: 100,
        end: 200,
        score: None,
        strand: Strand::Forward,
        phase: None,
        attributes: HashMap::new(),
    };

    let result = writer.write_record(&record);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("seqid cannot be empty"));
}

#[test]
fn test_gff3_writer_validation_invalid_interval() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();
    let mut writer = Gff3Writer::create(path).unwrap();

    let record = Gff3Record {
        seqid: "chr1".to_string(),
        source: "test".to_string(),
        feature_type: "gene".to_string(),
        start: 200,  // start >= end
        end: 100,
        score: None,
        strand: Strand::Forward,
        phase: None,
        attributes: HashMap::new(),
    };

    let result = writer.write_record(&record);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("Invalid interval"));
}

#[test]
fn test_gff3_writer_validation_invalid_phase() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();
    let mut writer = Gff3Writer::create(path).unwrap();

    let record = Gff3Record {
        seqid: "chr1".to_string(),
        source: "test".to_string(),
        feature_type: "CDS".to_string(),
        start: 100,
        end: 200,
        score: None,
        strand: Strand::Forward,
        phase: Some(3),  // Invalid phase (must be 0, 1, or 2)
        attributes: HashMap::new(),
    };

    let result = writer.write_record(&record);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("Invalid phase"));
}

#[test]
fn test_gff3_writer_multiple_records() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    // Write multiple records
    {
        let mut writer = Gff3Writer::create(path).unwrap();

        for i in 0..10 {
            let mut attributes = HashMap::new();
            attributes.insert("ID".to_string(), format!("gene{}", i));

            let record = Gff3Record {
                seqid: "chr1".to_string(),
                source: "test".to_string(),
                feature_type: "gene".to_string(),
                start: i * 1000,
                end: (i + 1) * 1000,
                score: Some(i as f64 * 10.0),
                strand: Strand::Forward,
                phase: None,
                attributes,
            };

            writer.write_record(&record).unwrap();
        }

        assert_eq!(writer.records_written(), 10);
        writer.finish().unwrap();
    }

    // Read back
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let parser = TabDelimitedParser::<_, Gff3Record>::new(reader);
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 10);
    for i in 0..10 {
        assert_eq!(records[i].start, i as u64 * 1000);
        assert_eq!(records[i].end, (i as u64 + 1) * 1000);
        assert_eq!(records[i].attributes.get("ID"), Some(&format!("gene{}", i)));
    }
}

#[test]
fn test_gff3_writer_records_written_counter() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();
    let mut writer = Gff3Writer::create(path).unwrap();

    assert_eq!(writer.records_written(), 0);

    for i in 0..50 {
        let record = Gff3Record {
            seqid: "chr1".to_string(),
            source: "test".to_string(),
            feature_type: "gene".to_string(),
            start: i * 1000,
            end: (i + 1) * 1000,
            score: None,
            strand: Strand::Forward,
            phase: None,
            attributes: HashMap::new(),
        };
        writer.write_record(&record).unwrap();
    }

    assert_eq!(writer.records_written(), 50);
    writer.finish().unwrap();
}
