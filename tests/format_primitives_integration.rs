//! Integration tests for format primitives.
//!
//! These tests demonstrate how the various primitives work together
//! to parse tab-delimited bioinformatics formats.

use biometal::formats::primitives::{
    fields::{parse_attributes, parse_optional, parse_required, split_fields},
    GenomicInterval, HeaderParser, Result, Strand, TabDelimitedParser, TabDelimitedRecord,
};
use std::str::FromStr;

/// Example GFF-like record using all primitives.
#[derive(Debug, PartialEq)]
struct GffLikeRecord {
    interval: GenomicInterval,
    source: String,
    feature_type: String,
    score: Option<f64>,
    strand: Strand,
    attributes: std::collections::HashMap<String, String>,
}

impl TabDelimitedRecord for GffLikeRecord {
    fn from_line(line: &str) -> Result<Self> {
        // GFF format: seqid source type start end score strand phase attributes
        let fields = split_fields(line, Some(9), 0)?;

        let chrom = fields[0].to_string();
        let start: u64 = parse_required(fields[3], "start", 0)?;
        let end: u64 = parse_required(fields[4], "end", 0)?;

        let interval = GenomicInterval::new(chrom, start - 1, end)?; // GFF is 1-based
        let source = fields[1].to_string();
        let feature_type = fields[2].to_string();
        let score = parse_optional(fields[5], "score", 0)?;
        let strand = Strand::from_str(fields[6])?;
        let attributes = parse_attributes(fields[8], 0)?;

        Ok(GffLikeRecord {
            interval,
            source,
            feature_type,
            score,
            strand,
            attributes,
        })
    }

    fn to_line(&self) -> String {
        let attrs_str = self
            .attributes
            .iter()
            .map(|(k, v)| {
                if v.is_empty() {
                    k.clone()
                } else {
                    format!("{}={}", k, v)
                }
            })
            .collect::<Vec<_>>()
            .join(";");

        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t.\t{}",
            self.interval.chrom,
            self.source,
            self.feature_type,
            self.interval.start + 1, // Convert back to 1-based
            self.interval.end,
            self.score
                .map(|s| s.to_string())
                .unwrap_or_else(|| ".".to_string()),
            self.strand,
            attrs_str
        )
    }
}

#[test]
fn test_full_integration() {
    let gff_data = "\
##gff-version 3
##sequence-region chr1 1 1000000
chr1\tEnsembl\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=ABC1;biotype=protein_coding
chr1\tEnsembl\texon\t1000\t1500\t100.5\t+\t.\tID=exon1;Parent=gene1
chr2\tEnsembl\tgene\t3000\t4000\t.\t-\t.\tID=gene2;Name=XYZ1
";

    // 1. Parse headers
    let header_parser = HeaderParser::new(gff_data.as_bytes());
    let (headers, first_line, reader) = header_parser.parse_headers_with_first_line().unwrap();

    assert_eq!(headers.len(), 2);
    assert!(headers[0].contains("gff-version"));
    assert!(headers[1].contains("sequence-region"));

    // 2. Parse first data line manually
    let first_record = GffLikeRecord::from_line(&first_line.unwrap()).unwrap();
    assert_eq!(first_record.interval.chrom, "chr1");
    assert_eq!(first_record.interval.start, 999); // 0-based
    assert_eq!(first_record.interval.end, 2000);
    assert_eq!(first_record.strand, Strand::Forward);
    assert_eq!(first_record.attributes.get("ID"), Some(&"gene1".to_string()));
    assert_eq!(
        first_record.attributes.get("Name"),
        Some(&"ABC1".to_string())
    );

    // 3. Parse remaining records with TabDelimitedParser
    let parser = TabDelimitedParser::<_, GffLikeRecord>::new(reader);
    let records: Vec<_> = parser.collect::<Result<_>>().unwrap();

    assert_eq!(records.len(), 2);

    // Verify exon record
    assert_eq!(records[0].feature_type, "exon");
    assert_eq!(records[0].score, Some(100.5));
    assert_eq!(
        records[0].attributes.get("Parent"),
        Some(&"gene1".to_string())
    );

    // Verify second gene
    assert_eq!(records[1].interval.chrom, "chr2");
    assert_eq!(records[1].strand, Strand::Reverse);
    assert_eq!(records[1].score, None);
}

#[test]
fn test_round_trip() {
    let original = GffLikeRecord {
        interval: GenomicInterval::new("chr1".to_string(), 999, 2000).unwrap(),
        source: "test".to_string(),
        feature_type: "gene".to_string(),
        score: Some(42.5),
        strand: Strand::Forward,
        attributes: {
            let mut map = std::collections::HashMap::new();
            map.insert("ID".to_string(), "gene1".to_string());
            map.insert("Name".to_string(), "TEST".to_string());
            map
        },
    };

    let line = original.to_line();
    let parsed = GffLikeRecord::from_line(&line).unwrap();

    assert_eq!(parsed.interval, original.interval);
    assert_eq!(parsed.source, original.source);
    assert_eq!(parsed.feature_type, original.feature_type);
    assert_eq!(parsed.score, original.score);
    assert_eq!(parsed.strand, original.strand);
    assert_eq!(parsed.attributes.get("ID"), original.attributes.get("ID"));
    assert_eq!(
        parsed.attributes.get("Name"),
        original.attributes.get("Name")
    );
}

#[test]
fn test_parse_optional_fields() {
    // Test record with missing score
    let line = "chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1";
    let record = GffLikeRecord::from_line(line).unwrap();

    assert_eq!(record.score, None);
    assert_eq!(record.strand, Strand::Forward);
}

#[test]
fn test_genomic_interval_operations() {
    let record1 = GffLikeRecord {
        interval: GenomicInterval::new("chr1".to_string(), 100, 200).unwrap(),
        source: "test".to_string(),
        feature_type: "gene".to_string(),
        score: None,
        strand: Strand::Forward,
        attributes: std::collections::HashMap::new(),
    };

    let record2 = GffLikeRecord {
        interval: GenomicInterval::new("chr1".to_string(), 150, 250).unwrap(),
        source: "test".to_string(),
        feature_type: "exon".to_string(),
        score: None,
        strand: Strand::Forward,
        attributes: std::collections::HashMap::new(),
    };

    // Test interval operations
    assert!(record1.interval.overlaps(&record2.interval));
    assert_eq!(record1.interval.length(), 100);
}
