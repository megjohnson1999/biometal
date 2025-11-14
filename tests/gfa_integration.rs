//! Integration tests for GFA format parsing.
//!
//! Tests GFA v1.0 parsing with real-world assembly graph scenarios.

use biometal::formats::gfa::{
    GfaHeader, GfaLink, GfaParser, GfaPath, GfaRecord, GfaSegment, Orientation,
};
use biometal::formats::TabDelimitedRecord;
use std::io::Cursor;

#[test]
fn test_gfa_segment_basic() {
    let line = "S\tcontig1\tACGTACGT";
    let seg = GfaSegment::from_line(line).unwrap();

    assert_eq!(seg.name, "contig1");
    assert_eq!(seg.sequence, "ACGTACGT");
    assert_eq!(seg.length(), 8);
}

#[test]
fn test_gfa_segment_with_tags() {
    let line = "S\tcontig1\tACGT\tLN:i:4\tRC:i:100\tFC:i:25";
    let seg = GfaSegment::from_line(line).unwrap();

    assert_eq!(seg.name, "contig1");
    assert_eq!(seg.sequence, "ACGT");
    assert!(seg.tags.contains_key("LN"));
    assert!(seg.tags.contains_key("RC"));
    assert!(seg.tags.contains_key("FC"));
}

#[test]
fn test_gfa_segment_absent_sequence() {
    let line = "S\tcontig1\t*\tLN:i:1000";
    let seg = GfaSegment::from_line(line).unwrap();

    assert_eq!(seg.sequence, "*");
    assert_eq!(seg.length(), 1000); // From LN tag
}

#[test]
fn test_gfa_link_forward_forward() {
    let line = "L\tcontig1\t+\tcontig2\t+\t4M";
    let link = GfaLink::from_line(line).unwrap();

    assert_eq!(link.from_segment, "contig1");
    assert_eq!(link.from_orient, Orientation::Forward);
    assert_eq!(link.to_segment, "contig2");
    assert_eq!(link.to_orient, Orientation::Forward);
    assert_eq!(link.overlap, "4M");
}

#[test]
fn test_gfa_link_reverse_orientations() {
    let line = "L\tcontig1\t-\tcontig2\t+\t3M";
    let link = GfaLink::from_line(line).unwrap();

    assert_eq!(link.from_orient, Orientation::Reverse);
    assert_eq!(link.to_orient, Orientation::Forward);
}

#[test]
fn test_gfa_link_no_overlap() {
    let line = "L\tcontig1\t+\tcontig2\t+\t*";
    let link = GfaLink::from_line(line).unwrap();

    assert_eq!(link.overlap, "*");
}

#[test]
fn test_gfa_path_basic() {
    let line = "P\tpath1\tcontig1+,contig2-,contig3+\t4M,5M";
    let path = GfaPath::from_line(line).unwrap();

    assert_eq!(path.name, "path1");
    assert_eq!(path.segments.len(), 3);
    assert_eq!(path.segments[0], "contig1+");
    assert_eq!(path.segments[1], "contig2-");
    assert_eq!(path.segments[2], "contig3+");
    assert_eq!(path.overlaps.len(), 2);
}

#[test]
fn test_gfa_path_no_overlaps() {
    let line = "P\tpath1\tcontig1+,contig2+\t*";
    let path = GfaPath::from_line(line).unwrap();

    assert_eq!(path.segments.len(), 2);
    assert_eq!(path.overlaps.len(), 0);
}

#[test]
fn test_gfa_header() {
    let line = "H\tVN:Z:1.0";
    let header = GfaHeader::from_line(line).unwrap();

    assert!(header.tags.contains_key("VN"));
}

#[test]
fn test_gfa_record_parsing() {
    let data = "\
H\tVN:Z:1.0
S\tcontig1\tACGT
L\tcontig1\t+\tcontig2\t+\t4M
P\tpath1\tcontig1+,contig2+\t4M
";

    let parser = GfaParser::new(Cursor::new(data.as_bytes()));
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 4);
    assert!(matches!(records[0], GfaRecord::Header(_)));
    assert!(matches!(records[1], GfaRecord::Segment(_)));
    assert!(matches!(records[2], GfaRecord::Link(_)));
    assert!(matches!(records[3], GfaRecord::Path(_)));
}

#[test]
fn test_gfa_simple_graph() {
    let data = "\
S\tA\tACGT
S\tB\tTGCA
L\tA\t+\tB\t+\t4M
";

    let parser = GfaParser::new(Cursor::new(data.as_bytes()));
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 3);

    let segments: Vec<_> = records
        .iter()
        .filter_map(|r| match r {
            GfaRecord::Segment(s) => Some(s),
            _ => None,
        })
        .collect();

    let links: Vec<_> = records
        .iter()
        .filter_map(|r| match r {
            GfaRecord::Link(l) => Some(l),
            _ => None,
        })
        .collect();

    assert_eq!(segments.len(), 2);
    assert_eq!(links.len(), 1);
    assert_eq!(links[0].from_segment, "A");
    assert_eq!(links[0].to_segment, "B");
}

#[test]
fn test_gfa_linear_path() {
    let data = "\
S\tctg1\tACGT
S\tctg2\tTGCA
S\tctg3\tGGAT
L\tctg1\t+\tctg2\t+\t1M
L\tctg2\t+\tctg3\t+\t1M
P\tlinear\tctg1+,ctg2+,ctg3+\t1M,1M
";

    let parser = GfaParser::new(Cursor::new(data.as_bytes()));
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    let paths: Vec<_> = records
        .iter()
        .filter_map(|r| match r {
            GfaRecord::Path(p) => Some(p),
            _ => None,
        })
        .collect();

    assert_eq!(paths.len(), 1);
    assert_eq!(paths[0].segments.len(), 3);
    assert_eq!(paths[0].overlaps.len(), 2);
}

#[test]
fn test_gfa_branching_graph() {
    let data = "\
S\tA\tACGT
S\tB\tTGCA
S\tC\tGGAT
S\tD\tCCAT
L\tA\t+\tB\t+\t1M
L\tA\t+\tC\t+\t1M
L\tB\t+\tD\t+\t1M
L\tC\t+\tD\t+\t1M
";

    let parser = GfaParser::new(Cursor::new(data.as_bytes()));
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    let links: Vec<_> = records
        .iter()
        .filter_map(|r| match r {
            GfaRecord::Link(l) => Some(l),
            _ => None,
        })
        .collect();

    // A has 2 outgoing edges (branch point)
    let a_outgoing = links
        .iter()
        .filter(|l| l.from_segment == "A")
        .count();
    assert_eq!(a_outgoing, 2);

    // D has 2 incoming edges (merge point)
    let d_incoming = links
        .iter()
        .filter(|l| l.to_segment == "D")
        .count();
    assert_eq!(d_incoming, 2);
}

#[test]
fn test_gfa_reverse_complement_edges() {
    let data = "\
S\tA\tACGT
S\tB\tTGCA
L\tA\t+\tB\t-\t2M
L\tB\t+\tA\t-\t2M
";

    let parser = GfaParser::new(Cursor::new(data.as_bytes()));
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    let links: Vec<_> = records
        .iter()
        .filter_map(|r| match r {
            GfaRecord::Link(l) => Some(l),
            _ => None,
        })
        .collect();

    assert_eq!(links.len(), 2);
    assert_eq!(links[0].to_orient, Orientation::Reverse);
    assert_eq!(links[1].to_orient, Orientation::Reverse);
}

#[test]
fn test_gfa_multiple_paths() {
    let data = "\
S\tA\tACGT
S\tB\tTGCA
S\tC\tGGAT
L\tA\t+\tB\t+\t1M
L\tA\t+\tC\t+\t1M
P\tpath1\tA+,B+\t1M
P\tpath2\tA+,C+\t1M
";

    let parser = GfaParser::new(Cursor::new(data.as_bytes()));
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    let paths: Vec<_> = records
        .iter()
        .filter_map(|r| match r {
            GfaRecord::Path(p) => Some(p),
            _ => None,
        })
        .collect();

    assert_eq!(paths.len(), 2);
    assert_eq!(paths[0].name, "path1");
    assert_eq!(paths[1].name, "path2");
}

#[test]
fn test_gfa_realistic_assembly() {
    // Simulated SPAdes-like assembly graph
    let data = "\
H\tVN:Z:1.0
S\tNODE_1\tACGTACGTACGTACGT\tLN:i:16\tKC:i:100
S\tNODE_2\tTGCAGCTGCAGCTGCA\tLN:i:16\tKC:i:150
S\tNODE_3\tGGATCCGATCCGATCC\tLN:i:16\tKC:i:80
S\tNODE_4\tCCATGGCATGCATGCA\tLN:i:16\tKC:i:120
L\tNODE_1\t+\tNODE_2\t+\t4M\tMQ:i:60
L\tNODE_2\t+\tNODE_3\t+\t3M\tMQ:i:55
L\tNODE_1\t+\tNODE_4\t-\t2M\tMQ:i:40
";

    let parser = GfaParser::new(Cursor::new(data.as_bytes()));
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    // Validate header
    let headers: Vec<_> = records
        .iter()
        .filter_map(|r| match r {
            GfaRecord::Header(h) => Some(h),
            _ => None,
        })
        .collect();
    assert_eq!(headers.len(), 1);

    // Validate segments have expected tags
    let segments: Vec<_> = records
        .iter()
        .filter_map(|r| match r {
            GfaRecord::Segment(s) => Some(s),
            _ => None,
        })
        .collect();

    for seg in &segments {
        assert!(seg.tags.contains_key("LN"));
        assert!(seg.tags.contains_key("KC"));
    }

    // Validate links have quality scores
    let links: Vec<_> = records
        .iter()
        .filter_map(|r| match r {
            GfaRecord::Link(l) => Some(l),
            _ => None,
        })
        .collect();

    for link in &links {
        assert!(link.tags.contains_key("MQ"));
    }
}

#[test]
fn test_gfa_round_trip_segment() {
    let original = "S\tcontig1\tACGT\tLN:i:4";
    let seg = GfaSegment::from_line(original).unwrap();
    let output = seg.to_line();
    assert_eq!(output, original);
}

#[test]
fn test_gfa_round_trip_link() {
    let original = "L\tA\t+\tB\t-\t4M";
    let link = GfaLink::from_line(original).unwrap();
    let output = link.to_line();
    assert_eq!(output, original);
}

#[test]
fn test_gfa_round_trip_path() {
    let original = "P\tpath1\tA+,B-,C+\t4M,3M";
    let path = GfaPath::from_line(original).unwrap();
    let output = path.to_line();
    assert_eq!(output, original);
}

#[test]
fn test_gfa_empty_graph() {
    let data = "H\tVN:Z:1.0\n";

    let parser = GfaParser::new(Cursor::new(data.as_bytes()));
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 1);
    assert!(matches!(records[0], GfaRecord::Header(_)));
}

#[test]
fn test_gfa_comments_ignored() {
    let data = "\
# This is a comment
H\tVN:Z:1.0
# Another comment
S\tA\tACGT
# More comments
L\tA\t+\tB\t+\t1M
";

    let parser = GfaParser::new(Cursor::new(data.as_bytes()));
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    // Comments should be filtered out
    assert_eq!(records.len(), 3);
}

#[test]
fn test_gfa_segment_name_parsing() {
    // Various naming conventions
    let names = vec!["contig1", "NODE_123", "utg000001l", "edge_42_length_1000"];

    for name in names {
        let line = format!("S\t{}\tACGT", name);
        let seg = GfaSegment::from_line(&line).unwrap();
        assert_eq!(seg.name, name);
    }
}

#[test]
fn test_gfa_long_sequence() {
    let long_seq = "A".repeat(10000);
    let line = format!("S\tlong_contig\t{}", long_seq);
    let seg = GfaSegment::from_line(&line).unwrap();

    assert_eq!(seg.name, "long_contig");
    assert_eq!(seg.length(), 10000);
}

#[test]
fn test_gfa_complex_overlap() {
    // Complex CIGAR strings
    let overlaps = vec!["10M", "5M1I5M", "3M1D7M", "*"];

    for overlap in overlaps {
        let line = format!("L\tA\t+\tB\t+\t{}", overlap);
        let link = GfaLink::from_line(&line).unwrap();
        assert_eq!(link.overlap, overlap);
    }
}
