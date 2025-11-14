//! Integration tests for VCF format parsing.
//!
//! Tests VCF v4.2 parsing with real-world variant calling scenarios.

use biometal::formats::vcf::{VcfParser, VcfRecord};
use biometal::formats::TabDelimitedRecord;
use std::io::Cursor;

#[test]
fn test_vcf_record_basic() {
    let line = "chr1\t12345\trs123\tA\tT\t30.0\tPASS\tDP=100";
    let record = VcfRecord::from_line(line).unwrap();

    assert_eq!(record.chrom, "chr1");
    assert_eq!(record.pos, 12345);
    assert_eq!(record.id, Some("rs123".to_string()));
    assert_eq!(record.reference, "A");
    assert_eq!(record.alternate, vec!["T"]);
    assert_eq!(record.quality, Some(30.0));
    assert_eq!(record.filter, Some("PASS".to_string()));
}

#[test]
fn test_vcf_record_missing_values() {
    let line = "chr1\t100\t.\tA\tT\t.\t.\t.";
    let record = VcfRecord::from_line(line).unwrap();

    assert_eq!(record.id, None);
    assert_eq!(record.quality, None);
    assert_eq!(record.filter, None);
    assert!(record.info.is_empty());
}

#[test]
fn test_vcf_multiple_alts() {
    let line = "chr1\t100\t.\tA\tT,G,C\t.\t.\t.";
    let record = VcfRecord::from_line(line).unwrap();

    assert_eq!(record.alternate.len(), 3);
    assert_eq!(record.alternate, vec!["T", "G", "C"]);
}

#[test]
fn test_vcf_info_parsing() {
    let line = "chr1\t100\t.\tA\tT\t.\t.\tDP=50;AF=0.25;DB";
    let record = VcfRecord::from_line(line).unwrap();

    assert_eq!(record.info.get("DP"), Some(&"50".to_string()));
    assert_eq!(record.info.get("AF"), Some(&"0.25".to_string()));
    assert_eq!(record.info.get("DB"), Some(&"".to_string()));
}

#[test]
fn test_vcf_info_complex() {
    let line = "chr1\t100\t.\tA\tT\t.\t.\tDP=100;AF=0.5,0.3;AN=2;AC=1;NS=3;DB";
    let record = VcfRecord::from_line(line).unwrap();

    assert_eq!(record.info.get("DP"), Some(&"100".to_string()));
    assert_eq!(record.info.get("AF"), Some(&"0.5,0.3".to_string()));
    assert_eq!(record.info.get("AN"), Some(&"2".to_string()));
    assert_eq!(record.info.get("AC"), Some(&"1".to_string()));
    assert_eq!(record.info.get("DB"), Some(&"".to_string()));
}

#[test]
fn test_vcf_with_genotypes() {
    let line = "chr1\t100\t.\tA\tT\t.\tPASS\tDP=50\tGT:GQ:DP\t0/1:99:30\t1/1:80:20";
    let record = VcfRecord::from_line(line).unwrap();

    assert_eq!(record.format, Some("GT:GQ:DP".to_string()));
    assert_eq!(record.samples.len(), 2);
    assert_eq!(record.samples[0], "0/1:99:30");
    assert_eq!(record.samples[1], "1/1:80:20");
}

#[test]
fn test_vcf_round_trip() {
    let original = "chr1\t12345\trs123\tA\tT\t30\tPASS\tAF=0.5;DP=100";
    let record = VcfRecord::from_line(original).unwrap();
    let output = record.to_line();

    let record2 = VcfRecord::from_line(&output).unwrap();
    assert_eq!(record, record2);
}

#[test]
fn test_vcf_header_parsing() {
    let data = "\
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FILTER=<ID=LowQual,Description=\"Low quality\">
##contig=<ID=chr1,length=248956422>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

    let mut parser = VcfParser::new(Cursor::new(data.as_bytes()));
    let header = parser.parse_header().unwrap();

    assert_eq!(header.fileformat, "VCFv4.2");
    assert_eq!(header.info_fields.len(), 1);
    assert_eq!(header.format_fields.len(), 1);
    assert_eq!(header.filters.len(), 1);
    assert_eq!(header.contigs.len(), 1);
}

#[test]
fn test_vcf_header_samples() {
    let data = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\tsample3
";

    let mut parser = VcfParser::new(Cursor::new(data.as_bytes()));
    let header = parser.parse_header().unwrap();

    assert_eq!(header.samples.len(), 3);
    assert_eq!(header.samples[0], "sample1");
    assert_eq!(header.samples[1], "sample2");
    assert_eq!(header.samples[2], "sample3");
}

#[test]
fn test_vcf_full_parsing() {
    let data = "\
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\trs1\tA\tT\t30\tPASS\tDP=50
chr1\t200\trs2\tG\tC\t40\tPASS\tDP=60
chr2\t300\trs3\tT\tA\t50\tPASS\tDP=70
";

    let mut parser = VcfParser::new(Cursor::new(data.as_bytes()));
    let header = parser.parse_header().unwrap();
    assert_eq!(header.fileformat, "VCFv4.2");

    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();
    assert_eq!(records.len(), 3);
    assert_eq!(records[0].pos, 100);
    assert_eq!(records[1].pos, 200);
    assert_eq!(records[2].pos, 300);
}

#[test]
fn test_vcf_snp_detection() {
    let snp = VcfRecord::from_line("chr1\t100\t.\tA\tT\t.\t.\t.").unwrap();
    assert_eq!(snp.reference.len(), 1);
    assert_eq!(snp.alternate[0].len(), 1);
}

#[test]
fn test_vcf_insertion() {
    let ins = VcfRecord::from_line("chr1\t100\t.\tA\tATG\t.\t.\t.").unwrap();
    assert!(ins.alternate[0].len() > ins.reference.len());
}

#[test]
fn test_vcf_deletion() {
    let del = VcfRecord::from_line("chr1\t100\t.\tATG\tA\t.\t.\t.").unwrap();
    assert!(del.reference.len() > del.alternate[0].len());
}

#[test]
fn test_vcf_multi_allelic() {
    let multi = VcfRecord::from_line("chr1\t100\t.\tA\tT,G,C\t.\t.\t.").unwrap();
    assert_eq!(multi.alternate.len(), 3);
}

#[test]
fn test_vcf_dbsnp_id() {
    let with_id = VcfRecord::from_line("chr1\t100\trs123456\tA\tT\t.\t.\t.").unwrap();
    assert_eq!(with_id.id, Some("rs123456".to_string()));

    let without_id = VcfRecord::from_line("chr1\t100\t.\tA\tT\t.\t.\t.").unwrap();
    assert_eq!(without_id.id, None);
}

#[test]
fn test_vcf_quality_scores() {
    let high_qual = VcfRecord::from_line("chr1\t100\t.\tA\tT\t99.9\t.\t.").unwrap();
    assert_eq!(high_qual.quality, Some(99.9));

    let low_qual = VcfRecord::from_line("chr1\t100\t.\tA\tT\t10.5\t.\t.").unwrap();
    assert_eq!(low_qual.quality, Some(10.5));

    let no_qual = VcfRecord::from_line("chr1\t100\t.\tA\tT\t.\t.\t.").unwrap();
    assert_eq!(no_qual.quality, None);
}

#[test]
fn test_vcf_filters() {
    let pass = VcfRecord::from_line("chr1\t100\t.\tA\tT\t.\tPASS\t.").unwrap();
    assert_eq!(pass.filter, Some("PASS".to_string()));

    let filtered = VcfRecord::from_line("chr1\t100\t.\tA\tT\t.\tLowQual\t.").unwrap();
    assert_eq!(filtered.filter, Some("LowQual".to_string()));

    let multi_filter = VcfRecord::from_line("chr1\t100\t.\tA\tT\t.\tLowQual;HighDP\t.").unwrap();
    assert_eq!(multi_filter.filter, Some("LowQual;HighDP".to_string()));
}

#[test]
fn test_vcf_realistic_1000g() {
    // Simulated 1000 Genomes Project variant
    let data = "\
##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00096\tHG00097
chr20\t14370\trs6054257\tG\tA\t29\tPASS\tAC=2;AF=0.5;AN=4;DP=14\tGT:GQ\t0|0:48\t1|0:48
chr20\t17330\t.\tT\tA\t3\tq10\tAC=1;AF=0.017;AN=6;DP=11\tGT:GQ\t0|0:49\t0|1:3
";

    let mut parser = VcfParser::new(Cursor::new(data.as_bytes()));
    let header = parser.parse_header().unwrap();

    assert_eq!(header.samples.len(), 2);
    assert_eq!(header.samples[0], "HG00096");

    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();
    assert_eq!(records.len(), 2);

    // First variant
    assert_eq!(records[0].id, Some("rs6054257".to_string()));
    assert_eq!(records[0].info.get("AC"), Some(&"2".to_string()));
    assert_eq!(records[0].info.get("AF"), Some(&"0.5".to_string()));

    // Second variant (novel)
    assert_eq!(records[1].id, None);
    assert_eq!(records[1].filter, Some("q10".to_string()));
}

#[test]
fn test_vcf_structural_variant() {
    // Large deletion
    let sv = VcfRecord::from_line(
        "chr1\t10000\t.\tNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\tN\t.\tPASS\tSVTYPE=DEL;SVLEN=-30"
    ).unwrap();

    assert!(sv.reference.len() > sv.alternate[0].len());
    assert_eq!(sv.info.get("SVTYPE"), Some(&"DEL".to_string()));
}

#[test]
fn test_vcf_empty_file() {
    let data = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

    let mut parser = VcfParser::new(Cursor::new(data.as_bytes()));
    let _header = parser.parse_header().unwrap();
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 0);
}

#[test]
fn test_vcf_chromosome_names() {
    let data = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t.\t.\t.
2\t200\t.\tG\tC\t.\t.\t.
chrX\t300\t.\tT\tA\t.\t.\t.
chrM\t400\t.\tC\tG\t.\t.\t.
";

    let mut parser = VcfParser::new(Cursor::new(data.as_bytes()));
    let _header = parser.parse_header().unwrap();
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 4);
    assert_eq!(records[0].chrom, "chr1");
    assert_eq!(records[1].chrom, "2");
    assert_eq!(records[2].chrom, "chrX");
    assert_eq!(records[3].chrom, "chrM");
}

#[test]
fn test_vcf_allele_frequencies() {
    let multi_af = VcfRecord::from_line(
        "chr1\t100\t.\tA\tT,G,C\t.\t.\tAF=0.5,0.3,0.2"
    ).unwrap();

    assert_eq!(multi_af.alternate.len(), 3);
    assert_eq!(multi_af.info.get("AF"), Some(&"0.5,0.3,0.2".to_string()));
}

#[test]
fn test_vcf_genotype_parsing() {
    let gt_line = "chr1\t100\t.\tA\tT\t.\tPASS\tDP=50\tGT:DP:GQ\t0/1:25:99\t1/1:25:80\t0/0:0:60";
    let record = VcfRecord::from_line(gt_line).unwrap();

    assert_eq!(record.samples.len(), 3);
    assert_eq!(record.samples[0], "0/1:25:99"); // Het
    assert_eq!(record.samples[1], "1/1:25:80"); // Hom alt
    assert_eq!(record.samples[2], "0/0:0:60");  // Hom ref
}

#[test]
fn test_vcf_phased_genotypes() {
    let phased = VcfRecord::from_line(
        "chr1\t100\t.\tA\tT\t.\tPASS\tDP=50\tGT\t0|1\t1|0"
    ).unwrap();

    assert_eq!(phased.samples[0], "0|1"); // Phased
    assert_eq!(phased.samples[1], "1|0"); // Phased
}

#[test]
fn test_vcf_missing_genotype() {
    let missing = VcfRecord::from_line(
        "chr1\t100\t.\tA\tT\t.\tPASS\tDP=50\tGT\t./.\t0/1"
    ).unwrap();

    assert_eq!(missing.samples[0], "./.");
    assert_eq!(missing.samples[1], "0/1");
}
