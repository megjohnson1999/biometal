//! Integration tests for GFF3 format parsing.
//!
//! Tests GFF3 gene annotation parsing with hierarchical features.

use biometal::formats::gff::{Gff3Parser, Gff3Record};
use biometal::formats::{Strand, TabDelimitedRecord};
use std::io::Cursor;

#[test]
fn test_gff3_basic() {
    let line = "chr1\tEnsembl\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=ABC1";
    let record = Gff3Record::from_line(line).unwrap();

    assert_eq!(record.seqid, "chr1");
    assert_eq!(record.source, "Ensembl");
    assert_eq!(record.feature_type, "gene");
    assert_eq!(record.start, 1000);
    assert_eq!(record.end, 2000);
    assert_eq!(record.strand, Strand::Forward);
    assert_eq!(record.get_id(), Some("gene1"));
}

#[test]
fn test_gff3_missing_values() {
    let line = "chr1\t.\tgene\t1000\t2000\t.\t.\t.\t.";
    let record = Gff3Record::from_line(line).unwrap();

    assert_eq!(record.source, ".");
    assert_eq!(record.score, None);
    assert_eq!(record.strand, Strand::Unknown);
    assert_eq!(record.phase, None);
    assert!(record.attributes.is_empty());
}

#[test]
fn test_gff3_cds_phase() {
    let phases = vec![
        ("chr1\t.\tCDS\t1000\t2000\t.\t+\t0\tID=cds1", Some(0)),
        ("chr1\t.\tCDS\t1000\t2000\t.\t+\t1\tID=cds1", Some(1)),
        ("chr1\t.\tCDS\t1000\t2000\t.\t+\t2\tID=cds1", Some(2)),
        ("chr1\t.\tCDS\t1000\t2000\t.\t+\t.\tID=cds1", None),
    ];

    for (line, expected_phase) in phases {
        let record = Gff3Record::from_line(line).unwrap();
        assert_eq!(record.phase, expected_phase);
    }
}

#[test]
fn test_gff3_attributes() {
    let line = "chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=ABC1;biotype=protein_coding;description=Test gene";
    let record = Gff3Record::from_line(line).unwrap();

    assert_eq!(record.get_id(), Some("gene1"));
    assert_eq!(record.attributes.get("Name"), Some(&"ABC1".to_string()));
    assert_eq!(record.attributes.get("biotype"), Some(&"protein_coding".to_string()));
    assert_eq!(record.attributes.get("description"), Some(&"Test gene".to_string()));
}

#[test]
fn test_gff3_parent_relationship() {
    let gene = Gff3Record::from_line("chr1\t.\tgene\t1000\t5000\t.\t+\t.\tID=gene1").unwrap();
    let mrna = Gff3Record::from_line("chr1\t.\tmRNA\t1000\t5000\t.\t+\t.\tID=mRNA1;Parent=gene1").unwrap();
    let exon = Gff3Record::from_line("chr1\t.\texon\t1000\t1500\t.\t+\t.\tID=exon1;Parent=mRNA1").unwrap();

    assert_eq!(gene.get_id(), Some("gene1"));
    assert_eq!(gene.get_parent(), None);

    assert_eq!(mrna.get_id(), Some("mRNA1"));
    assert_eq!(mrna.get_parent(), Some("gene1"));

    assert_eq!(exon.get_id(), Some("exon1"));
    assert_eq!(exon.get_parent(), Some("mRNA1"));
}

#[test]
fn test_gff3_interval_conversion() {
    let record = Gff3Record::from_line("chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1").unwrap();
    let interval = record.interval().unwrap();

    // GFF3: 1-based inclusive [1000, 2000]
    // Converted to 0-based half-open [999, 2000)
    assert_eq!(interval.start, 999);
    assert_eq!(interval.end, 2000);
    assert_eq!(interval.length(), 1001);
}

#[test]
fn test_gff3_feature_length() {
    let gene = Gff3Record::from_line("chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1").unwrap();
    assert_eq!(gene.length(), 1001); // 1000-2000 inclusive = 1001 bp

    let exon = Gff3Record::from_line("chr1\t.\texon\t1000\t1500\t.\t+\t.\tID=exon1").unwrap();
    assert_eq!(exon.length(), 501);
}

#[test]
fn test_gff3_full_parsing() {
    let data = "\
##gff-version 3
chr1\tEnsembl\tgene\t1000\t2000\t.\t+\t.\tID=gene1
chr1\tEnsembl\tmRNA\t1000\t2000\t.\t+\t.\tID=mRNA1;Parent=gene1
chr1\tEnsembl\texon\t1000\t1500\t.\t+\t.\tID=exon1;Parent=mRNA1
chr1\tEnsembl\texon\t1800\t2000\t.\t+\t.\tID=exon2;Parent=mRNA1
";

    let parser = Gff3Parser::new(Cursor::new(data.as_bytes()));
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 4);
    assert_eq!(records[0].feature_type, "gene");
    assert_eq!(records[1].feature_type, "mRNA");
    assert_eq!(records[2].feature_type, "exon");
    assert_eq!(records[3].feature_type, "exon");
}

#[test]
fn test_gff3_feature_types() {
    let types = vec![
        "gene", "mRNA", "exon", "CDS", "five_prime_UTR", "three_prime_UTR",
        "start_codon", "stop_codon", "tRNA", "ncRNA", "intron"
    ];

    for feature_type in types {
        let line = format!("chr1\t.\t{}\t1000\t2000\t.\t+\t.\tID=feat1", feature_type);
        let record = Gff3Record::from_line(&line).unwrap();
        assert_eq!(record.feature_type, feature_type);
    }
}

#[test]
fn test_gff3_strand_values() {
    let forward = Gff3Record::from_line("chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1").unwrap();
    assert_eq!(forward.strand, Strand::Forward);

    let reverse = Gff3Record::from_line("chr1\t.\tgene\t1000\t2000\t.\t-\t.\tID=gene1").unwrap();
    assert_eq!(reverse.strand, Strand::Reverse);

    let unknown = Gff3Record::from_line("chr1\t.\tgene\t1000\t2000\t.\t.\t.\tID=gene1").unwrap();
    assert_eq!(unknown.strand, Strand::Unknown);
}

#[test]
fn test_gff3_score_values() {
    let with_score = Gff3Record::from_line("chr1\t.\tgene\t1000\t2000\t100.5\t+\t.\tID=gene1").unwrap();
    assert_eq!(with_score.score, Some(100.5));

    let without_score = Gff3Record::from_line("chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1").unwrap();
    assert_eq!(without_score.score, None);
}

#[test]
fn test_gff3_hierarchical_gene() {
    let data = "\
##gff-version 3
chr1\t.\tgene\t1000\t5000\t.\t+\t.\tID=gene1;Name=ABC1
chr1\t.\tmRNA\t1000\t5000\t.\t+\t.\tID=mRNA1;Parent=gene1
chr1\t.\texon\t1000\t1500\t.\t+\t.\tID=exon1;Parent=mRNA1
chr1\t.\texon\t2000\t2500\t.\t+\t.\tID=exon2;Parent=mRNA1
chr1\t.\texon\t4500\t5000\t.\t+\t.\tID=exon3;Parent=mRNA1
chr1\t.\tCDS\t1200\t1500\t.\t+\t0\tID=cds1;Parent=mRNA1
chr1\t.\tCDS\t2000\t2500\t.\t+\t2\tID=cds2;Parent=mRNA1
chr1\t.\tCDS\t4500\t4800\t.\t+\t2\tID=cds3;Parent=mRNA1
";

    let parser = Gff3Parser::new(Cursor::new(data.as_bytes()));
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    // Count features by type
    let genes = records.iter().filter(|r| r.feature_type == "gene").count();
    let mrnas = records.iter().filter(|r| r.feature_type == "mRNA").count();
    let exons = records.iter().filter(|r| r.feature_type == "exon").count();
    let cds = records.iter().filter(|r| r.feature_type == "CDS").count();

    assert_eq!(genes, 1);
    assert_eq!(mrnas, 1);
    assert_eq!(exons, 3);
    assert_eq!(cds, 3);

    // Verify parent-child relationships
    let mrna_rec = records.iter().find(|r| r.feature_type == "mRNA").unwrap();
    assert_eq!(mrna_rec.get_parent(), Some("gene1"));

    let exon_recs: Vec<_> = records.iter().filter(|r| r.feature_type == "exon").collect();
    for exon in exon_recs {
        assert_eq!(exon.get_parent(), Some("mRNA1"));
    }
}

#[test]
fn test_gff3_multiple_transcripts() {
    let data = "\
##gff-version 3
chr1\t.\tgene\t1000\t5000\t.\t+\t.\tID=gene1
chr1\t.\tmRNA\t1000\t5000\t.\t+\t.\tID=mRNA1;Parent=gene1
chr1\t.\tmRNA\t1000\t4500\t.\t+\t.\tID=mRNA2;Parent=gene1
chr1\t.\texon\t1000\t1500\t.\t+\t.\tID=exon1;Parent=mRNA1
chr1\t.\texon\t4500\t5000\t.\t+\t.\tID=exon2;Parent=mRNA1
chr1\t.\texon\t1000\t1500\t.\t+\t.\tID=exon3;Parent=mRNA2
chr1\t.\texon\t4000\t4500\t.\t+\t.\tID=exon4;Parent=mRNA2
";

    let parser = Gff3Parser::new(Cursor::new(data.as_bytes()));
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    let mrnas = records.iter().filter(|r| r.feature_type == "mRNA").count();
    assert_eq!(mrnas, 2); // Two isoforms

    // Check both transcripts have gene1 as parent
    let mrna_recs: Vec<_> = records.iter().filter(|r| r.feature_type == "mRNA").collect();
    for mrna in mrna_recs {
        assert_eq!(mrna.get_parent(), Some("gene1"));
    }
}

#[test]
fn test_gff3_realistic_ensembl() {
    // Simulated Ensembl gene annotation
    let data = "\
##gff-version 3
##sequence-region chr7 1 159345973
chr7\tEnsembl\tgene\t127471196\t127495720\t.\t+\t.\tID=ENSG00000157764;Name=BRAF;biotype=protein_coding
chr7\tEnsembl\tmRNA\t127471196\t127495720\t.\t+\t.\tID=ENST00000288602;Parent=ENSG00000157764;Name=BRAF-201
chr7\tEnsembl\texon\t127471196\t127472090\t.\t+\t.\tID=ENSE00001484009;Parent=ENST00000288602
chr7\tEnsembl\texon\t127481592\t127481744\t.\t+\t.\tID=ENSE00003507205;Parent=ENST00000288602
chr7\tEnsembl\tCDS\t127471196\t127472090\t.\t+\t0\tID=ENSP00000288602;Parent=ENST00000288602
chr7\tEnsembl\tCDS\t127481592\t127481744\t.\t+\t1\tID=ENSP00000288602;Parent=ENST00000288602
";

    let parser = Gff3Parser::new(Cursor::new(data.as_bytes()));
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 6);

    // Validate gene
    let gene = &records[0];
    assert_eq!(gene.feature_type, "gene");
    assert_eq!(gene.get_id(), Some("ENSG00000157764"));
    assert_eq!(gene.attributes.get("Name"), Some(&"BRAF".to_string()));
    assert_eq!(gene.attributes.get("biotype"), Some(&"protein_coding".to_string()));

    // Validate mRNA
    let mrna = &records[1];
    assert_eq!(mrna.feature_type, "mRNA");
    assert_eq!(mrna.get_parent(), Some("ENSG00000157764"));
}

#[test]
fn test_gff3_round_trip() {
    let original = "chr1\tEnsembl\tgene\t1000\t2000\t50.5\t+\t.\tID=gene1;Name=ABC1";
    let record = Gff3Record::from_line(original).unwrap();
    let output = record.to_line();

    let record2 = Gff3Record::from_line(&output).unwrap();
    assert_eq!(record, record2);
}

#[test]
fn test_gff3_chromosome_names() {
    let chroms = vec!["chr1", "1", "chrX", "chrM", "scaffold_123"];

    for chrom in chroms {
        let line = format!("{}\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1", chrom);
        let record = Gff3Record::from_line(&line).unwrap();
        assert_eq!(record.seqid, chrom);
    }
}

#[test]
fn test_gff3_source_fields() {
    let sources = vec!["Ensembl", "NCBI", "RefSeq", "GENCODE", ".", "custom_pipeline"];

    for source in sources {
        let line = format!("chr1\t{}\tgene\t1000\t2000\t.\t+\t.\tID=gene1", source);
        let record = Gff3Record::from_line(&line).unwrap();
        assert_eq!(record.source, source);
    }
}

#[test]
fn test_gff3_empty_file() {
    let data = "##gff-version 3\n";

    let parser = Gff3Parser::new(Cursor::new(data.as_bytes()));
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 0);
}

#[test]
fn test_gff3_comments_ignored() {
    let data = "\
##gff-version 3
# This is a comment
chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1
# Another comment
chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=mRNA1;Parent=gene1
";

    let parser = Gff3Parser::new(Cursor::new(data.as_bytes()));
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 2);
}

#[test]
fn test_gff3_utr_features() {
    let data = "\
##gff-version 3
chr1\t.\tfive_prime_UTR\t1000\t1199\t.\t+\t.\tID=utr5_1;Parent=mRNA1
chr1\t.\tCDS\t1200\t2000\t.\t+\t0\tID=cds1;Parent=mRNA1
chr1\t.\tthree_prime_UTR\t2001\t2500\t.\t+\t.\tID=utr3_1;Parent=mRNA1
";

    let parser = Gff3Parser::new(Cursor::new(data.as_bytes()));
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 3);
    assert_eq!(records[0].feature_type, "five_prime_UTR");
    assert_eq!(records[1].feature_type, "CDS");
    assert_eq!(records[2].feature_type, "three_prime_UTR");
}

#[test]
fn test_gff3_noncoding_rna() {
    let data = "\
##gff-version 3
chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;biotype=lincRNA
chr1\t.\tncRNA\t1000\t2000\t.\t+\t.\tID=ncRNA1;Parent=gene1
chr1\t.\texon\t1000\t1500\t.\t+\t.\tID=exon1;Parent=ncRNA1
chr1\t.\texon\t1800\t2000\t.\t+\t.\tID=exon2;Parent=ncRNA1
";

    let parser = Gff3Parser::new(Cursor::new(data.as_bytes()));
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    let gene = &records[0];
    assert_eq!(gene.attributes.get("biotype"), Some(&"lincRNA".to_string()));

    let ncrna = &records[1];
    assert_eq!(ncrna.feature_type, "ncRNA");
}

#[test]
fn test_gff3_negative_strand_gene() {
    let data = "\
##gff-version 3
chr2\t.\tgene\t10000\t15000\t.\t-\t.\tID=gene2
chr2\t.\tmRNA\t10000\t15000\t.\t-\t.\tID=mRNA2;Parent=gene2
chr2\t.\texon\t10000\t11000\t.\t-\t.\tID=exon1;Parent=mRNA2
chr2\t.\tCDS\t10200\t11000\t.\t-\t0\tID=cds1;Parent=mRNA2
";

    let parser = Gff3Parser::new(Cursor::new(data.as_bytes()));
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    for record in &records {
        assert_eq!(record.strand, Strand::Reverse);
    }
}
