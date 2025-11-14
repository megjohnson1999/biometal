//! VCF format parsing examples.
//!
//! This example demonstrates:
//! - Parsing VCF headers and metadata
//! - Variant record parsing
//! - INFO field analysis
//! - Multi-sample genotype data

use biometal::formats::vcf::VcfParser;
use std::io::Cursor;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== VCF Format Parsing Demo ===\n");

    // Sample VCF data with header
    let vcf_data = "\
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership\">
##FILTER=<ID=LowQual,Description=\"Low quality variant\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2
chr1\t12345\trs123\tA\tT\t30.0\tPASS\tDP=100;AF=0.5;DB\tGT:GQ\t0/1:99\t1/1:60
chr1\t23456\t.\tG\tC\t50.0\tPASS\tDP=80;AF=0.25\tGT:GQ\t0/0:99\t0/1:80
chr2\t34567\trs456\tT\tA,G\t40.0\tPASS\tDP=120;AF=0.3,0.2\tGT:GQ\t0/1:90\t1/2:70
chr2\t45678\t.\tC\tG\t10.0\tLowQual\tDP=20;AF=0.1\tGT:GQ\t0/0:50\t0/1:30
";

    let mut parser = VcfParser::new(Cursor::new(vcf_data.as_bytes()));

    // Parse header
    println!("--- VCF Header ---");
    let header = parser.parse_header()?;
    println!("  File format: {}", header.fileformat);
    println!("  INFO fields: {}", header.info_fields.len());
    println!("  FORMAT fields: {}", header.format_fields.len());
    println!("  Filters: {}", header.filters.len());
    println!("  Contigs: {}", header.contigs.len());
    println!("  Samples: {:?}", header.samples);

    // Show INFO field definitions
    println!("\n--- INFO Field Definitions ---");
    for (id, desc) in &header.info_fields {
        println!("  {}: {}", id, desc);
    }

    // Show contig information
    println!("\n--- Contigs ---");
    for (id, length) in &header.contigs {
        if let Some(len) = length {
            println!("  {}: {} bp", id, len);
        } else {
            println!("  {}: length unknown", id);
        }
    }

    // Parse variants
    println!("\n--- Variant Records ---");
    let mut records = Vec::new();
    for result in parser {
        let record = result?;
        records.push(record);
    }

    println!("Total variants: {}\n", records.len());

    // Analyze each variant
    for (i, record) in records.iter().enumerate() {
        println!("Variant {}:", i + 1);
        println!("  Location: {}:{}", record.chrom, record.pos);
        println!(
            "  ID: {}",
            record.id.as_deref().unwrap_or("novel")
        );
        println!("  Ref: {} -> Alt: {}", record.reference, record.alternate.join(","));
        println!(
            "  Quality: {}",
            record
                .quality
                .map(|q| format!("{:.1}", q))
                .unwrap_or_else(|| "N/A".to_string())
        );
        println!(
            "  Filter: {}",
            record.filter.as_deref().unwrap_or("N/A")
        );

        // Parse INFO fields
        if !record.info.is_empty() {
            println!("  INFO:");
            for (key, value) in &record.info {
                if value.is_empty() {
                    println!("    {} (flag)", key);
                } else {
                    println!("    {}={}", key, value);
                }
            }
        }

        // Show genotypes
        if let Some(format) = &record.format {
            println!("  FORMAT: {}", format);
            for (j, genotype) in record.samples.iter().enumerate() {
                println!("    sample{}: {}", j + 1, genotype);
            }
        }
        println!();
    }

    // Variant type analysis
    println!("--- Variant Type Analysis ---");
    let mut snps = 0;
    let mut insertions = 0;
    let mut deletions = 0;
    let mut multi_allelic = 0;

    for record in &records {
        if record.alternate.len() > 1 {
            multi_allelic += 1;
        }

        for alt in &record.alternate {
            let ref_len = record.reference.len();
            let alt_len = alt.len();

            if ref_len == alt_len && ref_len == 1 {
                snps += 1;
            } else if alt_len > ref_len {
                insertions += 1;
            } else if alt_len < ref_len {
                deletions += 1;
            }
        }
    }

    println!("  SNPs: {}", snps);
    println!("  Insertions: {}", insertions);
    println!("  Deletions: {}", deletions);
    println!("  Multi-allelic sites: {}", multi_allelic);

    // Quality statistics
    println!("\n--- Quality Statistics ---");
    let qualities: Vec<f64> = records
        .iter()
        .filter_map(|r| r.quality)
        .collect();

    if !qualities.is_empty() {
        let min_qual = qualities.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_qual = qualities.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let avg_qual = qualities.iter().sum::<f64>() / qualities.len() as f64;

        println!("  Quality range: {:.1} - {:.1}", min_qual, max_qual);
        println!("  Average quality: {:.1}", avg_qual);
    }

    // Filter statistics
    println!("\n--- Filter Statistics ---");
    let mut pass_count = 0;
    let mut filtered_count = 0;

    for record in &records {
        match &record.filter {
            Some(f) if f == "PASS" => pass_count += 1,
            Some(_) => filtered_count += 1,
            None => {}
        }
    }

    println!("  PASS: {}", pass_count);
    println!("  Filtered: {}", filtered_count);

    // Allele frequency analysis
    println!("\n--- Allele Frequency Analysis ---");
    for record in &records {
        if let Some(af_str) = record.info.get("AF") {
            let afs: Vec<f64> = af_str
                .split(',')
                .filter_map(|s| s.parse().ok())
                .collect();

            for (i, af) in afs.iter().enumerate() {
                println!(
                    "  {}:{} {}>{} AF={:.3}",
                    record.chrom,
                    record.pos,
                    record.reference,
                    record.alternate.get(i).unwrap_or(&"?".to_string()),
                    af
                );
            }
        }
    }

    println!("\n✓ VCF parsing complete!");

    Ok(())
}
