//! GFF3 format parsing examples.
//!
//! This example demonstrates:
//! - Parsing GFF3 gene annotations
//! - Hierarchical feature relationships (gene -> mRNA -> exon)
//! - Attribute parsing (ID, Parent, Name)
//! - Gene structure analysis

use biometal::formats::gff::Gff3Parser;
use std::collections::HashMap;
use std::io::Cursor;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== GFF3 Format Parsing Demo ===\n");

    // Sample GFF3 data with hierarchical gene structure
    let gff_data = "\
##gff-version 3
##sequence-region chr1 1 248956422
chr1\tEnsembl\tgene\t1000\t5000\t.\t+\t.\tID=gene1;Name=ABC1;biotype=protein_coding
chr1\tEnsembl\tmRNA\t1000\t5000\t.\t+\t.\tID=mRNA1;Parent=gene1;Name=ABC1-201
chr1\tEnsembl\texon\t1000\t1500\t.\t+\t.\tID=exon1;Parent=mRNA1
chr1\tEnsembl\texon\t2000\t2500\t.\t+\t.\tID=exon2;Parent=mRNA1
chr1\tEnsembl\texon\t4500\t5000\t.\t+\t.\tID=exon3;Parent=mRNA1
chr1\tEnsembl\tCDS\t1200\t1500\t.\t+\t0\tID=cds1;Parent=mRNA1
chr1\tEnsembl\tCDS\t2000\t2500\t.\t+\t2\tID=cds2;Parent=mRNA1
chr1\tEnsembl\tCDS\t4500\t4800\t.\t+\t2\tID=cds3;Parent=mRNA1
chr2\tEnsembl\tgene\t10000\t15000\t.\t-\t.\tID=gene2;Name=XYZ1;biotype=protein_coding
chr2\tEnsembl\tmRNA\t10000\t15000\t.\t-\t.\tID=mRNA2;Parent=gene2;Name=XYZ1-201
chr2\tEnsembl\texon\t10000\t11000\t.\t-\t.\tID=exon4;Parent=mRNA2
chr2\tEnsembl\texon\t13000\t15000\t.\t-\t.\tID=exon5;Parent=mRNA2
chr2\tEnsembl\tCDS\t10200\t11000\t.\t-\t0\tID=cds4;Parent=mRNA2
chr2\tEnsembl\tCDS\t13000\t14500\t.\t-\t0\tID=cds5;Parent=mRNA2
";

    let parser = Gff3Parser::new(Cursor::new(gff_data.as_bytes()));

    // Parse all records
    let mut genes = Vec::new();
    let mut mrnas = Vec::new();
    let mut exons = Vec::new();
    let mut cds = Vec::new();

    for result in parser {
        let record = result?;
        match record.feature_type.as_str() {
            "gene" => genes.push(record),
            "mRNA" => mrnas.push(record),
            "exon" => exons.push(record),
            "CDS" => cds.push(record),
            _ => {}
        }
    }

    // Summary statistics
    println!("--- Feature Summary ---");
    println!("  Genes: {}", genes.len());
    println!("  mRNAs: {}", mrnas.len());
    println!("  Exons: {}", exons.len());
    println!("  CDS: {}", cds.len());

    // Analyze each gene
    println!("\n--- Gene Structure Analysis ---");
    for gene in &genes {
        let gene_id = gene.get_id().unwrap();
        let gene_name = gene.attributes.get("Name").unwrap();

        println!("\nGene: {} ({})", gene_name, gene_id);
        println!("  Location: {}:{}-{} ({})",
                 gene.seqid, gene.start, gene.end, gene.strand);
        println!("  Length: {} bp", gene.length());

        // Find mRNAs for this gene
        let gene_mrnas: Vec<_> = mrnas
            .iter()
            .filter(|m| m.get_parent() == Some(gene_id))
            .collect();

        println!("  Transcripts: {}", gene_mrnas.len());

        for mrna in &gene_mrnas {
            let mrna_id = mrna.get_id().unwrap();
            let mrna_name = mrna.attributes.get("Name").unwrap();

            println!("\n  Transcript: {} ({})", mrna_name, mrna_id);

            // Find exons for this mRNA
            let mrna_exons: Vec<_> = exons
                .iter()
                .filter(|e| e.get_parent() == Some(mrna_id))
                .collect();

            println!("    Exons: {}", mrna_exons.len());
            for (i, exon) in mrna_exons.iter().enumerate() {
                println!("      Exon {}: {}-{} ({} bp)",
                         i + 1, exon.start, exon.end, exon.length());
            }

            // Find CDS for this mRNA
            let mrna_cds: Vec<_> = cds
                .iter()
                .filter(|c| c.get_parent() == Some(mrna_id))
                .collect();

            println!("    CDS regions: {}", mrna_cds.len());
            let total_cds_length: u64 = mrna_cds.iter().map(|c| c.length()).sum();
            println!("    Total CDS length: {} bp", total_cds_length);
            println!("    Coding AA length: {} aa", total_cds_length / 3);

            for (i, cds_region) in mrna_cds.iter().enumerate() {
                println!("      CDS {}: {}-{} (phase: {})",
                         i + 1,
                         cds_region.start,
                         cds_region.end,
                         cds_region.phase.unwrap_or(0));
            }

            // Calculate UTRs
            let first_cds_start = mrna_cds.iter().map(|c| c.start).min().unwrap();
            let last_cds_end = mrna_cds.iter().map(|c| c.end).max().unwrap();

            if mrna.strand.to_string() == "+" {
                let utr5_length = first_cds_start - mrna.start;
                let utr3_length = mrna.end - last_cds_end;
                println!("    5' UTR: {} bp", utr5_length);
                println!("    3' UTR: {} bp", utr3_length);
            } else {
                let utr3_length = first_cds_start - mrna.start;
                let utr5_length = mrna.end - last_cds_end;
                println!("    5' UTR: {} bp", utr5_length);
                println!("    3' UTR: {} bp", utr3_length);
            }
        }
    }

    // Chromosome distribution
    println!("\n--- Chromosome Distribution ---");
    let mut chr_counts: HashMap<String, usize> = HashMap::new();
    for gene in &genes {
        *chr_counts.entry(gene.seqid.clone()).or_insert(0) += 1;
    }

    for (chr, count) in &chr_counts {
        println!("  {}: {} genes", chr, count);
    }

    // Strand analysis
    println!("\n--- Strand Distribution ---");
    let forward = genes.iter().filter(|g| g.strand.to_string() == "+").count();
    let reverse = genes.iter().filter(|g| g.strand.to_string() == "-").count();
    println!("  Forward (+): {}", forward);
    println!("  Reverse (-): {}", reverse);

    // Biotype analysis
    println!("\n--- Biotype Distribution ---");
    let mut biotype_counts: HashMap<String, usize> = HashMap::new();
    for gene in &genes {
        if let Some(biotype) = gene.attributes.get("biotype") {
            *biotype_counts.entry(biotype.clone()).or_insert(0) += 1;
        }
    }

    for (biotype, count) in &biotype_counts {
        println!("  {}: {}", biotype, count);
    }

    // Feature density
    println!("\n--- Feature Density ---");
    for gene in &genes {
        let gene_mrnas: Vec<_> = mrnas
            .iter()
            .filter(|m| {
                m.get_parent() == Some(gene.get_id().unwrap())
            })
            .collect();

        if let Some(mrna) = gene_mrnas.first() {
            let mrna_id = mrna.get_id().unwrap();
            let exon_count = exons
                .iter()
                .filter(|e| e.get_parent() == Some(mrna_id))
                .count();

            let gene_name = gene.attributes.get("Name").unwrap();
            println!("  {}: {} exons, {:.1} kb",
                     gene_name, exon_count, gene.length() as f64 / 1000.0);
        }
    }

    println!("\n✓ GFF3 analysis complete!");

    Ok(())
}
