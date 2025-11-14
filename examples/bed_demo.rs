//! BED format parsing examples.
//!
//! This example demonstrates parsing BED3, BED6, and BED12 formats
//! using biometal's streaming parsers.

use biometal::formats::bed::{Bed12Record, Bed3Record, Bed6Record};
use biometal::formats::{TabDelimitedParser, TabDelimitedRecord};
use std::io::Cursor;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== BED Format Parsing Demo ===\n");

    // BED3: Minimal genomic intervals
    println!("--- BED3 (chrom, start, end) ---");
    let bed3_data = "\
# BED3 example: genomic intervals
chr1\t1000\t2000
chr1\t3000\t4000
chr2\t5000\t6000
";

    let parser = TabDelimitedParser::<_, Bed3Record>::new(Cursor::new(bed3_data.as_bytes()));
    for result in parser {
        let record = result?;
        println!(
            "  {}: {}-{} (length: {})",
            record.interval.chrom,
            record.interval.start,
            record.interval.end,
            record.interval.length()
        );
    }

    // BED6: Standard annotations
    println!("\n--- BED6 (BED3 + name, score, strand) ---");
    let bed6_data = "\
# BED6 example: gene annotations
chr1\t1000\t2000\tgene1\t500\t+
chr1\t3000\t4000\tgene2\t800\t-
chr2\t5000\t6000\tgene3\t.\t.
";

    let parser = TabDelimitedParser::<_, Bed6Record>::new(Cursor::new(bed6_data.as_bytes()));
    for result in parser {
        let record = result?;
        println!(
            "  {}: {}-{} | name: {} | score: {} | strand: {}",
            record.bed3.interval.chrom,
            record.bed3.interval.start,
            record.bed3.interval.end,
            record.name.as_deref().unwrap_or("N/A"),
            record
                .score
                .map(|s| s.to_string())
                .unwrap_or_else(|| "N/A".to_string()),
            record
                .strand
                .map(|s| s.to_string())
                .unwrap_or_else(|| "N/A".to_string())
        );
    }

    // BED12: Full gene models
    println!("\n--- BED12 (BED6 + exon structure) ---");
    let bed12_data = "\
# BED12 example: gene models with exons
chr1\t1000\t5000\tTXN1\t100\t+\t1200\t4800\t0,0,255\t3\t400,600,400\t0,1000,4600
chr2\t2000\t8000\tTXN2\t200\t-\t2500\t7500\t255,0,0\t2\t1000,2000\t0,4000
";

    let parser = TabDelimitedParser::<_, Bed12Record>::new(Cursor::new(bed12_data.as_bytes()));
    for result in parser {
        let record = result?;
        println!(
            "  {}: {}-{} | {} | {} exons",
            record.bed6.bed3.interval.chrom,
            record.bed6.bed3.interval.start,
            record.bed6.bed3.interval.end,
            record.bed6.name.as_deref().unwrap_or("N/A"),
            record.block_count.unwrap_or(0)
        );

        // Show exon structure
        if let (Some(sizes), Some(starts)) = (&record.block_sizes, &record.block_starts) {
            for (i, (size, start)) in sizes.iter().zip(starts.iter()).enumerate() {
                let abs_start = record.bed6.bed3.interval.start + (*start as u64);
                let abs_end = abs_start + (*size as u64);
                println!("    Exon {}: {}-{} (size: {})", i + 1, abs_start, abs_end, size);
            }
        }

        // Show CDS region
        if let (Some(thick_start), Some(thick_end)) =
            (record.thick_start, record.thick_end)
        {
            println!("    CDS: {}-{}", thick_start, thick_end);
        }
    }

    // Demonstrate genomic operations
    println!("\n--- Genomic Interval Operations ---");
    let r1 = Bed3Record::from_line("chr1\t1000\t2000")?;
    let r2 = Bed3Record::from_line("chr1\t1500\t2500")?;
    let r3 = Bed3Record::from_line("chr1\t3000\t4000")?;

    println!("  r1: chr1:1000-2000");
    println!("  r2: chr1:1500-2500");
    println!("  r3: chr1:3000-4000");
    println!("  r1 overlaps r2? {}", r1.interval.overlaps(&r2.interval));
    println!("  r1 overlaps r3? {}", r1.interval.overlaps(&r3.interval));
    println!("  r1 contains r2? {}", r1.interval.contains(&r2.interval));

    // Round-trip serialization
    println!("\n--- Round-Trip Serialization ---");
    let original = Bed6Record::from_line("chr1\t1000\t2000\tgene1\t500\t+")?;
    let line = original.to_line();
    let parsed = Bed6Record::from_line(&line)?;
    println!("  Original: chr1:1000-2000, gene1, score=500, strand=+");
    println!("  Serialized: {}", line);
    println!("  Parsed matches original? {}", original == parsed);

    println!("\n✓ All BED format examples completed successfully!");

    Ok(())
}
