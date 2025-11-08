//! BAM to SAM conversion example (Phase 4 primitives).
//!
//! Demonstrates:
//! - Option B: Tag parsing with type-safe accessors
//! - Option C: SAM output writer for format conversion
//!
//! # Usage
//!
//! ```bash
//! cargo run --example bam_to_sam experiments/native-bam-implementation/test-data/synthetic_100000.bam output.sam
//! ```

use biometal::io::bam::{BamReader, SamWriter};
use std::env;
use std::fs::File;
use std::io::BufWriter;

fn main() -> biometal::Result<()> {
    // Get file paths from command line
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <input.bam> <output.sam>", args[0]);
        eprintln!("\nExample:");
        eprintln!("  cargo run --example bam_to_sam input.bam output.sam");
        std::process::exit(1);
    }

    let input_path = &args[1];
    let output_path = &args[2];

    println!("Converting {} to {}", input_path, output_path);
    println!();

    // Open BAM file with parallel BGZF decompression (Phase 2)
    let mut bam = BamReader::from_path(input_path)?;

    println!("=== BAM Header ===");
    println!("References: {}", bam.header().reference_count());
    println!();

    // Create SAM writer (Phase 4: Option C)
    let output_file = File::create(output_path)?;
    let writer = BufWriter::new(output_file);
    let mut sam = SamWriter::new(writer);

    // Write SAM header
    sam.write_header(bam.header())?;
    println!("✓ Header written");

    // Convert records: BAM (binary) → SAM (text)
    let mut record_count = 0;
    let mut tag_count = 0;

    for result in bam.records() {
        let record = result?;

        // Phase 4: Option B - Demonstrate tag parsing
        if !record.tags.is_empty() {
            // Parse tags (on-demand, no pre-parsing overhead)
            if let Ok(tags) = record.tags.iter() {
                tag_count += tags.len();

                // Example: Access common tags
                if record_count < 3 {
                    println!("\nRecord {}: {}", record_count + 1, record.name);
                    for tag in &tags {
                        println!("  Tag: {}", tag);
                    }
                }
            }
        }

        // Phase 4: Option C - Write to SAM format
        sam.write_record(&record)?;

        record_count += 1;

        if record_count % 10000 == 0 {
            println!("Converted {} records...", record_count);
        }
    }

    println!("\n=== Conversion Complete ===");
    println!("Total records: {}", record_count);
    println!("Total tags: {}", tag_count);
    println!("Output: {}", output_path);

    println!("\n=== Phase 4 Primitives Demonstrated ===");
    println!("✓ Option B: Tag parsing with type-safe accessors");
    println!("✓ Option C: SAM output writer (BAM→SAM conversion)");
    println!("\nThese are low-level primitives that enable tool building!");

    Ok(())
}
