//! Example: Display CAF file information.
//!
//! Usage:
//!   cargo run --example caf_info -- input.caf

use caf::io::CafFileReader;
use std::env;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();

    if args.len() != 2 {
        eprintln!("Usage: {} <input.caf>", args[0]);
        eprintln!("\nExample:");
        eprintln!("  cargo run --example caf_info -- output.caf");
        std::process::exit(1);
    }

    let caf_path = Path::new(&args[1]);

    if !caf_path.exists() {
        eprintln!("Error: CAF file not found: {:?}", caf_path);
        std::process::exit(1);
    }

    // Open CAF file
    let reader = CafFileReader::open(caf_path)?;

    // Display file information
    println!("CAF File Information");
    println!("===================");
    println!();

    let header = reader.header();
    let version_major = (header.version >> 8) & 0xFF;
    let version_minor = header.version & 0xFF;

    println!("Header:");
    println!("  Version:     {}.{}", version_major, version_minor);
    println!("  Block size:  {}", header.block_size);
    println!("  Total blocks: {}", reader.num_blocks());
    println!("  References:  {}", header.num_refs);
    println!();

    if header.num_refs > 0 {
        println!("References:");
        for (i, (name, length)) in header.ref_names.iter()
            .zip(header.ref_lengths.iter())
            .enumerate()
        {
            println!("  {}: {} ({} bp)", i, name, length);
        }
        println!();
    }

    if !header.sam_header.is_empty() {
        println!("SAM Header:");
        let sam_text = String::from_utf8_lossy(&header.sam_header);
        for line in sam_text.lines().take(5) {
            println!("  {}", line);
        }
        if sam_text.lines().count() > 5 {
            println!("  ... ({} more lines)", sam_text.lines().count() - 5);
        }
        println!();
    }

    println!("Summary:");
    println!("  Total records: {}", reader.total_records());
    println!("  Total blocks:  {}", reader.num_blocks());
    println!();

    // Show block offsets from index
    println!("Block Offsets (from index):");
    let index = reader.index();
    for i in 0..index.num_blocks.min(5) {
        let offset = index.block_offsets[i as usize];
        let meta = &index.block_metadata[i as usize];
        println!("  Block {}: offset={}, compressed_size={}, num_records={}",
                 i, offset, meta.compressed_size, meta.num_records);
    }
    if index.num_blocks > 5 {
        println!("  ... ({} more blocks)", index.num_blocks - 5);
    }
    println!();

    // Count records by iterating
    println!("Verifying records...");
    let mut count = 0;
    for result in reader.records()? {
        result?;
        count += 1;
        if count % 10_000 == 0 {
            print!("\r  Counted {} records...", count);
        }
    }
    println!("\r  Verified {} records total", count);

    Ok(())
}
