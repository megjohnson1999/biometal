//! Example: Convert BAM file to CAF format.
//!
//! Usage:
//!   cargo run --example bam_to_caf -- input.bam output.caf

use caf::bam_to_caf;
use std::env;
use std::path::Path;
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();

    let args: Vec<String> = env::args().collect();

    if args.len() != 3 {
        eprintln!("Usage: {} <input.bam> <output.caf>", args[0]);
        eprintln!("\nExample:");
        eprintln!("  cargo run --example bam_to_caf -- test.bam output.caf");
        std::process::exit(1);
    }

    let bam_path = Path::new(&args[1]);
    let caf_path = Path::new(&args[2]);

    if !bam_path.exists() {
        eprintln!("Error: Input BAM file not found: {:?}", bam_path);
        std::process::exit(1);
    }

    println!("Converting BAM → CAF");
    println!("  Input:  {:?}", bam_path);
    println!("  Output: {:?}", caf_path);
    println!();

    let start = Instant::now();
    bam_to_caf(bam_path, caf_path)?;
    let elapsed = start.elapsed();

    // Get file sizes
    let bam_size = std::fs::metadata(bam_path)?.len();
    let caf_size = std::fs::metadata(caf_path)?.len();

    println!();
    println!("Conversion complete!");
    println!("  Time:        {:?}", elapsed);
    println!("  BAM size:    {} bytes", bam_size);
    println!("  CAF size:    {} bytes", caf_size);
    println!("  Compression: {:.1}%", (caf_size as f64 / bam_size as f64) * 100.0);

    Ok(())
}
