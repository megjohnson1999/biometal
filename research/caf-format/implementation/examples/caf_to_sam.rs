//! Example: Convert CAF file to SAM format.
//!
//! Usage:
//!   cargo run --example caf_to_sam -- input.caf output.sam

use caf::caf_to_sam;
use std::env;
use std::path::Path;
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();

    let args: Vec<String> = env::args().collect();

    if args.len() != 3 {
        eprintln!("Usage: {} <input.caf> <output.sam>", args[0]);
        eprintln!("\nExample:");
        eprintln!("  cargo run --example caf_to_sam -- test.caf output.sam");
        std::process::exit(1);
    }

    let caf_path = Path::new(&args[1]);
    let sam_path = Path::new(&args[2]);

    if !caf_path.exists() {
        eprintln!("Error: Input CAF file not found: {:?}", caf_path);
        std::process::exit(1);
    }

    println!("Converting CAF → SAM");
    println!("  Input:  {:?}", caf_path);
    println!("  Output: {:?}", sam_path);
    println!();

    let start = Instant::now();
    caf_to_sam(caf_path, sam_path)?;
    let elapsed = start.elapsed();

    // Get file sizes
    let caf_size = std::fs::metadata(caf_path)?.len();
    let sam_size = std::fs::metadata(sam_path)?.len();

    println!();
    println!("Conversion complete!");
    println!("  Time:        {:?}", elapsed);
    println!("  CAF size:    {} bytes", caf_size);
    println!("  SAM size:    {} bytes", sam_size);
    println!("  Expansion:   {:.1}%", (sam_size as f64 / caf_size as f64) * 100.0);

    Ok(())
}
