//! Simple CAF file size analyzer.
//!
//! Usage:
//!   cargo run --example analyze_caf_simple -- input.caf

use std::env;
use std::fs;
use std::io::{BufReader, Read, Seek, SeekFrom};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();

    if args.len() != 2 {
        eprintln!("Usage: {} <input.caf>", args[0]);
        std::process::exit(1);
    }

    let caf_path = &args[1];
    let file_size = fs::metadata(caf_path)?.len();

    println!("CAF File Size Analysis");
    println!("======================");
    println!("File: {}", caf_path);
    println!("Total size: {:.2} MB ({} bytes)", file_size as f64 / (1024.0 * 1024.0), file_size);
    println!();

    // Read magic (4 bytes)
    let file = fs::File::open(caf_path)?;
    let mut reader = BufReader::new(file);

    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;
    println!("Magic: {:?} ({} bytes)", magic, 4);

    // Try to estimate header size by reading until we hit blocks
    // Header is bincode-serialized, so we need to deserialize it
    let header_start_pos = reader.stream_position()?;

    // Skip to end to read footer
    reader.seek(SeekFrom::End(-32))?;
    let footer_pos = reader.stream_position()?;

    println!();
    println!("Estimated Structure:");
    println!("  Magic: 4 bytes");
    println!("  Header: starts at offset {}", header_start_pos);
    println!("  Blocks + Index: {} - {} bytes",header_start_pos, footer_pos);
    println!("  Footer: 32 bytes at offset {}", footer_pos);
    println!();

    // For BAM comparison
    println!("Compression Analysis:");
    println!("  CAF: {:.2} MB", file_size as f64 / (1024.0 * 1024.0));
    println!("  Expected if similar to BAM: ~1 MB");
    println!("  Current expansion: ~{}×", file_size as f64 / (1024.0 * 1024.0));

    Ok(())
}
