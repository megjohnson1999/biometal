//! Analyze CAF file compression efficiency.
//!
//! Usage:
//!   cargo run --example analyze_caf_compression -- input.caf

use caf::io::CafFileReader;
use std::env;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();

    if args.len() != 2 {
        eprintln!("Usage: {} <input.caf>", args[0]);
        std::process::exit(1);
    }

    let caf_path = Path::new(&args[1]);

    println!("Analyzing CAF file: {:?}\n", caf_path);

    // Open CAF file
    let reader = CafFileReader::open(caf_path)?;
    let header = reader.header();

    println!("Header Information:");
    println!("  Version: {}.{}", header.version_major(), header.version_minor());
    println!("  Block size: {}", header.block_size);
    println!("  Num blocks: {}", header.num_blocks);
    println!("  Num refs: {}", header.num_refs);
    println!("  SAM header size: {} bytes", header.sam_header.len());
    println!();

    // Analyze blocks
    let mut total_uncompressed = 0u64;
    let mut total_compressed = 0u64;
    let mut total_records = 0u64;
    let num_blocks = header.num_blocks;

    println!("Per-Block Compression Analysis:");
    println!("Block | Records | Uncomp (KB) | Comp (KB) | Ratio");
    println!("------|---------|-------------|-----------|------");

    for block_id in 0..num_blocks {
        let block_reader = reader.read_block(block_id)?;
        let block = block_reader.block();

        let block_compressed = calculate_compressed_size(&block.columns);
        let block_uncompressed = block.uncompressed_size as u64;

        total_uncompressed += block_uncompressed;
        total_compressed += block_compressed;
        total_records += block.num_records as u64;

        let ratio = if block_compressed > 0 {
            block_uncompressed as f64 / block_compressed as f64
        } else {
            0.0
        };

        println!(
            "{:5} | {:7} | {:11.1} | {:9.1} | {:5.2}×",
            block_id,
            block.num_records,
            block_uncompressed as f64 / 1024.0,
            block_compressed as f64 / 1024.0,
            ratio
        );
    }

    println!();
    println!("Overall Statistics:");
    println!("  Total records: {}", total_records);
    println!("  Total uncompressed: {:.1} KB", total_uncompressed as f64 / 1024.0);
    println!("  Total compressed: {:.1} KB", total_compressed as f64 / 1024.0);
    println!(
        "  Overall ratio: {:.2}×",
        total_uncompressed as f64 / total_compressed as f64
    );
    println!();

    // File size analysis
    let file_size = std::fs::metadata(caf_path)?.len();
    let serialization_overhead = file_size as i64 - total_compressed as i64;
    let overhead_pct = (serialization_overhead as f64 / file_size as f64) * 100.0;

    println!("File Size Breakdown:");
    println!("  Actual file size: {:.1} KB ({} bytes)", file_size as f64 / 1024.0, file_size);
    println!("  Compressed data: {:.1} KB", total_compressed as f64 / 1024.0);
    println!("  Serialization overhead: {:.1} KB ({:.1}%)",
        serialization_overhead as f64 / 1024.0,
        overhead_pct
    );

    Ok(())
}

fn calculate_compressed_size(columns: &caf::types::ColumnData) -> u64 {
    columns.ref_ids.data.len() as u64
        + columns.positions.data.len() as u64
        + columns.mapq.data.len() as u64
        + columns.flags.data.len() as u64
        + columns.sequences.data.len() as u64
        + columns.seq_offsets.data.len() as u64
        + columns.qualities.data.len() as u64
        + columns.qual_offsets.data.len() as u64
        + columns.cigar_ops.data.len() as u64
        + columns.cigar_offsets.data.len() as u64
        + columns.read_names.data.len() as u64
        + columns.read_name_offsets.data.len() as u64
        + columns.mate_ref_ids.data.len() as u64
        + columns.mate_positions.data.len() as u64
        + columns.template_lengths.data.len() as u64
}
