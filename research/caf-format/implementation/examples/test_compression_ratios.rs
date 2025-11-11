//! Test actual compression ratios on real data.

use biometal::io::bam::BamReader;
use caf::block::{BlockBuilder, AlignmentRecord};
use caf::conversion::bam_record_to_alignment_record;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Testing Compression Ratios on Real Data\n");

    // Read first 10K records from BAM
    let bam_path = Path::new("/Users/scotthandley/Code/biometal/tests/data/synthetic_100k.bam");
    let mut bam_reader = BamReader::from_path(bam_path)?;

    let mut builder = BlockBuilder::new(0, 10_000);
    let mut count = 0;

    for result in bam_reader.records() {
        let record = result?;
        let caf_record = bam_record_to_alignment_record(&record)?;
        builder.add_record(caf_record)?;

        count += 1;
        if count >= 10_000 {
            break;
        }
    }

    println!("Built block with {} records", builder.len());

    let block = builder.build()?;

    println!("\nBlock Statistics:");
    println!("  Block ID: {}", block.block_id);
    println!("  Num records: {}", block.num_records);
    println!("  Uncompressed size: {:.1} KB", block.uncompressed_size as f64 / 1024.0);
    println!();

    // Analyze each column
    println!("Per-Column Compression:");
    println!("Column              | Uncomp | Comp  | Ratio");
    println!("--------------------|--------|-------|------");

    fn print_col<T>(name: &str, col: &caf::types::CompressedColumn<T>) {
        let ratio = col.compression_ratio();
        println!(
            "{:20}| {:6} | {:5} | {:5.2}×",
            name,
            col.uncompressed_len,
            col.compressed_len,
            ratio
        );
    }

    print_col("ref_ids", &block.columns.ref_ids);
    print_col("positions", &block.columns.positions);
    print_col("mapq", &block.columns.mapq);
    print_col("flags", &block.columns.flags);
    print_col("sequences", &block.columns.sequences);
    print_col("seq_offsets", &block.columns.seq_offsets);
    print_col("qualities", &block.columns.qualities);
    print_col("qual_offsets", &block.columns.qual_offsets);
    print_col("cigar_ops", &block.columns.cigar_ops);
    print_col("cigar_offsets", &block.columns.cigar_offsets);
    print_col("read_names", &block.columns.read_names);
    print_col("read_name_offsets", &block.columns.read_name_offsets);
    print_col("mate_ref_ids", &block.columns.mate_ref_ids);
    print_col("mate_positions", &block.columns.mate_positions);
    print_col("template_lengths", &block.columns.template_lengths);

    // Calculate total
    let total_compressed: u32 = block.columns.ref_ids.compressed_len
        + block.columns.positions.compressed_len
        + block.columns.mapq.compressed_len
        + block.columns.flags.compressed_len
        + block.columns.sequences.compressed_len
        + block.columns.seq_offsets.compressed_len
        + block.columns.qualities.compressed_len
        + block.columns.qual_offsets.compressed_len
        + block.columns.cigar_ops.compressed_len
        + block.columns.cigar_offsets.compressed_len
        + block.columns.read_names.compressed_len
        + block.columns.read_name_offsets.compressed_len
        + block.columns.mate_ref_ids.compressed_len
        + block.columns.mate_positions.compressed_len
        + block.columns.template_lengths.compressed_len;

    println!("--------------------|--------|-------|------");
    println!(
        "{:20}| {:6} | {:5} | {:5.2}×",
        "TOTAL",
        block.uncompressed_size,
        total_compressed,
        block.uncompressed_size as f64 / total_compressed as f64
    );

    println!("\nSerializedblock size: {} bytes", bincode::serialize(&block)?.len());

    Ok(())
}
