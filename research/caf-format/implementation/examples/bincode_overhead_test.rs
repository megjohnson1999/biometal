//! Test bincode serialization overhead for CAF structures.

use caf::types::{CompressedColumn, CompressionType};

fn main() {
    println!("Bincode Serialization Overhead Analysis\n");

    // Test 1: Empty CompressedColumn
    let empty: CompressedColumn<u8> = CompressedColumn::new(
        CompressionType::Raw,
        0,
        0,
        Vec::new(),
    );
    let empty_size = bincode::serialize(&empty).unwrap().len();
    println!("Empty CompressedColumn<u8>: {} bytes", empty_size);

    // Test 2: CompressedColumn with 1000 bytes of data
    let data_1k: CompressedColumn<u8> = CompressedColumn::new(
        CompressionType::Zstd,
        1000,
        2000,
        vec![0u8; 1000],
    );
    let data_1k_size = bincode::serialize(&data_1k).unwrap().len();
    let overhead_1k = data_1k_size - 1000;
    println!("CompressedColumn with 1000 bytes:");
    println!("  Total serialized: {} bytes", data_1k_size);
    println!("  Actual data: 1000 bytes");
    println!("  Overhead: {} bytes ({:.1}%)", overhead_1k, (overhead_1k as f64 / 1000.0) * 100.0);

    // Test 3: Estimate overhead for 15 columns per block
    println!("\nPer-block overhead (15 columns):");
    let columns_per_block = 15;
    let estimated_overhead_per_block = empty_size * columns_per_block;
    println!("  Metadata overhead: ~{} bytes", estimated_overhead_per_block);

    // Test 4: Estimate for 10 blocks (100K records)
    let num_blocks = 10;
    let total_overhead = estimated_overhead_per_block * num_blocks;
    println!("\nTotal overhead for {} blocks:", num_blocks);
    println!("  Metadata: ~{} bytes ({:.1} KB)", total_overhead, total_overhead as f64 / 1024.0);

    // Test 5: Vec<u8> length prefix
    let small_vec = vec![0u8; 10];
    let small_vec_size = bincode::serialize(&small_vec).unwrap().len();
    println!("\nVec<u8> overhead:");
    println!("  Vec with 10 bytes: {} bytes total ({} overhead)", small_vec_size, small_vec_size - 10);

    // Additional analysis
    println!("\nKey Findings:");
    println!("- Each CompressedColumn has {} bytes of metadata overhead", empty_size);
    println!("- 15 columns × {} blocks = {} CompressedColumn structures", num_blocks, columns_per_block * num_blocks);
    println!("- Total metadata overhead: ~{:.1} KB", total_overhead as f64 / 1024.0);
    println!("\nThis overhead accounts for the file size bloat.");
}
