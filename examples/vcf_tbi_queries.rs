//! Demonstrate TBI (Tabix Index) usage for VCF region queries
//!
//! This example shows:
//! - Loading a TBI index for VCF files
//! - Querying specific genomic regions
//! - Getting file chunks to read for efficient I/O
//! - Querying index metadata
//! - Practical workflow for targeted variant analysis
//!
//! Run with: cargo run --example vcf_tbi_queries

use biometal::formats::index::TbiIndex;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

/// Create a minimal TBI index for demonstration
fn create_demo_tbi() -> Vec<u8> {
    use std::io::Cursor;

    let mut data = Vec::new();
    let mut cursor = Cursor::new(&mut data);

    // Magic string "TBI\1"
    cursor.write_all(b"TBI\x01").unwrap();

    // n_ref = 2 (chr1, chr2)
    cursor.write_all(&2i32.to_le_bytes()).unwrap();

    // format = 2 (VCF)
    cursor.write_all(&2i32.to_le_bytes()).unwrap();

    // col_seq = 0 (CHROM column)
    cursor.write_all(&0i32.to_le_bytes()).unwrap();

    // col_beg = 1 (POS column)
    cursor.write_all(&1i32.to_le_bytes()).unwrap();

    // col_end = 0 (same as beg for VCF)
    cursor.write_all(&0i32.to_le_bytes()).unwrap();

    // meta = '#'
    cursor.write_all(&(b'#' as i32).to_le_bytes()).unwrap();

    // skip = 0
    cursor.write_all(&0i32.to_le_bytes()).unwrap();

    // l_nm = 10 ("chr1\0chr2\0")
    cursor.write_all(&10i32.to_le_bytes()).unwrap();
    cursor.write_all(b"chr1\0chr2\0").unwrap();

    // Index for chr1: 2 bins with chunks
    cursor.write_all(&2i32.to_le_bytes()).unwrap(); // n_bin

    // Bin 0 (covers entire chromosome)
    cursor.write_all(&0u32.to_le_bytes()).unwrap(); // bin_id
    cursor.write_all(&1i32.to_le_bytes()).unwrap(); // n_chunk
    cursor.write_all(&0x100000u64.to_le_bytes()).unwrap(); // chunk start
    cursor.write_all(&0x200000u64.to_le_bytes()).unwrap(); // chunk end

    // Bin 4681 (smallest bin, specific region)
    cursor.write_all(&4681u32.to_le_bytes()).unwrap(); // bin_id
    cursor.write_all(&1i32.to_le_bytes()).unwrap(); // n_chunk
    cursor.write_all(&0x150000u64.to_le_bytes()).unwrap(); // chunk start
    cursor.write_all(&0x180000u64.to_le_bytes()).unwrap(); // chunk end

    // Linear index for chr1: 16kb windows
    cursor.write_all(&3i32.to_le_bytes()).unwrap(); // n_intv
    cursor.write_all(&0x100000u64.to_le_bytes()).unwrap(); // interval 0
    cursor.write_all(&0x120000u64.to_le_bytes()).unwrap(); // interval 1
    cursor.write_all(&0x150000u64.to_le_bytes()).unwrap(); // interval 2

    // Index for chr2: 1 bin
    cursor.write_all(&1i32.to_le_bytes()).unwrap(); // n_bin

    // Bin 0
    cursor.write_all(&0u32.to_le_bytes()).unwrap(); // bin_id
    cursor.write_all(&1i32.to_le_bytes()).unwrap(); // n_chunk
    cursor.write_all(&0x300000u64.to_le_bytes()).unwrap(); // chunk start
    cursor.write_all(&0x400000u64.to_le_bytes()).unwrap(); // chunk end

    // Linear index for chr2
    cursor.write_all(&2i32.to_le_bytes()).unwrap(); // n_intv
    cursor.write_all(&0x300000u64.to_le_bytes()).unwrap(); // interval 0
    cursor.write_all(&0x350000u64.to_le_bytes()).unwrap(); // interval 1

    data
}

fn main() -> Result<(), Box<dyn Error>> {
    let test_dir = PathBuf::from("examples/data");
    std::fs::create_dir_all(&test_dir)?;

    let tbi_path = test_dir.join("demo.vcf.gz.tbi");

    // Create demo TBI file
    let tbi_data = create_demo_tbi();
    let mut tbi_file = File::create(&tbi_path)?;
    tbi_file.write_all(&tbi_data)?;
    drop(tbi_file);

    println!("=== VCF TBI Region Queries Demo ===\n");

    // === 1. Load TBI Index ===
    println!("1. Loading TBI index...");
    let index = TbiIndex::from_path(&tbi_path)?;
    println!("   ✓ Index loaded ({} references)", index.references().len());
    println!();

    // === 2. Query Index Metadata ===
    println!("2. Index metadata:");
    println!("   Format: {:?}", index.format());
    println!("   CHROM column: {}", index.col_seq());
    println!("   POS column: {}", index.col_beg());
    println!("   Comment character: '{}'", index.meta_char());
    println!("   Skip lines: {}", index.skip_lines());
    println!();

    // === 3. List Available References ===
    println!("3. Available references:");
    for (i, ref_info) in index.references().iter().enumerate() {
        println!("   [{}] {} ({} bins, {} linear intervals)",
            i, ref_info.name, ref_info.bins.len(), ref_info.intervals.len());
    }
    println!();

    // === 4. Check Reference Existence ===
    println!("4. Checking reference existence:");
    println!("   chr1 exists: {}", index.get_reference("chr1").is_some());
    println!("   chr2 exists: {}", index.get_reference("chr2").is_some());
    println!("   chr99 exists: {}", index.get_reference("chr99").is_some());
    println!();

    // === 5. Query Specific Regions ===
    println!("5. Querying genomic regions:");

    // Query chr1:0-100000
    println!("\n   Query: chr1:0-100000");
    let chunks = index.query("chr1", 0, 100000)?;
    println!("   Found {} chunks to read:", chunks.len());
    for (i, chunk) in chunks.iter().enumerate() {
        println!("     Chunk {}: offset 0x{:016x} - 0x{:016x}",
            i, chunk.start.as_raw(), chunk.end.as_raw());
        println!("               (compressed offset: 0x{:x}, uncompressed: 0x{:x})",
            chunk.start.compressed_offset(), chunk.start.uncompressed_offset());
    }

    // Query chr1:500000-600000
    println!("\n   Query: chr1:500000-600000");
    let chunks2 = index.query("chr1", 500000, 600000)?;
    println!("   Found {} chunks to read:", chunks2.len());
    for (i, chunk) in chunks2.iter().enumerate() {
        println!("     Chunk {}: offset 0x{:016x} - 0x{:016x}",
            i, chunk.start.as_raw(), chunk.end.as_raw());
    }

    // Query chr2:0-50000
    println!("\n   Query: chr2:0-50000");
    let chunks3 = index.query("chr2", 0, 50000)?;
    println!("   Found {} chunks to read:", chunks3.len());
    for (i, chunk) in chunks3.iter().enumerate() {
        println!("     Chunk {}: offset 0x{:016x} - 0x{:016x}",
            i, chunk.start.as_raw(), chunk.end.as_raw());
    }
    println!();

    // === 6. Error Handling ===
    println!("6. Error handling:");

    // Non-existent reference
    match index.query("chr99", 0, 100000) {
        Ok(_) => println!("   ❌ Should have failed for chr99"),
        Err(e) => println!("   ✓ Correctly rejected chr99: {}", e),
    }

    // Invalid range
    match index.query("chr1", 100000, 50000) {
        Ok(_) => println!("   ❌ Should have failed for invalid range"),
        Err(e) => println!("   ✓ Correctly rejected invalid range: {}", e),
    }
    println!();

    // === 7. Practical Workflow ===
    println!("7. Practical workflow - Targeted variant analysis:");
    println!("\n   Scenario: Find variants in BRCA2 gene region (chr13:32,889,611-32,973,805)");
    println!("   (Using chr1:32889000-32974000 as demo)\n");

    let gene_region = "chr1";
    let gene_start = 32889000u32;
    let gene_end = 32974000u32;

    println!("   Step 1: Query TBI index for region chunks");
    let gene_chunks = index.query(gene_region, gene_start, gene_end)?;
    println!("   → Found {} chunks covering this region", gene_chunks.len());

    println!("\n   Step 2: Seek to first chunk offset");
    if let Some(first_chunk) = gene_chunks.first() {
        println!("   → Seek to virtual offset: 0x{:016x}", first_chunk.start.as_raw());
        println!("   → BGZF block: {} @ offset {}",
            first_chunk.start.compressed_offset(),
            first_chunk.start.uncompressed_offset());
    }

    println!("\n   Step 3: Stream VCF records from chunk(s)");
    println!("   → Use biometal::formats::vcf::VcfStream with BGZF reader");
    println!("   → Filter records where POS >= {} AND POS < {}", gene_start, gene_end);

    println!("\n   Step 4: Process variants");
    println!("   → Parse genotypes, calculate allele frequencies");
    println!("   → Filter by quality, depth, etc.");
    println!();

    // === 8. Performance Characteristics ===
    println!("8. Performance characteristics:");
    println!("   - Index loading: O(n) where n = file size (~KB)");
    println!("   - Region query: O(log n) binary search in bins");
    println!("   - Chunk count: Typically 1-10 per query");
    println!("   - Memory usage: ~1-10 MB for typical genome index");
    println!("   - Speedup: 10-1000× vs full file scan (depends on region size)");
    println!();

    // === 9. Virtual Offset Explanation ===
    println!("9. Understanding virtual offsets:");
    println!("   BGZF virtual offset = (compressed_offset << 16) | uncompressed_offset");
    println!("   Example: 0x0000000000100ABC");
    println!("     - Compressed offset: 0x100 (256 bytes into file)");
    println!("     - Uncompressed offset: 0xABC (within decompressed block)");
    println!("   This allows seeking to exact positions in compressed files!");
    println!();

    println!("✅ All operations completed successfully!");
    println!("\n💡 Key takeaways:");
    println!("   • TBI enables random access to sorted, BGZF-compressed files");
    println!("   • Hierarchical binning provides O(log n) region queries");
    println!("   • Virtual offsets enable precise seeking in compressed data");
    println!("   • Ideal for: targeted variant calling, region-specific analysis");
    println!("   • Compatible with samtools tabix format");

    // Cleanup
    std::fs::remove_file(&tbi_path).ok();

    Ok(())
}
