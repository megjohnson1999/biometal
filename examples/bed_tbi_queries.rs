//! Demonstrate TBI usage for BED file peak lookups
//!
//! This example shows how to use TBI indices for efficient region queries
//! in BED files (ChIP-seq peaks, ATAC-seq peaks, genomic intervals)
//!
//! Run with: cargo run --example bed_tbi_queries

use biometal::formats::index::TbiIndex;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

fn create_bed_tbi() -> Vec<u8> {
    use std::io::Cursor;

    let mut data = Vec::new();
    let mut cursor = Cursor::new(&mut data);

    // Magic
    cursor.write_all(b"TBI\x01").unwrap();

    // n_ref = 3 (chr1, chr2, chrX)
    cursor.write_all(&3i32.to_le_bytes()).unwrap();

    // format = 0 (Generic BED)
    cursor.write_all(&0i32.to_le_bytes()).unwrap();

    // col_seq = 0, col_beg = 1, col_end = 2
    cursor.write_all(&0i32.to_le_bytes()).unwrap();
    cursor.write_all(&1i32.to_le_bytes()).unwrap();
    cursor.write_all(&2i32.to_le_bytes()).unwrap();

    // meta = '#'
    cursor.write_all(&(b'#' as i32).to_le_bytes()).unwrap();

    // skip = 1 (BED header line)
    cursor.write_all(&1i32.to_le_bytes()).unwrap();

    // l_nm = 15 ("chr1\0chr2\0chrX\0")
    cursor.write_all(&15i32.to_le_bytes()).unwrap();
    cursor.write_all(b"chr1\0chr2\0chrX\0").unwrap();

    // Index for chr1 (promoter regions)
    cursor.write_all(&2i32.to_le_bytes()).unwrap(); // n_bin
    cursor.write_all(&0u32.to_le_bytes()).unwrap(); // bin 0
    cursor.write_all(&1i32.to_le_bytes()).unwrap();
    cursor.write_all(&0x50000u64.to_le_bytes()).unwrap();
    cursor.write_all(&0xA0000u64.to_le_bytes()).unwrap();
    cursor.write_all(&4681u32.to_le_bytes()).unwrap(); // bin 4681
    cursor.write_all(&1i32.to_le_bytes()).unwrap();
    cursor.write_all(&0x60000u64.to_le_bytes()).unwrap();
    cursor.write_all(&0x80000u64.to_le_bytes()).unwrap();
    cursor.write_all(&2i32.to_le_bytes()).unwrap(); // n_intv
    cursor.write_all(&0x50000u64.to_le_bytes()).unwrap();
    cursor.write_all(&0x60000u64.to_le_bytes()).unwrap();

    // Index for chr2 (enhancer regions)
    cursor.write_all(&1i32.to_le_bytes()).unwrap();
    cursor.write_all(&0u32.to_le_bytes()).unwrap();
    cursor.write_all(&1i32.to_le_bytes()).unwrap();
    cursor.write_all(&0xB0000u64.to_le_bytes()).unwrap();
    cursor.write_all(&0xF0000u64.to_le_bytes()).unwrap();
    cursor.write_all(&1i32.to_le_bytes()).unwrap();
    cursor.write_all(&0xB0000u64.to_le_bytes()).unwrap();

    // Index for chrX (X-chromosome specific peaks)
    cursor.write_all(&1i32.to_le_bytes()).unwrap();
    cursor.write_all(&0u32.to_le_bytes()).unwrap();
    cursor.write_all(&1i32.to_le_bytes()).unwrap();
    cursor.write_all(&0x100000u64.to_le_bytes()).unwrap();
    cursor.write_all(&0x150000u64.to_le_bytes()).unwrap();
    cursor.write_all(&1i32.to_le_bytes()).unwrap();
    cursor.write_all(&0x100000u64.to_le_bytes()).unwrap();

    data
}

fn main() -> Result<(), Box<dyn Error>> {
    let test_dir = PathBuf::from("examples/data");
    std::fs::create_dir_all(&test_dir)?;

    let tbi_path = test_dir.join("peaks.bed.gz.tbi");

    // Create TBI index
    let tbi_data = create_bed_tbi();
    let mut tbi_file = File::create(&tbi_path)?;
    tbi_file.write_all(&tbi_data)?;
    drop(tbi_file);

    println!("=== BED Peak Lookup with TBI Demo ===\n");

    // === 1. Load Index ===
    println!("1. Loading BED TBI index...");
    let index = TbiIndex::from_path(&tbi_path)?;
    println!("   ✓ Loaded index for {} chromosomes\n", index.references().len());

    // === 2. Show Indexed Chromosomes ===
    println!("2. Indexed chromosomes:");
    for ref_info in index.references() {
        println!("   {} ({} bins, {} intervals)",
            ref_info.name, ref_info.bins.len(), ref_info.intervals.len());
    }
    println!();

    // === 3. Use Case 1: ChIP-seq Peak Overlap ===
    println!("3. Use Case: ChIP-seq peak overlap analysis");
    println!("\n   Scenario: Check if gene promoter overlaps with H3K27ac peaks");
    println!("   Gene: BRCA1 promoter at chr1:43,044,000-43,046,000\n");

    let promoter_chr = "chr1";
    let promoter_start = 43_044_000u32;
    let promoter_end = 43_046_000u32;

    println!("   Step 1: Query TBI for overlapping peaks");
    let chunks = index.query(promoter_chr, promoter_start, promoter_end)?;
    println!("   → Found {} chunks to scan", chunks.len());

    println!("\n   Step 2: Read chunks and filter peaks");
    println!("   → Read from offset: 0x{:016x}", chunks[0].start.as_raw());
    println!("   → Stream BED records until offset: 0x{:016x}", chunks[0].end.as_raw());

    println!("\n   Step 3: Check overlap");
    println!("   → For each peak: max(peak_start, {}) < min(peak_end, {})",
        promoter_start, promoter_end);
    println!();

    // === 4. Use Case 2: ATAC-seq Peak Enrichment ===
    println!("4. Use Case: ATAC-seq peak enrichment in specific region");
    println!("\n   Scenario: Count accessible chromatin peaks in 1 Mb window");
    println!("   Region: chr2:100,000,000-101,000,000\n");

    let window_chr = "chr2";
    let window_start = 100_000_000u32;
    let window_end = 101_000_000u32;

    println!("   Query TBI index...");
    let chunks = index.query(window_chr, window_start, window_end)?;
    println!("   → Found {} chunks", chunks.len());

    println!("\n   Workflow:");
    println!("   1. Stream BED records from chunks");
    println!("   2. Count peaks fully contained in [{}, {})", window_start, window_end);
    println!("   3. Calculate peak density = count / 1,000,000 bp");
    println!();

    // === 5. Use Case 3: Sex Chromosome Analysis ===
    println!("5. Use Case: X chromosome-specific peak calling");
    println!("\n   Scenario: Extract all peaks on chrX for dosage compensation study\n");

    println!("   Query entire chrX (0-156,000,000)");
    let chunks = index.query("chrX", 0, 156_000_000)?;
    println!("   → Found {} chunks covering chrX", chunks.len());

    println!("\n   Analysis pipeline:");
    println!("   1. Stream all chrX peaks from index");
    println!("   2. Calculate peak width distribution");
    println!("   3. Compare to autosomal chromosomes");
    println!("   4. Identify dosage compensation patterns");
    println!();

    // === 6. Multi-Region Query Optimization ===
    println!("6. Optimization: Batch region queries");
    println!("\n   Scenario: Check 1000 promoters for H3K4me3 marks\n");

    let test_promoters = vec![
        ("chr1", 1_000_000, 1_002_000),
        ("chr1", 5_000_000, 5_002_000),
        ("chr2", 10_000_000, 10_002_000),
    ];

    println!("   Querying {} promoters...", test_promoters.len());
    let mut total_chunks = 0;
    for (chr, start, end) in test_promoters {
        match index.query(chr, start, end) {
            Ok(chunks) => {
                total_chunks += chunks.len();
                println!("   {}:{}-{} → {} chunks", chr, start, end, chunks.len());
            }
            Err(_) => println!("   {}:{}-{} → no data", chr, start, end),
        }
    }
    println!("\n   Total chunks to read: {} (vs scanning entire file!)", total_chunks);
    println!();

    // === 7. Performance Comparison ===
    println!("7. Performance comparison:");
    println!("\n   Without TBI (sequential scan):");
    println!("   - Must decompress entire file");
    println!("   - Read all chromosomes even if querying chr22");
    println!("   - Time: O(n) where n = file size");
    println!("   - Example: 500 MB BED file = 10-30 seconds\n");

    println!("   With TBI (indexed access):");
    println!("   - Jump directly to relevant chunks");
    println!("   - Decompress only overlapping blocks");
    println!("   - Time: O(log n + m) where m = result size");
    println!("   - Example: Single gene = 0.01-0.1 seconds");
    println!("   - Speedup: 100-1000× for small regions!");
    println!();

    // === 8. Best Practices ===
    println!("8. Best practices for BED + TBI:");
    println!("   ✓ Sort BED file by chromosome and position");
    println!("   ✓ BGZF-compress with bgzip (not gzip)");
    println!("   ✓ Build TBI with: tabix -p bed file.bed.gz");
    println!("   ✓ Use TBI for queries, streaming for full scans");
    println!("   ✓ Batch small queries to amortize seek costs");
    println!();

    println!("✅ All operations completed successfully!");
    println!("\n💡 Key takeaways:");
    println!("   • TBI makes BED files queryable like databases");
    println!("   • Ideal for: peak overlap, region enrichment, multi-region queries");
    println!("   • 100-1000× faster than sequential scanning");
    println!("   • Essential for large-scale epigenomics analysis");
    println!("   • Compatible with samtools tabix");

    // Cleanup
    std::fs::remove_file(&tbi_path).ok();

    Ok(())
}
