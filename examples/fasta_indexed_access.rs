//! Demonstrate FAI (FASTA Index) usage for random access to sequences
//!
//! This example shows:
//! - Building a FASTA index from scratch
//! - Loading an existing index
//! - Fetching entire sequences by name
//! - Fetching specific regions (subsequences)
//! - Querying index metadata
//!
//! Run with: cargo run --example fasta_indexed_access

use biometal::io::fasta::index::FaiIndex;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn Error>> {
    // Create a temporary FASTA file for demonstration
    let test_dir = PathBuf::from("examples/data");
    std::fs::create_dir_all(&test_dir)?;

    let fasta_path = test_dir.join("demo.fasta");
    let fai_path = test_dir.join("demo.fasta.fai");

    // Write sample FASTA file
    let mut fasta = File::create(&fasta_path)?;
    writeln!(fasta, ">chr1 Human chromosome 1")?;
    writeln!(fasta, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")?;
    writeln!(fasta, "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA")?;
    writeln!(fasta, ">chr2 Human chromosome 2")?;
    writeln!(fasta, "GGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCC")?;
    writeln!(fasta, "AAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTX")?;
    writeln!(fasta, ">chrM Mitochondrial DNA")?;
    writeln!(fasta, "ATCGATCGATCGATCG")?;
    drop(fasta);

    println!("=== FASTA Indexed Access Demo ===\n");

    // === 1. Build Index ===
    println!("1. Building FASTA index...");
    let index = FaiIndex::build(&fasta_path)?;
    println!("   ✓ Index built with {} sequences\n", index.sequence_names.len());

    // === 2. Write Index to File ===
    println!("2. Writing index to {}", fai_path.display());
    index.write(&fai_path)?;
    println!("   ✓ Index written\n");

    // === 3. Load Index from File ===
    println!("3. Loading index from file...");
    let loaded_index = FaiIndex::from_path(&fai_path)?;
    println!("   ✓ Index loaded\n");

    // === 4. Query Sequence Metadata ===
    println!("4. Querying sequence metadata:");
    for seq_name in &loaded_index.sequence_names {
        if let Some(info) = loaded_index.sequences.get(seq_name) {
            println!("   {} (length={}, offset={}, line_bases={}, line_width={})",
                seq_name, info.length, info.offset, info.line_bases, info.line_width);
        }
    }
    println!();

    // === 5. Fetch Entire Sequences ===
    println!("5. Fetching entire sequences:");

    let chr1 = loaded_index.fetch("chr1", &fasta_path)?;
    println!("   chr1: {} bp", chr1.len());
    println!("   First 40 bp: {}", &chr1[..40]);

    let chrM = loaded_index.fetch("chrM", &fasta_path)?;
    println!("   chrM: {} bp", chrM.len());
    println!("   Full sequence: {}", chrM);
    println!();

    // === 6. Fetch Specific Regions ===
    println!("6. Fetching specific regions:");

    // Fetch chr1:0-20
    let region1 = loaded_index.fetch_region("chr1", 0, 20, &fasta_path)?;
    println!("   chr1:0-20   = {}", region1);

    // Fetch chr1:40-60 (crosses line boundary)
    let region2 = loaded_index.fetch_region("chr1", 40, 60, &fasta_path)?;
    println!("   chr1:40-60  = {}", region2);

    // Fetch chr2:20-40
    let region3 = loaded_index.fetch_region("chr2", 20, 40, &fasta_path)?;
    println!("   chr2:20-40  = {}", region3);
    println!();

    // === 7. Check Sequence Existence ===
    println!("7. Checking sequence existence:");
    println!("   chr1 exists: {}", loaded_index.sequences.contains_key("chr1"));
    println!("   chr99 exists: {}", loaded_index.sequences.contains_key("chr99"));
    println!();

    // === 8. Practical Use Case: Extract Gene Region ===
    println!("8. Practical example - Extract gene region:");
    println!("   Simulating gene at chr1:25-45");
    let gene_sequence = loaded_index.fetch_region("chr1", 25, 45, &fasta_path)?;
    println!("   Gene sequence: {}", gene_sequence);
    println!("   Length: {} bp", gene_sequence.len());

    // Calculate GC content
    let gc_count = gene_sequence.chars()
        .filter(|&c| c == 'G' || c == 'C')
        .count();
    let gc_percent = (gc_count as f64 / gene_sequence.len() as f64) * 100.0;
    println!("   GC content: {:.1}%", gc_percent);
    println!();

    // === 9. Performance Characteristics ===
    println!("9. Performance characteristics:");
    println!("   - Index loading: O(n) where n = number of sequences");
    println!("   - Sequence lookup: O(1) hash map access");
    println!("   - Region extraction: O(1) seek + O(m) read where m = region length");
    println!("   - Memory usage: ~200 bytes per sequence entry");
    println!();

    println!("✅ All operations completed successfully!");
    println!("\n💡 Key takeaways:");
    println!("   • FAI enables O(1) random access to any sequence");
    println!("   • No need to scan entire file for specific regions");
    println!("   • Ideal for: gene extraction, region queries, reference lookups");
    println!("   • Compatible with samtools faidx format");

    // Cleanup
    std::fs::remove_file(&fasta_path).ok();
    std::fs::remove_file(&fai_path).ok();

    Ok(())
}
