//! Simple BAM parsing example for profiling
//!
//! This example parses a BAM file and accesses all fields to generate
//! a realistic flamegraph profile.
//!
//! Run with: cargo flamegraph --example profile_bam

use biometal::io::bam::BamReader;
use std::env;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Allow file path to be specified via command line, default to large file
    let bam_path = env::args()
        .nth(1)
        .unwrap_or_else(|| "tests/data/large/large_1m.bam".to_string());

    eprintln!("Profiling BAM file: {}", bam_path);

    // Run parsing multiple times to get good profiling data
    for iteration in 0..10 {
        let mut bam = BamReader::from_path(&bam_path)?;

        let mut count = 0;
        let mut total_seq_len = 0;
        let mut total_cigar_ops = 0;
        let mut total_tags = 0;

        for record in bam.records() {
            let record = record?;

            // Access all fields to trigger all parsing paths
            count += 1;
            total_seq_len += record.sequence.len();
            total_cigar_ops += record.cigar.len();
            total_tags += record.tags.len();

            // Access some common fields
            let _name = &record.name;
            let _pos = record.position;
            let _mapq = record.mapq;
            let _flags = record.flags;
        }

        if iteration == 0 {
            eprintln!("Parsed {} records", count);
            eprintln!("Total sequence length: {} bases", total_seq_len);
            eprintln!("Total CIGAR ops: {}", total_cigar_ops);
            eprintln!("Total tags: {}", total_tags);
        }
    }

    Ok(())
}
