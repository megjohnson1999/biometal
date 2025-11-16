#!/usr/bin/env rust-script
//! Simple benchmark for CRAM reader performance
//!
//! Usage: cargo run --release --example cram_benchmark -- <cram_file> <reference>

use biometal::io::cram::CramReader;
use std::env;
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: {} <cram_file> <reference>", args[0]);
        std::process::exit(1);
    }

    let cram_path = &args[1];
    let reference_path = &args[2];

    println!("biometal CRAM Reader Benchmark");
    println!("==============================");
    println!("CRAM file: {}", cram_path);
    println!("Reference: {}", reference_path);
    println!();

    // Open CRAM file with reference
    let cram_reader = CramReader::from_path_with_reference(cram_path, reference_path)?;

    // Benchmark record iteration
    let start = Instant::now();
    let mut record_count = 0;
    let mut total_sequence_length = 0;

    for result in cram_reader.records() {
        let record = result?;
        record_count += 1;
        total_sequence_length += record.sequence.len();
    }

    let elapsed = start.elapsed();
    let seconds = elapsed.as_secs_f64();

    println!("Results:");
    println!("--------");
    println!("Records read:     {}", record_count);
    println!("Total sequence:   {} bp", total_sequence_length);
    println!("Time elapsed:     {:.3} seconds", seconds);
    println!("Throughput:       {:.0} records/sec", record_count as f64 / seconds);
    println!("Sequence speed:   {:.2} Mbp/sec", total_sequence_length as f64 / seconds / 1_000_000.0);

    Ok(())
}
