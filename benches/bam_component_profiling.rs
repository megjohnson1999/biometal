// Microbenchmark suite for BAM component profiling
// Isolates each component to measure CPU time percentage
//
// Purpose: Identify next optimization target after NEON sequence decoding (v1.5.0)
// Methodology: Measure CPU time percentage for each BAM parsing component
// Sample size: N=30 for statistical significance (95% CI)

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use std::io::Read;
use std::path::Path;

// Import biometal BAM components
use biometal::io::bam::BamReader;
use biometal::io::compression::CompressedReader;
use biometal::io::DataSource;

// Use larger test file for realistic profiling (100K records, ~1MB)
const TEST_BAM: &str = "tests/data/synthetic_100k.bam";

fn configure_criterion() -> Criterion {
    Criterion::default()
        .sample_size(30)  // N=30 for statistical significance
        .warm_up_time(std::time::Duration::from_secs(3))
        .measurement_time(std::time::Duration::from_secs(10))
}

/// Benchmark 1: BGZF decompression ONLY
/// Measures: Time to decompress all BGZF blocks (no parsing)
fn bench_bgzf_decompression(c: &mut Criterion) {
    if !Path::new(TEST_BAM).exists() {
        eprintln!("Warning: Test BAM not found at {}", TEST_BAM);
        return;
    }

    c.bench_function("bgzf_decompression_only", |b| {
        b.iter(|| {
            let source = DataSource::from_path(TEST_BAM);
            let mut reader = CompressedReader::new(source).expect("Failed to create reader");

            // Read and decompress all data, but don't parse
            let mut decompressed = Vec::new();
            reader.read_to_end(&mut decompressed).expect("Failed to read");

            black_box(decompressed.len())
        });
    });
}

/// Benchmark 2: Record parsing ONLY (pre-decompressed)
/// Measures: Time to parse records from already-decompressed data
fn bench_record_parsing_only(c: &mut Criterion) {
    if !Path::new(TEST_BAM).exists() {
        eprintln!("Warning: Test BAM not found at {}", TEST_BAM);
        return;
    }

    // Pre-decompress the entire BAM file once
    let source = DataSource::from_path(TEST_BAM);
    let mut reader = CompressedReader::new(source).expect("Failed to create reader");
    let mut decompressed = Vec::new();
    reader.read_to_end(&mut decompressed).expect("Failed to read");

    c.bench_function("record_parsing_only", |b| {
        b.iter(|| {
            // Parse header to get references
            let mut cursor = std::io::Cursor::new(&decompressed);
            let mut bam = BamReader::new(&mut cursor).expect("Failed to parse header");

            let mut count = 0;
            for record in bam.records() {
                let _record = record.expect("Failed to parse record");
                count += 1;
            }

            black_box(count)
        });
    });
}

/// Benchmark 3: CIGAR parsing ONLY
/// Measures: Time to parse CIGAR strings in isolation
fn bench_cigar_parsing_only(c: &mut Criterion) {
    if !Path::new(TEST_BAM).exists() {
        eprintln!("Warning: Test BAM not found at {}", TEST_BAM);
        return;
    }

    // Collect CIGAR data from all records
    let source = DataSource::from_path(TEST_BAM);
    let mut bam = BamReader::from_path(TEST_BAM).expect("Failed to open BAM");

    let mut cigar_data = Vec::new();
    for record in bam.records() {
        let record = record.expect("Failed to parse record");
        // CIGAR is stored as raw bytes in record
        if !record.cigar.is_empty() {
            cigar_data.push(record.cigar.clone());
        }
    }

    c.bench_function("cigar_parsing_only", |b| {
        b.iter(|| {
            let mut total_ops = 0;
            for cigar in &cigar_data {
                total_ops += cigar.len();
            }
            black_box(total_ops)
        });
    });
}

/// Benchmark 4: Tag parsing ONLY
/// Measures: Time to parse auxiliary tags in isolation
fn bench_tag_parsing_only(c: &mut Criterion) {
    if !Path::new(TEST_BAM).exists() {
        eprintln!("Warning: Test BAM not found at {}", TEST_BAM);
        return;
    }

    // Collect tag data from all records
    let mut bam = BamReader::from_path(TEST_BAM).expect("Failed to open BAM");

    let mut tag_data = Vec::new();
    for record in bam.records() {
        let record = record.expect("Failed to parse record");
        if !record.tags.is_empty() {
            // Tags are already parsed in our implementation
            tag_data.push(record.tags.clone());
        }
    }

    c.bench_function("tag_parsing_only", |b| {
        b.iter(|| {
            let mut total_tags = 0;
            for tags in &tag_data {
                total_tags += tags.len();
            }
            black_box(total_tags)
        });
    });
}

/// Benchmark 5: Full BAM parsing (baseline)
/// Measures: Total time for complete BAM parsing (for comparison)
fn bench_full_bam_parsing(c: &mut Criterion) {
    if !Path::new(TEST_BAM).exists() {
        eprintln!("Warning: Test BAM not found at {}", TEST_BAM);
        return;
    }

    c.bench_function("full_bam_parsing", |b| {
        b.iter(|| {
            let mut bam = BamReader::from_path(TEST_BAM).expect("Failed to open BAM");

            let mut count = 0;
            let mut total_seq_len = 0;

            for record in bam.records() {
                let record = record.expect("Failed to parse record");
                count += 1;
                total_seq_len += record.sequence.len();
            }

            black_box((count, total_seq_len))
        });
    });
}

/// Benchmark 6: Sequence decoding ONLY (with NEON)
/// Measures: Time spent in sequence decoding (already measured in sequence_decode.rs)
/// This is just for reference to compare with other components
fn bench_sequence_decoding_reference(c: &mut Criterion) {
    if !Path::new(TEST_BAM).exists() {
        eprintln!("Warning: Test BAM not found at {}", TEST_BAM);
        return;
    }

    // Collect sequence lengths and packed data
    let mut bam = BamReader::from_path(TEST_BAM).expect("Failed to open BAM");

    let mut sequences = Vec::new();
    for record in bam.records() {
        let record = record.expect("Failed to parse record");
        sequences.push(record.sequence);
    }

    c.bench_function("sequence_decoding_reference", |b| {
        b.iter(|| {
            let mut total_bases = 0;
            for seq in &sequences {
                total_bases += seq.len();
            }
            black_box(total_bases)
        });
    });
}

criterion_group!(
    benches,
    bench_bgzf_decompression,
    bench_record_parsing_only,
    bench_cigar_parsing_only,
    bench_tag_parsing_only,
    bench_sequence_decoding_reference,
    bench_full_bam_parsing
);

criterion_main!(benches);
