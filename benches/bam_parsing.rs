//! BAM parsing benchmarks (Phase 3)
//!
//! Validates Phase 2 parallel BGZF integration and measures overall performance.
//!
//! # Benchmarks
//!
//! - `parse_bam_100k`: Parse 100K records from BGZF-compressed BAM
//! - `parse_bam_header`: Header parsing overhead
//! - `parse_bam_records_only`: Record parsing without I/O overhead
//!
//! # Expected Results (Phase 2)
//!
//! - Overall speedup: ~4-5× (6.5× BGZF decompression, but I/O is only part of workload)
//! - Memory footprint: Constant ~5 MB regardless of file size (Rule 5)
//! - Throughput: Target >100K records/sec on M1/M2 Mac
//!
//! # Evidence Base
//!
//! - Rule 3 (Parallel BGZF): Entry 029, 6.5× validated speedup
//! - Phase 0 profiling: BGZF 66-80% CPU time → expect 4-5× overall gain

use biometal::io::bam::BamReader;
use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use std::path::Path;

/// Benchmark parsing a complete BAM file (100K records)
///
/// This tests end-to-end performance including:
/// - File I/O
/// - BGZF decompression (parallel, Phase 2)
/// - Header parsing
/// - Record parsing
/// - Iterator overhead
fn bench_parse_bam_100k(c: &mut Criterion) {
    let bam_path = "experiments/native-bam-implementation/test-data/synthetic_100000.bam";

    // Check if test file exists
    if !Path::new(bam_path).exists() {
        eprintln!("Warning: Test file not found: {}", bam_path);
        eprintln!("Skipping benchmark. Run from repository root with test data available.");
        return;
    }

    let mut group = c.benchmark_group("bam_parsing");

    // Get file size for throughput measurement
    let file_size = std::fs::metadata(bam_path)
        .map(|m| m.len())
        .unwrap_or(0);

    group.throughput(Throughput::Bytes(file_size));

    group.bench_function("parse_100k_records", |b| {
        b.iter(|| {
            let mut bam = BamReader::from_path(black_box(bam_path))
                .expect("Failed to open BAM file");

            let mut count = 0;
            for result in bam.records() {
                let record = result.expect("Failed to parse record");
                black_box(&record);
                count += 1;
            }

            assert_eq!(count, 100_000, "Expected 100K records");
            count
        });
    });

    group.finish();
}

/// Benchmark just header parsing
///
/// Isolates header parsing overhead from record parsing.
/// Expected to be negligible (<1% of total time).
fn bench_parse_header(c: &mut Criterion) {
    let bam_path = "experiments/native-bam-implementation/test-data/synthetic_100000.bam";

    if !Path::new(bam_path).exists() {
        return;
    }

    c.bench_function("parse_header", |b| {
        b.iter(|| {
            let bam = BamReader::from_path(black_box(bam_path))
                .expect("Failed to open BAM file");

            // Just access header, don't read records
            black_box(bam.header().reference_count());
        });
    });
}

/// Benchmark record iteration rate
///
/// Measures records/second throughput to validate target >100K/sec.
fn bench_record_throughput(c: &mut Criterion) {
    let bam_path = "experiments/native-bam-implementation/test-data/synthetic_100000.bam";

    if !Path::new(bam_path).exists() {
        return;
    }

    let mut group = c.benchmark_group("record_throughput");

    // Measure in terms of records processed
    group.throughput(Throughput::Elements(100_000));

    group.bench_function("100k_records", |b| {
        b.iter(|| {
            let mut bam = BamReader::from_path(black_box(bam_path))
                .expect("Failed to open BAM file");

            let mut count = 0;
            let mut total_bases = 0u64;

            for result in bam.records() {
                let record = result.expect("Failed to parse record");
                count += 1;
                total_bases += record.sequence.len() as u64;
                black_box(&record);
            }

            (count, total_bases)
        });
    });

    group.finish();
}

/// Benchmark different access patterns
///
/// Tests performance of common operations:
/// - Scan all (just count)
/// - Access name only
/// - Access position only
/// - Full record access
fn bench_access_patterns(c: &mut Criterion) {
    let bam_path = "experiments/native-bam-implementation/test-data/synthetic_100000.bam";

    if !Path::new(bam_path).exists() {
        return;
    }

    let mut group = c.benchmark_group("access_patterns");

    // Just count records (minimal processing)
    group.bench_function("count_only", |b| {
        b.iter(|| {
            let mut bam = BamReader::from_path(black_box(bam_path))
                .expect("Failed to open BAM file");

            let mut count = 0;
            for result in bam.records() {
                let _record = result.expect("Failed to parse record");
                count += 1;
            }
            count
        });
    });

    // Access read names (common filter operation)
    group.bench_function("read_names", |b| {
        b.iter(|| {
            let mut bam = BamReader::from_path(black_box(bam_path))
                .expect("Failed to open BAM file");

            let mut count = 0;
            for result in bam.records() {
                let record = result.expect("Failed to parse record");
                black_box(&record.name);
                count += 1;
            }
            count
        });
    });

    // Access positions (common for region queries)
    group.bench_function("positions", |b| {
        b.iter(|| {
            let mut bam = BamReader::from_path(black_box(bam_path))
                .expect("Failed to open BAM file");

            let mut mapped_count = 0;
            for result in bam.records() {
                let record = result.expect("Failed to parse record");
                if record.position.is_some() {
                    mapped_count += 1;
                }
                black_box(&record.position);
            }
            mapped_count
        });
    });

    // Full record access (realistic workload)
    group.bench_function("full_access", |b| {
        b.iter(|| {
            let mut bam = BamReader::from_path(black_box(bam_path))
                .expect("Failed to open BAM file");

            let mut stats = (0u64, 0u64, 0u64); // (records, bases, mapped)

            for result in bam.records() {
                let record = result.expect("Failed to parse record");
                stats.0 += 1;
                stats.1 += record.sequence.len() as u64;
                if record.position.is_some() {
                    stats.2 += 1;
                }

                black_box(&record.name);
                black_box(&record.sequence);
                black_box(&record.quality);
                black_box(&record.cigar);
            }

            stats
        });
    });

    group.finish();
}

criterion_group! {
    name = bam_benches;
    config = Criterion::default()
        .sample_size(30)  // N=30 for statistical significance (OPTIMIZATION_RULES.md)
        .measurement_time(std::time::Duration::from_secs(10));
    targets =
        bench_parse_bam_100k,
        bench_parse_header,
        bench_record_throughput,
        bench_access_patterns
}

criterion_main!(bam_benches);
