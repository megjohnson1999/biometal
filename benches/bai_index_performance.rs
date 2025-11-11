//! BAI Index Performance Benchmarks
//!
//! Validates O(log n) indexed region query performance vs O(n) full scan.
//!
//! # Benchmarks
//!
//! - `indexed_query_small`: Query chr1:1-1000 using BAI index (O(log n))
//! - `full_scan_small`: Full scan + filter chr1:1-1000 (O(n))
//! - `indexed_query_medium`: Query chr1:1-10000 using BAI index
//! - `full_scan_medium`: Full scan + filter chr1:1-10000
//! - `sequential_baseline`: Sequential read all records (baseline)
//!
//! # Expected Results
//!
//! - Indexed query should be significantly faster than full scan for small regions
//! - Speedup should increase as file size grows relative to query region
//! - For 100K record file, querying ~3% of records (small region):
//!   - Expected speedup: 10-30× (reads only relevant BGZF blocks)
//!
//! # Validation
//!
//! Sample size: N=30 (OPTIMIZATION_RULES.md requirement)
//! Measurement time: 10 seconds per benchmark

use biometal::io::bam::{BamReader, BaiIndex};
use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use std::path::Path;

const BAM_PATH: &str = "tests/data/synthetic_100k.bam";
const BAI_PATH: &str = "tests/data/synthetic_100k.bam.bai";

/// Benchmark indexed region query for small region (chr1:1-1000)
///
/// This uses the BAI index to jump directly to relevant BGZF blocks.
/// Expected: Fast, O(log n) complexity for finding blocks + O(m) for reading matching records.
fn bench_indexed_query_small(c: &mut Criterion) {
    if !Path::new(BAM_PATH).exists() || !Path::new(BAI_PATH).exists() {
        eprintln!("Warning: Test files not found ({})", BAM_PATH);
        eprintln!("Run 'cargo test' first to generate test data.");
        return;
    }

    let index = BaiIndex::from_path(BAI_PATH)
        .expect("Failed to load BAI index");

    c.bench_function("indexed_query_small_region", |b| {
        b.iter(|| {
            let mut query = BamReader::query(
                black_box(BAM_PATH),
                &index,
                "chr1",
                1,
                1000,
            ).expect("Failed to create query");

            let mut count = 0;
            for result in query {
                let record = result.expect("Failed to parse record");
                black_box(&record);
                count += 1;
            }

            count
        });
    });
}

/// Benchmark full scan for small region (chr1:1-1000)
///
/// This reads ALL records sequentially and filters for the target region.
/// Expected: Slow, O(n) complexity - must decompress and parse entire file.
fn bench_full_scan_small(c: &mut Criterion) {
    if !Path::new(BAM_PATH).exists() {
        return;
    }

    c.bench_function("full_scan_small_region", |b| {
        b.iter(|| {
            let mut bam = BamReader::from_path(black_box(BAM_PATH))
                .expect("Failed to open BAM file");

            let mut count = 0;

            // Manually filter for chr1:1-1000 (simulating no index)
            for result in bam.records() {
                let record = result.expect("Failed to parse record");

                // Filter by reference and position
                if record.reference_id == Some(0) {  // chr1
                    if let Some(pos) = record.position {
                        if pos >= 0 && pos < 1000 {
                            black_box(&record);
                            count += 1;
                        }
                    }
                }
            }

            count
        });
    });
}

/// Benchmark indexed region query for medium region (chr1:1-10000)
///
/// Larger region to test how speedup scales with query size.
fn bench_indexed_query_medium(c: &mut Criterion) {
    if !Path::new(BAM_PATH).exists() || !Path::new(BAI_PATH).exists() {
        return;
    }

    let index = BaiIndex::from_path(BAI_PATH)
        .expect("Failed to load BAI index");

    let mut group = c.benchmark_group("indexed_query");

    // Measure throughput in records processed
    group.throughput(Throughput::Elements(30_000)); // ~30K records expected

    group.bench_function("medium_region", |b| {
        b.iter(|| {
            let mut query = BamReader::query(
                black_box(BAM_PATH),
                &index,
                "chr1",
                1,
                10000,
            ).expect("Failed to create query");

            let mut count = 0;
            for result in query {
                let record = result.expect("Failed to parse record");
                black_box(&record);
                count += 1;
            }

            count
        });
    });

    group.finish();
}

/// Benchmark full scan for medium region (chr1:1-10000)
fn bench_full_scan_medium(c: &mut Criterion) {
    if !Path::new(BAM_PATH).exists() {
        return;
    }

    let mut group = c.benchmark_group("full_scan");
    group.throughput(Throughput::Elements(30_000));

    group.bench_function("medium_region", |b| {
        b.iter(|| {
            let mut bam = BamReader::from_path(black_box(BAM_PATH))
                .expect("Failed to open BAM file");

            let mut count = 0;

            for result in bam.records() {
                let record = result.expect("Failed to parse record");

                if record.reference_id == Some(0) {  // chr1
                    if let Some(pos) = record.position {
                        if pos >= 0 && pos < 10000 {
                            black_box(&record);
                            count += 1;
                        }
                    }
                }
            }

            count
        });
    });

    group.finish();
}

/// Benchmark sequential read baseline
///
/// Reads all records without filtering. This is the baseline for full scan.
fn bench_sequential_baseline(c: &mut Criterion) {
    if !Path::new(BAM_PATH).exists() {
        return;
    }

    let mut group = c.benchmark_group("baseline");

    // Get file size for throughput measurement
    let file_size = std::fs::metadata(BAM_PATH)
        .map(|m| m.len())
        .unwrap_or(0);

    group.throughput(Throughput::Bytes(file_size));

    group.bench_function("sequential_read_all", |b| {
        b.iter(|| {
            let mut bam = BamReader::from_path(black_box(BAM_PATH))
                .expect("Failed to open BAM file");

            let mut count = 0;
            for result in bam.records() {
                let record = result.expect("Failed to parse record");
                black_box(&record);
                count += 1;
            }

            count
        });
    });

    group.finish();
}

/// Benchmark index loading overhead
///
/// Measures how long it takes to load and parse the BAI index file.
/// Expected: Very fast (<1ms for typical indices).
fn bench_index_load(c: &mut Criterion) {
    if !Path::new(BAI_PATH).exists() {
        return;
    }

    c.bench_function("index_load", |b| {
        b.iter(|| {
            let index = BaiIndex::from_path(black_box(BAI_PATH))
                .expect("Failed to load BAI index");
            black_box(&index);
        });
    });
}

/// Benchmark query setup overhead
///
/// Measures the overhead of creating an indexed query (bin calculation, chunk lookup).
/// Expected: Negligible (<0.1ms).
fn bench_query_setup(c: &mut Criterion) {
    if !Path::new(BAM_PATH).exists() || !Path::new(BAI_PATH).exists() {
        return;
    }

    let index = BaiIndex::from_path(BAI_PATH)
        .expect("Failed to load BAI index");

    c.bench_function("query_setup", |b| {
        b.iter(|| {
            // Just create the query, don't iterate
            let query = BamReader::query(
                black_box(BAM_PATH),
                &index,
                "chr1",
                1,
                1000,
            ).expect("Failed to create query");

            black_box(query);
        });
    });
}

criterion_group! {
    name = bai_benches;
    config = Criterion::default()
        .sample_size(30)  // N=30 for statistical significance (OPTIMIZATION_RULES.md)
        .measurement_time(std::time::Duration::from_secs(10));
    targets =
        bench_index_load,
        bench_query_setup,
        bench_indexed_query_small,
        bench_full_scan_small,
        bench_indexed_query_medium,
        bench_full_scan_medium,
        bench_sequential_baseline
}

criterion_main!(bai_benches);
