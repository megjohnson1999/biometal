//! Benchmark for Rule 2: Block-based processing validation
//!
//! # Evidence Base
//!
//! Entry 027 predicts:
//! - Record-by-record: 82-86% NEON overhead
//! - Block-based: 4-8% overhead
//! - Expected speedup: ~14×
//!
//! This benchmark validates the prediction with N=30 samples.

use biometal::operations::base_counting::count_bases;
use biometal::operations::block::{count_bases_block, gc_content_block, mean_quality_block};
use biometal::operations::gc_content::gc_content;
use biometal::operations::quality_filter::mean_quality;
use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use std::hint::black_box as hint_black_box;

/// Generate test sequences of specified length
fn generate_sequences(count: usize, length: usize) -> Vec<Vec<u8>> {
    (0..count)
        .map(|i| {
            let mut seq = Vec::with_capacity(length);
            for j in 0..length {
                let base = match (i + j) % 4 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                };
                seq.push(base);
            }
            seq
        })
        .collect()
}

/// Generate test quality strings of specified length
fn generate_qualities(count: usize, length: usize) -> Vec<Vec<u8>> {
    (0..count)
        .map(|i| {
            let mut qual = Vec::with_capacity(length);
            for j in 0..length {
                // Phred scores 20-40 (ASCII 53-73)
                let q = 53 + ((i + j) % 21) as u8;
                qual.push(q);
            }
            qual
        })
        .collect()
}

/// Benchmark base counting: per-record vs block
fn bench_base_counting(c: &mut Criterion) {
    let mut group = c.benchmark_group("base_counting");

    // Test with different block sizes
    for size in [100, 1_000, 10_000].iter() {
        let sequences = generate_sequences(*size, 150);
        let seq_refs: Vec<&[u8]> = sequences.iter().map(|s| s.as_slice()).collect();

        // Per-record processing (current API)
        group.bench_with_input(
            BenchmarkId::new("per_record", size),
            size,
            |b, _| {
                b.iter(|| {
                    for seq in &sequences {
                        hint_black_box(count_bases(black_box(seq)));
                    }
                });
            },
        );

        // Block processing (Rule 2 optimization)
        group.bench_with_input(
            BenchmarkId::new("block", size),
            size,
            |b, _| {
                b.iter(|| {
                    hint_black_box(count_bases_block(black_box(&seq_refs)));
                });
            },
        );
    }

    group.finish();
}

/// Benchmark GC content: per-record vs block
fn bench_gc_content(c: &mut Criterion) {
    let mut group = c.benchmark_group("gc_content");

    for size in [100, 1_000, 10_000].iter() {
        let sequences = generate_sequences(*size, 150);
        let seq_refs: Vec<&[u8]> = sequences.iter().map(|s| s.as_slice()).collect();

        group.bench_with_input(
            BenchmarkId::new("per_record", size),
            size,
            |b, _| {
                b.iter(|| {
                    for seq in &sequences {
                        hint_black_box(gc_content(black_box(seq)));
                    }
                });
            },
        );

        group.bench_with_input(
            BenchmarkId::new("block", size),
            size,
            |b, _| {
                b.iter(|| {
                    hint_black_box(gc_content_block(black_box(&seq_refs)));
                });
            },
        );
    }

    group.finish();
}

/// Benchmark mean quality: per-record vs block
fn bench_mean_quality(c: &mut Criterion) {
    let mut group = c.benchmark_group("mean_quality");

    for size in [100, 1_000, 10_000].iter() {
        let qualities = generate_qualities(*size, 150);
        let qual_refs: Vec<&[u8]> = qualities.iter().map(|q| q.as_slice()).collect();

        group.bench_with_input(
            BenchmarkId::new("per_record", size),
            size,
            |b, _| {
                b.iter(|| {
                    for qual in &qualities {
                        hint_black_box(mean_quality(black_box(qual)));
                    }
                });
            },
        );

        group.bench_with_input(
            BenchmarkId::new("block", size),
            size,
            |b, _| {
                b.iter(|| {
                    hint_black_box(mean_quality_block(black_box(&qual_refs)));
                });
            },
        );
    }

    group.finish();
}

/// Benchmark realistic workflow: FASTQ QC pipeline
fn bench_realistic_qc_workflow(c: &mut Criterion) {
    let mut group = c.benchmark_group("qc_workflow");

    let size = 10_000;  // Typical block size
    let sequences = generate_sequences(size, 150);
    let qualities = generate_qualities(size, 150);
    let seq_refs: Vec<&[u8]> = sequences.iter().map(|s| s.as_slice()).collect();
    let qual_refs: Vec<&[u8]> = qualities.iter().map(|q| q.as_slice()).collect();

    // Per-record QC pipeline
    group.bench_function("per_record_pipeline", |b| {
        b.iter(|| {
            let mut passed = 0;
            for (seq, qual) in sequences.iter().zip(qualities.iter()) {
                let gc = gc_content(black_box(seq));
                let mean_q = mean_quality(black_box(qual));
                let counts = count_bases(black_box(seq));

                // Typical QC filters
                if gc >= 0.3 && gc <= 0.7 && mean_q >= 20.0 {
                    passed += 1;
                }

                hint_black_box((gc, mean_q, counts));
            }
            hint_black_box(passed);
        });
    });

    // Block-based QC pipeline (Rule 2)
    group.bench_function("block_pipeline", |b| {
        b.iter(|| {
            let gc_values = gc_content_block(black_box(&seq_refs));
            let mean_quals = mean_quality_block(black_box(&qual_refs));
            let counts = count_bases_block(black_box(&seq_refs));

            let mut passed = 0;
            for (gc, mean_q) in gc_values.iter().zip(mean_quals.iter()) {
                if *gc >= 0.3 && *gc <= 0.7 && *mean_q >= 20.0 {
                    passed += 1;
                }
            }

            hint_black_box((gc_values, mean_quals, counts, passed));
        });
    });

    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(30);  // N=30 for statistical rigor
    targets = bench_base_counting, bench_gc_content, bench_mean_quality, bench_realistic_qc_workflow
}

criterion_main!(benches);
