//! Quick validation benchmark: CAF column-selective vs full-record access
//!
//! This benchmark tests whether CAF's column-selective reading provides
//! speedups for quality filtering workflows.

use caf::block::{BlockBuilder, BlockReader, AlignmentRecord};
use caf::neon::filter_records_by_quality;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

/// Generate test CAF block with specified number of records
fn generate_test_caf_block(num_records: usize) -> caf::types::CafBlock {
    let mut builder = BlockBuilder::new(0, num_records as u32);

    for i in 0..num_records {
        let record = AlignmentRecord {
            ref_id: 0,
            position: (i * 100) as i32,
            mapq: 60,
            flags: 99,
            sequence: vec![b'A'; 100], // 100 base sequences
            qualities: vec![b'I'; 100], // Q=40 for all (high quality)
            cigar: vec![(0, 100)], // 100M
            read_name: format!("read{}", i).into_bytes(),
            mate_ref_id: 0,
            mate_position: (i * 100 + 100) as i32,
            template_length: 200,
        };

        builder.add_record(record).unwrap();
    }

    builder.build().unwrap()
}

/// CAF: Full record reconstruction + quality filtering (like BAM does)
fn bench_caf_full_record_access(block: &caf::types::CafBlock, min_quality: f64) -> usize {
    let reader = BlockReader::new(block.clone()).unwrap();

    let mut passing_count = 0;

    // Reconstruct full records (pays cost of decoding ALL columns)
    for i in 0..reader.len() {
        let record = reader.get_record(i).unwrap();

        // Calculate mean quality
        if !record.qualities.is_empty() {
            let sum: u32 = record.qualities.iter().map(|&q| (q - 33) as u32).sum();
            let mean = sum as f64 / record.qualities.len() as f64;

            if mean >= min_quality {
                passing_count += 1;
            }
        }
    }

    passing_count
}

/// CAF: Column-selective quality filtering (our hypothesis - should be faster)
fn bench_caf_column_selective_access(block: &caf::types::CafBlock, min_quality: f64) -> usize {
    let reader = BlockReader::new(block.clone()).unwrap();

    // Access ONLY quality column (column-selective reading)
    let (qualities_data, qual_offsets) = reader.qualities();

    // Filter using NEON-optimized function
    let passing_indices = filter_records_by_quality(qualities_data, qual_offsets, min_quality);

    passing_indices.len()
}

/// Manual quality filtering with column data (baseline)
fn bench_caf_column_manual_filter(block: &caf::types::CafBlock, min_quality: f64) -> usize {
    let reader = BlockReader::new(block.clone()).unwrap();

    // Access ONLY quality column
    let (qualities_data, qual_offsets) = reader.qualities();

    let mut passing_count = 0;

    // Manual filtering (no NEON, just scalar)
    for i in 0..qual_offsets.len() - 1 {
        let start = qual_offsets[i] as usize;
        let end = qual_offsets[i + 1] as usize;
        let quality_slice = &qualities_data[start..end];

        if !quality_slice.is_empty() {
            let sum: u32 = quality_slice.iter().map(|&q| (q - 33) as u32).sum();
            let mean = sum as f64 / quality_slice.len() as f64;

            if mean >= min_quality {
                passing_count += 1;
            }
        }
    }

    passing_count
}

fn bench_column_selective_advantage(c: &mut Criterion) {
    let mut group = c.benchmark_group("column_selective_vs_full_record");

    // Test with different dataset sizes
    for size in [1_000, 10_000].iter() {
        let caf_block = generate_test_caf_block(*size);

        // Full record reconstruction (decode ALL columns)
        group.bench_with_input(
            BenchmarkId::new("full_record_access", size),
            size,
            |b, _| b.iter(|| bench_caf_full_record_access(&caf_block, 30.0)),
        );

        // Column-selective manual (decode quality ONLY, scalar filter)
        group.bench_with_input(
            BenchmarkId::new("column_selective_manual", size),
            size,
            |b, _| b.iter(|| bench_caf_column_manual_filter(&caf_block, 30.0)),
        );

        // Column-selective NEON (decode quality ONLY, NEON filter)
        group.bench_with_input(
            BenchmarkId::new("column_selective_neon", size),
            size,
            |b, _| b.iter(|| bench_caf_column_selective_access(&caf_block, 30.0)),
        );
    }

    group.finish();
}

criterion_group!(benches, bench_column_selective_advantage);
criterion_main!(benches);
