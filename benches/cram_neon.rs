//! Benchmarks for CRAM ARM NEON optimizations (Phase 3).
//!
//! Validates NEON speedup claims:
//! - Reference comparison: 10-15× speedup target
//! - Base counting: 16-25× speedup target (Rule 1)
//! - Quality delta decoding: 15-20× speedup target

use biometal::io::cram::neon::{
    compare_to_reference_neon, count_bases_neon, decode_quality_deltas_neon,
};
use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use rand::Rng;

/// Generate random DNA sequence for benchmarking
fn generate_sequence(length: usize) -> Vec<u8> {
    let bases = b"ACGT";
    let mut rng = rand::thread_rng();
    (0..length)
        .map(|_| bases[rng.gen_range(0..4)])
        .collect()
}

/// Generate sequence with random mismatches
fn generate_sequence_with_mismatches(reference: &[u8], mismatch_rate: f64) -> Vec<u8> {
    let bases = b"ACGT";
    let mut rng = rand::thread_rng();
    reference
        .iter()
        .map(|&b| {
            if rng.gen::<f64>() < mismatch_rate {
                // Random mismatch
                bases[rng.gen_range(0..4)]
            } else {
                b
            }
        })
        .collect()
}

/// Benchmark reference comparison (NEON vs scalar)
fn bench_reference_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("cram_reference_comparison");

    for length in [100, 1_000, 10_000, 100_000].iter() {
        let reference = generate_sequence(*length);
        let read_5pct = generate_sequence_with_mismatches(&reference, 0.05);  // 5% mismatches (typical)
        let read_1pct = generate_sequence_with_mismatches(&reference, 0.01);  // 1% mismatches (good quality)

        group.bench_with_input(
            BenchmarkId::new("neon_5pct_mismatch", length),
            length,
            |b, _| {
                b.iter(|| {
                    black_box(compare_to_reference_neon(
                        black_box(&read_5pct),
                        black_box(&reference),
                    ))
                });
            },
        );

        group.bench_with_input(
            BenchmarkId::new("neon_1pct_mismatch", length),
            length,
            |b, _| {
                b.iter(|| {
                    black_box(compare_to_reference_neon(
                        black_box(&read_1pct),
                        black_box(&reference),
                    ))
                });
            },
        );

        // Scalar baseline
        group.bench_with_input(
            BenchmarkId::new("scalar_5pct_mismatch", length),
            length,
            |b, _| {
                b.iter(|| {
                    let mut mismatches = Vec::new();
                    for i in 0..black_box(&read_5pct).len() {
                        if read_5pct[i] != reference[i] {
                            mismatches.push((i, read_5pct[i], reference[i]));
                        }
                    }
                    black_box(mismatches)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark base counting (NEON vs scalar)
fn bench_base_counting(c: &mut Criterion) {
    let mut group = c.benchmark_group("cram_base_counting");

    for length in [100, 1_000, 10_000, 100_000].iter() {
        let sequence = generate_sequence(*length);

        group.bench_with_input(
            BenchmarkId::new("neon", length),
            length,
            |b, _| {
                b.iter(|| {
                    black_box(count_bases_neon(black_box(&sequence)))
                });
            },
        );

        // Scalar baseline
        group.bench_with_input(
            BenchmarkId::new("scalar", length),
            length,
            |b, _| {
                b.iter(|| {
                    let mut counts = (0, 0, 0, 0, 0);
                    for &base in black_box(&sequence) {
                        match base {
                            b'A' => counts.0 += 1,
                            b'C' => counts.1 += 1,
                            b'G' => counts.2 += 1,
                            b'T' => counts.3 += 1,
                            b'N' => counts.4 += 1,
                            _ => {}
                        }
                    }
                    black_box(counts)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark quality delta decoding (NEON vs scalar)
fn bench_quality_deltas(c: &mut Criterion) {
    let mut group = c.benchmark_group("cram_quality_deltas");

    for length in [100, 1_000, 10_000, 100_000].iter() {
        let mut rng = rand::thread_rng();
        let deltas: Vec<i8> = (0..*length)
            .map(|_| rng.gen_range(-5..=5))  // Typical delta range
            .collect();

        group.bench_with_input(
            BenchmarkId::new("neon", length),
            length,
            |b, _| {
                b.iter(|| {
                    black_box(decode_quality_deltas_neon(
                        black_box(&deltas),
                        30,  // Typical initial quality
                    ))
                });
            },
        );

        // Scalar baseline
        group.bench_with_input(
            BenchmarkId::new("scalar", length),
            length,
            |b, _| {
                b.iter(|| {
                    let mut qualities = Vec::with_capacity(deltas.len());
                    let mut current = 30i16;
                    for &delta in black_box(&deltas) {
                        current = (current + delta as i16).clamp(0, 93);
                        qualities.push((current as u8).saturating_add(33));
                    }
                    black_box(qualities)
                });
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_reference_comparison,
    bench_base_counting,
    bench_quality_deltas
);
criterion_main!(benches);
