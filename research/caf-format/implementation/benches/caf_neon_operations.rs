//! Benchmarks for CAF NEON operations
//!
//! Validates target speedups from OPTIMIZATION_RULES.md:
//! - Base counting: 25× speedup (Entry 020-025)
//! - Quality filtering: 25× speedup (Entry 020-025)
//! - MAPQ filtering: 16× speedup (Entry 020)
//!
//! Run with: cargo bench --bench caf_neon_operations

use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use caf::neon::*;

/// Generate random DNA sequence
fn generate_sequence(length: usize) -> Vec<u8> {
    use rand::Rng;
    let mut rng = rand::thread_rng();
    let bases = [b'A', b'C', b'G', b'T'];
    (0..length)
        .map(|_| bases[rng.gen_range(0..4)])
        .collect()
}

/// Generate random quality scores
fn generate_quality_scores(length: usize) -> Vec<u8> {
    use rand::Rng;
    let mut rng = rand::thread_rng();
    (0..length)
        .map(|_| rng.gen_range(33..73)) // Phred+33: Q0-40
        .collect()
}

/// Generate random MAPQ values
fn generate_mapq_values(length: usize) -> Vec<u8> {
    use rand::Rng;
    let mut rng = rand::thread_rng();
    (0..length).map(|_| rng.gen_range(0..60)).collect()
}

fn bench_base_counting(c: &mut Criterion) {
    let mut group = c.benchmark_group("base_counting");

    for size in [1_000, 10_000, 100_000].iter() {
        let sequence = generate_sequence(*size);

        // NEON benchmark (only on ARM)
        #[cfg(target_arch = "aarch64")]
        group.bench_with_input(
            BenchmarkId::new("neon", size),
            &sequence,
            |b, seq| {
                b.iter(|| unsafe { count_bases_neon(black_box(seq)) });
            },
        );

        // Scalar benchmark (always available)
        group.bench_with_input(
            BenchmarkId::new("scalar", size),
            &sequence,
            |b, seq| {
                b.iter(|| count_bases_scalar(black_box(seq)));
            },
        );

        // Public API benchmark (platform-appropriate)
        group.bench_with_input(
            BenchmarkId::new("auto", size),
            &sequence,
            |b, seq| {
                b.iter(|| count_bases(black_box(seq)));
            },
        );
    }

    group.finish();
}

fn bench_quality_filtering(c: &mut Criterion) {
    let mut group = c.benchmark_group("quality_filtering");

    for size in [1_000, 10_000, 100_000].iter() {
        let quality_scores = generate_quality_scores(*size);

        // NEON benchmark (only on ARM)
        #[cfg(target_arch = "aarch64")]
        group.bench_with_input(
            BenchmarkId::new("neon", size),
            &quality_scores,
            |b, qual| {
                b.iter(|| unsafe { mean_quality_neon(black_box(qual)) });
            },
        );

        // Scalar benchmark (always available)
        group.bench_with_input(
            BenchmarkId::new("scalar", size),
            &quality_scores,
            |b, qual| {
                b.iter(|| mean_quality_scalar(black_box(qual)));
            },
        );

        // Public API benchmark (platform-appropriate)
        group.bench_with_input(
            BenchmarkId::new("auto", size),
            &quality_scores,
            |b, qual| {
                b.iter(|| mean_quality(black_box(qual)));
            },
        );
    }

    group.finish();
}

fn bench_mapq_filtering(c: &mut Criterion) {
    let mut group = c.benchmark_group("mapq_filtering");

    for size in [1_000, 10_000, 100_000].iter() {
        let mapq_values = generate_mapq_values(*size);

        // NEON benchmark (only on ARM)
        #[cfg(target_arch = "aarch64")]
        group.bench_with_input(
            BenchmarkId::new("neon", size),
            &mapq_values,
            |b, mapqs| {
                b.iter(|| unsafe { count_high_mapq_neon(black_box(mapqs), 30) });
            },
        );

        // Scalar benchmark (always available)
        group.bench_with_input(
            BenchmarkId::new("scalar", size),
            &mapq_values,
            |b, mapqs| {
                b.iter(|| count_high_mapq_scalar(black_box(mapqs), 30));
            },
        );

        // Public API benchmark (platform-appropriate)
        group.bench_with_input(
            BenchmarkId::new("auto", size),
            &mapq_values,
            |b, mapqs| {
                b.iter(|| count_high_mapq(black_box(mapqs), 30));
            },
        );
    }

    group.finish();
}

fn bench_record_filtering(c: &mut Criterion) {
    let mut group = c.benchmark_group("record_filtering");

    for num_records in [1_000, 10_000].iter() {
        // Generate quality data: 100 bases per record
        let bases_per_record = 100;
        let total_bases = num_records * bases_per_record;
        let qualities = generate_quality_scores(total_bases);
        let offsets: Vec<u32> = (0..=*num_records)
            .map(|i| (i * bases_per_record) as u32)
            .collect();

        group.bench_with_input(
            BenchmarkId::new("filter_by_quality_Q30", num_records),
            &(&qualities, &offsets),
            |b, (qual, off)| {
                b.iter(|| {
                    filter_records_by_quality(black_box(qual), black_box(off), 30.0)
                });
            },
        );

        // MAPQ filtering
        let mapq_values = generate_mapq_values(*num_records);

        group.bench_with_input(
            BenchmarkId::new("filter_by_mapq_30", num_records),
            &mapq_values,
            |b, mapqs| {
                b.iter(|| filter_records_by_mapq(black_box(mapqs), 30));
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_base_counting,
    bench_quality_filtering,
    bench_mapq_filtering,
    bench_record_filtering
);
criterion_main!(benches);
