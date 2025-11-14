//! Smith-Waterman alignment benchmarks
//!
//! Compares CPU naive vs Metal GPU implementations across:
//! - Sequence lengths: 100bp, 500bp, 1000bp
//! - Batch sizes: 1, 10, 50, 100, 500, 1000 alignments
//! - Statistical rigor: N=30 samples
//!
//! Expected results (from CUDA literature):
//! - GPU speedup: 10-50× for batch sizes >100
//! - Dispatch overhead: ~1-3ms
//! - Sweet spot: 100-500 alignments per batch

use biometal::alignment::{smith_waterman_naive, ScoringMatrix};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rand::Rng;

/// Generate random DNA sequence of given length
fn generate_sequence(len: usize) -> Vec<u8> {
    let bases = b"ACGT";
    let mut rng = rand::thread_rng();
    (0..len).map(|_| bases[rng.gen_range(0..4)]).collect()
}

/// Benchmark CPU naive implementation - single alignment
fn bench_cpu_naive_single(c: &mut Criterion) {
    let mut group = c.benchmark_group("smith_waterman_cpu_single");
    group.sample_size(30); // N=30 for statistical rigor

    for seq_len in [100, 500, 1000].iter() {
        let query = generate_sequence(*seq_len);
        let reference = generate_sequence(*seq_len);
        let scoring = ScoringMatrix::default();

        group.throughput(Throughput::Elements(1));
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{}bp", seq_len)),
            seq_len,
            |b, _| {
                b.iter(|| {
                    black_box(smith_waterman_naive(
                        black_box(&query),
                        black_box(&reference),
                        black_box(&scoring),
                    ))
                })
            },
        );
    }

    group.finish();
}

/// Benchmark CPU naive implementation - batch processing
fn bench_cpu_naive_batch(c: &mut Criterion) {
    let mut group = c.benchmark_group("smith_waterman_cpu_batch");
    group.sample_size(30);

    let seq_len = 500; // Focus on realistic sequence length
    let scoring = ScoringMatrix::default();

    for batch_size in [10, 50, 100].iter() {
        let queries: Vec<_> = (0..*batch_size)
            .map(|_| generate_sequence(seq_len))
            .collect();
        let references: Vec<_> = (0..*batch_size)
            .map(|_| generate_sequence(seq_len))
            .collect();

        group.throughput(Throughput::Elements(*batch_size as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{}x500bp", batch_size)),
            batch_size,
            |b, _| {
                b.iter(|| {
                    for (q, r) in queries.iter().zip(references.iter()) {
                        black_box(smith_waterman_naive(
                            black_box(q),
                            black_box(r),
                            black_box(&scoring),
                        ));
                    }
                })
            },
        );
    }

    group.finish();
}

/// Benchmark GPU implementation - single alignment
#[cfg(feature = "gpu")]
fn bench_gpu_single(c: &mut Criterion) {
    use biometal::alignment::smith_waterman_gpu;

    let mut group = c.benchmark_group("smith_waterman_gpu_single");
    group.sample_size(30);

    for seq_len in [100, 500, 1000].iter() {
        let query = generate_sequence(*seq_len);
        let reference = generate_sequence(*seq_len);
        let scoring = ScoringMatrix::default();

        group.throughput(Throughput::Elements(1));
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{}bp", seq_len)),
            seq_len,
            |b, _| {
                b.iter(|| {
                    black_box(smith_waterman_gpu(
                        black_box(&query),
                        black_box(&reference),
                        black_box(&scoring),
                    ))
                })
            },
        );
    }

    group.finish();
}

/// Benchmark GPU implementation - batch processing
#[cfg(feature = "gpu")]
fn bench_gpu_batch(c: &mut Criterion) {
    use biometal::alignment::gpu::smith_waterman_batch_gpu;

    let mut group = c.benchmark_group("smith_waterman_gpu_batch");
    group.sample_size(30);

    let seq_len = 500;
    let scoring = ScoringMatrix::default();

    for batch_size in [10, 50, 100, 500, 1000].iter() {
        let queries: Vec<_> = (0..*batch_size)
            .map(|_| generate_sequence(seq_len))
            .collect();
        let references: Vec<_> = (0..*batch_size)
            .map(|_| generate_sequence(seq_len))
            .collect();

        let query_slices: Vec<&[u8]> = queries.iter().map(|s| s.as_slice()).collect();
        let ref_slices: Vec<&[u8]> = references.iter().map(|s| s.as_slice()).collect();

        group.throughput(Throughput::Elements(*batch_size as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{}x500bp", batch_size)),
            batch_size,
            |b, _| {
                b.iter(|| {
                    black_box(smith_waterman_batch_gpu(
                        black_box(&query_slices),
                        black_box(&ref_slices),
                        black_box(&scoring),
                    ))
                })
            },
        );
    }

    group.finish();
}

/// Benchmark dispatch overhead - GPU batch creation
#[cfg(feature = "gpu")]
fn bench_gpu_context_creation(c: &mut Criterion) {
    use biometal::alignment::gpu::GpuAlignmentBatch;

    let mut group = c.benchmark_group("smith_waterman_gpu_overhead");
    group.sample_size(30);

    group.bench_function("batch_creation", |b| {
        b.iter(|| {
            black_box(GpuAlignmentBatch::new())
        })
    });

    group.finish();
}

// Conditional compilation for GPU benchmarks
#[cfg(feature = "gpu")]
criterion_group!(
    benches,
    bench_cpu_naive_single,
    bench_cpu_naive_batch,
    bench_gpu_single,
    bench_gpu_batch,
    bench_gpu_context_creation
);

#[cfg(not(feature = "gpu"))]
criterion_group!(
    benches,
    bench_cpu_naive_single,
    bench_cpu_naive_batch
);

criterion_main!(benches);
