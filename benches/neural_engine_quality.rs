//! Benchmark: Neural Engine vs Traditional Quality Filtering
//!
//! Compares ML-powered quality prediction (Neural Engine) against
//! traditional average quality score filtering.

use biometal::operations::passes_quality_filter;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

#[cfg(feature = "neural-engine")]
use biometal::ml::quality::QualityPredictor;

// Generate test sequences with varying quality
fn generate_test_data(count: usize, avg_quality: u8) -> Vec<(Vec<u8>, Vec<u8>)> {
    use rand::Rng;
    let mut rng = rand::thread_rng();

    (0..count)
        .map(|_| {
            // Generate random 150bp sequence
            let sequence: Vec<u8> = (0..150)
                .map(|_| match rng.gen_range(0..4) {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                })
                .collect();

            // Generate quality scores around avg_quality
            let quality: Vec<u8> = (0..150)
                .map(|_| {
                    let variation = rng.gen_range(-10i16..10i16);
                    let q = (avg_quality as i16 + variation).max(0).min(93) as u8;
                    q
                })
                .collect();

            (sequence, quality)
        })
        .collect()
}

fn bench_traditional_quality_filter(c: &mut Criterion) {
    let mut group = c.benchmark_group("traditional_quality_filter");

    for &count in &[100, 1000, 10000] {
        group.throughput(Throughput::Elements(count));

        // High quality reads (Q30)
        group.bench_with_input(
            BenchmarkId::new("high_quality", count),
            &count,
            |b, &count| {
                let test_data = generate_test_data(count as usize, 30);
                b.iter(|| {
                    for (_sequence, quality) in &test_data {
                        black_box(passes_quality_filter(quality, 20.0));
                    }
                });
            },
        );

        // Low quality reads (Q10)
        group.bench_with_input(
            BenchmarkId::new("low_quality", count),
            &count,
            |b, &count| {
                let test_data = generate_test_data(count as usize, 10);
                b.iter(|| {
                    for (_sequence, quality) in &test_data {
                        black_box(passes_quality_filter(quality, 20.0));
                    }
                });
            },
        );
    }

    group.finish();
}

#[cfg(feature = "neural-engine")]
fn bench_neural_engine_quality_filter(c: &mut Criterion) {
    // Load model
    let model_path = "quality_model.onnx";
    let mut predictor = match QualityPredictor::new(model_path) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("Warning: Could not load ONNX model: {}", e);
            eprintln!("Skipping Neural Engine benchmarks");
            return;
        }
    };

    let mut group = c.benchmark_group("neural_engine_quality_filter");

    for &count in &[100, 1000, 10000] {
        group.throughput(Throughput::Elements(count));

        // High quality reads (Q30)
        group.bench_with_input(
            BenchmarkId::new("high_quality", count),
            &count,
            |b, &count| {
                let test_data = generate_test_data(count as usize, 30);
                b.iter(|| {
                    for (sequence, quality) in &test_data {
                        black_box(predictor.predict_quality(sequence, quality).unwrap());
                    }
                });
            },
        );

        // Low quality reads (Q10)
        group.bench_with_input(
            BenchmarkId::new("low_quality", count),
            &count,
            |b, &count| {
                let test_data = generate_test_data(count as usize, 10);
                b.iter(|| {
                    for (sequence, quality) in &test_data {
                        black_box(predictor.predict_quality(sequence, quality).unwrap());
                    }
                });
            },
        );
    }

    group.finish();
}

#[cfg(feature = "neural-engine")]
fn bench_quality_filter_comparison(c: &mut Criterion) {
    // Load model
    let model_path = "quality_model.onnx";
    let mut predictor = match QualityPredictor::new(model_path) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("Warning: Could not load ONNX model: {}", e);
            return;
        }
    };

    let mut group = c.benchmark_group("quality_filter_comparison");
    group.throughput(Throughput::Elements(1000));

    let test_data = generate_test_data(1000, 25); // Medium quality (Q25)

    group.bench_function("traditional", |b| {
        b.iter(|| {
            for (_sequence, quality) in &test_data {
                black_box(passes_quality_filter(quality, 20.0));
            }
        });
    });

    group.bench_function("neural_engine", |b| {
        b.iter(|| {
            for (sequence, quality) in &test_data {
                black_box(predictor.predict_quality(sequence, quality).unwrap());
            }
        });
    });

    group.finish();
}

// Traditional quality filter benchmarks (always available)
criterion_group!(traditional_benches, bench_traditional_quality_filter);

// Neural Engine benchmarks (only with feature flag)
#[cfg(feature = "neural-engine")]
criterion_group!(
    neural_engine_benches,
    bench_neural_engine_quality_filter,
    bench_quality_filter_comparison
);

#[cfg(feature = "neural-engine")]
criterion_main!(traditional_benches, neural_engine_benches);

#[cfg(not(feature = "neural-engine"))]
criterion_main!(traditional_benches);
