//! Minimizer Extraction Baseline Benchmark (ASBB Entry 035)
//!
//! This benchmark establishes the performance baseline for minimizer extraction
//! using the current FNV-1a hash + linear scan implementation (from Entry 034).
//!
//! # Purpose
//!
//! Pre-implementation baseline for comparing against ntHash + two stacks approach
//! (simd-minimizers-analysis experiment, GO decision Nov 6, 2025).
//!
//! # Expected Speedup (Entry 035-B comparison)
//!
//! - **Baseline** (Entry 035): ~50-100 Mbp/s (FNV-1a + linear scan)
//! - **Post-implementation** (Entry 035-B): ~400-600 Mbp/s (ntHash + two stacks + SIMD)
//! - **Target**: 4-8× speedup validation
//!
//! # Methodology
//!
//! - **N=30 repetitions** per configuration (statistical rigor, 95% CI)
//! - **4 sequence lengths**: 100bp, 1Kbp, 10Kbp, 100Kbp (scaling validation)
//! - **2 k-mer sizes**: k=21 (typical), k=31 (high-specificity)
//! - **2 window sizes**: w=11 (typical), w=19 (large)
//! - **Total**: 16 configurations × 30 = 480 measurements
//!
//! # Running the Benchmark
//!
//! ```bash
//! cargo bench --bench minimizer_baseline -- --measurement-time 60 --sample-size 30
//! ```
//!
//! # Evidence
//!
//! - **Entry 034** (pilot, N=3): 1.02-1.26× NEON speedup (scalar optimal)
//! - **Entry 035** (baseline, N=30): Full statistical rigor for comparison
//! - **Entry 035-B** (future): Post-implementation validation

use biometal::operations::kmer::extract_minimizers;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use std::collections::hash_map::RandomState;
use std::hash::{BuildHasher, Hash, Hasher};

/// Generate random DNA sequence for benchmarking
///
/// Uses seeded PRNG for reproducibility across runs
fn generate_random_dna(len: usize) -> Vec<u8> {
    let build_hasher = RandomState::new();
    (0..len)
        .map(|i| {
            let mut hasher = build_hasher.build_hasher();
            i.hash(&mut hasher);
            b"ACGT"[(hasher.finish() as usize) % 4]
        })
        .collect()
}

/// Benchmark minimizer extraction baseline (Entry 035)
///
/// Tests all 16 configurations:
/// - k ∈ {21, 31}
/// - w ∈ {11, 19}
/// - sequence length ∈ {100, 1000, 10000, 100000}
fn bench_minimizer_baseline(c: &mut Criterion) {
    let mut group = c.benchmark_group("minimizer_baseline");

    // Parameters (from Entry 035 specification)
    let k_values = [21, 31];
    let w_values = [11, 19];
    let seq_lengths = [100, 1_000, 10_000, 100_000];

    for &k in &k_values {
        for &w in &w_values {
            for &len in &seq_lengths {
                // Generate test sequence
                let seq = generate_random_dna(len);

                // Benchmark ID: k21_w11/100 (k-mer size, window size, sequence length)
                let id = BenchmarkId::new(format!("k{}_w{}", k, w), len);

                // Throughput: base pairs per second
                group.throughput(Throughput::Bytes(len as u64));

                // Benchmark: Extract minimizers
                group.bench_with_input(id, &seq, |b, seq| {
                    b.iter(|| extract_minimizers(black_box(seq), black_box(k), black_box(w)));
                });
            }
        }
    }

    group.finish();
}

/// Benchmark minimizer extraction at production scale
///
/// Additional benchmark for large sequences (1Mbp, 10Mbp) to validate
/// performance at realistic genomic scales
fn bench_minimizer_production_scale(c: &mut Criterion) {
    let mut group = c.benchmark_group("minimizer_production");

    // Typical genomics parameters
    let k = 21;
    let w = 11;
    let seq_lengths = [1_000_000, 10_000_000]; // 1Mbp, 10Mbp

    for &len in &seq_lengths {
        let seq = generate_random_dna(len);
        let id = BenchmarkId::new("k21_w11", len);

        group.throughput(Throughput::Bytes(len as u64));

        group.bench_with_input(id, &seq, |b, seq| {
            b.iter(|| extract_minimizers(black_box(seq), black_box(k), black_box(w)));
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_minimizer_baseline,
    bench_minimizer_production_scale
);
criterion_main!(benches);
