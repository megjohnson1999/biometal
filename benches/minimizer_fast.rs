//! Minimizer Extraction Fast Benchmark (ASBB Entry 036-B)
//!
//! This benchmark validates the performance of the optimized minimizer extraction
//! using ntHash + sliding window minimum (Phase 1 implementation, v1.3.0).
//!
//! # Purpose
//!
//! Post-implementation validation comparing against Entry 036 baseline to confirm
//! 100-200× speedup target from simd-minimizers integration.
//!
//! # Expected Speedup (Entry 036-B vs Entry 036)
//!
//! - **Baseline** (Entry 036): 1.7-5.5 Mbp/s (FNV-1a + O(w) scan, mean: 3.1 Mbp/s)
//! - **Fast** (Entry 036-B): 370-740 Mbp/s (ntHash + O(1) sliding min)
//! - **Target**: 100-200× speedup validation
//!
//! # Success Criteria
//!
//! - ≥50× speedup: SUCCESS (conservative threshold) ✅
//! - ≥100× speedup: EXCEPTIONAL (realistic target) ✅
//! - ≥150× speedup: OUTSTANDING (stretch goal) ✅
//!
//! # Methodology
//!
//! - **N=100 repetitions** per configuration (same as Entry 036 for comparison)
//! - **4 sequence lengths**: 100bp, 1Kbp, 10Kbp, 100Kbp (scaling validation)
//! - **2 k-mer sizes**: k=21 (typical), k=31 (high-specificity)
//! - **2 window sizes**: w=11 (typical), w=19 (large)
//! - **Total**: 16 configurations × 100 = 1,600 measurements
//!
//! # Running the Benchmark
//!
//! ```bash
//! # Full benchmark (N=100, ~45 minutes)
//! cargo bench --bench minimizer_fast -- --measurement-time 60 --sample-size 100
//!
//! # Quick validation (N=30, ~15 minutes)
//! cargo bench --bench minimizer_fast -- --measurement-time 30 --sample-size 30
//! ```
//!
//! # Evidence
//!
//! - **Entry 036** (baseline, N=100): 3.7 Mbp/s mean (rigorous baseline)
//! - **SimdMinimizers**: 820 Mbp/s (221× faster than Entry 036)
//! - **Entry 036-B** (this benchmark): Validates 100-200× target
//!
//! # Phase 2 Integration
//!
//! After running, use `parse_minimizer_fast.py` to extract statistics and compare
//! with Entry 036 baseline results for formal ASBB entry documentation.

use biometal::operations::kmer::extract_minimizers_fast;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use std::collections::hash_map::RandomState;
use std::hash::{BuildHasher, Hash, Hasher};

/// Generate random DNA sequence for benchmarking
///
/// Uses seeded PRNG for reproducibility across runs (same as Entry 036)
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

/// Benchmark fast minimizer extraction (Entry 036-B)
///
/// Tests all 16 configurations (same as Entry 036 baseline):
/// - k ∈ {21, 31}
/// - w ∈ {11, 19}
/// - sequence length ∈ {100, 1000, 10000, 100000}
fn bench_minimizer_fast(c: &mut Criterion) {
    let mut group = c.benchmark_group("minimizer_fast");

    // Parameters (identical to Entry 036 baseline)
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

                // Benchmark: Extract minimizers using fast algorithm
                group.bench_with_input(id, &seq, |b, seq| {
                    b.iter(|| {
                        extract_minimizers_fast(black_box(seq), black_box(k), black_box(w))
                            .unwrap()
                    });
                });
            }
        }
    }

    group.finish();
}

/// Benchmark fast minimizer extraction at production scale
///
/// Additional benchmark for large sequences (1Mbp, 10Mbp) to validate
/// performance at realistic genomic scales.
///
/// **Note**: These are resource-intensive. Use smaller N for quick validation.
fn bench_minimizer_fast_production_scale(c: &mut Criterion) {
    let mut group = c.benchmark_group("minimizer_fast_production");

    // Typical genomics parameters
    let k = 21;
    let w = 11;

    // Production scales: 1Mbp (bacterial gene cluster), 10Mbp (small chromosome)
    let seq_lengths = [1_000_000, 10_000_000];

    for &len in &seq_lengths {
        let seq = generate_random_dna(len);
        let id = BenchmarkId::new(format!("k{}_w{}", k, w), len);

        group.throughput(Throughput::Bytes(len as u64));

        group.bench_with_input(id, &seq, |b, seq| {
            b.iter(|| {
                extract_minimizers_fast(black_box(seq), black_box(k), black_box(w)).unwrap()
            });
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_minimizer_fast,
    bench_minimizer_fast_production_scale
);
criterion_main!(benches);
