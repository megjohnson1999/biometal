//! Neural Engine Read Quality Prediction Example
//!
//! This example demonstrates using Apple Neural Engine for ML-powered
//! FASTQ read quality prediction.
//!
//! # Requirements
//!
//! - macOS with Apple Silicon (M1/M2/M3/M4)
//! - Trained ONNX model (quality_model.onnx)
//! - neural-engine feature enabled
//!
//! # Usage
//!
//! ```bash
//! # Train model first
//! python research/neural-engine/train_quality_model.py \
//!     --fastq tests/data/sample.fastq.gz \
//!     --output quality_model.onnx
//!
//! # Run example
//! cargo run --example neural_quality --features neural-engine
//! ```

#[cfg(not(feature = "neural-engine"))]
fn main() {
    eprintln!("Error: This example requires the 'neural-engine' feature");
    eprintln!("Run with: cargo run --example neural_quality --features neural-engine");
    std::process::exit(1);
}

#[cfg(feature = "neural-engine")]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    use biometal::ml::quality::QualityPredictor;
    use biometal::ml::neural_engine::is_neural_engine_available;
    use std::time::Instant;

    println!("{}", "=".repeat(60));
    println!("Neural Engine Read Quality Prediction Example");
    println!("{}", "=".repeat(60));

    // Check Neural Engine availability
    if !is_neural_engine_available() {
        eprintln!("\nWarning: Neural Engine not available on this platform");
        eprintln!("This example requires macOS with Apple Silicon");
        return Ok(());
    }

    println!("\n✅ Neural Engine available!");

    // Load trained model
    let model_path = "quality_model.onnx";
    println!("\nLoading ONNX model: {}", model_path);

    let mut predictor = match QualityPredictor::new(model_path) {
        Ok(p) => {
            println!("✅ Model loaded successfully");
            p
        }
        Err(e) => {
            eprintln!("\n❌ Failed to load model: {}", e);
            eprintln!("\nPlease train the model first:");
            eprintln!("  python research/neural-engine/train_quality_model.py \\");
            eprintln!("      --fastq tests/data/sample.fastq.gz \\");
            eprintln!("      --output quality_model.onnx");
            return Err(e.into());
        }
    };

    // Test cases
    println!("\n{}", "=".repeat(60));
    println!("Running Test Cases");
    println!("{}", "=".repeat(60));

    // Create quality score vectors (need to outlive the test_cases vec)
    let high_quality = vec![40u8; 64];  // Q40
    let low_quality = vec![0u8; 64];    // Q0
    let medium_quality = vec![20u8; 64]; // Q20

    let test_cases = vec![
        (
            "High quality read (should PASS)",
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".as_slice(),
            &high_quality,
            true,
        ),
        (
            "Low quality read (should FAIL)",
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".as_slice(),
            &low_quality,
            false,
        ),
        (
            "Medium quality read (borderline)",
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".as_slice(),
            &medium_quality,
            true,
        ),
    ];

    for (i, (name, sequence, quality, expected)) in test_cases.iter().enumerate() {
        println!("\nTest {}: {}", i + 1, name);
        println!("Sequence: {}...", std::str::from_utf8(&sequence[..40]).unwrap());
        println!("Quality:  {} scores (avg: {:.1})", quality.len(), quality.iter().map(|&q| q as f32).sum::<f32>() / quality.len() as f32);

        // Predict with timing
        let start = Instant::now();
        let (passes, probability) = predictor.predict_quality_with_score(sequence, quality)?;
        let elapsed = start.elapsed();

        // Display results
        let status = if passes { "✅ PASS" } else { "❌ FAIL" };
        println!("Probability: {:.4}", probability);
        println!("Prediction:  {}", status);
        println!("Expected:    {}", if *expected { "✅ PASS" } else { "❌ FAIL" });
        println!("Latency:     {:?}", elapsed);

        if passes == *expected {
            println!("Result:      ✅ CORRECT");
        } else {
            println!("Result:      ⚠️  MISMATCH");
        }
    }

    // Batch inference benchmark
    println!("\n{}", "=".repeat(60));
    println!("Batch Inference Benchmark");
    println!("{}", "=".repeat(60));

    let num_reads = 1000;
    let test_sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let test_quality = vec![30u8; 64]; // Q30

    println!("\nProcessing {} reads...", num_reads);

    let start = Instant::now();
    let mut num_pass = 0;

    for _ in 0..num_reads {
        let passes = predictor.predict_quality(test_sequence, &test_quality)?;
        if passes {
            num_pass += 1;
        }
    }

    let elapsed = start.elapsed();
    let reads_per_sec = num_reads as f64 / elapsed.as_secs_f64();
    let latency_per_read = elapsed / num_reads;

    println!("\nResults:");
    println!("  Total reads:     {}", num_reads);
    println!("  Passed:          {} ({:.1}%)", num_pass, (num_pass as f64 / num_reads as f64) * 100.0);
    println!("  Total time:      {:?}", elapsed);
    println!("  Throughput:      {:.0} reads/sec", reads_per_sec);
    println!("  Latency/read:    {:?}", latency_per_read);

    println!("\n{}", "=".repeat(60));
    println!("Example Complete!");
    println!("{}", "=".repeat(60));

    println!("\nNext steps:");
    println!("1. Benchmark against traditional quality_filter:");
    println!("   cargo bench --features neural-engine");
    println!("2. Integrate into FASTQ processing pipeline");
    println!("3. Experiment with different quality thresholds");

    Ok(())
}
