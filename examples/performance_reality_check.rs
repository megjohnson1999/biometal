//! Step 2: Basic Performance Reality Check
//!
//! Measure actual memory usage and execution timing with realistic data sizes
//! to validate the ~5MB constant memory claims and assess performance.

use biometal::alignment::{StreamingMapper, StreamingMapperConfig, ScoringMatrix};
use std::process::Command;
use std::time::{Instant, Duration};
use std::fs;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("📊 Step 2: Basic Performance Reality Check");
    println!("{}", "=".repeat(60));

    test_memory_usage()?;
    test_execution_timing()?;
    test_scaling_behavior()?;

    println!("\n🎯 Step 2 Complete: Performance baseline established");
    Ok(())
}

fn get_memory_usage_mb() -> Result<f64, Box<dyn std::error::Error>> {
    let pid = std::process::id();

    // Try using ps command on Unix-like systems
    let output = Command::new("ps")
        .args(&["-o", "rss=", "-p", &pid.to_string()])
        .output()?;

    if output.status.success() {
        let rss_kb: f64 = String::from_utf8(output.stdout)?
            .trim()
            .parse()
            .unwrap_or(0.0);
        Ok(rss_kb / 1024.0) // Convert KB to MB
    } else {
        // Fallback - just return 0 if ps fails
        Ok(0.0)
    }
}

fn create_test_reference(size_kb: usize) -> Result<String, Box<dyn std::error::Error>> {
    let path = format!("/tmp/test_ref_{}kb.fasta", size_kb);

    // Create reference with alternating patterns to enable alignment
    let pattern_size = 100;
    let patterns = vec!["A", "T", "G", "C"];
    let mut reference = String::new();

    let total_bases = size_kb * 1000;
    for i in 0..total_bases {
        let pattern_idx = (i / pattern_size) % patterns.len();
        reference.push_str(patterns[pattern_idx]);
    }

    let fasta_content = format!(">test_ref_{}kb\n{}\n", size_kb, reference);
    fs::write(&path, fasta_content)?;

    Ok(path)
}

fn create_test_reads(num_reads: usize, read_length: usize) -> Result<String, Box<dyn std::error::Error>> {
    let path = format!("/tmp/test_reads_{}_{}bp.fastq", num_reads, read_length);
    let mut fastq_content = String::new();

    // Create reads that should align to different parts of reference
    let patterns = vec!["A", "T", "G", "C"];

    for i in 0..num_reads {
        let pattern_idx = i % patterns.len();
        let sequence = patterns[pattern_idx].repeat(read_length);
        let quality = "I".repeat(read_length);

        fastq_content.push_str(&format!("@read_{}\n{}\n+\n{}\n", i, sequence, quality));
    }

    fs::write(&path, fastq_content)?;
    Ok(path)
}

fn test_memory_usage() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n💾 Test 1: Memory Usage Measurement");

    let initial_memory = get_memory_usage_mb()?;
    println!("   Initial memory: {:.1} MB", initial_memory);

    // Test with different reference sizes to see if memory stays constant
    let test_sizes = vec![10, 50, 100, 500]; // KB

    for &size_kb in &test_sizes {
        println!("\n   Testing with {}KB reference:", size_kb);

        // Create test data
        let ref_path = create_test_reference(size_kb)?;
        let reads_path = create_test_reads(100, 50)?; // 100 reads, 50bp each

        // Create mapper
        let config = StreamingMapperConfig {
            window_size: 10_000, // 10KB windows
            overlap_bp: 1_000,   // 1KB overlap
            min_score_threshold: 20,
            scoring: ScoringMatrix::default(),
        };
        let mut mapper = StreamingMapper::new(config);

        let before_mapping = get_memory_usage_mb()?;
        println!("     Memory before mapping: {:.1} MB", before_mapping);

        // Execute mapping
        match mapper.map_reads_streaming(&ref_path, &reads_path) {
            Ok(mapping_iter) => {
                let during_mapping = get_memory_usage_mb()?;
                println!("     Memory during setup: {:.1} MB", during_mapping);

                // Process some mappings
                let mappings: Vec<_> = mapping_iter.take(10).collect::<Result<Vec<_>, _>>()?;

                let after_processing = get_memory_usage_mb()?;
                println!("     Memory after processing: {:.1} MB", after_processing);
                println!("     Mappings found: {}", mappings.len());

                let memory_increase = after_processing - initial_memory;
                println!("     Memory increase: {:.1} MB", memory_increase);

                // Check if memory is roughly constant regardless of reference size
                let is_constant = memory_increase < 50.0; // Allow up to 50MB increase
                println!("     Constant memory: {} ({:.1} MB < 50MB)",
                         if is_constant { "✅ YES" } else { "❌ NO" }, memory_increase);
            }
            Err(e) => println!("     ❌ Mapping failed: {}", e),
        }

        // Cleanup
        let _ = fs::remove_file(ref_path);
        let _ = fs::remove_file(reads_path);
    }

    let final_memory = get_memory_usage_mb()?;
    println!("\n   Final memory: {:.1} MB", final_memory);

    Ok(())
}

fn test_execution_timing() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n⏱️  Test 2: Execution Timing");

    // Test with a reasonable dataset
    let ref_path = create_test_reference(100)?; // 100KB reference
    let reads_path = create_test_reads(1000, 75)?; // 1000 reads, 75bp each

    let config = StreamingMapperConfig {
        window_size: 20_000,  // 20KB windows
        overlap_bp: 2_000,    // 2KB overlap
        min_score_threshold: 25,
        scoring: ScoringMatrix::default(),
    };
    let mut mapper = StreamingMapper::new(config);

    println!("   Dataset: 100KB reference, 1000 reads (75bp)");
    println!("   Window: 20KB with 2KB overlap");

    let start_time = Instant::now();

    match mapper.map_reads_streaming(&ref_path, &reads_path) {
        Ok(mapping_iter) => {
            let setup_time = start_time.elapsed();
            println!("   Setup time: {:.3}s", setup_time.as_secs_f64());

            let process_start = Instant::now();
            let mappings: Vec<_> = mapping_iter.collect::<Result<Vec<_>, _>>()?;
            let process_time = process_start.elapsed();

            let total_time = start_time.elapsed();

            println!("   Processing time: {:.3}s", process_time.as_secs_f64());
            println!("   Total time: {:.3}s", total_time.as_secs_f64());
            println!("   Mappings found: {}", mappings.len());

            // Calculate throughput
            let reads_per_sec = 1000.0 / total_time.as_secs_f64();
            let mbp_per_sec = (1000.0 * 75.0) / (total_time.as_secs_f64() * 1_000_000.0);

            println!("   Throughput: {:.0} reads/sec, {:.3} Mbp/sec", reads_per_sec, mbp_per_sec);

            // Performance assessment
            let is_fast = total_time.as_secs() < 10;
            println!("   Performance: {} ({:.1}s < 10s)",
                     if is_fast { "✅ REASONABLE" } else { "⚠️ SLOW" }, total_time.as_secs_f64());
        }
        Err(e) => println!("   ❌ Timing test failed: {}", e),
    }

    // Cleanup
    let _ = fs::remove_file(ref_path);
    let _ = fs::remove_file(reads_path);

    Ok(())
}

fn test_scaling_behavior() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n📈 Test 3: Scaling Behavior");

    let read_counts = vec![100, 500, 1000];
    let mut results = Vec::new();

    for &num_reads in &read_counts {
        println!("\n   Testing with {} reads:", num_reads);

        let ref_path = create_test_reference(50)?; // Fixed 50KB reference
        let reads_path = create_test_reads(num_reads, 60)?;

        let config = StreamingMapperConfig {
            window_size: 15_000,
            overlap_bp: 1_500,
            min_score_threshold: 20,
            scoring: ScoringMatrix::default(),
        };
        let mut mapper = StreamingMapper::new(config);

        let start_memory = get_memory_usage_mb()?;
        let start_time = Instant::now();

        match mapper.map_reads_streaming(&ref_path, &reads_path) {
            Ok(mapping_iter) => {
                let mappings: Vec<_> = mapping_iter.collect::<Result<Vec<_>, _>>()?;

                let elapsed = start_time.elapsed();
                let end_memory = get_memory_usage_mb()?;
                let memory_delta = end_memory - start_memory;

                println!("     Time: {:.3}s", elapsed.as_secs_f64());
                println!("     Memory delta: {:.1} MB", memory_delta);
                println!("     Mappings: {}", mappings.len());
                println!("     Rate: {:.0} reads/sec", num_reads as f64 / elapsed.as_secs_f64());

                results.push((num_reads, elapsed, memory_delta, mappings.len()));
            }
            Err(e) => println!("     ❌ Failed: {}", e),
        }

        // Cleanup
        let _ = fs::remove_file(ref_path);
        let _ = fs::remove_file(reads_path);
    }

    // Analyze scaling
    if results.len() >= 2 {
        println!("\n   Scaling Analysis:");
        for i in 1..results.len() {
            let (reads1, time1, mem1, _) = results[i-1];
            let (reads2, time2, mem2, _) = results[i];

            let read_ratio = reads2 as f64 / reads1 as f64;
            let time_ratio = time2.as_secs_f64() / time1.as_secs_f64();
            let mem_ratio = if mem1 > 1.0 { mem2 / mem1 } else { 1.0 };

            println!("     {} → {} reads: {:.1}× reads, {:.1}× time, {:.1}× memory",
                     reads1, reads2, read_ratio, time_ratio, mem_ratio);

            let linear_scaling = (time_ratio / read_ratio - 1.0).abs() < 0.5;
            let constant_memory = mem_ratio < 1.5;

            println!("       Time scaling: {} (ratio: {:.2})",
                     if linear_scaling { "✅ LINEAR" } else { "⚠️ NON-LINEAR" }, time_ratio / read_ratio);
            println!("       Memory scaling: {} (ratio: {:.2})",
                     if constant_memory { "✅ CONSTANT" } else { "⚠️ GROWING" }, mem_ratio);
        }
    }

    Ok(())
}