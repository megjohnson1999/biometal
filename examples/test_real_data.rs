//! Step 3: Real Data Testing - 10K FASTQ reads
//!
//! Test the streaming mapper with real sequencing data to validate
//! functionality and measure performance on authentic datasets.

use biometal::alignment::{StreamingMapper, StreamingMapperConfig, ScoringMatrix};
use std::process::Command;
use std::time::Instant;
use std::path::Path;

fn get_memory_usage_mb() -> Result<f64, Box<dyn std::error::Error>> {
    let pid = std::process::id();

    let output = Command::new("ps")
        .args(&["-o", "rss=", "-p", &pid.to_string()])
        .output()?;

    if output.status.success() {
        let rss_kb: f64 = String::from_utf8(output.stdout)?
            .trim()
            .parse()
            .unwrap_or(0.0);
        Ok(rss_kb / 1024.0)
    } else {
        Ok(0.0)
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🧬 Step 3: Real Data Testing - 10K FASTQ reads");
    println!("{}", "=".repeat(60));

    // Check if files exist
    let fastq_path = "/Users/Megan Johnson/Projects/biometal_validation/small_test_10k_reads.fastq.gz";
    let reference_path = "/Users/Megan Johnson/Projects/biometal_validation/test_reference.fasta";

    if !Path::new(fastq_path).exists() {
        println!("❌ FASTQ file not found: {}", fastq_path);
        return Ok(());
    }

    if !Path::new(reference_path).exists() {
        println!("❌ Reference file not found: {}", reference_path);
        return Ok(());
    }

    test_real_data_functionality(reference_path, fastq_path)?;
    test_real_data_performance(reference_path, fastq_path)?;
    test_memory_behavior_real_data(reference_path, fastq_path)?;

    println!("\n🎯 Step 3 Complete: Real data validation complete");
    Ok(())
}

fn test_real_data_functionality(reference_path: &str, fastq_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    println!("\n🔬 Test 1: Basic Functionality with Real Data");

    // Use permissive settings to ensure we find some mappings
    let config = StreamingMapperConfig {
        window_size: 50_000,   // 50KB windows for good coverage
        overlap_bp: 1_000,     // 1KB overlap
        min_score_threshold: 15, // Permissive threshold
        scoring: ScoringMatrix::default(),
    };

    let mut mapper = StreamingMapper::new(config);

    println!("   Dataset: 10K real FASTQ reads");
    println!("   Reference: Real bacterial sequences (2 chromosomes)");
    println!("   Settings: 50KB windows, 1KB overlap, score ≥ 15");

    let initial_memory = get_memory_usage_mb()?;
    println!("   Initial memory: {:.1} MB", initial_memory);

    match mapper.map_reads_streaming(reference_path, fastq_path) {
        Ok(mapping_iter) => {
            let mut total_mappings = 0;
            let mut score_sum = 0;
            let mut min_score = i32::MAX;
            let mut max_score = i32::MIN;
            let mut sample_mappings = Vec::new();

            println!("   Processing mappings...");

            for (i, mapping_result) in mapping_iter.enumerate() {
                match mapping_result {
                    Ok(mapping) => {
                        total_mappings += 1;
                        score_sum += mapping.alignment.score;
                        min_score = min_score.min(mapping.alignment.score);
                        max_score = max_score.max(mapping.alignment.score);

                        // Collect first 5 mappings as samples
                        if sample_mappings.len() < 5 {
                            sample_mappings.push(mapping);
                        }

                        // Progress indicator
                        if i % 1000 == 0 && i > 0 {
                            println!("     Processed {} mappings...", i);
                        }

                        // Limit processing for initial test
                        if total_mappings >= 5000 {
                            break;
                        }
                    }
                    Err(e) => {
                        println!("   ❌ Mapping error: {:?}", e);
                        break;
                    }
                }
            }

            let final_memory = get_memory_usage_mb()?;
            let memory_delta = final_memory - initial_memory;

            println!("\n   Results:");
            println!("     Total mappings: {}", total_mappings);
            println!("     Memory used: {:.1} MB", memory_delta);

            if total_mappings > 0 {
                let avg_score = score_sum as f64 / total_mappings as f64;
                println!("     Score range: {} - {} (avg: {:.1})", min_score, max_score, avg_score);

                println!("\n   Sample mappings:");
                for (i, mapping) in sample_mappings.iter().enumerate() {
                    println!("     {}: {} at {}-{}, score: {}",
                             i + 1,
                             mapping.query_id.chars().take(20).collect::<String>(),
                             mapping.global_ref_start,
                             mapping.global_ref_end,
                             mapping.alignment.score);
                }

                // Basic sanity checks
                let reasonable_mappings = total_mappings >= 100;
                let reasonable_scores = min_score > 0;
                let constant_memory = memory_delta < 50.0;

                println!("\n   Validation:");
                println!("     Sufficient mappings: {} ({} ≥ 100)",
                         if reasonable_mappings { "✅ YES" } else { "⚠️ NO" },
                         total_mappings);
                println!("     Positive scores: {} ({} > 0)",
                         if reasonable_scores { "✅ YES" } else { "⚠️ NO" },
                         min_score);
                println!("     Memory efficient: {} ({:.1} MB < 50MB)",
                         if constant_memory { "✅ YES" } else { "⚠️ NO" },
                         memory_delta);

            } else {
                println!("     ❌ No mappings found - possible issues:");
                println!("       - Threshold too high (try lower score threshold)");
                println!("       - Reference/query mismatch");
                println!("       - Window size issues");
            }
        }
        Err(e) => {
            println!("   ❌ Failed to create mapping iterator: {}", e);
            return Err(e.into());
        }
    }

    Ok(())
}

fn test_real_data_performance(reference_path: &str, fastq_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    println!("\n⚡ Test 2: Performance with Real Data");

    let config = StreamingMapperConfig {
        window_size: 100_000,  // 100KB windows
        overlap_bp: 2_000,     // 2KB overlap
        min_score_threshold: 20,
        scoring: ScoringMatrix::default(),
    };

    let mut mapper = StreamingMapper::new(config);

    println!("   Testing performance with 1000 reads (subset)");

    let start_time = Instant::now();

    match mapper.map_reads_streaming(reference_path, fastq_path) {
        Ok(mapping_iter) => {
            let mut processed_reads = 0;
            let mut mappings_found = 0;

            for mapping_result in mapping_iter {
                match mapping_result {
                    Ok(_mapping) => {
                        mappings_found += 1;
                        processed_reads += 1;

                        // Stop after processing 1000 reads worth of mappings
                        if processed_reads >= 1000 {
                            break;
                        }
                    }
                    Err(_) => {
                        // Count read attempts even if they fail
                        processed_reads += 1;
                        if processed_reads >= 1000 {
                            break;
                        }
                    }
                }
            }

            let elapsed = start_time.elapsed();
            let reads_per_sec = processed_reads as f64 / elapsed.as_secs_f64();

            println!("   Results:");
            println!("     Processed reads: {}", processed_reads);
            println!("     Mappings found: {}", mappings_found);
            println!("     Time: {:.3}s", elapsed.as_secs_f64());
            println!("     Rate: {:.0} reads/sec", reads_per_sec);

            let acceptable_speed = reads_per_sec >= 50.0;
            println!("     Speed assessment: {} ({:.0} reads/sec {} 50 reads/sec)",
                     if acceptable_speed { "✅ ACCEPTABLE" } else { "⚠️ SLOW" },
                     reads_per_sec,
                     if acceptable_speed { "≥" } else { "<" });
        }
        Err(e) => {
            println!("   ❌ Performance test failed: {}", e);
        }
    }

    Ok(())
}

fn test_memory_behavior_real_data(reference_path: &str, fastq_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    println!("\n💾 Test 3: Memory Behavior with Real Data");

    let config = StreamingMapperConfig {
        window_size: 75_000,   // 75KB windows
        overlap_bp: 1_500,     // 1.5KB overlap
        min_score_threshold: 18,
        scoring: ScoringMatrix::default(),
    };

    let mut mapper = StreamingMapper::new(config);

    let initial_memory = get_memory_usage_mb()?;
    println!("   Initial memory: {:.1} MB", initial_memory);

    match mapper.map_reads_streaming(reference_path, fastq_path) {
        Ok(mapping_iter) => {
            let mut checkpoints = vec![];
            let mut mapping_count = 0;

            for (i, mapping_result) in mapping_iter.enumerate() {
                if mapping_result.is_ok() {
                    mapping_count += 1;
                }

                // Memory checkpoint every 500 mappings processed
                if i % 500 == 0 && i > 0 {
                    let current_memory = get_memory_usage_mb().unwrap_or(0.0);
                    let delta = current_memory - initial_memory;
                    checkpoints.push((i, delta));

                    println!("     Checkpoint {}: {:.1} MB (+{:.1} MB)", i, current_memory, delta);
                }

                // Stop after reasonable sample
                if i >= 2500 {
                    break;
                }
            }

            let final_memory = get_memory_usage_mb()?;
            let final_delta = final_memory - initial_memory;

            println!("\n   Memory Analysis:");
            println!("     Final memory: {:.1} MB (+{:.1} MB)", final_memory, final_delta);
            println!("     Mappings processed: {}", mapping_count);

            if !checkpoints.is_empty() {
                let memory_growth = checkpoints.last().unwrap().1 - checkpoints.first().unwrap().1;
                println!("     Memory growth during processing: {:.1} MB", memory_growth);

                let constant_memory = memory_growth < 10.0;
                println!("     Memory stability: {} ({:.1} MB growth < 10MB)",
                         if constant_memory { "✅ STABLE" } else { "⚠️ GROWING" },
                         memory_growth);
            }

            let total_reasonable = final_delta < 100.0;
            println!("     Overall efficiency: {} ({:.1} MB < 100MB)",
                     if total_reasonable { "✅ EFFICIENT" } else { "⚠️ HIGH" },
                     final_delta);
        }
        Err(e) => {
            println!("   ❌ Memory test failed: {}", e);
        }
    }

    Ok(())
}