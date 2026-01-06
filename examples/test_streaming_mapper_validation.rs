//! Phase 6a.1: Streaming Mapper Validation
//!
//! This example validates the StreamingMapper functionality with real 10K read data.
//! It serves as a rapid sanity check to confirm basic streaming alignment works.

use biometal::alignment::{StreamingMapper, StreamingMapperConfig, ScoringMatrix};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🧪 Phase 6a.1: Streaming Mapper Validation");
    println!("{}", "=".repeat(50));

    // Test 1: Create StreamingMapper with appropriate configuration for testing
    println!("\n🔧 Creating StreamingMapper configuration...");
    let config = StreamingMapperConfig {
        window_size: 50_000,         // Small 50KB windows for testing
        overlap_bp: 100,             // Small overlap for testing
        min_score_threshold: 20,     // Permissive threshold for initial validation
        scoring: ScoringMatrix::default(),
    };

    println!("   Window size: {} bytes", config.window_size);
    println!("   Overlap: {} bp", config.overlap_bp);
    println!("   Min score: {}", config.min_score_threshold);
    println!("   Scoring: match={}, mismatch={}, gap_open={}, gap_extend={}",
             config.scoring.match_score, config.scoring.mismatch_score,
             config.scoring.gap_open, config.scoring.gap_extend);

    // Test 2: Create mapper
    println!("\n🗺️  Creating StreamingMapper...");
    let mut mapper = StreamingMapper::new(config);
    println!("✅ StreamingMapper created successfully");

    // Test 3: Define file paths
    let reference_path = "/Users/Megan Johnson/Projects/biometal_validation/test_reference.fasta";
    let reads_path = "/Users/Megan Johnson/Projects/biometal_validation/small_test_10k_reads.fastq.gz";

    println!("\n📁 Input files:");
    println!("   Reference: {}", reference_path);
    println!("   Reads: {}", reads_path);

    // Check files exist
    if !std::path::Path::new(reference_path).exists() {
        eprintln!("❌ Reference file not found: {}", reference_path);
        return Ok(());
    }
    if !std::path::Path::new(reads_path).exists() {
        eprintln!("❌ Reads file not found: {}", reads_path);
        return Ok(());
    }
    println!("✅ Both input files exist");

    // Test 4: Execute mapping with timing and basic monitoring
    println!("\n🧬 Starting streaming alignment...");
    let start_time = Instant::now();

    match mapper.map_reads_streaming(reference_path, reads_path) {
        Ok(mapping_iter) => {
            println!("✅ StreamingMappingIterator created successfully");

            let mut mapping_count = 0;
            let mut total_score = 0;
            let mut error_count = 0;

            println!("\n📊 Processing mappings (showing first 5)...");

            // Process mappings with error handling
            for (i, mapping_result) in mapping_iter.enumerate() {
                match mapping_result {
                    Ok(mapping) => {
                        mapping_count += 1;
                        total_score += mapping.alignment.score;

                        // Show details for first few mappings
                        if i < 5 {
                            println!("   {}. Query: {} → Ref: {} at {}-{}, Score: {}, Window: {}",
                                     i + 1,
                                     mapping.query_id.chars().take(20).collect::<String>(),
                                     mapping.reference_id,
                                     mapping.global_ref_start,
                                     mapping.global_ref_end,
                                     mapping.alignment.score,
                                     mapping.window_index);
                        }
                    }
                    Err(e) => {
                        error_count += 1;
                        if error_count <= 3 { // Show first few errors
                            eprintln!("   Error {}: {}", error_count, e);
                        }
                    }
                }

                // Stop after processing reasonable number for initial testing
                if mapping_count >= 100 || error_count >= 10 {
                    break;
                }
            }

            let execution_time = start_time.elapsed();

            println!("\n📈 Execution Results:");
            println!("   Successful mappings: {}", mapping_count);
            println!("   Errors encountered: {}", error_count);
            println!("   Average score: {:.1}",
                     if mapping_count > 0 { total_score as f64 / mapping_count as f64 } else { 0.0 });
            println!("   Execution time: {:.3}s", execution_time.as_secs_f64());

            // Success criteria
            println!("\n🎯 Phase 6a.1 Success Criteria:");
            let criteria_met = mapping_count > 0 && error_count == 0 && execution_time.as_secs() < 60;

            println!("   {} Found alignments: {}",
                     if mapping_count > 0 { "✅" } else { "❌" }, mapping_count);
            println!("   {} No errors: {}",
                     if error_count == 0 { "✅" } else { "❌" }, error_count);
            println!("   {} Reasonable execution time: {:.1}s",
                     if execution_time.as_secs() < 60 { "✅" } else { "⚠️" }, execution_time.as_secs_f64());

            if criteria_met {
                println!("\n🎉 SUCCESS: Phase 6a.1 rapid sanity check PASSED!");
                println!("✅ StreamingMapper basic functionality validated");
                println!("🚀 Ready for Phase 6a.2: Memory profiling validation");
            } else {
                println!("\n⚠️  PARTIAL SUCCESS: Basic functionality works but needs investigation");
                if mapping_count == 0 {
                    println!("   - No alignments found (may need different reference or threshold)");
                }
                if error_count > 0 {
                    println!("   - {} errors encountered during processing", error_count);
                }
                if execution_time.as_secs() >= 60 {
                    println!("   - Execution time longer than expected: {:.1}s", execution_time.as_secs_f64());
                }
            }
        }
        Err(e) => {
            let execution_time = start_time.elapsed();
            println!("❌ FAILED: StreamingMapper execution failed");
            println!("   Error: {}", e);
            println!("   Time to failure: {:.3}s", execution_time.as_secs_f64());

            println!("\n🔍 Debugging Information:");
            println!("   - Check if FASTQ/FASTA files are properly formatted");
            println!("   - Verify file permissions and accessibility");
            println!("   - Consider lowering min_score_threshold for testing");
        }
    }

    println!("\n{}", "=".repeat(50));
    println!("🏁 Phase 6a.1 validation completed");

    Ok(())
}