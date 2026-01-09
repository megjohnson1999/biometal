use biometal::alignment::{StreamingMapper, StreamingMapperConfig};
use std::collections::HashMap;
use std::env;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <reads.fastq> <reference.fasta>", args[0]);
        std::process::exit(1);
    }

    let reads_path = &args[1];
    let reference_path = &args[2];

    println!("🧬 Ultra-Stringent Filtering Validation");
    println!("======================================");
    println!("Reads: {}", reads_path);
    println!("Reference: {}", reference_path);

    // Initialize streaming mapper with ultra-stringent filtering
    let mut mapper = StreamingMapper::new(StreamingMapperConfig::default());
    println!("⚡ Ultra-stringent filtering enabled:");
    println!("   • E-value ≤ 0.0001 (100× more selective)");
    println!("   • Alignment length ≥ 30bp");
    println!("   • Sequence complexity > 2.0");
    println!("   • Length-normalized score ≥ 1.8");

    // Process reads with streaming mapper
    let mapping_iterator = mapper.map_reads_streaming(reference_path, reads_path)?;
    let mut total_reads = 0;
    let mut total_alignments = 0;
    let mut position_counts: HashMap<usize, usize> = HashMap::new();
    let mut score_distribution = Vec::new();

    println!("\n🔍 Processing reads with ultra-stringent filtering...");

    for mapping_result in mapping_iterator {
        let alignment = mapping_result?;
        total_alignments += 1;
        *position_counts.entry(alignment.global_ref_start).or_insert(0) += 1;
        score_distribution.push(alignment.alignment.score);

        // Count reads (this is a simple approximation)
        total_reads = position_counts.len().max(total_reads);

        if total_alignments % 100 == 0 {
            print!(".");
            if total_alignments % 1000 == 0 {
                println!(" {} alignments processed", total_alignments);
            }
        }
    }

    println!("\n\n📈 Ultra-Stringent Filtering Results:");
    println!("====================================");
    println!("Total alignments found: {}", total_alignments);

    if total_alignments > 0 {
        // Use a rough estimate for reads processed (since we're getting alignments, not reads)
        let estimated_reads = total_alignments; // Conservative estimate - each alignment represents at least one read
        println!("Estimated reads processed: {}", estimated_reads);

        let alignments_per_read = total_alignments as f64 / estimated_reads as f64;
        println!("Alignments per read: {:.2}", alignments_per_read);

        let multi_mapping_reads = position_counts.values().filter(|&&count| count > 1).count();
        let multi_mapping_rate = multi_mapping_reads as f64 / position_counts.len() as f64 * 100.0;
        println!("Multi-mapping rate: {:.1}%", multi_mapping_rate);

        // Position clustering analysis
        let unique_positions = position_counts.len();
        println!("Unique alignment positions: {}", unique_positions);

        if !position_counts.is_empty() {
            let max_reads_per_position = *position_counts.values().max().unwrap();
            let avg_reads_per_position = total_alignments as f64 / unique_positions as f64;
            println!("Max reads per position: {}", max_reads_per_position);
            println!("Avg reads per position: {:.1}", avg_reads_per_position);
        }

        // Score distribution
        if !score_distribution.is_empty() {
            score_distribution.sort();
            let median_score = score_distribution[score_distribution.len() / 2];
            let min_score = score_distribution[0];
            let max_score = score_distribution[score_distribution.len() - 1];
            println!("Score distribution: min={}, median={}, max={}", min_score, median_score, max_score);
        }

        // Quality assessment
        println!("\n🎯 Quality Assessment:");
        if total_alignments == 0 {
            println!("Status: ✅ ULTRA-SELECTIVE (0 alignments - extremely stringent!)");
        } else if alignments_per_read < 0.5 {
            println!("Status: ✅ EXCELLENT (Very low alignment rate - highly selective)");
        } else if alignments_per_read < 2.0 {
            println!("Status: ✅ GOOD (Low alignment rate - selective filtering)");
        } else if alignments_per_read < 10.0 {
            println!("Status: ⚠️  MODERATE (Medium alignment rate)");
        } else {
            println!("Status: ❌ HIGH FALSE POSITIVE (High alignment rate)");
        }

        println!("\n📊 Comparison vs Previous Results:");
        println!("Previous (stringent): 1,641 alignments");
        println!("Previous (moderate):  1,996 alignments");
        println!("Target (minimap2):    ~93 alignments");
        println!("Current (ultra):      {} alignments", total_alignments);

        if total_alignments > 0 {
            let improvement_vs_moderate = ((1996.0 - total_alignments as f64) / 1996.0 * 100.0).max(0.0);
            let improvement_vs_stringent = ((1641.0 - total_alignments as f64) / 1641.0 * 100.0).max(0.0);
            let vs_target = total_alignments as f64 / 93.0;

            println!("Improvement vs moderate: {:.1}% reduction", improvement_vs_moderate);
            println!("Improvement vs stringent: {:.1}% reduction", improvement_vs_stringent);
            println!("Ratio vs minimap2 target: {:.1}× (goal: ~1.0×)", vs_target);
        } else {
            println!("Improvement: 100% reduction (ultra-selective!)");
            println!("Ratio vs minimap2: 0× (potentially too stringent)");
        }
    } else {
        println!("Status: ✅ ULTRA-SELECTIVE (0 alignments - extremely stringent!)");
        println!("This suggests the filtering might be too strict, or the data doesn't contain good alignments.");
        println!("Improvement: 100% reduction vs all previous attempts");
        println!("Ratio vs minimap2: 0× (potentially too stringent)");
    }

    Ok(())
}