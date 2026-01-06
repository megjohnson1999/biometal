//! Alignment Quality Validation
//!
//! Assess whether our alignments are legitimate or false positives by examining:
//! - Alignment length and coverage
//! - Score distributions
//! - Mapping positions (are multiple reads piling up unrealistically?)
//! - Sequence similarity at reported positions

use biometal::alignment::{StreamingMapper, StreamingMapperConfig, ScoringMatrix};
use std::collections::HashMap;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔍 Alignment Quality Validation");
    println!("{}", "=".repeat(60));

    let fastq_path = "/Users/Megan Johnson/Projects/biometal_validation/small_test_10k_reads.fastq.gz";
    let reference_path = "/Users/Megan Johnson/Projects/biometal_validation/test_reference.fasta";

    if !Path::new(fastq_path).exists() || !Path::new(reference_path).exists() {
        println!("❌ Required files not found");
        return Ok(());
    }

    analyze_alignment_quality(reference_path, fastq_path)?;

    Ok(())
}

fn analyze_alignment_quality(reference_path: &str, fastq_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    println!("\n📊 Detailed Alignment Quality Analysis");

    let config = StreamingMapperConfig {
        window_size: 50_000,
        overlap_bp: 1_000,
        min_score_threshold: 15,
        scoring: ScoringMatrix::default(),
    };

    let mut mapper = StreamingMapper::new(config);

    match mapper.map_reads_streaming(reference_path, fastq_path) {
        Ok(mapping_iter) => {
            let mut alignments = Vec::new();
            let mut position_counts: HashMap<(String, usize), usize> = HashMap::new();
            let mut score_distribution = Vec::new();
            let mut length_distribution = Vec::new();
            let mut query_counts: HashMap<String, usize> = HashMap::new();

            println!("   Collecting alignment statistics...");

            for (i, mapping_result) in mapping_iter.enumerate() {
                match mapping_result {
                    Ok(mapping) => {
                        // Track mapping position (potential pileup indicator)
                        let pos_key = (mapping.reference_id.clone(), mapping.global_ref_start / 100 * 100); // 100bp bins
                        *position_counts.entry(pos_key).or_insert(0) += 1;

                        // Score distribution
                        score_distribution.push(mapping.alignment.score);

                        // Alignment length
                        let alignment_length = mapping.global_ref_end - mapping.global_ref_start;
                        length_distribution.push(alignment_length);

                        // Query mapping frequency
                        *query_counts.entry(mapping.query_id.clone()).or_insert(0) += 1;

                        alignments.push(mapping);

                        // Stop after collecting sufficient data
                        if alignments.len() >= 2000 {
                            break;
                        }
                    }
                    Err(_) => continue,
                }
            }

            analyze_statistics(&alignments, &position_counts, &score_distribution, &length_distribution, &query_counts);
        }
        Err(e) => {
            println!("   ❌ Analysis failed: {}", e);
        }
    }

    Ok(())
}

fn analyze_statistics(
    alignments: &[biometal::alignment::MappingResult],
    position_counts: &HashMap<(String, usize), usize>,
    scores: &[i32],
    lengths: &[usize],
    query_counts: &HashMap<String, usize>,
) {
    println!("\n📈 Quality Statistics:");
    println!("   Total alignments analyzed: {}", alignments.len());

    // Score analysis
    if !scores.is_empty() {
        let min_score = *scores.iter().min().unwrap();
        let max_score = *scores.iter().max().unwrap();
        let avg_score = scores.iter().sum::<i32>() as f64 / scores.len() as f64;
        let median_score = {
            let mut sorted_scores = scores.to_vec();
            sorted_scores.sort();
            sorted_scores[sorted_scores.len() / 2]
        };

        println!("\n🎯 Score Distribution:");
        println!("     Min score: {}", min_score);
        println!("     Max score: {}", max_score);
        println!("     Average: {:.1}", avg_score);
        println!("     Median: {}", median_score);

        // Score quality assessment
        let low_scores = scores.iter().filter(|&&s| s < 30).count();
        let high_scores = scores.iter().filter(|&&s| s >= 100).count();

        println!("     Low scores (<30): {} ({:.1}%)", low_scores,
                 low_scores as f64 / scores.len() as f64 * 100.0);
        println!("     High scores (≥100): {} ({:.1}%)", high_scores,
                 high_scores as f64 / scores.len() as f64 * 100.0);

        let quality_assessment = if avg_score >= 50.0 && low_scores < scores.len() / 2 {
            "✅ REASONABLE"
        } else if avg_score >= 25.0 {
            "⚠️ MARGINAL"
        } else {
            "❌ POOR"
        };
        println!("     Quality: {}", quality_assessment);
    }

    // Length analysis
    if !lengths.is_empty() {
        let min_length = *lengths.iter().min().unwrap();
        let max_length = *lengths.iter().max().unwrap();
        let avg_length = lengths.iter().sum::<usize>() as f64 / lengths.len() as f64;

        println!("\n📏 Alignment Length Distribution:");
        println!("     Min length: {} bp", min_length);
        println!("     Max length: {} bp", max_length);
        println!("     Average: {:.1} bp", avg_length);

        let short_alignments = lengths.iter().filter(|&&l| l < 20).count();
        let long_alignments = lengths.iter().filter(|&&l| l >= 50).count();

        println!("     Very short (<20bp): {} ({:.1}%)", short_alignments,
                 short_alignments as f64 / lengths.len() as f64 * 100.0);
        println!("     Long (≥50bp): {} ({:.1}%)", long_alignments,
                 long_alignments as f64 / lengths.len() as f64 * 100.0);

        let length_quality = if avg_length >= 30.0 && short_alignments < lengths.len() / 3 {
            "✅ REASONABLE"
        } else if avg_length >= 15.0 {
            "⚠️ SHORT"
        } else {
            "❌ TOO SHORT"
        };
        println!("     Quality: {}", length_quality);
    }

    // Position pileup analysis (false positive indicator)
    println!("\n📍 Position Distribution Analysis:");
    let mut pileup_counts = Vec::new();
    for count in position_counts.values() {
        pileup_counts.push(*count);
    }
    pileup_counts.sort();

    if !pileup_counts.is_empty() {
        let max_pileup = *pileup_counts.iter().max().unwrap();
        let avg_pileup = pileup_counts.iter().sum::<usize>() as f64 / pileup_counts.len() as f64;
        let high_pileup_positions = pileup_counts.iter().filter(|&&c| c >= 20).count();

        println!("     Unique positions: {}", position_counts.len());
        println!("     Max reads/position: {}", max_pileup);
        println!("     Avg reads/position: {:.1}", avg_pileup);
        println!("     High pileup positions (≥20): {} ({:.1}%)",
                 high_pileup_positions,
                 high_pileup_positions as f64 / position_counts.len() as f64 * 100.0);

        let pileup_assessment = if max_pileup > 100 || high_pileup_positions > position_counts.len() / 10 {
            "❌ SUSPICIOUS PILEUP"
        } else if max_pileup > 50 {
            "⚠️ MODERATE PILEUP"
        } else {
            "✅ REASONABLE"
        };
        println!("     Assessment: {}", pileup_assessment);
    }

    // Query mapping frequency (multi-mapping indicator)
    println!("\n🔄 Query Mapping Analysis:");
    let multi_mapping_queries = query_counts.values().filter(|&&c| c > 1).count();
    let total_queries = query_counts.len();

    if let Some(max_mappings) = query_counts.values().max() {
        println!("     Total unique queries: {}", total_queries);
        println!("     Multi-mapping queries: {} ({:.1}%)",
                 multi_mapping_queries,
                 multi_mapping_queries as f64 / total_queries as f64 * 100.0);
        println!("     Max mappings/query: {}", max_mappings);

        let multi_mapping_assessment = if multi_mapping_queries as f64 / total_queries as f64 > 0.8 {
            "❌ EXCESSIVE MULTI-MAPPING"
        } else if multi_mapping_queries as f64 / total_queries as f64 > 0.3 {
            "⚠️ HIGH MULTI-MAPPING"
        } else {
            "✅ REASONABLE"
        };
        println!("     Assessment: {}", multi_mapping_assessment);
    }

    // Overall false positive indicators
    println!("\n🚨 False Positive Risk Assessment:");

    let score_risk = if !scores.is_empty() {
        let avg_score = scores.iter().sum::<i32>() as f64 / scores.len() as f64;
        let low_score_pct = scores.iter().filter(|&&s| s < 30).count() as f64 / scores.len() as f64;
        avg_score < 30.0 || low_score_pct > 0.5
    } else { false };

    let length_risk = if !lengths.is_empty() {
        let avg_length = lengths.iter().sum::<usize>() as f64 / lengths.len() as f64;
        let short_pct = lengths.iter().filter(|&&l| l < 20).count() as f64 / lengths.len() as f64;
        avg_length < 20.0 || short_pct > 0.4
    } else { false };

    let pileup_risk = if !pileup_counts.is_empty() {
        let max_pileup = *pileup_counts.iter().max().unwrap();
        max_pileup > 100
    } else { false };

    let multi_mapping_risk = if !query_counts.is_empty() {
        let multi_pct = query_counts.values().filter(|&&c| c > 1).count() as f64 / query_counts.len() as f64;
        multi_pct > 0.7
    } else { false };

    let risk_factors = [score_risk, length_risk, pileup_risk, multi_mapping_risk].iter().filter(|&&r| r).count();

    println!("     Score quality risk: {}", if score_risk { "⚠️ HIGH" } else { "✅ LOW" });
    println!("     Length quality risk: {}", if length_risk { "⚠️ HIGH" } else { "✅ LOW" });
    println!("     Position pileup risk: {}", if pileup_risk { "⚠️ HIGH" } else { "✅ LOW" });
    println!("     Multi-mapping risk: {}", if multi_mapping_risk { "⚠️ HIGH" } else { "✅ LOW" });

    let overall_assessment = match risk_factors {
        0 => "✅ LOW FALSE POSITIVE RISK",
        1 => "⚠️ MODERATE FALSE POSITIVE RISK",
        2 => "⚠️ HIGH FALSE POSITIVE RISK",
        _ => "❌ VERY HIGH FALSE POSITIVE RISK - LIKELY SPURIOUS ALIGNMENTS"
    };

    println!("\n🏆 Overall Assessment: {}", overall_assessment);

    if risk_factors > 0 {
        println!("\n💡 Recommendations:");
        if score_risk {
            println!("   • Increase score threshold (try ≥30)");
        }
        if length_risk {
            println!("   • Filter short alignments (<20bp)");
        }
        if pileup_risk {
            println!("   • Investigate reference sequence quality");
        }
        if multi_mapping_risk {
            println!("   • Add proper multi-mapping filtering");
        }
        println!("   • Compare with ground truth alignments");
        println!("   • Validate against known sequence relationships");
    }
}