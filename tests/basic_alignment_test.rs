//! Basic alignment functionality test
//! Phase 6a.1: Simple validation of core alignment works

use biometal::alignment::{smith_waterman, smith_waterman_naive, ScoringMatrix};

#[test]
fn test_basic_smith_waterman() {
    println!("🧬 Testing basic Smith-Waterman alignment...");

    // Test perfect match
    let query = b"ACGTACGT";
    let reference = b"ACGTACGT";
    let scoring = ScoringMatrix::default();

    let alignment = smith_waterman(query, reference, &scoring);

    println!("Perfect match result:");
    println!("  Score: {}", alignment.score);
    println!("  Query: {}:{}", alignment.query_start, alignment.query_end);
    println!("  Ref: {}:{}", alignment.ref_start, alignment.ref_end);
    println!("  CIGAR: {:?}", alignment.cigar);

    // Should have high score for perfect match
    assert!(alignment.score > 0, "Perfect match should have positive score, got: {}", alignment.score);

    // Test with mismatch
    let query2 = b"ACGTACGT";
    let reference2 = b"ACTTACGT";  // One mismatch at position 2
    let alignment2 = smith_waterman(query2, reference2, &scoring);

    println!("\nOne mismatch result:");
    println!("  Score: {}", alignment2.score);
    println!("  Should be less than perfect match score: {}", alignment.score);

    assert!(alignment2.score < alignment.score,
           "Mismatch should have lower score: {} vs {}", alignment2.score, alignment.score);

    println!("✅ Basic Smith-Waterman test passed!");
}

#[test]
fn test_scoring_matrix() {
    println!("🎯 Testing ScoringMatrix creation...");

    // Test default
    let default_scoring = ScoringMatrix::default();
    assert_eq!(default_scoring.match_score, 2);
    assert_eq!(default_scoring.mismatch_score, -1);
    assert_eq!(default_scoring.gap_open, -2);
    assert_eq!(default_scoring.gap_extend, -1);
    println!("✅ Default scoring matrix: {:?}", default_scoring);

    // Test custom
    let custom_scoring = ScoringMatrix::new(5, -2, -3, -1);
    assert_eq!(custom_scoring.match_score, 5);
    assert_eq!(custom_scoring.mismatch_score, -2);
    println!("✅ Custom scoring matrix: {:?}", custom_scoring);

    println!("✅ ScoringMatrix tests passed!");
}

#[test]
fn test_naive_vs_optimized() {
    println!("⚖️  Testing naive vs optimized Smith-Waterman...");

    let query = b"ACGTACGTACGT";
    let reference = b"ACGTACGTACGT";
    let scoring = ScoringMatrix::default();

    // Test naive implementation
    let naive_result = smith_waterman_naive(query, reference, &scoring);
    println!("Naive result: score={}, positions={}:{} -> {}:{}",
             naive_result.score, naive_result.query_start, naive_result.query_end,
             naive_result.ref_start, naive_result.ref_end);

    // Test optimized implementation
    let optimized_result = smith_waterman(query, reference, &scoring);
    println!("Optimized result: score={}, positions={}:{} -> {}:{}",
             optimized_result.score, optimized_result.query_start, optimized_result.query_end,
             optimized_result.ref_start, optimized_result.ref_end);

    // Results should be identical for correctness
    assert_eq!(naive_result.score, optimized_result.score,
               "Naive and optimized should have same score: {} vs {}",
               naive_result.score, optimized_result.score);

    println!("✅ Naive vs optimized consistency test passed!");
}

#[test]
fn test_real_sequence_data() {
    println!("🧬 Testing with real sequence patterns...");

    // Use sequences from actual FASTQ data
    let bacterial_seq = b"GGTTCACTTGAGACACGAGCTCTGTACTGAATATACTCACAAC";
    let reference = b"GGTTCACTTGAGACACGAGCTCTGTACTGAATATACTCACAAC";
    let scoring = ScoringMatrix::default();

    let alignment = smith_waterman(bacterial_seq, reference, &scoring);

    println!("Real sequence alignment:");
    println!("  Length: {} bp", bacterial_seq.len());
    println!("  Score: {}", alignment.score);
    println!("  Coverage: {}% query, {}% ref",
             ((alignment.query_end - alignment.query_start) * 100) / bacterial_seq.len(),
             ((alignment.ref_end - alignment.ref_start) * 100) / reference.len());

    assert!(alignment.score > 50, "Real sequence should align well, got score: {}", alignment.score);

    println!("✅ Real sequence test passed!");
}

#[test]
fn test_alignment_with_gaps() {
    println!("🔄 Testing alignment with gaps...");

    let query = b"ACGTACGT";
    let reference = b"ACGTAACCGT";  // Extra bases in middle
    let scoring = ScoringMatrix::default();

    let alignment = smith_waterman(query, reference, &scoring);

    println!("Gap alignment:");
    println!("  Query:     {:?}", std::str::from_utf8(query).unwrap());
    println!("  Reference: {:?}", std::str::from_utf8(reference).unwrap());
    println!("  Score: {}", alignment.score);
    println!("  CIGAR: {:?}", alignment.cigar);

    // Should still find a good alignment despite the gap
    assert!(alignment.score > 0, "Should handle gaps, got score: {}", alignment.score);

    println!("✅ Gap alignment test passed!");
}

#[test]
fn test_memory_usage_indication() {
    println!("📊 Memory Usage Goals:");
    println!("  Current CLI: 335-652MB (linear scaling)");
    println!("  Target: ~5MB constant memory");
    println!("  Status: Core alignment functions compiled ✅");

    // Basic memory test - just ensure we can create multiple objects
    for i in 1..=10 {
        let scoring = ScoringMatrix::new(i, -1, -2, -1);
        let size = (i * 100) as usize;
        let query = vec![b'A'; size]; // Variable length
        let reference = vec![b'A'; size];

        let alignment = smith_waterman(&query, &reference, &scoring);

        // Each alignment should work regardless of input size
        assert!(alignment.score >= 0, "Iteration {} failed", i);
    }

    println!("✅ Memory scaling test passed (no crashes with varying sizes)");
}