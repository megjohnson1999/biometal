use biometal::alignment::{smith_waterman, ScoringMatrix};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🎯 Alignment Accuracy Test");
    println!("==========================");
    println!("Testing whether our scoring system finds CORRECT alignments");

    let scoring = ScoringMatrix::default();

    // Test Case 1: Perfect match
    println!("\n📍 Test 1: Perfect Match");
    let reference_seq = b"GGTTCACTTGAGACACGAGCTCTGTACTGAATATACTCACAACGCTTGTTAATAAGCTCT";
    let query_seq = reference_seq; // Identical

    let alignment = smith_waterman(query_seq, reference_seq, &scoring);
    println!("Perfect match score: {}", alignment.score);
    println!("Expected score: {}", reference_seq.len() as i32 * scoring.match_score);

    // Test quality control
    let passes_quality = scoring.passes_quality_threshold(
        alignment.score,
        alignment.len(),
        query_seq.len(),
        reference_seq.len(),
        query_seq,
    );

    if passes_quality {
        println!("✅ Perfect match passes ultra-stringent filtering!");
    } else {
        println!("❌ Perfect match REJECTED by ultra-stringent filtering!");

        // Analyze why it failed
        let evalue = scoring.calculate_evalue(alignment.score, query_seq.len(), reference_seq.len());
        let complexity = scoring.sequence_complexity(query_seq);
        let normalized_score = scoring.length_normalized_score(alignment.score, alignment.len());

        println!("   E-value: {} (threshold: ≤ 0.0001)", evalue);
        println!("   Alignment length: {} (threshold: ≥ 30)", alignment.len());
        println!("   Sequence complexity: {:.2} (threshold: > 1.8)", complexity);
        println!("   Normalized score: {:.2} (threshold: ≥ 1.8)", normalized_score);
    }

    // Test Case 2: Near-perfect match (1 mismatch)
    println!("\n📍 Test 2: Near-Perfect Match (1 mismatch)");
    let mut near_perfect = reference_seq.to_vec();
    near_perfect[20] = b'A'; // Change one base from original

    let alignment2 = smith_waterman(&near_perfect, reference_seq, &scoring);
    println!("Near-perfect score: {}", alignment2.score);

    let passes_quality2 = scoring.passes_quality_threshold(
        alignment2.score,
        alignment2.len(),
        near_perfect.len(),
        reference_seq.len(),
        &near_perfect,
    );

    if passes_quality2 {
        println!("✅ Near-perfect match passes filtering");
    } else {
        println!("⚠️ Near-perfect match rejected (may be acceptable if very stringent)");
    }

    // Test Case 3: Low complexity sequence (should be filtered)
    println!("\n📍 Test 3: Low Complexity Sequence (should be filtered)");
    let low_complexity = vec![b'A'; 60]; // All A's

    let alignment3 = smith_waterman(&low_complexity, reference_seq, &scoring);
    let passes_quality3 = scoring.passes_quality_threshold(
        alignment3.score,
        alignment3.len(),
        low_complexity.len(),
        reference_seq.len(),
        &low_complexity,
    );

    if passes_quality3 {
        println!("❌ Low complexity sequence passed filtering - should be rejected!");
    } else {
        println!("✅ Low complexity correctly filtered out");
        let complexity = scoring.sequence_complexity(&low_complexity);
        println!("   Complexity: {:.2} (threshold: > 1.8)", complexity);
    }

    // Test Case 4: Short read (should be filtered)
    println!("\n📍 Test 4: Short Read (<30bp, should be filtered)");
    let short_read = b"GGTTCACTTGAGACACGAGCTCTGT"; // 25bp

    let alignment4 = smith_waterman(short_read, reference_seq, &scoring);
    let passes_quality4 = scoring.passes_quality_threshold(
        alignment4.score,
        alignment4.len(),
        short_read.len(),
        reference_seq.len(),
        short_read,
    );

    if passes_quality4 {
        println!("❌ Short read passed filtering - should be rejected!");
    } else {
        println!("✅ Short read correctly filtered out");
        println!("   Alignment length: {} (threshold: ≥ 30)", alignment4.len());
    }

    // Test Case 5: High complexity, balanced sequence
    println!("\n📍 Test 5: High Complexity Sequence");
    let complex_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 60bp balanced

    let alignment5 = smith_waterman(complex_seq, complex_seq, &scoring); // Perfect match to itself
    let passes_quality5 = scoring.passes_quality_threshold(
        alignment5.score,
        alignment5.len(),
        complex_seq.len(),
        complex_seq.len(),
        complex_seq,
    );

    if passes_quality5 {
        println!("✅ High complexity sequence passes filtering");
        let complexity = scoring.sequence_complexity(complex_seq);
        println!("   Complexity: {:.2} (threshold: > 1.8)", complexity);
    } else {
        println!("❌ High complexity sequence rejected - may indicate over-filtering");
    }

    // Test Case 6: Test against wrong reference (should have low score)
    println!("\n📍 Test 6: Wrong Reference Alignment");
    let wrong_ref = vec![b'T'; 60]; // All T's

    let alignment6 = smith_waterman(reference_seq, &wrong_ref, &scoring);
    let passes_quality6 = scoring.passes_quality_threshold(
        alignment6.score,
        alignment6.len(),
        reference_seq.len(),
        wrong_ref.len(),
        reference_seq,
    );

    if passes_quality6 {
        println!("❌ Wrong reference alignment passed - false positive!");
        println!("   Score: {}", alignment6.score);
    } else {
        println!("✅ Wrong reference correctly rejected");
        println!("   Score: {} (low as expected)", alignment6.score);
    }

    println!("\n🏁 Ultra-Stringent Filtering Assessment:");
    println!("========================================");
    println!("✅ = Working correctly");
    println!("❌ = Problem detected");
    println!("⚠️ = Potentially over-strict (review if too many good alignments rejected)");
    println!("");
    println!("If perfect matches are being rejected, the thresholds may be too stringent.");
    println!("Consider adjusting E-value, complexity, or normalized score thresholds.");

    Ok(())
}