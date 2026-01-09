//! Scoring matrices for sequence alignment

/// Scoring matrix for sequence alignment
///
/// Defines scores for matches, mismatches, and gap penalties used in
/// Smith-Waterman and other alignment algorithms.
///
/// # Example
///
/// ```
/// use biometal::alignment::ScoringMatrix;
///
/// // Default scoring (match=2, mismatch=-1, gap_open=-2, gap_extend=-1)
/// let scoring = ScoringMatrix::default();
///
/// // Custom scoring
/// let custom = ScoringMatrix {
///     match_score: 5,
///     mismatch_score: -4,
///     gap_open: -10,
///     gap_extend: -1,
/// };
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScoringMatrix {
    /// Score for matching bases (positive)
    pub match_score: i32,
    /// Score for mismatching bases (negative)
    pub mismatch_score: i32,
    /// Penalty for opening a gap (negative)
    pub gap_open: i32,
    /// Penalty for extending a gap (negative)
    pub gap_extend: i32,
}

impl Default for ScoringMatrix {
    /// Default scoring parameters
    ///
    /// - Match: +2
    /// - Mismatch: -1
    /// - Gap open: -2
    /// - Gap extend: -1
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_score: -1,
            gap_open: -2,
            gap_extend: -1,
        }
    }
}

impl ScoringMatrix {
    /// Create a new scoring matrix
    pub fn new(match_score: i32, mismatch_score: i32, gap_open: i32, gap_extend: i32) -> Self {
        Self {
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
        }
    }

    /// Calculate the score for aligning two bases
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::alignment::ScoringMatrix;
    ///
    /// let scoring = ScoringMatrix::default();
    /// assert_eq!(scoring.score(b'A', b'A'), 2);  // Match
    /// assert_eq!(scoring.score(b'A', b'C'), -1); // Mismatch
    /// ```
    pub fn score(&self, a: u8, b: u8) -> i32 {
        if a == b {
            self.match_score
        } else {
            self.mismatch_score
        }
    }

    /// Calculate E-value for statistical significance of alignment
    ///
    /// Uses Karlin-Altschul statistics to compute the expected number of alignments
    /// with at least the given score occurring by chance in a database search.
    /// Lower E-values indicate more significant alignments.
    ///
    /// # Parameters
    /// - `score`: Raw alignment score
    /// - `query_len`: Length of query sequence
    /// - `ref_len`: Length of reference sequence (or database size)
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::alignment::ScoringMatrix;
    ///
    /// let scoring = ScoringMatrix::default();
    /// let evalue = scoring.calculate_evalue(50, 100, 1000000);
    /// assert!(evalue < 1.0); // Good alignment should have E-value < 1
    /// ```
    pub fn calculate_evalue(&self, score: i32, query_len: usize, ref_len: usize) -> f64 {
        // Karlin-Altschul parameters for nucleotide sequences with match=+2, mismatch=-1
        // Adjusted parameters based on empirical studies for biometal's scoring scheme
        const LAMBDA: f64 = 0.825; // Information content parameter (higher for +2/-1 scoring)
        const K: f64 = 0.046;      // Search space parameter

        // Apply effective length corrections (standard in BLAST)
        // Account for edge effects in short sequences
        let effective_query_len = (query_len as f64).max(10.0) - 8.0;
        let effective_ref_len = (ref_len as f64).max(100.0) - 40.0;

        let search_space = effective_query_len * effective_ref_len;

        // E-value = K * m * n * exp(-lambda * S)
        let evalue = K * search_space * (-LAMBDA * score as f64).exp();
        evalue.max(1e-100) // Floor to prevent underflow
    }

    /// Calculate length-normalized alignment score
    ///
    /// Normalizes raw scores by alignment length to prevent bias toward longer
    /// alignments. This helps distinguish high-quality short alignments from
    /// low-quality long alignments.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::alignment::ScoringMatrix;
    ///
    /// let scoring = ScoringMatrix::default();
    /// let normalized = scoring.length_normalized_score(20, 10); // Score 20, length 10
    /// assert_eq!(normalized, 2.0); // 20/10 = 2.0 per base
    /// ```
    pub fn length_normalized_score(&self, score: i32, alignment_len: usize) -> f64 {
        if alignment_len == 0 {
            return 0.0;
        }
        score as f64 / alignment_len as f64
    }

    /// Calculate sequence complexity score
    ///
    /// Computes the Shannon entropy of a sequence to identify low-complexity
    /// regions (e.g., homopolymers, repeats) that may cause false positive
    /// alignments. Returns a value between 0 (no complexity) and 2 (maximum complexity).
    ///
    /// Sequences with complexity < 1.0 are typically repetitive and should be
    /// filtered to reduce false positives.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::alignment::ScoringMatrix;
    ///
    /// let scoring = ScoringMatrix::default();
    /// let complex_seq = b"ACGTACGTACGT"; // High complexity
    /// let simple_seq = b"AAAAAAAAAAAA"; // Low complexity (homopolymer)
    ///
    /// assert!(scoring.sequence_complexity(complex_seq) > 1.5);
    /// assert!(scoring.sequence_complexity(simple_seq) < 0.1);
    /// ```
    pub fn sequence_complexity(&self, sequence: &[u8]) -> f64 {
        if sequence.is_empty() {
            return 0.0;
        }

        // Count nucleotide frequencies
        let mut counts = [0u32; 256];
        for &base in sequence {
            counts[base as usize] += 1;
        }

        // Calculate Shannon entropy: H = -Σ(p_i * log2(p_i))
        let seq_len = sequence.len() as f64;
        let mut entropy = 0.0;

        for count in &counts {
            if *count > 0 {
                let probability = *count as f64 / seq_len;
                entropy -= probability * probability.log2();
            }
        }

        entropy
    }

    /// Check if an alignment passes stringent quality thresholds
    ///
    /// Combines E-value significance, length requirements, and complexity filtering
    /// to provide a single quality assessment. This replaces the simple score
    /// threshold with statistical rigor tuned to match minimap2 stringency.
    ///
    /// # Quality Criteria (Ultra-Stringent for minimap2 Selectivity)
    /// - E-value ≤ 0.0001 (ultra-high significance, 100× more selective)
    /// - Alignment length ≥ 30 bp (substantial, reliable overlap)
    /// - Sequence complexity > 1.8 (aggressive filtering while preserving biological diversity)
    /// - Length-normalized score ≥ 1.8 (very high per-base quality requirement)
    ///
    /// These ultra-stringent thresholds are designed to produce alignment counts
    /// approaching minimap2's selectivity while preserving high-quality biological sequences.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::alignment::ScoringMatrix;
    ///
    /// let scoring = ScoringMatrix::default();
    /// let query = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 35bp, high complexity
    ///
    /// let passes = scoring.passes_quality_threshold(
    ///     80,    // very high score for ultra-stringent E-value
    ///     32,    // good length (≥30bp required)
    ///     100,   // query_len
    ///     50000, // smaller ref for better E-value
    ///     query
    /// );
    /// assert!(passes);
    /// ```
    pub fn passes_quality_threshold(
        &self,
        score: i32,
        alignment_len: usize,
        query_len: usize,
        ref_len: usize,
        sequence: &[u8],
    ) -> bool {
        // E-value significance test (ultra-stringent threshold - 100x more selective)
        let evalue = self.calculate_evalue(score, query_len, ref_len);
        if evalue > 0.0001 {
            return false;
        }

        // Minimum alignment length (ultra-stringent for minimap2-like selectivity)
        if alignment_len < 30 {
            return false;
        }

        // Sequence complexity filter (aggressive filtering while preserving biological diversity)
        let complexity = self.sequence_complexity(sequence);
        if complexity < 1.8 {
            return false;
        }

        // Length-normalized quality (ultra-stringent per-base requirement)
        let normalized_score = self.length_normalized_score(score, alignment_len);
        if normalized_score < 1.8 {
            return false;
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_scoring() {
        let scoring = ScoringMatrix::default();
        assert_eq!(scoring.match_score, 2);
        assert_eq!(scoring.mismatch_score, -1);
        assert_eq!(scoring.gap_open, -2);
        assert_eq!(scoring.gap_extend, -1);
    }

    #[test]
    fn test_custom_scoring() {
        let scoring = ScoringMatrix::new(5, -4, -10, -1);
        assert_eq!(scoring.match_score, 5);
        assert_eq!(scoring.mismatch_score, -4);
        assert_eq!(scoring.gap_open, -10);
        assert_eq!(scoring.gap_extend, -1);
    }

    #[test]
    fn test_score_match() {
        let scoring = ScoringMatrix::default();
        assert_eq!(scoring.score(b'A', b'A'), 2);
        assert_eq!(scoring.score(b'C', b'C'), 2);
        assert_eq!(scoring.score(b'G', b'G'), 2);
        assert_eq!(scoring.score(b'T', b'T'), 2);
    }

    #[test]
    fn test_score_mismatch() {
        let scoring = ScoringMatrix::default();
        assert_eq!(scoring.score(b'A', b'C'), -1);
        assert_eq!(scoring.score(b'A', b'G'), -1);
        assert_eq!(scoring.score(b'A', b'T'), -1);
        assert_eq!(scoring.score(b'C', b'G'), -1);
    }

    #[test]
    fn test_evalue_calculation() {
        let scoring = ScoringMatrix::default();

        // High-scoring alignment should have low E-value
        let high_score_evalue = scoring.calculate_evalue(100, 100, 1000000);
        assert!(high_score_evalue < 1e-25, "High score E-value: {}", high_score_evalue);

        // Low-scoring alignment should have high E-value
        let low_score_evalue = scoring.calculate_evalue(10, 100, 1000000);
        assert!(low_score_evalue > 10.0, "Low score E-value: {}", low_score_evalue);

        // Longer search space should have higher E-values for same score
        let short_evalue = scoring.calculate_evalue(50, 50, 1000);
        let long_evalue = scoring.calculate_evalue(50, 100, 1000000);
        assert!(long_evalue > short_evalue, "E-value should increase with search space");

        // Test significance threshold
        let moderate_score = scoring.calculate_evalue(40, 100, 1000000);
        assert!(moderate_score < 1.0, "Moderate score should be significant: {}", moderate_score);
    }

    #[test]
    fn test_length_normalized_score() {
        let scoring = ScoringMatrix::default();

        // Perfect alignment: score = 2 * length
        assert_eq!(scoring.length_normalized_score(20, 10), 2.0);
        assert_eq!(scoring.length_normalized_score(40, 20), 2.0);

        // Poor alignment: low score per base
        assert_eq!(scoring.length_normalized_score(10, 20), 0.5);

        // Edge case: zero length
        assert_eq!(scoring.length_normalized_score(50, 0), 0.0);
    }

    #[test]
    fn test_sequence_complexity() {
        let scoring = ScoringMatrix::default();

        // High complexity: balanced nucleotides
        let complex_seq = b"ACGTACGTACGTACGT";
        let complexity = scoring.sequence_complexity(complex_seq);
        assert!(complexity > 1.8, "High complexity sequence: {}", complexity);

        // Low complexity: homopolymer
        let simple_seq = b"AAAAAAAAAAAAAAAA";
        let simple_complexity = scoring.sequence_complexity(simple_seq);
        assert!(simple_complexity < 0.1, "Homopolymer complexity: {}", simple_complexity);

        // Medium complexity: two nucleotides
        let medium_seq = b"ATATATATATATATAT";
        let medium_complexity = scoring.sequence_complexity(medium_seq);
        assert!(medium_complexity > 0.8 && medium_complexity < 1.2,
                "Medium complexity: {}", medium_complexity);

        // Edge case: empty sequence
        assert_eq!(scoring.sequence_complexity(b""), 0.0);

        // Edge case: single nucleotide
        assert_eq!(scoring.sequence_complexity(b"A"), 0.0);
    }

    #[test]
    fn test_passes_quality_threshold() {
        let scoring = ScoringMatrix::default();

        // High-quality alignment should pass ultra-stringent thresholds
        let good_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 39bp, high complexity
        let passes = scoring.passes_quality_threshold(
            100,       // very high score for ultra-stringent E-value (0.0001)
            35,        // good length (≥30bp required)
            100,       // query length
            10000,     // smaller reference for better E-value
            good_seq
        );
        assert!(passes, "High-quality alignment should pass ultra-stringent filtering");

        // Low-scoring alignment should fail
        let fails_score = scoring.passes_quality_threshold(
            30,        // poor score (won't meet E-value 0.0001)
            35,        // good length
            100,       // query length
            10000,     // reference length
            good_seq
        );
        assert!(!fails_score, "Low-scoring alignment should fail E-value test");

        // Short alignment should fail (new 30bp minimum)
        let fails_length = scoring.passes_quality_threshold(
            100,       // good score
            25,        // too short (<30bp required)
            100,       // query length
            10000,     // reference length
            good_seq
        );
        assert!(!fails_length, "Short alignment should fail length test");

        // Low-complexity sequence should fail (1.8 complexity threshold)
        let repeat_seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 39bp homopolymer
        let fails_complexity = scoring.passes_quality_threshold(
            100,       // good score
            35,        // good length
            100,       // query length
            10000,     // reference length
            repeat_seq
        );
        assert!(!fails_complexity, "Low-complexity alignment should fail complexity test");

        // Poor per-base quality should fail (1.8 normalized score requirement)
        let mixed_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 39bp, good complexity
        let fails_quality = scoring.passes_quality_threshold(
            50,        // score too low (50/35 = 1.43 < 1.8)
            35,        // good length
            100,       // query length
            5000,      // small reference for good E-value
            mixed_seq
        );
        assert!(!fails_quality, "Poor per-base quality should fail normalized score test");
    }

    #[test]
    fn test_quality_threshold_edge_cases() {
        let scoring = ScoringMatrix::default();
        let good_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 41bp, high complexity

        // Borderline cases with ultra-stringent parameters
        let borderline_score = scoring.passes_quality_threshold(
            100,       // very high score for ultra-stringent E-value
            30,        // exactly minimum length
            100,       // query length
            3000,      // small reference for better E-value
            good_seq
        );
        // Should pass with good parameters under ultra-stringent filtering
        assert!(borderline_score, "Borderline case should pass with ultra-stringent parameters");

        // Test minimum length boundary (now 30bp)
        let min_length = scoring.passes_quality_threshold(
            100,       // very high score
            30,        // exactly minimum length (ultra-stringent: 30bp)
            100,       // query length
            3000,      // reference length
            good_seq
        );
        assert!(min_length, "Alignment at minimum length (30bp) should pass");

        // One below minimum length
        let below_min = scoring.passes_quality_threshold(
            100,       // very high score
            29,        // one below minimum (ultra-stringent: 30bp)
            100,       // query length
            3000,      // reference length
            good_seq
        );
        assert!(!below_min, "Alignment below minimum length (30bp) should fail");

        // Test complexity boundary (1.8 threshold)
        let moderate_complexity = b"ACGTACGTACGTACGTACGTACGTACGTACGTAAAAAAAAAA"; // 41bp, moderate complexity
        let complexity_score = scoring.sequence_complexity(moderate_complexity);
        let complexity_test = scoring.passes_quality_threshold(
            120,       // very high score for ultra-stringent requirements
            35,        // good length
            100,       // query length
            3000,      // reference length
            moderate_complexity
        );

        // This should depend on the actual complexity score
        if complexity_score >= 1.8 {
            assert!(complexity_test, "Sequences with complexity ≥1.8 should pass");
        } else {
            assert!(!complexity_test, "Sequences with complexity <1.8 should fail");
        }
    }
}
