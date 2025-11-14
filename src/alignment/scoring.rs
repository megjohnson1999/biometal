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
}
