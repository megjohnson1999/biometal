//! CIGAR (Compact Idiosyncratic Gapped Alignment Report) operations
//!
//! CIGAR strings represent sequence alignments compactly using operation codes.

/// CIGAR operation types
///
/// Represents the operations in a sequence alignment:
/// - Match: Aligned bases (may be match or mismatch)
/// - Insertion: Bases in query not in reference
/// - Deletion: Bases in reference not in query
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    /// M: Alignment match (length) - could be match or mismatch
    Match(usize),
    /// I: Insertion to reference (length)
    Insertion(usize),
    /// D: Deletion from reference (length)
    Deletion(usize),
}

impl CigarOp {
    /// Get the length of this operation
    pub fn len(&self) -> usize {
        match self {
            CigarOp::Match(n) | CigarOp::Insertion(n) | CigarOp::Deletion(n) => *n,
        }
    }

    /// Check if this operation has zero length (should not happen in valid CIGAR)
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Get the operation code as a character
    pub fn code(&self) -> char {
        match self {
            CigarOp::Match(_) => 'M',
            CigarOp::Insertion(_) => 'I',
            CigarOp::Deletion(_) => 'D',
        }
    }
}

impl std::fmt::Display for CigarOp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}", self.len(), self.code())
    }
}

/// Compress consecutive CIGAR operations of the same type
///
/// Converts a vector like `[M, M, M, I, I, D]` into `[3M, 2I, 1D]`
///
/// # Example
///
/// ```
/// use biometal::alignment::{CigarOp, compress_cigar};
///
/// let cigar = vec![
///     CigarOp::Match(1),
///     CigarOp::Match(1),
///     CigarOp::Match(1),
///     CigarOp::Insertion(1),
///     CigarOp::Insertion(1),
/// ];
///
/// let compressed = compress_cigar(cigar);
/// assert_eq!(compressed, vec![CigarOp::Match(3), CigarOp::Insertion(2)]);
/// ```
pub fn compress_cigar(cigar: Vec<CigarOp>) -> Vec<CigarOp> {
    if cigar.is_empty() {
        return cigar;
    }

    let mut compressed = Vec::new();
    let mut current = cigar[0];
    let mut count = current.len();

    for &op in &cigar[1..] {
        if std::mem::discriminant(&op) == std::mem::discriminant(&current) {
            // Same operation type, accumulate count
            count += op.len();
        } else {
            // Different operation, push accumulated operation
            compressed.push(match current {
                CigarOp::Match(_) => CigarOp::Match(count),
                CigarOp::Insertion(_) => CigarOp::Insertion(count),
                CigarOp::Deletion(_) => CigarOp::Deletion(count),
            });
            current = op;
            count = op.len();
        }
    }

    // Push final operation
    compressed.push(match current {
        CigarOp::Match(_) => CigarOp::Match(count),
        CigarOp::Insertion(_) => CigarOp::Insertion(count),
        CigarOp::Deletion(_) => CigarOp::Deletion(count),
    });

    compressed
}

/// Format CIGAR string for display
///
/// # Example
///
/// ```
/// use biometal::alignment::CigarOp;
///
/// let cigar = vec![CigarOp::Match(4), CigarOp::Insertion(2), CigarOp::Match(3)];
/// let cigar_str: String = cigar.iter().map(|op| op.to_string()).collect();
/// assert_eq!(cigar_str, "4M2I3M");
/// ```
pub fn format_cigar(cigar: &[CigarOp]) -> String {
    cigar.iter().map(|op| op.to_string()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cigar_op_len() {
        assert_eq!(CigarOp::Match(5).len(), 5);
        assert_eq!(CigarOp::Insertion(3).len(), 3);
        assert_eq!(CigarOp::Deletion(7).len(), 7);
    }

    #[test]
    fn test_cigar_op_code() {
        assert_eq!(CigarOp::Match(5).code(), 'M');
        assert_eq!(CigarOp::Insertion(3).code(), 'I');
        assert_eq!(CigarOp::Deletion(7).code(), 'D');
    }

    #[test]
    fn test_cigar_op_display() {
        assert_eq!(CigarOp::Match(5).to_string(), "5M");
        assert_eq!(CigarOp::Insertion(3).to_string(), "3I");
        assert_eq!(CigarOp::Deletion(7).to_string(), "7D");
    }

    #[test]
    fn test_compress_cigar_empty() {
        let cigar = vec![];
        let compressed = compress_cigar(cigar);
        assert_eq!(compressed, vec![]);
    }

    #[test]
    fn test_compress_cigar_single() {
        let cigar = vec![CigarOp::Match(5)];
        let compressed = compress_cigar(cigar);
        assert_eq!(compressed, vec![CigarOp::Match(5)]);
    }

    #[test]
    fn test_compress_cigar_consecutive() {
        let cigar = vec![
            CigarOp::Match(1),
            CigarOp::Match(1),
            CigarOp::Match(1),
        ];
        let compressed = compress_cigar(cigar);
        assert_eq!(compressed, vec![CigarOp::Match(3)]);
    }

    #[test]
    fn test_compress_cigar_mixed() {
        let cigar = vec![
            CigarOp::Match(1),
            CigarOp::Match(1),
            CigarOp::Insertion(1),
            CigarOp::Insertion(1),
            CigarOp::Deletion(1),
            CigarOp::Match(1),
        ];
        let compressed = compress_cigar(cigar);
        assert_eq!(
            compressed,
            vec![
                CigarOp::Match(2),
                CigarOp::Insertion(2),
                CigarOp::Deletion(1),
                CigarOp::Match(1),
            ]
        );
    }

    #[test]
    fn test_format_cigar() {
        let cigar = vec![CigarOp::Match(4), CigarOp::Insertion(2), CigarOp::Match(3)];
        assert_eq!(format_cigar(&cigar), "4M2I3M");
    }
}
