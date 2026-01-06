//! Streaming pattern matching workflows with alignment integration
//!
//! # Overview
//!
//! This module provides streaming pattern matching workflows that integrate
//! Smith-Waterman alignment with biometal's format parsers for real-world
//! bioinformatics tasks.
//!
//! # Use Cases
//!
//! - **Motif finding**: Locate regulatory sequences in streaming fashion
//! - **Primer design**: Find optimal primer binding sites
//! - **Adapter detection**: Identify sequencing adapters for trimming
//! - **Variant calling**: Local realignment around variant sites
//! - **Quality control**: Pattern-based sequence validation
//!
//! # Architecture
//!
//! Integrates with biometal's streaming parsers (FASTQ, FASTA, BAM) to provide
//! constant-memory pattern matching regardless of dataset size. Uses Smith-Waterman
//! for flexible, gapped alignment rather than exact string matching.
//!
//! # Example
//!
//! ```no_run
//! use biometal::alignment::pattern_matcher::{MotifFinder, MotifPattern};
//! use std::path::Path;
//!
//! let motifs = vec![
//!     MotifPattern::new("TATAAA", "TATA box"),
//!     MotifPattern::new("CAAT", "CAAT box"),
//! ];
//!
//! let mut finder = MotifFinder::new(motifs, 40); // min score = 40
//! for result in finder.find_in_fastq(Path::new("sequences.fastq"))? {
//!     let match_result = result?;
//!     println!("Found {} at position {}", match_result.motif_name, match_result.position);
//! }
//! ```

use crate::alignment::{smith_waterman, Alignment, ScoringMatrix};
use crate::{FastaStream, FastqStream, Result};
use std::io::BufRead;
use std::path::Path;

/// Pattern to search for in sequences
#[derive(Debug, Clone)]
pub struct MotifPattern {
    /// Pattern sequence (DNA/RNA)
    pub sequence: Vec<u8>,
    /// Human-readable name or description
    pub name: String,
    /// Reverse complement of the pattern (for bidirectional search)
    pub reverse_complement: Vec<u8>,
}

impl MotifPattern {
    /// Create a new motif pattern
    ///
    /// # Arguments
    ///
    /// * `sequence` - Pattern sequence (e.g., "TATAAA")
    /// * `name` - Descriptive name (e.g., "TATA box")
    ///
    /// # Returns
    ///
    /// New MotifPattern with reverse complement computed
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::alignment::pattern_matcher::MotifPattern;
    ///
    /// let motif = MotifPattern::new("TATAAA", "TATA box");
    /// assert_eq!(motif.sequence, b"TATAAA");
    /// assert_eq!(motif.reverse_complement, b"TTTATTAT"); // Reverse complement
    /// ```
    pub fn new(sequence: &str, name: &str) -> Self {
        let seq_bytes = sequence.as_bytes().to_vec();
        let reverse_comp = Self::reverse_complement(&seq_bytes);

        Self {
            sequence: seq_bytes,
            name: name.to_string(),
            reverse_complement: reverse_comp,
        }
    }

    /// Compute reverse complement of DNA sequence
    fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .rev()
            .map(|&base| match base {
                b'A' | b'a' => b'T',
                b'T' | b't' => b'A',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                _ => base, // Unknown bases unchanged
            })
            .collect()
    }
}

/// Result of motif matching
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MotifMatch {
    /// Query sequence identifier
    pub sequence_id: String,
    /// Matched motif name
    pub motif_name: String,
    /// Alignment score
    pub score: i32,
    /// Position in sequence (0-indexed)
    pub position: usize,
    /// Length of match
    pub length: usize,
    /// Strand (+ for forward, - for reverse complement)
    pub strand: char,
    /// Detailed alignment information
    pub alignment: Alignment,
}

/// Streaming motif finder with constant memory
///
/// Finds motifs in streaming fashion using Smith-Waterman alignment.
/// Supports both forward and reverse complement matching.
pub struct MotifFinder {
    /// Patterns to search for
    patterns: Vec<MotifPattern>,
    /// Minimum alignment score threshold
    min_score: i32,
    /// Scoring matrix for alignment
    scoring: ScoringMatrix,
}

impl MotifFinder {
    /// Create a new motif finder
    ///
    /// # Arguments
    ///
    /// * `patterns` - Motif patterns to search for
    /// * `min_score` - Minimum alignment score to report
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::alignment::pattern_matcher::{MotifFinder, MotifPattern};
    ///
    /// let patterns = vec![
    ///     MotifPattern::new("TATAAA", "TATA box"),
    ///     MotifPattern::new("CAAT", "CAAT box"),
    /// ];
    ///
    /// let finder = MotifFinder::new(patterns, 30);
    /// ```
    pub fn new(patterns: Vec<MotifPattern>, min_score: i32) -> Self {
        Self {
            patterns,
            min_score,
            scoring: ScoringMatrix::default(),
        }
    }

    /// Create motif finder with custom scoring matrix
    ///
    /// # Arguments
    ///
    /// * `patterns` - Motif patterns to search for
    /// * `min_score` - Minimum alignment score threshold
    /// * `scoring` - Custom scoring matrix
    pub fn with_scoring(patterns: Vec<MotifPattern>, min_score: i32, scoring: ScoringMatrix) -> Self {
        Self {
            patterns,
            min_score,
            scoring,
        }
    }

    /// Find motifs in FASTQ file with streaming
    ///
    /// # Arguments
    ///
    /// * `fastq_path` - Path to FASTQ file
    ///
    /// # Returns
    ///
    /// Iterator over motif matches found in the file
    ///
    /// # Memory Usage
    ///
    /// - Constant ~5MB regardless of file size
    /// - Processes one read at a time
    /// - No accumulation of results in memory
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::alignment::pattern_matcher::{MotifFinder, MotifPattern};
    /// use std::path::Path;
    ///
    /// let patterns = vec![MotifPattern::new("TATAAA", "TATA box")];
    /// let mut finder = MotifFinder::new(patterns, 40);
    ///
    /// for result in finder.find_in_fastq(Path::new("reads.fastq"))? {
    ///     let motif_match = result?;
    ///     println!("Found {} in {}", motif_match.motif_name, motif_match.sequence_id);
    /// }
    /// ```
    pub fn find_in_fastq<P: AsRef<Path>>(
        &self,
        fastq_path: P,
    ) -> Result<MotifMatchIterator<impl BufRead>> {
        let fastq_stream = FastqStream::from_path(fastq_path)?;
        Ok(MotifMatchIterator::new_fastq(
            fastq_stream,
            self.patterns.clone(),
            self.min_score,
            self.scoring.clone(),
        ))
    }

    /// Find motifs in FASTA file with streaming
    ///
    /// # Arguments
    ///
    /// * `fasta_path` - Path to FASTA file
    ///
    /// # Returns
    ///
    /// Iterator over motif matches found in the file
    pub fn find_in_fasta<P: AsRef<Path>>(
        &self,
        fasta_path: P,
    ) -> Result<MotifMatchIterator<impl BufRead>> {
        let fasta_stream = FastaStream::from_path(fasta_path)?;
        Ok(MotifMatchIterator::new_fasta(
            fasta_stream,
            self.patterns.clone(),
            self.min_score,
            self.scoring.clone(),
        ))
    }

    /// Find motifs in a single sequence
    ///
    /// # Arguments
    ///
    /// * `sequence_id` - Identifier for the sequence
    /// * `sequence` - DNA/RNA sequence bytes
    ///
    /// # Returns
    ///
    /// Vector of all motif matches found in the sequence
    pub fn find_in_sequence(&self, sequence_id: &str, sequence: &[u8]) -> Vec<MotifMatch> {
        let mut matches = Vec::new();

        for pattern in &self.patterns {
            // Search forward strand
            if let Some(match_result) = self.find_pattern_in_sequence(
                sequence_id,
                sequence,
                &pattern.sequence,
                &pattern.name,
                '+',
            ) {
                matches.push(match_result);
            }

            // Search reverse complement strand
            if let Some(match_result) = self.find_pattern_in_sequence(
                sequence_id,
                sequence,
                &pattern.reverse_complement,
                &pattern.name,
                '-',
            ) {
                matches.push(match_result);
            }
        }

        matches
    }

    /// Find a specific pattern in sequence using Smith-Waterman
    fn find_pattern_in_sequence(
        &self,
        sequence_id: &str,
        sequence: &[u8],
        pattern: &[u8],
        pattern_name: &str,
        strand: char,
    ) -> Option<MotifMatch> {
        let alignment = smith_waterman(pattern, sequence, &self.scoring);

        if alignment.score >= self.min_score {
            Some(MotifMatch {
                sequence_id: sequence_id.to_string(),
                motif_name: pattern_name.to_string(),
                score: alignment.score,
                position: alignment.ref_start,
                length: alignment.ref_end - alignment.ref_start,
                strand,
                alignment,
            })
        } else {
            None
        }
    }
}

/// Iterator for motif matches from streaming sources
pub struct MotifMatchIterator<R: BufRead> {
    source: MotifMatchSource<R>,
    patterns: Vec<MotifPattern>,
    min_score: i32,
    scoring: ScoringMatrix,
    current_matches: Vec<MotifMatch>,
    match_index: usize,
}

/// Source of sequences for motif matching
enum MotifMatchSource<R: BufRead> {
    Fastq(FastqStream<R>),
    Fasta(FastaStream<R>),
}

impl<R: BufRead> MotifMatchIterator<R> {
    fn new_fastq(
        fastq_stream: FastqStream<R>,
        patterns: Vec<MotifPattern>,
        min_score: i32,
        scoring: ScoringMatrix,
    ) -> Self {
        Self {
            source: MotifMatchSource::Fastq(fastq_stream),
            patterns,
            min_score,
            scoring,
            current_matches: Vec::new(),
            match_index: 0,
        }
    }

    fn new_fasta(
        fasta_stream: FastaStream<R>,
        patterns: Vec<MotifPattern>,
        min_score: i32,
        scoring: ScoringMatrix,
    ) -> Self {
        Self {
            source: MotifMatchSource::Fasta(fasta_stream),
            patterns,
            min_score,
            scoring,
            current_matches: Vec::new(),
            match_index: 0,
        }
    }

    /// Load motif matches from next sequence record
    fn load_next_sequence(&mut self) -> Result<bool> {
        self.current_matches.clear();
        self.match_index = 0;

        match &mut self.source {
            MotifMatchSource::Fastq(stream) => {
                if let Some(record_result) = stream.next() {
                    let record = record_result?;
                    let finder = MotifFinder {
                        patterns: self.patterns.clone(),
                        min_score: self.min_score,
                        scoring: self.scoring.clone(),
                    };
                    self.current_matches = finder.find_in_sequence(&record.id, &record.sequence);
                    Ok(true)
                } else {
                    Ok(false) // No more records
                }
            }
            MotifMatchSource::Fasta(stream) => {
                if let Some(record_result) = stream.next() {
                    let record = record_result?;
                    let finder = MotifFinder {
                        patterns: self.patterns.clone(),
                        min_score: self.min_score,
                        scoring: self.scoring.clone(),
                    };
                    self.current_matches = finder.find_in_sequence(&record.id, &record.sequence);
                    Ok(true)
                } else {
                    Ok(false) // No more records
                }
            }
        }
    }
}

impl<R: BufRead> Iterator for MotifMatchIterator<R> {
    type Item = Result<MotifMatch>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // Return next match from current sequence
            if self.match_index < self.current_matches.len() {
                let match_result = self.current_matches[self.match_index].clone();
                self.match_index += 1;
                return Some(Ok(match_result));
            }

            // Load next sequence
            match self.load_next_sequence() {
                Ok(true) => continue,  // Successfully loaded, check for matches
                Ok(false) => return None, // No more sequences
                Err(e) => return Some(Err(e)), // Error loading sequence
            }
        }
    }
}

/// Primer finding functionality
pub struct PrimerFinder {
    motif_finder: MotifFinder,
}

impl PrimerFinder {
    /// Create primer finder for common primer sites
    ///
    /// Searches for common primer binding sites used in PCR and sequencing.
    pub fn new_standard() -> Self {
        let patterns = vec![
            MotifPattern::new("GTTTCCCAGTCACGAC", "M13 Forward"),
            MotifPattern::new("CAGGAAACAGCTATGAC", "M13 Reverse"),
            MotifPattern::new("TGTAAAACGACGGCCAGT", "M13F (-20)"),
            MotifPattern::new("AACAGCTATGACCATG", "M13R (-24)"),
        ];

        Self {
            motif_finder: MotifFinder::new(patterns, 50), // High stringency for primers
        }
    }

    /// Find primers in FASTQ file
    pub fn find_primers_fastq<P: AsRef<Path>>(
        &self,
        fastq_path: P,
    ) -> Result<MotifMatchIterator<impl BufRead>> {
        self.motif_finder.find_in_fastq(fastq_path)
    }
}

/// Adapter detection functionality
pub struct AdapterDetector {
    motif_finder: MotifFinder,
}

impl AdapterDetector {
    /// Create adapter detector for Illumina adapters
    ///
    /// Detects common Illumina sequencing adapters for quality control.
    pub fn new_illumina() -> Self {
        let patterns = vec![
            MotifPattern::new("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "Illumina Universal"),
            MotifPattern::new("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", "Illumina Small RNA 3'"),
            MotifPattern::new("TGGAATTCTCGGGTGCCAAGG", "Illumina Small RNA 5'"),
        ];

        Self {
            motif_finder: MotifFinder::new(patterns, 60), // Very high stringency
        }
    }

    /// Detect adapters in FASTQ file
    pub fn detect_adapters_fastq<P: AsRef<Path>>(
        &self,
        fastq_path: P,
    ) -> Result<MotifMatchIterator<impl BufRead>> {
        self.motif_finder.find_in_fastq(fastq_path)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_motif_pattern_creation() {
        let motif = MotifPattern::new("TATAAA", "TATA box");
        assert_eq!(motif.sequence, b"TATAAA");
        assert_eq!(motif.name, "TATA box");
        assert_eq!(motif.reverse_complement, b"TTTATA");
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(
            MotifPattern::reverse_complement(b"ACGT"),
            b"ACGT"
        );
        assert_eq!(
            MotifPattern::reverse_complement(b"AAAA"),
            b"TTTT"
        );
        assert_eq!(
            MotifPattern::reverse_complement(b"ATCG"),
            b"CGAT"
        );
    }

    #[test]
    fn test_motif_finder_creation() {
        let patterns = vec![
            MotifPattern::new("TATAAA", "TATA box"),
            MotifPattern::new("CAAT", "CAAT box"),
        ];

        let finder = MotifFinder::new(patterns, 30);
        assert_eq!(finder.patterns.len(), 2);
        assert_eq!(finder.min_score, 30);
    }

    #[test]
    fn test_find_in_sequence() {
        let patterns = vec![MotifPattern::new("TATA", "TATA box")];
        let finder = MotifFinder::new(patterns, 6); // 4 matches × 2 = 8 > 6

        let sequence = b"ACGTTATAACGT";
        let matches = finder.find_in_sequence("test", sequence);

        assert!(!matches.is_empty());
        assert_eq!(matches[0].motif_name, "TATA box");
        assert_eq!(matches[0].strand, '+');
        assert!(matches[0].position > 0);
    }

    #[test]
    fn test_primer_finder() {
        let finder = PrimerFinder::new_standard();
        assert_eq!(finder.motif_finder.patterns.len(), 4);
    }

    #[test]
    fn test_adapter_detector() {
        let detector = AdapterDetector::new_illumina();
        assert_eq!(detector.motif_finder.patterns.len(), 3);
    }

    #[test]
    fn test_custom_scoring() {
        let patterns = vec![MotifPattern::new("ACGT", "test")];
        let scoring = ScoringMatrix {
            match_score: 3,
            mismatch_score: -2,
            gap_open: -3,
            gap_extend: -1,
        };

        let finder = MotifFinder::with_scoring(patterns, 10, scoring);
        assert_eq!(finder.scoring.match_score, 3);
    }
}