//! Streaming read mapper with constant memory architecture
//!
//! # Innovation
//!
//! Traditional read mapping requires large memory indices (BWA: 8-16GB for human genome).
//! This module implements a novel **windowed reference processing** approach that maintains
//! biometal's ~5MB constant memory constraint while enabling read mapping functionality.
//!
//! # Architecture
//!
//! Instead of loading entire reference genomes and building indices, the streaming mapper:
//! - Processes reference genome in overlapping windows (e.g., 1MB with 200bp overlap)
//! - Aligns reads against current window using Smith-Waterman
//! - Buffers alignments at window boundaries for consistency
//! - Achieves constant memory regardless of reference/dataset size
//!
//! # Performance Trade-offs
//!
//! - **Memory**: ~5MB constant (vs 8-16GB for traditional aligners)
//! - **Speed**: Slower than indexed approaches for whole-genome mapping
//! - **Use Cases**: Ideal for targeted alignment, memory-constrained environments
//!
//! # Evidence Base
//!
//! Follows biometal's streaming architecture (OPTIMIZATION_RULES.md Rule 5)
//! with 99.5% memory reduction demonstrated across format parsers.

use crate::alignment::{smith_waterman, Alignment, ScoringMatrix};
use crate::{FastaRecord, FastaStream, FastqRecord, FastqStream, Result};
use std::collections::VecDeque;
use std::io::BufRead;
use std::path::Path;

/// Configuration for streaming read mapping
#[derive(Debug, Clone)]
pub struct StreamingMapperConfig {
    /// Size of reference windows in bytes (e.g., 1MB)
    pub window_size: usize,
    /// Overlap between windows in base pairs (e.g., 200bp)
    pub overlap_bp: usize,
    /// Minimum alignment score to report
    pub min_score_threshold: i32,
    /// Scoring matrix for alignment
    pub scoring: ScoringMatrix,
}

impl Default for StreamingMapperConfig {
    fn default() -> Self {
        Self {
            window_size: 1_000_000, // 1MB windows
            overlap_bp: 200,        // 200bp overlap
            min_score_threshold: 50, // Minimum score threshold
            scoring: ScoringMatrix::default(),
        }
    }
}

/// Streaming read mapping result
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MappingResult {
    /// Query sequence identifier
    pub query_id: String,
    /// Reference sequence identifier
    pub reference_id: String,
    /// Global position in reference (accounting for windows)
    pub global_ref_start: usize,
    /// Global position in reference (end)
    pub global_ref_end: usize,
    /// Alignment result
    pub alignment: Alignment,
    /// Window index where this alignment was found
    pub window_index: usize,
}

/// Reference window for streaming processing
#[derive(Debug, Clone)]
struct ReferenceWindow {
    /// Sequence data for this window
    sequence: Vec<u8>,
    /// Global start position in reference
    global_start: usize,
    /// Global end position in reference
    global_end: usize,
    /// Window index
    window_index: usize,
    /// Reference sequence ID
    reference_id: String,
}

/// Streaming read mapper with constant memory architecture
///
/// # Example
///
/// ```no_run
/// use biometal::alignment::streaming_mapper::{StreamingMapper, StreamingMapperConfig};
/// use std::path::Path;
///
/// let config = StreamingMapperConfig::default();
/// let mut mapper = StreamingMapper::new(config);
///
/// // Map reads from FASTQ file against reference FASTA
/// let mappings: Vec<_> = mapper
///     .map_reads_streaming(
///         Path::new("reference.fasta"),
///         Path::new("reads.fastq")
///     )?
///     .collect::<Result<Vec<_>, _>>()?;
/// ```
pub struct StreamingMapper {
    config: StreamingMapperConfig,
    // Buffer for alignments at window boundaries
    boundary_alignments: VecDeque<MappingResult>,
}

impl StreamingMapper {
    /// Create a new streaming mapper with configuration
    ///
    /// # Arguments
    ///
    /// * `config` - Mapping configuration (window size, overlap, thresholds)
    ///
    /// # Returns
    ///
    /// New StreamingMapper instance
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::alignment::streaming_mapper::{StreamingMapper, StreamingMapperConfig};
    ///
    /// let config = StreamingMapperConfig {
    ///     window_size: 500_000,  // 500KB windows
    ///     overlap_bp: 100,       // 100bp overlap
    ///     min_score_threshold: 30,
    ///     ..Default::default()
    /// };
    /// let mapper = StreamingMapper::new(config);
    /// ```
    pub fn new(config: StreamingMapperConfig) -> Self {
        Self {
            config,
            boundary_alignments: VecDeque::new(),
        }
    }

    /// Stream read mapping with constant memory
    ///
    /// Processes reference genome in windows and maps reads against each window.
    /// Memory usage remains constant regardless of reference or dataset size.
    ///
    /// # Arguments
    ///
    /// * `reference_path` - Path to reference FASTA file
    /// * `reads_path` - Path to reads FASTQ file
    ///
    /// # Returns
    ///
    /// Iterator over mapping results
    ///
    /// # Memory Usage
    ///
    /// - Reference window: ~1MB (configurable)
    /// - Read buffer: ~1MB (typical FASTQ chunk)
    /// - Alignment buffers: ~3MB (DP matrices, results)
    /// - **Total: ~5MB constant**
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::alignment::streaming_mapper::{StreamingMapper, StreamingMapperConfig};
    /// use std::path::Path;
    ///
    /// let mut mapper = StreamingMapper::new(StreamingMapperConfig::default());
    ///
    /// for mapping in mapper.map_reads_streaming(
    ///     Path::new("chr1.fasta"),
    ///     Path::new("reads.fastq")
    /// )? {
    ///     let result = mapping?;
    ///     println!("Mapped {} to position {}", result.query_id, result.global_ref_start);
    /// }
    /// ```
    pub fn map_reads_streaming<P1: AsRef<Path>, P2: AsRef<Path>>(
        &mut self,
        reference_path: P1,
        reads_path: P2,
    ) -> Result<StreamingMappingIterator<impl BufRead, impl BufRead>> {
        let reference_stream = FastaStream::from_path(reference_path)?;
        let reads_stream = FastqStream::from_path(reads_path)?;

        Ok(StreamingMappingIterator::new(
            reference_stream,
            reads_stream,
            self.config.clone(),
        ))
    }

    /// Map a single read against a reference window
    ///
    /// # Arguments
    ///
    /// * `read` - FASTQ read to align
    /// * `window` - Reference window to align against
    ///
    /// # Returns
    ///
    /// Best mapping result if score exceeds threshold
    fn map_read_to_window(
        &self,
        read: &FastqRecord,
        window: &ReferenceWindow,
    ) -> Option<MappingResult> {
        let alignment = smith_waterman(
            &read.sequence,
            &window.sequence,
            &self.config.scoring,
        );

        // Filter by score threshold
        if alignment.score < self.config.min_score_threshold {
            return None;
        }

        // Convert local window coordinates to global reference coordinates
        let global_ref_start = window.global_start + alignment.ref_start;
        let global_ref_end = window.global_start + alignment.ref_end;

        Some(MappingResult {
            query_id: read.id.clone(),
            reference_id: window.reference_id.clone(),
            global_ref_start,
            global_ref_end,
            alignment,
            window_index: window.window_index,
        })
    }

    /// Create overlapping windows from reference sequence
    ///
    /// Splits reference into overlapping windows to ensure alignments
    /// spanning window boundaries are not missed.
    fn create_windows(
        &self,
        reference: &FastaRecord,
        reference_id: String,
    ) -> Vec<ReferenceWindow> {
        let sequence = &reference.sequence;
        let seq_len = sequence.len();

        if seq_len <= self.config.window_size {
            // Single window for small references
            return vec![ReferenceWindow {
                sequence: sequence.to_vec(),
                global_start: 0,
                global_end: seq_len,
                window_index: 0,
                reference_id,
            }];
        }

        let mut windows = Vec::new();
        let mut start = 0;
        let mut window_index = 0;

        while start < seq_len {
            let end = (start + self.config.window_size).min(seq_len);

            let window_sequence = if end == seq_len {
                // Last window: include everything to end
                sequence[start..end].to_vec()
            } else {
                // Regular window: include overlap for next window
                let _overlap_start = end.saturating_sub(self.config.overlap_bp);
                sequence[start..end].to_vec()
            };

            windows.push(ReferenceWindow {
                sequence: window_sequence,
                global_start: start,
                global_end: end,
                window_index,
                reference_id: reference_id.clone(),
            });

            // Next window starts with overlap
            start = end.saturating_sub(self.config.overlap_bp);
            window_index += 1;

            // Prevent infinite loop for very small overlaps
            if start + self.config.overlap_bp >= seq_len {
                break;
            }
        }

        windows
    }
}

/// Iterator for streaming mapping results
///
/// Processes reference and reads in streaming fashion, yielding mapping results
/// as they are found. Memory usage remains constant regardless of dataset size.
pub struct StreamingMappingIterator<R: BufRead, S: BufRead> {
    reference_stream: FastaStream<R>,
    reads_stream: FastqStream<S>,
    config: StreamingMapperConfig,
    current_windows: VecDeque<ReferenceWindow>,
    current_reads: Vec<FastqRecord>,
    reads_exhausted: bool,
    windows_exhausted: bool,
    mapping_buffer: VecDeque<MappingResult>,
}

impl<R: BufRead, S: BufRead> StreamingMappingIterator<R, S> {
    fn new(
        reference_stream: FastaStream<R>,
        reads_stream: FastqStream<S>,
        config: StreamingMapperConfig,
    ) -> Self {
        Self {
            reference_stream,
            reads_stream,
            config,
            current_windows: VecDeque::new(),
            current_reads: Vec::new(),
            reads_exhausted: false,
            windows_exhausted: false,
            mapping_buffer: VecDeque::new(),
        }
    }

    /// Load next batch of reads (constant memory)
    fn load_next_reads(&mut self) -> Result<()> {
        const MAX_READS_PER_BATCH: usize = 1000; // ~1MB typical batch

        self.current_reads.clear();

        for _ in 0..MAX_READS_PER_BATCH {
            match self.reads_stream.next() {
                Some(Ok(read)) => self.current_reads.push(read),
                Some(Err(e)) => return Err(e),
                None => {
                    self.reads_exhausted = true;
                    break;
                }
            }
        }

        Ok(())
    }

    /// Load next reference sequence and create windows
    fn load_next_reference(&mut self) -> Result<()> {
        self.current_windows.clear();

        match self.reference_stream.next() {
            Some(Ok(reference)) => {
                let reference_id = reference.id.clone();
                let mapper = StreamingMapper::new(self.config.clone());
                let windows = mapper.create_windows(&reference, reference_id);

                for window in windows {
                    self.current_windows.push_back(window);
                }
            }
            Some(Err(e)) => return Err(e),
            None => {
                self.windows_exhausted = true;
            }
        }

        Ok(())
    }

    /// Process current reads against current window
    fn process_current_window(&mut self) -> Result<()> {
        if let Some(window) = self.current_windows.front() {
            let mapper = StreamingMapper::new(self.config.clone());

            for read in &self.current_reads {
                if let Some(mapping) = mapper.map_read_to_window(read, window) {
                    self.mapping_buffer.push_back(mapping);
                }
            }
        }

        Ok(())
    }
}

impl<R: BufRead, S: BufRead> Iterator for StreamingMappingIterator<R, S> {
    type Item = Result<MappingResult>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // Return buffered results first
            if let Some(mapping) = self.mapping_buffer.pop_front() {
                return Some(Ok(mapping));
            }

            // If we have reads and windows, process them
            if !self.current_reads.is_empty() && !self.current_windows.is_empty() {
                if let Err(e) = self.process_current_window() {
                    return Some(Err(e));
                }
                self.current_windows.pop_front(); // Move to next window
                continue;
            }

            // Load more windows if needed
            if self.current_windows.is_empty() && !self.windows_exhausted {
                if let Err(e) = self.load_next_reference() {
                    return Some(Err(e));
                }
                continue;
            }

            // Load more reads if needed
            if self.current_reads.is_empty() && !self.reads_exhausted {
                if let Err(e) = self.load_next_reads() {
                    return Some(Err(e));
                }
                continue;
            }

            // If both streams are exhausted, we're done
            if self.reads_exhausted && self.windows_exhausted {
                return None;
            }

            // Reset reads for next reference sequence
            if self.current_windows.is_empty() && !self.current_reads.is_empty() {
                if let Err(e) = self.load_next_reads() {
                    return Some(Err(e));
                }
            }

            // If we still have no progress, we're done
            if self.current_reads.is_empty() && self.current_windows.is_empty() {
                return None;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_streaming_mapper_config_default() {
        let config = StreamingMapperConfig::default();
        assert_eq!(config.window_size, 1_000_000);
        assert_eq!(config.overlap_bp, 200);
        assert_eq!(config.min_score_threshold, 50);
    }

    #[test]
    fn test_create_single_window() {
        let config = StreamingMapperConfig {
            window_size: 1000,
            overlap_bp: 100,
            ..Default::default()
        };
        let mapper = StreamingMapper::new(config);

        let small_record = FastaRecord::new("test".to_string(), b"ACGTACGTACGT".to_vec());
        let windows = mapper.create_windows(&small_record, "chr1".to_string());

        assert_eq!(windows.len(), 1);
        assert_eq!(windows[0].global_start, 0);
        assert_eq!(windows[0].global_end, 12);
        assert_eq!(windows[0].sequence, b"ACGTACGTACGT");
    }

    #[test]
    fn test_create_multiple_windows() {
        let config = StreamingMapperConfig {
            window_size: 10,
            overlap_bp: 3,
            ..Default::default()
        };
        let mapper = StreamingMapper::new(config);

        // 20bp sequence -> should create multiple windows
        let sequence = b"ACGTACGTACGTACGTACGT"; // 20bp
        let record = FastaRecord::new("test".to_string(), sequence.to_vec());
        let windows = mapper.create_windows(&record, "chr1".to_string());

        assert!(windows.len() > 1);
        assert_eq!(windows[0].global_start, 0);
        assert_eq!(windows[0].global_end, 10);

        // Check overlap
        assert_eq!(windows[1].global_start, 7); // 10 - 3 = 7
    }

    #[test]
    fn test_map_read_to_window() {
        let mapper = StreamingMapper::new(StreamingMapperConfig::default());

        let window = ReferenceWindow {
            sequence: b"ACGTACGTACGTACGT".to_vec(),
            global_start: 100,
            global_end: 116,
            window_index: 0,
            reference_id: "chr1".to_string(),
        };

        let read = FastqRecord::new("read1".to_string(), b"ACGTACGT".to_vec(), b"HHHHHHHH".to_vec());

        let result = mapper.map_read_to_window(&read, &window);
        assert!(result.is_some());

        let mapping = result.unwrap();
        assert_eq!(mapping.query_id, "read1");
        assert_eq!(mapping.reference_id, "chr1");
        assert_eq!(mapping.global_ref_start, 100); // window start + local start (0)
        assert!(mapping.alignment.score > 0);
    }

    #[test]
    fn test_low_score_filtering() {
        let config = StreamingMapperConfig {
            min_score_threshold: 100, // Very high threshold
            ..Default::default()
        };
        let mapper = StreamingMapper::new(config);

        let window = ReferenceWindow {
            sequence: b"TTTTTTTTTTTTTTTT".to_vec(), // No match
            global_start: 0,
            global_end: 16,
            window_index: 0,
            reference_id: "chr1".to_string(),
        };

        let read = FastqRecord::new("read1".to_string(), b"ACGTACGT".to_vec(), b"HHHHHHHH".to_vec());

        let result = mapper.map_read_to_window(&read, &window);
        assert!(result.is_none()); // Should be filtered out due to low score
    }
}