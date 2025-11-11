//! BAM streaming reader.
//!
//! Provides a streaming interface for reading BAM files with constant memory.
//! Follows biometal's Rule 5 (streaming architecture, ~5 MB memory).
//!
//! # Design
//!
//! - Iterator-based interface (constant memory)
//! - Header read separately from records
//! - Records are not accumulated (streaming)
//! - **Phase 2**: Parallel BGZF decompression (6.5× speedup, Rule 3)
//!
//! # Usage
//!
//! ```no_run
//! use biometal::io::bam::BamReader;
//!
//! # fn main() -> biometal::Result<()> {
//! // Automatically uses parallel BGZF decompression for .bam files
//! let mut bam = BamReader::from_path("alignments.bam")?;
//!
//! println!("Header: {} references", bam.header().reference_count());
//!
//! for result in bam.records() {
//!     let record = result?;
//!     println!("{} at {}", record.name, record.position.unwrap_or(-1));
//! }
//! # Ok(())
//! # }
//! ```

use super::header::{read_header, Header};
use super::index::{BaiIndex, Chunk};
use super::record::{parse_record, Record};
use crate::io::compression::{CompressedReader, DataSource};
use std::fs::File;
use std::io::{self, BufRead, Read};
use std::path::Path;

/// BAM file reader with streaming interface.
///
/// Reads BAM files with constant memory (Rule 5: streaming architecture).
/// The header is read once during construction, then records are streamed.
///
/// # Buffer Reuse
///
/// Maintains an internal buffer that's reused across record reads to avoid
/// repeated allocations. This buffer grows to accommodate the largest record
/// seen, then stays at that size for subsequent reads.
pub struct BamReader<R> {
    /// Underlying reader (typically BufReader<File>)
    reader: R,
    /// BAM header (read during construction)
    header: Header,
    /// Reusable buffer for reading record data (avoids repeated allocations)
    buffer: Vec<u8>,
}

impl<R: BufRead> BamReader<R> {
    /// Create a new BAM reader.
    ///
    /// Reads and validates the BAM header immediately.
    ///
    /// # Errors
    ///
    /// Returns error if:
    /// - Cannot read header
    /// - Invalid magic bytes
    /// - Header is malformed
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::bam::BamReader;
    /// use std::fs::File;
    /// use std::io::BufReader;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let file = File::open("alignments.bam")?;
    /// let reader = BufReader::new(file);
    /// let bam = BamReader::new(reader)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(mut reader: R) -> io::Result<Self> {
        let header = read_header(&mut reader)?;
        // Initialize buffer with reasonable capacity (typical record ~500 bytes)
        let buffer = Vec::with_capacity(512);
        Ok(Self {
            reader,
            header,
            buffer,
        })
    }

    /// Get a reference to the BAM header.
    ///
    /// The header contains SAM header text and reference sequence information.
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Create an iterator over BAM records.
    ///
    /// Records are streamed with constant memory (not accumulated).
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::bam::BamReader;
    /// use std::fs::File;
    /// use std::io::BufReader;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let file = File::open("alignments.bam")?;
    /// let reader = BufReader::new(file);
    /// let mut bam = BamReader::new(reader)?;
    ///
    /// for result in bam.records() {
    ///     let record = result?;
    ///     // Process one record at a time (constant memory)
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn records(&mut self) -> Records<'_, R> {
        Records { reader: self }
    }

    /// Read a single record.
    ///
    /// Returns `Ok(None)` when EOF is reached.
    ///
    /// # Buffer Reuse
    ///
    /// This method reuses an internal buffer to avoid repeated allocations.
    /// The buffer grows to accommodate larger records but is not shrunk,
    /// making subsequent reads more efficient.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::bam::BamReader;
    /// use std::fs::File;
    /// use std::io::BufReader;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let file = File::open("alignments.bam")?;
    /// let reader = BufReader::new(file);
    /// let mut bam = BamReader::new(reader)?;
    ///
    /// while let Some(record) = bam.read_record()? {
    ///     println!("Read: {}", record.name);
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn read_record(&mut self) -> io::Result<Option<Record>> {
        // Read block size (4 bytes, little-endian)
        let mut size_buf = [0u8; 4];
        match self.reader.read_exact(&mut size_buf) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
            Err(e) => return Err(e),
        }

        let block_size = i32::from_le_bytes(size_buf);
        if block_size < 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Invalid block size: {}", block_size),
            ));
        }

        let block_size = block_size as usize;

        // Reuse buffer, resizing if needed to fit block_size + 4 bytes for the size itself
        self.buffer.clear();
        self.buffer.reserve(block_size + 4);

        // Put block size at start (parse_record expects it)
        self.buffer.extend_from_slice(&size_buf);

        unsafe {
            // SAFETY: We're immediately reading into this memory
            let start_len = self.buffer.len();
            self.buffer.set_len(start_len + block_size);
        }

        // Read record data into buffer (after the block size)
        self.reader.read_exact(&mut self.buffer[4..])?;

        // Parse record from buffer (includes block size at start)
        let record = parse_record(&self.buffer)?;
        Ok(Some(record))
    }
}

impl BamReader<CompressedReader> {
    /// Open a BAM file from a path with automatic parallel BGZF decompression.
    ///
    /// **Phase 2**: Integrates parallel BGZF decompression (Rule 3: 6.5× speedup).
    ///
    /// # Performance (Rule 3)
    ///
    /// - Automatically detects BGZF compression (peeks at magic bytes)
    /// - If compressed: Uses parallel decompression (8 blocks, 6.5× speedup)
    /// - If uncompressed: Direct passthrough (no overhead)
    /// - Maintains constant ~1 MB memory (Rule 5)
    ///
    /// # Evidence
    ///
    /// - Rule 3 (Parallel BGZF): Entry 029, validated 6.5× speedup
    /// - Expected overall speedup: ~4-5× for typical BAM parsing workloads
    ///
    /// # Errors
    ///
    /// Returns error if:
    /// - File cannot be opened
    /// - Header is invalid
    /// - Decompression fails
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::bam::BamReader;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// // Automatically uses parallel BGZF for .bam files (6.5× faster)
    /// let bam = BamReader::from_path("alignments.bam")?;
    /// println!("Opened BAM with {} references", bam.header().reference_count());
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_path<P: AsRef<Path>>(path: P) -> crate::Result<Self> {
        // Create data source (supports local files, Rule 6 future: HTTP/SRA)
        let source = DataSource::from_path(path);

        // Create compressed reader with automatic BGZF detection
        // - Peeks at magic bytes (31, 139) to detect gzip
        // - If BGZF: Parallel decompression (Rule 3: 6.5× speedup)
        // - If uncompressed: Direct passthrough
        let reader = CompressedReader::new(source)?;

        // Create BAM reader from compressed reader
        // - io::Error automatically converts to BiometalError via From trait
        Ok(Self::new(reader)?)
    }

    /// Plan a region query using an index.
    ///
    /// Returns the chunks (virtual file offset ranges) that need to be read
    /// to retrieve all alignments overlapping the specified region.
    ///
    /// # Phase 1 API
    ///
    /// This is a Phase 1 implementation that returns the query plan (chunks to read)
    /// rather than automatically executing the query. Full automatic region queries
    /// with seeking will be added in Phase 2.
    ///
    /// # Arguments
    ///
    /// * `index` - BAI index for this BAM file
    /// * `reference_name` - Reference sequence name (e.g., "chr1", "1")
    /// * `start` - Start position (0-based, inclusive)
    /// * `end` - End position (0-based, exclusive)
    ///
    /// # Returns
    ///
    /// * `Some(chunks)` - List of chunks to read for this region
    /// * `None` - Reference not found or no data for region
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use biometal::io::bam::{BamReader, BaiIndex};
    /// # fn main() -> biometal::Result<()> {
    /// let bam = BamReader::from_path("alignments.bam")?;
    /// let index = BaiIndex::from_path("alignments.bam.bai")?;
    ///
    /// // Get chunks for region chr1:1000-2000
    /// if let Some(chunks) = bam.query_chunks(&index, "chr1", 1000, 2000) {
    ///     println!("Need to read {} chunks for this region", chunks.len());
    ///     for chunk in &chunks {
    ///         println!("  Chunk: compressed offset {} to {}",
    ///                  chunk.start.compressed_offset(),
    ///                  chunk.end.compressed_offset());
    ///     }
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn query_chunks(
        &self,
        index: &BaiIndex,
        reference_name: &str,
        start: i32,
        end: i32,
    ) -> Option<Vec<Chunk>> {
        // Look up reference ID from name
        let ref_id = self.header.reference_id(reference_name)?;

        // Query index for chunks
        index.query_chunks(ref_id, start, end)
    }

    /// Plan a region query using reference ID.
    ///
    /// Similar to `query_chunks` but takes a reference ID directly instead of name.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use biometal::io::bam::{BamReader, BaiIndex};
    /// # fn main() -> biometal::Result<()> {
    /// let bam = BamReader::from_path("alignments.bam")?;
    /// let index = BaiIndex::from_path("alignments.bam.bai")?;
    ///
    /// // Query first reference (ID 0), region 1000-2000
    /// if let Some(chunks) = bam.query_chunks_by_id(&index, 0, 1000, 2000) {
    ///     println!("Found {} chunks", chunks.len());
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn query_chunks_by_id(
        &self,
        index: &BaiIndex,
        ref_id: usize,
        start: i32,
        end: i32,
    ) -> Option<Vec<Chunk>> {
        index.query_chunks(ref_id, start, end)
    }

    /// Execute a region query and return an iterator over overlapping records.
    ///
    /// **Phase 2**: This is the complete region query implementation with automatic
    /// seeking and filtering. Records are streamed with constant memory.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the BAM file (must be seekable, local file)
    /// * `index` - BAI index for this BAM file
    /// * `reference_name` - Reference sequence name (e.g., "chr1", "1")
    /// * `start` - Start position (0-based, inclusive)
    /// * `end` - End position (0-based, exclusive)
    ///
    /// # Returns
    ///
    /// An iterator that yields records overlapping the specified region.
    /// Records are filtered to only include those that overlap [start, end).
    ///
    /// # Memory Footprint (Rule 5)
    ///
    /// - Seekable reader: ~130 KB (one decompressed block)
    /// - Record buffer: ~500 bytes (reused)
    /// - Total: Constant memory regardless of region size
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use biometal::io::bam::{BamReader, BaiIndex};
    /// # fn main() -> biometal::Result<()> {
    /// let index = BaiIndex::from_path("alignments.bam.bai")?;
    ///
    /// // Query chr1:1000-2000
    /// let query = BamReader::query(
    ///     "alignments.bam",
    ///     &index,
    ///     "chr1",
    ///     1000,
    ///     2000
    /// )?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     println!("Read {} at position {}",
    ///              record.name,
    ///              record.position.unwrap_or(-1));
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn query<P: AsRef<Path>>(
        path: P,
        index: &BaiIndex,
        reference_name: &str,
        start: i32,
        end: i32,
    ) -> crate::Result<RegionQuery> {
        // Open BAM file to get header and reference ID
        let bam = BamReader::from_path(&path)?;

        // Look up reference ID from name
        let ref_id = bam.header.reference_id(reference_name)
            .ok_or_else(|| io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Reference '{}' not found in BAM header", reference_name)
            ))?;

        // Get chunks for this region
        let chunks = index.query_chunks(ref_id, start, end)
            .ok_or_else(|| io::Error::new(
                io::ErrorKind::NotFound,
                format!("No data found for reference {} region {}:{}",
                         reference_name, start, end)
            ))?;

        // Open seekable reader
        let file = File::open(path.as_ref())?;
        let mut reader = crate::io::compression::SeekableBgzfReader::new(file)?;

        // Seek to first chunk if chunks exist and get its end offset
        let current_chunk_end = if let Some(first_chunk) = chunks.first() {
            reader.seek_to_virtual_offset(first_chunk.start.as_raw())?;
            first_chunk.end.as_raw()
        } else {
            0
        };

        Ok(RegionQuery {
            reader,
            header: bam.header,
            chunks,
            current_chunk_idx: 0,
            ref_id,
            start,
            end,
            buffer: Vec::with_capacity(512),
            current_chunk_end,
        })
    }
}

/// Iterator over BAM records.
///
/// Created by [`BamReader::records()`]. Streams records with constant memory.
pub struct Records<'a, R> {
    reader: &'a mut BamReader<R>,
}

impl<'a, R: BufRead> Iterator for Records<'a, R> {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// Iterator over BAM records in a specific genomic region.
///
/// Created by [`BamReader::query()`]. Seeks through indexed chunks and yields
/// only records that overlap the query region.
///
/// # Memory Footprint (Rule 5)
///
/// - Seekable reader: ~130 KB (one decompressed block, cached)
/// - Record buffer: ~500 bytes (reused)
/// - Chunks vector: Typically <1 KB (small number of chunks per region)
/// - Total: Constant memory regardless of region size or file size
///
/// # Filtering
///
/// Records are filtered to ensure they overlap [start, end):
/// - Reference ID must match query reference
/// - Record's aligned region must overlap query region
/// - Unmapped records are excluded
pub struct RegionQuery {
    /// Seekable BGZF reader
    reader: crate::io::compression::SeekableBgzfReader,
    /// BAM header (for reference validation)
    header: Header,
    /// Chunks to read (from BAI index query)
    chunks: Vec<Chunk>,
    /// Current chunk index
    current_chunk_idx: usize,
    /// Query reference ID
    ref_id: usize,
    /// Query start position (0-based, inclusive)
    start: i32,
    /// Query end position (0-based, exclusive)
    end: i32,
    /// Reusable buffer for reading records
    buffer: Vec<u8>,
    /// Virtual offset marking end of current chunk
    current_chunk_end: u64,
}

impl RegionQuery {
    /// Check if a record overlaps the query region.
    ///
    /// A record overlaps if:
    /// - Reference ID matches
    /// - Record is mapped
    /// - Aligned region [pos, pos + alignment_length) overlaps [start, end)
    fn record_overlaps(&self, record: &Record) -> bool {
        // Check reference ID
        if record.reference_id != Some(self.ref_id) {
            return false;
        }

        // Check if mapped
        let pos = match record.position {
            Some(p) => p,
            None => return false, // Unmapped
        };

        // Calculate alignment end position from CIGAR
        // For region overlap, we need: record_start < query_end && record_end > query_start
        let alignment_length = record.cigar.iter()
            .map(|op| op.reference_length())
            .sum::<i32>();

        let record_end = pos + alignment_length;

        // Check overlap: record must start before query end and end after query start
        pos < self.end && record_end > self.start
    }

    /// Read next record from current chunk.
    ///
    /// Returns Ok(Some(record)) if record read successfully,
    /// Ok(None) if end of chunk reached,
    /// Err if read/parse error.
    fn read_next_record(&mut self) -> io::Result<Option<Record>> {
        // Read block size (4 bytes, little-endian)
        let mut size_buf = [0u8; 4];
        match self.reader.read(&mut size_buf)? {
            0 => return Ok(None), // EOF or end of chunk
            n if n < 4 => return Ok(None), // Incomplete read, end of chunk
            _ => {}
        }

        let block_size = i32::from_le_bytes(size_buf);
        if block_size < 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Invalid block size: {}", block_size),
            ));
        }

        let block_size = block_size as usize;

        // Reuse buffer for record data
        self.buffer.clear();
        self.buffer.reserve(block_size + 4);

        // Put block size at start
        self.buffer.extend_from_slice(&size_buf);

        // Resize to hold block data
        self.buffer.resize(4 + block_size, 0);

        // Read record data
        self.reader.read_exact(&mut self.buffer[4..])?;

        // Parse record
        let record = parse_record(&self.buffer)?;
        Ok(Some(record))
    }
}

impl Iterator for RegionQuery {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // If we've processed all chunks, we're done
            if self.current_chunk_idx >= self.chunks.len() {
                return None;
            }

            // Check if we've reached the end of the current chunk
            if self.reader.virtual_offset() >= self.current_chunk_end {
                // Move to next chunk
                self.current_chunk_idx += 1;

                if self.current_chunk_idx >= self.chunks.len() {
                    return None; // No more chunks
                }

                // Seek to next chunk
                let next_chunk = &self.chunks[self.current_chunk_idx];
                match self.reader.seek_to_virtual_offset(next_chunk.start.as_raw()) {
                    Ok(_) => {
                        self.current_chunk_end = next_chunk.end.as_raw();
                    }
                    Err(e) => return Some(Err(e)),
                }
            }

            // Try to read next record
            match self.read_next_record() {
                Ok(Some(record)) => {
                    // Check if record overlaps query region
                    if self.record_overlaps(&record) {
                        return Some(Ok(record));
                    }
                    // Record doesn't overlap, try next
                    continue;
                }
                Ok(None) => {
                    // EOF reached, move to next chunk
                    self.current_chunk_idx += 1;

                    if self.current_chunk_idx >= self.chunks.len() {
                        return None; // No more chunks
                    }

                    // Seek to next chunk
                    let next_chunk = &self.chunks[self.current_chunk_idx];
                    match self.reader.seek_to_virtual_offset(next_chunk.start.as_raw()) {
                        Ok(_) => {
                            self.current_chunk_end = next_chunk.end.as_raw();
                        }
                        Err(e) => return Some(Err(e)),
                    }
                    continue;
                }
                Err(e) => {
                    return Some(Err(e));
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn create_minimal_bam() -> Vec<u8> {
        let mut data = Vec::new();

        // Magic
        data.extend_from_slice(b"BAM\x01");

        // SAM header text (empty)
        data.extend_from_slice(&0i32.to_le_bytes());

        // References (0)
        data.extend_from_slice(&0i32.to_le_bytes());

        // One minimal record
        let mut record_data = Vec::new();

        // Block size placeholder
        let block_size_pos = record_data.len();
        record_data.extend_from_slice(&0i32.to_le_bytes());

        // Reference ID (-1)
        record_data.extend_from_slice(&(-1i32).to_le_bytes());

        // Position (-1)
        record_data.extend_from_slice(&(-1i32).to_le_bytes());

        // l_read_name (5)
        record_data.push(5);

        // MAPQ (255)
        record_data.push(255);

        // BAI bin (0)
        record_data.extend_from_slice(&0u16.to_le_bytes());

        // n_cigar_op (0)
        record_data.extend_from_slice(&0u16.to_le_bytes());

        // FLAGS (4 = unmapped)
        record_data.extend_from_slice(&4u16.to_le_bytes());

        // l_seq (0)
        record_data.extend_from_slice(&0i32.to_le_bytes());

        // next_refID (-1)
        record_data.extend_from_slice(&(-1i32).to_le_bytes());

        // next_pos (-1)
        record_data.extend_from_slice(&(-1i32).to_le_bytes());

        // tlen (0)
        record_data.extend_from_slice(&0i32.to_le_bytes());

        // read_name ("read\0")
        record_data.extend_from_slice(b"read\0");

        // Update block size
        let block_size = (record_data.len() - 4) as i32;
        record_data[block_size_pos..block_size_pos + 4]
            .copy_from_slice(&block_size.to_le_bytes());

        data.extend_from_slice(&record_data);

        data
    }

    #[test]
    fn test_bam_reader_new() {
        let bam_data = create_minimal_bam();
        let cursor = Cursor::new(bam_data);
        let bam = BamReader::new(cursor).unwrap();

        assert_eq!(bam.header().reference_count(), 0);
        assert_eq!(bam.header().text, "");
    }

    #[test]
    fn test_bam_reader_read_record() {
        let bam_data = create_minimal_bam();
        let cursor = Cursor::new(bam_data);
        let mut bam = BamReader::new(cursor).unwrap();

        // Read first record
        let record = bam.read_record().unwrap();
        assert!(record.is_some());
        let record = record.unwrap();
        assert_eq!(record.name, "read");

        // EOF
        let record = bam.read_record().unwrap();
        assert!(record.is_none());
    }

    #[test]
    fn test_bam_reader_records_iterator() {
        let bam_data = create_minimal_bam();
        let cursor = Cursor::new(bam_data);
        let mut bam = BamReader::new(cursor).unwrap();

        let records: Vec<_> = bam.records().collect::<io::Result<Vec<_>>>().unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, "read");
    }

    #[test]
    fn test_invalid_magic() {
        let data = b"INVALID";
        let cursor = Cursor::new(data);
        let result = BamReader::new(cursor);
        assert!(result.is_err());
    }
}
