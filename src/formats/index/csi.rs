//! CSI (Coordinate-Sorted Index) format support
//!
//! This module implements parsing and querying of CSI index files (.csi),
//! which are the successor to BAI/TBI indexes with support for larger
//! reference sequences and configurable binning parameters.
//!
//! # Format Specification
//!
//! CSI is an extended index format that overcomes BAI/TBI limitations:
//! - **BAI/TBI limitation**: 512 Mbp maximum reference size
//! - **CSI advantage**: Configurable binning supports larger references
//!
//! ## Header
//! - Magic: "CSI\1" (4 bytes)
//! - min_shift: Binning depth (int32, typically 14 = 16kb)
//! - depth: Number of binning levels (int32, typically 5)
//! - aux_size: Size of auxiliary data (int32)
//! - aux_data: Optional auxiliary metadata (aux_size bytes)
//! - n_ref: Number of reference sequences (int32)
//!
//! ## Index Data (per reference)
//! - Binning index: Hierarchical bins with chunks and loffset
//! - Linear index: Intervals at min_shift resolution
//!
//! # Binning Scheme
//!
//! CSI uses configurable binning with parameters:
//! - **min_shift**: Base bin size (e.g., 14 = 16kb bins)
//! - **depth**: Number of hierarchical levels
//!
//! Default (min_shift=14, depth=5):
//! - Level 0: 1 bin (entire sequence)
//! - Level 1: 8 bins
//! - Level 2: 64 bins
//! - Level 3: 512 bins
//! - Level 4: 4096 bins
//! - Level 5: 32768 bins (16kb each)
//!
//! # Virtual File Offsets
//!
//! BGZF virtual offsets (64-bit):
//! - High 48 bits: Compressed file offset
//! - Low 16 bits: Uncompressed offset within block
//!
//! # Example
//!
//! ```no_run
//! use biometal::formats::index::CsiIndex;
//!
//! # fn main() -> biometal::Result<()> {
//! // Load CSI index
//! let index = CsiIndex::from_path("alignments.bam.csi")?;
//!
//! // Get metadata
//! println!("Min shift: {}", index.min_shift());
//! println!("Depth: {}", index.depth());
//! println!("References: {}", index.references().len());
//!
//! // Query region (if reference names available)
//! if let Some(chunks) = index.query("chr1", 1000000, 2000000)? {
//!     println!("Found {} chunks", chunks.len());
//! }
//! # Ok(())
//! # }
//! ```

use crate::error::{BiometalError, Result};
use crate::io::bam::index::{Chunk, VirtualOffset};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

/// CSI file format magic string
const CSI_MAGIC: &[u8; 4] = b"CSI\x01";

/// A bin in the hierarchical binning index
#[derive(Debug, Clone)]
pub struct CsiBin {
    /// Bin number
    pub bin_id: u32,
    /// Left-most file offset (optimization hint)
    pub loffset: VirtualOffset,
    /// Chunks of data in this bin
    pub chunks: Vec<Chunk>,
}

impl CsiBin {
    /// Create a new bin
    pub fn new(bin_id: u32, loffset: VirtualOffset) -> Self {
        CsiBin {
            bin_id,
            loffset,
            chunks: Vec::new(),
        }
    }
}

/// Reference sequence index data
#[derive(Debug, Clone)]
pub struct CsiReference {
    /// Reference sequence name (if available from aux data)
    pub name: Option<String>,
    /// Bins for this reference (hierarchical spatial index)
    pub bins: Vec<CsiBin>,
    // Note: CSI format does NOT have linear intervals like BAI/TBI
}

impl CsiReference {
    /// Create a new reference
    pub fn new(name: Option<String>) -> Self {
        CsiReference {
            name,
            bins: Vec::new(),
        }
    }
}

/// CSI (Coordinate-Sorted Index)
///
/// Provides fast random access to sorted genomic files with configurable
/// binning parameters that support larger reference sequences than BAI/TBI.
#[derive(Debug, Clone)]
pub struct CsiIndex {
    /// Binning depth (typically 14 = 16kb)
    min_shift: i32,
    /// Number of binning levels (typically 5)
    depth: i32,
    /// Auxiliary data (metadata, reference names, etc.)
    aux_data: Vec<u8>,
    /// Reference sequences
    references: Vec<CsiReference>,
    /// Reference name to index mapping (if names available)
    ref_map: HashMap<String, usize>,
    /// Optional unmapped read count
    n_no_coor: Option<u64>,
}

impl CsiIndex {
    /// Load CSI index from a file
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::index::CsiIndex;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let index = CsiIndex::from_path("alignments.bam.csi")?;
    /// println!("Loaded {} references", index.references().len());
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path.as_ref())?;
        let mut reader = BufReader::new(file);

        Self::parse(&mut reader)
    }

    /// Parse CSI index from a reader
    pub fn parse<R: Read>(reader: &mut R) -> Result<Self> {
        // Read and verify magic string
        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)?;
        if &magic != CSI_MAGIC {
            return Err(BiometalError::InvalidInput {
                msg: format!(
                    "Invalid CSI magic: expected {:?}, got {:?}",
                    CSI_MAGIC, magic
                ),
            });
        }

        // Read header fields
        let min_shift = read_i32(reader)?;
        let depth = read_i32(reader)?;
        let aux_size = read_i32(reader)?;

        if min_shift < 0 || depth < 0 || aux_size < 0 {
            return Err(BiometalError::InvalidInput {
                msg: format!(
                    "Invalid CSI header: min_shift={}, depth={}, aux_size={}",
                    min_shift, depth, aux_size
                ),
            });
        }

        // Read auxiliary data
        let mut aux_data = vec![0u8; aux_size as usize];
        if aux_size > 0 {
            reader.read_exact(&mut aux_data)?;
        }

        // Parse reference names from aux data if present
        let ref_names = parse_reference_names(&aux_data)?;

        // Read number of references
        let n_ref = read_i32(reader)?;
        if n_ref < 0 {
            return Err(BiometalError::InvalidInput {
                msg: format!("Invalid reference count: {}", n_ref),
            });
        }

        // Parse index data for each reference
        let mut references = Vec::new();
        let mut ref_map = HashMap::new();

        for idx in 0..n_ref as usize {
            let name = ref_names.get(idx).cloned();
            let mut reference = CsiReference::new(name.clone());

            // Read binning index
            let n_bin = read_i32(reader)?;
            if n_bin < 0 {
                return Err(BiometalError::InvalidInput {
                    msg: format!("Invalid bin count: {}", n_bin),
                });
            }

            // Calculate metadata bin ID
            // Formula from CSIv1 spec: bin_limit = (1 << ((depth + 1) * 3)) / 7
            // metadata_id = bin_limit + 1
            // For depth=5: (1 << 18) / 7 + 1 = 262144 / 7 + 1 = 37450
            let bin_limit = (1u32 << ((depth + 1) * 3)) / 7;
            let metadata_bin_id = bin_limit + 1;

            for _ in 0..n_bin {
                let bin_id = read_u32(reader)?;
                let loffset = VirtualOffset::from_raw(read_u64(reader)?);
                let n_chunk = read_i32(reader)?;

                if n_chunk < 0 {
                    return Err(BiometalError::InvalidInput {
                        msg: format!("Invalid chunk count: {}", n_chunk),
                    });
                }

                // Bin 37450 (for depth=5) is a special metadata bin
                // It contains ref_beg, ref_end, n_mapped, n_unmapped instead of chunks
                if bin_id == metadata_bin_id {
                    // Read metadata fields but don't create chunks
                    // n_chunk should be 2 (ref_beg+ref_end, n_mapped+n_unmapped)
                    for _ in 0..n_chunk {
                        let _field1 = read_u64(reader)?;  // ref_beg or n_mapped
                        let _field2 = read_u64(reader)?;  // ref_end or n_unmapped
                    }
                    // Don't add metadata bin to reference bins (we don't need it for queries)
                } else {
                    // Regular bin with chunks
                    let mut bin = CsiBin::new(bin_id, loffset);
                    for _ in 0..n_chunk {
                        let chunk_beg = read_u64(reader)?;
                        let chunk_end = read_u64(reader)?;
                        bin.chunks.push(Chunk::new(
                            VirtualOffset::from_raw(chunk_beg),
                            VirtualOffset::from_raw(chunk_end),
                        ));
                    }
                    reference.bins.push(bin);
                }
            }

            // Note: CSI format does NOT have linear intervals (n_intv) like BAI/TBI
            // Each reference only contains n_bin + bin data

            if let Some(ref name) = name {
                ref_map.insert(name.clone(), idx);
            }
            references.push(reference);
        }

        // Read optional n_no_coor (unmapped read count) - u64
        // This is optional; if we hit EOF, just continue without it
        let n_no_coor = match read_u64(reader) {
            Ok(n) => Some(n),
            Err(BiometalError::Io(ref e)) if e.kind() == std::io::ErrorKind::UnexpectedEof => None,
            Err(e) => return Err(e),
        };

        Ok(CsiIndex {
            min_shift,
            depth,
            aux_data,
            references,
            ref_map,
            n_no_coor,
        })
    }

    /// Get binning depth parameter
    pub fn min_shift(&self) -> i32 {
        self.min_shift
    }

    /// Get number of binning levels
    pub fn depth(&self) -> i32 {
        self.depth
    }

    /// Get auxiliary data
    pub fn aux_data(&self) -> &[u8] {
        &self.aux_data
    }

    /// Get all references
    pub fn references(&self) -> &[CsiReference] {
        &self.references
    }

    /// Get reference by name
    pub fn get_reference(&self, name: &str) -> Option<&CsiReference> {
        self.ref_map.get(name).map(|&idx| &self.references[idx])
    }

    /// Get reference by index
    pub fn get_reference_by_index(&self, idx: usize) -> Option<&CsiReference> {
        self.references.get(idx)
    }

    /// Query region by name and get chunks to read
    ///
    /// Returns None if reference name is not found in the index.
    ///
    /// # Arguments
    ///
    /// * `ref_name` - Reference sequence name
    /// * `start` - Start position (0-based, inclusive)
    /// * `end` - End position (0-based, exclusive)
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::index::CsiIndex;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let index = CsiIndex::from_path("data.bam.csi")?;
    /// if let Some(chunks) = index.query("chr1", 1000000, 2000000)? {
    ///     println!("Found {} chunks", chunks.len());
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn query(&self, ref_name: &str, start: u32, end: u32) -> Result<Option<Vec<Chunk>>> {
        let reference = match self.get_reference(ref_name) {
            Some(r) => r,
            None => return Ok(None),
        };

        self.query_reference(reference, start, end)
    }

    /// Query region by index and get chunks to read
    ///
    /// # Arguments
    ///
    /// * `ref_idx` - Reference sequence index (0-based)
    /// * `start` - Start position (0-based, inclusive)
    /// * `end` - End position (0-based, exclusive)
    pub fn query_by_index(&self, ref_idx: usize, start: u32, end: u32) -> Result<Option<Vec<Chunk>>> {
        let reference = match self.get_reference_by_index(ref_idx) {
            Some(r) => r,
            None => return Ok(None),
        };

        self.query_reference(reference, start, end)
    }

    /// Internal query implementation
    fn query_reference(&self, reference: &CsiReference, start: u32, end: u32) -> Result<Option<Vec<Chunk>>> {
        if start >= end {
            return Err(BiometalError::InvalidRange(format!(
                "Invalid range: start ({}) >= end ({})",
                start, end
            )));
        }

        // Get candidate bins that overlap [start, end)
        let bins = reg2bins_csi(start, end, self.min_shift, self.depth);

        // Collect chunks from overlapping bins
        let mut chunks = Vec::new();
        for bin_id in bins {
            if let Some(bin) = reference.bins.iter().find(|b| b.bin_id == bin_id) {
                chunks.extend_from_slice(&bin.chunks);
            }
        }

        // Note: CSI does not have linear intervals for optimization
        // We rely purely on the hierarchical binning index

        // Sort and merge overlapping chunks
        chunks.sort_by_key(|c| c.start.as_raw());
        let merged = merge_chunks(&chunks);

        Ok(Some(merged))
    }
}

/// Parse reference names from auxiliary data
///
/// CSI aux data typically contains tab-delimited reference names
fn parse_reference_names(aux: &[u8]) -> Result<Vec<String>> {
    if aux.is_empty() {
        return Ok(Vec::new());
    }

    // Try parsing as null-terminated strings
    let mut names = Vec::new();
    let mut start = 0;

    for (i, &byte) in aux.iter().enumerate() {
        if byte == 0 || byte == b'\t' || byte == b'\n' {
            if i > start {
                if let Ok(name) = std::str::from_utf8(&aux[start..i]) {
                    if !name.is_empty() {
                        names.push(name.to_string());
                    }
                }
            }
            start = i + 1;
        }
    }

    // Handle last name if no trailing delimiter
    if start < aux.len() {
        if let Ok(name) = std::str::from_utf8(&aux[start..]) {
            if !name.is_empty() {
                names.push(name.to_string());
            }
        }
    }

    Ok(names)
}

/// Calculate bin IDs that overlap a region with configurable binning
///
/// Uses CSI binning scheme with min_shift and depth parameters.
///
/// # Arguments
///
/// * `beg` - Start position (0-based, inclusive)
/// * `end` - End position (0-based, exclusive)
/// * `min_shift` - Base bin size (e.g., 14 = 16kb)
/// * `depth` - Number of binning levels
fn reg2bins_csi(beg: u32, end: u32, min_shift: i32, depth: i32) -> Vec<u32> {
    let mut bins = Vec::new();
    let end = end.saturating_sub(1); // Make end inclusive for calculation

    // Add bin 0 (entire sequence)
    bins.push(0);

    // Calculate bins for each level
    for level in 0..depth {
        let shift = min_shift + 3 * level;
        let offset = ((1u32 << (3 * (level + 1))) - 1) / 7;
        let beg_bin = offset + (beg >> shift);
        let end_bin = offset + (end >> shift);

        for bin in beg_bin..=end_bin {
            bins.push(bin);
        }
    }

    bins
}

/// Merge overlapping chunks
fn merge_chunks(chunks: &[Chunk]) -> Vec<Chunk> {
    if chunks.is_empty() {
        return Vec::new();
    }

    let mut merged = Vec::new();
    let mut current = chunks[0].clone();

    for chunk in &chunks[1..] {
        if chunk.start.as_raw() <= current.end.as_raw() {
            // Overlapping or adjacent - merge
            if chunk.end.as_raw() > current.end.as_raw() {
                current.end = chunk.end;
            }
        } else {
            // Non-overlapping - save current and start new
            merged.push(current.clone());
            current = chunk.clone();
        }
    }
    merged.push(current);

    merged
}

// Helper functions for reading binary data (little-endian)

fn read_i32<R: Read>(reader: &mut R) -> Result<i32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(i32::from_le_bytes(buf))
}

fn read_u32<R: Read>(reader: &mut R) -> Result<u32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_u64<R: Read>(reader: &mut R) -> Result<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reg2bins_csi_default() {
        // Default CSI parameters (min_shift=14, depth=5)
        let bins = reg2bins_csi(1000, 2000, 14, 5);
        assert!(bins.contains(&0)); // Level 0
        assert!(!bins.is_empty());
    }

    #[test]
    fn test_reg2bins_csi_large_shift() {
        // Larger min_shift for larger bin sizes
        let bins = reg2bins_csi(1000000, 2000000, 18, 5);
        assert!(bins.contains(&0));
        assert!(!bins.is_empty());
    }

    #[test]
    fn test_reg2bins_csi_single_position() {
        let bins = reg2bins_csi(1000, 1001, 14, 5);
        assert!(bins.contains(&0));
        // Should have bins from multiple levels
        assert!(bins.len() >= 2);
    }

    #[test]
    fn test_parse_reference_names_empty() {
        let names = parse_reference_names(&[]).unwrap();
        assert!(names.is_empty());
    }

    #[test]
    fn test_parse_reference_names_tab_delimited() {
        let data = b"chr1\tchr2\tchr3";
        let names = parse_reference_names(data).unwrap();
        assert_eq!(names, vec!["chr1", "chr2", "chr3"]);
    }

    #[test]
    fn test_parse_reference_names_null_terminated() {
        let data = b"chr1\0chr2\0chr3\0";
        let names = parse_reference_names(data).unwrap();
        assert_eq!(names, vec!["chr1", "chr2", "chr3"]);
    }

    #[test]
    fn test_parse_reference_names_newline_delimited() {
        let data = b"chr1\nchr2\nchr3\n";
        let names = parse_reference_names(data).unwrap();
        assert_eq!(names, vec!["chr1", "chr2", "chr3"]);
    }

    #[test]
    fn test_merge_chunks() {
        let chunks = vec![
            Chunk::new(VirtualOffset::from_raw(100), VirtualOffset::from_raw(200)),
            Chunk::new(VirtualOffset::from_raw(150), VirtualOffset::from_raw(250)),
            Chunk::new(VirtualOffset::from_raw(300), VirtualOffset::from_raw(400)),
        ];

        let merged = merge_chunks(&chunks);
        assert_eq!(merged.len(), 2);
        assert_eq!(merged[0].start.as_raw(), 100);
        assert_eq!(merged[0].end.as_raw(), 250);
        assert_eq!(merged[1].start.as_raw(), 300);
        assert_eq!(merged[1].end.as_raw(), 400);
    }

    #[test]
    fn test_merge_chunks_non_overlapping() {
        let chunks = vec![
            Chunk::new(VirtualOffset::from_raw(100), VirtualOffset::from_raw(200)),
            Chunk::new(VirtualOffset::from_raw(300), VirtualOffset::from_raw(400)),
        ];

        let merged = merge_chunks(&chunks);
        assert_eq!(merged.len(), 2);
    }

    #[test]
    fn test_merge_chunks_all_overlapping() {
        let chunks = vec![
            Chunk::new(VirtualOffset::from_raw(100), VirtualOffset::from_raw(300)),
            Chunk::new(VirtualOffset::from_raw(150), VirtualOffset::from_raw(350)),
            Chunk::new(VirtualOffset::from_raw(200), VirtualOffset::from_raw(400)),
        ];

        let merged = merge_chunks(&chunks);
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].start.as_raw(), 100);
        assert_eq!(merged[0].end.as_raw(), 400);
    }

    #[test]
    fn test_csi_bin_creation() {
        let bin = CsiBin::new(42, VirtualOffset::from_raw(1000));
        assert_eq!(bin.bin_id, 42);
        assert_eq!(bin.loffset.as_raw(), 1000);
        assert!(bin.chunks.is_empty());
    }
}
