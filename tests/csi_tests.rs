//! Integration tests for CSI (Coordinate-Sorted Index) format
//!
//! These tests validate the CSI index parser against various scenarios
//! including edge cases and real-world usage patterns.

use biometal::formats::index::CsiIndex;
use biometal::io::bam::index::{Chunk, VirtualOffset};
use std::io::Write;

/// Helper to create a minimal valid CSI index for testing
fn create_minimal_csi() -> Vec<u8> {
    let mut data = Vec::new();

    // Magic
    data.extend_from_slice(b"CSI\x01");

    // min_shift (14 = 16kb)
    data.extend_from_slice(&14i32.to_le_bytes());

    // depth (5 levels)
    data.extend_from_slice(&5i32.to_le_bytes());

    // aux_size (0 = no aux data)
    data.extend_from_slice(&0i32.to_le_bytes());

    // n_ref (1 reference)
    data.extend_from_slice(&1i32.to_le_bytes());

    // Reference 0:
    // n_bin (1 bin)
    data.extend_from_slice(&1i32.to_le_bytes());

    // Bin 0 (entire sequence):
    data.extend_from_slice(&0u32.to_le_bytes()); // bin_id
    data.extend_from_slice(&0u64.to_le_bytes()); // loffset
    data.extend_from_slice(&1i32.to_le_bytes()); // n_chunk

    // Chunk:
    data.extend_from_slice(&1000u64.to_le_bytes()); // chunk_beg
    data.extend_from_slice(&2000u64.to_le_bytes()); // chunk_end

    // Note: CSI does not have n_intv like BAI/TBI

    data
}

/// Helper to create CSI with auxiliary reference names
fn create_csi_with_names() -> Vec<u8> {
    let mut data = Vec::new();

    // Magic
    data.extend_from_slice(b"CSI\x01");

    // min_shift (14)
    data.extend_from_slice(&14i32.to_le_bytes());

    // depth (5)
    data.extend_from_slice(&5i32.to_le_bytes());

    // aux_data with tab-delimited reference names
    let aux = b"chr1\tchr2\tchr3";
    data.extend_from_slice(&(aux.len() as i32).to_le_bytes());
    data.extend_from_slice(aux);

    // n_ref (3 references)
    data.extend_from_slice(&3i32.to_le_bytes());

    // References (empty for simplicity)
    for _ in 0..3 {
        data.extend_from_slice(&0i32.to_le_bytes()); // n_bin
        // Note: CSI does not have n_intv
    }

    data
}

/// Helper to create CSI with multiple bins
fn create_csi_with_multiple_bins() -> Vec<u8> {
    let mut data = Vec::new();

    // Magic
    data.extend_from_slice(b"CSI\x01");

    // min_shift (14)
    data.extend_from_slice(&14i32.to_le_bytes());

    // depth (5)
    data.extend_from_slice(&5i32.to_le_bytes());

    // aux_size (0)
    data.extend_from_slice(&0i32.to_le_bytes());

    // n_ref (1)
    data.extend_from_slice(&1i32.to_le_bytes());

    // Reference 0:
    // n_bin (3 bins)
    data.extend_from_slice(&3i32.to_le_bytes());

    // Bin 0:
    data.extend_from_slice(&0u32.to_le_bytes());
    data.extend_from_slice(&0u64.to_le_bytes());
    data.extend_from_slice(&1i32.to_le_bytes());
    data.extend_from_slice(&1000u64.to_le_bytes());
    data.extend_from_slice(&2000u64.to_le_bytes());

    // Bin 1:
    data.extend_from_slice(&1u32.to_le_bytes());
    data.extend_from_slice(&1000u64.to_le_bytes());
    data.extend_from_slice(&2i32.to_le_bytes());
    data.extend_from_slice(&1000u64.to_le_bytes());
    data.extend_from_slice(&1500u64.to_le_bytes());
    data.extend_from_slice(&1500u64.to_le_bytes());
    data.extend_from_slice(&2000u64.to_le_bytes());

    // Bin 2:
    data.extend_from_slice(&2u32.to_le_bytes());
    data.extend_from_slice(&2000u64.to_le_bytes());
    data.extend_from_slice(&1i32.to_le_bytes());
    data.extend_from_slice(&2000u64.to_le_bytes());
    data.extend_from_slice(&3000u64.to_le_bytes());

    // Note: CSI does not have n_intv

    data
}

#[test]
fn test_csi_parse_minimal() {
    let data = create_minimal_csi();
    let mut cursor = std::io::Cursor::new(data);

    let index = CsiIndex::parse(&mut cursor).expect("Failed to parse CSI");

    assert_eq!(index.min_shift(), 14);
    assert_eq!(index.depth(), 5);
    assert_eq!(index.aux_data().len(), 0);
    assert_eq!(index.references().len(), 1);

    let ref0 = &index.references()[0];
    assert_eq!(ref0.bins.len(), 1);
    // Note: CSI does not have intervals

    let bin0 = &ref0.bins[0];
    assert_eq!(bin0.bin_id, 0);
    assert_eq!(bin0.loffset.as_raw(), 0);
    assert_eq!(bin0.chunks.len(), 1);
    assert_eq!(bin0.chunks[0].start.as_raw(), 1000);
    assert_eq!(bin0.chunks[0].end.as_raw(), 2000);
}

#[test]
fn test_csi_parse_with_names() {
    let data = create_csi_with_names();
    let mut cursor = std::io::Cursor::new(data);

    let index = CsiIndex::parse(&mut cursor).expect("Failed to parse CSI");

    assert_eq!(index.references().len(), 3);
    assert!(index.aux_data().len() > 0);

    // Check reference names were parsed
    assert!(index.get_reference("chr1").is_some());
    assert!(index.get_reference("chr2").is_some());
    assert!(index.get_reference("chr3").is_some());
    assert!(index.get_reference("chr4").is_none());
}

// Note: Test removed - CSI does not have linear intervals like BAI/TBI

#[test]
fn test_csi_parse_with_multiple_bins() {
    let data = create_csi_with_multiple_bins();
    let mut cursor = std::io::Cursor::new(data);

    let index = CsiIndex::parse(&mut cursor).expect("Failed to parse CSI");

    let ref0 = &index.references()[0];
    assert_eq!(ref0.bins.len(), 3);

    // Bin 0
    assert_eq!(ref0.bins[0].bin_id, 0);
    assert_eq!(ref0.bins[0].chunks.len(), 1);

    // Bin 1 (2 chunks)
    assert_eq!(ref0.bins[1].bin_id, 1);
    assert_eq!(ref0.bins[1].chunks.len(), 2);
    assert_eq!(ref0.bins[1].loffset.as_raw(), 1000);

    // Bin 2
    assert_eq!(ref0.bins[2].bin_id, 2);
    assert_eq!(ref0.bins[2].chunks.len(), 1);
    assert_eq!(ref0.bins[2].loffset.as_raw(), 2000);
}

#[test]
fn test_csi_invalid_magic() {
    let mut data = create_minimal_csi();
    data[0] = b'X'; // Corrupt magic

    let mut cursor = std::io::Cursor::new(data);
    let result = CsiIndex::parse(&mut cursor);

    assert!(result.is_err());
    assert!(result
        .unwrap_err()
        .to_string()
        .contains("Invalid CSI magic"));
}

#[test]
fn test_csi_query_by_index() {
    let data = create_csi_with_multiple_bins();
    let mut cursor = std::io::Cursor::new(data);
    let index = CsiIndex::parse(&mut cursor).expect("Failed to parse CSI");

    // Query region
    let chunks = index
        .query_by_index(0, 1000, 2000)
        .expect("Query failed")
        .expect("No chunks");

    assert!(!chunks.is_empty());
}

#[test]
fn test_csi_query_by_name() {
    let data = create_csi_with_names();
    let mut cursor = std::io::Cursor::new(data);
    let index = CsiIndex::parse(&mut cursor).expect("Failed to parse CSI");

    // Query existing reference
    let result = index.query("chr1", 1000, 2000);
    assert!(result.is_ok());

    // Query non-existent reference
    let result = index.query("chr99", 1000, 2000);
    assert!(result.unwrap().is_none());
}

#[test]
fn test_csi_query_invalid_range() {
    let data = create_minimal_csi();
    let mut cursor = std::io::Cursor::new(data);
    let index = CsiIndex::parse(&mut cursor).expect("Failed to parse CSI");

    // Invalid range: start >= end
    let result = index.query_by_index(0, 2000, 1000);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("Invalid range"));
}

#[test]
fn test_csi_different_min_shift() {
    let mut data = Vec::new();

    // Create CSI with min_shift=18 (256kb bins instead of 16kb)
    data.extend_from_slice(b"CSI\x01");
    data.extend_from_slice(&18i32.to_le_bytes()); // min_shift=18
    data.extend_from_slice(&5i32.to_le_bytes());
    data.extend_from_slice(&0i32.to_le_bytes());
    data.extend_from_slice(&1i32.to_le_bytes());
    data.extend_from_slice(&0i32.to_le_bytes()); // n_bin
    data.extend_from_slice(&0i32.to_le_bytes()); // n_intv

    let mut cursor = std::io::Cursor::new(data);
    let index = CsiIndex::parse(&mut cursor).expect("Failed to parse CSI");

    assert_eq!(index.min_shift(), 18);
}

#[test]
fn test_csi_different_depth() {
    let mut data = Vec::new();

    // Create CSI with depth=3 (fewer binning levels)
    data.extend_from_slice(b"CSI\x01");
    data.extend_from_slice(&14i32.to_le_bytes());
    data.extend_from_slice(&3i32.to_le_bytes()); // depth=3
    data.extend_from_slice(&0i32.to_le_bytes());
    data.extend_from_slice(&1i32.to_le_bytes());
    data.extend_from_slice(&0i32.to_le_bytes()); // n_bin
    data.extend_from_slice(&0i32.to_le_bytes()); // n_intv

    let mut cursor = std::io::Cursor::new(data);
    let index = CsiIndex::parse(&mut cursor).expect("Failed to parse CSI");

    assert_eq!(index.depth(), 3);
}

#[test]
fn test_csi_get_reference_by_index() {
    let data = create_csi_with_names();
    let mut cursor = std::io::Cursor::new(data);
    let index = CsiIndex::parse(&mut cursor).expect("Failed to parse CSI");

    // Get by valid index
    let ref0 = index.get_reference_by_index(0);
    assert!(ref0.is_some());
    assert_eq!(ref0.unwrap().name.as_deref(), Some("chr1"));

    let ref1 = index.get_reference_by_index(1);
    assert!(ref1.is_some());
    assert_eq!(ref1.unwrap().name.as_deref(), Some("chr2"));

    // Get by invalid index
    let ref_invalid = index.get_reference_by_index(99);
    assert!(ref_invalid.is_none());
}

#[test]
fn test_csi_empty_references() {
    let mut data = Vec::new();

    // Create CSI with 0 references
    data.extend_from_slice(b"CSI\x01");
    data.extend_from_slice(&14i32.to_le_bytes());
    data.extend_from_slice(&5i32.to_le_bytes());
    data.extend_from_slice(&0i32.to_le_bytes());
    data.extend_from_slice(&0i32.to_le_bytes()); // n_ref=0

    let mut cursor = std::io::Cursor::new(data);
    let index = CsiIndex::parse(&mut cursor).expect("Failed to parse CSI");

    assert_eq!(index.references().len(), 0);
}

#[test]
fn test_csi_aux_data_preservation() {
    let data = create_csi_with_names();
    let mut cursor = std::io::Cursor::new(data);
    let index = CsiIndex::parse(&mut cursor).expect("Failed to parse CSI");

    // Check aux data is preserved
    let aux = index.aux_data();
    assert!(!aux.is_empty());
    assert!(std::str::from_utf8(aux).is_ok());
}

#[test]
fn test_csi_chunk_loffset_tracking() {
    let data = create_csi_with_multiple_bins();
    let mut cursor = std::io::Cursor::new(data);
    let index = CsiIndex::parse(&mut cursor).expect("Failed to parse CSI");

    let ref0 = &index.references()[0];

    // Check loffsets are correctly tracked
    assert_eq!(ref0.bins[0].loffset.as_raw(), 0);
    assert_eq!(ref0.bins[1].loffset.as_raw(), 1000);
    assert_eq!(ref0.bins[2].loffset.as_raw(), 2000);
}
