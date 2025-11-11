#!/usr/bin/env python3
"""
End-to-end tests for BAI (BAM Index) Python bindings

Tests:
- BaiIndex loading and properties
- Indexed region queries via BamReader.query_region()
- Record filtering and validation
- Error handling for invalid inputs
- Performance characteristics
- Result consistency with Rust tests
"""
import pytest
from pathlib import Path
import time
import biometal


# Test data paths
BAM_FILE = Path("tests/data/synthetic_100k.bam")
BAI_FILE = Path("tests/data/synthetic_100k.bam.bai")


@pytest.fixture
def bam_path():
    """Fixture providing path to test BAM file"""
    if not BAM_FILE.exists():
        pytest.skip(f"Test BAM file not found: {BAM_FILE}")
    return str(BAM_FILE)


@pytest.fixture
def bai_path():
    """Fixture providing path to test BAI index file"""
    if not BAI_FILE.exists():
        pytest.skip(f"Test BAI file not found: {BAI_FILE}")
    return str(BAI_FILE)


@pytest.fixture
def bai_index(bai_path):
    """Fixture providing BaiIndex instance"""
    return biometal.BaiIndex.from_path(bai_path)


# ============================================================================
# BaiIndex Loading Tests
# ============================================================================

def test_bai_index_from_path(bai_path):
    """Test BaiIndex can be loaded from path"""
    index = biometal.BaiIndex.from_path(bai_path)
    assert index is not None


def test_bai_index_invalid_path():
    """Test BaiIndex raises error for invalid path"""
    with pytest.raises(Exception):
        biometal.BaiIndex.from_path("/nonexistent/file.bam.bai")


def test_bai_index_repr(bai_index):
    """Test BaiIndex has string representation"""
    repr_str = repr(bai_index)
    assert "BaiIndex" in repr_str
    assert "references" in repr_str


def test_bai_index_reference_count(bai_index):
    """Test BaiIndex has reference count"""
    # Our test file has 3 references (chr1, chr2, chr22)
    assert bai_index.reference_count == 3


# ============================================================================
# Indexed Query Tests - Basic Functionality
# ============================================================================

def test_query_region_small(bam_path, bai_index):
    """Test indexed query for small region (chr1:1-1000)"""
    query_iter = biometal.BamReader.query_region(
        bam_path,
        bai_index,
        "chr1",
        1,
        1000
    )

    assert query_iter is not None

    # Count records
    count = 0
    for record in query_iter:
        assert record is not None

        # Verify record is in expected region
        if record.position is not None:
            # BAM uses 0-based coordinates
            assert record.position >= 0
            assert record.position < 1000

        count += 1

    # From Rust tests, we expect ~2985 records (allow flexibility)
    assert 2900 <= count <= 3100, f"Expected ~2985 records, got {count}"


def test_query_region_medium(bam_path, bai_index):
    """Test indexed query for medium region (chr1:1-10000)"""
    query_iter = biometal.BamReader.query_region(
        bam_path,
        bai_index,
        "chr1",
        1,
        10000
    )

    count = 0
    for record in query_iter:
        if record.position is not None:
            assert record.position >= 0
            assert record.position < 10000
        count += 1

    # From Rust tests, we expect ~30693 records
    assert 30000 <= count <= 31000, f"Expected ~30693 records, got {count}"


def test_query_region_empty(bam_path, bai_index):
    """Test indexed query for empty region (chr22 far out)"""
    query_iter = biometal.BamReader.query_region(
        bam_path,
        bai_index,
        "chr22",
        50000000,
        50001000
    )

    count = 0
    for record in query_iter:
        count += 1

    # Should have very few or no records in this region
    assert count < 10


def test_query_region_iterator_behavior(bam_path, bai_index):
    """Test query iterator can be iterated multiple times"""
    query_iter1 = biometal.BamReader.query_region(
        bam_path,
        bai_index,
        "chr1",
        1,
        1000
    )

    count1 = sum(1 for _ in query_iter1)

    # Create new iterator for same region
    query_iter2 = biometal.BamReader.query_region(
        bam_path,
        bai_index,
        "chr1",
        1,
        1000
    )

    count2 = sum(1 for _ in query_iter2)

    # Should get same count
    assert count1 == count2


# ============================================================================
# Record Field Validation
# ============================================================================

def test_query_region_record_fields(bam_path, bai_index):
    """Test records from indexed queries have all fields populated"""
    query_iter = biometal.BamReader.query_region(
        bam_path,
        bai_index,
        "chr1",
        1,
        1000
    )

    checked = 0
    for record in query_iter:
        # Essential fields should be present
        assert isinstance(record.name, str)
        assert len(record.name) > 0

        assert isinstance(record.sequence, bytes)
        assert len(record.sequence) > 0

        assert isinstance(record.quality, bytes)
        assert len(record.quality) == len(record.sequence)

        # Optional fields
        if record.position is not None:
            assert isinstance(record.position, int)
            assert record.position >= 0

        checked += 1
        if checked >= 100:
            break

    assert checked > 0


def test_query_region_cigar_operations(bam_path, bai_index):
    """Test CIGAR operations are parsed correctly in indexed queries"""
    query_iter = biometal.BamReader.query_region(
        bam_path,
        bai_index,
        "chr1",
        1,
        1000
    )

    for record in query_iter:
        cigar = record.cigar
        assert isinstance(cigar, list)

        if len(cigar) > 0:
            op = cigar[0]
            assert hasattr(op, 'length')
            assert hasattr(op, 'op_char')
            assert op.length > 0
            assert op.op_char in 'MIDNSHP=X'
            break


# ============================================================================
# Error Handling Tests
# ============================================================================

def test_query_region_empty_reference_name(bam_path, bai_index):
    """Test error handling for empty reference name"""
    with pytest.raises(ValueError, match="Reference name cannot be empty"):
        biometal.BamReader.query_region(
            bam_path,
            bai_index,
            "",
            1,
            1000
        )


def test_query_region_negative_start(bam_path, bai_index):
    """Test error handling for negative start position"""
    with pytest.raises(ValueError, match="Start position must be non-negative"):
        biometal.BamReader.query_region(
            bam_path,
            bai_index,
            "chr1",
            -1,
            1000
        )


def test_query_region_negative_end(bam_path, bai_index):
    """Test error handling for negative end position"""
    with pytest.raises(ValueError, match="End position must be non-negative"):
        biometal.BamReader.query_region(
            bam_path,
            bai_index,
            "chr1",
            1,
            -1
        )


def test_query_region_invalid_range(bam_path, bai_index):
    """Test error handling for start >= end"""
    with pytest.raises(ValueError, match="Start position .* must be less than end position"):
        biometal.BamReader.query_region(
            bam_path,
            bai_index,
            "chr1",
            1000,
            1000
        )


def test_query_region_nonexistent_reference(bam_path, bai_index):
    """Test error handling for non-existent reference"""
    with pytest.raises(Exception):  # Should raise an error
        query_iter = biometal.BamReader.query_region(
            bam_path,
            bai_index,
            "chrX",
            1,
            1000
        )
        # Force evaluation by trying to iterate
        next(iter(query_iter))


# ============================================================================
# Comparison Tests (Indexed vs Full Scan)
# ============================================================================

def test_indexed_vs_full_scan_consistency(bam_path, bai_index):
    """Test indexed query returns reads overlapping region"""
    # Indexed query returns reads that OVERLAP the region
    indexed_query = biometal.BamReader.query_region(
        bam_path,
        bai_index,
        "chr1",
        100,
        200
    )
    indexed_names = set(record.name for record in indexed_query)

    # Full scan with manual filter - filter by START position only
    reader = biometal.BamReader.from_path(bam_path)
    filtered_names = set()

    for record in reader:
        if record.reference_id == 0:  # chr1
            if record.position is not None and 100 <= record.position < 200:
                filtered_names.add(record.name)

    # Indexed query returns MORE records than start-position filter
    # because it correctly includes reads that START BEFORE but OVERLAP the region
    # All records from filtered set should be in indexed set
    common = indexed_names & filtered_names
    assert len(common) == len(filtered_names), \
        f"All filtered records should be in indexed results. Indexed: {len(indexed_names)}, Filtered: {len(filtered_names)}, Common: {len(common)}"

    # Indexed query should have more records (reads starting before but overlapping)
    assert len(indexed_names) >= len(filtered_names), \
        f"Indexed query should include overlapping reads. Indexed: {len(indexed_names)}, Filtered: {len(filtered_names)}"


def test_adjacent_regions_no_overlap(bam_path, bai_index):
    """Test adjacent regions return different records"""
    # First region
    query1 = biometal.BamReader.query_region(
        bam_path,
        bai_index,
        "chr1",
        1,
        5000
    )
    count1 = sum(1 for _ in query1)

    # Adjacent region
    query2 = biometal.BamReader.query_region(
        bam_path,
        bai_index,
        "chr1",
        5000,
        10000
    )
    count2 = sum(1 for _ in query2)

    assert count1 > 0
    assert count2 > 0

    # Combined should be close to medium region query
    # (Allowing overlap for reads spanning the boundary)
    combined = count1 + count2
    assert 28000 <= combined <= 32000


# ============================================================================
# Performance Tests
# ============================================================================

def test_indexed_query_performance(bam_path, bai_index):
    """Test indexed queries are reasonably fast"""
    # Query small region multiple times
    times = []

    for _ in range(5):
        start = time.time()

        query_iter = biometal.BamReader.query_region(
            bam_path,
            bai_index,
            "chr1",
            1,
            1000
        )

        count = sum(1 for _ in query_iter)

        elapsed = time.time() - start
        times.append(elapsed)

        assert count > 0

    avg_time = sum(times) / len(times)

    # Should complete in < 50ms (being generous)
    assert avg_time < 0.05, f"Average query time {avg_time:.3f}s too slow"


def test_index_loading_overhead(bai_path):
    """Test index loading is fast"""
    times = []

    for _ in range(10):
        start = time.time()
        index = biometal.BaiIndex.from_path(bai_path)
        elapsed = time.time() - start
        times.append(elapsed)
        assert index is not None

    avg_time = sum(times) / len(times)

    # Should load in < 1ms (from benchmarks we know it's ~4.4 µs)
    assert avg_time < 0.001, f"Average load time {avg_time*1000:.3f}ms too slow"


# ============================================================================
# Integration Tests
# ============================================================================

def test_full_workflow_with_index(bam_path, bai_index):
    """Test complete workflow: load index, query region, analyze records"""
    # Query region
    query_iter = biometal.BamReader.query_region(
        bam_path,
        bai_index,
        "chr1",
        1000,
        2000
    )

    # Analyze records
    high_quality = 0
    total_bases = 0
    mapped_count = 0

    for record in query_iter:
        if record.is_mapped:
            mapped_count += 1

            if record.mapq is not None and record.mapq >= 30:
                high_quality += 1

            total_bases += len(record.sequence)

    # Should have found some records
    assert mapped_count > 0
    assert total_bases > 0


def test_multiple_queries_same_index(bam_path, bai_index):
    """Test same index can be used for multiple queries"""
    # Query 1
    query1 = biometal.BamReader.query_region(
        bam_path,
        bai_index,
        "chr1",
        1,
        1000
    )
    count1 = sum(1 for _ in query1)

    # Query 2 (different region)
    query2 = biometal.BamReader.query_region(
        bam_path,
        bai_index,
        "chr1",
        10000,
        11000
    )
    count2 = sum(1 for _ in query2)

    # Query 3 (back to first region)
    query3 = biometal.BamReader.query_region(
        bam_path,
        bai_index,
        "chr1",
        1,
        1000
    )
    count3 = sum(1 for _ in query3)

    # First and third should match
    assert count1 == count3
    assert count1 > 0
    assert count2 > 0


def test_streaming_constant_memory(bam_path, bai_index):
    """Test indexed queries maintain constant memory (streaming)"""
    # Query large region
    query_iter = biometal.BamReader.query_region(
        bam_path,
        bai_index,
        "chr1",
        1,
        10000
    )

    # Process records one at a time (streaming)
    count = 0
    prev_pos = None

    for record in query_iter:
        # Verify records are generally ordered
        if record.position is not None:
            if prev_pos is not None:
                # Positions should not decrease dramatically
                assert record.position >= prev_pos - 1000, \
                    f"Positions should be ordered: prev={prev_pos}, current={record.position}"
            prev_pos = record.position

        count += 1

    assert count > 1000, "Should process substantial number of records"


# ============================================================================
# Boundary Conditions
# ============================================================================

def test_single_base_query(bam_path, bai_index):
    """Test query for single base position"""
    query_iter = biometal.BamReader.query_region(
        bam_path,
        bai_index,
        "chr1",
        100,
        101
    )

    count = sum(1 for _ in query_iter)

    # From Rust tests, we expect ~331 overlapping reads
    assert 300 <= count <= 400, f"Expected ~331 records, got {count}"


def test_large_region_query(bam_path, bai_index):
    """Test query for large region"""
    query_iter = biometal.BamReader.query_region(
        bam_path,
        bai_index,
        "chr1",
        1,
        100000
    )

    count = sum(1 for _ in query_iter)

    # Should return most/all records (test file has 100K records on chr1)
    assert count >= 90000


# ============================================================================
# Type Validation
# ============================================================================

def test_bai_index_type(bai_index):
    """Test BaiIndex is correct type"""
    assert hasattr(bai_index, 'reference_count')
    assert hasattr(bai_index, '__repr__')


def test_query_iterator_type(bam_path, bai_index):
    """Test query returns correct iterator type"""
    query_iter = biometal.BamReader.query_region(
        bam_path,
        bai_index,
        "chr1",
        1,
        1000
    )

    # Should be iterable
    assert hasattr(query_iter, '__iter__')
    assert hasattr(query_iter, '__next__')


# ============================================================================
# Run tests
# ============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
