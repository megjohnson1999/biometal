#!/usr/bin/env python3
"""
Test SAM reader functionality.
"""

import biometal
import tempfile
import os


def test_sam_reader_basic():
    """Test SamReader basic functionality"""
    print("Testing SamReader...")

    # Create a temporary SAM file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.sam', delete=False) as f:
        temp_path = f.name
        f.write("@HD\tVN:1.6\tSO:unsorted\n")
        f.write("@SQ\tSN:chr1\tLN:1000\n")
        f.write("@SQ\tSN:chr2\tLN:2000\n")
        f.write("read1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\t**********\n")
        f.write("read2\t16\tchr2\t200\t30\t5M2I3M\t*\t0\t0\tACGTACGTAC\t**********\n")
        f.write("read3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t****\n")

    try:
        # Read SAM file
        sam = biometal.SamReader.from_path(temp_path)

        # Check header
        assert sam.reference_count == 2
        header = sam.header
        assert header.reference_count == 2
        assert header.reference_name(0) == "chr1"
        assert header.reference_name(1) == "chr2"
        assert header.reference_length(0) == 1000
        assert header.reference_length(1) == 2000

        # Check records
        records = list(sam)
        assert len(records) == 3

        # Check first record (mapped)
        assert records[0].name == "read1"
        assert records[0].reference_id == 0  # chr1
        assert records[0].position == 99  # 0-based (SAM is 1-based)
        assert records[0].mapq == 60
        assert records[0].is_mapped
        assert records[0].is_forward
        assert records[0].sequence_str == "ACGTACGTAC"
        assert len(records[0].cigar) == 1
        assert records[0].cigar[0].length == 10
        assert records[0].cigar[0].is_match()

        # Check second record (reverse strand, with insertion)
        assert records[1].name == "read2"
        assert records[1].reference_id == 1  # chr2
        assert records[1].position == 199  # 0-based
        assert records[1].is_reverse
        assert len(records[1].cigar) == 3
        assert records[1].cigar[0].is_match()
        assert records[1].cigar[1].is_insertion()
        assert records[1].cigar[2].is_match()

        # Check third record (unmapped)
        assert records[2].name == "read3"
        assert records[2].reference_id is None
        assert records[2].position is None
        assert not records[2].is_mapped

        print("  ✓ SAM reader basic test passed!")

    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


def test_sam_reader_empty_header():
    """Test SamReader with minimal header"""
    print("\nTesting SamReader with empty header...")

    with tempfile.NamedTemporaryFile(mode='w', suffix='.sam', delete=False) as f:
        temp_path = f.name
        f.write("@HD\tVN:1.6\n")
        f.write("read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t****\n")

    try:
        sam = biometal.SamReader.from_path(temp_path)

        # No references
        assert sam.reference_count == 0

        # One unmapped record
        records = list(sam)
        assert len(records) == 1
        assert records[0].name == "read1"

        print("  ✓ SAM empty header test passed!")

    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


if __name__ == '__main__':
    print("SAM Reader Tests")
    print("=" * 60)

    try:
        test_sam_reader_basic()
        test_sam_reader_empty_header()

        print("\n" + "=" * 60)
        print("✓ All SAM reader tests passed!")
        print("=" * 60)

    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
