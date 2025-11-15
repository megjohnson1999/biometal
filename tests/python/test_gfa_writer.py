#!/usr/bin/env python3
"""
Test GFA writer functionality.
"""

import biometal
import tempfile
import os


def test_gfa_writer_basic():
    """Test GfaWriter basic functionality"""
    print("Testing GfaWriter...")

    with tempfile.NamedTempFile(suffix='.gfa', delete=False) as f:
        temp_path = f.name

    try:
        # Create writer
        writer = biometal.GfaWriter.create(temp_path)

        # Write header
        writer.write_header({"VN": "Z:1.0"})

        # Write segment
        segment = biometal.GfaSegment.from_line("S\tcontig1\tACGTACGT")
        writer.write_segment(segment)

        # Write link
        link = biometal.GfaLink.from_line("L\tcontig1\t+\tcontig2\t+\t4M")
        writer.write_link(link)

        # Write path
        path = biometal.GfaPath.from_line("P\tpath1\tcontig1+,contig2+\t4M")
        writer.write_path(path)

        assert writer.records_written() == 4
        writer.finish()

        # Read back and verify
        with open(temp_path) as f:
            lines = f.read().strip().split('\n')

        assert len(lines) == 4
        assert lines[0].startswith('H')
        assert lines[1].startswith('S')
        assert lines[2].startswith('L')
        assert lines[3].startswith('P')

        print("  ✓ Basic writer test passed")

    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


def test_gfa_writer_round_trip():
    """Test GFA writer round-trip"""
    print("\nTesting GFA writer round-trip...")

    with tempfile.NamedTempFile(suffix='.gfa', delete=False) as f:
        temp_path = f.name

    try:
        # Write data
        writer = biometal.GfaWriter.create(temp_path)

        segment1 = biometal.GfaSegment.from_line("S\tctg1\tACGT")
        segment2 = biometal.GfaSegment.from_line("S\tctg2\tTGCA")
        link = biometal.GfaLink.from_line("L\tctg1\t+\tctg2\t-\t4M")

        writer.write_segment(segment1)
        writer.write_segment(segment2)
        writer.write_link(link)
        writer.finish()

        # Read back
        stream = biometal.GfaStream.from_path(temp_path)
        records = list(stream)

        assert len(records) == 3

        # Verify records
        assert isinstance(records[0], biometal.GfaSegment)
        assert records[0].name == "ctg1"
        assert records[0].sequence == "ACGT"

        assert isinstance(records[1], biometal.GfaSegment)
        assert records[1].name == "ctg2"
        assert records[1].sequence == "TGCA"

        assert isinstance(records[2], biometal.GfaLink)
        assert records[2].from_segment == "ctg1"
        assert records[2].to_segment == "ctg2"

        print("  ✓ Round-trip test passed")

    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


def test_gfa_writer_compressed():
    """Test GFA writer with compression"""
    print("\nTesting GFA writer with compression...")

    with tempfile.NamedTempFile(suffix='.gfa.gz', delete=False) as f:
        temp_path = f.name

    try:
        writer = biometal.GfaWriter.create(temp_path)

        segment = biometal.GfaSegment.from_line("S\tctg1\tACGT")
        writer.write_segment(segment)
        writer.finish()

        # Read back compressed file
        stream = biometal.GfaStream.from_path(temp_path)
        records = list(stream)

        assert len(records) == 1
        assert records[0].name == "ctg1"

        print("  ✓ Compressed writer test passed")

    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


if __name__ == '__main__':
    print("GFA Writer Tests")
    print("=" * 60)

    try:
        test_gfa_writer_basic()
        test_gfa_writer_round_trip()
        test_gfa_writer_compressed()

        print("\n" + "=" * 60)
        print("✓ All GFA writer tests passed!")
        print("=" * 60)

    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
