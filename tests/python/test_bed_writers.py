#!/usr/bin/env python3
"""
Test BED writer Python bindings with round-trip verification.
"""

import biometal
import tempfile
import os


def test_bed3_writer_round_trip():
    """Test Bed3Writer round-trip"""
    print("Testing Bed3Writer...")

    with tempfile.NamedTemporaryFile(suffix='.bed.gz', delete=False) as f:
        temp_path = f.name

    try:
        # Write records
        writer = biometal.Bed3Writer.create(temp_path)

        records = [
            biometal.Bed3Record.from_line("chr1\t1000\t2000"),
            biometal.Bed3Record.from_line("chr2\t3000\t4000"),
            biometal.Bed3Record.from_line("chr3\t5000\t6000"),
        ]

        for record in records:
            writer.write_record(record)

        assert writer.records_written() == 3, f"Expected 3 records, got {writer.records_written()}"
        writer.finish()

        # Read back
        stream = biometal.Bed3Stream.from_path(temp_path)
        read_back = []
        for record in stream:
            read_back.append(record)

        # Verify
        assert len(read_back) == 3, f"Expected 3 records, got {len(read_back)}"
        for i, (orig, read) in enumerate(zip(records, read_back)):
            assert read.chrom == orig.chrom, f"Record {i}: chrom mismatch"
            assert read.start == orig.start, f"Record {i}: start mismatch"
            assert read.end == orig.end, f"Record {i}: end mismatch"

        print(f"  ✓ Wrote and verified {len(records)} BED3 records")

    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


def test_bed6_writer_round_trip():
    """Test Bed6Writer round-trip"""
    print("\nTesting Bed6Writer...")

    with tempfile.NamedTemporaryFile(suffix='.bed.gz', delete=False) as f:
        temp_path = f.name

    try:
        # Write records
        writer = biometal.Bed6Writer.create(temp_path)

        records = [
            biometal.Bed6Record.from_line("chr1\t1000\t2000\tgene1\t100\t+"),
            biometal.Bed6Record.from_line("chr2\t3000\t4000\tgene2\t200\t-"),
            biometal.Bed6Record.from_line("chr3\t5000\t6000\tgene3\t300\t."),
        ]

        for record in records:
            writer.write_record(record)

        assert writer.records_written() == 3
        writer.finish()

        # Read back
        stream = biometal.Bed6Stream.from_path(temp_path)
        read_back = list(stream)

        # Verify
        assert len(read_back) == 3
        for i, (orig, read) in enumerate(zip(records, read_back)):
            assert read.chrom == orig.chrom, f"Record {i}: chrom mismatch"
            assert read.start == orig.start, f"Record {i}: start mismatch"
            assert read.end == orig.end, f"Record {i}: end mismatch"
            assert read.name == orig.name, f"Record {i}: name mismatch"
            assert read.score == orig.score, f"Record {i}: score mismatch"
            assert read.strand == orig.strand, f"Record {i}: strand mismatch"

        print(f"  ✓ Wrote and verified {len(records)} BED6 records")

    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


def test_bed12_writer_round_trip():
    """Test Bed12Writer round-trip"""
    print("\nTesting Bed12Writer...")

    # Use a simpler BED12 example with required fields
    line = "chr1\t1000\t5000\tgene1\t500\t+\t1200\t4800\t255,0,0\t3\t500,500,500\t0,2000,4000"

    with tempfile.NamedTemporaryFile(suffix='.bed.gz', delete=False) as f:
        temp_path = f.name

    try:
        # Write record
        writer = biometal.Bed12Writer.create(temp_path)
        record = biometal.Bed12Record.from_line(line)
        writer.write_record(record)
        assert writer.records_written() == 1
        writer.finish()

        # Read back
        stream = biometal.Bed12Stream.from_path(temp_path)
        read_back = list(stream)

        # Verify
        assert len(read_back) == 1
        read = read_back[0]
        assert read.chrom == record.chrom
        assert read.start == record.start
        assert read.end == record.end
        assert read.name == record.name
        assert read.block_count == record.block_count

        print(f"  ✓ Wrote and verified 1 BED12 record with {record.block_count} blocks")

    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


def test_compression():
    """Test that .gz compression works"""
    print("\nTesting compression...")

    with tempfile.NamedTemporaryFile(suffix='.bed.gz', delete=False) as f:
        temp_path = f.name

    try:
        writer = biometal.Bed3Writer.create(temp_path)
        for i in range(100):
            record = biometal.Bed3Record.from_line(f"chr{i}\t{i*1000}\t{i*1000+500}")
            writer.write_record(record)
        writer.finish()

        # File should exist and be compressed
        assert os.path.exists(temp_path)
        file_size = os.path.getsize(temp_path)

        # Read back to verify compression worked
        stream = biometal.Bed3Stream.from_path(temp_path)
        count = sum(1 for _ in stream)
        assert count == 100

        print(f"  ✓ Wrote 100 records to compressed file ({file_size} bytes)")

    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


if __name__ == '__main__':
    print("BED Writer Python Bindings Round-Trip Tests")
    print("=" * 60)

    try:
        test_bed3_writer_round_trip()
        test_bed6_writer_round_trip()
        test_bed12_writer_round_trip()
        test_compression()

        print("\n" + "=" * 60)
        print("✓ All BED writer round-trip tests passed!")
        print("=" * 60)

    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
