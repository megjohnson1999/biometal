"""Test BAM writer Python bindings."""

import biometal
import tempfile
import os


def test_bam_writer_basic():
    """Test basic BAM writing from Python."""
    # Read input BAM
    reader = biometal.BamReader.from_path("tests/data/synthetic_100k.bam")
    header = reader.header

    # Create output BAM
    with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
        output_path = f.name

    try:
        # Create writer
        writer = biometal.BamWriter.create(output_path, header)

        # Write first 10 records
        count = 0
        for record in reader:
            writer.write_record(record)
            count += 1
            if count >= 10:
                break

        # Finish writing
        writer.finish()

        # Verify
        assert writer.records_written() == 10

        # Read back and verify
        reader2 = biometal.BamReader.from_path(output_path)
        read_count = 0
        for _ in reader2:
            read_count += 1

        assert read_count == 10, f"Expected 10 records, got {read_count}"

    finally:
        # Clean up
        if os.path.exists(output_path):
            os.remove(output_path)


def test_bam_writer_filtering():
    """Test BAM filtering workflow."""
    # Read input BAM
    reader = biometal.BamReader.from_path("tests/data/synthetic_100k.bam")
    header = reader.header

    # Create filtered BAM
    with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
        filtered_path = f.name

    try:
        writer = biometal.BamWriter.create(filtered_path, header)

        # Filter by MAPQ ≥ 30
        for record in reader:
            if record.mapq and record.mapq >= 30:
                writer.write_record(record)

        writer.finish()
        filtered_count = writer.records_written()

        # Verify all records meet threshold
        reader2 = biometal.BamReader.from_path(filtered_path)
        for record in reader2:
            assert record.mapq is None or record.mapq >= 30, \
                f"Found record with MAPQ {record.mapq} < 30"

        print(f"Filtered to {filtered_count} high-quality reads")

    finally:
        if os.path.exists(filtered_path):
            os.remove(filtered_path)


if __name__ == "__main__":
    # Check if BamWriter is available
    if hasattr(biometal, 'BamWriter'):
        print("✅ BamWriter is available")
        test_bam_writer_basic()
        test_bam_writer_filtering()
        print("✅ All tests passed")
    else:
        print("❌ BamWriter is not available in Python bindings")
        print("This is a known PyO3 registration issue.")
        print("See KNOWN_ISSUES.md for details.")
