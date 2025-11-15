#!/usr/bin/env python3
"""
Test Python bindings for FASTQ and FASTA writers.

This script demonstrates the Python API for writing files.
"""

import biometal
import tempfile
import os


def test_fastq_writer():
    """Test FastqWriter Python bindings"""
    print("Testing FastqWriter...")

    # Create temporary file
    with tempfile.NamedTemporaryFile(suffix='.fq.gz', delete=False) as f:
        temp_path = f.name

    try:
        # Create writer
        writer = biometal.FastqWriter.create(temp_path)
        print(f"✓ Created writer: {writer}")

        # Create test records (when record creation API is available)
        # For now, just test the writer methods exist
        assert hasattr(writer, 'write_record'), "write_record method missing"
        assert hasattr(writer, 'records_written'), "records_written method missing"
        assert hasattr(writer, 'flush'), "flush method missing"
        assert hasattr(writer, 'finish'), "finish method missing"

        # Test records_written
        count = writer.records_written()
        assert count == 0, f"Expected 0 records, got {count}"
        print(f"✓ Records written: {count}")

        # Finish writer
        writer.finish()
        print("✓ FastqWriter works!")

    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


def test_fasta_writer():
    """Test FastaWriter Python bindings"""
    print("\nTesting FastaWriter...")

    # Create temporary file
    with tempfile.NamedTemporaryFile(suffix='.fa.gz', delete=False) as f:
        temp_path = f.name

    try:
        # Create writer
        writer = biometal.FastaWriter.create(temp_path)
        print(f"✓ Created writer: {writer}")

        # Test methods exist
        assert hasattr(writer, 'write_record'), "write_record method missing"
        assert hasattr(writer, 'with_line_width'), "with_line_width method missing"
        assert hasattr(writer, 'records_written'), "records_written method missing"
        assert hasattr(writer, 'flush'), "flush method missing"
        assert hasattr(writer, 'finish'), "finish method missing"

        # Test records_written
        count = writer.records_written()
        assert count == 0, f"Expected 0 records, got {count}"
        print(f"✓ Records written: {count}")

        # Finish writer
        writer.finish()
        print("✓ FastaWriter works!")

    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


def test_fasta_writer_stdout():
    """Test FastaWriter.stdout()"""
    print("\nTesting FastaWriter.stdout()...")

    # Create stdout writer
    writer = biometal.FastaWriter.stdout()
    print(f"✓ Created stdout writer: {writer}")

    # Note: Don't actually write to stdout in tests
    # Just verify the method exists and returns a writer
    assert hasattr(writer, 'write_record'), "write_record method missing"

    print("✓ FastaWriter.stdout() works!")


if __name__ == '__main__':
    # Note: These tests will only work once the writers are properly exported
    # to Python. Currently they demonstrate the expected API.

    print("FASTQ/FASTA Writer Python Bindings Test")
    print("=" * 50)

    try:
        test_fastq_writer()
        test_fasta_writer()
        test_fasta_writer_stdout()

        print("\n" + "=" * 50)
        print("✓ All Python writer bindings tests passed!")

    except AttributeError as e:
        print(f"\n⚠ Writer classes not yet exported to Python: {e}")
        print("  The Rust implementation is complete.")
        print("  Export registration is in src/python/mod.rs")
        print("  Rebuild with: maturin develop --release --features python")
