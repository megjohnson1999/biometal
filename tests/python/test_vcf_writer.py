#!/usr/bin/env python3
"""
Test VCF writer round-trip functionality.
"""

import biometal
import tempfile
import os


def test_vcf_writer_basic():
    """Test VcfWriter basic functionality"""
    print("Testing VcfWriter...")

    with tempfile.NamedTemporaryFile(suffix='.vcf', delete=False) as f:
        temp_path = f.name

    try:
        # Create header
        header = biometal.VcfHeader("VCFv4.2")
        header.add_info("DP", "Total Depth")
        header.add_info("AF", "Allele Frequency")

        # Write VCF
        writer = biometal.VcfWriter.create(temp_path, header)

        line1 = "chr1\t12345\trs123\tA\tT\t30.0\tPASS\tDP=100;AF=0.5"
        record1 = biometal.VcfRecord.from_line(line1)
        writer.write_record(record1)

        line2 = "chr1\t12346\trs124\tG\tC\t40.0\tPASS\tDP=150;AF=0.75"
        record2 = biometal.VcfRecord.from_line(line2)
        writer.write_record(record2)

        assert writer.records_written() == 2
        writer.finish()

        # Read back
        stream = biometal.VcfStream.from_path(temp_path)
        read_header = stream.header()

        # Verify header
        assert read_header.fileformat == "VCFv4.2"
        assert "DP" in read_header.info_fields
        assert "AF" in read_header.info_fields

        # Verify records
        read_back = list(stream)
        assert len(read_back) == 2

        assert read_back[0].chrom == "chr1"
        assert read_back[0].pos == 12345
        assert read_back[0].id == "rs123"
        assert read_back[0].reference == "A"
        assert read_back[0].alternate == ["T"]

        assert read_back[1].chrom == "chr1"
        assert read_back[1].pos == 12346

        print("  ✓ VCF basic writer successful")

    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


def test_vcf_writer_compressed():
    """Test VcfWriter with gzip compression"""
    print("\nTesting VcfWriter with compression...")

    with tempfile.NamedTemporaryFile(suffix='.vcf.gz', delete=False) as f:
        temp_path = f.name

    try:
        # Create header with samples
        header = biometal.VcfHeader("VCFv4.2")
        header.add_info("DP", "Total Depth")
        header.samples = ["sample1", "sample2"]

        # Write compressed VCF
        writer = biometal.VcfWriter.create(temp_path, header)

        line = "chr1\t100\t.\tA\tT\t.\t.\tDP=50"
        record = biometal.VcfRecord.from_line(line)
        writer.write_record(record)
        writer.finish()

        # Read back compressed file
        stream = biometal.VcfStream.from_path(temp_path)
        read_header = stream.header()

        # Verify header includes samples
        assert read_header.fileformat == "VCFv4.2"
        assert read_header.samples == ["sample1", "sample2"]

        # Verify record
        read_back = list(stream)
        assert len(read_back) == 1
        assert read_back[0].chrom == "chr1"
        assert read_back[0].pos == 100

        print("  ✓ VCF compressed writer successful")

    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


def test_vcf_writer_with_samples():
    """Test VcfWriter with sample genotypes"""
    print("\nTesting VcfWriter with samples...")

    with tempfile.NamedTemporaryFile(suffix='.vcf', delete=False) as f:
        temp_path = f.name

    try:
        # Create header with samples
        header = biometal.VcfHeader("VCFv4.2")
        header.add_info("DP", "Total Depth")
        header.add_format("GT", "Genotype")
        header.add_format("DP", "Read Depth")
        header.samples = ["NA12878", "NA12879"]

        # Create VCF with genotypes
        writer = biometal.VcfWriter.create(temp_path, header)

        # Note: We're using from_line to create records, which should preserve samples
        line = "chr1\t100\t.\tA\tT\t30\tPASS\tDP=50\tGT:DP\t0/1:25\t1/1:25"
        record = biometal.VcfRecord.from_line(line)
        writer.write_record(record)
        writer.finish()

        # Read back
        stream = biometal.VcfStream.from_path(temp_path)
        read_header = stream.header()

        # Verify header
        assert read_header.samples == ["NA12878", "NA12879"]
        assert "GT" in read_header.format_fields

        # Verify record
        read_back = list(stream)
        assert len(read_back) == 1
        assert read_back[0].format == "GT:DP"
        assert len(read_back[0].samples) == 2

        print("  ✓ VCF with samples successful")

    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


if __name__ == '__main__':
    print("VCF Writer Tests")
    print("=" * 60)

    try:
        test_vcf_writer_basic()
        test_vcf_writer_compressed()
        test_vcf_writer_with_samples()

        print("\n" + "=" * 60)
        print("✓ All VCF writer tests passed!")
        print("=" * 60)

    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
