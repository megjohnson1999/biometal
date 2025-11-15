#!/usr/bin/env python3
"""
Quick test for tab-delimited format writers (GFF3, GTF, PAF, narrowPeak).
"""

import biometal
import tempfile
import os


def test_gff3_writer():
    """Test Gff3Writer round-trip"""
    print("Testing Gff3Writer...")

    line = "chr1\tEnsembl\tgene\t1000\t5000\t.\t+\t.\tID=gene1;Name=ABC1"

    with tempfile.NamedTemporaryFile(suffix='.gff3.gz', delete=False) as f:
        temp_path = f.name

    try:
        # Write
        writer = biometal.Gff3Writer.create(temp_path)
        record = biometal.Gff3Record.from_line(line)
        writer.write_record(record)
        assert writer.records_written() == 1
        writer.finish()

        # Read back
        stream = biometal.Gff3Stream.from_path(temp_path)
        read_back = list(stream)

        # Verify
        assert len(read_back) == 1
        assert read_back[0].seqid == "chr1"
        assert read_back[0].feature_type == "gene"
        assert read_back[0].start == 1000

        print("  ✓ GFF3 round-trip successful")

    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


def test_gtf_writer():
    """Test GtfWriter round-trip"""
    print("\nTesting GtfWriter...")

    line = 'chr1\tENSEMBL\texon\t1000\t2000\t.\t+\t.\tgene_id "ENSG00001"; transcript_id "ENST00001";'

    with tempfile.NamedTemporaryFile(suffix='.gtf.gz', delete=False) as f:
        temp_path = f.name

    try:
        # Write
        writer = biometal.GtfWriter.create(temp_path)
        record = biometal.GtfRecord.from_line(line)
        writer.write_record(record)
        assert writer.records_written() == 1
        writer.finish()

        # Read back
        stream = biometal.GtfStream.from_path(temp_path)
        read_back = list(stream)

        # Verify
        assert len(read_back) == 1
        assert read_back[0].seqname == "chr1"
        assert read_back[0].feature == "exon"
        assert read_back[0].start == 1000

        print("  ✓ GTF round-trip successful")

    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


def test_paf_writer():
    """Test PafWriter round-trip"""
    print("\nTesting PafWriter...")

    # Minimal PAF line
    line = "query1\t100\t10\t90\t+\ttarget1\t200\t50\t150\t75\t80\t60"

    with tempfile.NamedTemporaryFile(suffix='.paf.gz', delete=False) as f:
        temp_path = f.name

    try:
        # Write
        writer = biometal.PafWriter.create(temp_path)
        record = biometal.PafRecord.from_line(line)
        writer.write_record(record)
        assert writer.records_written() == 1
        writer.finish()

        # Read back
        stream = biometal.PafStream.from_path(temp_path)
        read_back = list(stream)

        # Verify
        assert len(read_back) == 1
        assert read_back[0].query_name == "query1"
        assert read_back[0].target_name == "target1"

        print("  ✓ PAF round-trip successful")

    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


if __name__ == '__main__':
    print("Tab-Delimited Writer Tests (GFF3, GTF, PAF)")
    print("=" * 60)

    try:
        test_gff3_writer()
        test_gtf_writer()
        test_paf_writer()

        print("\n" + "=" * 60)
        print("✓ All tab-delimited writer tests passed!")
        print("=" * 60)

    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
