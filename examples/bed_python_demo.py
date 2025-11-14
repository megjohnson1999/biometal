#!/usr/bin/env python3
"""
BED format parsing examples using biometal.

Demonstrates:
- Parsing BED3, BED6, and BED12 formats
- Streaming architecture (constant memory)
- Coordinate manipulation
- Feature filtering
"""

import biometal
from pathlib import Path
import tempfile


def demo_bed3():
    """BED3 format - minimal genomic intervals"""
    print("=== BED3 Format Demo ===\n")

    # Create sample BED3 data
    bed3_data = """chr1\t1000\t2000
chr1\t5000\t6000
chr2\t10000\t15000
chr2\t20000\t25000
chrX\t100000\t150000"""

    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write(bed3_data)
        bed_path = f.name

    try:
        # Parse BED3 file
        stream = biometal.Bed3Stream.from_path(bed_path)

        intervals_by_chr = {}
        for record in stream:
            chrom = record.chrom
            if chrom not in intervals_by_chr:
                intervals_by_chr[chrom] = []
            intervals_by_chr[chrom].append((record.start, record.end))

        # Report statistics
        print("Chromosome distribution:")
        for chrom, intervals in sorted(intervals_by_chr.items()):
            total_length = sum(end - start for start, end in intervals)
            print(f"  {chrom}: {len(intervals)} intervals, {total_length:,} bp")

    finally:
        Path(bed_path).unlink()


def demo_bed6():
    """BED6 format - standard annotations"""
    print("\n=== BED6 Format Demo ===\n")

    # Create sample BED6 data (ENCODE-style peaks)
    bed6_data = """chr1\t1000\t2000\tpeak1\t100\t+
chr1\t5000\t6000\tpeak2\t500\t+
chr2\t10000\t15000\tpeak3\t800\t-
chr2\t20000\t25000\tpeak4\t200\t-
chrX\t100000\t150000\tpeak5\t1000\t+"""

    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write(bed6_data)
        bed_path = f.name

    try:
        # Parse BED6 file
        stream = biometal.Bed6Stream.from_path(bed_path)

        # Filter high-quality peaks (score > 400)
        high_quality_peaks = []
        for record in stream:
            if record.score and record.score > 400:
                high_quality_peaks.append({
                    'name': record.name,
                    'location': f"{record.chrom}:{record.start}-{record.end}",
                    'score': record.score,
                    'strand': record.strand,
                    'length': record.end - record.start
                })

        print(f"Found {len(high_quality_peaks)} high-quality peaks (score > 400):")
        for peak in high_quality_peaks:
            print(f"  {peak['name']}: {peak['location']} ({peak['strand']}) - "
                  f"{peak['length']:,} bp, score={peak['score']}")

    finally:
        Path(bed_path).unlink()


def demo_bed12():
    """BED12 format - full gene structure"""
    print("\n=== BED12 Format Demo ===\n")

    # Create sample BED12 data (gene with 3 exons)
    bed12_data = """chr1\t1000\t5000\tgene1\t500\t+\t1200\t4800\t255,0,0\t3\t500,500,500\t0,2000,4000
chr2\t10000\t20000\tgene2\t800\t-\t11000\t19000\t0,255,0\t2\t2000,3000\t0,7000"""

    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write(bed12_data)
        bed_path = f.name

    try:
        # Parse BED12 file
        stream = biometal.Bed12Stream.from_path(bed_path)

        print("Gene structure analysis:")
        for record in stream:
            print(f"\nGene: {record.name}")
            print(f"  Location: {record.chrom}:{record.start}-{record.end} ({record.strand})")
            print(f"  Total length: {record.length():,} bp")

            if record.block_count:
                print(f"  Exons: {record.block_count}")

                if record.block_sizes and record.block_starts:
                    for i, (size, start) in enumerate(zip(record.block_sizes, record.block_starts)):
                        exon_start = record.start + start
                        exon_end = exon_start + size
                        print(f"    Exon {i+1}: {exon_start}-{exon_end} ({size} bp)")

            if record.thick_start and record.thick_end:
                cds_length = record.thick_end - record.thick_start
                print(f"  CDS: {record.thick_start}-{record.thick_end} ({cds_length:,} bp)")
                print(f"  Coding potential: {cds_length // 3} amino acids")

    finally:
        Path(bed_path).unlink()


def demo_coordinate_operations():
    """Demonstrate coordinate operations"""
    print("\n=== Coordinate Operations ===\n")

    # Parse a BED record
    record = biometal.Bed6Record.from_line("chr1\t1000\t2000\tfeature1\t100\t+")

    print(f"Original: {record.chrom}:{record.start}-{record.end}")
    print(f"Length: {record.length()} bp")
    print(f"Midpoint: {(record.start + record.end) // 2}")

    # Check overlap with another interval
    other_start, other_end = 1500, 2500
    overlap_start = max(record.start, other_start)
    overlap_end = min(record.end, other_end)

    if overlap_start < overlap_end:
        overlap_length = overlap_end - overlap_start
        print(f"\nOverlap with [{other_start}, {other_end}):")
        print(f"  Region: [{overlap_start}, {overlap_end})")
        print(f"  Length: {overlap_length} bp")


if __name__ == "__main__":
    print("BED Format Parsing with biometal\n")
    print("Constant memory streaming architecture")
    print("=" * 50 + "\n")

    demo_bed3()
    demo_bed6()
    demo_bed12()
    demo_coordinate_operations()

    print("\n✓ All BED format demos complete!")
