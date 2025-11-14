#!/usr/bin/env python3
"""
VCF (Variant Call Format) parsing with biometal.

Demonstrates:
- Parsing VCF v4.2 format
- Header and metadata extraction
- Variant classification (SNP, indel, etc.)
- Quality filtering
- INFO field parsing
"""

import biometal
from pathlib import Path
import tempfile
from collections import defaultdict


def demo_basic_parsing():
    """Parse basic VCF file"""
    print("=== VCF Format Parsing ===\n")

    # Create sample VCF data
    vcf_data = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t12345\trs123\tA\tT\t30.0\tPASS\tDP=100;AF=0.25
chr1\t23456\t.\tG\tC\t50.0\tPASS\tDP=150;AF=0.5
chr2\t34567\trs456\tATG\tA\t45.0\tPASS\tDP=80;AF=0.3"""

    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
        f.write(vcf_data)
        vcf_path = f.name

    try:
        # Parse VCF file
        stream = biometal.VcfStream.from_path(vcf_path)

        # Get header
        header = stream.header()
        print(f"VCF version: {header.fileformat}")
        print(f"Contigs: {len(header.contigs)}")
        print(f"INFO fields: {len(header.info_fields)}")
        print(f"FORMAT fields: {len(header.format_fields)}\n")

        # Parse variants
        variants = []
        for record in stream:
            variants.append(record)

        print(f"Total variants: {len(variants)}\n")

        # Show variant details
        print("Variant details:")
        for var in variants:
            var_id = var.id or "."
            quality = f"{var.quality:.1f}" if var.quality else "."
            print(f"  {var.chrom}:{var.pos} {var.reference}->{','.join(var.alternate)}")
            print(f"    ID: {var_id}, Quality: {quality}, Filter: {var.filter}")
            if var.info:
                info_str = "; ".join(f"{k}={v}" for k, v in var.info.items())
                print(f"    INFO: {info_str}")

    finally:
        Path(vcf_path).unlink()


def demo_variant_classification():
    """Classify variants by type"""
    print("\n=== Variant Classification ===\n")

    vcf_data = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1000\t.\tA\tT\t.\tPASS\t.
chr1\t2000\t.\tG\tC\t.\tPASS\t.
chr1\t3000\t.\tA\tATG\t.\tPASS\t.
chr1\t4000\t.\tATG\tA\t.\tPASS\t.
chr1\t5000\t.\tA\tT,G,C\t.\tPASS\t."""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
        f.write(vcf_data)
        vcf_path = f.name

    try:
        stream = biometal.VcfStream.from_path(vcf_path)
        _ = stream.header()  # Skip header

        # Classify variants
        snps = []
        insertions = []
        deletions = []
        multi_allelic = []

        for record in stream:
            if len(record.alternate) > 1:
                multi_allelic.append(record)
            elif record.is_snp():
                snps.append(record)
            elif record.is_insertion():
                insertions.append(record)
            elif record.is_deletion():
                deletions.append(record)

        print(f"SNPs: {len(snps)}")
        print(f"Insertions: {len(insertions)}")
        print(f"Deletions: {len(deletions)}")
        print(f"Multi-allelic: {len(multi_allelic)}\n")

        # Show examples
        if snps:
            print("Example SNP:")
            snp = snps[0]
            print(f"  {snp.chrom}:{snp.pos} {snp.reference}->{snp.alternate[0]}")

        if insertions:
            print("\nExample insertion:")
            ins = insertions[0]
            ins_size = len(ins.alternate[0]) - len(ins.reference)
            print(f"  {ins.chrom}:{ins.pos} +{ins_size}bp ({ins.alternate[0]})")

        if deletions:
            print("\nExample deletion:")
            del_var = deletions[0]
            del_size = len(del_var.reference) - len(del_var.alternate[0])
            print(f"  {del_var.chrom}:{del_var.pos} -{del_size}bp")

        if multi_allelic:
            print("\nExample multi-allelic:")
            multi = multi_allelic[0]
            alts = ','.join(multi.alternate)
            print(f"  {multi.chrom}:{multi.pos} {multi.reference}->{alts}")

    finally:
        Path(vcf_path).unlink()


def demo_quality_filtering():
    """Filter variants by quality"""
    print("\n=== Quality Filtering ===\n")

    vcf_data = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1000\t.\tA\tT\t10.0\tLowQual\tDP=10
chr1\t2000\t.\tG\tC\t30.0\tPASS\tDP=50
chr1\t3000\t.\tT\tA\t50.0\tPASS\tDP=100
chr1\t4000\t.\tC\tG\t70.0\tPASS\tDP=150
chr1\t5000\t.\tA\tT\t5.0\tLowQual\tDP=5"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
        f.write(vcf_data)
        vcf_path = f.name

    try:
        stream = biometal.VcfStream.from_path(vcf_path)
        _ = stream.header()

        # Filter by quality
        quality_threshold = 30.0
        high_quality = []

        for record in stream:
            if record.quality and record.quality >= quality_threshold:
                if record.filter == "PASS":
                    high_quality.append(record)

        print(f"High-quality variants (QUAL >= {quality_threshold}):")
        print(f"  Found: {len(high_quality)}\n")

        # Quality distribution
        stream2 = biometal.VcfStream.from_path(vcf_path)
        _ = stream2.header()

        quality_bins = defaultdict(int)
        for record in stream2:
            if record.quality:
                bin_val = int(record.quality // 10) * 10
                quality_bins[bin_val] += 1

        print("Quality score distribution:")
        for bin_val in sorted(quality_bins.keys()):
            count = quality_bins[bin_val]
            bar = "█" * count
            print(f"  {bin_val:3d}-{bin_val+9:3d}: {bar} ({count})")

    finally:
        Path(vcf_path).unlink()


def demo_info_parsing():
    """Parse and analyze INFO fields"""
    print("\n=== INFO Field Analysis ===\n")

    vcf_data = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele Count">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1000\t.\tA\tT\t30.0\tPASS\tDP=50;AF=0.1;AC=5
chr1\t2000\t.\tG\tC\t40.0\tPASS\tDP=100;AF=0.25;AC=25
chr1\t3000\t.\tT\tA\t50.0\tPASS\tDP=200;AF=0.5;AC=100
chr1\t4000\t.\tC\tG\t60.0\tPASS\tDP=80;AF=0.05;AC=4"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
        f.write(vcf_data)
        vcf_path = f.name

    try:
        stream = biometal.VcfStream.from_path(vcf_path)
        _ = stream.header()

        # Collect INFO statistics
        depths = []
        allele_freqs = []

        for record in stream:
            if 'DP' in record.info:
                depths.append(int(record.info['DP']))

            if 'AF' in record.info:
                allele_freqs.append(float(record.info['AF']))

        # Calculate statistics
        if depths:
            avg_depth = sum(depths) / len(depths)
            print(f"Coverage statistics:")
            print(f"  Mean depth: {avg_depth:.1f}x")
            print(f"  Min depth: {min(depths)}x")
            print(f"  Max depth: {max(depths)}x")

        if allele_freqs:
            avg_af = sum(allele_freqs) / len(allele_freqs)
            print(f"\nAllele frequency statistics:")
            print(f"  Mean AF: {avg_af:.3f}")
            print(f"  Min AF: {min(allele_freqs):.3f}")
            print(f"  Max AF: {max(allele_freqs):.3f}")

            # Classify by frequency
            rare = sum(1 for af in allele_freqs if af < 0.05)
            common = sum(1 for af in allele_freqs if af >= 0.05)
            print(f"\n  Rare variants (AF < 0.05): {rare}")
            print(f"  Common variants (AF >= 0.05): {common}")

    finally:
        Path(vcf_path).unlink()


if __name__ == "__main__":
    print("VCF Variant Calling Format Parsing with biometal\n")
    print("Streaming architecture for population genetics")
    print("=" * 60 + "\n")

    demo_basic_parsing()
    demo_variant_classification()
    demo_quality_filtering()
    demo_info_parsing()

    print("\n✓ All VCF format demos complete!")
