#!/usr/bin/env python3
"""
GFF3 (General Feature Format) parsing with biometal.

Demonstrates:
- Parsing gene annotations
- Hierarchical feature relationships (gene -> mRNA -> exon)
- Attribute extraction
- Gene structure analysis
"""

import biometal
from pathlib import Path
import tempfile
from collections import defaultdict


def demo_basic_parsing():
    """Parse basic GFF3 annotations"""
    print("=== GFF3 Format Parsing ===\n")

    # Create sample GFF3 data
    gff_data = """##gff-version 3
##sequence-region chr1 1 248956422
chr1\tEnsembl\tgene\t1000\t5000\t.\t+\t.\tID=gene1;Name=ABC1;biotype=protein_coding
chr1\tEnsembl\tmRNA\t1000\t5000\t.\t+\t.\tID=mRNA1;Parent=gene1;Name=ABC1-201
chr1\tEnsembl\texon\t1000\t1500\t.\t+\t.\tID=exon1;Parent=mRNA1
chr1\tEnsembl\texon\t2000\t2500\t.\t+\t.\tID=exon2;Parent=mRNA1
chr1\tEnsembl\texon\t4500\t5000\t.\t+\t.\tID=exon3;Parent=mRNA1
chr2\tEnsembl\tgene\t10000\t15000\t.\t-\t.\tID=gene2;Name=XYZ1;biotype=protein_coding"""

    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.gff3', delete=False) as f:
        f.write(gff_data)
        gff_path = f.name

    try:
        # Parse GFF3 file
        stream = biometal.Gff3Stream.from_path(gff_path)

        # Collect features by type
        features_by_type = defaultdict(list)

        for record in stream:
            features_by_type[record.feature_type].append(record)

        print("Feature summary:")
        for feature_type, features in sorted(features_by_type.items()):
            print(f"  {feature_type}: {len(features)}")

        print("\nGene details:")
        for gene in features_by_type['gene']:
            gene_id = gene.get_id()
            gene_name = gene.get_name()
            print(f"  {gene_name} ({gene_id})")
            print(f"    Location: {gene.seqid}:{gene.start}-{gene.end} ({gene.strand})")
            print(f"    Length: {gene.length():,} bp")

    finally:
        Path(gff_path).unlink()


def demo_hierarchical_features():
    """Demonstrate parent-child relationships"""
    print("\n=== Hierarchical Feature Analysis ===\n")

    gff_data = """##gff-version 3
chr1\tEnsembl\tgene\t1000\t5000\t.\t+\t.\tID=gene1;Name=ABC1
chr1\tEnsembl\tmRNA\t1000\t5000\t.\t+\t.\tID=mRNA1;Parent=gene1;Name=ABC1-201
chr1\tEnsembl\texon\t1000\t1500\t.\t+\t.\tID=exon1;Parent=mRNA1
chr1\tEnsembl\texon\t2000\t2500\t.\t+\t.\tID=exon2;Parent=mRNA1
chr1\tEnsembl\texon\t4500\t5000\t.\t+\t.\tID=exon3;Parent=mRNA1
chr1\tEnsembl\tCDS\t1200\t1500\t.\t+\t0\tID=cds1;Parent=mRNA1
chr1\tEnsembl\tCDS\t2000\t2500\t.\t+\t2\tID=cds2;Parent=mRNA1
chr1\tEnsembl\tCDS\t4500\t4800\t.\t+\t2\tID=cds3;Parent=mRNA1"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.gff3', delete=False) as f:
        f.write(gff_data)
        gff_path = f.name

    try:
        # Parse all features
        stream = biometal.Gff3Stream.from_path(gff_path)
        all_features = list(stream)

        # Build feature index
        features_by_id = {}
        for feature in all_features:
            feature_id = feature.get_id()
            if feature_id:
                features_by_id[feature_id] = feature

        # Analyze gene structure
        genes = [f for f in all_features if f.feature_type == 'gene']

        for gene in genes:
            gene_id = gene.get_id()
            gene_name = gene.get_name()

            print(f"Gene: {gene_name} ({gene_id})")
            print(f"  Location: {gene.seqid}:{gene.start}-{gene.end} ({gene.strand})")
            print(f"  Length: {gene.length():,} bp")

            # Find mRNAs for this gene
            mrnas = [f for f in all_features
                    if f.feature_type == 'mRNA' and f.get_parent() == gene_id]

            print(f"  Transcripts: {len(mrnas)}")

            for mrna in mrnas:
                mrna_id = mrna.get_id()
                mrna_name = mrna.get_name()

                print(f"\n  Transcript: {mrna_name} ({mrna_id})")

                # Find exons for this mRNA
                exons = [f for f in all_features
                        if f.feature_type == 'exon' and f.get_parent() == mrna_id]

                print(f"    Exons: {len(exons)}")
                for i, exon in enumerate(sorted(exons, key=lambda e: e.start)):
                    print(f"      Exon {i+1}: {exon.start}-{exon.end} ({exon.length()} bp)")

                # Find CDS for this mRNA
                cds_regions = [f for f in all_features
                              if f.feature_type == 'CDS' and f.get_parent() == mrna_id]

                if cds_regions:
                    total_cds = sum(cds.length() for cds in cds_regions)
                    print(f"    CDS regions: {len(cds_regions)}")
                    print(f"    Total CDS length: {total_cds} bp")
                    print(f"    Coding potential: {total_cds // 3} amino acids")

                    # Calculate UTRs
                    first_cds_start = min(cds.start for cds in cds_regions)
                    last_cds_end = max(cds.end for cds in cds_regions)

                    if mrna.strand == '+':
                        utr5_length = first_cds_start - mrna.start
                        utr3_length = mrna.end - last_cds_end
                        print(f"    5' UTR: {utr5_length} bp")
                        print(f"    3' UTR: {utr3_length} bp")
                    else:
                        utr3_length = first_cds_start - mrna.start
                        utr5_length = mrna.end - last_cds_end
                        print(f"    5' UTR: {utr5_length} bp")
                        print(f"    3' UTR: {utr3_length} bp")

    finally:
        Path(gff_path).unlink()


def demo_coordinate_conversion():
    """Demonstrate GFF3 to BED coordinate conversion"""
    print("\n=== Coordinate Conversion (GFF3 -> BED) ===\n")

    gff_data = """##gff-version 3
chr1\tEnsembl\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=ABC1"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.gff3', delete=False) as f:
        f.write(gff_data)
        gff_path = f.name

    try:
        stream = biometal.Gff3Stream.from_path(gff_path)

        for record in stream:
            # GFF3 coordinates (1-based, inclusive)
            print(f"GFF3 format (1-based, inclusive):")
            print(f"  {record.seqid}:{record.start}-{record.end}")
            print(f"  Length: {record.length()} bp")

            # Convert to BED coordinates (0-based, half-open)
            bed_start, bed_end = record.to_0based()
            print(f"\nBED format (0-based, half-open):")
            print(f"  {record.seqid}:{bed_start}-{bed_end}")
            print(f"  Length: {bed_end - bed_start} bp")

            print(f"\nConversion:")
            print(f"  GFF3 start {record.start} -> BED start {bed_start} (subtract 1)")
            print(f"  GFF3 end {record.end} -> BED end {bed_end} (no change)")

    finally:
        Path(gff_path).unlink()


def demo_chromosome_distribution():
    """Analyze gene distribution across chromosomes"""
    print("\n=== Chromosome Distribution ===\n")

    gff_data = """##gff-version 3
chr1\tEnsembl\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=ABC1
chr1\tEnsembl\tgene\t5000\t6000\t.\t-\t.\tID=gene2;Name=ABC2
chr2\tEnsembl\tgene\t10000\t11000\t.\t+\t.\tID=gene3;Name=XYZ1
chr2\tEnsembl\tgene\t20000\t21000\t.\t+\t.\tID=gene4;Name=XYZ2
chr2\tEnsembl\tgene\t30000\t31000\t.\t-\t.\tID=gene5;Name=XYZ3
chrX\tEnsembl\tgene\t100000\t101000\t.\t+\t.\tID=gene6;Name=XIST"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.gff3', delete=False) as f:
        f.write(gff_data)
        gff_path = f.name

    try:
        stream = biometal.Gff3Stream.from_path(gff_path)

        # Count genes per chromosome
        chr_counts = defaultdict(int)
        strand_counts = defaultdict(lambda: {'+': 0, '-': 0})

        for record in stream:
            if record.feature_type == 'gene':
                chr_counts[record.seqid] += 1
                strand_counts[record.seqid][record.strand] += 1

        print("Genes per chromosome:")
        for chrom in sorted(chr_counts.keys()):
            count = chr_counts[chrom]
            forward = strand_counts[chrom]['+']
            reverse = strand_counts[chrom]['-']
            print(f"  {chrom}: {count} genes (+ {forward}, - {reverse})")

    finally:
        Path(gff_path).unlink()


if __name__ == "__main__":
    print("GFF3 Gene Annotation Parsing with biometal\n")
    print("Streaming architecture for genomic features")
    print("=" * 60 + "\n")

    demo_basic_parsing()
    demo_hierarchical_features()
    demo_coordinate_conversion()
    demo_chromosome_distribution()

    print("\n✓ All GFF3 format demos complete!")
