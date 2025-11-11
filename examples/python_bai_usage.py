#!/usr/bin/env python3
"""
Example: Using BAI index for efficient BAM region queries in Python

This example demonstrates the performance difference between:
1. Legacy query() - O(n) full file scan
2. Indexed query_region() - O(log n) with BAI index
"""

import biometal
import time

def example_indexed_query():
    """Example of efficient indexed region query."""
    print("=" * 60)
    print("Example: Indexed BAM Region Query (O(log n))")
    print("=" * 60)

    # Load BAI index once (reusable for multiple queries)
    print("\n1. Loading BAI index...")
    index = biometal.BaiIndex.from_path("alignments.bam.bai")
    print(f"   Loaded index with {index.reference_count} references")

    # Query a specific region using index
    print("\n2. Querying chr1:1000-2000 with index...")
    start_time = time.time()

    record_count = 0
    for record in biometal.BamReader.query_region(
        "alignments.bam",
        index,
        "chr1",        # Reference name
        1000,          # Start position (0-based, inclusive)
        2000           # End position (0-based, exclusive)
    ):
        record_count += 1
        if record_count == 1:
            print(f"   First record: {record.name} at position {record.position}")

    elapsed = time.time() - start_time
    print(f"   Found {record_count} records in {elapsed:.3f}s")
    print(f"   Performance: O(log n) - only reads overlapping chunks")

def example_multiple_regions():
    """Example of querying multiple regions efficiently."""
    print("\n" + "=" * 60)
    print("Example: Multiple Region Queries (reusing index)")
    print("=" * 60)

    # Load index once, query multiple regions
    index = biometal.BaiIndex.from_path("alignments.bam.bai")

    regions = [
        ("chr1", 1000, 2000),
        ("chr1", 5000, 6000),
        ("chr2", 10000, 11000),
    ]

    for ref_name, start, end in regions:
        count = sum(1 for _ in biometal.BamReader.query_region(
            "alignments.bam",
            index,
            ref_name,
            start,
            end
        ))
        print(f"   {ref_name}:{start}-{end}: {count} records")

def example_legacy_comparison():
    """Show legacy query for comparison."""
    print("\n" + "=" * 60)
    print("Example: Legacy Query (O(n) - for comparison)")
    print("=" * 60)

    print("\n⚠ Warning: This scans the entire BAM file!")

    start_time = time.time()
    record_count = 0

    # Legacy query - scans entire file to find matching records
    for record in biometal.BamReader.query(
        "alignments.bam",
        reference_id=0,  # Numeric ID, not name
        start=1000,
        end=2000
    ):
        record_count += 1

    elapsed = time.time() - start_time
    print(f"   Found {record_count} records in {elapsed:.3f}s")
    print(f"   Performance: O(n) - reads entire file")

def main():
    """Run all examples."""
    try:
        example_indexed_query()
        example_multiple_regions()
        example_legacy_comparison()

        print("\n" + "=" * 60)
        print("Summary")
        print("=" * 60)
        print("✓ Indexed queries (query_region) provide O(log n) performance")
        print("✓ Load index once, reuse for multiple queries")
        print("✓ Legacy queries (query) scan entire file - use only when needed")

    except FileNotFoundError as e:
        print(f"\nℹ Example files not found: {e}")
        print("\nTo run this example:")
        print("  1. Create a BAM file: alignments.bam")
        print("  2. Index it: samtools index alignments.bam")
        print("  3. Run this script")

if __name__ == "__main__":
    main()
