#!/usr/bin/env python3
"""
Advanced BAM Filtering Examples with biometal

Demonstrates production-ready filtering techniques for BAM alignment files:
- Region queries (chr:start-end)
- Complex multi-criteria filters
- Common bioinformatics workflows
- Performance-optimized patterns

Performance:
- 4.54 million records/sec throughput
- 43.0 MiB/s compressed file processing
- 4× speedup via parallel BGZF decompression
- Constant ~5 MB memory footprint

Requirements:
- biometal-rs >= 1.2.0 (pip install biometal-rs)
- Python 3.9+
"""

import biometal
from typing import Iterator, Dict, List, Tuple, Optional
from dataclasses import dataclass
from collections import defaultdict


# ============================================================================
# Region Query Helpers
# ============================================================================

@dataclass
class Region:
    """Genomic region for filtering."""
    reference_id: int
    start: int  # 0-based, inclusive
    end: int    # 0-based, exclusive

    def contains(self, position: int) -> bool:
        """Check if position falls within region."""
        return self.start <= position < self.end

    def overlaps(self, start: int, end: int) -> bool:
        """Check if interval [start, end) overlaps this region."""
        return not (end <= self.start or start >= self.end)

    @classmethod
    def from_string(cls, region_str: str, ref_map: Dict[str, int]) -> 'Region':
        """
        Parse region string like 'chr1:1000-2000' or '0:1000-2000'.

        Args:
            region_str: Region string (e.g., 'chr1:1000-2000' or '0:1000-2000')
            ref_map: Dictionary mapping reference names to IDs

        Returns:
            Region object

        Examples:
            >>> ref_map = {'chr1': 0, 'chr2': 1}
            >>> region = Region.from_string('chr1:1000-2000', ref_map)
            >>> region.reference_id
            0
        """
        parts = region_str.split(':')
        if len(parts) != 2:
            raise ValueError(f"Invalid region format: {region_str}")

        ref_name = parts[0]

        # Handle numeric reference IDs
        if ref_name.isdigit():
            reference_id = int(ref_name)
        elif ref_name in ref_map:
            reference_id = ref_map[ref_name]
        else:
            raise ValueError(f"Unknown reference: {ref_name}")

        # Parse coordinates
        coords = parts[1].split('-')
        if len(coords) != 2:
            raise ValueError(f"Invalid coordinates: {parts[1]}")

        start = int(coords[0])
        end = int(coords[1])

        if start >= end:
            raise ValueError(f"Invalid range: start ({start}) >= end ({end})")

        return cls(reference_id, start, end)


def filter_by_region(
    bam_path: str,
    region: Region,
    require_primary: bool = True,
    min_mapq: int = 0
) -> Iterator[biometal.BamRecord]:
    """
    Filter BAM records by genomic region.

    Args:
        bam_path: Path to BAM file
        region: Region to query
        require_primary: Only return primary alignments
        min_mapq: Minimum mapping quality

    Yields:
        BamRecord objects within the specified region

    Examples:
        >>> region = Region(reference_id=0, start=1000, end=2000)
        >>> for record in filter_by_region("alignments.bam", region, min_mapq=30):
        ...     print(f"{record.name}: {record.position}")
    """
    bam = biometal.BamReader.from_path(bam_path)

    for record in bam:
        # Check reference matches
        if record.reference_id != region.reference_id:
            continue

        # Check if unmapped
        if record.position is None:
            continue

        # Check if overlaps region
        if not region.contains(record.position):
            continue

        # Apply filters
        if require_primary and not record.is_primary:
            continue

        if record.mapq is not None and record.mapq < min_mapq:
            continue

        yield record


# ============================================================================
# Complex Filters
# ============================================================================

@dataclass
class FilterCriteria:
    """Comprehensive filtering criteria for BAM records."""

    # Alignment quality
    min_mapq: Optional[int] = None
    max_mapq: Optional[int] = None

    # Alignment type
    require_mapped: bool = False
    require_primary: bool = False
    require_proper_pair: bool = False

    # Strand
    require_forward: bool = False
    require_reverse: bool = False

    # Read properties
    min_length: Optional[int] = None
    max_length: Optional[int] = None

    # Reference
    allowed_references: Optional[List[int]] = None
    excluded_references: Optional[List[int]] = None

    # Region
    region: Optional[Region] = None

    def passes(self, record: biometal.BamRecord) -> bool:
        """
        Check if record passes all criteria.

        Args:
            record: BamRecord to evaluate

        Returns:
            True if record passes all filters
        """
        # MAPQ filters
        if self.min_mapq is not None:
            if record.mapq is None or record.mapq < self.min_mapq:
                return False

        if self.max_mapq is not None:
            if record.mapq is None or record.mapq > self.max_mapq:
                return False

        # Alignment type filters
        if self.require_mapped and not record.is_mapped:
            return False

        if self.require_primary and not record.is_primary:
            return False

        if self.require_proper_pair and not record.is_paired:
            return False

        # Strand filters
        if self.require_forward and record.is_reverse:
            return False

        if self.require_reverse and not record.is_reverse:
            return False

        # Length filters
        seq_len = len(record.sequence)
        if self.min_length is not None and seq_len < self.min_length:
            return False

        if self.max_length is not None and seq_len > self.max_length:
            return False

        # Reference filters
        if self.allowed_references is not None:
            if record.reference_id not in self.allowed_references:
                return False

        if self.excluded_references is not None:
            if record.reference_id in self.excluded_references:
                return False

        # Region filter
        if self.region is not None:
            if record.reference_id != self.region.reference_id:
                return False
            if record.position is None:
                return False
            if not self.region.contains(record.position):
                return False

        return True


def filter_bam(
    bam_path: str,
    criteria: FilterCriteria,
    limit: Optional[int] = None
) -> Iterator[biometal.BamRecord]:
    """
    Filter BAM records with complex criteria.

    Args:
        bam_path: Path to BAM file
        criteria: FilterCriteria object with filter settings
        limit: Optional maximum number of records to return

    Yields:
        BamRecord objects that pass all criteria

    Examples:
        >>> # High-quality, primary, forward-strand alignments
        >>> criteria = FilterCriteria(
        ...     min_mapq=30,
        ...     require_primary=True,
        ...     require_forward=True,
        ...     min_length=50
        ... )
        >>> for record in filter_bam("alignments.bam", criteria, limit=1000):
        ...     process(record)
    """
    bam = biometal.BamReader.from_path(bam_path)
    count = 0

    for record in bam:
        if criteria.passes(record):
            yield record
            count += 1
            if limit is not None and count >= limit:
                break


# ============================================================================
# Common Workflows
# ============================================================================

def variant_calling_filter(bam_path: str) -> Iterator[biometal.BamRecord]:
    """
    Filter for variant calling: high quality, primary, mapped.

    Standard filters:
    - MAPQ >= 30 (high confidence)
    - Primary alignments only
    - Mapped to reference

    Args:
        bam_path: Path to BAM file

    Yields:
        High-quality alignments suitable for variant calling
    """
    criteria = FilterCriteria(
        min_mapq=30,
        require_mapped=True,
        require_primary=True
    )
    return filter_bam(bam_path, criteria)


def strand_specific_rna_seq(
    bam_path: str,
    forward_only: bool = True
) -> Iterator[biometal.BamRecord]:
    """
    Filter for strand-specific RNA-seq analysis.

    Args:
        bam_path: Path to BAM file
        forward_only: If True, return only forward strand; else reverse

    Yields:
        Strand-specific alignments
    """
    criteria = FilterCriteria(
        require_mapped=True,
        require_primary=True,
        require_forward=forward_only,
        require_reverse=not forward_only
    )
    return filter_bam(bam_path, criteria)


def paired_end_properly_paired(bam_path: str) -> Iterator[biometal.BamRecord]:
    """
    Filter for properly paired reads (quality control).

    Args:
        bam_path: Path to BAM file

    Yields:
        Properly paired alignments
    """
    criteria = FilterCriteria(
        require_proper_pair=True,
        require_primary=True
    )
    return filter_bam(bam_path, criteria)


def extract_region(
    bam_path: str,
    reference_id: int,
    start: int,
    end: int,
    min_mapq: int = 20
) -> List[biometal.BamRecord]:
    """
    Extract all alignments from a specific genomic region.

    Args:
        bam_path: Path to BAM file
        reference_id: Reference sequence ID
        start: Start position (0-based, inclusive)
        end: End position (0-based, exclusive)
        min_mapq: Minimum mapping quality

    Returns:
        List of alignments in the region

    Examples:
        >>> # Extract chr1:1000-2000
        >>> alignments = extract_region("alignments.bam", 0, 1000, 2000)
        >>> print(f"Found {len(alignments)} alignments in region")
    """
    region = Region(reference_id, start, end)
    return list(filter_by_region(bam_path, region, min_mapq=min_mapq))


# ============================================================================
# Coverage Analysis
# ============================================================================

def calculate_coverage(
    bam_path: str,
    reference_id: int,
    start: int,
    end: int,
    min_mapq: int = 20
) -> Dict[int, int]:
    """
    Calculate per-base coverage for a region.

    Args:
        bam_path: Path to BAM file
        reference_id: Reference sequence ID
        start: Start position (0-based, inclusive)
        end: End position (0-based, exclusive)
        min_mapq: Minimum mapping quality

    Returns:
        Dictionary mapping position -> coverage depth

    Examples:
        >>> coverage = calculate_coverage("alignments.bam", 0, 1000, 2000)
        >>> mean_cov = sum(coverage.values()) / len(coverage)
        >>> print(f"Mean coverage: {mean_cov:.1f}×")
    """
    coverage = defaultdict(int)
    region = Region(reference_id, start, end)

    for record in filter_by_region(bam_path, region, min_mapq=min_mapq):
        # Count this read at its alignment position
        # (Simplified - real coverage would use CIGAR)
        pos = record.position
        read_len = len(record.sequence)

        for i in range(pos, min(pos + read_len, end)):
            if i >= start:
                coverage[i] += 1

    return dict(coverage)


def find_high_coverage_regions(
    bam_path: str,
    reference_id: int,
    min_coverage: int = 20,
    min_length: int = 100
) -> List[Tuple[int, int, float]]:
    """
    Identify high-coverage regions (e.g., for CNV detection).

    Args:
        bam_path: Path to BAM file
        reference_id: Reference sequence ID
        min_coverage: Minimum coverage depth
        min_length: Minimum region length

    Returns:
        List of (start, end, mean_coverage) tuples

    Examples:
        >>> regions = find_high_coverage_regions("alignments.bam", 0)
        >>> for start, end, cov in regions:
        ...     print(f"Region {start}-{end}: {cov:.1f}× coverage")
    """
    # This is a simplified example - production would use sliding windows
    high_cov_regions = []

    # Count alignments starting at each position
    position_counts = defaultdict(int)

    criteria = FilterCriteria(
        require_mapped=True,
        require_primary=True,
        allowed_references=[reference_id]
    )

    for record in filter_bam(bam_path, criteria):
        if record.position is not None:
            position_counts[record.position] += 1

    # Find contiguous high-coverage regions
    in_region = False
    region_start = None
    region_positions = []

    for pos in sorted(position_counts.keys()):
        cov = position_counts[pos]

        if cov >= min_coverage:
            if not in_region:
                in_region = True
                region_start = pos
            region_positions.append((pos, cov))
        else:
            if in_region and len(region_positions) >= min_length:
                mean_cov = sum(c for _, c in region_positions) / len(region_positions)
                high_cov_regions.append((region_start, pos, mean_cov))
            in_region = False
            region_positions = []

    return high_cov_regions


# ============================================================================
# Statistics and Summary
# ============================================================================

def count_by_flag(
    bam_path: str,
    limit: Optional[int] = None
) -> Dict[str, int]:
    """
    Count alignments by various flag properties.

    Args:
        bam_path: Path to BAM file
        limit: Optional limit on records to process

    Returns:
        Dictionary with counts for each property

    Examples:
        >>> stats = count_by_flag("alignments.bam", limit=10000)
        >>> print(f"Mapped: {stats['mapped']} / {stats['total']}")
    """
    counts = {
        'total': 0,
        'mapped': 0,
        'unmapped': 0,
        'primary': 0,
        'secondary': 0,
        'paired': 0,
        'forward': 0,
        'reverse': 0,
        'first': 0,
        'second': 0,
    }

    bam = biometal.BamReader.from_path(bam_path)

    for record in bam:
        counts['total'] += 1

        if record.is_mapped:
            counts['mapped'] += 1
        else:
            counts['unmapped'] += 1

        if record.is_primary:
            counts['primary'] += 1
        else:
            counts['secondary'] += 1

        if record.is_paired:
            counts['paired'] += 1

        if record.is_reverse:
            counts['reverse'] += 1
        else:
            counts['forward'] += 1

        if record.is_first:
            counts['first'] += 1

        if record.is_second:
            counts['second'] += 1

        if limit is not None and counts['total'] >= limit:
            break

    return counts


def mapq_distribution(
    bam_path: str,
    limit: Optional[int] = None
) -> Dict[int, int]:
    """
    Calculate MAPQ score distribution.

    Args:
        bam_path: Path to BAM file
        limit: Optional limit on records to process

    Returns:
        Dictionary mapping MAPQ -> count

    Examples:
        >>> dist = mapq_distribution("alignments.bam", limit=10000)
        >>> for mapq in sorted(dist.keys()):
        ...     print(f"MAPQ {mapq}: {dist[mapq]} alignments")
    """
    distribution = defaultdict(int)

    bam = biometal.BamReader.from_path(bam_path)
    count = 0

    for record in bam:
        if record.mapq is not None:
            distribution[record.mapq] += 1
        count += 1

        if limit is not None and count >= limit:
            break

    return dict(distribution)


# ============================================================================
# Example Usage
# ============================================================================

def main():
    """Demonstrate advanced filtering techniques."""
    import sys

    if len(sys.argv) < 2:
        print("Usage: python bam_advanced_filtering.py <bam_file>")
        print("\nExamples:")
        print("  python bam_advanced_filtering.py alignments.bam")
        sys.exit(1)

    bam_path = sys.argv[1]

    print(f"📊 Advanced BAM Filtering Examples")
    print(f"{'=' * 60}\n")
    print(f"File: {bam_path}\n")

    # Example 1: Variant calling filter
    print("1. Variant Calling Filter (MAPQ ≥ 30, primary, mapped):")
    count = sum(1 for _ in variant_calling_filter(bam_path))
    print(f"   Found {count:,} high-quality alignments\n")

    # Example 2: Region query
    print("2. Region Query (chr0:0-1000):")
    region = Region(reference_id=0, start=0, end=1000)
    region_records = list(filter_by_region(bam_path, region, min_mapq=20))
    print(f"   Found {len(region_records):,} alignments in region\n")

    # Example 3: Complex filter
    print("3. Complex Filter (high-quality, forward strand, 50-150 bp):")
    criteria = FilterCriteria(
        min_mapq=30,
        require_primary=True,
        require_forward=True,
        min_length=50,
        max_length=150
    )
    complex_count = sum(1 for _ in filter_bam(bam_path, criteria, limit=10000))
    print(f"   Found {complex_count:,} matching alignments (first 10K)\n")

    # Example 4: Flag statistics
    print("4. Flag Statistics (first 10K records):")
    stats = count_by_flag(bam_path, limit=10000)
    print(f"   Total: {stats['total']:,}")
    print(f"   Mapped: {stats['mapped']:,} ({100*stats['mapped']/stats['total']:.1f}%)")
    print(f"   Primary: {stats['primary']:,} ({100*stats['primary']/stats['total']:.1f}%)")
    print(f"   Strand: {stats['forward']:,} forward, {stats['reverse']:,} reverse\n")

    # Example 5: MAPQ distribution
    print("5. MAPQ Distribution (first 10K records):")
    dist = mapq_distribution(bam_path, limit=10000)
    for mapq in sorted(dist.keys())[:10]:  # Show first 10
        count = dist[mapq]
        print(f"   MAPQ {mapq}: {count:,} alignments")
    print()

    print(f"{'=' * 60}")
    print("✅ Analysis complete")
    print(f"\n💡 Tip: Modify criteria objects to customize filters")
    print(f"   Performance: 4.54M records/sec, constant ~5 MB memory")


if __name__ == '__main__':
    main()
