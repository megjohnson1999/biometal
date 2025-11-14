#!/usr/bin/env python3
"""
GFA (Graphical Fragment Assembly) format parsing with biometal.

Demonstrates:
- Parsing assembly graphs
- Segment, link, and path records
- Graph connectivity analysis
- Path extraction
"""

import biometal
from pathlib import Path
import tempfile
from collections import defaultdict


def demo_basic_parsing():
    """Parse basic GFA assembly graph"""
    print("=== GFA Format Parsing ===\n")

    # Create sample GFA data (simple assembly graph)
    gfa_data = """H\tVN:Z:1.0
S\tcontig1\tACGTACGT\tLN:i:8\tKC:i:100
S\tcontig2\tTGCATGCA\tLN:i:8\tKC:i:150
S\tcontig3\tGGAAGGAA\tLN:i:8\tKC:i:200
L\tcontig1\t+\tcontig2\t+\t4M
L\tcontig2\t+\tcontig3\t+\t4M
P\tpath1\tcontig1+,contig2+,contig3+\t4M,4M"""

    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.gfa', delete=False) as f:
        f.write(gfa_data)
        gfa_path = f.name

    try:
        # Parse GFA file
        stream = biometal.GfaStream.from_path(gfa_path)

        segments = []
        links = []
        paths = []

        for record in stream:
            if isinstance(record, biometal.GfaSegment):
                segments.append(record)
            elif isinstance(record, biometal.GfaLink):
                links.append(record)
            elif isinstance(record, biometal.GfaPath):
                paths.append(record)

        print(f"Segments: {len(segments)}")
        print(f"Links: {len(links)}")
        print(f"Paths: {len(paths)}\n")

        # Show segments
        print("Segment details:")
        for seg in segments:
            length = seg.length() or 0
            coverage = seg.coverage()
            cov_str = f", coverage={coverage:.1f}x" if coverage else ""
            print(f"  {seg.name}: {length} bp{cov_str}")

        # Show links
        print("\nGraph connections:")
        for link in links:
            print(f"  {link.from_segment}{link.from_orient} -> "
                  f"{link.to_segment}{link.to_orient} (overlap: {link.overlap})")

        # Show paths
        print("\nPaths through graph:")
        for path in paths:
            print(f"  {path.name}: {' -> '.join(path.segments)}")

    finally:
        Path(gfa_path).unlink()


def demo_graph_analysis():
    """Analyze graph connectivity and structure"""
    print("\n=== Graph Analysis ===\n")

    # Create branching assembly graph
    gfa_data = """H\tVN:Z:1.0
S\tstart\tACGT\tLN:i:4
S\tbranch1\tTGCA\tLN:i:4
S\tbranch2\tGGAA\tLN:i:4
S\tend\tTTCC\tLN:i:4
L\tstart\t+\tbranch1\t+\t2M
L\tstart\t+\tbranch2\t+\t2M
L\tbranch1\t+\tend\t+\t2M
L\tbranch2\t+\tend\t+\t2M"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.gfa', delete=False) as f:
        f.write(gfa_data)
        gfa_path = f.name

    try:
        stream = biometal.GfaStream.from_path(gfa_path)

        # Build adjacency lists
        outgoing = defaultdict(list)
        incoming = defaultdict(list)

        segments_by_name = {}

        for record in stream:
            if isinstance(record, biometal.GfaSegment):
                segments_by_name[record.name] = record
            elif isinstance(record, biometal.GfaLink):
                outgoing[record.from_segment].append(
                    (record.to_segment, record.from_orient, record.to_orient)
                )
                incoming[record.to_segment].append(
                    (record.from_segment, record.from_orient, record.to_orient)
                )

        # Analyze graph structure
        print("Node degree analysis:")
        for name in sorted(segments_by_name.keys()):
            in_degree = len(incoming[name])
            out_degree = len(outgoing[name])
            total_degree = in_degree + out_degree

            node_type = "intermediate"
            if in_degree == 0:
                node_type = "source"
            elif out_degree == 0:
                node_type = "sink"
            elif in_degree > 1 or out_degree > 1:
                node_type = "branch point"

            print(f"  {name}: in={in_degree}, out={out_degree} ({node_type})")

        # Find source and sink nodes
        sources = [n for n in segments_by_name if len(incoming[n]) == 0]
        sinks = [n for n in segments_by_name if len(outgoing[n]) == 0]

        print(f"\nSource nodes (in-degree 0): {', '.join(sources)}")
        print(f"Sink nodes (out-degree 0): {', '.join(sinks)}")

        # Find branch points
        branches = [n for n in segments_by_name
                   if len(incoming[n]) > 1 or len(outgoing[n]) > 1]
        if branches:
            print(f"Branch points: {', '.join(branches)}")

    finally:
        Path(gfa_path).unlink()


def demo_path_extraction():
    """Extract paths through assembly graph"""
    print("\n=== Path Extraction ===\n")

    # Create graph with multiple paths
    gfa_data = """H\tVN:Z:1.0
S\ts1\tACGT\tLN:i:4
S\ts2\tTGCA\tLN:i:4
S\ts3\tGGAA\tLN:i:4
S\ts4\tTTCC\tLN:i:4
L\ts1\t+\ts2\t+\t2M
L\ts2\t+\ts3\t+\t2M
L\ts3\t+\ts4\t+\t2M
P\tpath_main\ts1+,s2+,s3+,s4+\t2M,2M,2M
P\tpath_sub\ts2+,s3+\t2M"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.gfa', delete=False) as f:
        f.write(gfa_data)
        gfa_path = f.name

    try:
        stream = biometal.GfaStream.from_path(gfa_path)

        segments = {}
        paths = []

        for record in stream:
            if isinstance(record, biometal.GfaSegment):
                segments[record.name] = record
            elif isinstance(record, biometal.GfaPath):
                paths.append(record)

        print("Paths through assembly:")
        for path in paths:
            # Calculate total path length
            total_length = sum(
                segments[seg.rstrip('+-')].length() or 0
                for seg in path.segments
                if seg.rstrip('+-') in segments
            )

            print(f"\n{path.name}:")
            print(f"  Segments: {len(path.segments)}")
            print(f"  Total length: {total_length} bp")
            print(f"  Path: {' -> '.join(path.segments)}")

            # Show segment details
            for i, seg_name in enumerate(path.segments):
                clean_name = seg_name.rstrip('+-')
                if clean_name in segments:
                    seg = segments[clean_name]
                    length = seg.length() or 0
                    print(f"    {i+1}. {clean_name}: {length} bp")

    finally:
        Path(gfa_path).unlink()


if __name__ == "__main__":
    print("GFA Assembly Graph Parsing with biometal\n")
    print("Streaming architecture for pangenomes and assemblies")
    print("=" * 60 + "\n")

    demo_basic_parsing()
    demo_graph_analysis()
    demo_path_extraction()

    print("\n✓ All GFA format demos complete!")
