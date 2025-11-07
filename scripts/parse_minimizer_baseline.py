#!/usr/bin/env python3
"""
Parse Criterion benchmark results for Entry 035 baseline analysis.

Extracts mean, 95% CI, CV, and throughput for all minimizer configurations.
"""

import json
import os
from pathlib import Path
from typing import Dict, List, Tuple

def parse_estimates(estimates_path: Path) -> Dict:
    """Parse estimates.json from criterion output."""
    with open(estimates_path) as f:
        data = json.load(f)

    mean = data['mean']['point_estimate']
    mean_lower = data['mean']['confidence_interval']['lower_bound']
    mean_upper = data['mean']['confidence_interval']['upper_bound']
    std_dev = data['std_dev']['point_estimate']

    # Convert from nanoseconds to seconds
    mean_s = mean / 1e9
    mean_lower_s = mean_lower / 1e9
    mean_upper_s = mean_upper / 1e9

    # Calculate CV (coefficient of variation)
    cv = (std_dev / mean) * 100

    return {
        'mean_ns': mean,
        'mean_s': mean_s,
        'mean_lower_s': mean_lower_s,
        'mean_upper_s': mean_upper_s,
        'std_dev_ns': std_dev,
        'cv_percent': cv
    }

def calculate_throughput(mean_s: float, seq_length: int) -> float:
    """Calculate throughput in Mbp/s."""
    return (seq_length / mean_s) / 1e6

def parse_all_results(criterion_dir: Path) -> List[Dict]:
    """Parse all minimizer baseline results."""
    results = []

    # Expected configurations
    k_values = [21, 31]
    w_values = [11, 19]
    lengths = [100, 1000, 10000, 100000]

    for k in k_values:
        for w in w_values:
            for length in lengths:
                # Path: target/criterion/minimizer_baseline/k21_w11/100/base/estimates.json
                config_path = criterion_dir / f"k{k}_w{w}" / str(length) / "base" / "estimates.json"

                if not config_path.exists():
                    print(f"WARNING: Missing {config_path}")
                    continue

                stats = parse_estimates(config_path)
                throughput = calculate_throughput(stats['mean_s'], length)
                throughput_lower = calculate_throughput(stats['mean_upper_s'], length)  # inverted
                throughput_upper = calculate_throughput(stats['mean_lower_s'], length)  # inverted

                results.append({
                    'k': k,
                    'w': w,
                    'length': length,
                    'mean_s': stats['mean_s'],
                    'mean_lower_s': stats['mean_lower_s'],
                    'mean_upper_s': stats['mean_upper_s'],
                    'cv_percent': stats['cv_percent'],
                    'throughput_mbps': throughput,
                    'throughput_ci_lower': throughput_lower,
                    'throughput_ci_upper': throughput_upper
                })

    return results

def format_results_table(results: List[Dict]) -> str:
    """Format results as markdown table."""
    lines = []
    lines.append("| Configuration | Seq Length | Mean Time (ms) | 95% CI | CV (%) | Throughput (Mbp/s) | 95% CI |")
    lines.append("|---------------|------------|----------------|--------|--------|--------------------|----|")

    for r in sorted(results, key=lambda x: (x['k'], x['w'], x['length'])):
        config = f"k={r['k']}, w={r['w']}"
        length_str = format_length(r['length'])
        mean_ms = r['mean_s'] * 1000
        ci_str = f"[{r['mean_lower_s']*1000:.2f}, {r['mean_upper_s']*1000:.2f}]"
        cv_str = f"{r['cv_percent']:.1f}"
        throughput_str = f"{r['throughput_mbps']:.1f}"
        throughput_ci_str = f"[{r['throughput_ci_lower']:.1f}, {r['throughput_ci_upper']:.1f}]"

        lines.append(f"| {config} | {length_str} | {mean_ms:.2f} | {ci_str} | {cv_str} | {throughput_str} | {throughput_ci_str} |")

    return '\n'.join(lines)

def format_length(length: int) -> str:
    """Format sequence length (100, 1K, 10K, 100K)."""
    if length >= 1000:
        return f"{length//1000}K"
    return str(length)

def main():
    criterion_dir = Path("target/criterion/minimizer_baseline")

    if not criterion_dir.exists():
        print(f"ERROR: {criterion_dir} does not exist")
        print("Run: cargo bench --bench minimizer_baseline")
        return

    print("Parsing Entry 035 baseline results...\n")
    results = parse_all_results(criterion_dir)

    print(f"Found {len(results)} configurations")
    print(f"Expected: 16 configurations (4 k/w combinations × 4 lengths)\n")

    if len(results) < 16:
        print(f"WARNING: Only {len(results)}/16 configurations complete")
        print("Benchmark may still be running or was interrupted\n")

    print("=" * 80)
    print("ENTRY 035: Minimizer Extraction Baseline Results")
    print("=" * 80)
    print()
    print(format_results_table(results))
    print()

    # Summary statistics
    if results:
        throughputs = [r['throughput_mbps'] for r in results]
        cvs = [r['cv_percent'] for r in results]

        print("=" * 80)
        print("SUMMARY STATISTICS")
        print("=" * 80)
        print(f"Throughput range: {min(throughputs):.1f} - {max(throughputs):.1f} Mbp/s")
        print(f"Mean throughput: {sum(throughputs)/len(throughputs):.1f} Mbp/s")
        print(f"CV range: {min(cvs):.1f}% - {max(cvs):.1f}%")
        print(f"Mean CV: {sum(cvs)/len(cvs):.1f}%")
        print()

        # Scaling analysis
        print("=" * 80)
        print("SCALING ANALYSIS (k=21, w=11)")
        print("=" * 80)
        k21_w11 = [r for r in results if r['k'] == 21 and r['w'] == 11]
        for r in sorted(k21_w11, key=lambda x: x['length']):
            print(f"{format_length(r['length']):>6}: {r['throughput_mbps']:>6.1f} Mbp/s (CV: {r['cv_percent']:>4.1f}%)")
        print()

if __name__ == "__main__":
    main()
