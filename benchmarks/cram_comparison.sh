#!/bin/bash
# CRAM reader performance comparison: biometal vs samtools
# Usage: ./benchmarks/cram_comparison.sh <cram_file> <reference>

set -e

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <cram_file> <reference>"
    exit 1
fi

CRAM_FILE="$1"
REFERENCE="$2"

echo "================================================================"
echo "CRAM Reader Performance Comparison"
echo "================================================================"
echo "CRAM file:  $CRAM_FILE"
echo "Reference:  $REFERENCE"
echo "Date:       $(date)"
echo ""

# Check if samtools is installed
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools not found. Please install samtools first."
    exit 1
fi

echo "samtools version:"
samtools --version | head -1
echo ""

# Benchmark 1: samtools view (count records)
echo "================================================================"
echo "Benchmark 1: samtools view (count records)"
echo "================================================================"
for i in {1..3}; do
    echo "Run $i:"
    time samtools view -c --reference "$REFERENCE" "$CRAM_FILE"
    echo ""
done

# Benchmark 2: biometal CRAM reader
echo "================================================================"
echo "Benchmark 2: biometal CRAM reader"
echo "================================================================"
for i in {1..3}; do
    echo "Run $i:"
    cargo run --release --example cram_benchmark -- "$CRAM_FILE" "$REFERENCE"
    echo ""
done

echo "================================================================"
echo "Comparison Complete"
echo "================================================================"
