#!/usr/bin/env bash
#
# Comprehensive benchmark comparing samtools and biometal
# Tests: sequential reading, indexed queries, memory usage
#
# Usage: ./samtools_vs_biometal.sh <bam_file> <bai_file>

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

BAM_FILE="${1:-tests/data/synthetic_100k.bam}"
BAI_FILE="${2:-tests/data/synthetic_100k.bam.bai}"

# Verify files exist
if [[ ! -f "$BAM_FILE" ]]; then
    echo -e "${RED}Error: BAM file not found: $BAM_FILE${NC}"
    exit 1
fi

if [[ ! -f "$BAI_FILE" ]]; then
    echo -e "${RED}Error: BAI file not found: $BAI_FILE${NC}"
    echo "Generate with: samtools index $BAM_FILE"
    exit 1
fi

echo "=========================================="
echo "samtools vs biometal Benchmark Suite"
echo "=========================================="
echo "BAM file: $BAM_FILE"
echo "BAI file: $BAI_FILE"
echo "Date: $(date)"
echo "Machine: $(uname -m)"
echo "samtools: $(samtools --version | head -1)"
echo "=========================================="
echo ""

# Create output directory
OUTPUT_DIR="benchmarks/comparison/results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"

# Helper function to run benchmark
run_benchmark() {
    local name="$1"
    local command="$2"
    local iterations="${3:-5}"

    echo -e "${BLUE}Benchmark: $name${NC}"
    echo "Command: $command"
    echo "Iterations: $iterations"

    local times=()
    local total=0

    for i in $(seq 1 $iterations); do
        echo -n "  Run $i/$iterations... "
        local start=$(perl -MTime::HiRes=time -e 'print time')
        eval "$command" > /dev/null 2>&1
        local end=$(perl -MTime::HiRes=time -e 'print time')
        local elapsed=$(echo "$end - $start" | bc)
        times+=($elapsed)
        total=$(echo "$total + $elapsed" | bc)
        echo "${elapsed}s"
    done

    # Calculate statistics
    local mean=$(echo "scale=3; $total / $iterations" | bc)

    # Find min and max
    local min=${times[0]}
    local max=${times[0]}
    for t in "${times[@]}"; do
        if (( $(echo "$t < $min" | bc -l) )); then
            min=$t
        fi
        if (( $(echo "$t > $max" | bc -l) )); then
            max=$t
        fi
    done

    echo -e "${GREEN}  Mean: ${mean}s  Min: ${min}s  Max: ${max}s${NC}"
    echo "$name,$mean,$min,$max" >> "$OUTPUT_DIR/results.csv"
    echo ""
}

# Initialize results CSV
echo "benchmark,mean_seconds,min_seconds,max_seconds" > "$OUTPUT_DIR/results.csv"

echo "=========================================="
echo "Test 1: Sequential BAM Reading (View All)"
echo "=========================================="
echo ""

# samtools view (all records)
run_benchmark "samtools_view_all" "samtools view -c $BAM_FILE"

# biometal (would need a CLI tool, so we'll use a Python script)
cat > "$OUTPUT_DIR/biometal_view_all.py" << 'EOF'
#!/usr/bin/env python3
import sys
import biometal

bam_path = sys.argv[1]
count = 0
for record in biometal.BamReader.from_path(bam_path):
    count += 1
EOF

chmod +x "$OUTPUT_DIR/biometal_view_all.py"
run_benchmark "biometal_view_all" "python3 $OUTPUT_DIR/biometal_view_all.py $BAM_FILE"

echo "=========================================="
echo "Test 2: Indexed Region Query (chr1:1-1000)"
echo "=========================================="
echo ""

# samtools view with region
run_benchmark "samtools_query_small" "samtools view -c $BAM_FILE chr1:1-1000"

# biometal indexed query
cat > "$OUTPUT_DIR/biometal_query_small.py" << 'EOF'
#!/usr/bin/env python3
import sys
import biometal

bam_path = sys.argv[1]
bai_path = sys.argv[2]

index = biometal.BaiIndex.from_path(bai_path)
count = 0
for record in biometal.BamReader.query_region(bam_path, index, "chr1", 1, 1000):
    count += 1
EOF

chmod +x "$OUTPUT_DIR/biometal_query_small.py"
run_benchmark "biometal_query_small" "python3 $OUTPUT_DIR/biometal_query_small.py $BAM_FILE $BAI_FILE"

echo "=========================================="
echo "Test 3: Indexed Region Query (chr1:1-10000)"
echo "=========================================="
echo ""

# samtools view with larger region
run_benchmark "samtools_query_medium" "samtools view -c $BAM_FILE chr1:1-10000"

# biometal indexed query (medium)
cat > "$OUTPUT_DIR/biometal_query_medium.py" << 'EOF'
#!/usr/bin/env python3
import sys
import biometal

bam_path = sys.argv[1]
bai_path = sys.argv[2]

index = biometal.BaiIndex.from_path(bai_path)
count = 0
for record in biometal.BamReader.query_region(bam_path, index, "chr1", 1, 10000):
    count += 1
EOF

chmod +x "$OUTPUT_DIR/biometal_query_medium.py"
run_benchmark "biometal_query_medium" "python3 $OUTPUT_DIR/biometal_query_medium.py $BAM_FILE $BAI_FILE"

echo "=========================================="
echo "Test 4: MAPQ Filtering (Q>=30)"
echo "=========================================="
echo ""

# samtools view with MAPQ filter
run_benchmark "samtools_mapq_filter" "samtools view -c -q 30 $BAM_FILE"

# biometal MAPQ filter
cat > "$OUTPUT_DIR/biometal_mapq_filter.py" << 'EOF'
#!/usr/bin/env python3
import sys
import biometal

bam_path = sys.argv[1]
count = 0
for record in biometal.BamReader.from_path(bam_path):
    if record.is_mapped and record.mapq and record.mapq >= 30:
        count += 1
EOF

chmod +x "$OUTPUT_DIR/biometal_mapq_filter.py"
run_benchmark "biometal_mapq_filter" "python3 $OUTPUT_DIR/biometal_mapq_filter.py $BAM_FILE"

echo "=========================================="
echo "Test 5: Memory Usage"
echo "=========================================="
echo ""

# Memory usage for samtools view
echo -e "${BLUE}samtools view memory:${NC}"
/usr/bin/time -l samtools view $BAM_FILE > /dev/null 2> "$OUTPUT_DIR/samtools_memory.txt"
grep "maximum resident set size" "$OUTPUT_DIR/samtools_memory.txt" | awk '{print "  Peak RSS: " $1/1024/1024 " MB"}'

# Memory usage for biometal
echo -e "${BLUE}biometal memory:${NC}"
/usr/bin/time -l python3 "$OUTPUT_DIR/biometal_view_all.py" $BAM_FILE 2> "$OUTPUT_DIR/biometal_memory.txt"
grep "maximum resident set size" "$OUTPUT_DIR/biometal_memory.txt" | awk '{print "  Peak RSS: " $1/1024/1024 " MB"}'

echo ""
echo "=========================================="
echo "Results Summary"
echo "=========================================="
echo ""

# Parse results and create comparison
python3 - <<EOF
import csv
import sys

with open('$OUTPUT_DIR/results.csv', 'r') as f:
    reader = csv.DictReader(f)
    results = list(reader)

# Group by test type
samtools_results = {r['benchmark']: float(r['mean_seconds']) for r in results if 'samtools' in r['benchmark']}
biometal_results = {r['benchmark']: float(r['mean_seconds']) for r in results if 'biometal' in r['benchmark']}

print("Performance Comparison (Mean Times):")
print("=" * 70)
print(f"{'Test':<30} {'samtools':<12} {'biometal':<12} {'Speedup':>10}")
print("-" * 70)

tests = [
    ("Sequential read (all)", "samtools_view_all", "biometal_view_all"),
    ("Query small region", "samtools_query_small", "biometal_query_small"),
    ("Query medium region", "samtools_query_medium", "biometal_query_medium"),
    ("MAPQ filter (Q>=30)", "samtools_mapq_filter", "biometal_mapq_filter"),
]

for name, st_key, bm_key in tests:
    if st_key in samtools_results and bm_key in biometal_results:
        st_time = samtools_results[st_key]
        bm_time = biometal_results[bm_key]
        speedup = st_time / bm_time if bm_time > 0 else 0

        print(f"{name:<30} {st_time:>10.3f}s {bm_time:>10.3f}s {speedup:>9.2f}×")

print("=" * 70)
EOF

echo ""
echo "Full results saved to: $OUTPUT_DIR/"
echo "  - results.csv: Raw benchmark data"
echo "  - *_memory.txt: Memory usage data"
echo "  - *.py: Benchmark scripts"
echo ""
echo -e "${GREEN}Benchmark complete!${NC}"
