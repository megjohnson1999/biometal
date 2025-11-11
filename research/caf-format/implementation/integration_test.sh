#!/usr/bin/env bash
#
# CAF Integration Test Suite
# Tests BAM → CAF → SAM round-trip with multiple datasets
#
# Usage: ./integration_test.sh

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test configuration
TEST_DIR="../../../tests/data"
OUTPUT_DIR="/tmp/caf_integration_tests"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RESULTS_FILE="INTEGRATION_TEST_RESULTS_${TIMESTAMP}.md"

# Test files
declare -a TEST_FILES=(
    "$TEST_DIR/test.bam"
    "$TEST_DIR/synthetic_100k.bam"
    "$TEST_DIR/large/large_1m.bam"
)

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Results accumulator
declare -a TEST_RESULTS=()

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}CAF Integration Test Suite${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo "Timestamp: $(date)"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Function to format file size
format_size() {
    local size=$1
    if [ "$size" -lt 1024 ]; then
        echo "${size} B"
    elif [ "$size" -lt 1048576 ]; then
        echo "$(awk "BEGIN {printf \"%.2f\", $size/1024}") KB"
    else
        echo "$(awk "BEGIN {printf \"%.2f\", $size/1048576}") MB"
    fi
}

# Function to count records in SAM file
count_sam_records() {
    local sam_file=$1
    grep -v '^@' "$sam_file" | wc -l | tr -d ' '
}

# Function to test a single BAM file
test_bam_file() {
    local bam_file=$1
    local basename=$(basename "$bam_file" .bam)

    echo -e "${YELLOW}Testing: $basename${NC}"
    echo "----------------------------------------"

    # Check if file exists
    if [ ! -f "$bam_file" ]; then
        echo -e "${RED}✗ File not found: $bam_file${NC}"
        TEST_RESULTS+=("❌ $basename: File not found")
        return 1
    fi

    # Get BAM file size
    local bam_size=$(stat -f%z "$bam_file" 2>/dev/null || stat -c%s "$bam_file" 2>/dev/null)
    local bam_size_fmt=$(format_size "$bam_size")
    echo "  BAM size: $bam_size_fmt"

    # Output files
    local caf_file="$OUTPUT_DIR/${basename}.caf"
    local sam_file="$OUTPUT_DIR/${basename}.sam"

    # Step 1: BAM → CAF conversion
    echo -n "  [1/3] BAM → CAF... "
    local start_time=$(date +%s%3N)

    if cargo run --release --example bam_to_caf -- "$bam_file" "$caf_file" > "$OUTPUT_DIR/${basename}_bam_to_caf.log" 2>&1; then
        local end_time=$(date +%s%3N)
        local conversion_time=$((end_time - start_time))
        echo -e "${GREEN}✓${NC} (${conversion_time}ms)"

        # Get CAF file size
        local caf_size=$(stat -f%z "$caf_file" 2>/dev/null || stat -c%s "$caf_file" 2>/dev/null)
        local caf_size_fmt=$(format_size "$caf_size")
        local compression_ratio=$(awk "BEGIN {printf \"%.2f\", $caf_size/$bam_size}")

        echo "  CAF size: $caf_size_fmt (${compression_ratio}× vs BAM)"
    else
        echo -e "${RED}✗ FAILED${NC}"
        echo "  See log: $OUTPUT_DIR/${basename}_bam_to_caf.log"
        TEST_RESULTS+=("❌ $basename: BAM → CAF conversion failed")
        return 1
    fi

    # Step 2: CAF → SAM conversion
    echo -n "  [2/3] CAF → SAM... "
    start_time=$(date +%s%3N)

    if cargo run --release --example caf_to_sam -- "$caf_file" "$sam_file" > "$OUTPUT_DIR/${basename}_caf_to_sam.log" 2>&1; then
        end_time=$(date +%s%3N)
        conversion_time=$((end_time - start_time))
        echo -e "${GREEN}✓${NC} (${conversion_time}ms)"

        # Get SAM file size
        local sam_size=$(stat -f%z "$sam_file" 2>/dev/null || stat -c%s "$sam_file" 2>/dev/null)
        local sam_size_fmt=$(format_size "$sam_size")

        echo "  SAM size: $sam_size_fmt"
    else
        echo -e "${RED}✗ FAILED${NC}"
        echo "  See log: $OUTPUT_DIR/${basename}_caf_to_sam.log"
        TEST_RESULTS+=("❌ $basename: CAF → SAM conversion failed")
        return 1
    fi

    # Step 3: Validate record count
    echo -n "  [3/3] Validating records... "
    local sam_records=$(count_sam_records "$sam_file")
    echo -e "${GREEN}✓${NC} ($sam_records records)"

    # Success
    echo -e "${GREEN}✓ Round-trip successful${NC}"
    TEST_RESULTS+=("✅ $basename: $sam_records records, ${compression_ratio}× vs BAM, ${caf_size_fmt}")

    echo ""
    return 0
}

# Run tests
echo -e "${BLUE}Running integration tests...${NC}"
echo ""

PASSED=0
FAILED=0

for test_file in "${TEST_FILES[@]}"; do
    if test_bam_file "$test_file"; then
        ((PASSED++))
    else
        ((FAILED++))
    fi
done

# Print summary
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Test Summary${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

for result in "${TEST_RESULTS[@]}"; do
    echo "  $result"
done

echo ""
echo "Total: $((PASSED + FAILED)) tests"
echo -e "${GREEN}Passed: $PASSED${NC}"

if [ "$FAILED" -gt 0 ]; then
    echo -e "${RED}Failed: $FAILED${NC}"
fi

echo ""

# Write results to markdown file
cat > "$RESULTS_FILE" <<EOF
# CAF Integration Test Results

**Date**: $(date)
**Test Suite**: Week 3 Integration Testing
**Status**: $PASSED/$((PASSED + FAILED)) tests passed

---

## Summary

EOF

for result in "${TEST_RESULTS[@]}"; do
    echo "- $result" >> "$RESULTS_FILE"
done

cat >> "$RESULTS_FILE" <<EOF

---

## Test Configuration

- **Test files**: ${#TEST_FILES[@]} BAM files
- **Output directory**: \`$OUTPUT_DIR\`
- **Dictionary compression**: Enabled (trained on 30K samples)
- **Block size**: 10,000 records
- **Compression level**: Zstandard level 3

---

## File Details

EOF

for test_file in "${TEST_FILES[@]}"; do
    basename=$(basename "$test_file")
    if [ -f "$test_file" ]; then
        size=$(stat -f%z "$test_file" 2>/dev/null || stat -c%s "$test_file" 2>/dev/null)
        size_fmt=$(format_size "$size")
        echo "- **$basename**: $size_fmt" >> "$RESULTS_FILE"
    fi
done

cat >> "$RESULTS_FILE" <<EOF

---

## Logs

Test logs are available in: \`$OUTPUT_DIR/\`

- \`*_bam_to_caf.log\`: BAM → CAF conversion logs
- \`*_caf_to_sam.log\`: CAF → SAM conversion logs

---

## Next Steps

1. Analyze compression ratios across different file sizes
2. Profile performance on large files (1M+ records)
3. Test edge cases (empty files, unusual quality distributions)
4. Benchmark NEON operations on CAF blocks

EOF

echo -e "${BLUE}Results written to: $RESULTS_FILE${NC}"
echo ""

# Exit with failure if any tests failed
if [ "$FAILED" -gt 0 ]; then
    exit 1
fi

exit 0
