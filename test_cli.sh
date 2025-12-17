#!/bin/bash
# Biometal CLI Test Suite
# Tests all implemented CLI commands for functionality and error handling

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test counters
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Test data files
TEST_FASTA="/tmp/test_cli.fa"
TEST_FASTQ="/tmp/test_cli.fq"
TEST_MIXED="/tmp/test_mixed.fq"

# Create test data
setup_test_data() {
    echo -e "${YELLOW}📁 Setting up test data...${NC}"

    # Create FASTA test file
    cat > $TEST_FASTA << 'EOF'
>seq1 Simple sequence
ACGTACGTACGT
>seq2 GC-rich sequence
GGGGCCCCGGGGCCCC
>seq3 AT-rich sequence
AAAATTTTAAAATTTT
EOF

    # Create FASTQ test file
    cat > $TEST_FASTQ << 'EOF'
@read1 High quality
ACGTACGTACGT
+
IIIIIIIIIIII
@read2 Low quality
TTTTTTTTTTTT
+
""""""""""""
@read3 Mixed sequence
AAAACCCCGGGGTTTT
+
IIIIHHHHGGGGFFFF
EOF

    # Create FASTQ with mixed quality
    cat > $TEST_MIXED << 'EOF'
@mixed1
ACGTACGTACGTNNNN
+
IIIIHHHHGGGG####
@mixed2
TTTTTTTTTTTTAAAA
+
GGGGFFFFEEEE!!!!
EOF

    echo "✅ Test data created"
}

echo -e "${YELLOW}🧬 Biometal CLI Test Suite${NC}"
echo "=============================="

# Setup test data
setup_test_data

# Test function
run_test() {
    local test_name="$1"
    local command="$2"
    local expected_pattern="$3"
    local should_fail="${4:-false}"

    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    echo -n "Testing $test_name... "

    if [ "$should_fail" = "true" ]; then
        # Test should fail
        if $command >/dev/null 2>&1; then
            echo -e "${RED}FAIL${NC} (expected failure but command succeeded)"
            FAILED_TESTS=$((FAILED_TESTS + 1))
            return 1
        else
            echo -e "${GREEN}PASS${NC}"
            PASSED_TESTS=$((PASSED_TESTS + 1))
            return 0
        fi
    else
        # Test should succeed
        if output=$($command 2>&1); then
            if [[ -z "$expected_pattern" ]] || echo "$output" | grep -q "$expected_pattern"; then
                echo -e "${GREEN}PASS${NC}"
                PASSED_TESTS=$((PASSED_TESTS + 1))
                return 0
            else
                echo -e "${RED}FAIL${NC} (pattern '$expected_pattern' not found)"
                echo "  Output: $output"
                FAILED_TESTS=$((FAILED_TESTS + 1))
                return 1
            fi
        else
            echo -e "${RED}FAIL${NC} (command failed)"
            echo "  Error: $output"
            FAILED_TESTS=$((FAILED_TESTS + 1))
            return 1
        fi
    fi
}

# Core Statistics Tests
echo -e "${YELLOW}📊 Testing Core Statistics (Tier 1)${NC}"
run_test "count-bases FASTA" "./target/release/biometal count-bases $TEST_FASTA" "NEON acceleration"
run_test "count-bases FASTQ" "./target/release/biometal count-bases $TEST_FASTQ" "Total"
run_test "gc-content calculation" "./target/release/biometal gc-content $TEST_FASTA" "GC content"
run_test "mean-quality calculation" "./target/release/biometal mean-quality $TEST_FASTQ" "Mean quality"
run_test "complexity-score calculation" "./target/release/biometal complexity-score $TEST_FASTA" "Mean complexity"

# Sequence Operations Tests
echo -e "${YELLOW}🔄 Testing Sequence Operations (Tier 1)${NC}"
run_test "reverse-complement" "./target/release/biometal reverse-complement $TEST_FASTA" "seq1"
run_test "complement only" "./target/release/biometal complement $TEST_FASTA" "seq1"
run_test "reverse only" "./target/release/biometal reverse $TEST_FASTA" "seq1"
run_test "validate-dna" "./target/release/biometal validate-dna $TEST_FASTA" "valid DNA"

# Quality Management Tests
echo -e "${YELLOW}⚡ Testing Quality Management (Tier 2A)${NC}"
run_test "trim-quality default" "./target/release/biometal trim-quality $TEST_FASTQ" "@read"
run_test "trim-quality with threshold" "./target/release/biometal trim-quality --threshold 30 $TEST_FASTQ" "@read"
run_test "mask-low-quality" "./target/release/biometal mask-low-quality --threshold 20 $TEST_FASTQ" "NNNN"
run_test "extract-region" "./target/release/biometal extract-region --start 2 --end 8 $TEST_FASTA" "GTACGT"

# Format Conversion Tests
echo -e "${YELLOW}🔄 Testing Format Conversion (Tier 2B)${NC}"
run_test "fastq-to-fasta conversion" "./target/release/biometal fastq-to-fasta $TEST_FASTQ" ">read1"
run_test "count-reads FASTQ" "./target/release/biometal count-reads $TEST_FASTQ" "3 reads"
run_test "count-reads FASTA" "./target/release/biometal count-reads $TEST_FASTA" "3 reads"

# Pattern Matching Tests
echo -e "${YELLOW}🔍 Testing Pattern Matching (Tier 2C)${NC}"
run_test "find-pattern basic" "./target/release/biometal find-pattern --pattern ACGT $TEST_FASTQ" "read1"
run_test "find-pattern all occurrences" "./target/release/biometal find-pattern --pattern TTTT --all $TEST_FASTQ" "read"
run_test "count-pattern" "./target/release/biometal count-pattern --pattern AAAA $TEST_FASTQ" "NEON acceleration"
run_test "find-adapters" "./target/release/biometal find-adapters $TEST_FASTQ" "Adapter Contamination Summary"

# Error Handling Tests
echo -e "${YELLOW}❌ Testing Error Handling${NC}"
run_test "nonexistent file" "./target/release/biometal count-bases /tmp/nonexistent.fa" "" true
run_test "invalid command" "./target/release/biometal invalid-command" "" true
run_test "missing pattern argument" "./target/release/biometal find-pattern $TEST_FASTQ" "" true

# Help and Version Tests
echo -e "${YELLOW}❓ Testing Help and Version${NC}"
run_test "main help" "./target/release/biometal --help" "Core Statistics"
run_test "version" "./target/release/biometal --version" "biometal"
run_test "command help" "./target/release/biometal count-bases --help" "Count base frequencies"

# Output format tests
echo -e "${YELLOW}📄 Testing Output Formats${NC}"
run_test "find-pattern positions" "./target/release/biometal find-pattern --pattern ACGT --format positions $TEST_FASTQ" "read1"
run_test "find-pattern sequences" "./target/release/biometal find-pattern --pattern TTTT --format sequences $TEST_FASTQ" "TTTT"

# Edge case tests
echo -e "${YELLOW}🔬 Testing Edge Cases${NC}"
run_test "empty pattern" "./target/release/biometal find-pattern --pattern \"\" $TEST_FASTQ" ""
run_test "very short region extract" "./target/release/biometal extract-region --start 1 --end 2 $TEST_FASTA" "C"

# Summary
echo ""
echo -e "${YELLOW}📊 Test Results Summary${NC}"
echo "======================="
echo -e "Total tests: $TOTAL_TESTS"
echo -e "${GREEN}Passed: $PASSED_TESTS${NC}"
if [ $FAILED_TESTS -gt 0 ]; then
    echo -e "${RED}Failed: $FAILED_TESTS${NC}"
else
    echo -e "${GREEN}Failed: $FAILED_TESTS${NC}"
fi

if [ $FAILED_TESTS -eq 0 ]; then
    echo -e "${GREEN}🎉 All tests passed!${NC}"
    echo -e "${GREEN}✅ CLI implementation is ready for production${NC}"
    exit 0
else
    echo -e "${RED}❌ Some tests failed${NC}"
    exit 1
fi

# Cleanup
cleanup() {
    rm -f $TEST_FASTA $TEST_FASTQ $TEST_MIXED
}
trap cleanup EXIT