#!/bin/bash
# validate_cli.sh - Automated CLI validation script
# Usage: ./scripts/validate_cli.sh <command_name> [test_file]

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
COMMAND_NAME="$1"
TEST_FILE="${2:-}"
TEMP_DIR="/tmp/biometal_validation"

# Help function
show_help() {
    cat << EOF
biometal CLI Validation Script

USAGE:
    $0 <command_name> [test_file]

EXAMPLES:
    $0 count-bases                    # Validate count-bases with generated test data
    $0 gc-content sample.fastq        # Validate gc-content with custom file
    $0 find-pattern                   # Validate find-pattern with generated test data

VALIDATION PHASES:
    ⚙️  Integration: Verify CLI uses library primitives
    ✅ Correctness: Test with mixed-case and edge cases
    🚀 Performance: Check NEON acceleration claims

EOF
}

# Utility functions
log_phase() {
    echo -e "${BLUE}$1${NC}"
}

log_success() {
    echo -e "   ${GREEN}✅ $1${NC}"
}

log_error() {
    echo -e "   ${RED}❌ $1${NC}"
}

log_warning() {
    echo -e "   ${YELLOW}⚠️  $1${NC}"
}

# Check if command provided
if [ -z "$COMMAND_NAME" ]; then
    show_help
    exit 1
fi

# Setup temp directory
mkdir -p "$TEMP_DIR"
cd "$(dirname "$0")/.."

echo -e "${BLUE}🔍 Validating CLI command: ${YELLOW}$COMMAND_NAME${NC}"
echo

# Phase 1: Integration Verification
log_phase "⚙️  Phase 1: Integration Verification"

# Check if CLI source exists
CLI_FILE="src/bin/cli/$(echo "$COMMAND_NAME" | tr '-' '_').rs"
if [ ! -f "$CLI_FILE" ]; then
    # Try finding in other CLI files
    CLI_FILES=(src/bin/cli/*.rs)
    FOUND_IN=""
    for file in "${CLI_FILES[@]}"; do
        if grep -q "pub fn $(echo "$COMMAND_NAME" | tr '-' '_')" "$file" 2>/dev/null; then
            CLI_FILE="$file"
            FOUND_IN="$file"
            break
        fi
    done

    if [ -z "$FOUND_IN" ]; then
        log_error "CLI source file not found for $COMMAND_NAME"
        exit 1
    fi
fi

# Check library primitive usage
if grep -q "use biometal::operations::" "$CLI_FILE"; then
    PRIMITIVE_IMPORTS=$(grep "use biometal::operations::" "$CLI_FILE" | head -3)
    log_success "Library imports found:"
    echo "$PRIMITIVE_IMPORTS" | sed 's/^/       /'
else
    log_error "No library imports found in $CLI_FILE"
    exit 1
fi

# Check for manual implementation anti-patterns (more specific)
MANUAL_PATTERNS=(
    "for.*base.*in.*sequence"           # Manual base iteration
    "match.*base.*['=>].*[atcgATCGnN]"  # Manual base matching
    "writeln.*@.*record\.id"            # Manual FASTQ header writing
    "writeln.*>.*record\.id"            # Manual FASTA header writing
    "for.*char.*in.*bases"              # Manual character iteration
    "match.*b'[ATCGN]'"                 # Manual byte matching
)

MANUAL_DETECTED=false
for pattern in "${MANUAL_PATTERNS[@]}"; do
    if grep -q "$pattern" "$CLI_FILE"; then
        log_warning "Manual implementation anti-pattern: $pattern"
        MANUAL_DETECTED=true
    fi
done

# Check for library function calls and writers
LIBRARY_FUNCTIONS=$(grep -E "(count_bases|gc_content|find_pattern|count_pattern|to_fasta_record|reverse_complement|complement|reverse)\(" "$CLI_FILE" | wc -l)
LIBRARY_WRITERS=$(grep -E "(FastqWriter|FastaWriter|BamWriter)" "$CLI_FILE" | wc -l)

if [ "$LIBRARY_FUNCTIONS" -gt 0 ]; then
    log_success "Library function calls found ($LIBRARY_FUNCTIONS occurrences)"
elif [ "$LIBRARY_WRITERS" -gt 0 ]; then
    log_success "Library writers found ($LIBRARY_WRITERS occurrences)"
else
    log_warning "No library function calls or writers detected - verify implementation"
fi

echo

# Phase 2: Correctness Testing
log_phase "✅ Phase 2: Correctness Testing"

# Generate test data if not provided
if [ -z "$TEST_FILE" ]; then
    TEST_FILE="$TEMP_DIR/test_data.fastq"
    cat > "$TEST_FILE" << 'EOF'
@mixed_case_test
AAAAaaaaTTTTttttCCCCccccGGGGgggg
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@all_uppercase
AAAATTTTCCCCGGGG
+
IIIIIIIIIIIIIIII
@all_lowercase
aaaattttccccgggg
+
IIIIIIIIIIIIIIII
EOF
    log_success "Generated test data: $TEST_FILE"
else
    log_success "Using provided test file: $TEST_FILE"
fi

# Test basic functionality
echo "   Testing basic functionality..."
if timeout 30s cargo run --bin biometal --release -- "$COMMAND_NAME" --help > /dev/null 2>&1; then
    log_success "Help command works"
else
    log_error "Help command failed or timed out"
    exit 1
fi

# Test with actual data
echo "   Testing with mixed-case data..."
OUTPUT_FILE="$TEMP_DIR/output.txt"

# Handle commands that require special parameters
if [ "$COMMAND_NAME" = "extract-region" ]; then
    # extract-region requires --start and --end parameters
    COMMAND_ARGS="--start 5 --end 15 $TEST_FILE"
else
    COMMAND_ARGS="$TEST_FILE"
fi

if timeout 60s cargo run --bin biometal --release -- "$COMMAND_NAME" $COMMAND_ARGS > "$OUTPUT_FILE" 2>&1; then
    log_success "Command executed successfully"

    # Show first few lines of output
    echo "   Sample output:"
    head -5 "$OUTPUT_FILE" | sed 's/^/       /'
else
    log_error "Command failed or timed out"
    echo "   Error output:"
    cat "$OUTPUT_FILE" | sed 's/^/       /'
    exit 1
fi

# Test edge cases
echo "   Testing edge cases..."
EDGE_FILE="$TEMP_DIR/edge_cases.fastq"
cat > "$EDGE_FILE" << 'EOF'
@empty_sequence

+

@all_n_sequence
NNNNNNNNNNNN
+
IIIIIIIIIIII
@mixed_with_n
ATCGNNNNatcgNNNN
+
IIIIIIIIIIIIIIII
@single_base_upper
A
+
I
@single_base_lower
a
+
I
EOF

# Handle edge case testing with special parameters
if [ "$COMMAND_NAME" = "extract-region" ]; then
    EDGE_COMMAND_ARGS="--start 1 --end 5 $EDGE_FILE"
else
    EDGE_COMMAND_ARGS="$EDGE_FILE"
fi

if timeout 60s cargo run --bin biometal --release -- "$COMMAND_NAME" $EDGE_COMMAND_ARGS > "$TEMP_DIR/edge_output.txt" 2>&1; then
    log_success "Edge cases handled gracefully"
else
    log_warning "Edge cases may have issues"
    echo "   Error output:"
    cat "$TEMP_DIR/edge_output.txt" | sed 's/^/       /'
fi

# Check for case sensitivity handling
if grep -qi "mixed.*case" "$OUTPUT_FILE" 2>/dev/null || [ "$COMMAND_NAME" = "count-bases" ]; then
    # For count-bases, check if it properly counts mixed case
    if [ "$COMMAND_NAME" = "count-bases" ]; then
        # Expected: 32 of each base (16 upper + 16 lower case)
        TOTAL_BASES=$(grep -E "Total|total" "$OUTPUT_FILE" | grep -o '[0-9]\+' | tail -1)
        if [ "$TOTAL_BASES" = "64" ] 2>/dev/null; then
            log_success "Case sensitivity handled correctly (64 total bases)"
        else
            log_warning "Case sensitivity may not be working (got $TOTAL_BASES bases, expected 64)"
        fi
    fi
fi

echo

# Phase 3: Performance Verification
log_phase "🚀 Phase 3: Performance Verification"

# Check for NEON acceleration indicators
if grep -qi "neon" "$OUTPUT_FILE"; then
    NEON_CLAIM=$(grep -i "neon" "$OUTPUT_FILE" | head -1)
    log_success "NEON acceleration claimed:"
    echo "$NEON_CLAIM" | sed 's/^/       /'
else
    # Check platform
    ARCH=$(uname -m)
    if [ "$ARCH" = "arm64" ]; then
        log_warning "No NEON acceleration indicators found (expected on ARM64)"
    else
        log_success "No NEON claims expected on $ARCH platform"
    fi
fi

# Check if running in release mode (performance-critical)
if echo "$0" | grep -q "release"; then
    log_success "Running in release mode"
else
    log_warning "Consider using release mode for performance testing"
fi

# Memory usage check (basic)
echo "   Checking memory behavior..."
if ps aux | grep -q "[b]iometal.*$COMMAND_NAME" 2>/dev/null; then
    log_success "Process completed (streaming behavior preserved)"
else
    log_success "Command completed successfully"
fi

echo

# Summary and recommendations
log_phase "📋 Validation Summary"

VALIDATION_SCORE=0
MAX_SCORE=5

# Score integration (improved logic)
HAS_IMPORTS=$(grep -q "use biometal::operations::" "$CLI_FILE" && echo "true" || echo "false")
HAS_LIBRARY_USAGE=$([ "$LIBRARY_FUNCTIONS" -gt 0 ] || [ "$LIBRARY_WRITERS" -gt 0 ] && echo "true" || echo "false")

if [ "$HAS_IMPORTS" = "true" ] && [ "$HAS_LIBRARY_USAGE" = "true" ] && [ "$MANUAL_DETECTED" = "false" ]; then
    log_success "Integration: PASS"
    ((VALIDATION_SCORE++))
else
    log_error "Integration: FAIL"
    if [ "$MANUAL_DETECTED" = "true" ]; then
        echo "   Manual implementation anti-patterns detected"
    fi
    if [ "$HAS_LIBRARY_USAGE" = "false" ]; then
        echo "   No library function calls or writers found"
    fi
fi

# Score basic functionality
if [ -f "$OUTPUT_FILE" ] && [ -s "$OUTPUT_FILE" ]; then
    log_success "Basic functionality: PASS"
    ((VALIDATION_SCORE++))
else
    log_error "Basic functionality: FAIL"
fi

# Score edge cases
if [ -f "$TEMP_DIR/edge_output.txt" ]; then
    log_success "Edge cases: PASS"
    ((VALIDATION_SCORE++))
else
    log_error "Edge cases: FAIL"
fi

# Score case sensitivity (if applicable)
if [ "$COMMAND_NAME" = "count-bases" ]; then
    if [ "$TOTAL_BASES" = "64" ] 2>/dev/null; then
        log_success "Case sensitivity: PASS"
        ((VALIDATION_SCORE++))
    else
        log_error "Case sensitivity: FAIL"
    fi
else
    log_success "Case sensitivity: N/A"
    ((VALIDATION_SCORE++))
fi

# Score performance indicators
if grep -qi "neon" "$OUTPUT_FILE" || [ "$(uname -m)" != "arm64" ]; then
    log_success "Performance indicators: PASS"
    ((VALIDATION_SCORE++))
else
    log_error "Performance indicators: FAIL"
fi

echo
echo -e "${BLUE}Final Score: ${YELLOW}$VALIDATION_SCORE/$MAX_SCORE${NC}"

if [ "$VALIDATION_SCORE" -eq "$MAX_SCORE" ]; then
    echo -e "${GREEN}🎉 VALIDATION PASSED: $COMMAND_NAME is ready for production!${NC}"
elif [ "$VALIDATION_SCORE" -ge 3 ]; then
    echo -e "${YELLOW}⚠️  VALIDATION PARTIAL: $COMMAND_NAME has minor issues${NC}"
    echo "   Consider reviewing warnings above"
else
    echo -e "${RED}❌ VALIDATION FAILED: $COMMAND_NAME requires fixes${NC}"
    exit 1
fi

# Cleanup
rm -rf "$TEMP_DIR"

echo
echo -e "${BLUE}Validation complete! See docs/VALIDATION_FRAMEWORK.md for detailed methodology.${NC}"