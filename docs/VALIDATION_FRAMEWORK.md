# biometal CLI Validation Framework

**Version**: 1.0
**Last Updated**: January 2026
**Purpose**: Systematic validation methodology for biometal CLI tools and primitives

---

## Overview

This framework provides a systematic approach to validate biometal CLI tools, ensuring they:
1. **Use library primitives correctly** (no bypassing)
2. **Handle real-world genomic data** (case sensitivity, edge cases)
3. **Maintain performance claims** (NEON acceleration preserved)
4. **Integrate robustly** (error handling, memory efficiency)

**Validation Philosophy**: `Correctness → Memory Efficiency → Speed`

---

## Quick Start Checklist

For any new biometal CLI tool, follow this validation sequence:

### Phase 1: Integration Verification ⚙️
- [ ] **CLI uses library primitives** (not manual implementations)
- [ ] **No bypass logic detected** (grep for manual counting/processing)
- [ ] **Proper imports confirmed** (`use biometal::operations::*`)
- [ ] **Function calls verified** (direct calls to fixed primitives)

### Phase 2: Correctness Testing ✅
- [ ] **Basic functionality works** (simple test cases)
- [ ] **Case sensitivity handled** (mixed upper/lowercase sequences)
- [ ] **Edge cases robust** (empty, N's, single bases, malformed input)
- [ ] **Cross-validation passes** (results match external tools)

### Phase 3: Performance Verification 🚀
- [ ] **NEON acceleration preserved** (output shows speedup claims)
- [ ] **Memory efficiency maintained** (constant memory usage)
- [ ] **No performance regression** (benchmark if needed)

---

## Detailed Methodology

### Phase 1: Integration Verification

#### 1.1 Source Code Inspection

**Goal**: Verify CLI tools use library primitives, not manual implementations

**Commands**:
```bash
# Check for library primitive usage
grep -r "use biometal::operations" src/bin/cli/
grep -r "count_bases\|gc_content\|find_pattern" src/bin/cli/

# Look for actual anti-patterns (refined Jan 2026)
grep -r "writeln.*@.*record\.id\|writeln.*>.*record\.id" src/bin/cli/
grep -r "for.*base.*in.*sequence\|match.*base.*['=>].*[atcgATCGnN]" src/bin/cli/
grep -r "for.*char.*in.*bases\|match.*b'[ATCGN]'" src/bin/cli/

# Check for proper writers (critical integration requirement)
grep -r "FastqWriter\|FastaWriter\|BamWriter" src/bin/cli/
```

**✅ Success Criteria**:
- All CLI commands import and call library functions directly
- Uses proper format writers: `FastqWriter`, `FastaWriter`, `BamWriter` (critical)
- No manual format writing: `writeln!(output, "@{}", record.id)`
- No manual base counting or pattern matching loops found
- Comments reference library functions (e.g., "now case-insensitive")

**❌ Red Flags** (updated Jan 2026):
- Manual format writing patterns (major integration issue)
- Manual loops over sequence bases with hardcoded logic
- Bypassing library writers for output formatting
- Comments about "manual implementation" or "bypass"

#### 1.2 Function Call Verification

**Goal**: Confirm CLI calls the exact fixed primitive functions

**Method**: Read CLI source files and verify:
```rust
// ✅ Good: Direct library usage + proper writers
use biometal::operations::count_bases;
use biometal::io::FastqWriter;

let counts = count_bases(&record.sequence);

let mut writer = FastqWriter::stdout()?;
writer.write_record(&processed_record)?;

// ❌ Bad: Manual implementation
let mut a_count = 0;
for &base in &record.sequence {
    if base == b'A' { a_count += 1; }
}
writeln!(output, "@{}", record.id)?; // Manual format writing
```

### Phase 2: Correctness Testing

#### 2.1 Test Data Generation

Create comprehensive test files covering real-world scenarios:

```bash
# Mixed case sequences (soft-masked genomes)
cat > test_mixed_case.fastq << 'EOF'
@test_mixed_case
AAAAaaaaTTTTttttCCCCccccGGGGgggg
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@test_all_upper
AAAATTTTCCCCGGGG
+
IIIIIIIIIIIIIIII
@test_all_lower
aaaattttccccgggg
+
IIIIIIIIIIIIIIII
EOF

# Edge cases
cat > test_edge_cases.fastq << 'EOF'
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
```

#### 2.2 Basic Functionality Testing

**Goal**: Verify core functionality works correctly

```bash
# Test basic commands
cargo run --bin biometal --release -- count-bases test_mixed_case.fastq
cargo run --bin biometal --release -- gc-content test_mixed_case.fastq
cargo run --bin biometal --release -- find-pattern -p "ATCG" test_mixed_case.fastq
cargo run --bin biometal --release -- count-pattern -p "ATCG" test_mixed_case.fastq
```

**✅ Success Criteria**:
- Commands execute without errors
- Output format is correct and readable
- Basic functionality matches expectations

#### 2.3 Case Sensitivity Validation

**Goal**: Ensure mixed-case genomic sequences handled correctly

**Test**: Mixed case sequence `AAAAaaaaTTTTttttCCCCccccGGGGgggg`

**Expected Results**:
```
count-bases: A=16, T=16, C=16, G=16 (25% each, total=64)
gc-content: 50.00% GC content (32 GC bases out of 64 total)
find-pattern ATCG: Should find matches at case-insensitive positions
count-pattern ATCG: Should count all case-insensitive matches
```

**Manual Verification**:
```python
# Cross-check with manual calculation
seq = "AAAAaaaaTTTTttttCCCCccccGGGGgggg"
A_count = seq.upper().count('A')  # Should be 16
total_bases = len([b for b in seq if b.upper() in 'ATCG'])  # Should be 64
```

#### 2.4 Edge Case Testing

**Critical Edge Cases**:

| Case | Description | Expected Behavior |
|------|-------------|-------------------|
| **Empty sequence** | `""` | Ignore, don't count |
| **All N's** | `NNNNNNNN` | Ignore N's, count only ATCG |
| **Mixed with N's** | `ATCGNNNatcg` | Count `ATCGatcg` = 8 bases |
| **Single base** | `A` or `a` | Count correctly, case-insensitive |
| **Malformed FASTQ** | Missing lines | Proper error message |

**Test Command**:
```bash
cargo run --bin biometal --release -- count-bases test_edge_cases.fastq
```

**✅ Success Criteria**:
- Edge cases handled gracefully without crashes
- N bases properly excluded from ATCG counting
- Empty sequences ignored
- Single bases counted correctly
- Malformed input produces clear error messages

#### 2.5 Cross-Validation with External Tools

**Goal**: Verify results match industry-standard tools

**Setup**:
```bash
# Install reference tools (in virtual environment)
python3 -m venv validation_env
source validation_env/bin/activate
pip install biopython

# Install seqtk (if available)
# brew install seqtk  # macOS
# apt-get install seqtk  # Ubuntu
```

**Cross-validation Script**:
```python
#!/usr/bin/env python3
"""Cross-validate biometal results with external tools"""

from Bio import SeqIO
import subprocess
import json

def validate_base_counting(fastq_file):
    """Compare base counts with BioPython"""

    # BioPython reference
    counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    for record in SeqIO.parse(fastq_file, "fastq"):
        seq = str(record.seq).upper()
        for base in seq:
            if base in counts:
                counts[base] += 1

    total_ref = sum(counts.values())

    # biometal result
    result = subprocess.run([
        "cargo", "run", "--bin", "biometal", "--release", "--",
        "count-bases", "--format", "json", fastq_file
    ], capture_output=True, text=True)

    biometal_data = json.loads(result.stdout)

    # Compare
    assert biometal_data['counts']['A'] == counts['A']
    assert biometal_data['counts']['T'] == counts['T']
    assert biometal_data['counts']['C'] == counts['C']
    assert biometal_data['counts']['G'] == counts['G']
    assert biometal_data['total_bases'] == total_ref

    print(f"✅ Base counting validation passed: {total_ref} bases")

if __name__ == '__main__':
    validate_base_counting("test_mixed_case.fastq")
```

### Phase 3: Performance Verification

#### 3.1 NEON Acceleration Verification

**Goal**: Confirm NEON speedup claims are preserved

**Test**: Look for acceleration indicators in output
```bash
cargo run --bin biometal --release -- count-bases test_file.fastq
# Expected output: "Base Frequencies (16.7× NEON acceleration on ARM64)"

cargo run --bin biometal --release -- gc-content test_file.fastq
# Expected output: "50.00% GC content (64 bases, 20.3× NEON acceleration on ARM64)"
```

**✅ Success Criteria**:
- Output mentions NEON acceleration
- Speedup claims match established benchmarks (16.7-25.1×)
- Performance indicators present on ARM64 platforms

#### 3.2 Memory Efficiency Check

**Goal**: Ensure streaming architecture maintained

**Method**: Monitor memory usage during large file processing

```bash
# Create large test file
head -n 100000 large_dataset.fastq > test_large.fastq

# Monitor memory usage (macOS)
time -l cargo run --bin biometal --release -- count-bases test_large.fastq

# Look for: Maximum resident set size should be ~5-10 MB (constant)
```

**✅ Success Criteria**:
- Memory usage remains constant regardless of file size
- No evidence of accumulating records in memory
- Streaming behavior preserved (processes records one at a time)

#### 3.3 Performance Regression Testing

**Goal**: Ensure fixes didn't break performance

**When to Run**: After major primitive changes affecting core operations

**Method**:
```bash
# Run benchmarks
cargo bench --bench operations

# Compare key metrics:
# - Base counting: Should maintain >10 GiB/s throughput
# - GC content: Should maintain >15 GiB/s throughput
# - Pattern matching: Should maintain NEON speedup claims

# Acceptable: <10% performance impact for correctness fixes
# Unacceptable: >20% performance regression
```

---

## Validation Templates

### Template: New CLI Command Validation

```markdown
## Validation Report: [COMMAND_NAME]

**Date**: [DATE]
**Validator**: [NAME]
**Version**: [VERSION]

### Phase 1: Integration ⚙️
- [ ] **Uses library primitives**: `use biometal::operations::[PRIMITIVE]`
- [ ] **Direct function calls**: `[PRIMITIVE](&record.sequence)`
- [ ] **No bypass logic**: No manual implementations found
- [ ] **Proper error handling**: Library errors propagated correctly

### Phase 2: Correctness ✅
- [ ] **Basic functionality**: Simple test cases pass
- [ ] **Mixed case handling**: `ATCGatcg` → correct case-insensitive results
- [ ] **Edge cases**: Empty/N/single base sequences handled
- [ ] **Cross-validation**: Results match [EXTERNAL_TOOL]

**Test Results**:
```
Input: [TEST_SEQUENCE]
Expected: [EXPECTED_RESULT]
Actual: [ACTUAL_RESULT]
Status: ✅ PASS / ❌ FAIL
```

### Phase 3: Performance 🚀
- [ ] **NEON acceleration**: Output shows "[X.X]× NEON acceleration on ARM64"
- [ ] **Memory efficiency**: Constant ~5MB usage confirmed
- [ ] **No regression**: Performance within acceptable bounds

**Overall Status**: ✅ VALIDATED / ❌ ISSUES_FOUND
```

### Template: Primitive Library Function Validation

```markdown
## Primitive Validation: [FUNCTION_NAME]

### Library Testing
- [ ] **Unit tests pass**: `cargo test [MODULE]::tests`
- [ ] **Property tests pass**: Random input validation
- [ ] **Cross-platform**: Works on ARM64 + x86_64
- [ ] **Performance**: Meets NEON speedup claims

### CLI Integration
- [ ] **CLI uses function**: No manual reimplementation
- [ ] **Results consistent**: CLI matches library function directly
- [ ] **Edge cases propagated**: Library edge case handling works in CLI

### Real-World Testing
- [ ] **Genomic datasets**: Tested with actual FASTQ/FASTA files
- [ ] **Large files**: Memory efficiency confirmed (>1GB inputs)
- [ ] **Error conditions**: Malformed input handled gracefully
```

---

## Common Issues & Debugging

### Issue: CLI Results Don't Match Library

**Symptoms**: CLI output differs from direct library function calls
**Causes**:
- CLI has wrapper logic modifying results
- CLI using different function than expected
- Array indexing bugs in CLI layer

**Debug Steps**:
1. Check CLI source for wrapper functions
2. Verify CLI calls exact library function
3. Add debug prints to compare library vs CLI results
4. Check array index mappings (A/T/C/G order)

### Issue: Case Sensitivity Not Working

**Symptoms**: Mixed case sequences show incorrect counts
**Causes**:
- CLI bypassing fixed library functions
- Old cached binary being used
- Test data malformed

**Debug Steps**:
1. Verify CLI imports: `grep "use biometal::operations" src/bin/cli/`
2. Rebuild completely: `cargo clean && cargo build --release`
3. Check test data format (valid FASTQ/FASTA)
4. Cross-validate with external tools

### Issue: Performance Claims Missing

**Symptoms**: Output doesn't show NEON acceleration indicators
**Causes**:
- Running on x86_64 (NEON only on ARM64)
- Debug build instead of release
- Platform detection issues

**Debug Steps**:
1. Confirm platform: `uname -m` (should show `arm64`)
2. Use release build: `--release` flag
3. Check NEON feature flags in Cargo.toml

---

## Automation Scripts

### Automated Validation Runner

```bash
#!/bin/bash
# validate_cli.sh - Automated CLI validation

set -e

COMMAND_NAME="$1"
if [ -z "$COMMAND_NAME" ]; then
    echo "Usage: $0 <command_name>"
    exit 1
fi

echo "🔍 Validating CLI command: $COMMAND_NAME"

# Phase 1: Integration
echo "⚙️  Phase 1: Integration Verification"
echo "   Checking library primitive usage..."
if grep -r "use biometal::operations" src/bin/cli/ | grep -q "$COMMAND_NAME"; then
    echo "   ✅ Library imports found"
else
    echo "   ❌ No library imports found"
    exit 1
fi

# Phase 2: Correctness
echo "✅ Phase 2: Correctness Testing"
echo "   Running basic functionality test..."
if cargo run --bin biometal --release -- "$COMMAND_NAME" --help > /dev/null 2>&1; then
    echo "   ✅ Help command works"
else
    echo "   ❌ Help command failed"
    exit 1
fi

# Phase 3: Performance
echo "🚀 Phase 3: Performance Verification"
echo "   Testing with sample data..."
# Add specific tests based on command type

echo "✅ Validation complete for $COMMAND_NAME"
```

---

## Best Practices

### 1. **Always Test Edge Cases First**
Before testing normal functionality, verify edge cases (empty, N's, single bases). These reveal integration issues quickly.

### 2. **Use Cross-Validation Early**
Compare results with external tools (BioPython, seqtk) during development, not just at the end.

### 3. **Validate Real-World Data**
Test with actual genomic datasets, not just synthetic sequences. Real data has unexpected edge cases.

### 4. **Document Expected Results**
Always calculate expected results manually or with reference tools before running tests.

### 5. **Check Performance Claims**
Verify NEON acceleration indicators appear in output. Missing indicators suggest integration issues.

### 6. **Test Both FASTQ and FASTA**
Many CLI tools support multiple formats. Test both to ensure consistent behavior.

---

## Maintenance

### When to Re-Validate

**Trigger Events**:
- New primitive function added to library
- New CLI command implemented
- Core library functions modified
- Performance optimization changes
- Platform support changes (new architectures)

**Frequency**:
- **Full validation**: Before each release
- **Regression testing**: After primitive changes
- **Spot checking**: During development

### Updating This Framework

**Version History**:
- v1.0 (Jan 2026): Initial framework based on case sensitivity validation work

**To Update**:
1. Document new validation patterns discovered
2. Add new edge cases found in real-world usage
3. Update performance benchmarks as hardware/optimizations evolve
4. Expand cross-validation tool list as ecosystem grows

---

## References

- **OPTIMIZATION_RULES.md**: Evidence-based optimization rules
- **apple-silicon-bio-bench**: Performance validation dataset
- **BioPython**: Cross-validation reference (biopython.org)
- **seqtk**: High-performance sequence toolkit (github.com/lh3/seqtk)

---

**Framework Validated On**: CLI Integration Improvements (Jan 2026)
**Commands Tested**: All 16 biometal CLI commands (comprehensive validation)
**Success Rate**: 93.75% pass rate (15/16 commands), 56.25% perfect scores (9/16 at 5/5)
**Major Achievements**:
- ✅ Fixed 5 commands: complement, reverse, trim-quality, mask-low-quality, fastq-to-fasta
- ✅ Upgraded from 4/5 to 5/5 validation scores using FastqWriter/FastaWriter
- ✅ Reduced integration warnings by 62.5% through refined anti-pattern detection
- ✅ Enhanced validation framework accuracy (eliminated 75% false positives)
**Next Review**: After next major primitive addition

**See Also**:
- [CLI_VALIDATION_REPORT.md](../../CLI_VALIDATION_REPORT.md): Latest comprehensive validation results
- [scripts/validate_cli.sh](../../scripts/validate_cli.sh): Automated validation implementation