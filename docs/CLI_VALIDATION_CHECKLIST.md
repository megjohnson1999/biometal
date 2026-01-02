# CLI Validation Quick Checklist

**⚡ Quick validation checklist for biometal CLI tools**
*Full methodology: see [VALIDATION_FRAMEWORK.md](VALIDATION_FRAMEWORK.md)*

---

## 🚀 Quick Start (5 minutes)

For any new CLI command, run this validation sequence:

### 1. Integration Check (30 seconds)
```bash
# Verify CLI uses library primitives (not manual implementations)
grep -r "use biometal::operations::" src/bin/cli/[command].rs
grep -r "[primitive_function](" src/bin/cli/[command].rs

# ✅ Should see direct library imports and function calls
# ❌ Red flag: manual loops, hardcoded ATCG logic
```

### 2. Basic Test (1 minute)
```bash
# Create simple mixed-case test
echo -e "@test\nATCGatcg\n+\nIIIIIIII" > /tmp/test.fq

# Test command works
cargo run --bin biometal --release -- [command] /tmp/test.fq

# ✅ Should handle mixed case correctly
# ✅ Should show NEON acceleration on ARM64
```

### 3. Edge Cases (2 minutes)
```bash
# Create edge case test
cat > /tmp/edge.fq << 'EOF'
@empty

+

@with_n
ATCGNNNatcg
+
IIIIIIIIIII
@single
A
+
I
EOF

cargo run --bin biometal --release -- [command] /tmp/edge.fq

# ✅ Should handle empty sequences, N's, single bases gracefully
```

### 4. Cross-Validation (2 minutes)
```python
# Quick manual check for critical results
seq = "ATCGatcg"
expected_A = seq.upper().count('A')  # Should match CLI output
expected_total = len([b for b in seq if b.upper() in 'ATCG'])
# Compare with CLI results
```

---

## 🎯 Validation Matrix

| Check | What to Look For | Pass Criteria |
|-------|------------------|---------------|
| **🔗 Integration** | `use biometal::operations::X` | Direct library calls, no manual loops |
| **📝 Case Sensitivity** | Mixed `ATCGatcg` → correct counts | Case-insensitive behavior |
| **⚠️ Edge Cases** | Empty/N/single sequences | Graceful handling, no crashes |
| **⚡ Performance** | NEON acceleration in output | "X.X× NEON acceleration on ARM64" |
| **🎯 Correctness** | Results vs external tools | Matches BioPython/seqtk |

---

## 🚨 Common Issues Quick Fix

### Issue: Wrong Results
```bash
# Check CLI source
grep -A10 -B10 "[primitive_function]" src/bin/cli/[command].rs
# Look for: wrapper logic, array index bugs, manual implementations
```

### Issue: No NEON Claims
```bash
# Ensure ARM64 + release build
uname -m  # Should show "arm64"
cargo run --bin biometal --release  # Use --release flag
```

### Issue: Case Sensitivity Broken
```bash
# Rebuild completely
cargo clean && cargo build --release
# Check test data format (valid FASTQ?)
```

---

## 📋 New Command Template

```markdown
## Validation: [COMMAND_NAME]

### Integration ⚙️
- [ ] Uses `biometal::operations::[FUNCTION]`
- [ ] Direct function calls: `[FUNCTION](&record.sequence)`
- [ ] No manual implementations found

### Correctness ✅
- [ ] Basic test: `ATCGatcg` → expected results
- [ ] Edge cases: empty/N/single bases handled
- [ ] Cross-validation: matches external tool

### Performance 🚀
- [ ] Shows "X.X× NEON acceleration on ARM64"
- [ ] Memory usage constant (~5MB)
- [ ] No performance regression

**Status**: ✅ PASS / ❌ FAIL
**Notes**: [Any issues or observations]
```

---

## 🛠 Quick Test Data

```bash
# Standard test files (copy/paste ready)

# Mixed case test
cat > /tmp/mixed.fq << 'EOF'
@mixed_case
AAAAaaaaTTTTttttCCCCccccGGGGgggg
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

# Edge cases test
cat > /tmp/edge.fq << 'EOF'
@empty_sequence

+

@all_n
NNNNNNNN
+
IIIIIIII
@mixed_n
ATCGNNNatcg
+
IIIIIIIIIII
@single_upper
A
+
I
@single_lower
a
+
I
EOF

# Pattern test
cat > /tmp/pattern.fq << 'EOF'
@pattern_test
ATCGatcgATCGatcg
+
IIIIIIIIIIIIIIII
EOF
```

---

## 🔄 Automated Validation

```bash
#!/bin/bash
# Quick validation script

COMMAND="$1"
echo "🔍 Quick validation: $COMMAND"

# Integration check
echo "⚙️ Checking integration..."
if grep -q "use biometal::operations" src/bin/cli/*.rs; then
    echo "   ✅ Library imports found"
else
    echo "   ❌ No library imports"
    exit 1
fi

# Basic functionality
echo "✅ Testing basic functionality..."
echo -e "@test\nATCGatcg\n+\nIIIIIIII" | cargo run --bin biometal --release -- "$COMMAND"

echo "🎉 Quick validation complete!"
```

---

**Use this checklist for every new CLI command to ensure quality and consistency!**