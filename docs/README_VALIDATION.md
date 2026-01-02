# biometal CLI Validation System

**Systematic validation framework for CLI tools and primitives**

## Quick Start

### Validate a CLI Command (30 seconds)
```bash
# Automated validation
./scripts/validate_cli.sh count-bases

# Manual quick check
./scripts/validate_cli.sh gc-content my_data.fastq
```

### Expected Output
```
🔍 Validating CLI command: count-bases

⚙️  Phase 1: Integration Verification
   ✅ Library imports found:
       use biometal::operations::count_bases;
   ✅ Library function calls found (3 occurrences)

✅ Phase 2: Correctness Testing
   ✅ Generated test data: /tmp/biometal_validation/test_data.fastq
   ✅ Help command works
   ✅ Command executed successfully
   ✅ Edge cases handled gracefully
   ✅ Case sensitivity handled correctly (64 total bases)

🚀 Phase 3: Performance Verification
   ✅ NEON acceleration claimed:
       Base Frequencies (16.7× NEON acceleration on ARM64)
   ✅ Performance indicators: PASS

📋 Validation Summary
   ✅ Integration: PASS
   ✅ Basic functionality: PASS
   ✅ Edge cases: PASS
   ✅ Case sensitivity: PASS
   ✅ Performance indicators: PASS

Final Score: 5/5
🎉 VALIDATION PASSED: count-bases is ready for production!
```

## Documentation

| Document | Purpose | When to Use |
|----------|---------|-------------|
| **[VALIDATION_FRAMEWORK.md](VALIDATION_FRAMEWORK.md)** | Complete methodology | Implementing new primitives |
| **[CLI_VALIDATION_CHECKLIST.md](CLI_VALIDATION_CHECKLIST.md)** | Quick reference | Daily development |
| **[validate_cli.sh](../scripts/validate_cli.sh)** | Automation script | Continuous validation |

## Validation Philosophy

**Correctness → Memory Efficiency → Speed**

1. **Correctness First**: Ensure accurate results with real-world genomic data
2. **Memory Efficiency**: Maintain constant ~5MB usage (streaming architecture)
3. **Speed Last**: Verify NEON acceleration claims are preserved

## Common Use Cases

### New CLI Command
```bash
# After implementing new command
./scripts/validate_cli.sh my-new-command

# Check specific aspects
grep "use biometal::operations" src/bin/cli/my_command.rs
```

### After Library Changes
```bash
# Validate all commands still work
for cmd in count-bases gc-content find-pattern count-pattern; do
    ./scripts/validate_cli.sh "$cmd"
done
```

### Before Release
```bash
# Full regression testing
cargo test
./scripts/validate_cli.sh count-bases
./scripts/validate_cli.sh gc-content
./scripts/validate_cli.sh find-pattern
./scripts/validate_cli.sh count-pattern
```

## Integration with Development

### Git Hooks (Optional)
Add to `.git/hooks/pre-commit`:
```bash
#!/bin/bash
# Quick validation before commit
./scripts/validate_cli.sh count-bases > /dev/null
echo "✅ CLI validation passed"
```

### CI/CD Integration
```yaml
# .github/workflows/validate.yml
- name: Validate CLI Tools
  run: |
    ./scripts/validate_cli.sh count-bases
    ./scripts/validate_cli.sh gc-content
```

## Troubleshooting

### Common Issues

| Issue | Cause | Solution |
|-------|-------|----------|
| "No library imports found" | Manual implementation | Use `biometal::operations::*` functions |
| "Case sensitivity: FAIL" | Old binary cached | `cargo clean && cargo build --release` |
| "No NEON claims" | Wrong platform/build | Use ARM64 + `--release` flag |
| "Command failed" | Malformed input | Check FASTQ/FASTA format |

### Debug Mode
```bash
# Run with debug output
RUST_LOG=debug ./scripts/validate_cli.sh count-bases

# Manual testing
echo -e "@test\nATCGatcg\n+\nIIIIIIII" | cargo run --bin biometal --release -- count-bases
```

## Contributing

### Adding New Validation Patterns
1. Update `validate_cli.sh` with new checks
2. Document in `VALIDATION_FRAMEWORK.md`
3. Test with existing commands
4. Update checklist

### Reporting Issues
Include validation output:
```bash
./scripts/validate_cli.sh [command] > validation_report.txt 2>&1
# Attach validation_report.txt to issue
```

---

**This validation system ensures biometal CLI tools are robust, correct, and performant for production genomic workflows.** 🧬✨