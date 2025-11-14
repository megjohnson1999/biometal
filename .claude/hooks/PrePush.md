# PrePush Hook: Cross-Platform Validation

Validate cross-platform compatibility before pushing code.

## Validation Steps

1. **Test Suite**: Run full test suite (582 tests must pass)
   ```bash
   cargo test --all-features
   ```

2. **Platform-Specific Code Audit**:
   - Check for `#[cfg(target_arch = "aarch64")]` without scalar fallback
   - Check for Mac-specific APIs (Metal, CoreML) without feature flags
   - Verify `#[cfg(not(target_arch = "aarch64"))]` fallback exists

3. **Clippy Lints**: Ensure no warnings
   ```bash
   cargo clippy --all-features -- -D warnings
   ```

4. **Documentation Tests**: Verify examples still work
   ```bash
   cargo test --doc
   ```

5. **Feature Flag Validation**:
   - Test with no features: `cargo test --no-default-features`
   - Test each feature independently
   - Test all features: `cargo test --all-features`

6. **Python Bindings** (if modified):
   ```bash
   cd python_bindings
   maturin develop
   pytest tests/
   ```

## Platform-Specific Checks

### Mac ARM (Primary Platform)
- ✅ NEON intrinsics compile
- ✅ Metal shaders compile (if Phase 1 code)
- ✅ CoreML models load (if Phase 2 code)

### Linux ARM (Graviton - Portability Target)
- ⚠️ Simulated via feature flags (actual testing requires CI/hardware)
- ✅ NEON code compiles with `target_arch = "aarch64"`
- ✅ No Mac-specific APIs (Metal, CoreML) without feature guards

### x86_64 (Fallback Target)
- ✅ Scalar fallbacks compile
- ✅ Tests pass without NEON optimizations
- ✅ Performance degrades gracefully (no panics)

## Failure Actions

If any validation fails:
1. **Report failure** with specific error
2. **Block push** (user must fix before pushing)
3. **Suggest fix**:
   - Missing tests: `cargo test` output
   - Platform issues: Point to cfg patterns in existing code
   - Python issues: Check PyO3 bindings consistency

## Success Actions

If all validations pass:
1. **Report success** with test count
2. **Allow push** to proceed
3. **Suggest next steps**:
   - Update CHANGELOG.md if significant changes
   - Consider updating documentation
   - Run benchmarks if performance-sensitive code changed

## Example Output

```
🔍 PrePush Validation: Cross-Platform Testing

✅ Test suite: 582 tests passed
✅ Clippy: No warnings
✅ Documentation tests: 121 tests passed
✅ Platform audit: All ARM code has scalar fallback
⚠️ Python bindings: Skipped (no changes detected)

✨ Cross-platform validation PASSED

📋 Optional next steps:
  - CHANGELOG.md updated? (if user-facing changes)
  - Benchmarks run? (if performance code changed)
  - Documentation updated? (if API changes)

✅ Safe to push
```

## Configuration

Validation can be customized via environment variables:
- `SKIP_PYTHON_TESTS=1`: Skip Python binding tests (faster for Rust-only changes)
- `SKIP_CLIPPY=1`: Skip Clippy (not recommended)
- `VERBOSE=1`: Show full test output

## Integration with Strategic Pivot

As GPU/Metal/Neural Engine code is added (Phase 1-3):

1. **Metal Shaders** (Phase 1):
   - Verify `.metal` files compile on Mac
   - Verify feature flag `metal` guards Metal APIs
   - Ensure fallback exists (CPU implementation)

2. **CoreML Models** (Phase 2):
   - Verify `.mlmodel` files load on Mac
   - Verify feature flag `coreml` guards CoreML APIs
   - Ensure fallback exists (CPU/ONNX implementation)

3. **CUDA/OpenCL** (if added for Linux portability):
   - Verify CUDA compiles with feature flag
   - Verify OpenCL fallback exists
   - Ensure CPU fallback always available

## Anti-Patterns to Prevent

- ❌ Pushing ARM-only code without x86 fallback
- ❌ Pushing Mac-specific APIs without feature guards
- ❌ Pushing failing tests
- ❌ Pushing code that breaks Python bindings
- ❌ Pushing without running Clippy

---

**Purpose**: Catch cross-platform issues before they reach CI/production. Maintain portability as complexity increases (GPU/Metal/Neural Engine).
