---
description: Automated release workflow (version bump, changelog, git tag, publish)
---

Execute complete release workflow following SemVer and best practices.

Steps:
1. **Version Analysis**:
   - Analyze changes since last release
   - Determine version bump (patch/minor/major)
   - Suggest version based on SemVer rules

2. **Pre-Release Validation**:
   - Run full test suite (582 tests)
   - Run all benchmarks (N=30)
   - Build release binary
   - Test Python bindings

3. **Changelog Update**:
   - Generate changelog from git log
   - Categorize changes (features, fixes, performance, docs)
   - Format according to Keep a Changelog

4. **Version Bump**:
   - Update Cargo.toml
   - Update python_bindings/Cargo.toml
   - Update CLAUDE.md version references
   - Update documentation

5. **Git Tag**:
   - Create annotated tag (v1.X.Y)
   - Push tag to remote

6. **Publish**:
   - Publish to crates.io: cargo publish
   - Publish to PyPI: maturin publish (from python_bindings/)

7. **Post-Release**:
   - Create GitHub release with changelog
   - Update documentation site (if applicable)
   - Announce on social media (draft posts)

Arguments: <version> (e.g., "patch", "minor", "major", or explicit "1.8.0")

Example: /release minor
Example: /release 1.8.0

Expected time: 15-30 minutes (including builds and validation)

Generates: Git tag, crates.io release, PyPI release, GitHub release
