# AfterEdits Hook: Documentation Sync Check

Detect when code changes might require documentation updates.

## Trigger Conditions

After code edits, check if documentation might be affected:

### Public API Changes
If changes touch:
- `pub fn`, `pub struct`, `pub enum` (API surface)
- `src/lib.rs` (public exports)
- Python bindings (`python_bindings/src/lib.rs`)

→ **Action**: Flag that API documentation may need updates

### Performance-Critical Code
If changes touch:
- Functions with `#[inline]` or `#[target_feature]`
- NEON intrinsics (`vld1q_u8`, `vaddq_u8`, etc.)
- BGZF decompression
- BAM/FASTQ/FASTA parsers

→ **Action**: Flag that performance documentation may need updates

### Configuration/Architecture
If changes touch:
- `Cargo.toml` (dependencies, features)
- `.claude/` (agent/hook configuration)
- `benchmarks/` (benchmark methodology)

→ **Action**: Flag that CLAUDE.md or guides may need updates

### Examples/Tutorials
If changes touch:
- `examples/` (example code)
- `notebooks/` (Jupyter tutorials)

→ **Action**: Flag that examples may need validation

## Detection Logic

Analyze edited files and report:

```markdown
📝 Documentation Sync Check

Edited files: <list of files changed>

Potential documentation impacts:
- [ ] API documentation (public API changed)
- [ ] Performance benchmarks (performance code changed)
- [ ] User Guide (usage patterns changed)
- [ ] CLAUDE.md (architecture/configuration changed)
- [ ] Examples (example code may be outdated)
- [ ] Python docs (Python bindings changed)

Recommended actions:
1. <specific doc file to review>
2. <specific section to update>
3. <validation to perform>

Use `/update-docs` to systematically update affected documentation.
```

## Documentation Map

Map code locations to documentation:

### API Documentation
- `src/io/fastq.rs` → `docs/USER_GUIDE.md` (FASTQ section)
- `src/io/fasta.rs` → `docs/USER_GUIDE.md` (FASTA section)
- `src/io/bam/` → `docs/USER_GUIDE.md` (BAM section), `docs/BAM_API.md`
- `src/operations/` → `docs/USER_GUIDE.md` (Operations section)
- `python_bindings/` → `docs/PYTHON.md`

### Performance Documentation
- `src/operations/simd/` → `docs/PERFORMANCE_OPTIMIZATION_GUIDE.md` (NEON section)
- `benchmarks/` → `OPTIMIZATION_RULES.md`, benchmark reports
- BGZF code → `docs/PERFORMANCE_OPTIMIZATION_GUIDE.md` (Decompression section)

### Architecture Documentation
- `src/lib.rs` → `docs/ARCHITECTURE.md`
- `.claude/` → `CLAUDE.md` (Workflow section)
- `Cargo.toml` → `CLAUDE.md` (Dependencies section)

### Tutorial/Example Documentation
- `examples/` → Validate all examples compile and run
- `notebooks/` → Validate all cells execute

## Smart Filtering

Don't flag documentation for:
- Internal implementation details (private functions)
- Test code (`#[cfg(test)]`, `tests/`)
- Documentation files themselves (`.md` files)
- Minor refactoring (variable renames, whitespace)

DO flag documentation for:
- Public API signature changes
- Performance characteristic changes
- Breaking changes
- New features
- Deprecated features

## Example Output

```markdown
📝 Documentation Sync Check

Edited files:
- src/io/bam/parser.rs (public API changed)
- src/operations/simd/base_count.rs (performance code changed)

Potential documentation impacts:
- [x] API documentation (pub fn parse_bam_record signature changed)
- [x] Performance benchmarks (base_count implementation changed)
- [ ] User Guide (usage patterns unchanged)
- [ ] CLAUDE.md (no architectural changes)
- [ ] Examples (example code still valid)
- [ ] Python docs (Python bindings unchanged)

Recommended actions:
1. Update docs/BAM_API.md: Document new parse_bam_record parameters
2. Re-run benchmarks: cargo bench --bench base_counting (N=30)
3. Update OPTIMIZATION_RULES.md: If base_count performance changed significantly
4. Validate examples/bam_parser.rs: Ensure still compiles with new API

Use `/update-docs` to systematically address these items.
```

## Integration with Weekly Progress

If changes are significant (new features, breaking changes):
- **Suggest updating PROJECT_TODOS.md**: Mark tasks complete or adjust timelines
- **Suggest updating CHANGELOG.md**: Document user-facing changes
- **Suggest updating version**: Follow SemVer (patch/minor/major)

## Workflow Integration

After hook detects documentation drift:
1. **Review suggestions**: User decides which docs need updating
2. **Use /update-docs**: Systematic documentation update workflow
3. **Use doc-specialist agent**: For comprehensive documentation rewrites
4. **Validate examples**: Run examples/tests to confirm docs match reality

## Anti-Patterns to Prevent

- ❌ Merging public API changes without updating docs
- ❌ Changing performance characteristics without updating benchmarks
- ❌ Adding features without updating User Guide
- ❌ Breaking examples without validation

---

**Purpose**: Maintain documentation quality as codebase evolves. Prevent documentation drift, especially critical given 142K+ word documentation base.
