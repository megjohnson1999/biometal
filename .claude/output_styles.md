# Custom Output Styles for biometal

## Benchmark Results Style

**Use Case**: Structured benchmark output for N=30 runs

**Format**:
```markdown
## Benchmark: {operation_name}

**Methodology**: N=30 runs, {platform}, {date}

| Implementation | Mean ± StdDev | Median | 95% CI | Speedup |
|----------------|---------------|--------|--------|---------|
| {baseline}     | {mean} ± {std} | {median} | [{ci_low}, {ci_high}] | 1.0× |
| {optimized}    | {mean} ± {std} | {median} | [{ci_low}, {ci_high}] | {speedup}× |

**Statistical Significance**: p={p_value}, effect size={cohen_d}

**ASBB Reference**: Rule {rule_num}, Entry {entry_num}

**Conclusion**: {interpretation}
```

**Trigger**: When reporting benchmark results

**Example**:
```
Output benchmark results for base_counting in benchmark_results_style
```

---

## Test Summary Style

**Use Case**: Concise test execution summary

**Format**:
```
✅ Test Suite: {passed}/{total} passed ({percent}%)

Library Tests: {lib_passed}/{lib_total}
Integration Tests: {integration_passed}/{integration_total}
Documentation Tests: {doc_passed}/{doc_total}
Python Tests: {python_passed}/{python_total}

Execution Time: {duration}

{failures_section if any}
```

**Trigger**: After cargo test or pytest execution

**Example**:
```
Output test results in test_summary_style
```

---

## Progress Report Style

**Use Case**: Weekly progress updates (used by progress-tracker agent)

**Format**:
```markdown
# Week {N} Progress Report ({date_range})

## Summary
{1-2 sentence overview}

**Schedule Status**: {on_track|ahead|behind}

## Completed Tasks ✅
- {task_1}
- {task_2}

## In Progress 🔄
- {task_1} (~{percent}%)

## Git Activity
- Commits: {count}
- Files Changed: {count}
- Lines: +{added} / -{removed}

## Testing
- Tests: {count} ({change})
- Benchmarks: {Y/N}

## Next Week Focus
- {task_1}
- {task_2}

## Blockers/Risks
- {blocker if any}

---
**Week {N} of 24** | **Phase {phase}** | **Overall: {percent}%**
```

**Trigger**: /weekly-review command

---

## Evidence Validation Style

**Use Case**: Structured output from evidence-validator agent

**Format**:
```markdown
## Evidence Validation: {optimization_name}

**Proposed Change**: {description}

**Performance Claim**: {expected_speedup}

**Evidence Check**:
- Rule 1 (NEON): {✅|⚠️|❌} {comment}
- Rule 2 (Block-based): {✅|⚠️|❌} {comment}
- Rule 3 (Parallel BGZF): {✅|⚠️|❌} {comment}
- Rule 4 (mmap): {✅|⚠️|❌} {comment}
- Rule 5 (Streaming): {✅|⚠️|❌} {comment}
- Rule 6 (Network): {✅|⚠️|❌} {comment}

**ASBB References**:
- {entry_citation}

**Risk Assessment**:
- Performance risk: {LOW|MEDIUM|HIGH}
- Correctness risk: {LOW|MEDIUM|HIGH}
- Portability risk: {LOW|MEDIUM|HIGH}

**Recommendation**: {GO|GO_WITH_CAUTION|NO_GO|VALIDATE_FIRST}

**Rationale**: {explanation}

{next_steps_if_validate_first}
```

**Trigger**: /check-evidence command

---

## Performance Comparison Style

**Use Case**: Comparing biometal vs samtools/pysam/etc.

**Format**:
```markdown
## Performance Comparison: {operation}

| Tool | Throughput | Memory | Platform | Notes |
|------|------------|--------|----------|-------|
| biometal | {value} | {value} | {platform} | {notes} |
| {competitor} | {value} | {value} | {platform} | {notes} |

**biometal Advantage**: {speedup}× throughput, {memory_reduction}× memory

**Use Cases**:
- ✅ {where_biometal_excels}
- ⚠️ {where_competitive}
- ❌ {where_competitor_better}
```

**Trigger**: When comparing with other tools

---

## Documentation Update Style

**Use Case**: Structured documentation sync recommendations

**Format**:
```markdown
📝 Documentation Sync Check

**Edited Files**: {list}

**Potential Impacts**:
- [ ] API documentation
- [ ] Performance benchmarks
- [ ] User Guide
- [ ] CLAUDE.md
- [ ] Examples

**Recommended Actions**:
1. {specific_file}: {specific_section} ({reason})
2. {specific_file}: {specific_section} ({reason})

**Validation Required**:
- [ ] Rerun benchmarks
- [ ] Validate examples
- [ ] Update changelog

Use `/update-docs` to proceed.
```

**Trigger**: AfterEdits hook

---

## Git Commit Message Style

**Use Case**: Structured commit messages following project conventions

**Format**:
```
{type}({scope}): {description}

{detailed_explanation}

{benchmarks_if_performance}

{breaking_changes_if_any}

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>
```

**Types**: feat, fix, perf, docs, test, refactor, chore

**Examples**:
- `perf(bam): Switch to cloudflare_zlib backend (+67% throughput)`
- `feat(index): Add BAI indexed queries (1.68-500× speedup)`
- `docs(guide): Update USER_GUIDE.md with BAI tutorial`

---

## Implementation Notes

These output styles should be used by:
- **Agents**: reference these styles in their responses
- **Slash Commands**: specify output style in command description
- **Hooks**: use appropriate style for hook output

**Priority**:
1. Benchmark Results Style (High - used frequently)
2. Evidence Validation Style (High - critical for decision making)
3. Progress Report Style (Medium - weekly usage)
4. Test Summary Style (Medium - regular usage)
5. Others (Low - nice-to-have)

---

**Note**: Claude Code output styles may evolve. Check official docs at https://code.claude.com/docs/en/output-styles for latest features.
