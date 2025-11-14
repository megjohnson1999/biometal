---
description: Run all benchmarks and generate comprehensive performance report
---

Execute all project benchmarks with N=30 statistical rigor and generate comprehensive report.

Steps:
1. Discover all benchmarks in benchmarks/ directory
2. Run each benchmark via cargo bench
3. Collect and analyze all results
4. Generate comprehensive performance matrix:
   - Operation-by-operation breakdown
   - Speedup vs scalar baseline
   - Memory usage analysis
   - Cross-platform comparison (if applicable)
5. Identify performance regressions vs previous version
6. Update CHANGELOG.md if significant changes detected
7. Archive results with timestamp

Use the benchmark-specialist agent for execution.

Expected time: 15-30 minutes depending on benchmark count.

Generates: benchmarks/results/<timestamp>_full_report.md
