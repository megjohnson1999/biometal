---
description: Run statistically rigorous benchmark (N=30) and generate formatted report
---

Run a statistically rigorous benchmark for the specified operation following ASBB methodology (N=30 samples).

Steps:
1. Identify the benchmark in benchmarks/ directory
2. Execute with cargo bench (Criterion handles N=30)
3. Extract results from target/criterion/
4. Calculate statistical measures (mean, median, stddev, 95% CI)
5. Format results according to benchmark report template
6. Cross-reference with OPTIMIZATION_RULES.md
7. Provide speedup analysis vs baseline

Use the benchmark-specialist agent for execution.

Arguments expected: <benchmark_name> (e.g., "base_counting", "bam_parsing")

Generate a markdown report ready for inclusion in documentation.
