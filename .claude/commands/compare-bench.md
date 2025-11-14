---
description: Compare benchmark results across implementations/versions with statistical analysis
---

Compare benchmark results between two implementations or versions with full statistical analysis.

Steps:
1. Run benchmark for baseline implementation (N=30)
2. Run benchmark for comparison implementation (N=30)
3. Perform statistical tests:
   - Two-sample t-test for significance
   - Effect size (Cohen's d)
   - Confidence intervals for difference
4. Generate comparison table with speedup ratios
5. Visualize results if appropriate
6. Provide interpretation in context of OPTIMIZATION_RULES.md

Use the benchmark-specialist agent for execution.

Arguments expected: <baseline> <comparison> (e.g., "scalar neon", "v1.6.0 v1.7.0")

Output includes statistical significance and practical significance assessment.
