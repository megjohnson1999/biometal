---
description: Validate optimization proposal against ASBB evidence (OPTIMIZATION_RULES.md)
---

Validate a proposed optimization against experimental evidence from OPTIMIZATION_RULES.md.

Use the evidence-validator agent to:
1. Parse the optimization proposal
2. Check against Rules 1-6
3. Identify relevant ASBB entries
4. Assess risks (performance, correctness, portability)
5. Provide GO / GO_WITH_CAUTION / NO_GO / VALIDATE_FIRST recommendation
6. Suggest next steps if validation needed

Arguments expected: <optimization_description>

Example: /check-evidence "Use NEON SIMD for CIGAR parsing"
Example: /check-evidence "Parallelize BAM record parsing"
Example: /check-evidence "Load records into Vec for faster access"

Output: Structured evidence validation report with recommendation and rationale.
