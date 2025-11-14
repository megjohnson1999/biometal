---
description: Systematically update documentation after code changes
---

Systematically update documentation to reflect recent code changes.

Steps:
1. Analyze recent git changes (git diff, git log)
2. Identify affected documentation:
   - API docs (if public API changed)
   - Performance docs (if benchmarks/optimizations changed)
   - User Guide (if usage patterns changed)
   - CLAUDE.md (if architecture changed)
   - Examples (if example code affected)
3. For each affected doc:
   - Show current content
   - Propose updates based on code changes
   - Validate consistency with code
4. Update CHANGELOG.md if user-facing changes
5. Suggest version bump if appropriate (SemVer)

Use doc-specialist agent for comprehensive updates.

Optional arguments: <scope> (e.g., "api", "performance", "user-guide", "all")

Example: /update-docs api
Example: /update-docs all

Generates: Updated documentation files + CHANGELOG.md entry
