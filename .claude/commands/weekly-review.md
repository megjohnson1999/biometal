---
description: Generate weekly progress report and update PROJECT_TODOS.md
---

Generate comprehensive weekly progress report and update 24-week plan (PROJECT_TODOS.md).

Steps:
1. Analyze git activity for past week:
   - Commits (git log --since="1 week ago" --oneline)
   - Files changed (git diff --stat HEAD~7..HEAD)
   - Lines changed (additions/deletions)
2. Review PROJECT_TODOS.md:
   - Which tasks were completed?
   - Which tasks are in progress?
   - Are we on schedule?
3. Check test/benchmark status:
   - Test count (cargo test --list)
   - Recent benchmark runs (benchmarks/results/)
4. Documentation updates:
   - New/modified docs (git diff --name-only HEAD~7..HEAD '*.md')
5. Generate report:
   ```markdown
   # Week N Progress Report (<date range>)

   ## Completed Tasks
   - <task from PROJECT_TODOS.md> ✅
   - <task from PROJECT_TODOS.md> ✅

   ## In Progress
   - <task from PROJECT_TODOS.md> 🔄 (<percent complete>)

   ## Git Activity
   - Commits: <count>
   - Files changed: <count>
   - Lines added: <count>, removed: <count>

   ## Testing
   - Tests: <count> (<change from last week>)
   - Benchmarks: <run this week? Y/N>

   ## Documentation
   - <new/updated docs>

   ## Next Week Focus
   - <tasks from PROJECT_TODOS.md>

   ## Blockers/Risks
   - <any issues encountered>

   ## Schedule Status
   - <on track / ahead / behind> (current week vs planned)
   ```
6. Update PROJECT_TODOS.md:
   - Mark completed tasks ✅
   - Update in-progress percentages
   - Adjust timeline if needed
7. Archive report: `planning/weekly_reports/week_N_YYYY_MM_DD.md`

Optional: Generate comparison with STRATEGIC_PIVOT_PLAN.md to validate alignment.

Expected time: 10-15 minutes

Generates: Weekly report + updated PROJECT_TODOS.md
