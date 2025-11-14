# MCP Integration Recommendations for biometal

## Recommended MCP Servers

### 1. Git Integration (High Priority)
**Use Case**: Enhanced git operations (blame, history, complex diffs)
**Server**: `@modelcontextprotocol/server-git`
**Benefits**:
- Deeper git history analysis for weekly reviews
- Better attribution tracking
- Enhanced diff visualization

**Setup**:
```bash
npm install -g @modelcontextprotocol/server-git
```

### 2. GitHub Integration (High Priority)
**Use Case**: Issue tracking, PR management, release automation
**Server**: `@modelcontextprotocol/server-github`
**Benefits**:
- Create issues for bugs discovered during benchmarking
- Automate PR creation for weekly milestones
- Track community feedback
- Release note generation

**Setup**:
```bash
npm install -g @modelcontextprotocol/server-github
# Requires GITHUB_TOKEN environment variable
```

### 3. Database for Benchmark Results (Medium Priority)
**Use Case**: Track benchmark history, detect regressions, generate performance reports
**Server**: `@modelcontextprotocol/server-sqlite`
**Benefits**:
- Store N=30 benchmark results in queryable format
- Track performance trends over time
- Generate regression reports
- Compare across versions/platforms

**Setup**:
```bash
npm install -g @modelcontextprotocol/server-sqlite
# Create benchmarks.db for historical benchmark data
```

**Schema**:
```sql
CREATE TABLE benchmarks (
    id INTEGER PRIMARY KEY,
    date TEXT,
    operation TEXT,
    implementation TEXT,
    mean REAL,
    stddev REAL,
    median REAL,
    min REAL,
    max REAL,
    platform TEXT,
    version TEXT
);
```

### 4. Slack/Discord Integration (Low Priority)
**Use Case**: Community communication, milestone announcements
**Server**: `@modelcontextprotocol/server-slack` or Discord equivalent
**Benefits**:
- Announce weekly progress to community
- Post benchmark results
- Coordinate with contributors

**Deferred**: Until community grows

### 5. Documentation Search (Medium Priority)
**Use Case**: Quick reference to 142K+ words of documentation
**Server**: Custom MCP server for local markdown search
**Benefits**:
- Semantic search across CLAUDE.md, USER_GUIDE.md, etc.
- Find relevant sections during development
- Ensure consistency across docs

**Implementation**:
Could use `@modelcontextprotocol/server-filesystem` with grep-like capabilities, or build custom server with semantic search.

## Integration Priority

**Phase 1 (Immediate)**:
1. Git integration - Enhance weekly reviews
2. GitHub integration - Issue tracking and releases

**Phase 2 (Next Month)**:
3. SQLite for benchmarks - Historical tracking
4. Documentation search - Maintain consistency

**Phase 3 (Future)**:
5. Community tools (Slack/Discord) - When community reaches critical mass

## Configuration Example

`.claude/mcp.json`:
```json
{
  "servers": {
    "git": {
      "command": "npx",
      "args": ["-y", "@modelcontextprotocol/server-git"],
      "cwd": "/Users/scotthandley/Code/biometal"
    },
    "github": {
      "command": "npx",
      "args": ["-y", "@modelcontextprotocol/server-github"],
      "env": {
        "GITHUB_TOKEN": "${GITHUB_TOKEN}"
      }
    },
    "benchmarks": {
      "command": "npx",
      "args": ["-y", "@modelcontextprotocol/server-sqlite", "benchmarks/benchmark_history.db"]
    }
  }
}
```

## Usage Examples

### Git History Analysis
```
"Use the git MCP server to analyze commits in the past 2 months that touched src/io/bam/"
```

### GitHub Issue Creation
```
"Create a GitHub issue for the NEON optimization regression detected in benchmark N=30 run"
```

### Benchmark Tracking
```
"Query the benchmarks database to compare base_counting performance across the last 5 versions"
```

## Custom MCP Server: Bioinformatics Data

**Future Consideration**: Build custom MCP server for:
- NCBI/SRA API integration (fetch datasets)
- UniProt/RefSeq reference data
- Genomic coordinate utilities
- BLAST/BLAT alignment services

**Deferred**: Until core functionality complete (Phase 3+)

---

**Recommendation**: Start with Git + GitHub integration (high ROI, low setup cost). Add SQLite benchmarking database once historical data accumulates (Week 4+).
