# Repository Size Investigation Report

**Date**: November 13, 2025
**Total Repository Size**: 9.0 GB
**Git History Size**: 669 MB
**Working Tree Size**: ~8.3 GB

---

## Executive Summary

**Root Cause**: Rust build artifacts (`target/` directories) were accidentally committed to git in research and experiment subdirectories.

**Impact**:
- 17,226 unnecessary files committed to git history
- 669 MB bloated .git directory
- 2.3+ GB of committed build artifacts
- Slow clone times for contributors

**Solution**: Remove build artifacts from git history, update .gitignore

---

## Detailed Breakdown

### 1. Working Tree Size (8.3 GB)

| Directory | Size | Status | Action |
|-----------|------|--------|--------|
| **target/** | 5.8 GB | ✅ In .gitignore | Keep (local builds) |
| **research/** | 1.7 GB | ❌ **COMMITTED** | Remove from git |
| **experiments/** | 756 MB | ❌ **COMMITTED** | Remove from git |
| **venv/** | 65 MB | ✅ In .gitignore | Keep (local venv) |
| **tests/** | 58 MB | ⚠️ Partially committed | Review test data |
| Other | ~200 MB | ✅ Legitimate code/docs | Keep |

---

### 2. Git Repository Analysis

**Total tracked files**: 17,589

**Breakdown by directory**:
```
16,542 files  research/          (94.0% of tracked files) ❌
   833 files  experiments/        (4.7% of tracked files) ❌
    46 files  src/               (0.3% of tracked files) ✅
    22 files  docs/              (0.1% of tracked files) ✅
    19 files  examples/          (0.1% of tracked files) ✅
    18 files  benches/           (0.1% of tracked files) ✅
    17 files  tests/             (0.1% of tracked files) ✅
```

**Problem**: 17,375 files (98.8%) are from research/experiments, mostly build artifacts!

---

### 3. Specific Culprits

#### A. Research CAF Target Directory ❌

**Path**: `research/caf-format/implementation/target/`
**Size**: 1.7 GB
**Tracked Files**: 16,476
**Status**: COMMITTED TO GIT (should be in .gitignore)

**Largest files in git history**:
```
27.5 MB  libreqwest-8534ade147fc52d7.rlib
27.5 MB  libreqwest-dec41876a454f0f4.rlib
20.7 MB  libcriterion-563bf8bff0cb5cf0.rlib
20.4 MB  bam_to_caf (binary)
18.1 MB  libtokio-27b60af946c5e6d9.rlib
16.9 MB  libproptest-391cdd06ed87615d.rlib
14.8 MB  test_compression_ratios (binary)
12.8 MB  dep-graph.bin (incremental)
```

**Why this happened**: The CAF research project (Nov 4-11, 2025) had its own Cargo project with a `target/` directory that was NOT in the .gitignore at that level.

#### B. Experiments Target Directories ❌

**Paths**:
- `experiments/archive/simd-minimizers-analysis/external/simd-minimizers/target/` (434 MB)
- `experiments/native-bam-implementation/profiling/target/` (115 MB)
- `experiments/archive/sra-decoder/target/` (8.7 MB)

**Tracked Files**: 750
**Status**: COMMITTED TO GIT

#### C. Test Data ⚠️

**Path**: `tests/data/large/large_1m.bam`
**Size**: 9.45 MB
**Status**: Committed to git (borderline - test data can be regenerated)

---

### 4. Current .gitignore Analysis

**What's correctly ignored** ✅:
```gitignore
/target/           # Main project target (works)
.venv/             # Python venv (works)
*.fq, *.fq.gz      # FASTQ files (works)
*.fa, *.fa.gz      # FASTA files (works)
```

**What's MISSING** ❌:
```gitignore
# Nested target/ directories (research, experiments)
**/target/         # Would catch ALL target/ directories

# BAM/BAI files (test data, large)
*.bam
*.bam.bai

# Criterion benchmark outputs
target/criterion/

# Virtual environments (broader)
venv/
.venv/
env/

# Research project artifacts
research/*/target/
experiments/*/target/
```

---

## Recommendations

### Priority 1: Remove Build Artifacts from Git History (CRITICAL)

**Impact**: Reduce git history from 669 MB to ~50 MB (92% reduction)

**Steps**:

1. **Create backup** (safety first):
   ```bash
   cd /Users/scotthandley/Code/biometal
   git clone --mirror . ../biometal-backup.git
   ```

2. **Install git-filter-repo** (safe alternative to filter-branch):
   ```bash
   brew install git-filter-repo  # macOS
   ```

3. **Remove target/ directories from git history**:
   ```bash
   git filter-repo --path-glob '**/target/**' --invert-paths --force
   ```

4. **Remove large BAM test files** (optional):
   ```bash
   git filter-repo --path-glob 'tests/data/large/*.bam' --invert-paths --force
   ```

5. **Force push** (⚠️ BREAKING CHANGE - coordinate with any collaborators):
   ```bash
   git push origin --force --all
   git push origin --force --tags
   ```

**Expected result**:
- Git history: 669 MB → ~50-100 MB
- Clone time: Minutes → Seconds
- Tracked files: 17,589 → ~200 (99% reduction)

**Risk**: Force push required. If others have cloned, they need to re-clone.

---

### Priority 2: Update .gitignore (IMMEDIATE)

**Add these patterns**:

```gitignore
# Rust build artifacts (catch ALL target directories)
**/target/
target/
Cargo.lock  # Already present

# Python virtual environments
venv/
.venv/
env/
.env/

# Test data (large files)
*.bam
*.bam.bai
tests/data/large/

# Benchmark results (criterion generates large files)
target/criterion/
bench_results/

# Jupyter notebook outputs
.ipynb_checkpoints/
**/.ipynb_checkpoints/

# Research/experiment build artifacts (belt-and-suspenders)
research/*/target/
research/*/Cargo.lock
experiments/*/target/
experiments/*/Cargo.lock

# OS/IDE files (already present, but complete list)
.DS_Store
.vscode/
.idea/
*.swp
*.swo
*~
```

**Apply immediately**:
```bash
# Update .gitignore (use Edit tool or manual edit)
git add .gitignore
git commit -m "chore: Expand .gitignore to catch nested target/ directories"
```

---

### Priority 3: Clean Working Tree (OPTIONAL)

**After** git history is cleaned, optionally remove local build artifacts to save disk space:

```bash
# Remove CAF research target/
rm -rf research/caf-format/implementation/target/

# Remove experiment targets
find experiments/ -name "target" -type d -exec rm -rf {} +

# Remove main target/ (will be rebuilt on next cargo build)
rm -rf target/

# Remove venv (will be rebuilt with pip install)
rm -rf venv/
```

**Savings**: ~8 GB → ~500 MB working tree

**Trade-off**: Next build will take longer (rebuilding from scratch)

---

### Priority 4: Consider Git LFS for Large Test Data (OPTIONAL)

If you need large test files (BAM, FASTQ):

1. **Install Git LFS**:
   ```bash
   brew install git-lfs
   git lfs install
   ```

2. **Track large file types**:
   ```bash
   git lfs track "*.bam"
   git lfs track "*.bam.bai"
   git lfs track "*.fq.gz"
   git lfs track "*.fa.gz"
   ```

3. **Update .gitattributes**:
   ```gitattributes
   *.bam filter=lfs diff=lfs merge=lfs -text
   *.bam.bai filter=lfs diff=lfs merge=lfs -text
   *.fq.gz filter=lfs diff=lfs merge=lfs -text
   *.fa.gz filter=lfs diff=lfs merge=lfs -text
   ```

**Benefits**:
- Large files stored externally (LFS server)
- Git repository stays small
- Files downloaded on-demand

**Costs**:
- GitHub LFS: 1 GB free, then $5/mo per 50 GB
- Alternative: Host test data externally (S3, Zenodo, Figshare)

---

## Alternative: External Test Data

**Instead of Git LFS**, use external hosting:

1. **Upload large test data to Zenodo** (free, permanent DOI):
   - Create dataset: "biometal test data"
   - Upload BAM/FASTQ files
   - Get permanent URL

2. **Add download script**:
   ```bash
   # scripts/download-test-data.sh
   #!/bin/bash
   wget https://zenodo.org/record/XXXXX/large_1m.bam -O tests/data/large/large_1m.bam
   ```

3. **Update documentation**:
   ```markdown
   ## Running Tests

   Large test data is hosted externally. Download with:
   ```bash
   ./scripts/download-test-data.sh
   ```
   ```

**Benefits**:
- Free (Zenodo, Figshare)
- Permanent DOI (citeable)
- No git bloat

---

## Implementation Plan

### Phase 1: Immediate (15 minutes)

```bash
# 1. Update .gitignore
cat >> .gitignore <<'EOF'

# Nested Rust build artifacts
**/target/

# Large test data
*.bam
*.bam.bai
tests/data/large/

# Additional Python environments
venv/
env/
EOF

# 2. Commit
git add .gitignore
git commit -m "chore: Prevent future target/ commits"
```

### Phase 2: Safe Cleanup (1 hour)

```bash
# 1. Create backup
cd /Users/scotthandley/Code/biometal
git clone --mirror . ../biometal-backup.git

# 2. Install git-filter-repo
brew install git-filter-repo

# 3. Remove build artifacts from history
git filter-repo --path-glob '**/target/**' --invert-paths --force

# 4. Verify size reduction
du -sh .git  # Should be ~50-100 MB (down from 669 MB)
git count-objects -vH

# 5. Force push (⚠️ COORDINATE WITH TEAM)
git push origin --force --all
git push origin --force --tags
```

### Phase 3: Optional Cleanup (30 minutes)

```bash
# Remove local build artifacts (if disk space needed)
rm -rf target/
rm -rf research/caf-format/implementation/target/
find experiments/ -name "target" -type d -exec rm -rf {} +
rm -rf venv/

# Reclaim space
du -sh .  # Should be ~500 MB (down from 9.0 GB)
```

---

## Expected Outcomes

| Metric | Before | After | Improvement |
|--------|--------|-------|-------------|
| **Total size** | 9.0 GB | ~500 MB | **95% reduction** |
| **Git history** | 669 MB | ~50-100 MB | **85-92% reduction** |
| **Tracked files** | 17,589 | ~200 | **99% reduction** |
| **Clone time** | ~5-10 minutes | ~10-30 seconds | **20-30× faster** |
| **Contributor friction** | HIGH (large clone) | LOW (small clone) | **Much better** |

---

## Verification Checklist

After cleanup:

- [ ] `.gitignore` updated with `**/target/` pattern
- [ ] Git history cleaned (`.git/` < 100 MB)
- [ ] Total repository < 1 GB
- [ ] Test clone works: `git clone <repo> /tmp/test-clone`
- [ ] Builds work: `cd /tmp/test-clone && cargo build`
- [ ] Tests pass: `cargo test --all-features`
- [ ] CI updated (if using GitHub Actions)

---

## Risks and Mitigation

### Risk 1: Force Push Breaks Collaborators

**Probability**: LOW (solo dev currently)
**Impact**: MEDIUM (collaborators need to re-clone)
**Mitigation**:
1. Check for collaborators: `git remote show origin`
2. Coordinate before force push
3. Document migration: "Re-clone required after Nov 13"

### Risk 2: Lost Important Files

**Probability**: VERY LOW (only removing build artifacts)
**Impact**: LOW (build artifacts regenerable)
**Mitigation**:
1. Create mirror backup: `git clone --mirror`
2. Verify backup before force push
3. Keep backup for 30 days

### Risk 3: CI/CD Breakage

**Probability**: MEDIUM (if CI caches .git directory)
**Impact**: MEDIUM (CI fails until cache cleared)
**Mitigation**:
1. Clear GitHub Actions cache after force push
2. Monitor first CI run after cleanup
3. Update CI if needed

---

## Post-Cleanup Monitoring

**After cleanup, track**:

1. **Repository size**:
   ```bash
   du -sh .git  # Should stay < 100 MB
   ```

2. **New large files**:
   ```bash
   git ls-files | xargs du -h | sort -hr | head -20
   ```

3. **Tracked file count**:
   ```bash
   git ls-files | wc -l  # Should be ~200, not 17,000
   ```

**Set monthly reminder**: Review `.gitignore` and repository size

---

## Conclusion

**Root Cause**: Nested `target/` directories not ignored
**Impact**: 9.0 GB repository (669 MB git history)
**Solution**: Update `.gitignore` + clean git history
**Effort**: 1-2 hours
**Result**: 95% size reduction, 20× faster clones

**Recommendation**: Execute Phase 1 (update .gitignore) IMMEDIATELY, then Phase 2 (cleanup) when ready for force push.

---

**Investigation Date**: November 13, 2025
**Investigator**: git-github-manager agent
**Status**: Ready for cleanup
