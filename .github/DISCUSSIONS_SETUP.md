# GitHub Discussions Setup Guide

This guide outlines how to set up and manage GitHub Discussions for the biometal community.

## Enabling Discussions

1. Go to **Settings** → **General** → **Features**
2. Check **Discussions**
3. Click **Set up discussions**

## Recommended Categories

### 1. 📣 Announcements
- **Purpose**: Official project announcements (releases, roadmap updates)
- **Format**: Announcement
- **Who can post**: Maintainers only
- **Description**: "Official biometal announcements, releases, and roadmap updates"

### 2. 💡 Ideas
- **Purpose**: Feature requests and enhancement proposals
- **Format**: Idea (with voting)
- **Who can post**: Anyone
- **Description**: "Suggest new features or enhancements. Vote on ideas you'd like to see!"

### 3. 🙏 Q&A
- **Purpose**: Questions about using biometal
- **Format**: Q&A (with accepted answers)
- **Who can post**: Anyone
- **Description**: "Ask questions about installation, usage, performance, or troubleshooting"

### 4. 🚀 Show and Tell
- **Purpose**: Community showcases (papers, pipelines, use cases)
- **Format**: Open-ended discussion
- **Who can post**: Anyone
- **Description**: "Share your biometal projects, papers, pipelines, or success stories!"

### 5. 📊 Performance
- **Purpose**: Performance discussions, benchmarks, optimization tips
- **Format**: Open-ended discussion
- **Who can post**: Anyone
- **Description**: "Discuss performance, share benchmarks, and optimization strategies"

### 6. 🧬 Bioinformatics Workflows
- **Purpose**: Share and discuss bioinformatics workflows using biometal
- **Format**: Open-ended discussion
- **Who can post**: Anyone
- **Description**: "Share complete workflows, pipelines, and analysis patterns"

### 7. 🛠️ Development
- **Purpose**: Contributing, development setup, architecture discussions
- **Format**: Open-ended discussion
- **Who can post**: Anyone
- **Description**: "Discuss contributions, development setup, and architecture decisions"

### 8. 🌍 General
- **Purpose**: Anything that doesn't fit other categories
- **Format**: Open-ended discussion
- **Who can post**: Anyone
- **Description**: "General discussion about biometal and ARM-native bioinformatics"

## Welcome Post Template

Create a pinned welcome post in the **General** category:

---

**Title**: Welcome to biometal Discussions! 👋

**Body**:
```markdown
# Welcome to biometal Discussions! 👋

We're excited to have you here! This is the place for:

- ❓ **Questions** about installation, usage, or troubleshooting
- 💡 **Ideas** for new features or enhancements
- 🚀 **Showcasing** your projects, papers, or pipelines
- 📊 **Discussing** performance and optimization
- 🧬 **Sharing** bioinformatics workflows
- 🛠️ **Contributing** to development

## Quick Links

- 📘 [User Guide](https://github.com/scotthandley/biometal/blob/main/docs/USER_GUIDE.md) - Comprehensive documentation
- ⚡ [Performance Guide](https://github.com/scotthandley/biometal/blob/main/docs/PERFORMANCE_OPTIMIZATION_GUIDE.md) - Optimization tips
- 📊 [Benchmark Comparison](https://github.com/scotthandley/biometal/blob/main/benchmarks/comparison/BENCHMARK_COMPARISON.md) - vs samtools/pysam
- 🔧 [Contributing Guide](https://github.com/scotthandley/biometal/blob/main/CONTRIBUTING.md) - How to contribute
- 📖 [API Docs](https://docs.rs/biometal) - Complete API reference

## Getting Started

**Installation:**
```bash
# Python
pip install biometal-rs

# Rust
cargo add biometal
```

**Quick Example:**
```python
import biometal

# BAI-indexed query (500× faster on large files)
index = biometal.BaiIndex.from_path("alignments.bam.bai")
for record in biometal.BamReader.query_region(
    "alignments.bam", index, "chr1", 1_000_000, 2_000_000
):
    print(f"Read {record.qname} at {record.pos}")
```

## Community Guidelines

- **Be respectful** and inclusive
- **Search first** before asking questions
- **Share your knowledge** - help others learn
- **Provide context** - include code, errors, platform info
- **Give credit** - cite papers, acknowledge contributors

## What's New in v1.6.0

- ✨ **BAI index support** - 1.68-500× speedup for targeted analysis
- 📊 **Comprehensive benchmarks** - vs samtools/pysam (validated)
- 💾 **10-200× lower memory** - Constant 5 MB usage
- 🚀 **4-25× ARM NEON speedup** - Apple Silicon & Graviton optimized
- 📚 **40,000+ words of docs** - User guide, tutorials, optimization guide

[Read the full announcement →](blog/v1.6.0_announcement.md)

## Join the Conversation

Choose a category above to get started, or introduce yourself below! We'd love to hear:

- 🌍 Where are you from?
- 🔬 What do you work on?
- 💻 What platform do you use (Mac ARM, Graviton, x86_64)?
- 🧬 What bioinformatics workflows are you running?

Looking forward to building this community together! 🚀

— biometal team
```

---

## Moderation Guidelines

### Response Priority

1. **High priority** (respond within 24 hours):
   - Questions tagged with `urgent`
   - Bug reports
   - Breaking issues

2. **Medium priority** (respond within 48-72 hours):
   - General questions
   - Feature discussions
   - Performance questions

3. **Low priority** (respond when available):
   - General discussions
   - Showcases (congratulate!)

### Response Templates

**For common questions:**

> Thanks for your question! This is covered in our [User Guide](link-to-section). Here's the relevant section:
>
> [quote relevant section]
>
> Let me know if that helps!

**For feature requests:**

> Thanks for the suggestion! This aligns with our [roadmap/doesn't align with our current focus].
>
> [If aligned]: Would you be interested in contributing to this? Check our [Contributing Guide](CONTRIBUTING.md).
>
> [If not aligned]: We're focusing on X right now, but we'll keep this in mind for future phases.

**For showcases:**

> Amazing work! 🎉 Thanks for sharing. Would you be interested in writing this up as a tutorial for the community?

### Weekly Maintenance

**Every Monday:**
- Review unanswered questions
- Close stale discussions (30+ days inactive)
- Pin important discussions
- Update announcements category

**Every release:**
- Post announcement in Announcements category
- Update pinned welcome post
- Highlight community contributions

## Automation (Optional)

Consider setting up GitHub Actions for:

1. **Auto-label** - Automatically label discussions by category
2. **Stale discussions** - Mark inactive discussions
3. **Welcome bot** - Welcome first-time contributors

---

**Last Updated**: November 10, 2025
