# Contributing to biometal

Thank you for your interest in contributing to biometal! This guide will help you get started.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Project Structure](#project-structure)
- [Development Workflow](#development-workflow)
- [Coding Standards](#coding-standards)
- [Testing](#testing)
- [Documentation](#documentation)
- [Submitting Changes](#submitting-changes)
- [Community](#community)

---

## Code of Conduct

biometal is committed to providing a welcoming and inclusive environment for all contributors. Please be respectful, considerate, and constructive in all interactions.

**Core principles:**
- Be respectful and inclusive
- Focus on constructive feedback
- Welcome newcomers and help them learn
- Assume good intentions

---

## Getting Started

### Ways to Contribute

1. **Report bugs** - Found a bug? [Open an issue](https://github.com/scotthandley/biometal/issues/new/choose)
2. **Suggest features** - Have an idea? [Create a feature request](https://github.com/scotthandley/biometal/issues/new/choose)
3. **Improve documentation** - Docs can always be better
4. **Fix bugs** - Check [open issues](https://github.com/scotthandley/biometal/issues?q=is%3Aissue+is%3Aopen+label%3Abug)
5. **Add features** - Implement something from the roadmap
6. **Optimize performance** - Help make biometal faster
7. **Write tutorials** - Share your use cases

### Good First Issues

Look for issues tagged with `good-first-issue` for beginner-friendly contributions:
[Good First Issues](https://github.com/scotthandley/biometal/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22)

---

## Development Setup

### Prerequisites

**For Rust development:**
- Rust 1.75.0 or later
- Cargo (comes with Rust)
- Git

**For Python development:**
- Python 3.8 or later
- maturin 1.0+ (`pip install maturin`)
- pytest (`pip install pytest`)

**Platform notes:**
- **macOS**: Xcode command line tools
- **Linux**: build-essential or equivalent
- **Windows**: MSVC build tools

### Clone the Repository

```bash
git clone https://github.com/scotthandley/biometal.git
cd biometal
```

### Build and Test

**Rust library:**
```bash
# Build
cargo build

# Run tests
cargo test

# Run benchmarks
cargo bench
```

**Python bindings:**
```bash
# Build and install locally (development mode)
maturin develop --release

# Run Python tests
pytest tests/python/

# Build wheels
maturin build --release
```

### Development Tools

**Recommended:**
```bash
# Install clippy (linter)
rustup component add clippy

# Install rustfmt (formatter)
rustup component add rustfmt

# Check code
cargo clippy --all-targets --all-features

# Format code
cargo fmt

# Generate documentation
cargo doc --no-deps --open
```

---

## Project Structure

```
biometal/
├── src/
│   ├── lib.rs              # Public API
│   ├── io/                 # Streaming parsers
│   │   ├── fastq.rs
│   │   ├── fasta.rs
│   │   ├── bam/            # BAM/SAM parser
│   │   └── ...
│   ├── operations/         # Analysis operations
│   ├── python/             # Python bindings (PyO3)
│   └── ...
├── tests/                  # Integration tests
├── benches/                # Criterion benchmarks
├── examples/               # Usage examples
├── docs/                   # Documentation
├── experiments/            # Research experiments
└── .github/                # GitHub templates

Key files:
- OPTIMIZATION_RULES.md     # Evidence-based optimization rules
- CLAUDE.md                 # Development guide (for Claude sessions)
- CHANGELOG.md              # Version history
- CONTRIBUTING.md           # This file
```

---

## Development Workflow

### 1. Create a Branch

```bash
git checkout -b feature/your-feature-name
# or
git checkout -b fix/bug-description
```

**Branch naming:**
- `feature/` - New features
- `fix/` - Bug fixes
- `docs/` - Documentation changes
- `perf/` - Performance improvements
- `refactor/` - Code refactoring

### 2. Make Changes

Follow the [Coding Standards](#coding-standards) below.

### 3. Test Your Changes

```bash
# Run all tests
cargo test

# Run specific test
cargo test test_name

# Run benchmarks (if performance-related)
cargo bench
```

### 4. Update Documentation

- Add/update doc comments for public APIs
- Update relevant documentation files
- Add examples if introducing new functionality

### 5. Commit Your Changes

**Commit message format:**
```
<type>(<scope>): <subject>

<body>

<footer>
```

**Types:**
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `perf`: Performance improvement
- `refactor`: Code refactoring
- `test`: Test additions/changes
- `chore`: Build/tooling changes

**Examples:**
```bash
git commit -m "feat(bam): Add CRAM format support"
git commit -m "fix(fastq): Handle empty sequence lines correctly"
git commit -m "docs: Update BAI index tutorial with generation instructions"
git commit -m "perf(neon): Optimize base counting for M3 chips"
```

### 6. Push and Create PR

```bash
git push origin feature/your-feature-name
```

Then create a pull request on GitHub.

---

## Coding Standards

### Rust Code Style

**Follow Rust conventions:**
- Use `rustfmt` for consistent formatting
- Follow [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/)
- Use meaningful names (avoid abbreviations)

**Error handling:**
```rust
// ✅ GOOD: Use Result types
pub fn parse_record(data: &[u8]) -> Result<Record> {
    validate(data)?;
    Ok(Record::from_bytes(data))
}

// ❌ BAD: Don't use unwrap() or expect() in library code
pub fn parse_record(data: &[u8]) -> Record {
    Record::from_bytes(data).unwrap()  // Panics on error!
}
```

**Memory management:**
```rust
// ✅ GOOD: Streaming, constant memory
pub fn process_stream(path: &str) -> Result<()> {
    for record in FastqStream::from_path(path)? {
        process(record?);
    }
    Ok(())
}

// ❌ BAD: Accumulates in memory
pub fn process_all(path: &str) -> Result<Vec<Record>> {
    let mut records = Vec::new();
    for record in FastqStream::from_path(path)? {
        records.push(record?);  // Grows without bound!
    }
    Ok(records)
}
```

**Platform-specific code:**
```rust
// ✅ GOOD: Provide fallback for all platforms
#[cfg(target_arch = "aarch64")]
pub fn operation(data: &[u8]) -> Result {
    operation_neon(data)
}

#[cfg(not(target_arch = "aarch64"))]
pub fn operation(data: &[u8]) -> Result {
    operation_scalar(data)
}
```

### Evidence-Based Optimization

biometal follows evidence-based optimization from [apple-silicon-bio-bench](https://github.com/scotthandley/apple-silicon-bio-bench).

**Before implementing optimizations:**
1. Check `OPTIMIZATION_RULES.md` for validated rules
2. Reference specific experiments/evidence
3. Don't guess parameters - use validated values

**Example:**
```rust
// ✅ GOOD: Evidence-based
// Rule 2: Block size from Entry 027 (1,440 measurements)
const BLOCK_SIZE: usize = 10_000;

// ❌ BAD: Arbitrary choice
const BLOCK_SIZE: usize = 8_192;  // No evidence!
```

### Documentation

**All public APIs must have doc comments:**

```rust
/// Parse a FASTQ record from a byte slice.
///
/// This function validates the FASTQ format and returns a structured
/// `FastqRecord` with quality scores decoded.
///
/// # Arguments
///
/// * `data` - Raw FASTQ record bytes (4 lines)
///
/// # Returns
///
/// * `Ok(FastqRecord)` - Successfully parsed record
/// * `Err(BiometalError)` - Invalid format or I/O error
///
/// # Example
///
/// ```
/// use biometal::io::fastq::parse_record;
///
/// let data = b"@read1\nACGT\n+\nIIII\n";
/// let record = parse_record(data)?;
/// assert_eq!(record.id, "read1");
/// # Ok::<(), biometal::BiometalError>(())
/// ```
///
/// # Performance
///
/// - Constant memory: ~200 bytes per record
/// - Throughput: ~500 Krecords/s (scalar), ~8 Mrecords/s (NEON)
pub fn parse_record(data: &[u8]) -> Result<FastqRecord> {
    // implementation
}
```

---

## Testing

### Test Requirements

**All contributions must include tests:**

1. **Unit tests** - Test individual functions
2. **Integration tests** - Test end-to-end workflows
3. **Property-based tests** - For complex logic (use proptest)
4. **Benchmarks** - For performance-critical code

### Writing Tests

**Unit tests:**
```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_valid_record() {
        let data = b"@read1\nACGT\n+\nIIII\n";
        let record = parse_record(data).unwrap();
        assert_eq!(record.id, "read1");
        assert_eq!(record.sequence, b"ACGT");
    }

    #[test]
    fn test_parse_invalid_format() {
        let data = b"invalid";
        assert!(parse_record(data).is_err());
    }
}
```

**Property-based tests:**
```rust
use proptest::prelude::*;

proptest! {
    #[test]
    fn test_roundtrip(seq in "[ACGT]{1,1000}") {
        let encoded = encode(&seq);
        let decoded = decode(&encoded).unwrap();
        prop_assert_eq!(seq, decoded);
    }
}
```

**Benchmarks:**
```rust
use criterion::{criterion_group, criterion_main, Criterion};

fn bench_parse(c: &mut Criterion) {
    let data = generate_fastq_data(100_000);
    c.bench_function("parse_fastq_100k", |b| {
        b.iter(|| parse_fastq(&data))
    });
}

criterion_group!(benches, bench_parse);
criterion_main!(benches);
```

### Running Tests

```bash
# All tests
cargo test

# Specific test
cargo test test_parse_valid_record

# Integration tests only
cargo test --test '*'

# With output
cargo test -- --nocapture

# Benchmarks
cargo bench
```

---

## Documentation

### Types of Documentation

1. **API docs** - Rust doc comments (/// or //!)
2. **User guide** - `docs/USER_GUIDE.md`
3. **Tutorials** - Jupyter notebooks in `notebooks/`
4. **Examples** - Runnable code in `examples/`
5. **README** - Project overview

### Documentation Standards

**Doc comments should include:**
- Clear description
- Arguments and return values
- Example code (runnable with doc tests)
- Performance characteristics (if relevant)
- Error conditions

**Example format:**
```rust
/// [One-line summary]
///
/// [Detailed description]
///
/// # Arguments
///
/// * `arg1` - Description
///
/// # Returns
///
/// * `Ok(T)` - Success case
/// * `Err(E)` - Error cases
///
/// # Example
///
/// ```
/// [runnable example]
/// ```
///
/// # Performance
///
/// [Memory usage, throughput, etc.]
```

### Building Documentation

```bash
# Generate and open docs
cargo doc --no-deps --open

# Test documentation examples
cargo test --doc
```

---

## Submitting Changes

### Pull Request Checklist

Before submitting a PR, ensure:

- [ ] Code compiles without warnings (`cargo build`)
- [ ] All tests pass (`cargo test`)
- [ ] Code is formatted (`cargo fmt`)
- [ ] Clippy passes (`cargo clippy`)
- [ ] Documentation is updated
- [ ] Examples are added (if new feature)
- [ ] Benchmarks are added (if performance-related)
- [ ] Commit messages follow convention
- [ ] PR description is clear and complete

### PR Description Template

```markdown
## Description
[Clear description of what this PR does]

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Performance improvement
- [ ] Documentation update
- [ ] Refactoring

## Related Issue
Fixes #[issue number]

## Testing
[How did you test these changes?]

## Performance Impact
[If applicable, benchmark results before/after]

## Checklist
- [ ] Tests pass
- [ ] Documentation updated
- [ ] No new warnings
```

### Review Process

1. **Automated checks** - CI runs tests, clippy, fmt
2. **Code review** - Maintainer reviews code
3. **Discussion** - Address feedback, make changes
4. **Approval** - Once approved, PR is merged

**Expected timeline:**
- Initial review: 1-3 days
- Back-and-forth: As needed
- Merge: After approval and passing CI

---

## Community

### Communication Channels

- **GitHub Issues** - Bug reports, feature requests
- **GitHub Discussions** - Questions, ideas, showcase
- **Pull Requests** - Code contributions

### Getting Help

**Have questions?**
1. Check the [User Guide](docs/USER_GUIDE.md)
2. Search [existing issues](https://github.com/scotthandley/biometal/issues)
3. Ask in [Discussions](https://github.com/scotthandley/biometal/discussions)
4. Open a new issue with the `question` label

### Recognition

Contributors are recognized in:
- `CHANGELOG.md` - For each release
- GitHub contributor graph
- Project README (for significant contributions)

---

## Development Resources

### Project-Specific

- **OPTIMIZATION_RULES.md** - Evidence-based optimization rules
- **CLAUDE.md** - Development guide (for Claude AI sessions)
- **experiments/** - Research experiments and prototypes
- **apple-silicon-bio-bench** - [Evidence base](https://github.com/scotthandley/apple-silicon-bio-bench)

### Rust Resources

- [Rust Book](https://doc.rust-lang.org/book/)
- [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/)
- [Rust Performance Book](https://nnethercote.github.io/perf-book/)

### Bioinformatics Resources

- [SAM/BAM Specification](https://samtools.github.io/hts-specs/)
- [FASTQ Format](https://en.wikipedia.org/wiki/FASTQ_format)
- [BGZF Format](https://samtools.github.io/hts-specs/SAMv1.pdf)

---

## Questions?

If you have questions about contributing, please:
1. Check this guide first
2. Search [existing discussions](https://github.com/scotthandley/biometal/discussions)
3. Open a new discussion or issue

Thank you for contributing to biometal! 🚀

---

**Last Updated**: November 10, 2025
