# Known Issues

## Python Bindings

### ✅ RESOLVED: PyO3 Registration Bug (SAM Reader, GFA Writer, BAM Writer)

**Status**: ✅ RESOLVED
**Resolution Date**: November 16, 2025
**Original Date**: November 14-15, 2025

#### Resolution Summary
All three previously failing Python bindings (PySamReader, PyGfaWriter, PyBamWriter) are now fully functional after fixing a type error in `src/python/gtf.rs`.

#### Root Cause
A compilation error in `src/python/gtf.rs` line 126 was preventing successful PyO3 linking:
```rust
// BEFORE (incorrect):
fn transcript_id(&self) -> String {
    self.inner.transcript_id().to_string()  // Error: Option<&str> has no to_string()
}

// AFTER (correct):
fn transcript_id(&self) -> Option<String> {
    self.inner.transcript_id().map(|s| s.to_string())
}
```

#### Verification
All three classes tested and verified functional (November 16, 2025):

✅ **BamWriter**: Complete write/read roundtrip verified
```python
reader = biometal.BamReader.from_path("input.bam")
writer = biometal.BamWriter.create("output.bam", reader.header)
for record in reader:
    writer.write_record(record)
writer.finish()
```

✅ **SamReader**: Can read SAM files, access headers, iterate records
```python
reader = biometal.SamReader.from_path("input.sam")
print(reader.header.reference_count)
for record in reader:
    print(record.name)
```

✅ **GfaWriter**: Can create and write GFA files
```python
writer = biometal.GfaWriter.create("output.gfa.gz")
writer.finish()
```

#### Lessons Learned
1. PyO3 compilation errors can cause linking failures that prevent class registration
2. The error was not in the failing classes themselves, but in an unrelated file (gtf.rs)
3. Using `maturin build` instead of `cargo build --features python` properly handles PyO3 linking
4. Always verify .so file is updated after build (see CLAUDE.md Python Bindings section)

---

## Historical Documentation (RESOLVED ISSUES)

### SAM Reader Python Binding Registration Issue (RESOLVED)

**Status**: ✅ RESOLVED
**Severity**: Low (Rust implementation fully functional)
**Date**: November 14, 2025

#### Description
The SAM reader implementation is fully functional in Rust with all tests passing (4/4), but the Python bindings fail to register despite correct PyO3 annotations.

#### Details
- **Rust Implementation**: ✅ Complete (`src/io/bam/sam_reader.rs`)
  - All 4 unit tests passing
  - Streaming reader with constant memory
  - Full SAM format support
- **Python Bindings**: ❌ Not appearing in module
  - `PySamReader` struct properly annotated with `#[pyclass(name = "SamReader")]`
  - Registered in `src/python/mod.rs` with `m.add_class::<PySamReader>()?`
  - Compiles without errors or warnings
  - **Symbols not present in compiled .so file** (confirmed with `nm`)

#### Investigation Steps Taken
1. ✅ Verified `#[pyclass]` annotation
2. ✅ Verified module registration
3. ✅ Clean rebuild (cargo clean && maturin develop)
4. ✅ Checked for compilation errors (none)
5. ✅ Verified other similar classes (SamWriter works fine)
6. ✅ Attempted CompressedReader type change
7. ✅ Checked symbol table with `nm` (symbols missing)

#### Hypothesis
Likely a PyO3 macro expansion issue with generic type `SamReader<CompressedReader>` where:
- Trait bounds on generic types may not be compatible with PyO3's requirements
- Macro silently fails without producing compiler diagnostics
- Similar pattern works for `BamReader<CompressedReader>` but not for `SamReader<CompressedReader>`

#### Workaround
Use Rust API directly. SAM reading functionality is fully available via:
```rust
use biometal::io::bam::SamReader;

let sam = SamReader::from_path("input.sam")?;
for record in sam.records() {
    // process records
}
```

#### Next Steps
1. Investigate PyO3 trait bounds requirements
2. Consider using `SamReader<File>` instead of `SamReader<CompressedReader>` for Python bindings
3. Open issue with PyO3 project if needed
4. Consult PyO3 community/Discord

#### References
- PyO3 documentation: https://pyo3.rs/
- Implementation: `src/python/bam.rs:1851-1926`
- Rust tests: `src/io/bam/sam_reader.rs:399-462`

---

### GFA Writer Python Binding Registration Issue

**Status**: Under Investigation
**Severity**: Low (Rust implementation fully functional)
**Date**: November 14, 2025

#### Description
The GFA writer implementation is fully functional in Rust with all tests passing (2/2), but the Python bindings fail to register despite correct PyO3 annotations. **This exhibits the exact same pattern as the SAM reader issue above.**

#### Details
- **Rust Implementation**: ✅ Complete (`src/formats/gfa.rs:657-918`)
  - All 2 unit tests passing
  - Full GFA writing support (header, segment, link, path)
  - Automatic compression detection (.gfa vs .gfa.gz)
- **Python Bindings**: ❌ Not appearing in module
  - `PyGfaWriter` struct properly annotated with `#[pyclass(name = "GfaWriter", unsendable)]`
  - Registered in `src/python/mod.rs` with `m.add_class::<PyGfaWriter>()?`
  - Compiles without errors or warnings
  - **Symbols not present in compiled .so file** (confirmed with `nm`)
  - Full rebuild (cargo clean) does not resolve the issue

#### Pattern Analysis
This is the **second occurrence** of this PyO3 registration mystery:
1. ❌ `PySamReader` (SAM reading)
2. ❌ `PyGfaWriter` (GFA writing)

Both follow identical pattern:
- Correct `#[pyclass]` annotations
- Successful compilation
- Module registration in mod.rs
- Symbols missing from .so file
- Other classes in same file work fine (e.g., `PyGfaSegment`, `PyGfaLink`, `PyGfaPath` all work)

#### Hypothesis
Common factor between SAM reader and GFA writer that differs from working classes:
- Both are **I/O classes** (reader/writer)
- Both may have **generic constraints** or **trait bounds** that PyO3 macros can't handle
- Working classes are simpler data structures with straightforward PyO3 mappings

#### Workaround
Use Rust API directly. GFA writing functionality is fully available via:
```rust
use biometal::formats::gfa::GfaWriter;

let mut writer = GfaWriter::create("output.gfa.gz")?;
writer.write_header(tags)?;
writer.write_segment(&segment)?;
writer.finish()?;
```

#### Next Steps
1. Compare PyGfaWriter with working writers (PyVcfWriter, PyBed3Writer)
2. Investigate if `Option<GfaWriter>` pattern causes issues
3. Test with simplified struct (remove `path: String` field)
4. Coordinate investigation with SAM reader issue (likely same root cause)

#### References
- PyO3 documentation: https://pyo3.rs/
- Implementation: `src/python/gfa.rs:374-575`
- Rust tests: `src/formats/gfa.rs:1034-1151`

---

### BAM Writer Python Binding Registration Issue

**Status**: Under Investigation
**Severity**: Medium (Rust implementation fully functional, but blocks Python workflows)
**Date**: November 15, 2025

#### Description
The BAM writer implementation is fully functional in Rust with all tests passing (14/14: 6 unit + 8 integration), but the Python bindings fail to register despite correct PyO3 annotations. **This is the third occurrence of this exact pattern.**

#### Details
- **Rust Implementation**: ✅ Complete (`src/io/bam/writer.rs`)
  - All 14 tests passing (6 unit + 8 integration)
  - Full BAM writing support (header, records, BGZF compression)
  - Round-trip verification with real-world 100K+ record files
  - Production-ready quality
- **Python Bindings**: ❌ Not appearing in module
  - `PyBamWriter` struct properly annotated with `#[pyclass(name = "BamWriter", unsendable)]`
  - Registered in `src/python/mod.rs` with `m.add_class::<PyBamWriter>()?`
  - Compiles without errors or warnings
  - **Symbols not present in compiled .so file** (confirmed with `nm`)
  - Full rebuild (maturin develop --release) does not resolve the issue

#### Pattern Analysis
This is the **third occurrence** of this PyO3 registration mystery:
1. ❌ `PySamReader` (SAM reading)
2. ❌ `PyGfaWriter` (GFA writing)
3. ❌ `PyBamWriter` (BAM writing) **NEW**

All three follow identical pattern:
- Correct `#[pyclass]` annotations (with `unsendable` flag)
- Successful compilation with no errors
- Module registration in mod.rs
- Symbols missing from .so file
- Other classes in same file work fine

#### Common Factors
Analyzing what's different between working and failing classes:
1. **All are I/O writer/reader classes** with complex internal state
2. **All use `Option<T>` wrapping**:
   - `PySamReader`: `Option<SamReader<CompressedReader>>`
   - `PyGfaWriter`: `Option<GfaWriter>`
   - `PyBamWriter`: `Option<BamWriter>`
3. **All contain enum types internally**:
   - `BamWriter` contains `CompressedWriter` (enum)
   - `GfaWriter` likely similar
   - `SamReader<CompressedReader>` uses `CompressedReader` (enum)
4. **Working classes** (PyBamReader, PySamWriter) also use `Option<T>` with enums, so this is not the distinguishing factor

#### Hypothesis
Based on three occurrences, the common pattern suggests:
- Issue may be related to **specific enum types** or **trait bounds** on generic parameters
- PyO3 macro expansion fails silently without compiler diagnostics
- May be a PyO3 version-specific bug (using PyO3 0.27)
- Could be related to `unsendable` flag combined with specific type patterns

#### Workaround
Use Rust API directly. BAM writing functionality is fully available via:
```rust
use biometal::io::bam::BamWriter;

let mut writer = BamWriter::create("output.bam", header)?;
writer.write_record(&record)?;
writer.finish()?;
```

#### Impact
- **Rust users**: No impact, full functionality available
- **Python users**: Cannot write BAM files from Python
  - Blocks filtering workflows (read BAM → filter → write BAM)
  - Blocks format conversion (SAM → BAM, etc.)
  - Workaround: Use Rust API or call biometal CLI tools

#### Next Steps
1. Test with different PyO3 versions (0.26, 0.28) to isolate version-specific behavior
2. Compare working `PyBamReader` and `PySamWriter` implementations line-by-line
3. Try removing `unsendable` flag to see if that's a contributing factor
4. Investigate if `CompressedWriter` enum needs special PyO3 treatment
5. Coordinate investigation across all three issues (likely same root cause)
6. Consider reporting to PyO3 project with minimal reproduction case

#### References
- PyO3 documentation: https://pyo3.rs/
- Implementation: `src/python/bam.rs:1829-1987`
- Rust tests: `src/io/bam/writer.rs:500-590`, `tests/bam_writer_integration.rs:1-422`
- Python test: `tests/python/test_bam_writer.py` (demonstrates issue)
