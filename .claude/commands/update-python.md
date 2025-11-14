---
description: Update Python bindings to reflect Rust API changes
---

Update Python bindings (PyO3) to reflect recent Rust API changes.

Steps:
1. Analyze recent Rust API changes (public functions/structs)
2. Identify missing or outdated Python bindings
3. Generate/update PyO3 wrappers:
   - Function signatures
   - Error handling (Result → PyResult)
   - Type conversions
   - Docstrings
4. Generate Python tests for new/updated bindings
5. Build and test:
   ```bash
   cd python_bindings
   maturin develop
   pytest tests/
   ```
6. Update docs/PYTHON.md with new APIs
7. Validate performance (benchmark if critical)

Use python-bindings agent for implementation.

Optional arguments: <scope> (e.g., "fastq", "bam", "all")

Example: /update-python bam
Example: /update-python all

Generates: Updated python_bindings/src/lib.rs + tests + documentation
