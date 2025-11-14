# Norovirus Receptor Discovery Pipeline - Project Summary

## Overview

Complete computational pipeline for discovering the unknown cellular receptor for human norovirus GII.17 using AlphaFold3-based structure prediction and systematic screening of human intestinal membrane proteins.

## What Was Built

### 1. Core Pipeline Modules

**Candidate Library Builder** (`src/candidates/library_builder.py`)
- Queries UniProt REST API for human intestinal membrane proteins
- Filters by tissue specificity, topology, and localization
- Extracts metadata: expression, localization, protein class
- Exports: JSON, FASTA, statistics

**Structure Preparation** (`src/structures/structure_prep.py`)
- Extracts ectodomain regions from full-length proteins
- Topology prediction (UniProt annotations + hydrophobicity heuristic)
- VP1 P-domain sequence retrieval
- AlphaFold3 JSON input generation
- Optional: DeepTMHMM/TMHMM integration

**AlphaFold3 Runner** (`src/prediction/alphafold_runner.py`)
- Multi-backend support: local, AlphaFold3 Server, Google Colab
- Batch processing with configurable parallelization
- Generates 5 models per complex
- Extracts ipTM, pTM, pLDDT scores

**Interface Analyzer** (`src/analysis/interface_analyzer.py`)
- Identifies interface residues (distance-based)
- Counts interactions: H-bonds, salt bridges, hydrophobic, aromatic
- Calculates interface area and geometry scores
- Extracts pLDDT at binding interface
- Uses BioPython for PDB parsing

**Scoring & Ranking** (`src/analysis/scoring.py`)
- Integrates 7 scoring components with configurable weights
- Calculates composite scores: structural confidence, biological relevance
- Ranks all candidates by overall score
- Assigns confidence tiers: high/medium/low

**Visualization** (`src/validation/visualization.py`)
- Score distribution plots
- Confidence tier pie charts
- ipTM vs expression scatter plots
- Score comparison bar charts

### 2. Pipeline Orchestrator

**Main Pipeline** (`pipeline.py`)
- CLI interface for full pipeline or individual stages
- Configuration management (JSON)
- Progress logging and error handling
- Automatic report generation

### 3. Configuration & Documentation

**Configuration Files**
- `config/default.json`: Default parameters
- Configurable: backends, weights, filters, thresholds

**Documentation**
- `README.md`: Complete user guide (6000+ words)
  - Installation instructions
  - Stage-by-stage explanations
  - Examples and troubleshooting
  - Validation recommendations
- `PROJECT_SUMMARY.md`: This file

**Dependencies**
- `requirements.txt`: All Python dependencies
- Core: BioPython, NumPy, pandas, matplotlib, seaborn
- Optional: PyMOL, freesasa

### 4. Analysis Tools

**Jupyter Notebook** (`notebooks/01_analyze_results.ipynb`)
- Load and explore results
- Visualize score distributions
- Identify top candidates
- Export validation lists
- Summary statistics

### 5. Automation

**Quick Start Script** (`quickstart.sh`)
- Automated setup: venv, dependencies, directories
- Python version checking
- Usage instructions

## Pipeline Architecture

```
Input: GII.17 VP1 P-domain + Human proteome
  ↓
[1] Candidate Library Builder
  → 500-2000 intestinal membrane proteins
  ↓
[2] Structure Preparation
  → Ectodomain extraction + AF3 inputs
  ↓
[3] AlphaFold3 Predictions
  → VP1-receptor complex structures (5 models each)
  ↓
[4] Interface Analysis
  → Binding interactions + interface metrics
  ↓
[5] Scoring & Ranking
  → Integrated scores → Top candidates
  ↓
Output: Ranked list with validation recommendations
```

## Key Features

### Scalability
- Handles 500-2000 candidates
- Parallelizable across multiple GPUs
- Batch processing with progress tracking

### Flexibility
- Configurable scoring weights
- Multiple AlphaFold3 backends
- Customizable filters and thresholds
- Modular design (run individual stages)

### Robustness
- Error handling throughout
- Multiple fallbacks (topology prediction)
- Comprehensive logging
- Input validation

### Analysis
- Multi-metric scoring (7 components)
- Statistical analysis tools
- Visualization suite
- Jupyter notebook integration

## Computational Requirements

### Minimum
- Python 3.9+
- 8 GB RAM
- 50 GB disk space

### Recommended
- Python 3.9+
- 32 GB RAM
- 1x GPU (NVIDIA A100 or similar)
- 200 GB disk space
- AlphaFold3 installed locally

### For 500 Candidates
- CPU-only: ~10-25 days
- 1x A100 GPU: ~25-40 hours
- 10x A100 GPUs: ~3-4 hours

## Expected Outputs

### Data Files
- `data/processed/candidates/candidate_library.json` (~5-10 MB)
- `data/processed/structures/af3_inputs/*.json` (500-2000 files)
- `results/predictions/*/ranked_*.pdb` (~50 GB total)
- `results/interface_analysis.json` (~10-20 MB)
- `results/ranked_candidates.json` (~2-5 MB)

### Reports
- `results/PIPELINE_REPORT.md`: Full pipeline summary
- Console output: Top 10 candidates
- Plots: Distribution, correlation, confidence tiers

### Typical Results
- High confidence (ipTM > 0.8): 5-20 candidates
- Medium confidence (ipTM 0.6-0.8): 50-100 candidates
- Top candidate should have:
  - ipTM > 0.8
  - Interface area > 800 Å²
  - Multiple interaction types
  - High intestinal expression

## Validation Path

### Computational (Immediate)
1. Review top 10 candidates
2. Check interface details (H-bonds, contacts)
3. Verify biological plausibility
4. Run MD simulations on top 3-5 (optional)

### Experimental (Next Steps)
1. **Biochemical**: VLP binding assays (ELISA, SPR)
2. **Cell-based**: CRISPR knockout/overexpression
3. **Structural**: Co-crystallization
4. **In vivo**: Animal models

## Success Criteria

**Pipeline Success**:
- ✓ Process ≥500 candidates
- ✓ ≥5 candidates with ipTM > 0.8
- ✓ Clear ranking with justification
- ✓ Reproducible results

**Biological Success** (Experimental):
- True receptor binds VP1 with Kd < 1 μM
- Knockout abolishes infection
- Overexpression confers susceptibility
- Crystal structure validates prediction

## Future Enhancements

### Short-term
- DeepTMHMM/Phobius integration
- Conservation analysis across GII strains
- Population genetics (gnomAD integration)
- PyMOL visualization scripts

### Medium-term
- Multi-strain screening (GII.4, GII.17, GI.1)
- RoseTTAFold-AA validation
- DiffDock orthogonal docking
- Automated MD simulations

### Long-term
- Integration with experimental data (RNA-seq, proteomics)
- Machine learning for score optimization
- Web interface for interactive analysis
- Real-time prediction monitoring

## File Structure

```
norovirus-receptor-discovery/
├── pipeline.py                 # Main orchestrator
├── quickstart.sh              # Setup script
├── requirements.txt           # Dependencies
├── README.md                  # User guide
├── LICENSE                    # MIT license
├── PROJECT_SUMMARY.md         # This file
│
├── config/
│   └── default.json          # Default configuration
│
├── src/
│   ├── candidates/
│   │   └── library_builder.py
│   ├── structures/
│   │   └── structure_prep.py
│   ├── prediction/
│   │   └── alphafold_runner.py
│   ├── analysis/
│   │   ├── interface_analyzer.py
│   │   └── scoring.py
│   └── validation/
│       └── visualization.py
│
├── notebooks/
│   └── 01_analyze_results.ipynb
│
├── data/
│   ├── raw/
│   ├── processed/
│   │   ├── candidates/
│   │   └── structures/
│   └── results/
│
└── results/
    ├── predictions/
    ├── plots/
    └── PIPELINE_REPORT.md
```

## Module Sizes

- `library_builder.py`: ~450 lines
- `structure_prep.py`: ~550 lines
- `alphafold_runner.py`: ~400 lines
- `interface_analyzer.py`: ~550 lines
- `scoring.py`: ~500 lines
- `visualization.py`: ~350 lines
- `pipeline.py`: ~350 lines
- **Total**: ~3,150 lines of production Python code

## Testing Recommendations

### Unit Tests
- Topology prediction accuracy
- Ectodomain extraction correctness
- Scoring function validation
- Interface analysis metrics

### Integration Tests
- Full pipeline on 10 test proteins
- Config loading and validation
- Error handling (missing files, bad inputs)

### Performance Tests
- Memory usage monitoring
- Prediction timing benchmarks
- Scalability testing (100, 500, 1000 candidates)

## Citation

When using this pipeline, please cite:
- AlphaFold3: Abramson et al. (2024) Nature
- UniProt: UniProt Consortium (2023) Nucleic Acids Research
- BioPython: Cock et al. (2009) Bioinformatics

## Acknowledgments

This pipeline integrates multiple state-of-the-art tools and databases:
- AlphaFold3 (DeepMind/Google)
- UniProt (EMBL-EBI/SIB/PIR)
- BioPython (Open Bioinformatics Foundation)

## Contact & Support

- Documentation: See README.md
- Issues: GitHub issue tracker
- Questions: Open a discussion

---

**Version**: 1.0.0
**Date**: November 2025
**Status**: Production-ready, requires AlphaFold3 access
**License**: MIT
