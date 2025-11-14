# Norovirus Receptor Discovery Pipeline

Computational pipeline for discovering the unknown cellular receptor for human norovirus (GII.17) using AlphaFold3 and structure-based screening.

## Overview

While histo-blood group antigens (HBGAs) are known **attachment factors** for norovirus, the true **cellular receptor** that mediates viral entry remains unknown. This pipeline uses state-of-the-art structural prediction (AlphaFold3) to systematically screen human intestinal membrane proteins and identify high-confidence receptor candidates.

### Approach

1. **Build Candidate Library**: Query UniProt for human intestinal epithelial membrane proteins with extracellular domains
2. **Prepare Structures**: Extract ectodomains and prepare AlphaFold3 inputs
3. **Predict Binding**: Use AlphaFold3 to predict VP1 P-domain binding to each candidate
4. **Analyze Interfaces**: Detailed analysis of binding interfaces (interactions, geometry, quality)
5. **Score and Rank**: Integrate structural and biological data to rank candidates

## Installation

### Requirements

- Python 3.9+
- AlphaFold3 (see installation options below)
- GPU recommended (but not required)

### Setup

```bash
# Clone repository
git clone <repository-url>
cd norovirus-receptor-discovery

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### AlphaFold3 Installation Options

**Option 1: Local Installation** (Recommended for large-scale screening)
- Follow [AlphaFold3 installation guide](https://github.com/deepmind/alphafold)
- Requires GPU with CUDA support
- ~3-5 minutes per prediction on A100 GPU

**Option 2: AlphaFold3 Server** (Free tier: 20 predictions/day)
- Create account at [AlphaFold Server](https://alphafoldserver.com)
- Set API key in config
- Best for initial pilot studies

**Option 3: Google Colab** (Free GPU, manual upload)
- Use AlphaFold3 Colab notebook
- Manual file upload/download required
- Good for small batches (10-50 predictions)

## Quick Start

### Run Complete Pipeline

```bash
# Run with default configuration
python pipeline.py --full

# Run with custom config
python pipeline.py --full --config config/custom.json
```

### Run Individual Stages

```bash
# Stage 1: Build candidate library
python pipeline.py --stage candidates

# Stage 2: Prepare structures
python pipeline.py --stage structures

# Stage 3: Run predictions (requires AlphaFold3)
python pipeline.py --stage predictions

# Stage 4: Analyze interfaces
python pipeline.py --stage analysis

# Stage 5: Score and rank
python pipeline.py --stage scoring
```

## Pipeline Stages in Detail

### Stage 1: Candidate Library Builder

**What it does:**
- Queries UniProt for human intestinal membrane proteins
- Filters for proteins with extracellular domains
- Extracts metadata (expression, localization, function)

**Query criteria:**
- Organism: Human (taxonomy 9606)
- Location: Membrane (cell membrane)
- Tissue: Intestine/colon/duodenum/jejunum/ileum
- Quality: Reviewed (Swiss-Prot)
- Topology: Has transmembrane or GPI anchor

**Output:**
- `data/processed/candidates/candidate_library.json`: Full metadata
- `data/processed/candidates/candidate_library.fasta`: Sequences
- `data/processed/candidates/library_statistics.json`: Summary stats

**Example:**
```bash
python -m src.candidates.library_builder \
  --output-dir data/processed/candidates
```

**Expected:** ~500-2000 candidates

### Stage 2: Structure Preparation

**What it does:**
- Extracts ectodomain regions from full-length proteins
- Predicts topology (if not known from UniProt)
- Retrieves GII.17 VP1 P-domain sequence
- Prepares AlphaFold3 JSON inputs

**Topology prediction:**
- Priority: UniProt annotations (most reliable)
- Fallback: Hydrophobicity-based heuristic
- Optional: DeepTMHMM/TMHMM (if installed)

**Output:**
- `data/processed/structures/af3_inputs/*.json`: AF3 input files
- `data/processed/structures/ectodomains.json`: Ectodomain info
- `data/processed/structures/preparation_summary.json`: Summary

**Example:**
```bash
python -m src.structures.structure_prep \
  --library data/processed/candidates/candidate_library.json \
  --output-dir data/processed/structures
```

### Stage 3: AlphaFold3 Predictions

**What it does:**
- Runs AlphaFold3 on VP1 + ectodomain complexes
- Generates 5 models per complex
- Extracts ipTM, pTM, and pLDDT scores

**Computational requirements:**
- GPU: 1× A100 (~3-5 min/prediction)
- CPU: ~30-60 min/prediction (much slower)
- Disk: ~100 MB per complex

**For 500 candidates:**
- A100 GPU: ~25-40 hours
- CPU: ~10-25 days (not recommended)

**Output:**
- `results/predictions/<complex_id>/ranked_*.pdb`: Predicted structures
- `results/predictions/<complex_id>/scores.json`: Metrics
- `results/predictions/prediction_summary.json`: Overall summary

**Example:**
```bash
# Local installation
python -m src.prediction.alphafold_runner \
  --input-dir data/processed/structures/af3_inputs \
  --output-dir results/predictions \
  --backend local

# Using AlphaFold3 Server
python -m src.prediction.alphafold_runner \
  --input-dir data/processed/structures/af3_inputs \
  --output-dir results/predictions \
  --backend server
```

### Stage 4: Interface Analysis

**What it does:**
- Identifies interface residues (distance < 5 Å)
- Counts interactions (H-bonds, salt bridges, hydrophobic, aromatic)
- Calculates interface area and geometry scores
- Extracts pLDDT at interface

**Metrics computed:**
- Interface area (Å²)
- Interface pLDDT (confidence at binding site)
- Hydrogen bonds
- Salt bridges
- Hydrophobic contacts
- Aromatic interactions
- Interface gap score (tightness)
- Shape complementarity

**Output:**
- `results/interface_analysis.json`: Complete analyses

**Example:**
```bash
python -m src.analysis.interface_analyzer \
  --prediction-dir results/predictions \
  --output results/interface_analysis.json \
  --distance-cutoff 5.0
```

### Stage 5: Scoring and Ranking

**What it does:**
- Integrates structural + biological metrics
- Calculates composite scores
- Ranks all candidates
- Identifies high-confidence hits

**Scoring components:**
| Component | Weight | Description |
|-----------|--------|-------------|
| ipTM | 30% | Interface prediction confidence |
| Interface quality | 20% | H-bonds, salt bridges, contacts |
| Interface area | 10% | Binding surface size |
| Consistency | 10% | Agreement across models |
| Expression | 10% | Intestinal expression level |
| Localization | 10% | Apical membrane targeting |
| Conservation | 10% | Binding site conservation |

**Confidence tiers:**
- **High**: ipTM > 0.8, pLDDT > 80, overall > 0.75
- **Medium**: ipTM > 0.6, pLDDT > 70, overall > 0.6
- **Low**: Below medium thresholds

**Output:**
- `results/ranked_candidates.json`: Ranked list with scores
- Console: Top 10 candidates

**Example:**
```bash
python -m src.analysis.scoring \
  --interface-analyses results/interface_analysis.json \
  --candidate-library data/processed/candidates/candidate_library.json \
  --output results/ranked_candidates.json
```

## Configuration

Create custom configuration in `config/custom.json`:

```json
{
  "data_dir": "data",
  "results_dir": "results",

  "af3_backend": "local",
  "max_candidates": 500,
  "distance_cutoff": 5.0,

  "scoring_weights": {
    "ipTM_weight": 0.30,
    "interface_quality_weight": 0.20,
    "expression_weight": 0.15
  },

  "filters": {
    "min_ectodomain_length": 50,
    "min_interface_area": 500,
    "min_ipTM": 0.5
  }
}
```

## Interpreting Results

### Key Metrics

**ipTM (Interface Predicted TM-score):**
- Range: 0-1 (higher is better)
- > 0.8: High confidence binding
- 0.6-0.8: Medium confidence
- < 0.6: Low confidence

**Interface pLDDT:**
- Range: 0-100 (higher is better)
- > 80: High confidence structure
- 70-80: Medium confidence
- < 70: Low confidence

**Overall Score:**
- Range: 0-1 (weighted composite)
- > 0.75: Strong candidate for validation
- 0.6-0.75: Worth investigating
- < 0.6: Lower priority

### Expected Results

**Typical screening outcomes:**
- High confidence candidates (ipTM > 0.8): 5-20
- Medium confidence (ipTM 0.6-0.8): 50-100
- Low confidence (ipTM < 0.6): Remaining

**Top candidates should show:**
- High ipTM (> 0.8)
- Large interface area (> 800 Å²)
- Multiple interaction types
- High intestinal expression
- Apical membrane localization

## Validation Recommendations

### Top 3-5 Candidates

**1. Biochemical Validation**
- Express ectodomain as Fc-fusion protein
- VLP binding assays (ELISA, SPR, BLI)
- Measure Kd (expect nM-μM range)

**2. Cell-Based Validation**
- CRISPR knockout in permissive cells → loss of infection
- Overexpression in non-permissive cells → gain of susceptibility
- Antibody blocking experiments

**3. Structural Validation**
- Co-crystallize VP1 + receptor ectodomain
- Compare with AF3 prediction (RMSD < 2 Å = excellent)

**4. In Vivo Validation**
- Tissue expression pattern (matches norovirus tropism?)
- Humanized mouse models
- Clinical correlation (receptor polymorphisms vs susceptibility)

## Advanced Usage

### Parallel Predictions (10 GPUs)

```bash
# Split inputs into 10 batches
python scripts/split_inputs.py \
  --input-dir data/processed/structures/af3_inputs \
  --num-splits 10

# Run on each GPU (separate terminals or job scheduler)
for i in {0..9}; do
  CUDA_VISIBLE_DEVICES=$i python pipeline.py --stage predictions \
    --config config/gpu_$i.json &
done
```

### Custom Scoring Weights

Emphasize structural confidence over biological data:

```json
{
  "scoring_weights": {
    "ipTM_weight": 0.40,
    "interface_quality_weight": 0.30,
    "interface_area_weight": 0.15,
    "consistency_weight": 0.10,
    "expression_weight": 0.03,
    "localization_weight": 0.02,
    "conservation_weight": 0.00
  }
}
```

### Alternative VP1 Strains

Modify `config/vp1_strains.json` to screen against multiple strains:

```json
{
  "strains": [
    {"name": "GII.17", "p_domain_start": 225, "p_domain_end": 530},
    {"name": "GII.4", "p_domain_start": 225, "p_domain_end": 530},
    {"name": "GI.1", "p_domain_start": 225, "p_domain_end": 530}
  ]
}
```

## Troubleshooting

### Issue: UniProt query returns too few candidates

**Solution:**
- Broaden tissue specificity filter
- Include "provisional" (TrEMBL) proteins (lower quality)
- Check UniProt API status

### Issue: AlphaFold3 predictions fail

**Common causes:**
- Out of GPU memory → Reduce batch size or use CPU
- Sequence too long → Check ectodomain extraction
- Missing dependencies → Check AF3 installation

### Issue: Low ipTM scores across all candidates

**Possible reasons:**
- Ectodomain boundaries incorrect → Review topology predictions
- VP1 sequence wrong → Verify strain and region
- True receptor not in candidate library → Expand search criteria

## Citation

If you use this pipeline in your research, please cite:

- AlphaFold3: [citation pending]
- UniProt: The UniProt Consortium (2023) Nucleic Acids Res.
- BioPython: Cock et al. (2009) Bioinformatics

## Contributing

Contributions welcome! Areas for improvement:
- Additional topology prediction tools (DeepTMHMM, Phobius)
- Conservation analysis across norovirus strains
- Integration with experimental data (RNA-seq, proteomics)
- Visualization tools (PyMOL scripts, interactive dashboards)

## License

MIT License - see LICENSE file for details

## Contact

For questions or issues, please open a GitHub issue or contact [your contact info].

---

**Disclaimer**: This is a computational prediction pipeline. All high-confidence candidates require experimental validation. Computational predictions, even with high confidence scores, do not guarantee biological function.
