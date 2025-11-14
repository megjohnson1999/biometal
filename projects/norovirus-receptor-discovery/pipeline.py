#!/usr/bin/env python3
"""
Norovirus Receptor Discovery Pipeline

Main orchestrator for the computational receptor discovery workflow.
"""

import logging
import sys
import argparse
from pathlib import Path
from typing import Optional
import json

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from candidates.library_builder import CandidateLibraryBuilder
from structures.structure_prep import StructurePreparationPipeline
from prediction.alphafold_runner import AlphaFold3Runner
from analysis.interface_analyzer import analyze_batch
from analysis.scoring import score_all_candidates, ScoringWeights

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class NororvirusReceptorPipeline:
    """Main pipeline orchestrator"""

    def __init__(self, config_file: Optional[Path] = None):
        """
        Args:
            config_file: Optional path to configuration JSON
        """
        self.config = self._load_config(config_file)
        self.data_dir = Path(self.config.get('data_dir', 'data'))
        self.results_dir = Path(self.config.get('results_dir', 'results'))

    def _load_config(self, config_file: Optional[Path]) -> dict:
        """Load configuration from file or use defaults"""
        default_config = {
            'data_dir': 'data',
            'results_dir': 'results',
            'af3_backend': 'local',
            'max_candidates': 500,
            'distance_cutoff': 5.0,
            'scoring_weights': {}
        }

        if config_file and config_file.exists():
            with open(config_file) as f:
                user_config = json.load(f)
            default_config.update(user_config)

        return default_config

    def run_full_pipeline(self):
        """Execute complete pipeline"""
        logger.info("=" * 60)
        logger.info("NOROVIRUS RECEPTOR DISCOVERY PIPELINE")
        logger.info("=" * 60)

        # Stage 1: Build candidate library
        logger.info("\n[Stage 1/5] Building candidate library...")
        candidates = self.build_candidate_library()

        # Stage 2: Prepare structures
        logger.info("\n[Stage 2/5] Preparing structures for AlphaFold3...")
        prepared = self.prepare_structures()

        # Stage 3: Run predictions
        logger.info("\n[Stage 3/5] Running AlphaFold3 predictions...")
        predictions = self.run_predictions()

        # Stage 4: Analyze interfaces
        logger.info("\n[Stage 4/5] Analyzing binding interfaces...")
        analyses = self.analyze_interfaces()

        # Stage 5: Score and rank
        logger.info("\n[Stage 5/5] Scoring and ranking candidates...")
        ranked = self.score_and_rank()

        # Generate summary report
        self.generate_report(candidates, prepared, predictions, analyses, ranked)

        logger.info("\n" + "=" * 60)
        logger.info("PIPELINE COMPLETE")
        logger.info("=" * 60)

    def build_candidate_library(self):
        """Stage 1: Build candidate library"""
        output_dir = self.data_dir / "processed" / "candidates"

        builder = CandidateLibraryBuilder(output_dir)
        candidates = builder.build_library()

        # Optional: Filter to top N candidates
        max_candidates = self.config.get('max_candidates', 500)
        if len(candidates) > max_candidates:
            logger.info(f"Filtering to top {max_candidates} candidates")
            candidates = self._filter_candidates(candidates, max_candidates)

        logger.info(f"✓ Built library with {len(candidates)} candidates")
        return candidates

    def prepare_structures(self):
        """Stage 2: Prepare structures"""
        library_file = self.data_dir / "processed" / "candidates" / "candidate_library.json"
        output_dir = self.data_dir / "processed" / "structures"

        pipeline = StructurePreparationPipeline(output_dir)
        prepared = pipeline.prepare_from_library(library_file)

        logger.info(f"✓ Prepared {len(prepared)} AlphaFold3 inputs")
        return prepared

    def run_predictions(self):
        """Stage 3: Run AlphaFold3 predictions"""
        input_dir = self.data_dir / "processed" / "structures" / "af3_inputs"
        output_dir = self.results_dir / "predictions"
        backend = self.config.get('af3_backend', 'local')

        input_files = list(input_dir.glob("*.json"))
        logger.info(f"Found {len(input_files)} complexes to predict")

        runner = AlphaFold3Runner(backend=backend, output_dir=output_dir)
        predictions = runner.run_batch(input_files, max_parallel=1)

        successful = sum(1 for p in predictions.values() if p)
        logger.info(f"✓ Completed {successful}/{len(input_files)} predictions")
        return predictions

    def analyze_interfaces(self):
        """Stage 4: Analyze interfaces"""
        prediction_dir = self.results_dir / "predictions"
        output_file = self.results_dir / "interface_analysis.json"
        distance_cutoff = self.config.get('distance_cutoff', 5.0)

        analyses = analyze_batch(prediction_dir, output_file, distance_cutoff)

        logger.info(f"✓ Analyzed {len(analyses)} interfaces")
        return analyses

    def score_and_rank(self):
        """Stage 5: Score and rank candidates"""
        interface_file = self.results_dir / "interface_analysis.json"
        library_file = self.data_dir / "processed" / "candidates" / "candidate_library.json"
        output_file = self.results_dir / "ranked_candidates.json"

        # Load custom weights if provided
        weights_config = self.config.get('scoring_weights', {})
        weights = ScoringWeights(**weights_config) if weights_config else None

        ranked = score_all_candidates(
            interface_file,
            library_file,
            output_file,
            weights
        )

        logger.info(f"✓ Ranked {len(ranked)} candidates")
        return ranked

    def _filter_candidates(self, candidates, max_count: int):
        """Filter candidates by priority criteria"""
        # Prioritize by protein class
        priority_classes = [
            'receptor', 'integrin', 'lectin', 'cd_antigen',
            'tight_junction', 'adhesion'
        ]

        prioritized = []
        for pclass in priority_classes:
            matching = [c for c in candidates if c.protein_class == pclass]
            prioritized.extend(matching)

        # Add remaining candidates
        remaining = [c for c in candidates if c not in prioritized]
        prioritized.extend(remaining)

        return prioritized[:max_count]

    def generate_report(self, candidates, prepared, predictions, analyses, ranked):
        """Generate final summary report"""
        report_file = self.results_dir / "PIPELINE_REPORT.md"

        report = f"""# Norovirus Receptor Discovery Pipeline Report

## Pipeline Summary

### Stage 1: Candidate Library
- Total candidates identified: {len(candidates)}
- Data source: UniProt (reviewed, human, intestinal membrane proteins)

### Stage 2: Structure Preparation
- Successfully prepared: {len(prepared)}
- VP1 P-domain: GII.17 strain (residues 225-530)

### Stage 3: AlphaFold3 Predictions
- Total complexes: {len(predictions)}
- Successful predictions: {sum(1 for p in predictions.values() if p)}
- Backend: {self.config.get('af3_backend', 'local')}

### Stage 4: Interface Analysis
- Interfaces analyzed: {len(analyses)}
- Distance cutoff: {self.config.get('distance_cutoff', 5.0)} Å

### Stage 5: Scoring and Ranking
- Total scored: {len(ranked)}
- High confidence: {sum(1 for r in ranked if r.confidence_tier == 'high')}
- Medium confidence: {sum(1 for r in ranked if r.confidence_tier == 'medium')}
- Low confidence: {sum(1 for r in ranked if r.confidence_tier == 'low')}

## Top 20 Candidates

| Rank | Gene | Protein ID | Overall Score | ipTM | Confidence |
|------|------|------------|---------------|------|------------|
"""
        for i, candidate in enumerate(ranked[:20], 1):
            report += f"| {i} | {candidate.gene_name} | {candidate.protein_id} | {candidate.overall_score:.3f} | {candidate.ipTM_score:.3f} | {candidate.confidence_tier} |\n"

        report += f"""
## High Confidence Candidates (ipTM > 0.8)

"""
        high_conf = [r for r in ranked if r.ipTM_score > 0.8]
        if high_conf:
            for candidate in high_conf[:10]:
                report += f"""### {candidate.gene_name} ({candidate.protein_id})
- **Overall Score**: {candidate.overall_score:.3f}
- **ipTM**: {candidate.ipTM_score:.3f}
- **Interface pLDDT**: {candidate.interface_pLDDT_score * 100:.1f}
- **Structural Confidence**: {candidate.structural_confidence:.3f}
- **Biological Relevance**: {candidate.biological_relevance:.3f}

"""
        else:
            report += "No candidates with ipTM > 0.8\n"

        report += f"""
## Recommended Next Steps

### Experimental Validation (Top 3-5 Candidates)
1. **Ectodomain Expression**: Express as Fc-fusion proteins
2. **VLP Binding Assays**: ELISA, SPR, or BLI
3. **Cell-Based Validation**:
   - CRISPR knockout in permissive cells
   - Overexpression in non-permissive cells
4. **Structure Determination**: Co-crystal of VP1-receptor complex

### Computational Follow-up
1. **Molecular Dynamics**: 100-500 ns MD on top 5 candidates
2. **Alternative Tools**: Validate with RoseTTAFold, DiffDock
3. **Conservation Analysis**: Check binding site conservation across GII strains
4. **Population Genetics**: Check for resistance-associated variants

## Files Generated

- Candidate library: `data/processed/candidates/candidate_library.json`
- AlphaFold3 inputs: `data/processed/structures/af3_inputs/`
- Predictions: `results/predictions/`
- Interface analyses: `results/interface_analysis.json`
- Ranked candidates: `results/ranked_candidates.json`
- This report: `results/PIPELINE_REPORT.md`

---
*Report generated by Norovirus Receptor Discovery Pipeline*
"""

        report_file.parent.mkdir(parents=True, exist_ok=True)
        with open(report_file, 'w') as f:
            f.write(report)

        logger.info(f"✓ Report saved to {report_file}")


def main():
    """Main CLI entry point"""
    parser = argparse.ArgumentParser(
        description="Norovirus Receptor Discovery Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run full pipeline
  python pipeline.py --full

  # Run specific stages
  python pipeline.py --stage candidates
  python pipeline.py --stage structures
  python pipeline.py --stage predictions
  python pipeline.py --stage analysis
  python pipeline.py --stage scoring

  # With custom config
  python pipeline.py --full --config config/custom.json
        """
    )

    parser.add_argument(
        '--full',
        action='store_true',
        help='Run complete pipeline'
    )
    parser.add_argument(
        '--stage',
        choices=['candidates', 'structures', 'predictions', 'analysis', 'scoring'],
        help='Run specific pipeline stage'
    )
    parser.add_argument(
        '--config',
        type=Path,
        help='Path to configuration JSON file'
    )

    args = parser.parse_args()

    if not args.full and not args.stage:
        parser.print_help()
        sys.exit(1)

    # Initialize pipeline
    pipeline = NororvirusReceptorPipeline(args.config)

    try:
        if args.full:
            pipeline.run_full_pipeline()
        elif args.stage == 'candidates':
            pipeline.build_candidate_library()
        elif args.stage == 'structures':
            pipeline.prepare_structures()
        elif args.stage == 'predictions':
            pipeline.run_predictions()
        elif args.stage == 'analysis':
            pipeline.analyze_interfaces()
        elif args.stage == 'scoring':
            pipeline.score_and_rank()

        logger.info("\n✓ Pipeline completed successfully")

    except Exception as e:
        logger.error(f"\n✗ Pipeline failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
