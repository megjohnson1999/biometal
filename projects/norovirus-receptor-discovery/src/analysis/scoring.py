"""
Scoring and Ranking Module

Integrates multiple metrics to score and rank receptor candidates.
"""

import logging
import json
import numpy as np
from dataclasses import dataclass, asdict
from typing import List, Dict, Optional
from pathlib import Path

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class CandidateScore:
    """Complete scoring for a receptor candidate"""
    complex_id: str
    protein_id: str
    gene_name: str

    # Structural scores (from AlphaFold3)
    ipTM_score: float  # 0-1, higher = better
    pTM_score: float
    interface_pLDDT_score: float  # 0-100
    interface_area_score: float  # Normalized 0-1

    # Interface quality scores
    interaction_score: float  # Based on H-bonds, salt bridges, etc.
    geometry_score: float  # Gap score + shape complementarity
    consistency_score: float  # Across multiple models

    # Biological plausibility scores
    expression_score: float  # Intestinal expression
    localization_score: float  # Apical membrane
    conservation_score: float  # Binding site conservation

    # Composite scores
    structural_confidence: float  # Weighted structural metrics
    biological_relevance: float  # Weighted biological metrics
    overall_score: float  # Final integrated score

    # Ranking
    rank: Optional[int] = None
    confidence_tier: Optional[str] = None  # 'high', 'medium', 'low'

    # Supporting data
    num_models: int = 5
    best_model: int = 1
    notes: str = ""


class ScoringWeights:
    """Configurable weights for scoring components"""

    def __init__(self, **kwargs):
        # Structural weights
        self.ipTM_weight = kwargs.get('ipTM_weight', 0.30)
        self.interface_quality_weight = kwargs.get('interface_quality_weight', 0.20)
        self.interface_area_weight = kwargs.get('interface_area_weight', 0.10)
        self.consistency_weight = kwargs.get('consistency_weight', 0.10)

        # Biological weights
        self.expression_weight = kwargs.get('expression_weight', 0.10)
        self.localization_weight = kwargs.get('localization_weight', 0.10)
        self.conservation_weight = kwargs.get('conservation_weight', 0.10)

        # Validate weights sum to 1.0
        total = (
            self.ipTM_weight +
            self.interface_quality_weight +
            self.interface_area_weight +
            self.consistency_weight +
            self.expression_weight +
            self.localization_weight +
            self.conservation_weight
        )

        if abs(total - 1.0) > 0.01:
            logger.warning(f"Weights sum to {total}, not 1.0. Normalizing...")
            norm = 1.0 / total
            self.ipTM_weight *= norm
            self.interface_quality_weight *= norm
            self.interface_area_weight *= norm
            self.consistency_weight *= norm
            self.expression_weight *= norm
            self.localization_weight *= norm
            self.conservation_weight *= norm


class CandidateScorer:
    """Score and rank receptor candidates"""

    def __init__(self, weights: Optional[ScoringWeights] = None):
        self.weights = weights or ScoringWeights()

    def score_candidate(
        self,
        complex_id: str,
        protein_id: str,
        gene_name: str,
        interface_analyses: List[Dict],
        candidate_data: Dict,
        conservation_data: Optional[Dict] = None
    ) -> CandidateScore:
        """
        Score a single candidate based on all available data

        Args:
            complex_id: Complex identifier
            protein_id: UniProt ID
            gene_name: Gene name
            interface_analyses: List of interface analysis results (one per model)
            candidate_data: Candidate protein metadata
            conservation_data: Optional conservation analysis results

        Returns:
            CandidateScore object
        """
        # Extract best model metrics
        best_model = self._select_best_model(interface_analyses)

        # Structural scores
        ipTM_score = best_model.get('ipTM', 0.0)
        pTM_score = best_model.get('pTM', 0.0)
        interface_pLDDT = best_model.get('interface_pLDDT', 0.0)
        interface_area = best_model.get('interface_area', 0.0)

        # Normalize scores
        interface_pLDDT_score = interface_pLDDT / 100.0  # Convert to 0-1
        interface_area_score = self._normalize_area(interface_area)

        # Calculate interaction score
        interaction_score = self._calculate_interaction_score(best_model)

        # Calculate geometry score
        geometry_score = self._calculate_geometry_score(best_model)

        # Calculate consistency across models
        consistency_score = self._calculate_consistency(interface_analyses)

        # Biological scores
        expression_score = self._score_expression(candidate_data)
        localization_score = self._score_localization(candidate_data)
        conservation_score = self._score_conservation(conservation_data)

        # Composite scores
        structural_confidence = self._calculate_structural_confidence(
            ipTM_score, interface_pLDDT_score, interface_area_score,
            interaction_score, geometry_score, consistency_score
        )

        biological_relevance = self._calculate_biological_relevance(
            expression_score, localization_score, conservation_score
        )

        # Overall score
        overall_score = self._calculate_overall_score(
            ipTM_score,
            interaction_score,
            interface_area_score,
            consistency_score,
            expression_score,
            localization_score,
            conservation_score
        )

        # Determine confidence tier
        confidence_tier = self._determine_confidence_tier(
            ipTM_score, interface_pLDDT_score, overall_score
        )

        return CandidateScore(
            complex_id=complex_id,
            protein_id=protein_id,
            gene_name=gene_name,
            ipTM_score=ipTM_score,
            pTM_score=pTM_score,
            interface_pLDDT_score=interface_pLDDT_score,
            interface_area_score=interface_area_score,
            interaction_score=interaction_score,
            geometry_score=geometry_score,
            consistency_score=consistency_score,
            expression_score=expression_score,
            localization_score=localization_score,
            conservation_score=conservation_score,
            structural_confidence=structural_confidence,
            biological_relevance=biological_relevance,
            overall_score=overall_score,
            num_models=len(interface_analyses),
            best_model=best_model.get('model_number', 1),
            confidence_tier=confidence_tier
        )

    def _select_best_model(self, analyses: List[Dict]) -> Dict:
        """Select best model based on ipTM and interface quality"""
        if not analyses:
            return {}

        # Sort by ipTM (primary) and interface_pLDDT (secondary)
        sorted_analyses = sorted(
            analyses,
            key=lambda x: (x.get('ipTM', 0), x.get('interface_pLDDT', 0)),
            reverse=True
        )

        return sorted_analyses[0]

    def _normalize_area(self, area: float) -> float:
        """
        Normalize interface area to 0-1 score

        Typical protein-protein interfaces: 800-2000 Å²
        """
        min_area = 500.0
        max_area = 2000.0

        if area < min_area:
            return area / min_area * 0.5  # Penalize small interfaces
        elif area > max_area:
            return 1.0
        else:
            return 0.5 + (area - min_area) / (max_area - min_area) * 0.5

    def _calculate_interaction_score(self, analysis: Dict) -> float:
        """
        Score based on interaction types

        H-bonds, salt bridges, hydrophobic contacts
        """
        h_bonds = analysis.get('hydrogen_bonds', 0)
        salt_bridges = analysis.get('salt_bridges', 0)
        hydrophobic = analysis.get('hydrophobic_contacts', 0)

        # Typical ranges (approximate)
        h_bond_score = min(1.0, h_bonds / 10.0)
        salt_bridge_score = min(1.0, salt_bridges / 5.0)
        hydrophobic_score = min(1.0, hydrophobic / 20.0)

        # Weighted combination
        interaction_score = (
            0.4 * h_bond_score +
            0.3 * salt_bridge_score +
            0.3 * hydrophobic_score
        )

        return interaction_score

    def _calculate_geometry_score(self, analysis: Dict) -> float:
        """Score based on interface geometry"""
        gap_score = analysis.get('interface_gap_score', 0.0)
        shape_comp = analysis.get('shape_complementarity', 0.0)

        return (gap_score + shape_comp) / 2.0

    def _calculate_consistency(self, analyses: List[Dict]) -> float:
        """
        Calculate consistency across multiple models

        High consistency = same interface residues across models
        """
        if len(analyses) < 2:
            return 0.5  # Neutral score if only one model

        # Compare ipTM scores across models
        ipTMs = [a.get('ipTM', 0) for a in analyses]
        ipTM_std = np.std(ipTMs)

        # Low standard deviation = high consistency
        # ipTM std < 0.1 = very consistent
        consistency = max(0.0, 1.0 - ipTM_std / 0.2)

        return consistency

    def _score_expression(self, candidate_data: Dict) -> float:
        """Score based on intestinal expression"""
        expression = candidate_data.get('intestinal_expression')

        if expression is None:
            return 0.5  # Neutral score if unknown

        # Assuming TPM values
        # High expression: > 50 TPM
        # Medium: 10-50 TPM
        # Low: < 10 TPM

        if expression > 50:
            return 1.0
        elif expression > 10:
            return 0.7
        elif expression > 1:
            return 0.4
        else:
            return 0.2

    def _score_localization(self, candidate_data: Dict) -> float:
        """Score based on subcellular localization"""
        localization = set(candidate_data.get('localization', []))

        if 'apical' in localization:
            return 1.0  # Ideal - apical membrane
        elif 'membrane' in localization:
            return 0.7  # Good - membrane, but not specifically apical
        elif 'basolateral' in localization:
            return 0.3  # Less likely - wrong side
        else:
            return 0.5  # Unknown

    def _score_conservation(self, conservation_data: Optional[Dict]) -> float:
        """Score based on binding site conservation"""
        if conservation_data is None:
            return 0.5  # Neutral if no data

        # Would analyze conservation across viral strains
        # For now, return placeholder
        return 0.7

    def _calculate_structural_confidence(
        self,
        ipTM: float,
        pLDDT: float,
        area: float,
        interactions: float,
        geometry: float,
        consistency: float
    ) -> float:
        """Calculate structural confidence composite score"""
        return (
            0.35 * ipTM +
            0.20 * pLDDT +
            0.15 * area +
            0.15 * interactions +
            0.10 * geometry +
            0.05 * consistency
        )

    def _calculate_biological_relevance(
        self,
        expression: float,
        localization: float,
        conservation: float
    ) -> float:
        """Calculate biological relevance composite score"""
        return (
            0.4 * expression +
            0.4 * localization +
            0.2 * conservation
        )

    def _calculate_overall_score(
        self,
        ipTM: float,
        interactions: float,
        area: float,
        consistency: float,
        expression: float,
        localization: float,
        conservation: float
    ) -> float:
        """Calculate overall score using configured weights"""
        w = self.weights

        score = (
            w.ipTM_weight * ipTM +
            w.interface_quality_weight * interactions +
            w.interface_area_weight * area +
            w.consistency_weight * consistency +
            w.expression_weight * expression +
            w.localization_weight * localization +
            w.conservation_weight * conservation
        )

        return score

    def _determine_confidence_tier(
        self,
        ipTM: float,
        pLDDT: float,
        overall: float
    ) -> str:
        """Determine confidence tier"""
        if ipTM > 0.8 and pLDDT > 0.8 and overall > 0.75:
            return 'high'
        elif ipTM > 0.6 and pLDDT > 0.7 and overall > 0.6:
            return 'medium'
        else:
            return 'low'

    def rank_candidates(self, scores: List[CandidateScore]) -> List[CandidateScore]:
        """
        Rank candidates by overall score

        Returns sorted list with rank assigned
        """
        # Sort by overall score (descending)
        sorted_scores = sorted(
            scores,
            key=lambda x: x.overall_score,
            reverse=True
        )

        # Assign ranks
        for i, score in enumerate(sorted_scores, 1):
            score.rank = i

        return sorted_scores


def score_all_candidates(
    interface_analyses_file: Path,
    candidate_library_file: Path,
    output_file: Path,
    weights: Optional[ScoringWeights] = None
) -> List[CandidateScore]:
    """
    Score all candidates and generate ranked list

    Args:
        interface_analyses_file: JSON file with interface analyses
        candidate_library_file: JSON file with candidate metadata
        output_file: Where to save scored results
        weights: Optional custom scoring weights

    Returns:
        Ranked list of CandidateScore objects
    """
    logger.info("Scoring all candidates...")

    # Load data
    with open(interface_analyses_file) as f:
        analyses = json.load(f)

    with open(candidate_library_file) as f:
        candidates = json.load(f)

    # Create mapping of protein_id to candidate data
    candidate_map = {c['uniprot_id']: c for c in candidates}

    # Group analyses by complex_id
    analyses_by_complex = {}
    for analysis in analyses:
        complex_id = analysis['complex_id']
        if complex_id not in analyses_by_complex:
            analyses_by_complex[complex_id] = []
        analyses_by_complex[complex_id].append(analysis)

    # Score each candidate
    scorer = CandidateScorer(weights)
    scores = []

    for complex_id, complex_analyses in analyses_by_complex.items():
        # Extract protein_id from complex_id (format: GII17_VP1_PROTEIN_ID)
        parts = complex_id.split('_')
        if len(parts) >= 3:
            protein_id = '_'.join(parts[2:])
        else:
            protein_id = complex_id

        if protein_id not in candidate_map:
            logger.warning(f"No candidate data for {protein_id}")
            continue

        candidate_data = candidate_map[protein_id]
        gene_name = candidate_data.get('gene_name', 'Unknown')

        score = scorer.score_candidate(
            complex_id=complex_id,
            protein_id=protein_id,
            gene_name=gene_name,
            interface_analyses=complex_analyses,
            candidate_data=candidate_data,
            conservation_data=None  # Would add if available
        )

        scores.append(score)

    # Rank candidates
    ranked_scores = scorer.rank_candidates(scores)

    # Save results
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(
            [asdict(s) for s in ranked_scores],
            f,
            indent=2
        )

    logger.info(f"Scored and ranked {len(ranked_scores)} candidates")
    logger.info(f"Results saved to {output_file}")

    # Print top 10
    print("\n=== Top 10 Candidates ===\n")
    for score in ranked_scores[:10]:
        print(f"{score.rank}. {score.gene_name} ({score.protein_id})")
        print(f"   Overall Score: {score.overall_score:.3f}")
        print(f"   ipTM: {score.ipTM_score:.3f}")
        print(f"   Confidence: {score.confidence_tier}")
        print()

    return ranked_scores


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Score and rank receptor candidates"
    )
    parser.add_argument(
        "--interface-analyses",
        type=Path,
        required=True,
        help="JSON file with interface analyses"
    )
    parser.add_argument(
        "--candidate-library",
        type=Path,
        required=True,
        help="JSON file with candidate library"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default="data/results/ranked_candidates.json",
        help="Output file for ranked candidates"
    )

    args = parser.parse_args()

    ranked = score_all_candidates(
        args.interface_analyses,
        args.candidate_library,
        args.output
    )

    print(f"\n✓ Scored and ranked {len(ranked)} candidates")
    print(f"✓ Results saved to {args.output}")


if __name__ == "__main__":
    main()
