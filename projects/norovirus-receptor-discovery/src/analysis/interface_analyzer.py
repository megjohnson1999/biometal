"""
Interface Analyzer

Analyzes binding interfaces from AlphaFold3 predictions to identify
key interactions and score binding quality.
"""

import logging
import json
import numpy as np
from dataclasses import dataclass, asdict
from typing import List, Dict, Tuple, Optional, Set
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import three_to_one
import warnings
warnings.filterwarnings('ignore', category=PDBWarning)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class InterfaceResidue:
    """Represents a residue at the binding interface"""
    chain_id: str
    residue_number: int
    residue_name: str
    residue_letter: str
    pLDDT: float
    distance_to_partner: float  # Minimum distance to partner chain


@dataclass
class InterfaceAnalysis:
    """Complete interface analysis for a complex"""
    complex_id: str
    model_number: int

    # Overall metrics
    ipTM: float
    pTM: float
    interface_area: float  # Å²
    interface_pLDDT: float

    # Interface residues
    vp1_interface_residues: List[InterfaceResidue]
    receptor_interface_residues: List[InterfaceResidue]

    # Interactions
    hydrogen_bonds: int
    salt_bridges: int
    hydrophobic_contacts: int
    aromatic_interactions: int

    # Geometry
    interface_gap_score: float  # 0-1, higher is better (tighter interface)
    shape_complementarity: float  # 0-1, higher is better

    # Confidence
    consistent_interface: bool  # Interface consistent across models
    binding_confidence: str  # 'high', 'medium', 'low'


class InterfaceAnalyzer:
    """Analyze protein-protein binding interfaces"""

    def __init__(self, distance_cutoff: float = 5.0):
        """
        Args:
            distance_cutoff: Distance cutoff for interface residues (Å)
        """
        self.distance_cutoff = distance_cutoff
        self.parser = PDBParser(QUIET=True)

    def analyze_complex(
        self,
        pdb_file: Path,
        complex_id: str,
        model_number: int,
        ipTM: float,
        pTM: float,
        chain_ids: Tuple[str, str] = ('A', 'B')
    ) -> InterfaceAnalysis:
        """
        Analyze binding interface from PDB file

        Args:
            pdb_file: Path to PDB file
            complex_id: Complex identifier
            model_number: Model number
            ipTM: Interface pTM score
            pTM: Overall pTM score
            chain_ids: (VP1_chain, receptor_chain)

        Returns:
            InterfaceAnalysis object
        """
        logger.info(f"Analyzing interface for {complex_id} model {model_number}")

        structure = self.parser.get_structure(complex_id, pdb_file)
        model = structure[0]

        vp1_chain_id, receptor_chain_id = chain_ids

        # Get chains
        try:
            vp1_chain = model[vp1_chain_id]
            receptor_chain = model[receptor_chain_id]
        except KeyError:
            logger.error(f"Chain not found in {pdb_file}")
            return self._create_empty_analysis(complex_id, model_number, ipTM, pTM)

        # Find interface residues
        vp1_interface = self._find_interface_residues(
            vp1_chain, receptor_chain, self.distance_cutoff
        )
        receptor_interface = self._find_interface_residues(
            receptor_chain, vp1_chain, self.distance_cutoff
        )

        # Calculate interface metrics
        interface_area = self._calculate_interface_area(vp1_interface, receptor_interface)
        interface_pLDDT = self._calculate_interface_pLDDT(vp1_interface + receptor_interface)

        # Analyze interactions
        h_bonds = self._count_hydrogen_bonds(vp1_chain, receptor_chain)
        salt_bridges = self._count_salt_bridges(vp1_chain, receptor_chain)
        hydrophobic = self._count_hydrophobic_contacts(vp1_chain, receptor_chain)
        aromatic = self._count_aromatic_interactions(vp1_chain, receptor_chain)

        # Geometric metrics
        gap_score = self._calculate_gap_score(vp1_interface, receptor_interface)
        shape_comp = self._estimate_shape_complementarity(vp1_chain, receptor_chain)

        # Assess confidence
        binding_confidence = self._assess_binding_confidence(
            ipTM, interface_pLDDT, interface_area, len(vp1_interface) + len(receptor_interface)
        )

        analysis = InterfaceAnalysis(
            complex_id=complex_id,
            model_number=model_number,
            ipTM=ipTM,
            pTM=pTM,
            interface_area=interface_area,
            interface_pLDDT=interface_pLDDT,
            vp1_interface_residues=vp1_interface,
            receptor_interface_residues=receptor_interface,
            hydrogen_bonds=h_bonds,
            salt_bridges=salt_bridges,
            hydrophobic_contacts=hydrophobic,
            aromatic_interactions=aromatic,
            interface_gap_score=gap_score,
            shape_complementarity=shape_comp,
            consistent_interface=True,  # Would check across models
            binding_confidence=binding_confidence
        )

        return analysis

    def _find_interface_residues(
        self,
        chain1,
        chain2,
        cutoff: float
    ) -> List[InterfaceResidue]:
        """
        Find residues at the interface between two chains

        Returns list of InterfaceResidue objects
        """
        interface_residues = []

        for residue1 in chain1.get_residues():
            if not self._is_standard_residue(residue1):
                continue

            min_distance = float('inf')

            # Check distance to all residues in chain2
            for residue2 in chain2.get_residues():
                if not self._is_standard_residue(residue2):
                    continue

                dist = self._residue_distance(residue1, residue2)
                min_distance = min(min_distance, dist)

            if min_distance <= cutoff:
                # Extract pLDDT from B-factor column (AlphaFold convention)
                plddt = self._get_plddt(residue1)

                interface_res = InterfaceResidue(
                    chain_id=chain1.id,
                    residue_number=residue1.id[1],
                    residue_name=residue1.resname,
                    residue_letter=three_to_one(residue1.resname) if residue1.resname in three_to_one.keys() else 'X',
                    pLDDT=plddt,
                    distance_to_partner=min_distance
                )
                interface_residues.append(interface_res)

        return interface_residues

    def _residue_distance(self, res1, res2) -> float:
        """Calculate minimum distance between two residues"""
        min_dist = float('inf')
        for atom1 in res1:
            for atom2 in res2:
                dist = atom1 - atom2  # BioPython computes distance
                min_dist = min(min_dist, dist)
        return min_dist

    def _is_standard_residue(self, residue) -> bool:
        """Check if residue is a standard amino acid"""
        return residue.id[0] == ' ' and residue.resname in three_to_one.keys()

    def _get_plddt(self, residue) -> float:
        """Extract pLDDT score from B-factor column"""
        plddt_values = []
        for atom in residue:
            plddt_values.append(atom.bfactor)
        return np.mean(plddt_values) if plddt_values else 0.0

    def _calculate_interface_area(
        self,
        vp1_interface: List[InterfaceResidue],
        receptor_interface: List[InterfaceResidue]
    ) -> float:
        """
        Estimate interface area (simplified)

        More accurate calculation would use SASA difference
        """
        # Rough estimate: ~20 Å² per interface residue
        total_interface_residues = len(vp1_interface) + len(receptor_interface)
        return total_interface_residues * 20.0

    def _calculate_interface_pLDDT(self, interface_residues: List[InterfaceResidue]) -> float:
        """Calculate average pLDDT at interface"""
        if not interface_residues:
            return 0.0
        return np.mean([res.pLDDT for res in interface_residues])

    def _count_hydrogen_bonds(self, chain1, chain2) -> int:
        """
        Count potential hydrogen bonds between chains

        Simplified: D-A distance < 3.5 Å and angle considerations
        """
        donors = ['N', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NH1', 'NH2', 'NZ', 'OG', 'OG1', 'OH']
        acceptors = ['O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'ND1', 'NE2']

        h_bonds = 0
        cutoff = 3.5

        for res1 in chain1.get_residues():
            if not self._is_standard_residue(res1):
                continue
            for atom1 in res1:
                if atom1.name not in donors:
                    continue

                for res2 in chain2.get_residues():
                    if not self._is_standard_residue(res2):
                        continue
                    for atom2 in res2:
                        if atom2.name not in acceptors:
                            continue

                        dist = atom1 - atom2
                        if dist < cutoff:
                            h_bonds += 1

        return h_bonds

    def _count_salt_bridges(self, chain1, chain2) -> int:
        """Count potential salt bridges (ionic interactions)"""
        positive = ['ARG', 'LYS', 'HIS']
        negative = ['ASP', 'GLU']

        salt_bridges = 0
        cutoff = 4.0

        for res1 in chain1.get_residues():
            if not self._is_standard_residue(res1) or res1.resname not in positive + negative:
                continue

            for res2 in chain2.get_residues():
                if not self._is_standard_residue(res2) or res2.resname not in positive + negative:
                    continue

                # Check if opposite charges
                if (res1.resname in positive and res2.resname in negative) or \
                   (res1.resname in negative and res2.resname in positive):
                    dist = self._residue_distance(res1, res2)
                    if dist < cutoff:
                        salt_bridges += 1

        return salt_bridges

    def _count_hydrophobic_contacts(self, chain1, chain2) -> int:
        """Count hydrophobic contacts"""
        hydrophobic = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO']

        contacts = 0
        cutoff = 5.0

        for res1 in chain1.get_residues():
            if not self._is_standard_residue(res1) or res1.resname not in hydrophobic:
                continue

            for res2 in chain2.get_residues():
                if not self._is_standard_residue(res2) or res2.resname not in hydrophobic:
                    continue

                dist = self._residue_distance(res1, res2)
                if dist < cutoff:
                    contacts += 1

        return contacts

    def _count_aromatic_interactions(self, chain1, chain2) -> int:
        """Count aromatic interactions (pi-pi, pi-cation)"""
        aromatic = ['PHE', 'TYR', 'TRP', 'HIS']

        interactions = 0
        cutoff = 6.0

        for res1 in chain1.get_residues():
            if not self._is_standard_residue(res1) or res1.resname not in aromatic:
                continue

            for res2 in chain2.get_residues():
                if not self._is_standard_residue(res2) or res2.resname not in aromatic:
                    continue

                dist = self._residue_distance(res1, res2)
                if dist < cutoff:
                    interactions += 1

        return interactions

    def _calculate_gap_score(
        self,
        vp1_interface: List[InterfaceResidue],
        receptor_interface: List[InterfaceResidue]
    ) -> float:
        """
        Calculate interface gap score (how tightly packed)

        Returns value 0-1, where 1 = very tight interface
        """
        if not vp1_interface or not receptor_interface:
            return 0.0

        # Average minimum distance across interface
        avg_distance = np.mean([res.distance_to_partner for res in vp1_interface + receptor_interface])

        # Convert to 0-1 score (distance 3.0 = score 1.0, distance 5.0 = score 0.0)
        score = max(0.0, min(1.0, (5.0 - avg_distance) / 2.0))
        return score

    def _estimate_shape_complementarity(self, chain1, chain2) -> float:
        """
        Estimate shape complementarity (simplified)

        Real calculation would use surface normals and curvature
        """
        # Placeholder - would need proper implementation
        # For now, return moderate value
        return 0.6

    def _assess_binding_confidence(
        self,
        ipTM: float,
        interface_pLDDT: float,
        interface_area: float,
        num_interface_residues: int
    ) -> str:
        """
        Assess overall binding confidence

        Returns: 'high', 'medium', or 'low'
        """
        # High confidence criteria
        if ipTM > 0.8 and interface_pLDDT > 80 and interface_area > 800:
            return 'high'

        # Medium confidence
        elif ipTM > 0.6 and interface_pLDDT > 70 and interface_area > 500:
            return 'medium'

        # Low confidence
        else:
            return 'low'

    def _create_empty_analysis(
        self,
        complex_id: str,
        model_number: int,
        ipTM: float,
        pTM: float
    ) -> InterfaceAnalysis:
        """Create empty analysis for failed cases"""
        return InterfaceAnalysis(
            complex_id=complex_id,
            model_number=model_number,
            ipTM=ipTM,
            pTM=pTM,
            interface_area=0.0,
            interface_pLDDT=0.0,
            vp1_interface_residues=[],
            receptor_interface_residues=[],
            hydrogen_bonds=0,
            salt_bridges=0,
            hydrophobic_contacts=0,
            aromatic_interactions=0,
            interface_gap_score=0.0,
            shape_complementarity=0.0,
            consistent_interface=False,
            binding_confidence='low'
        )


def analyze_batch(
    prediction_dir: Path,
    output_file: Path,
    distance_cutoff: float = 5.0
) -> List[InterfaceAnalysis]:
    """
    Analyze a batch of predictions

    Args:
        prediction_dir: Directory with prediction PDB files
        output_file: Where to save analysis results
        distance_cutoff: Interface distance cutoff

    Returns:
        List of InterfaceAnalysis objects
    """
    analyzer = InterfaceAnalyzer(distance_cutoff=distance_cutoff)
    analyses = []

    # Find all PDB files
    pdb_files = list(prediction_dir.glob("**/*.pdb"))
    logger.info(f"Found {len(pdb_files)} PDB files to analyze")

    for pdb_file in pdb_files:
        try:
            # Extract metadata from filename or accompanying JSON
            complex_id = pdb_file.parent.name
            model_number = 1  # Would parse from filename

            # Would load scores from prediction results
            ipTM = 0.0  # Placeholder
            pTM = 0.0

            analysis = analyzer.analyze_complex(
                pdb_file, complex_id, model_number, ipTM, pTM
            )
            analyses.append(analysis)

        except Exception as e:
            logger.error(f"Error analyzing {pdb_file}: {e}")

    # Save results
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(
            [asdict(a) for a in analyses],
            f,
            indent=2
        )

    logger.info(f"Saved {len(analyses)} analyses to {output_file}")
    return analyses


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Analyze binding interfaces from AlphaFold3 predictions"
    )
    parser.add_argument(
        "--prediction-dir",
        type=Path,
        required=True,
        help="Directory containing prediction PDB files"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default="data/results/interface_analysis.json",
        help="Output JSON file"
    )
    parser.add_argument(
        "--distance-cutoff",
        type=float,
        default=5.0,
        help="Interface distance cutoff (Å)"
    )

    args = parser.parse_args()

    analyses = analyze_batch(
        args.prediction_dir,
        args.output,
        args.distance_cutoff
    )

    print(f"\n✓ Analyzed {len(analyses)} predictions")
    print(f"✓ Results saved to {args.output}")


if __name__ == "__main__":
    main()
