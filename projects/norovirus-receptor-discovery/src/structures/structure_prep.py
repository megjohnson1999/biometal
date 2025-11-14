"""
Structure Preparation for AlphaFold3 Predictions

Extracts ectodomains, predicts topology, and prepares sequences for structure prediction.
"""

import logging
import json
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
from pathlib import Path
import subprocess
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class EctodomainSequence:
    """Represents an extracted ectodomain sequence"""
    protein_id: str
    gene_name: str
    full_sequence: str
    ectodomain_sequence: str
    ectodomain_regions: List[Tuple[int, int]]  # In 1-based coordinates
    ectodomain_start: int
    ectodomain_end: int
    topology_method: str  # 'uniprot', 'deeptmhmm', 'tmhmm'
    confidence: float

    @property
    def ectodomain_length(self) -> int:
        return len(self.ectodomain_sequence)


class TopologyPredictor:
    """Predict membrane topology using various tools"""

    def __init__(self, tool: str = "deeptmhmm"):
        """
        Args:
            tool: Topology prediction tool ('deeptmhmm', 'tmhmm', 'phobius')
        """
        self.tool = tool

    def predict_topology(self, sequence: str) -> Dict:
        """
        Predict membrane topology

        Returns:
            Dictionary with:
            - tm_regions: List of (start, end) for TM regions
            - ectodomain_regions: List of (start, end) for extracellular regions
            - topology_type: 'single_pass', 'multi_pass', etc.
            - confidence: Prediction confidence
        """
        if self.tool == "deeptmhmm":
            return self._predict_deeptmhmm(sequence)
        elif self.tool == "tmhmm":
            return self._predict_tmhmm(sequence)
        else:
            raise ValueError(f"Unknown tool: {self.tool}")

    def _predict_deeptmhmm(self, sequence: str) -> Dict:
        """
        Predict using DeepTMHMM (biolib or local installation)

        Note: This is a placeholder. In practice, you would either:
        1. Use biolib API
        2. Install DeepTMHMM locally
        3. Use a web service
        """
        logger.warning("DeepTMHMM prediction not implemented - using simple heuristic")
        return self._simple_topology_heuristic(sequence)

    def _predict_tmhmm(self, sequence: str) -> Dict:
        """
        Predict using TMHMM (requires local installation)

        Note: This requires TMHMM to be installed and in PATH
        """
        try:
            # Write sequence to temp file
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
                f.write(f">seq\n{sequence}\n")
                temp_file = f.name

            # Run TMHMM
            result = subprocess.run(
                ['tmhmm', temp_file],
                capture_output=True,
                text=True,
                timeout=60
            )

            # Parse output (simplified)
            return self._parse_tmhmm_output(result.stdout)

        except FileNotFoundError:
            logger.warning("TMHMM not found - using simple heuristic")
            return self._simple_topology_heuristic(sequence)
        except Exception as e:
            logger.warning(f"TMHMM prediction failed: {e}")
            return self._simple_topology_heuristic(sequence)

    def _simple_topology_heuristic(self, sequence: str) -> Dict:
        """
        Simple hydrophobicity-based topology prediction
        This is a fallback when specialized tools aren't available
        """
        # Kyte-Doolittle hydrophobicity scale
        hydrophobicity = {
            'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
            'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
            'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
            'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
        }

        window_size = 19
        threshold = 1.6

        tm_regions = []
        ectodomain_regions = []

        # Sliding window hydrophobicity
        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i + window_size]
            avg_hydro = sum(hydrophobicity.get(aa, 0) for aa in window) / window_size

            if avg_hydro > threshold:
                # Potential TM region
                tm_regions.append((i + 1, i + window_size))  # 1-based

        # Merge overlapping TM regions
        if tm_regions:
            merged_tm = [tm_regions[0]]
            for start, end in tm_regions[1:]:
                if start <= merged_tm[-1][1] + 5:  # Allow 5 aa gap
                    merged_tm[-1] = (merged_tm[-1][0], max(merged_tm[-1][1], end))
                else:
                    merged_tm.append((start, end))
            tm_regions = merged_tm

        # Identify extracellular regions (before first TM and between TMs)
        if tm_regions:
            # N-terminal region
            if tm_regions[0][0] > 30:  # At least 30 aa N-terminal
                ectodomain_regions.append((1, tm_regions[0][0] - 1))

            # Between TM regions (assume alternating topology)
            for i in range(len(tm_regions) - 1):
                region_start = tm_regions[i][1] + 1
                region_end = tm_regions[i + 1][0] - 1
                if region_end - region_start > 30:  # At least 30 aa
                    # Odd-numbered loops are typically extracellular
                    if i % 2 == 0:
                        ectodomain_regions.append((region_start, region_end))

        topology_type = (
            "single_pass" if len(tm_regions) == 1
            else "multi_pass" if len(tm_regions) > 1
            else "soluble"
        )

        return {
            "tm_regions": tm_regions,
            "ectodomain_regions": ectodomain_regions,
            "topology_type": topology_type,
            "confidence": 0.6,  # Low confidence for heuristic
            "method": "hydrophobicity_heuristic"
        }

    def _parse_tmhmm_output(self, output: str) -> Dict:
        """Parse TMHMM output"""
        # Simplified parser - would need to be more robust
        lines = output.strip().split('\n')
        tm_regions = []
        ectodomain_regions = []

        for line in lines:
            if 'TMhelix' in line:
                parts = line.split()
                start = int(parts[2])
                end = int(parts[3])
                tm_regions.append((start, end))
            elif 'outside' in line.lower():
                parts = line.split()
                start = int(parts[2])
                end = int(parts[3])
                if end - start > 30:
                    ectodomain_regions.append((start, end))

        topology_type = (
            "single_pass" if len(tm_regions) == 1
            else "multi_pass" if len(tm_regions) > 1
            else "soluble"
        )

        return {
            "tm_regions": tm_regions,
            "ectodomain_regions": ectodomain_regions,
            "topology_type": topology_type,
            "confidence": 0.85,
            "method": "tmhmm"
        }


class EctodomainExtractor:
    """Extract ectodomain sequences from full-length proteins"""

    def __init__(self, topology_predictor: Optional[TopologyPredictor] = None):
        self.topology_predictor = topology_predictor or TopologyPredictor()

    def extract_ectodomain(
        self,
        protein_id: str,
        gene_name: str,
        sequence: str,
        known_regions: Optional[List[Tuple[int, int]]] = None
    ) -> Optional[EctodomainSequence]:
        """
        Extract ectodomain from protein sequence

        Args:
            protein_id: UniProt ID or identifier
            gene_name: Gene name
            sequence: Full-length protein sequence
            known_regions: Optional known ectodomain regions from UniProt

        Returns:
            EctodomainSequence object or None if no suitable ectodomain
        """
        if known_regions and len(known_regions) > 0:
            # Use known regions
            method = "uniprot"
            confidence = 0.95
            ectodomain_regions = known_regions
        else:
            # Predict topology
            topology = self.topology_predictor.predict_topology(sequence)
            method = topology["method"]
            confidence = topology["confidence"]
            ectodomain_regions = topology["ectodomain_regions"]

        if not ectodomain_regions:
            logger.warning(f"No ectodomain found for {protein_id}")
            return None

        # Extract longest ectodomain region
        longest_region = max(ectodomain_regions, key=lambda x: x[1] - x[0])
        start, end = longest_region

        # Extract sequence (convert to 0-based indexing)
        ectodomain_seq = sequence[start - 1:end]

        if len(ectodomain_seq) < 50:
            logger.warning(f"Ectodomain too short for {protein_id}: {len(ectodomain_seq)} aa")
            return None

        return EctodomainSequence(
            protein_id=protein_id,
            gene_name=gene_name,
            full_sequence=sequence,
            ectodomain_sequence=ectodomain_seq,
            ectodomain_regions=ectodomain_regions,
            ectodomain_start=start,
            ectodomain_end=end,
            topology_method=method,
            confidence=confidence
        )


class VP1StructureRetriever:
    """Retrieve or prepare GII.17 VP1 P-domain structure"""

    def __init__(self):
        self.pdb_base = "https://files.rcsb.org/download"

    def get_vp1_sequence(self, strain: str = "GII.17") -> str:
        """
        Get VP1 P-domain sequence for specified strain

        For GII.17, the P-domain is approximately residues 225-530
        """
        # GII.17 VP1 sequence (example - would need actual sequence)
        # This is a placeholder - in practice, retrieve from GenBank
        gii17_vp1_full = """
MKMASNDAAPSNDGAAGLVPEINNELKVAGQPLQSTVANGSIYAGLAVKAGFDFRLRDE
VVDGDLLGTTQLSPVNICTFRGDVTHIAGRQLAIPHATGVAQFQGNSSDVTTGFTPDFN
KFLVPPTVESKLKPFTLPILSSQSGALTRRFPVPHLDFPFLTQTNPVHDFRTTMKPDAS
ALSVPGTNDDLGQNQIHQPGFSKSAWALAPFPDPGPPGSGDVFFRLDNNWLQQVDPSVQ
TQIGGDLISVGVPFNDPSQSQGWFSGLDYLPQSGHVYGQNFSPKVRVMIVDGGQPTPSE
IIDGFTNGWFTWVSARNSNRFTYPVGQTHNGGLLNPTTSSLISGFTGALRQGPFQFPHD
PFPTWSQPYNYGLASPEYSPNNPDEWVQTLNLSIAQNLGFGSGFFSTLGFPNQLPYIKP
SPTTGLLVEGRDLLPVHLDLEYTFAELAPNGVYFDIGAQDGQFSVPIEIDITGQNSDQL
KLHTPFQGPQGTWQNWQGLPQDPGDMFGVDVRNLTPPGQS
"""
        # P-domain is approximately residues 225-530
        p_domain_start = 224  # 0-based
        p_domain_end = 530

        full_seq = gii17_vp1_full.replace("\n", "").strip()
        p_domain_seq = full_seq[p_domain_start:p_domain_end]

        return p_domain_seq

    def fetch_pdb_structure(self, pdb_id: str, output_path: Path) -> bool:
        """
        Fetch PDB structure from RCSB

        Args:
            pdb_id: PDB identifier (e.g., '5J4K')
            output_path: Where to save the PDB file

        Returns:
            True if successful
        """
        url = f"{self.pdb_base}/{pdb_id}.pdb"

        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()

            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, 'w') as f:
                f.write(response.text)

            logger.info(f"Downloaded PDB {pdb_id} to {output_path}")
            return True

        except Exception as e:
            logger.error(f"Failed to download PDB {pdb_id}: {e}")
            return False


class AlphaFold3InputPrep:
    """Prepare input files for AlphaFold3 predictions"""

    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def prepare_complex_input(
        self,
        vp1_sequence: str,
        ectodomain: EctodomainSequence,
        complex_id: Optional[str] = None
    ) -> Dict:
        """
        Prepare AlphaFold3 input for VP1 + ectodomain complex

        Args:
            vp1_sequence: VP1 P-domain sequence
            ectodomain: Ectodomain sequence object
            complex_id: Optional identifier for this complex

        Returns:
            Dictionary with input information
        """
        if complex_id is None:
            complex_id = f"GII17_VP1_{ectodomain.protein_id}"

        # Prepare JSON input (AlphaFold3 format)
        af3_input = {
            "name": complex_id,
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": vp1_sequence,
                        "id": "VP1_P_domain"
                    }
                },
                {
                    "proteinChain": {
                        "sequence": ectodomain.ectodomain_sequence,
                        "id": f"{ectodomain.gene_name}_ectodomain"
                    }
                }
            ],
            "modelSeeds": [1, 2, 3, 4, 5],  # Generate 5 models
            "dialect": "alphafold3",
            "version": 1
        }

        # Save JSON file
        json_file = self.output_dir / f"{complex_id}.json"
        with open(json_file, 'w') as f:
            json.dump(af3_input, f, indent=2)

        # Also save as FASTA (for other tools)
        fasta_file = self.output_dir / f"{complex_id}.fasta"
        records = [
            SeqRecord(Seq(vp1_sequence), id="VP1_P_domain", description="GII.17 VP1 P-domain"),
            SeqRecord(
                Seq(ectodomain.ectodomain_sequence),
                id=f"{ectodomain.gene_name}_ectodomain",
                description=f"{ectodomain.protein_id} ectodomain {ectodomain.ectodomain_start}-{ectodomain.ectodomain_end}"
            )
        ]
        SeqIO.write(records, fasta_file, "fasta")

        return {
            "complex_id": complex_id,
            "json_input": str(json_file),
            "fasta_input": str(fasta_file),
            "vp1_length": len(vp1_sequence),
            "ectodomain_length": ectodomain.ectodomain_length
        }


class StructurePreparationPipeline:
    """Main pipeline for structure preparation"""

    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.extractor = EctodomainExtractor()
        self.vp1_retriever = VP1StructureRetriever()
        self.input_prep = AlphaFold3InputPrep(output_dir / "af3_inputs")

    def prepare_from_library(self, library_file: Path) -> List[Dict]:
        """
        Prepare AlphaFold3 inputs from candidate library

        Args:
            library_file: Path to candidate_library.json

        Returns:
            List of prepared input dictionaries
        """
        logger.info(f"Loading candidate library from {library_file}")

        with open(library_file) as f:
            library = json.load(f)

        # Get VP1 sequence
        vp1_sequence = self.vp1_retriever.get_vp1_sequence("GII.17")
        logger.info(f"VP1 P-domain length: {len(vp1_sequence)} aa")

        # Process each candidate
        prepared_inputs = []
        ectodomains = []

        for candidate in library:
            try:
                # Extract ectodomain
                ectodomain = self.extractor.extract_ectodomain(
                    protein_id=candidate["uniprot_id"],
                    gene_name=candidate["gene_name"],
                    sequence=candidate["sequence"],
                    known_regions=[(r[0], r[1]) for r in candidate.get("ectodomain_regions", [])]
                )

                if ectodomain:
                    ectodomains.append(ectodomain)

                    # Prepare AF3 input
                    input_info = self.input_prep.prepare_complex_input(
                        vp1_sequence, ectodomain
                    )
                    prepared_inputs.append(input_info)

            except Exception as e:
                logger.error(f"Error processing {candidate['uniprot_id']}: {e}")

        logger.info(f"Prepared {len(prepared_inputs)} AlphaFold3 inputs")

        # Save ectodomain information
        ectodomain_file = self.output_dir / "ectodomains.json"
        with open(ectodomain_file, 'w') as f:
            json.dump(
                [vars(e) for e in ectodomains],
                f,
                indent=2
            )

        # Save summary
        summary_file = self.output_dir / "preparation_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(
                {
                    "total_candidates": len(library),
                    "successful_preparations": len(prepared_inputs),
                    "vp1_length": len(vp1_sequence),
                    "inputs": prepared_inputs
                },
                f,
                indent=2
            )

        return prepared_inputs


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Prepare structures for AlphaFold3 predictions"
    )
    parser.add_argument(
        "--library",
        type=Path,
        required=True,
        help="Path to candidate_library.json"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default="data/processed/structures",
        help="Output directory"
    )

    args = parser.parse_args()

    pipeline = StructurePreparationPipeline(args.output_dir)
    prepared = pipeline.prepare_from_library(args.library)

    print(f"\n✓ Prepared {len(prepared)} AlphaFold3 inputs")
    print(f"✓ Saved to {args.output_dir}")


if __name__ == "__main__":
    main()
