"""
Candidate Library Builder for Norovirus Receptor Discovery

Queries UniProt, Human Protein Atlas, and other databases to build a library
of human intestinal epithelial membrane proteins with extracellular domains.
"""

import logging
import json
import time
from dataclasses import dataclass, asdict
from typing import List, Dict, Optional, Set
from pathlib import Path
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class CandidateProtein:
    """Represents a candidate receptor protein"""
    uniprot_id: str
    gene_name: str
    protein_name: str
    sequence: str
    length: int
    membrane_type: str  # 'single_pass', 'multi_pass', 'gpi_anchored'
    topology: Optional[str]  # Predicted topology
    ectodomain_regions: List[tuple]  # List of (start, end) tuples
    intestinal_expression: Optional[float]  # TPM or similar
    tissue_specificity: Optional[str]
    localization: Set[str]  # 'apical', 'basolateral', 'membrane'
    protein_class: Optional[str]  # 'receptor', 'transporter', 'adhesion', etc.
    go_terms: List[str]
    references: List[str]

    def to_dict(self):
        """Convert to dictionary with sets as lists"""
        d = asdict(self)
        d['localization'] = list(self.localization)
        return d


class UniProtQuerier:
    """Query UniProt for human membrane proteins"""

    BASE_URL = "https://rest.uniprot.org/uniprotkb"

    def __init__(self, rate_limit: float = 0.5):
        """
        Args:
            rate_limit: Seconds to wait between requests
        """
        self.rate_limit = rate_limit
        self.session = requests.Session()

    def query_intestinal_membrane_proteins(
        self,
        output_format: str = "json",
        max_results: int = 2000
    ) -> List[Dict]:
        """
        Query UniProt for human intestinal membrane proteins

        Query strategy:
        - Organism: Human (9606)
        - Location: Membrane
        - Tissue: Intestine/colon/duodenum/jejunum/ileum
        - Topology: Has transmembrane region or GPI anchor
        - Quality: Reviewed (Swiss-Prot)

        Returns:
            List of protein entry dictionaries
        """
        logger.info("Querying UniProt for intestinal membrane proteins...")

        # Build query string
        query_parts = [
            "(organism_id:9606)",  # Human
            "AND (reviewed:true)",  # Swiss-Prot only (high quality)
            "AND (cc_subcellular_location:membrane)",  # Membrane proteins
            "AND (cc_subcellular_location:\"Cell membrane\")",  # Cell membrane
            "AND (cc_tissue_specificity:intestin* OR cc_tissue_specificity:colon OR cc_tissue_specificity:duoden* OR cc_tissue_specificity:jejun* OR cc_tissue_specificity:ile*)",  # Intestinal
            "AND (ft_transmem:* OR ft_lipid:\"GPI-anchor\")",  # Has TM or GPI
        ]

        query = " ".join(query_parts)

        params = {
            "query": query,
            "format": output_format,
            "size": max_results,
            "fields": ",".join([
                "accession",
                "gene_names",
                "protein_name",
                "sequence",
                "length",
                "ft_transmem",
                "ft_topo_dom",
                "ft_lipid",
                "cc_subcellular_location",
                "cc_tissue_specificity",
                "go",
                "protein_families",
            ])
        }

        results = []
        url = f"{self.BASE_URL}/search"

        try:
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()

            data = response.json()
            results = data.get("results", [])

            logger.info(f"Retrieved {len(results)} proteins from UniProt")

            # Handle pagination if needed
            next_link = data.get("next")
            while next_link and len(results) < max_results:
                time.sleep(self.rate_limit)
                response = self.session.get(next_link, timeout=30)
                response.raise_for_status()
                data = response.json()
                results.extend(data.get("results", []))
                next_link = data.get("next")
                logger.info(f"Retrieved {len(results)} proteins so far...")

        except requests.exceptions.RequestException as e:
            logger.error(f"Error querying UniProt: {e}")

        return results

    def parse_topology(self, features: List[Dict]) -> tuple:
        """
        Parse transmembrane topology from UniProt features

        Returns:
            (membrane_type, topology_string, ectodomain_regions)
        """
        tm_regions = []
        topo_domains = []
        gpi_anchor = False

        for feature in features:
            ftype = feature.get("type")
            location = feature.get("location", {})
            start = location.get("start", {}).get("value")
            end = location.get("end", {}).get("value")
            description = feature.get("description", "")

            if ftype == "Transmembrane":
                if start and end:
                    tm_regions.append((start, end))
            elif ftype == "Topological domain":
                if start and end:
                    topo_domains.append({
                        "start": start,
                        "end": end,
                        "description": description
                    })
            elif ftype == "Lipidation" and "GPI-anchor" in description:
                gpi_anchor = True

        # Determine membrane type
        if gpi_anchor:
            membrane_type = "gpi_anchored"
        elif len(tm_regions) == 1:
            membrane_type = "single_pass"
        elif len(tm_regions) > 1:
            membrane_type = "multi_pass"
        else:
            membrane_type = "unknown"

        # Extract extracellular/ectodomain regions
        ectodomain_regions = []
        for domain in topo_domains:
            desc = domain["description"].lower()
            if "extracellular" in desc or "lumenal" in desc:
                ectodomain_regions.append((domain["start"], domain["end"]))

        # Build topology string
        topology = f"{len(tm_regions)}TM"
        if gpi_anchor:
            topology += "+GPI"

        return membrane_type, topology, ectodomain_regions


class HumanProteinAtlasQuerier:
    """Query Human Protein Atlas for expression data"""

    BASE_URL = "https://www.proteinatlas.org"

    def __init__(self, rate_limit: float = 1.0):
        self.rate_limit = rate_limit
        self.session = requests.Session()

    def get_tissue_expression(self, uniprot_id: str) -> Optional[Dict]:
        """
        Get tissue-specific expression data

        Returns:
            Dictionary with tissue expression levels
        """
        time.sleep(self.rate_limit)

        # HPA API endpoint (may need adjustment based on actual API)
        url = f"{self.BASE_URL}/{uniprot_id}.json"

        try:
            response = self.session.get(url, timeout=10)
            if response.status_code == 200:
                return response.json()
        except:
            pass

        return None


class CandidateLibraryBuilder:
    """Main builder class for candidate library"""

    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.uniprot = UniProtQuerier()
        self.hpa = HumanProteinAtlasQuerier()

    def build_library(self) -> List[CandidateProtein]:
        """
        Build complete candidate library

        Returns:
            List of CandidateProtein objects
        """
        logger.info("Starting candidate library build...")

        # Query UniProt
        uniprot_results = self.uniprot.query_intestinal_membrane_proteins()

        candidates = []
        for entry in uniprot_results:
            try:
                candidate = self._parse_uniprot_entry(entry)
                if candidate and self._filter_candidate(candidate):
                    candidates.append(candidate)
            except Exception as e:
                logger.warning(f"Error parsing entry {entry.get('primaryAccession')}: {e}")

        logger.info(f"Built library with {len(candidates)} candidates")

        # Save results
        self._save_library(candidates)

        return candidates

    def _parse_uniprot_entry(self, entry: Dict) -> Optional[CandidateProtein]:
        """Parse a UniProt entry into CandidateProtein"""

        # Basic info
        uniprot_id = entry.get("primaryAccession")
        gene_names = entry.get("genes", [])
        gene_name = gene_names[0].get("geneName", {}).get("value") if gene_names else "Unknown"

        protein_name_data = entry.get("proteinDescription", {})
        protein_name = protein_name_data.get("recommendedName", {}).get("fullName", {}).get("value", "Unknown")

        # Sequence
        sequence_data = entry.get("sequence", {})
        sequence = sequence_data.get("value", "")
        length = sequence_data.get("length", 0)

        # Parse topology
        features = entry.get("features", [])
        membrane_type, topology, ectodomain_regions = self.uniprot.parse_topology(features)

        # Subcellular location
        locations = entry.get("comments", [])
        localization = set()
        for comment in locations:
            if comment.get("commentType") == "SUBCELLULAR LOCATION":
                for loc in comment.get("subcellularLocations", []):
                    location_name = loc.get("location", {}).get("value", "").lower()
                    if "apical" in location_name:
                        localization.add("apical")
                    elif "basolateral" in location_name:
                        localization.add("basolateral")
                    elif "membrane" in location_name:
                        localization.add("membrane")

        # GO terms
        go_terms = []
        for go in entry.get("uniProtKBCrossReferences", []):
            if go.get("database") == "GO":
                go_terms.append(go.get("id"))

        # Protein class (simplified)
        protein_class = self._classify_protein(protein_name, go_terms)

        candidate = CandidateProtein(
            uniprot_id=uniprot_id,
            gene_name=gene_name,
            protein_name=protein_name,
            sequence=sequence,
            length=length,
            membrane_type=membrane_type,
            topology=topology,
            ectodomain_regions=ectodomain_regions,
            intestinal_expression=None,  # Would get from HPA
            tissue_specificity=None,
            localization=localization,
            protein_class=protein_class,
            go_terms=go_terms,
            references=[]
        )

        return candidate

    def _classify_protein(self, name: str, go_terms: List[str]) -> str:
        """Classify protein into functional category"""
        name_lower = name.lower()

        if "receptor" in name_lower:
            return "receptor"
        elif "integrin" in name_lower:
            return "integrin"
        elif "claudin" in name_lower or "occludin" in name_lower:
            return "tight_junction"
        elif "cadherin" in name_lower:
            return "adhesion"
        elif "transporter" in name_lower or "channel" in name_lower:
            return "transporter"
        elif "cd" in name_lower[:3]:  # CD antigens
            return "cd_antigen"
        elif "lectin" in name_lower or "selectin" in name_lower:
            return "lectin"
        else:
            return "other"

    def _filter_candidate(self, candidate: CandidateProtein) -> bool:
        """
        Filter candidates based on quality criteria

        Criteria:
        - Must have ectodomain > 50 amino acids
        - Must be membrane protein
        - Should be apical or general membrane (not strictly basolateral)
        """
        # Must have ectodomain
        if not candidate.ectodomain_regions:
            return False

        # Check ectodomain size
        total_ectodomain_length = sum(
            end - start + 1
            for start, end in candidate.ectodomain_regions
        )
        if total_ectodomain_length < 50:
            return False

        # Filter basolateral-only proteins
        if candidate.localization == {"basolateral"}:
            return False

        return True

    def _save_library(self, candidates: List[CandidateProtein]):
        """Save candidate library to files"""

        # Save as JSON
        json_file = self.output_dir / "candidate_library.json"
        with open(json_file, 'w') as f:
            json.dump(
                [c.to_dict() for c in candidates],
                f,
                indent=2
            )
        logger.info(f"Saved library to {json_file}")

        # Save as FASTA (full sequences)
        fasta_file = self.output_dir / "candidate_library.fasta"
        records = []
        for candidate in candidates:
            record = SeqRecord(
                Seq(candidate.sequence),
                id=candidate.uniprot_id,
                description=f"{candidate.gene_name}|{candidate.protein_name}|{candidate.topology}"
            )
            records.append(record)

        SeqIO.write(records, fasta_file, "fasta")
        logger.info(f"Saved sequences to {fasta_file}")

        # Save summary statistics
        stats = self._generate_statistics(candidates)
        stats_file = self.output_dir / "library_statistics.json"
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        logger.info(f"Saved statistics to {stats_file}")

    def _generate_statistics(self, candidates: List[CandidateProtein]) -> Dict:
        """Generate summary statistics"""

        stats = {
            "total_candidates": len(candidates),
            "by_membrane_type": {},
            "by_protein_class": {},
            "length_stats": {
                "min": min(c.length for c in candidates),
                "max": max(c.length for c in candidates),
                "mean": sum(c.length for c in candidates) / len(candidates)
            }
        }

        # Count by membrane type
        for candidate in candidates:
            mtype = candidate.membrane_type
            stats["by_membrane_type"][mtype] = stats["by_membrane_type"].get(mtype, 0) + 1

        # Count by protein class
        for candidate in candidates:
            pclass = candidate.protein_class
            stats["by_protein_class"][pclass] = stats["by_protein_class"].get(pclass, 0) + 1

        return stats


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Build candidate library for norovirus receptor discovery"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default="data/processed/candidates",
        help="Output directory for candidate library"
    )

    args = parser.parse_args()

    builder = CandidateLibraryBuilder(args.output_dir)
    candidates = builder.build_library()

    print(f"\n✓ Built library with {len(candidates)} candidates")
    print(f"✓ Saved to {args.output_dir}")


if __name__ == "__main__":
    main()
