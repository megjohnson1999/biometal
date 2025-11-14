"""
Norovirus Receptor Discovery Pipeline

Computational pipeline for discovering the cellular receptor for human norovirus
using AlphaFold3 structure prediction and systematic screening.
"""

__version__ = "1.0.0"
__author__ = "Your Name"

from pathlib import Path

# Project root
PROJECT_ROOT = Path(__file__).parent.parent

# Data directories
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
