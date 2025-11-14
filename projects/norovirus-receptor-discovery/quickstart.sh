#!/bin/bash
# Quick start script for Norovirus Receptor Discovery Pipeline

set -e

echo "========================================"
echo "Norovirus Receptor Discovery Pipeline"
echo "Quick Start Setup"
echo "========================================"
echo

# Check Python version
echo "[1/5] Checking Python version..."
python_version=$(python3 --version 2>&1 | awk '{print $2}' | cut -d. -f1,2)
echo "Found Python $python_version"

if [ $(echo "$python_version < 3.9" | bc -l) -eq 1 ]; then
    echo "ERROR: Python 3.9+ required"
    exit 1
fi

# Create virtual environment
echo
echo "[2/5] Creating virtual environment..."
if [ ! -d "venv" ]; then
    python3 -m venv venv
    echo "✓ Virtual environment created"
else
    echo "✓ Virtual environment already exists"
fi

# Activate virtual environment
echo
echo "[3/5] Activating virtual environment..."
source venv/bin/activate
echo "✓ Virtual environment activated"

# Install dependencies
echo
echo "[4/5] Installing dependencies..."
pip install --upgrade pip > /dev/null 2>&1
pip install -r requirements.txt
echo "✓ Dependencies installed"

# Create directory structure
echo
echo "[5/5] Setting up directory structure..."
mkdir -p data/{raw,processed/{candidates,structures},results}
mkdir -p results/{predictions,plots}
mkdir -p logs
echo "✓ Directories created"

echo
echo "========================================"
echo "Setup Complete!"
echo "========================================"
echo
echo "Next steps:"
echo
echo "1. Activate the virtual environment:"
echo "   source venv/bin/activate"
echo
echo "2. Configure AlphaFold3 backend in config/default.json"
echo
echo "3. Run the pipeline:"
echo "   python pipeline.py --full"
echo
echo "4. Or run individual stages:"
echo "   python pipeline.py --stage candidates"
echo "   python pipeline.py --stage structures"
echo "   python pipeline.py --stage predictions"
echo "   python pipeline.py --stage analysis"
echo "   python pipeline.py --stage scoring"
echo
echo "5. Analyze results:"
echo "   jupyter notebook notebooks/01_analyze_results.ipynb"
echo
echo "For more information, see README.md"
echo
