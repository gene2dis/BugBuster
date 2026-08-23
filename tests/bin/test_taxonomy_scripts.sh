#!/bin/bash
#
# Test suite for taxonomy Python scripts
# Tests taxonomy_report.py and taxonomy_phyloseq.py functionality
#

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BIN_DIR="${SCRIPT_DIR}/../../bin"

echo "=== Taxonomy Scripts Test Suite ==="
echo ""

# Check if scripts exist
echo "Checking script availability..."
if [ ! -f "${BIN_DIR}/taxonomy_report.py" ]; then
    echo "ERROR: taxonomy_report.py not found"
    exit 1
fi

if [ ! -f "${BIN_DIR}/taxonomy_phyloseq.py" ]; then
    echo "ERROR: taxonomy_phyloseq.py not found"
    exit 1
fi

if [ ! -f "${BIN_DIR}/tables_to_phyloseq.R" ]; then
    echo "ERROR: tables_to_phyloseq.R not found"
    exit 1
fi

echo "✓ All scripts found"
echo ""

# Check Python dependencies
echo "Checking Python dependencies..."
if python3 -c "import pandas, matplotlib, numpy, h5py" 2>/dev/null; then
    echo "✓ Python dependencies available"
else
    echo "WARNING: Some Python dependencies missing (pandas, matplotlib, numpy, h5py)"
    echo "Install with: pip install pandas matplotlib numpy h5py"
fi
echo ""

# Test taxonomy_report.py help
echo "Testing taxonomy_report.py --help..."
if python3 "${BIN_DIR}/taxonomy_report.py" --help > /dev/null 2>&1; then
    echo "✓ taxonomy_report.py help works"
else
    echo "ERROR: taxonomy_report.py help failed"
    exit 1
fi
echo ""

# Test taxonomy_phyloseq.py help
echo "Testing taxonomy_phyloseq.py --help..."
if python3 "${BIN_DIR}/taxonomy_phyloseq.py" --help > /dev/null 2>&1; then
    echo "✓ taxonomy_phyloseq.py help works"
else
    echo "ERROR: taxonomy_phyloseq.py help failed"
    exit 1
fi
echo ""

# Test tables_to_phyloseq.R help
echo "Testing tables_to_phyloseq.R --help..."
if Rscript "${BIN_DIR}/tables_to_phyloseq.R" --help > /dev/null 2>&1; then
    echo "✓ tables_to_phyloseq.R help works"
else
    echo "WARNING: tables_to_phyloseq.R help failed (R or phyloseq may not be installed)"
fi
echo ""

echo "=== Basic Tests Complete ==="
