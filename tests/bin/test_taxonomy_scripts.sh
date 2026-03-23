#!/bin/bash
#
# Test suite for taxonomy Python scripts
# Tests taxonomy_report.py and taxonomy_phyloseq.py functionality
#

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BIN_DIR="${SCRIPT_DIR}/../../bin"
TEST_DATA_DIR="${SCRIPT_DIR}/../data/taxonomy"

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
python3 -c "import pandas, matplotlib, numpy, h5py" 2>/dev/null
if [ $? -eq 0 ]; then
    echo "✓ Python dependencies available"
else
    echo "WARNING: Some Python dependencies missing (pandas, matplotlib, numpy, h5py)"
    echo "Install with: pip install pandas matplotlib numpy h5py"
fi
echo ""

# Test taxonomy_report.py help
echo "Testing taxonomy_report.py --help..."
python3 "${BIN_DIR}/taxonomy_report.py" --help > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "✓ taxonomy_report.py help works"
else
    echo "ERROR: taxonomy_report.py help failed"
    exit 1
fi
echo ""

# Test taxonomy_phyloseq.py help
echo "Testing taxonomy_phyloseq.py --help..."
python3 "${BIN_DIR}/taxonomy_phyloseq.py" --help > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "✓ taxonomy_phyloseq.py help works"
else
    echo "ERROR: taxonomy_phyloseq.py help failed"
    exit 1
fi
echo ""

# Test tables_to_phyloseq.R help
echo "Testing tables_to_phyloseq.R --help..."
Rscript "${BIN_DIR}/tables_to_phyloseq.R" --help > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "✓ tables_to_phyloseq.R help works"
else
    echo "WARNING: tables_to_phyloseq.R help failed (R or phyloseq may not be installed)"
fi
echo ""

echo "=== Basic Tests Complete ==="
echo ""
echo "To run integration tests with real data:"
echo "  1. Prepare test data in ${TEST_DATA_DIR}"
echo "  2. Run: bash tests/bin/test_taxonomy_integration.sh"
echo ""
