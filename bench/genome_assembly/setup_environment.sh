#!/bin/bash
# Set up the benchmark environment with required tools
# Usage: ./setup_environment.sh

set -e

echo "=== Setting up Genome Assembly Benchmark Environment ==="

# Check for conda
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found. Please install miniconda or anaconda first."
    exit 1
fi

ENV_NAME="raptor_bench"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_YML="${SCRIPT_DIR}/environment.yml"

if [ ! -f "${ENV_YML}" ]; then
    echo "ERROR: Missing environment file: ${ENV_YML}"
    exit 1
fi

if conda env list | grep -q "^${ENV_NAME} "; then
    echo "Updating existing environment ${ENV_NAME} from ${ENV_YML}"
    conda env update -n "${ENV_NAME}" -f "${ENV_YML}" --prune
else
    echo "Creating environment ${ENV_NAME} from ${ENV_YML}"
    conda env create -f "${ENV_YML}"
fi

eval "$(conda shell.bash hook)"
conda activate "${ENV_NAME}"

echo ""
echo "=== Setup Complete ==="
echo ""
echo "To use the benchmark environment, run:"
echo "  conda activate ${ENV_NAME}"
echo ""
echo "Then run benchmarks with:"
echo "  ./run_benchmark.sh simulated 8"
echo ""
echo "Environment file:"
echo "  ${ENV_YML}"
echo ""

# Show versions
echo "Installed versions:"
spades.py --version 2>&1 | head -1 || echo "SPAdes: not installed"
quast.py --version 2>&1 | head -1 || echo "QUAST: not installed"
