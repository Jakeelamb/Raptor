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

# Create conda environment for benchmarking
ENV_NAME="raptor_bench"

if conda env list | grep -q "^${ENV_NAME} "; then
    echo "Environment ${ENV_NAME} already exists"
    echo "Activating..."
else
    echo "Creating conda environment: ${ENV_NAME}"
    conda create -n "${ENV_NAME}" -y python=3.10
fi

# Activate and install tools
echo "Installing benchmark tools..."
eval "$(conda shell.bash hook)"
conda activate "${ENV_NAME}"

# Install SPAdes
if ! command -v spades.py &> /dev/null; then
    echo "Installing SPAdes..."
    conda install -c bioconda spades -y
fi

# Install QUAST for evaluation
if ! command -v quast.py &> /dev/null; then
    echo "Installing QUAST..."
    conda install -c bioconda quast -y
fi

# Install other useful tools
conda install -c bioconda seqkit samtools -y 2>/dev/null || true

echo ""
echo "=== Setup Complete ==="
echo ""
echo "To use the benchmark environment, run:"
echo "  conda activate ${ENV_NAME}"
echo ""
echo "Then run benchmarks with:"
echo "  ./run_benchmark.sh simulated 8"
echo ""

# Show versions
echo "Installed versions:"
spades.py --version 2>&1 | head -1 || echo "SPAdes: not installed"
quast.py --version 2>&1 | head -1 || echo "QUAST: not installed"
