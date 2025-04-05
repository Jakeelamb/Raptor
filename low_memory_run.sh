#!/bin/bash

# Ultra-low memory execution script for Raptor RNA-Seq pipeline
# Use this script when working on systems with limited memory (<16GB)

set -e

# Configuration
THREADS=2  # Absolute minimum threads for low memory usage
MAX_READS=5000000  # Maximum number of reads to process (subsampling)
COVERAGE_TARGET=500  # Higher coverage target for better normalization
MIN_CONTIG_LEN=50  # Lower minimum contig length for low-coverage data

# Print header
echo "============================================================"
echo "RAPTOR LOW-MEMORY EXECUTION MODE"
echo "============================================================"
echo "This script uses extreme memory optimization settings to run"
echo "Raptor on systems with less than 16GB of available RAM."
echo "Thread count: $THREADS"
echo "Max reads: $MAX_READS (subsampling enabled)"
echo "============================================================"

# Step 1: Compile with optimizations for memory efficiency
echo "[1/5] Compiling Raptor with memory optimizations..."

# Set memory-efficient compilation flags
export RUSTFLAGS="-C target-cpu=native -C opt-level=2 -C codegen-units=1"
export RAYON_NUM_THREADS=$THREADS

# Compile binaries
cargo build --release --bin normalize_paired_reads
cargo build --release --bin assemble

echo "Compilation complete."
echo "============================================================"

# Step 2: Check if input files exist
echo "[2/5] Checking input files..."

# For testing, use the smaller test files
# Comment these lines and uncomment the ones below for real data
TEST_MODE=false

if [ "$TEST_MODE" = true ]; then
    # Use test data for validation
    TEST_DATA_DIR="test_data"
    if [ ! -f "$TEST_DATA_DIR/test_R1.fastq.gz" ] || [ ! -f "$TEST_DATA_DIR/test_R2.fastq.gz" ]; then
        echo "Test data not found. Generating..."
        ./create_test_data.sh
    fi
    INPUT_R1="$TEST_DATA_DIR/test_R1.fastq.gz"
    INPUT_R2="$TEST_DATA_DIR/test_R2.fastq.gz"
    OUTPUT_DIR="low_mem_test_output"
else
    # Use real data for processing
    INPUT_R1="Cxt-r2-38_R1_001.fastq.gz"
    INPUT_R2="Cxt-r2-38_R2_001.fastq.gz"
    OUTPUT_DIR="low_mem_output"
fi

mkdir -p $OUTPUT_DIR

if [ ! -f "$INPUT_R1" ] || [ ! -f "$INPUT_R2" ]; then
    echo "Error: Input files not found"
    echo "  R1: $INPUT_R1"
    echo "  R2: $INPUT_R2"
    exit 1
fi

echo "Input files OK."
echo "Using: $INPUT_R1 and $INPUT_R2"
echo "Output will be saved to: $OUTPUT_DIR"
echo "============================================================"

# Step 3: Subsample the data if not in test mode
echo "[3/5] Subsampling input data..."

if [ "$TEST_MODE" = false ]; then
    SUBSAMPLED_R1="$OUTPUT_DIR/subsampled_R1.fastq.gz"
    SUBSAMPLED_R2="$OUTPUT_DIR/subsampled_R2.fastq.gz"
    
    echo "Subsampling $MAX_READS reads from input files..."
    zcat "$INPUT_R1" | head -n $(($MAX_READS * 4)) | gzip > "$SUBSAMPLED_R1"
    zcat "$INPUT_R2" | head -n $(($MAX_READS * 4)) | gzip > "$SUBSAMPLED_R2"
    
    # Update input paths to use subsampled data
    INPUT_R1="$SUBSAMPLED_R1"
    INPUT_R2="$SUBSAMPLED_R2"
    
    echo "Subsampling complete."
    echo "Subsampled files:"
    echo "  $SUBSAMPLED_R1"
    echo "  $SUBSAMPLED_R2"
else
    echo "Using test data, skipping subsampling."
fi
echo "============================================================"

# Step 4: Run normalization with memory-efficient settings
echo "[4/5] Running normalization with GPU acceleration..."

# Run normalizer with GPU acceleration for optimal memory usage
target/release/normalize_paired_reads \
    "$INPUT_R1" \
    "$INPUT_R2" \
    "$OUTPUT_DIR/normalized" \
    25 $COVERAGE_TARGET 1 --gpu

echo "Normalization complete."
echo "Output files:"
echo "  $OUTPUT_DIR/normalized_R1.norm.fastq.gz"
echo "  $OUTPUT_DIR/normalized_R2.norm.fastq.gz"
echo "============================================================"

# Step 5: Run assembly with minimal memory footprint
echo "[5/5] Running assembly with memory-efficient settings..."

target/release/assemble \
    "$OUTPUT_DIR/normalized_R1.norm.fastq.gz" \
    "$OUTPUT_DIR/assembly.fasta"

echo "Assembly complete."
echo "Output file: $OUTPUT_DIR/assembly.fasta"
echo "============================================================"

echo "Low-memory pipeline completed successfully!"
echo ""

if [ "$TEST_MODE" = true ]; then
    echo "Test mode: Completed successfully with test data."
    echo "To run with your real data, edit this script and set TEST_MODE=false"
else
    echo "Pipeline completed with subsampled data ($MAX_READS reads)."
    echo "If you need more comprehensive results, increase MAX_READS value."
    echo ""
    echo "For full pipeline with enhanced low-coverage processing:"
    echo ""
    echo "RUSTFLAGS=\"-C target-cpu=native\" RAYON_NUM_THREADS=$THREADS \\"
    echo "cargo run --release -- assemble \\"
    echo "    -i \"$OUTPUT_DIR/normalized_R1.norm.fastq.gz\" \\"
    echo "    -o \"$OUTPUT_DIR/assembly_full\" \\"
    echo "    --threads $THREADS \\"
    echo "    --isoforms \\"
    echo "    --gfa \\"
    echo "    --gtf \"$OUTPUT_DIR/transcripts.gtf\" \\"
    echo "    --min-len $MIN_CONTIG_LEN \\"
    echo "    --min-confidence 0.9 \\"
    echo "    --min-path-len 50 \\"
    echo "    --dev-mode"
fi

echo "============================================================"

# Make the script executable
chmod +x low_memory_run.sh 