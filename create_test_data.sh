#!/bin/bash

# Create valid small FASTQ test files for testing the pipeline
set -e

OUTPUT_DIR="test_data"
mkdir -p $OUTPUT_DIR

echo "Generating sample test data..."

# Use the sample generator that's already in the codebase
# First, compile it specifically
RUSTFLAGS="-C target-cpu=native" RAYON_NUM_THREADS=2 \
cargo build --release --bin generate_sample

# Run the binary directly with simpler parameters
target/release/generate_sample "$OUTPUT_DIR/test_R1.fastq"
# Duplicate it for paired-end testing
cp "$OUTPUT_DIR/test_R1.fastq" "$OUTPUT_DIR/test_R2.fastq"

# Compress them with gzip
gzip "$OUTPUT_DIR/test_R1.fastq"
gzip "$OUTPUT_DIR/test_R2.fastq"

if [ -f "$OUTPUT_DIR/test_R1.fastq.gz" ] && [ -f "$OUTPUT_DIR/test_R2.fastq.gz" ]; then
    echo "✓ Test data successfully generated"
    ls -lah $OUTPUT_DIR
else
    echo "✗ Failed to generate test data"
    exit 1
fi 