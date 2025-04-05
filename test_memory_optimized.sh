#!/bin/bash

# Test script for memory-optimized Raptor pipeline
# Uses small test files to validate functionality

set -e

# Configuration
THREADS=2
OUTPUT_DIR="test_output"
mkdir -p $OUTPUT_DIR

# First generate test data if it doesn't exist
TEST_DATA_DIR="test_data"
if [ ! -f "$TEST_DATA_DIR/test_R1.fastq.gz" ] || [ ! -f "$TEST_DATA_DIR/test_R2.fastq.gz" ]; then
    echo "Test data not found. Generating..."
    ./create_test_data.sh
fi

echo "==== Testing memory-optimized normalize_paired_reads ===="

# Use memory-efficient settings for testing
RUSTFLAGS="-C target-cpu=native" RAYON_NUM_THREADS=$THREADS \
cargo run --release --bin normalize_paired_reads -- \
    "$TEST_DATA_DIR/test_R1.fastq.gz" \
    "$TEST_DATA_DIR/test_R2.fastq.gz" \
    "$OUTPUT_DIR/test" \
    15 3 1 --streaming

if [ -f "$OUTPUT_DIR/test_R1.norm.fastq.gz" ] && [ -f "$OUTPUT_DIR/test_R2.norm.fastq.gz" ]; then
    echo "✓ Normalization successful"
else
    echo "✗ Normalization failed"
    exit 1
fi

echo "==== Testing memory-optimized assembly ===="

# Test assembly with the normalized reads
RUSTFLAGS="-C target-cpu=native" RAYON_NUM_THREADS=$THREADS \
cargo run --release --bin assemble -- \
    "$OUTPUT_DIR/test_R1.norm.fastq.gz" \
    "$OUTPUT_DIR/test_assembly.fasta"

if [ -f "$OUTPUT_DIR/test_assembly.fasta" ]; then
    echo "✓ Assembly successful"
else
    echo "✗ Assembly failed"
    exit 1
fi

echo "==== Checking file sizes ===="
ls -lah $OUTPUT_DIR

echo "==== All tests passed ===="
echo "The memory-optimized pipeline is functioning correctly." 