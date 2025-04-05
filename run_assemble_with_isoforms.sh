#!/bin/bash

# This script demonstrates how to run the Raptor with isoform detection

# Set variables for input and output
INPUT_FASTQ="sample_large.fastq"
OUTPUT_PREFIX="assembled"
MIN_CONFIDENCE=0.25
THREADS=$(nproc)

# Check if input file exists
if [ ! -f "$INPUT_FASTQ" ]; then
    echo "Error: Input file $INPUT_FASTQ not found"
    echo "Generating sample data..."
    cargo run --bin generate_sample
    
    if [ ! -f "$INPUT_FASTQ" ]; then
        echo "Failed to generate sample data. Exiting."
        exit 1
    fi
    
    echo "Sample data generated successfully."
fi

# Function to check if compilation is needed
need_compile() {
    if [ ! -f "target/release/raptor" ]; then
        return 0  # Need to compile
    fi
    
    # Check if any source files are newer than the binary
    if find src -name "*.rs" -newer "target/release/raptor" | grep -q .; then
        return 0  # Need to compile
    fi
    
    return 1  # No need to compile
}

# Display header
echo "======================================================"
echo "Raptor RNA-Seq Assembly with Isoform Detection"
echo "======================================================"
echo "Input: $INPUT_FASTQ"
echo "Output: ${OUTPUT_PREFIX}_*.fasta"
echo "Threads: $THREADS"
echo "Min confidence: $MIN_CONFIDENCE"
echo "======================================================"

# Check if we need to compile
if need_compile; then
    echo "Run Raptor with isoform detection enabled"
    cargo build --release
fi

echo "======================================================"
# Run with isoform detection disabled (basic mode)
echo "Running Raptor in basic mode (without isoform detection):"
target/release/raptor assemble \
    -i "$INPUT_FASTQ" \
    -o "${OUTPUT_PREFIX}_basic" \
    --threads "$THREADS"

echo "======================================================"
# Run with isoform detection enabled
echo "Running Raptor with isoform detection:"
target/release/raptor assemble \
    -i "$INPUT_FASTQ" \
    -o "$OUTPUT_PREFIX" \
    --isoforms \
    --threads "$THREADS" \
    --min-confidence "$MIN_CONFIDENCE" \
    --gfa --gfa2 --gtf "${OUTPUT_PREFIX}.gtf" --gff3 "${OUTPUT_PREFIX}.gff3" \
    --json-metadata "${OUTPUT_PREFIX}_meta.json" --tsv-metadata "${OUTPUT_PREFIX}_meta.tsv"

# If needed, compile with custom features
if [ "$NEED_CUSTOM_COMPILE" = "true" ]; then
    echo "Compiling Raptor..."
    cargo build --release --features "gpu,avx2"
fi

# Run the Raptor with isoform detection
echo "Running Raptor..."

# Check if running succeeded
if [ $? -eq 0 ]; then
    echo "======================================================"
    echo "Assembly completed successfully!"
    echo "Output files:"
    echo "  - ${OUTPUT_PREFIX}.fasta        (Main assembly)"
    echo "  - ${OUTPUT_PREFIX}_isoforms.fasta   (Isoform sequences)"
    echo "  - ${OUTPUT_PREFIX}.gfa          (Assembly graph)"
    echo "  - ${OUTPUT_PREFIX}.gtf          (Gene annotations)"
    echo "  - ${OUTPUT_PREFIX}.gff3         (Gene annotations, GFF3 format)"
    echo "  - ${OUTPUT_PREFIX}_meta.json    (Metadata, JSON format)"
    echo "  - ${OUTPUT_PREFIX}_meta.tsv     (Metadata, TSV format)"
    echo "======================================================"
else
    echo "Assembly failed."
    exit 1
fi
