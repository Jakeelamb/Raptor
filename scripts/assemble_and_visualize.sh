#!/bin/bash
set -e

# RNA-Seq Assembly Pipeline
# Usage: ./assemble_and_visualize.sh read1.fastq.gz read2.fastq.gz threads use_gpu

# Check if required arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 read1.fastq.gz read2.fastq.gz [threads=4] [use_gpu=false]"
    exit 1
fi

# Parse arguments
READ1="$1"
READ2="$2"
THREADS="${3:-4}"
USE_GPU="${4:-false}"

echo "RNA-Seq Assembly Pipeline"
echo "==============================================="
echo "Input files:"
echo "  Read 1: $READ1"
echo "  Read 2: $READ2"
echo "Settings:"
echo "  Threads: $THREADS"
echo "  GPU acceleration: $USE_GPU"
echo ""

# Create output directory
OUTPUT_DIR="assembly_output"
mkdir -p "$OUTPUT_DIR"

# Step 1: Normalize reads
echo "[1/3] Normalizing reads..."
if [ "$USE_GPU" = "true" ] || [ "$USE_GPU" = "--gpu" ]; then
    cargo run --release -- normalize \
        --input1 "$READ1" \
        --input2 "$READ2" \
        --output "$OUTPUT_DIR/normalized" \
        --threads "$THREADS" \
        --gpu \
        --streaming
else
    cargo run --release -- normalize \
        --input1 "$READ1" \
        --input2 "$READ2" \
        --output "$OUTPUT_DIR/normalized" \
        --threads "$THREADS" \
        --streaming
fi

# Step 2: Assemble normalized reads
echo "[2/3] Assembling normalized reads..."
cargo run --release -- assemble \
    --input "$OUTPUT_DIR/normalized_R1.fastq.gz" \
    --paired-end "$OUTPUT_DIR/normalized_R2.fastq.gz" \
    --output "$OUTPUT_DIR/assembly" \
    --threads "$THREADS" \
    --use-gpu "$USE_GPU" \
    --k 25 \
    --min-count 2 \
    --output-gfa true \
    --adaptive-k true \
    --isoforms true \
    --gff3 true

# Step 3: Generate annotated GFA file for visualization
echo "[3/3] Generating annotated GFA file..."
cargo run --release -- traverse \
    --input "$OUTPUT_DIR/assembly.gfa" \
    --segments "$OUTPUT_DIR/assembly.segments.tsv" \
    --output "$OUTPUT_DIR/annotated" \
    --formats "fasta,json" \
    --include-edges true \
    --visualize true \
    --metadata true

echo ""
echo "Assembly complete! Output files are in the $OUTPUT_DIR directory."
echo "  - Normalized reads: $OUTPUT_DIR/normalized_R[1,2].fastq.gz"
echo "  - Assembly graph: $OUTPUT_DIR/assembly.gfa"
echo "  - Transcript sequences: $OUTPUT_DIR/annotated.fasta"
echo "  - Path metadata: $OUTPUT_DIR/annotated.json" 