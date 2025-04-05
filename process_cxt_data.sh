#!/bin/bash

# Process Cxt RNA-Seq data through the Raptor pipeline
# Uses GPU-accelerated normalization and parallel assembly
# Memory-efficient version: limits threads and enables streaming mode

set -e

# Configuration
INPUT_R1="Cxt-r2-38_R1_001.fastq.gz"
INPUT_R2="Cxt-r2-38_R2_001.fastq.gz"
OUTPUT_DIR="assembly_output"
THREADS=2  # Reduced to 2 for even less memory usage
USE_GPU=true

# Check input files
if [ ! -f "$INPUT_R1" ] || [ ! -f "$INPUT_R2" ]; then
    echo "Error: Input files not found"
    echo "  R1: $INPUT_R1"
    echo "  R2: $INPUT_R2"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "======================================================"
echo "Raptor RNA-Seq Assembly Pipeline (Memory-Efficient Mode)"
echo "======================================================"
echo "Input R1: $INPUT_R1"
echo "Input R2: $INPUT_R2"
echo "Output directory: $OUTPUT_DIR"
echo "Threads: $THREADS (reduced for memory efficiency)"
echo "GPU-accelerated: $USE_GPU"
echo "======================================================"

# Step 1: Normalize reads (with streaming flag to reduce memory usage)
echo "Step 1: Normalizing paired-end reads (streaming mode)..."
if [ "$USE_GPU" = true ]; then
    echo "Using GPU acceleration"
    GPU_FLAG="--gpu"
else
    GPU_FLAG=""
fi

RUSTFLAGS="-C target-cpu=native" RAYON_NUM_THREADS=$THREADS \
cargo run --release --bin normalize_paired_reads -- \
    "$INPUT_R1" \
    "$INPUT_R2" \
    "$OUTPUT_DIR/normalized" \
    15 5 2 $GPU_FLAG --streaming

echo "Normalization complete. Output files:"
echo "  $OUTPUT_DIR/normalized_R1.norm.fastq.gz"
echo "  $OUTPUT_DIR/normalized_R2.norm.fastq.gz"
echo "======================================================"

# Step 2: Assemble normalized reads
echo "Step 2: Assembling normalized reads (reduced memory)..."

RUSTFLAGS="-C target-cpu=native" RAYON_NUM_THREADS=$THREADS \
cargo run --release --bin assemble -- \
    "$OUTPUT_DIR/normalized_R1.norm.fastq.gz" \
    "$OUTPUT_DIR/cxt_contigs.fasta"

echo "Assembly complete. Output files:"
echo "  $OUTPUT_DIR/cxt_contigs.fasta"
echo "======================================================"

# Step 3: Run the full Raptor pipeline (including isoform detection)
echo "Step 3: Running full Raptor pipeline with isoform detection..."

RUSTFLAGS="-C target-cpu=native" RAYON_NUM_THREADS=$THREADS \
cargo run --release -- assemble \
    -i "$OUTPUT_DIR/normalized_R1.norm.fastq.gz" \
    -o "$OUTPUT_DIR/cxt_assembly" \
    --threads "$THREADS" \
    --isoforms \
    --gfa \
    --gtf "$OUTPUT_DIR/cxt_transcripts.gtf" \
    --gff3 "$OUTPUT_DIR/cxt_transcripts.gff3" \
    --compute-tpm \
    --streaming

echo "Full pipeline complete. Output files:"
echo "  $OUTPUT_DIR/cxt_assembly.fasta (Main assembly)"
echo "  $OUTPUT_DIR/cxt_assembly_isoforms.fasta (Isoform sequences)"
echo "  $OUTPUT_DIR/cxt_assembly.gfa (Assembly graph)"
echo "  $OUTPUT_DIR/cxt_transcripts.gtf (Gene annotations)"
echo "  $OUTPUT_DIR/cxt_transcripts.gff3 (Gene annotations, GFF3 format)"
echo "======================================================"

echo "All processing complete!" 