#!/bin/bash
#SBATCH --job-name=raptor_assembly
#SBATCH --output=raptor_%j.out
#SBATCH --error=raptor_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

# Print job info
echo "Starting Raptor Pipeline on HPC"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURM_NODELIST"
echo "Cores allocated: $SLURM_CPUS_PER_TASK"
echo "Memory allocated: 128GB"
echo "GPU allocated: 1"

# Configuration
THREADS=$SLURM_CPUS_PER_TASK
K_SIZE=25
COVERAGE_TARGET=500
MIN_CONTIG_LEN=50
MIN_CONFIDENCE=0.9
MIN_PATH_LEN=50
RAPTOR_BIN="$PWD/target/release/raptor"
OUTPUT_DIR="raptor_output_${SLURM_JOB_ID}"

# Input files - replace with your actual paths
INPUT_R1="$1"  # First argument is Read 1
INPUT_R2="$2"  # Second argument is Read 2

# Check if input files exist
if [ ! -f "$INPUT_R1" ] || [ ! -f "$INPUT_R2" ]; then
    echo "Error: Input files not found!"
    echo "Usage: sbatch raptor_hpc.sh <read1.fastq.gz> <read2.fastq.gz>"
    exit 1
fi

echo "Input files:"
echo "  R1: $INPUT_R1"
echo "  R2: $INPUT_R2"

# Create output directory
mkdir -p $OUTPUT_DIR
echo "Output will be saved to: $OUTPUT_DIR"

# 1. Normalization step with GPU acceleration
echo "============================================================"
echo "[1/3] Running normalization with GPU acceleration"
echo "============================================================"
echo "Parameters:"
echo "  K-mer size: $K_SIZE"
echo "  Coverage target: $COVERAGE_TARGET"
echo "  Threads: $THREADS"

$RAPTOR_BIN normalize \
    -i $INPUT_R1 \
    -I $INPUT_R2 \
    -o $OUTPUT_DIR/normalized \
    --gpu \
    --threads $THREADS \
    --coverage-target $COVERAGE_TARGET

echo "Normalization complete."
echo "Output files:"
echo "  $OUTPUT_DIR/normalized_1.fq.gz"
echo "  $OUTPUT_DIR/normalized_2.fq.gz"

# 2. Assembly step
echo "============================================================"
echo "[2/3] Running assembly with optimized parameters"
echo "============================================================"
echo "Parameters:"
echo "  Min contig length: $MIN_CONTIG_LEN"
echo "  Threads: $THREADS"

$RAPTOR_BIN assemble \
    -i $OUTPUT_DIR/normalized_1.fq.gz \
    -o $OUTPUT_DIR/assembly \
    --min-len $MIN_CONTIG_LEN \
    --threads $THREADS \
    --gfa \
    --export-metadata \
    --json-metadata $OUTPUT_DIR/assembly_meta.json \
    --polish

echo "Assembly complete."
echo "Output files:"
echo "  $OUTPUT_DIR/assembly.fasta"
echo "  $OUTPUT_DIR/assembly.gfa"
echo "  $OUTPUT_DIR/assembly_meta.json"

# 3. Isoform detection and transcript assembly
echo "============================================================"
echo "[3/3] Running isoform detection and transcript assembly"
echo "============================================================"
echo "Parameters:"
echo "  Min confidence: $MIN_CONFIDENCE"
echo "  Min path length: $MIN_PATH_LEN"
echo "  Threads: $THREADS"

$RAPTOR_BIN assemble \
    -i $OUTPUT_DIR/normalized_1.fq.gz \
    -o $OUTPUT_DIR/transcripts \
    --threads $THREADS \
    --isoforms \
    --gfa \
    --compute-tpm \
    --gtf $OUTPUT_DIR/transcripts.gtf \
    --gff3 $OUTPUT_DIR/transcripts.gff3 \
    --min-len $MIN_CONTIG_LEN \
    --min-confidence $MIN_CONFIDENCE \
    --min-path-len $MIN_PATH_LEN \
    --export-metadata \
    --json-metadata $OUTPUT_DIR/transcripts_meta.json \
    --counts-matrix

echo "Isoform assembly complete."
echo "Output files:"
echo "  $OUTPUT_DIR/transcripts.fasta"
echo "  $OUTPUT_DIR/transcripts.gtf"
echo "  $OUTPUT_DIR/transcripts.gff3"
echo "  $OUTPUT_DIR/transcripts_meta.json"
echo "  $OUTPUT_DIR/transcripts.counts.matrix"

# 4. Visualize expression data
echo "============================================================"
echo "[4/4] Generating visualizations"
echo "============================================================"

$RAPTOR_BIN visualize \
    --matrix $OUTPUT_DIR/transcripts.counts.matrix \
    --output $OUTPUT_DIR/visualization.svg \
    --heatmap $OUTPUT_DIR/expression_heatmap.png \
    --pca $OUTPUT_DIR/pca_plot.png

echo "Visualization complete."
echo "Output files:"
echo "  $OUTPUT_DIR/visualization.svg"
echo "  $OUTPUT_DIR/expression_heatmap.png"
echo "  $OUTPUT_DIR/pca_plot.png"

echo "============================================================"
echo "Raptor pipeline completed successfully!"
echo "All results saved to: $OUTPUT_DIR"
echo "============================================================" 