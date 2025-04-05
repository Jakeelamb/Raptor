#!/bin/bash

# This script demonstrates how to run the assembler with isoform detection
# and custom GTF output

# Make script exit on any error
set -e

# Default values
INPUT_FILE="reads.fastq.gz"
OUTPUT_PREFIX="output"
GTF_FILE="output.gtf"
MIN_CONFIDENCE=0.25
MAX_PATH_DEPTH=20

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--input)
      INPUT_FILE="$2"
      shift 2
      ;;
    -o|--output)
      OUTPUT_PREFIX="$2"
      shift 2
      ;;
    --gtf)
      GTF_FILE="$2"
      shift 2
      ;;
    --min-confidence)
      MIN_CONFIDENCE="$2"
      shift 2
      ;;
    --max-path-depth)
      MAX_PATH_DEPTH="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: $0 [-i INPUT_FILE] [-o OUTPUT_PREFIX] [--gtf GTF_FILE] [--min-confidence CONFIDENCE] [--max-path-depth DEPTH]"
      echo "Run assembler with isoform detection enabled"
      echo ""
      echo "Options:"
      echo "  -i, --input FILE          Input FASTQ(.gz) file (default: reads.fastq.gz)"
      echo "  -o, --output PREFIX       Output file prefix (default: output)"
      echo "  --gtf FILE                GTF output file (default: output.gtf)"
      echo "  --min-confidence VALUE    Minimum confidence for keeping isoform paths (default: 0.25)"
      echo "  --max-path-depth VALUE    Maximum path depth for isoform traversal (default: 20)"
      echo "  -h, --help                Show this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
  echo "Error: Input file $INPUT_FILE does not exist"
  exit 1
fi

echo "Running assembler with isoform detection:"
echo "  Input: $INPUT_FILE"
echo "  Output prefix: $OUTPUT_PREFIX"
echo "  GTF file: $GTF_FILE"
echo "  Min confidence: $MIN_CONFIDENCE"
echo "  Max path depth: $MAX_PATH_DEPTH"
echo ""

# Compile if needed (comment out if already compiled)
echo "Compiling assembler..."
cargo build --release

# Run the assembler with isoform detection
echo "Running assembler..."
cargo run --release -- assemble \
  -i "$INPUT_FILE" \
  -o "$OUTPUT_PREFIX" \
  --isoforms \
  --gtf "$GTF_FILE" \
  --min-confidence "$MIN_CONFIDENCE" \
  --max-path-depth "$MAX_PATH_DEPTH"

# Check results
if [ -f "$GTF_FILE" ]; then
  echo "Success! GTF file created: $GTF_FILE"
  echo "Transcript count: $(grep -c 'transcript' "$GTF_FILE")"
else
  echo "Error: GTF file was not created"
  exit 1
fi

echo "Isoform FASTA file: ${OUTPUT_PREFIX}.isoforms.fasta"
echo "Isoform GFA file: ${OUTPUT_PREFIX}.isoforms.gfa"
echo "Isoform JSON file: ${OUTPUT_PREFIX}.isoforms.json"
echo "Statistics file: ${OUTPUT_PREFIX}.isoforms.stats.json"

echo "Done."
