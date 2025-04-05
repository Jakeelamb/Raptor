#!/bin/bash

# Example workflow for yeast transcriptome assembly
# This script demonstrates a complete RNA-Seq assembly pipeline using Raptor
# Modified to use less memory (max 15GB)

set -e  # Exit on error

# Set a reasonable thread count for systems with limited memory
THREADS=4

echo "=== Raptor Yeast RNA-Seq Workflow ==="
echo "Limited to $THREADS threads to reduce memory usage"
echo "Downloading test data..."

# Download SRA yeast RNA-Seq data from S. cerevisiae
wget -c -q --show-progress ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR/SRR948/SRR948778/SRR948778_1.fastq.gz

echo "Running normalization..."
# Normalize reads to reduce redundancy while preserving coverage
RUSTFLAGS="-C target-cpu=native" RAYON_NUM_THREADS=$THREADS \
raptor normalize -i SRR948778_1.fastq.gz -o norm.fastq.gz --streaming --threads $THREADS

echo "Assembling transcripts..."
# Run assembly with isoform detection, TPM calculation, and GTF export
RUSTFLAGS="-C target-cpu=native" RAYON_NUM_THREADS=$THREADS \
raptor assemble -i norm.fastq.gz -o yeast --isoforms --compute-tpm \
  --threads $THREADS --streaming \
  --json-metadata yeast_meta.json --gtf yeast.gtf --gff3 yeast.gff3

echo "Generating PCA visualization..."
# Generate PCA plot for samples
raptor stats --input yeast_isoform.counts.matrix --pca yeast_pca.svg

echo "=== Workflow complete ==="
echo "Output files:"
echo "  yeast.fasta              - Assembled transcripts"
echo "  yeast_isoform.tpm.tsv    - Transcript expression values"
echo "  yeast_meta.json          - Transcript metrics in JSON format"
echo "  yeast.gtf                - GTF annotation file"
echo "  yeast.gff3               - GFF3 annotation file"
echo "  yeast_pca.svg            - PCA visualization" 