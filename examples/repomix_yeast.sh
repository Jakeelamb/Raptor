#!/bin/bash

# Example workflow for yeast transcriptome assembly
# This script demonstrates a complete RNA-Seq assembly pipeline using Repomix

set -e  # Exit on error

echo "=== Repomix Yeast RNA-Seq Workflow ==="
echo "Downloading test data..."

# Download SRA yeast RNA-Seq data from S. cerevisiae
wget -c -q --show-progress ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR/SRR948/SRR948778/SRR948778_1.fastq.gz

echo "Running normalization..."
# Normalize reads to reduce redundancy while preserving coverage
repomix normalize -i SRR948778_1.fastq.gz -o norm.fastq.gz --streaming

echo "Assembling transcripts..."
# Run assembly with isoform detection, TPM calculation, and GTF export
repomix assemble -i norm.fastq.gz -o yeast --isoforms --compute-tpm \
  --json-metadata yeast_meta.json --gtf yeast.gtf --gff3 yeast.gff3

echo "Generating PCA visualization..."
# Generate PCA plot for samples
repomix visualize --matrix yeast_isoform.counts.matrix --output yeast_pca.svg

echo "=== Workflow complete ==="
echo "Output files:"
echo "  yeast.fasta              - Assembled transcripts"
echo "  yeast_isoform.tpm.tsv    - Transcript expression values"
echo "  yeast_meta.json          - Transcript metrics in JSON format"
echo "  yeast.gtf                - GTF annotation file"
echo "  yeast.gff3               - GFF3 annotation file"
echo "  yeast_pca.svg            - PCA visualization" 