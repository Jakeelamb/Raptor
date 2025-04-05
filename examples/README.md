# Raptor Example Workflows

This directory contains example workflows for common RNA-Seq analysis scenarios.

## Yeast RNA-Seq Assembly

A complete workflow for assembling transcripts from _Saccharomyces cerevisiae_ RNA-Seq data.

```bash
bash examples/raptor_yeast.sh
```

### Outputs:

- `yeast.fasta` - Assembled transcripts in FASTA format
- `yeast_isoform.tpm.tsv` - Transcript expression values (TPM)
- `yeast_meta.json` - Detailed transcript metrics in JSON format
- `yeast.gtf` - Gene Transfer Format annotation file
- `yeast.gff3` - GFF3 annotation file for genome browsers
- `yeast_pca.svg` - Principal Component Analysis visualization

## Custom Workflow For Your Data

Create your own workflow script based on these examples:

```bash
#!/bin/bash

# Step 1: Normalize your reads
raptor normalize -i your_reads.fastq.gz -o norm.fastq.gz

# Step 2: Assemble transcripts with isoform detection
raptor assemble -i norm.fastq.gz -o output_prefix \
  --isoforms \
  --compute-tpm \
  --json-metadata metrics.json \
  --gtf transcripts.gtf
  
# Step 3: Filter transcripts by expression
raptor assemble -i norm.fastq.gz -o filtered \
  --isoforms \
  --compute-tpm \
  --min-tpm 1.0
  
# Step 4: Polish with long reads (if available)
minimap2 -ax map-ont output_prefix.fasta nanopore_reads.fastq > alignments.sam
raptor assemble -i norm.fastq.gz -o polished \
  --isoforms \
  --polish-reads alignments.sam
```

## Reproducibility

All examples are designed to be reproducible and deterministic. When reporting results, please include the version of Raptor used:

```bash
raptor --version
``` 