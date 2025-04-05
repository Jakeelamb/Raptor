# Repomix

A high-performance Rust-based RNA-seq assembler with advanced isoform detection and quantification.

## Features

- **K-mer normalization**: Efficient preprocessing with both GPU and CPU support
- **Greedy de novo assembly**: Fast contig generation from normalized reads
- **Isoform inference**: Sophisticated path traversal for transcript isoform detection
- **RLE compression**: Repeat-aware compression for efficient storage and processing
- **Quantification**: TPM calculation and expression matrix generation
- **Differential expression**: Statistical analysis of expression differences
- **Visualization**: PCA plots of expression data
- **Multiple outputs**: GFA1, GFA2, GFF3, GTF export formats
- **Hybrid assembly**: Long-read polishing and hybrid assembly capabilities
- **Benchmarking**: Evaluation against reference transcriptomes

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/repomix.git
cd repomix

# Build and install
cargo build --release
cargo install --path .
```

## Quick Start

```bash
# Normalize reads
repomix normalize -i reads.fastq.gz -o norm.fastq.gz

# Assemble transcripts
repomix assemble -i norm.fastq.gz -o assembled --isoforms

# Analyze differential expression
repomix diffexp --matrix expression.matrix --group-a sample1,sample2 --group-b sample3,sample4 --output diff_results.tsv
```

## Assembly Pipeline

1. **Normalization**: Reduce read redundancy while preserving coverage
   ```bash
   repomix normalize --input1 reads.fastq.gz --output normalized
   ```

2. **Assembly**: Reconstruct transcripts from normalized reads
   ```bash
   repomix assemble --input normalized_norm.fastq.gz --output assembled --isoforms --gtf transcripts.gtf
   ```

3. **Quantification**: Calculate expression levels
   ```bash
   # Uses the --compute-tpm flag during assembly
   repomix assemble --input normalized_norm.fastq.gz --output assembled --isoforms --compute-tpm
   ```

4. **Differential Expression**: Compare expression between conditions
   ```bash
   repomix diffexp --matrix assembled_isoform.counts.matrix --group-a sample1,sample2 --group-b sample3,sample4 --output diff_results.tsv
   ```

## Advanced Features

### Transcript Filtering by Expression

Filter out low-abundance transcripts:
```bash
repomix assemble --input reads.fastq.gz --output assembled --isoforms --compute-tpm --min-tpm 1.0
```

### Transcript Polishing with Long Reads

Use long-read alignments to improve transcript accuracy:
```bash
# First align long reads to transcripts with minimap2
minimap2 -ax map-ont assembled.fasta nanopore_reads.fastq > alignments.sam

# Then use the alignments to polish transcripts
repomix assemble --input reads.fastq.gz --output polished --isoforms --polish-reads alignments.sam
```

### Multi-Sample Analysis

Analyze multiple samples together:
```bash
# Create sample file (CSV format: sample_name,alignment_file)
echo "sample1,alignments1.sam" > samples.csv
echo "sample2,alignments2.sam" >> samples.csv

# Run assembly with multi-sample support
repomix assemble --input reads.fastq.gz --output multi --isoforms --compute-tpm --samples samples.csv
```

## GFA Path Traversal

This tool supports traversing paths defined in GFA files and reconstructing the full sequences. This is particularly useful for isoform analysis and visualization.

### Features

- **Edge-aware segment navigation**: Navigate through graph segments while respecting edges
- **Orientation-aware traversal**: Handle both forward (+) and reverse (-) segment orientation
- **Export options**: Output formats include FASTA, DOT (for Graphviz), TSV metadata, and more
- **Path visualization**: Generate visualizations for tools like Bandage or ODGI

### Usage

```bash
# Basic usage
cargo run --release -- traverse -i input.gfa -s segments.tsv -o output

# Advanced usage with all options
cargo run --release -- traverse \
    -i input.gfa \
    -s segments.tsv \
    -o output \
    --formats fasta,dot,json,odgi \
    --include-edges \
    --visualize \
    --metadata
```

### Input File Formats

- **GFA file**: Standard GFA format containing segments (S lines) and paths (P lines)
- **Segment sequences**: Tab-separated file with `segment_id\tsequence` format

### Output Files

- `output.fasta`: FASTA sequences for each path
- `output.dot`: DOT graph visualization (can be rendered using Graphviz)
- `output.json`: JSON metadata about paths
- `output.path_stats.tsv`: TSV file with path statistics
- `output.paths.gfa`: ODGI-compatible path definitions

### Example Script

For a complete example, see the `scripts/run_path_traversal.sh` script, which demonstrates:

1. Creating sample GFA and segment files
2. Running the traversal command
3. Visualizing the output paths

## Output Files

- `*.fasta`: Assembled sequences
- `*.gfa`: Graph Fragment Assembly format for visualizing the assembly graph
- `*.gtf`: Gene Transfer Format for genomic feature annotation
- `*.gff3`: General Feature Format version 3 for genomic features
- `*.tpm.tsv`: Transcript expression in Transcripts Per Million
- `*.counts.matrix`: Expression matrix for multiple samples
- `*_pca.svg`: PCA visualization of sample relationships

## Authors

- Your name and contributors

## License

This project is licensed under the MIT License - see the LICENSE file for details. 