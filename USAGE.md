# Raptor Command-Line Usage

This document provides a detailed guide to using the Raptor assembler and its various subcommands.

## Table of Contents

*   [Installation](#installation)
*   [Global Options](#global-options)
*   [Subcommands](#subcommands)
    *   [normalize](#normalize)
    *   [assemble](#assemble)
    *   [stats](#stats)
    *   [benchmark](#benchmark)
    *   [isoform](#isoform)
    *   [diffexp](#diffexp)
    *   [gtfcompare](#gtfcompare)
    *   [eval](#eval)
    *   [visualize](#visualize)
    *   [traverse](#traverse)
*   [Pipeline Examples](#pipeline-examples)
    *   [Basic Assembly](#basic-assembly)
    *   [Assembly with Isoform Detection & Quantification](#assembly-with-isoform-detection--quantification)
    *   [GPU-Accelerated Normalization & Assembly](#gpu-accelerated-normalization--assembly)

## Installation

Please refer to the main [README.md](./README.md) for installation instructions, including options for GPU and MPI support.

## Global Options

These options apply to the main `raptor` command:

*   `--version`: Display Raptor version information.
*   `--help`: Display the help message.

## Subcommands

Raptor operates through several subcommands, each targeting a specific stage of the assembly or analysis process.

---

### `normalize`

Normalize input reads using k-mer counting, with optional GPU acceleration. This step helps to remove sequencing errors and reduce data complexity before assembly.

**Arguments:**

| Argument                 | Short | Long                 | Type     | Description                                      | Default           |
| :----------------------- | :---- | :------------------- | :------- | :----------------------------------------------- | :---------------- |
| Input Read 1           | `-i`  | `--input1`           | `String` | Input read 1 FASTQ(.gz)                          | **Required**      |
| Input Read 2           | `-I`  | `--input2`           | `String` | Optional input read 2 for paired-end           | `None`            |
| Output Prefix            | `-o`  | `--output`           | `String` | Output prefix for normalized files             | **Required**      |
| GPU Acceleration         |       | `--gpu`              | `bool`   | Enable GPU-based k-mer counting                | `false`           |
| Threads                  |       | `--threads`          | `usize`  | Number of threads to use                       | System CPUs       |
| Streaming Mode           |       | `--streaming`        | `bool`   | Enable streaming mode for large datasets       | `false`           |
| Coverage Target          |       | `--coverage-target`  | `usize`  | Target coverage threshold for normalization    | `500`             |
| Max Reads                |       | `--max-reads`        | `usize`  | Maximum number of reads to process             | `5_000_000`       |

**Examples:**

*   Normalize single-end reads:
    ```bash
    raptor normalize -i reads.fastq.gz -o norm_reads --threads 8
    ```
*   Normalize paired-end reads using streaming and GPU:
    ```bash
    raptor normalize -i reads_R1.fastq.gz -I reads_R2


---

### `assemble`

Assemble normalized reads into contigs and optionally infer isoforms.

**Arguments:**

| Argument                   | Short | Long                 | Type     | Description                                                       | Default           |
| :------------------------- | :---- | :------------------- | :------- | :---------------------------------------------------------------- | :---------------- |
| Input FASTQ                | `-i`  | `--input`            | `String` | Input FASTQ(.gz) file (normalized)                                | **Required**      |
| Output Prefix              | `-o`  | `--output`           | `String` | Output prefix for FASTA/GFA/GTF etc. files                      | **Required**      |
| Minimum Contig Length      |       | `--min-len`          | `usize`  | Minimum contig length to output                                     | `50`              |
| Threads                    |       | `--threads`          | `usize`  | Number of threads                                                 | System CPUs       |
| Output GFA1                |       | `--gfa`              | `bool`   | Output GFA1 format as well                                        | `false`           |
| Output GFA2                |       | `--gfa2`             | `bool`   | Output GFA2 format as well                                        | `false`           |
| Adaptive K-mer             |       | `--adaptive-k`       | `bool`   | Enable adaptive k-mer selection                                   | `false`           |
| Run-Length Encoding        |       | `--rle`              | `bool`   | Enable run-length encoding for compression (affects assembly)     | `false`           |
| Collapse Repeats (RLE)     |       | `--collapse-repeats` | `bool`   | Enable repeat collapsing using RLE                                | `false`           |
| Min Repeat Length (RLE)    |       | `--min-repeat-len`   | `usize`  | Minimum repeat length to collapse (in RLE tuples)                 | `20`              |
| Contig Polishing           |       | `--polish`           | `bool`   | Enable contig polishing using input reads                         | `false`           |
| Polishing Window Size      |       | `--polish-window`    | `usize`  | Window size for polishing                                         | `25`              |
| Streaming Mode             |       | `--streaming`        | `bool`   | Enable streaming mode for large datasets (affects assembly)     | `false`           |
| Export Metadata            |       | `--export-metadata`  | `bool`   | Export metadata in JSON format (`<output>.contig_meta.json`)    | `false`           |
| JSON Metadata Path         |       | `--json-metadata`    | `String` | Optional path to write contig metadata as JSON                  | `None`            |
| TSV Metadata Path          |       | `--tsv-metadata`     | `String` | Optional path to write contig metadata as TSV                   | `None`            |
| Infer Isoforms             |       | `--isoforms`         | `bool`   | Enable isoform inference and transcript path export             | `false`           |
| GTF Output Path            |       | `--gtf`              | `String` | Optional path to write isoform GTF                              | `None`            |
| GFF3 Output Path           |       | `--gff3`             | `String` | Optional path to write isoform GFF3                             | `None`            |
| Max Path Depth (Isoform)   |       | `--max-path-depth`   | `usize`  | Maximum path depth for isoform traversal                          | `20`              |
| Min Confidence (Isoform)   |       | `--min-confidence`   | `f64`    | Minimum confidence for keeping isoform paths (0.0-1.0)            | `0.9`             |
| Min Path Length (Isoform)  |       | `--min-path-len`     | `usize`  | Minimum path length for isoform traversal                       | `50`              |
| Development Mode           |       | `--dev-mode`         | `bool`   | Bypasses filters and outputs debug info                         | `false`           |
| Compute TPM                |       | `--compute-tpm`      | `bool`   | Compute TPM expression values for transcripts                   | `false`           |
| Counts Matrix Output       |       | `--counts-matrix`    | `bool`   | Export transcript counts matrix (`<output>_isoform.counts.matrix`) | `false`           |
| Polish Isoforms            |       | `--polish-isoforms`  | `bool`   | Polish isoform sequences using aligned reads                    | `false`           |
| Samples CSV (TPM/Polish) |       | `--samples`          | `String` | CSV file with `sample_name,alignment_path` (SAM/BAM)           | `None`            |
| Min TPM (Isoform)          |       | `--min-tpm`          | `f64`    | Minimum TPM value for keeping transcripts                       | `0.1`             |
| Long Reads (Polish)        |       | `--polish-reads`     | `String` | SAM/BAM file with long reads mapped for polishing               | `None`            |
| Distributed Assembly       |       | `--distributed`      | `bool`   | Enable distributed assembly (Experimental)                      | `false`           |
| Buckets (Distributed)      |       | `--buckets`          | `usize`  | Number of buckets for distributed assembly                      | `16`              |

**Examples:**

*   Basic assembly from normalized reads:
    ```bash
    raptor assemble -i norm_reads.fastq.gz -o basic_assembly --threads 16
    ```
*   Assemble and output GFA graph:
    ```bash
    raptor assemble -i norm_reads.fastq.gz -o assembly_with_gfa --gfa
    ```
*   Assemble, infer isoforms, compute TPM, and output GTF:
    ```bash
    raptor assemble -i norm_reads.fastq.gz -o full_isoform_assembly \
      --isoforms --compute-tpm --gtf full_isoform_assembly.gtf \
      --counts-matrix --samples sample_alignments.csv --threads 16
    ```
*   Assemble with polishing using long reads:
    ```bash
    # Assume long reads are aligned to an initial assembly (e.g., basic_assembly.fasta)
    # minimap2 -ax map-ont basic_assembly.fasta long_reads.fastq.gz > long_read_alignments.sam

    raptor assemble -i norm_reads.fastq.gz -o polished_assembly \
      --isoforms --polish-isoforms --polish-reads long_read_alignments.sam
    ```

---

### `stats`

Calculate statistics for assembly output (FASTA or GFA).

**Arguments:**

| Argument       | Short | Long     | Type     | Description                            | Default |
| :------------- | :---- | :------- | :------- | :------------------------------------- | :------ |
| Input File     | `-i`  | `--input`  | `String` | Input file (FASTA or GFA)            | **Required** |
| Output Format  |       | `--format` | `String` | Output format (`json` or `tsv`)      | `json`  |
| Graph Stats    |       | `--graph`  | `bool`   | Enable graph stats (branchiness, etc.) | `false` |

**Examples:**

*   Calculate basic FASTA statistics (N50, lengths) and output as JSON:
    ```bash
    raptor stats -i assembly.fasta
    ```
*   Calculate statistics including graph metrics from a GFA file and output as TSV:
    ```bash
    raptor stats -i assembly.gfa --graph --format tsv
    ```

---

### `benchmark`

Benchmark k-mer counting performance.

**Arguments:**

| Argument   | Short | Long      | Type     | Description         | Default      |
| :--------- | :---- | :-------- | :------- | :------------------ | :----------- |
| Input FASTQ| `-i`  | `--input` | `String` | Input FASTQ file    | **Required** |
| K-mer Size |       | `-k`      | `usize`  | K-mer size to test  | `25`         |
| Threads    |       | `--threads` | `usize`  | Number of threads | System CPUs  |

**Example:**

```bash
raptor benchmark -i reads.fastq.gz -k 31 --threads 8
```

---

### `isoform`

Reconstruct isoforms from a GFA graph and expression data. (Note: Isoform reconstruction is typically done within the `assemble` command using the `--isoforms` flag. This subcommand might be for specialized workflows.)

**Arguments:**

| Argument               | Short | Long                   | Type     | Description                                                | Default        |
| :--------------------- | :---- | :--------------------- | :------- | :--------------------------------------------------------- | :------------- |
| Input GFA              | `-i`  | `--input`              | `String` | Input GFA file with contigs/overlaps                     | **Required**   |
| Expression Data (TSV)  | `-e`  | `--expression`         | `String` | TSV file (`contig_id`, `coverage`)                       | **Required**   |
| Output Prefix          | `-o`  | `--output`             | `String` | Output prefix for generated files                          | **Required**   |
| Minimum Confidence     |       | `--min-confidence`     | `f64`    | Minimum confidence score for transcript paths (0.0-1.0)  | `0.25`         |
| Maximum Depth          |       | `--max-depth`          | `usize`  | Maximum depth for graph traversal                          | `50`           |
| Output Formats         |       | `--formats`            | `String` | Comma-separated: `fasta`, `gfa`, `gtf`                   | `fasta,gfa`    |
| Threads                |       | `--threads`            | `usize`  | Number of threads to use                                   | System CPUs    |
| Output Stats           |       | `--stats`              | `bool`   | Output transcript statistics                               | `false`        |
| Filter Similar         |       | `--filter-similar`     | `bool`   | Enable similarity filtering to remove redundant transcripts | `false`        |
| Similarity Threshold   |       | `--similarity-threshold` | `f64`    | Similarity threshold for filtering (0.0-1.0)               | `0.8`          |
| Merge Similar          |       | `--merge-similar`      | `bool`   | Merge similar transcripts instead of filtering them          | `false`        |

**Example:**

```bash
raptor isoform -i assembly.gfa -e contig_coverage.tsv -o reconstructed_isoforms \
  --min-confidence 0.5 --threads 8 --formats fasta,gtf --stats
```

---

### `diffexp`

Perform differential expression analysis on a transcript counts matrix.

**Arguments:**

| Argument        | Short | Long          | Type     | Description                                             | Default     |
| :-------------- | :---- | :------------ | :------- | :------------------------------------------------------ | :---------- |
| Input Matrix    | `-m`  | `--matrix`    | `String` | Input counts matrix (`<output>_isoform.counts.matrix`) | **Required**|
| Group A Samples |       | `--group-a`   | `String` | Comma-separated list of sample names for group A        | **Required**|
| Group B Samples |       | `--group-b`   | `String` | Comma-separated list of sample names for group B        | **Required**|
| Output File     | `-o`  | `--output`    | `String` | Output file for differential expression results         | **Required**|
| P-value Threshold|      | `--p-value`   | `f64`    | P-value threshold for significance                    | `0.05`      |
| Fold Change Threshold| | `--fold-change`| `f64`    | Log2 fold-change threshold for significance           | `1.0`       |

**Example:**

```bash
raptor diffexp --matrix my_assembly_isoform.counts.matrix \
  --group-a sample1,sample2,sample3 \
  --group-b sample4,sample5,sample6 \
  -o diffexp_results.tsv --p-value 0.01
```

---

### `gtfcompare`

Compare predicted and truth GTF files to evaluate transcript accuracy.

**Arguments:**

| Argument        | Short | Long         | Type     | Description                              | Default     |
| :-------------- | :---- | :----------- | :------- | :--------------------------------------- | :---------- |
| Truth GTF       | `-t`  | `--truth`    | `String` | Truth/reference GTF file                 | **Required**|
| Predicted GTF   | `-p`  | `--predicted`| `String` | Predicted GTF file                       | **Required**|
| Output Metrics  | `-o`  | `--output`   | `String` | Output file for comparison metrics       | `None`      |

**Example:**

```bash
raptor gtfcompare --truth reference_annotation.gtf --predicted raptor_assembly.gtf -o comparison_metrics.txt
```

---

### `eval`

Evaluate assembly results against ground truth GTF. (Similar to `gtfcompare` but potentially with different metrics or focus).

**Arguments:**

| Argument        | Short | Long      | Type     | Description                           | Default     |
| :-------------- | :---- | :-------- | :------- | :------------------------------------ | :---------- |
| Truth GTF       |       | `--truth` | `String` | Truth/reference GTF file              | **Required**|
| Predicted GTF   |       | `--pred`  | `String` | Predicted GTF file                    | **Required**|
| Output Metrics  |       | `--output`| `String` | Output file for evaluation metrics    | `None`      |

**Example:**

```bash
raptor eval --truth reference_annotation.gtf --pred raptor_assembly.gtf --output eval_results.json
```

---

### `visualize`

Visualize TPM matrix with PCA plot and heatmap.

**Arguments:**

| Argument          | Short | Long         | Type     | Description                           | Default     |
| :---------------- | :---- | :----------- | :------- | :------------------------------------ | :---------- |
| Input TPM Matrix  |       | `--matrix`   | `String` | Path to TPM matrix file             | **Required**|
| Output PCA Plot   | `-o`  | `--output`   | `String` | Output PCA plot file (SVG or PNG)   | **Required**|
| Output Heatmap    |       | `--heatmap`  | `String` | Output heatmap file (PNG)           | `None`      |
| Output PCA File   |       | `--pca`      | `String` | Output PCA file (PNG)             | `None`      |
| PCA Components    |       | `--components` | `usize`  | Number of components for PCA        | `2`         |

**Example:**

```bash
raptor visualize --matrix my_assembly_isoform.counts.matrix \
  --output pca_plot.svg --heatmap heatmap.png --pca pca_data.png
```

---

### `traverse`

Traverse paths defined in a GFA file and export the corresponding sequences.

**Arguments:**

| Argument            | Short | Long           | Type     | Description                                       | Default     |
| :------------------ | :---- | :------------- | :------- | :------------------------------------------------ | :---------- |
| Input GFA           | `-i`  | `--input`      | `String` | Input GFA file with path definitions (P lines)    | **Required**|
| Segments TSV        | `-s`  | `--segments`   | `String` | TSV file: `segment_id\tsequence`                 | **Required**|
| Output Prefix       | `-o`  | `--output`     | `String` | Output file prefix                                | **Required**|
| Output Formats      |       | `--formats`    | `String` | Comma-separated: `fasta`, `dot`, `json`           | `fasta`     |
| Include Edges       |       | `--include-edges`| `bool`   | Include edge information in path (Not common)   | `false`     |
| Visualize (DOT)     |       | `--visualize`  | `bool`   | Generate DOT graph visualization                | `false`     |
| Export Metadata     |       | `--metadata`   | `bool`   | Export path metadata                              | `false`     |

**Example:**

*   Export path sequences to FASTA:
    ```bash
    raptor traverse -i assembly_with_paths.gfa -s assembly_segments.tsv -o traversed_paths
    ```
*   Export paths to FASTA and generate a DOT graph:
    ```bash
    raptor traverse -i assembly_with_paths.gfa -s assembly_segments.tsv -o traversed_paths_viz --formats fasta,dot --visualize
    ```

---

## Pipeline Examples

Here are examples combining multiple Raptor steps for common workflows.

### Basic Assembly

Normalize reads and assemble them into contigs.

```bash
#!/bin/bash

# Input reads
READS_R1="input_R1.fastq.gz"
READS_R2="input_R2.fastq.gz" # Optional: Leave empty for single-end

# Output prefix
PREFIX="my_basic_assembly"

# Parameters
THREADS=16

# 1. Normalize Reads
echo "Normalizing reads..."
if [ -z "$READS_R2" ]; then
  raptor normalize -i "$READS_R1" -o "${PREFIX}_norm" --threads "$THREADS"
  NORMALIZED_READS="${PREFIX}_norm.fastq.gz"
else
  raptor normalize -i "$READS_R1" -I "$READS_R2" -o "${PREFIX}_norm" --threads "$THREADS"
  NORMALIZED_READS="${PREFIX}_norm_R1.fastq.gz" # Assuming paired output naming convention
  # Check if normalization produced paired or single output based on actual implementation
fi

# 2. Assemble Contigs
echo "Assembling contigs..."
raptor assemble -i "$NORMALIZED_READS" -o "$PREFIX" --threads "$THREADS" --gfa

echo "Basic assembly finished. Output prefix: $PREFIX"
echo "Contigs: ${PREFIX}.fasta"
echo "Graph: ${PREFIX}.gfa"

```

### Assembly with Isoform Detection & Quantification

Normalize reads, assemble contigs, infer isoforms, compute TPM values, and generate annotation files.

```bash
#!/bin/bash

# Input reads
READS_R1="input_R1.fastq.gz"
READS_R2="input_R2.fastq.gz"

# Alignment files (one SAM/BAM per sample, listed in CSV)
SAMPLE_CSV="sample_alignments.csv"
# Format of sample_alignments.csv:
# sample1,path/to/sample1.bam
# sample2,path/to/sample2.bam
# ...

# Output prefix
PREFIX="my_isoform_assembly"

# Parameters
THREADS=16
MIN_CONFIDENCE=0.85
MIN_TPM=1.0

# 1. Normalize Reads
echo "Normalizing reads..."
raptor normalize -i "$READS_R1" -I "$READS_R2" -o "${PREFIX}_norm" --threads "$THREADS"
# Assuming paired output, adjust if needed
NORMALIZED_READS_R1="${PREFIX}_norm_R1.fastq.gz"
# Check normalization output for actual file names if needed for assembly input

# 2. Assemble, Infer Isoforms, and Quantify
echo "Assembling and processing isoforms..."
raptor assemble -i "$NORMALIZED_READS_R1" -o "$PREFIX" \
  --threads "$THREADS" \
  --isoforms \
  --compute-tpm \
  --counts-matrix \
  --samples "$SAMPLE_CSV" \
  --min-confidence "$MIN_CONFIDENCE" \
  --min-tpm "$MIN_TPM" \
  --gtf "${PREFIX}_isoforms.gtf" \
  --gff3 "${PREFIX}_isoforms.gff3" \
  --json-metadata "${PREFIX}_metadata.json"

echo "Isoform assembly finished. Output prefix: $PREFIX"
echo "Contigs: ${PREFIX}.fasta"
echo "Isoforms (FASTA): ${PREFIX}_isoforms.fasta"
echo "Isoforms (GTF): ${PREFIX}_isoforms.gtf"
echo "Isoforms (GFF3): ${PREFIX}_isoforms.gff3"
echo "TPM Matrix: ${PREFIX}_isoform.counts.matrix"
echo "Metadata: ${PREFIX}_metadata.json"

# 3. Optional: Visualize TPMs
echo "Generating PCA plot..."
raptor visualize --matrix "${PREFIX}_isoform.counts.matrix" --output "${PREFIX}_pca.svg"

```

### GPU-Accelerated Normalization & Assembly

Leverage GPU for the normalization step.

```bash
#!/bin/bash

# Ensure Raptor was compiled with GPU support!
# cargo build --release --features "gpu"

# Input reads
READS_R1="input_R1.fastq.gz"
READS_R2="input_R2.fastq.gz"

# Output prefix
PREFIX="my_gpu_assembly"

# Parameters
THREADS=16

# 1. Normalize Reads (with GPU)
echo "Normalizing reads using GPU..."
raptor normalize -i "$READS_R1" -I "$READS_R2" -o "${PREFIX}_norm" --gpu --threads "$THREADS" --streaming
NORMALIZED_READS_R1="${PREFIX}_norm_R1.fastq.gz"

# 2. Assemble Contigs
echo "Assembling contigs..."
raptor assemble -i "$NORMALIZED_READS_R1" -o "$PREFIX" --threads "$THREADS"

echo "GPU-accelerated assembly finished. Output prefix: $PREFIX"
```
