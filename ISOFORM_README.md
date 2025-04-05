# Isoform Detection and Transcript Assembly

This feature adds Trinity-class RNA-Seq isoform assembly capabilities to the assembler. It enables the detection and reconstruction of transcript isoforms from assembled contigs using parallel graph traversal algorithms.

## Command Line Usage

To enable isoform detection during assembly, use the `--isoforms` flag:

```bash
cargo run -- assemble -i reads.fastq.gz -o output --isoforms
```

### Optional Parameters

You can customize the isoform detection process with these additional parameters:

- `--gtf <file>`: Specify a custom output path for the GTF file (default: output.isoforms.gtf)
- `--max-path-depth <value>`: Maximum depth for transcript path traversal (default: 20)
- `--min-confidence <value>`: Minimum confidence threshold for keeping transcript paths (default: 0.25)

Example with all options:

```bash
cargo run -- assemble -i reads.fastq.gz -o output --isoforms --gtf output.gtf --max-path-depth 30 --min-confidence 0.3
```

## Using the Helper Script

For convenience, a helper script `run_assemble_with_isoforms.sh` is provided to run the assembler with isoform detection:

```bash
./run_assemble_with_isoforms.sh -i reads.fastq.gz -o output --gtf output.gtf
```

The script accepts the following parameters:

- `-i, --input FILE`: Input FASTQ(.gz) file (default: reads.fastq.gz)
- `-o, --output PREFIX`: Output file prefix (default: output)
- `--gtf FILE`: GTF output file (default: output.gtf)
- `--min-confidence VALUE`: Minimum confidence for keeping isoform paths (default: 0.25)
- `--max-path-depth VALUE`: Maximum path depth for isoform traversal (default: 20)
- `-h, --help`: Show help message

## Expression Quantification

The assembler now supports several methods for quantifying transcript expression:

### TPM Quantification from Read Support

To compute TPM (Transcripts Per Million) values based on the reads used for assembly:

```bash
cargo run -- assemble -i reads.fastq.gz -o output --isoforms --compute-tpm
```

This will generate an additional file `output_isoform.tpm.tsv` with transcript IDs, lengths, and TPM values.

### Multi-Sample Expression Quantification

To quantify expression across multiple samples using SAM/BAM alignments:

1. Create a CSV file containing sample names and paths to SAM alignment files:
   ```
   sample1,path/to/sample1.sam
   sample2,path/to/sample2.sam
   sample3,path/to/sample3.sam
   ```

2. Run the assembler with both `--compute-tpm` and `--samples` options:
   ```bash
   cargo run -- assemble -i reads.fastq.gz -o output --isoforms --compute-tpm --samples samples.csv
   ```

This will generate:
- `output_isoform.tpm.tsv` - TPM values from the assembly reads
- `output_isoform.counts.matrix` - A Trinity-style expression matrix with TPM values for all samples

### Alignment-based Quantification

To generate transcript quantification using external aligners:

1. Assemble transcripts using the normal isoform pipeline
2. Align reads to the assembled transcripts using your favorite aligner (e.g., HISAT2, minimap2)
3. Use the resulting SAM/BAM files with the `--samples` option as described above

Example workflow with minimap2:
```bash
# First, assemble transcripts
cargo run -- assemble -i reads.fastq.gz -o transcripts --isoforms

# Align sample reads to the transcripts
minimap2 -ax sr transcripts.isoforms.fasta sample1.fastq > sample1.sam
minimap2 -ax sr transcripts.isoforms.fasta sample2.fastq > sample2.sam

# Create a samples.csv file
echo "sample1,sample1.sam" > samples.csv
echo "sample2,sample2.sam" >> samples.csv

# Run quantification
cargo run -- assemble -i reads.fastq.gz -o transcripts --isoforms --compute-tpm --samples samples.csv
```

## Additional Features

### Transcript Polishing

Improve transcript sequence quality by aligning and polishing with the input reads:

```bash
cargo run -- assemble -i reads.fastq.gz -o output --isoforms --polish-isoforms
```

### Transcript Statistics

When running with `--isoforms`, the assembler automatically provides transcript statistics including:
- Total number of transcripts
- Average transcript length
- N50 (weighted median transcript length)

## Output Files

When using the isoform features, the following files are generated:

| File | Description |
|------|-------------|
| `output.isoforms.fasta` | FASTA format transcript sequences |
| `output.isoforms.gtf` | GTF format transcript annotations |
| `output.isoforms.gfa` | GFA format transcript paths |
| `output.isoforms.json` | JSON format transcript details |
| `output.isoforms.stats.json` | Transcript statistics in JSON format |
| `output_isoform.tpm.tsv` | TPM expression values (when using `--compute-tpm`) |
| `output_isoform.counts.matrix` | Multi-sample expression matrix (when using `--samples`) |

## Technical Details

The isoform detection process works in several stages:

1. **Graph Building**: Constructs an isoform graph from the assembled contigs
2. **Parallel Path Discovery**: Finds potential transcript paths through the graph using multi-threaded traversal
3. **Confidence Scoring**: Evaluates paths based on k-mer coverage and graph connectivity
4. **Path Filtering**: Removes low-confidence paths that likely don't represent real transcripts
5. **Similarity Filtering**: Collapses redundant isoforms to reduce output complexity
6. **Transcript Assembly**: Generates final transcript sequences and annotations

All stages use parallel processing for optimal performance on multi-core systems.

## Differential Expression Analysis

The assembler now supports differential expression analysis between groups of samples:

### Running Differential Expression Analysis

After generating a counts matrix with multiple samples, you can identify differentially expressed transcripts:

```bash
cargo run -- diffexp --matrix output_isoform.counts.matrix \
  --group-a sample1,sample2 \
  --group-b sample3,sample4 \
  --output de_results.tsv
```

Optional parameters:
- `--p-value VALUE`: P-value threshold for significance (default: 0.05)
- `--fold-change VALUE`: Log2 fold-change threshold (default: 1.0)

### Output Format

The differential expression output file contains the following columns:
- `transcript_id`: Identifier of the transcript
- `log2FC`: Log2 fold change (a positive value means higher expression in group B)
- `p_value`: Statistical significance (smaller values indicate higher significance)

### Workflow Example

```bash
# 1. Assemble transcripts with isoform detection
cargo run -- assemble -i reads.fastq.gz -o transcripts --isoforms

# 2. Align reads from multiple samples
minimap2 -ax sr transcripts.isoforms.fasta sample1.fastq > sample1.sam
minimap2 -ax sr transcripts.isoforms.fasta sample2.fastq > sample2.sam
minimap2 -ax sr transcripts.isoforms.fasta sample3.fastq > sample3.sam
minimap2 -ax sr transcripts.isoforms.fasta sample4.fastq > sample4.sam

# 3. Create a samples file
echo "sample1,sample1.sam" > samples.csv
echo "sample2,sample2.sam" >> samples.csv
echo "sample3,sample3.sam" >> samples.csv
echo "sample4,sample4.sam" >> samples.csv

# 4. Generate the expression matrix
cargo run -- assemble -i reads.fastq.gz -o transcripts --isoforms \
  --compute-tpm --samples samples.csv

# 5. Perform differential expression analysis
cargo run -- diffexp --matrix transcripts_isoform.counts.matrix \
  --group-a sample1,sample2 --group-b sample3,sample4 \
  --output de_results.tsv
```

## GTF Comparison

The assembler now supports comparing predicted GTF files with reference/truth GTF files to evaluate assembly accuracy:

### Running GTF Comparison

```bash
cargo run -- gtf-compare --truth reference.gtf --predicted output.isoforms.gtf
```

Optional parameters:
- `--output PATH`: Write detailed metrics to a TSV file

### Output Metrics

The comparison outputs the following metrics:
- `True Positives`: Number of correctly predicted transcript features
- `False Positives`: Number of predicted features that don't match the reference
- `False Negatives`: Number of reference features that weren't predicted
- `Precision`: Proportion of predicted features that are correct
- `Recall`: Proportion of reference features that were predicted
- `F1 Score`: Harmonic mean of precision and recall

### Use Cases

This feature is particularly useful for:
- Evaluating assembler performance on simulated data where the ground truth is known
- Comparing the output of different assembly approaches
- Benchmarking against known transcript sets (e.g., Ensembl annotations)

## Summary Table of Features

| Feature | Description | Command |
|---------|-------------|---------|
| Isoform Assembly | Reconstruct transcript isoforms | `assemble --isoforms` |
| TPM Quantification | Calculate TPM expression values | `assemble --isoforms --compute-tpm` |
| Multi-sample Analysis | Process multiple sample alignments | `assemble --isoforms --samples samples.csv` |
| Differential Expression | Compare expression between sample groups | `diffexp --matrix --group-a --group-b` |
| GTF Comparison | Evaluate assembly accuracy | `gtf-compare --truth --predicted` |
| Transcript Polishing | Improve sequence quality | `assemble --isoforms --polish-isoforms` | 