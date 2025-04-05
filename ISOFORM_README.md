# Isoform Detection in Raptor

This document provides detailed information about the isoform detection functionality in Raptor, including how it works, the available parameters, and the output formats generated.

## Overview

RNA-seq assembly aims to reconstruct transcripts from short reads. However, alternative splicing creates multiple transcript isoforms from the same gene, which presents computational challenges. Raptor implements a path traversal-based approach to detect and quantify these alternative isoforms.

## How Isoform Detection Works

Raptor performs isoform detection through the following steps:

1. **De Bruijn Graph Construction**: After normalization, reads are decomposed into k-mers and assembled into a de Bruijn graph.

2. **Graph Simplification**: The graph is simplified by removing tips, bubbles, and other artifacts.

3. **Contig Generation**: Initial contigs are generated using a greedy path traversal algorithm.

4. **Isoform Graph Construction**: Contigs are connected to form an isoform graph representing potential splice junctions.

5. **Path Traversal**: The graph is traversed to identify all possible paths, which represent potential isoforms.

6. **Filtering and Confidence Scoring**: Paths are filtered based on coverage, length, and confidence scores.

7. **Splice Junction Detection**: Splice junctions are identified by analyzing the patterns in the isoform graph.

8. **Strand Determination**: The strand of each transcript is determined when possible.

9. **Quantification**: Expression levels are calculated for each isoform using the read coverage information.

### Splicing Detection

Raptor can detect several types of alternative splicing events:

- **Exon Skipping**: An exon that appears in some isoforms is skipped in others.
- **Alternative 5' or 3' Splice Sites**: Different exon boundaries are used in different isoforms.
- **Intron Retention**: An intron that is spliced out in some isoforms is retained in others.
- **Mutually Exclusive Exons**: Different exons are used in different isoforms.

## Command-Line Usage

### Basic Usage

```bash
raptor assemble --input normalized_reads.fastq.gz --output results --isoforms
```

### Complete Usage with All Isoform Options

```bash
raptor assemble \
  --input normalized_reads.fastq.gz \
  --output results \
  --threads 12 \
  --isoforms \
  --compute-tpm \
  --min-confidence 0.3 \
  --max-path-depth 50 \
  --min-isoform-length 200 \
  --polish-isoforms \
  --polish-reads long_reads.sam \
  --counts-matrix \
  --gtf transcripts.gtf \
  --gff3 transcripts.gff3
```

## Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--input` | Input normalized FASTQ file |
| `--output` | Prefix for output files |
| `--isoforms` | Enable isoform detection (flag) |

### Optional Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--threads` | Number of CPU threads to use | 4 |
| `--compute-tpm` | Calculate TPM values (flag) | false |
| `--min-confidence` | Minimum confidence score to keep an isoform (0.0-1.0) | 0.2 |
| `--max-path-depth` | Maximum depth for path traversal in the graph | 100 |
| `--min-isoform-length` | Minimum length (bp) to keep an isoform | 150 |
| `--polish-isoforms` | Enable isoform polishing (flag) | false |
| `--polish-reads` | Path to long reads SAM/BAM file for polishing | - |
| `--counts-matrix` | Export transcript counts matrix (flag) | false |
| `--gtf` | Output GTF file path | - |
| `--gff3` | Output GFF3 file path | - |
| `--samples` | CSV file with sample names and alignment files | - |
| `--strand-specific` | Input data is strand-specific (flag) | false |
| `--min-tpm` | Minimum TPM to keep an isoform | 0.1 |
| `--cluster-similar` | Cluster similar isoforms (flag) | false |
| `--similarity-threshold` | Similarity threshold for clustering (0.0-1.0) | 0.9 |

## Output Files

When running with the `--isoforms` flag, Raptor produces the following output files:

| File | Description |
|------|-------------|
| `[prefix].isoforms.fasta` | FASTA file containing the sequences of all detected isoforms |
| `[prefix].isoforms.gfa` | GFA format representation of the isoform graph |
| `[prefix].isoform.tpm.tsv` | TSV file with TPM values for each isoform (if `--compute-tpm` is enabled) |
| `[prefix]_isoform.counts.matrix` | Matrix of isoform counts across samples (if `--counts-matrix` is enabled) |
| `[prefix].isoforms.gtf` | GTF annotation file for the isoforms (if `--gtf` path is provided) |
| `[prefix].isoforms.gff3` | GFF3 annotation file for the isoforms (if `--gff3` path is provided) |
| `[prefix].splice_junctions.bed` | BED file containing all detected splice junctions |
| `[prefix].isoform_clusters.tsv` | TSV file with clustering information (if `--cluster-similar` is enabled) |

## Output Format Details

### FASTA Output

The FASTA output contains the sequence of each isoform with headers in the following format:

```
>transcript_1 length=1245 confidence=0.87 tpm=12.3 splicing=SE,IR strand=+
ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG...
```

The header includes:
- Transcript ID
- Length in base pairs
- Confidence score (0.0-1.0)
- TPM value (if `--compute-tpm` is enabled)
- Splicing types detected (SE: Skipped Exon, IR: Intron Retention, etc.)
- Strand information (+ or -)

### GTF/GFF3 Output

The GTF/GFF3 output follows standard format requirements and includes the following features:

- Transcript records
- Exon records with parent-child relationships
- CDS predictions (when possible)
- Splice junction annotations
- Expression information in the attributes field

Example GTF line:
```
scaffold_1  raptor  transcript  1245  2468  .  +  .  transcript_id "transcript_1"; gene_id "gene_1"; confidence "0.87"; TPM "12.3";
scaffold_1  raptor  exon        1245  1500  .  +  .  transcript_id "transcript_1"; gene_id "gene_1"; exon_number "1";
```

### TPM Values File

The TPM (Transcripts Per Million) file is a tab-separated file with the following columns:

1. Transcript ID
2. Length (bp)
3. TPM value

### Count Matrix

The count matrix file format has transcripts as rows and samples as columns. The first line is a header with sample names:

```
transcript_id sample1 sample2 sample3 ...
transcript_1  124     98      156     ...
transcript_2  56      43      87      ...
...
```

### Splice Junction BED

The splice junction BED file uses the BED12 format to represent splice junctions:

```
scaffold_1  1244  2468  transcript_1  87  +  1245  2467  0,0,0  2  255,223  0,1000
```

This format allows visualization of splice junctions in genome browsers.

## Advanced Features

### Transcript Clustering

When `--cluster-similar` is enabled, Raptor clusters similar isoforms to reduce redundancy. The similarity threshold controls how similar transcripts must be to be clustered together. The clustering information is output to `[prefix].isoform_clusters.tsv`.

### Transcript Polishing

When `--polish-isoforms` is enabled, Raptor uses information from aligned long reads (provided via `--polish-reads`) to improve the accuracy of isoform sequences. This can correct errors in the assembly and improve the detection of splice junctions.

### Multi-Sample TPM Calculation

When a samples file is provided via `--samples`, Raptor calculates TPM values for each sample. The samples file format is a CSV with two columns:

```
sample_name,alignment_file
```

Where `alignment_file` is a SAM/BAM file with reads aligned to the assembled transcripts.

## Example Workflows

### Basic Isoform Detection

```bash
# Normalize reads
raptor normalize --input1 reads_1.fastq.gz --input2 reads_2.fastq.gz --output normalized

# Detect isoforms
raptor assemble --input normalized_norm.fastq.gz --output results --isoforms
```

### Isoform Detection with Expression Quantification

```bash
# Normalize reads
raptor normalize --input1 reads_1.fastq.gz --input2 reads_2.fastq.gz --output normalized

# Detect isoforms and calculate TPM
raptor assemble --input normalized_norm.fastq.gz --output results --isoforms --compute-tpm
```

### Isoform Detection with Long-Read Polishing

```bash
# Normalize short reads
raptor normalize --input1 reads_1.fastq.gz --input2 reads_2.fastq.gz --output normalized

# Align long reads to assembled transcripts
# (first run a basic assembly)
raptor assemble --input normalized_norm.fastq.gz --output initial

# Align long reads to the initial assembly
minimap2 -ax map-ont initial.isoforms.fasta long_reads.fastq > long_alignments.sam

# Run full assembly with polishing
raptor assemble --input normalized_norm.fastq.gz --output final \
  --isoforms --compute-tpm --polish-isoforms --polish-reads long_alignments.sam
```

### Multi-Sample Analysis

```bash
# Normalize reads
raptor normalize --input1 reads_1.fastq.gz --input2 reads_2.fastq.gz --output normalized

# Create sample information file
echo "sample1,sample1.sam" > samples.csv
echo "sample2,sample2.sam" >> samples.csv

# Detect isoforms with multi-sample quantification
raptor assemble --input normalized_norm.fastq.gz --output results \
  --isoforms --compute-tpm --samples samples.csv --counts-matrix
```

## Troubleshooting

### Common Issues

1. **Too many or too few isoforms detected**
   - Adjust `--min-confidence` (lower for more isoforms, higher for fewer)
   - Adjust `--max-path-depth` (higher for more complex graphs)
   - Verify input read quality and coverage

2. **Missing splice junctions**
   - Check that read coverage is sufficient across potential junction sites
   - Consider using longer reads if available

3. **Poor quantification accuracy**
   - Ensure sufficient read depth for reliable quantification
   - Use biological replicates when possible
   - Consider validating key isoforms with targeted methods (RT-PCR, etc.)

## Performance Considerations

- Isoform detection is computationally intensive; use more threads when available
- Memory usage increases with graph complexity; for very large datasets, consider:
  - Running on a high-memory machine
  - Splitting the analysis by chromosome or scaffold
  - Using a higher confidence threshold to reduce the number of paths

## References

For more details on the algorithmic approaches used in Raptor's isoform detection, consult the following papers:

1. Grabherr MG, et al. (2011). "Full-length transcriptome assembly from RNA-Seq data without a reference genome."
2. Pertea M, et al. (2015). "StringTie enables improved reconstruction of a transcriptome from RNA-seq reads."
3. Li B & Dewey CN (2011). "RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome."

## Citation

If you use Raptor's isoform detection in your research, please cite:

```
[Citation information will be added upon publication]
``` 