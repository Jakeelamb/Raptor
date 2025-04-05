# Benchmarking Setup

This directory contains tools and datasets for benchmarking the RNA-Seq assembler against reference transcripts.

## Required Files

- `truth.gtf`: Reference GTF with known transcript annotations
- `reads.fastq.gz`: Simulated or real reads matching the reference
- `assembled.gtf`: Assembler output (generated with `--isoforms --gtf` options)

## Running the Evaluation

To compare predicted transcripts against the reference transcripts:

```bash
cargo run -- gtf-compare --truth truth.gtf --pred assembled.gtf
```

Optional parameters:
- `--output <file>`: Write detailed metrics to a TSV file

## Sample Benchmark Datasets

### Simulated Data

The `simulate_data.py` script can generate a synthetic dataset with known transcripts:

```bash
python simulate_data.py --transcripts 500 --reads 100000 --output benchmark_data
```

This will create:
- `benchmark_data/truth.gtf` - Reference transcript annotations
- `benchmark_data/truth.fasta` - Reference transcript sequences
- `benchmark_data/reads.fastq.gz` - Simulated RNA-Seq reads

### Real Data with References

We also provide benchmark datasets from:

1. **Human (Ensembl):** Selected chromosomes with reference transcripts
2. **Yeast:** Complete transcriptome with reference annotations
3. **Drosophila:** Selected genes with alternative splicing

## Metrics

The evaluation computes standard transcript recovery metrics:

- **Precision:** Proportion of predicted transcripts that match reference
- **Recall:** Proportion of reference transcripts that were predicted
- **F1 Score:** Harmonic mean of precision and recall

## Complete Benchmarking Pipeline

```bash
# 1. Generate or download benchmark data
python simulate_data.py --transcripts 500 --reads 100000 --output benchmark_data

# 2. Run the assembler
cargo run -- assemble -i benchmark_data/reads.fastq.gz -o assembled --isoforms --gtf assembled.gtf

# 3. Evaluate the results
cargo run -- gtf-compare --truth benchmark_data/truth.gtf --pred assembled.gtf --output metrics.tsv

# 4. Analyze the results
python analyze_results.py metrics.tsv
``` 