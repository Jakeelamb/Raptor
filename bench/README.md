# Benchmarking Setup

This directory currently contains three benchmark tracks:

- Trinity parity benchmarking in `bench/trinity_parity/`
- transcript/isoform benchmarking for GTF-aware evaluation
- genome assembly benchmarking in `bench/genome_assembly/`

For whole-genome assembler comparisons, use the dedicated workflow documented in `bench/genome_assembly/README.md`.
For Trinity replacement work, use `bench/trinity_parity/README.md`. That harness is the authority for Raptor-vs-Trinity transcriptome assembly evidence.

The older transcript/isoform section below describes the intended GTF-aware evaluation shape. Treat it as historical until the files it names are replaced by the Trinity parity harness.

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

Historical intended workflow. The old `simulate_data.py` script is not currently present in this directory; use `trinity_parity/run_tiny_fixture.py` for the active deterministic fixture.

```bash
python3 bench/trinity_parity/run_tiny_fixture.py
```

This writes generated truth, reads, Raptor output, and `report.json` under `target/trinity_parity/tiny_alt_isoform/`.

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
# Active tiny fixture harness
python3 bench/trinity_parity/run_tiny_fixture.py

# Existing GTF comparison command, once truth/predicted GTFs exist
cargo run -- gtf-compare --truth truth.gtf --pred assembled.gtf --output metrics.tsv
```
