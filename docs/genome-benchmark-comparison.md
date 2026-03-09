# Genome Benchmark Comparison

This document is generated from `bench/genome_assembly/benchmark_summary.csv`.

It is the publication-facing comparison view for Raptor's genome assembly benchmarks.

## Current Status

- Raptor baseline runs are available in the current repo state.
- Comparator conclusions should only be made for datasets where the same summary contains both Raptor and comparator rows with non-`N/A` metrics.

## Dataset: `drosophila`

Tools present in summary: spades

| Tool | Time (s) | Peak RSS (KB) | Contigs | N50 | Genome Fraction % | Misassemblies | NGA50 | Duplication Ratio |
|------|----------|---------------|---------|-----|-------------------|---------------|-------|-------------------|
| spades | N/A | N/A | N/A | N/A | N/A | N/A | N/A | N/A |

Comparator metrics are not yet complete for this dataset; treat this as an execution record, not a winner table.

## Dataset: `quick_test`

Tools present in summary: raptor, raptor, raptor, spades

| Tool | Time (s) | Peak RSS (KB) | Contigs | N50 | Genome Fraction % | Misassemblies | NGA50 | Duplication Ratio |
|------|----------|---------------|---------|-----|-------------------|---------------|-------|-------------------|
| raptor | N/A | N/A | N/A | N/A | N/A | N/A | N/A | N/A |
| raptor | N/A | N/A | N/A | N/A | N/A | N/A | N/A | N/A |
| raptor | N/A | N/A | N/A | N/A | N/A | N/A | N/A | N/A |
| spades | N/A | N/A | N/A | N/A | N/A | N/A | N/A | N/A |

Comparator metrics are not yet complete for this dataset; treat this as an execution record, not a winner table.

## Update Workflow

```bash
python3 ./bench/genome_assembly/summarize_results.py
python3 ./bench/genome_assembly/render_comparison_report.py
```

