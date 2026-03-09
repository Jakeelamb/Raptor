# Genome Benchmark Comparison

This document is generated from `bench/genome_assembly/benchmark_summary.csv`.

It is the publication-facing comparison view for Raptor's genome assembly benchmarks.

## Current Status

- Raptor baseline runs are available in the current repo state.
- Comparator conclusions should only be made for datasets where the same summary contains both Raptor and comparator rows with non-`N/A` metrics.

## Dataset: `drosophila`

Latest complete rows by tool: spades

| Tool | Source Run | Time (s) | Peak RSS (KB) | Contigs | N50 | Genome Fraction % | Misassemblies | NGA50 | Duplication Ratio |
|------|------------|----------|---------------|---------|-----|-------------------|---------------|-------|-------------------|
| spades | drosophila_20260130_153531 | N/A | N/A | N/A | N/A | N/A | N/A | N/A | N/A |

Comparator metrics are not yet complete for this dataset; treat this as an execution record, not a winner table.

## Dataset: `quick_test`

Latest complete rows by tool: raptor, spades

| Tool | Source Run | Time (s) | Peak RSS (KB) | Contigs | N50 | Genome Fraction % | Misassemblies | NGA50 | Duplication Ratio |
|------|------------|----------|---------------|---------|-----|-------------------|---------------|-------|-------------------|
| raptor | quick_test_20260309_222334 | 575.48 | N/A | 23 | 255072 | 90.723 | 2 | 105008 | 1.001 |
| spades | quick_test_20260309_140001 | 232.424 | N/A | 37 | 130153 | 88.601 | 0 | 120152 | 1.000 |

Comparator metrics are present for this dataset.

## Update Workflow

```bash
python3 ./bench/genome_assembly/summarize_results.py
python3 ./bench/genome_assembly/render_comparison_report.py
```

