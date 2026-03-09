# Genome Benchmark Report

Date: 2026-03-09

This document records the current checked-in benchmark evidence for Raptor's large-genome assembly pipeline.

## Environment

- Host OS: Linux 6.18.9-arch1-2 x86_64
- CPU threads available: 24
- System memory: 98,130,588 kB (about 93.6 GiB)
- Build: `cargo build --release`

## Dataset

- Benchmark dataset: `quick_test`
- Reference size: 2,000,000 bp
- Input reads: paired-end FASTQ
- Reads processed: 666,666
- Bases processed: 99,999,900

The `quick_test` dataset is intended as a reproducible local baseline, not as a publication-grade benchmark on its own.

## Command

```bash
target/release/raptor assemble-large \
  -i bench/genome_assembly/data/quick_test/reads_1.fastq.gz \
  --input2 bench/genome_assembly/data/quick_test/reads_2.fastq.gz \
  -o artifacts/benchmarks/quick_test_20260309_raptor/contigs.fa \
  -t 8 \
  -k 31 \
  --min-count 0 \
  --scaffold \
  --polish
```

## Baseline Results

### Assembly

- End-to-end runtime: 247.54 s
- Disk used during k-mer bucketing: 0.64 GB
- Adaptive min-count selected: 3
- Total k-mers: 79,999,920
- Unique k-mers: 20,298,322
- Filtered k-mers: 1,841,744
- Error-corrected k-mers: 14,981,404
- Tips removed: 2,508
- Bubbles popped: 2

### Raw Contigs

- Contigs: 24
- Total assembled bases: 1,817,131 bp
- Reference recovery by assembled length: 90.86%
- QUAST genome fraction: 90.651%
- N50: 130,057 bp
- Largest contig: 420,010 bp

### Scaffolds

- Scaffolds: 21
- Total scaffolded bases: 1,817,131 bp
- Scaffold N50: 360,081 bp
- Largest scaffold: 575,626 bp

### Polishing

- Corrections made: 1,082
- Correction types: 1,082 substitutions, 0 insertions, 0 deletions

## Interpretation

This run shows that the checked-in `assemble-large` pipeline is operational end-to-end on a nontrivial paired-end dataset and produces multi-hundred-kilobase scaffold continuity from a local benchmark.

## Comparator Run

On 2026-03-09, the same `quick_test` dataset was also compared against SPAdes 4.2.0 with QUAST 5.3.0 evaluation after enabling adaptive cutoff selection in `assemble-large`.

### Comparator Summary

| Tool | Runtime (s) | Contigs | N50 | Genome Fraction (%) | NGA50 | Misassemblies |
|------|-------------|---------|-----|---------------------|-------|---------------|
| SPAdes | 232.424 | 37 | 130,153 | 88.601 | 120,152 | 0 |
| Raptor | 247.54 | 24 | 130,057 | 90.651 | 111,000 | 2 |

### Comparator Interpretation

The large regression on `quick_test` turned out to be a cutoff-policy problem, not just a graph-traversal problem. With adaptive min-count selection enabled, Raptor recovers from 1,681 fragmented contigs to 24 contigs and reaches practical parity with SPAdes on raw contig N50 while exceeding SPAdes on genome fraction for this dataset.

This is not a clean win yet. SPAdes is still faster, still has the better NGA50, and still avoids the two misassemblies QUAST reports for the current Raptor output.

The current evidence is good enough for:

- regression tracking
- performance profiling
- release-note claims about implemented functionality

It is not yet strong enough for:

- broad superiority claims against established assemblers
- publication-grade comparative benchmarking

## What Is Missing

To support stronger external claims, the benchmark suite still needs:

- at least one larger public real-world benchmark
- peak RSS memory capture on the benchmark host
- repeated runs or variance estimates
- alignment-aware improvements that reduce the remaining misassemblies and close the NGA50 gap

## Next Benchmark Priorities

1. Improve branch resolution and polishing so the raw contigs keep this continuity while removing the two QUAST misassemblies.
2. Add at least one larger public dataset beyond `quick_test`.
3. Capture peak RSS memory reliably for both tools on this host.
4. Repeat each benchmark enough times to report variance instead of a single run.
