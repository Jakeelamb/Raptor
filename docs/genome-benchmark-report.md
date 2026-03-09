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
  --min-count 3 \
  --scaffold \
  --polish
```

## Results

### Assembly

- End-to-end runtime: 289 s
- Disk used during k-mer bucketing: 0.64 GB
- Total k-mers: 79,999,920
- Unique k-mers: 20,298,322
- Filtered k-mers: 1,841,744
- Error-corrected k-mers: 14,981,045
- Tips removed: 2,508
- Bubbles popped: 2

### Raw Contigs

- Contigs: 24
- Total assembled bases: 1,817,131 bp
- Reference recovery by assembled length: 90.86%
- N50: 120,054 bp
- Largest contig: 420,010 bp

### Scaffolds

- Scaffolds: 21
- Total scaffolded bases: 1,817,134 bp
- Scaffold N50: 360,081 bp
- Largest scaffold: 575,626 bp

### Polishing

- Corrections made: 986
- Correction types: 986 substitutions, 0 insertions, 0 deletions

## Interpretation

This run shows that the checked-in `assemble-large` pipeline is operational end-to-end on a nontrivial paired-end dataset and produces multi-hundred-kilobase scaffold continuity from a local benchmark.

The current evidence is good enough for:

- regression tracking
- performance profiling
- release-note claims about implemented functionality

It is not yet strong enough for:

- broad superiority claims against established assemblers
- publication-grade comparative benchmarking

## What Is Missing

To support stronger external claims, the benchmark suite still needs:

- comparator runs from the same scripts and hardware
- at least one larger public real-world benchmark
- QUAST-style evaluation metrics
- peak RSS memory capture on the benchmark host
- repeated runs or variance estimates

## Next Benchmark Priorities

1. Install and benchmark SPAdes with the same datasets and thread counts.
2. Add QUAST output parsing into the benchmark summary pipeline.
3. Add at least one larger public dataset beyond `quick_test`.
4. Aggregate results into a single cross-tool summary table for README and release notes.
