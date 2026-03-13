# Contig-Extraction Component Bench

This component isolates phase-5 contig extraction in
[large_genome_assembler.rs](/home/jake/Projects/Raptor/src/pipeline/large_genome_assembler.rs).

It evaluates `build_contigs_from_graph` on fixed synthetic graph fixtures
instead of paying full assembly cost.

## Commands

Generate one synthetic task:

```bash
cargo run -- component-bench prepare-contig-extraction \
  --output artifacts/autoresearch_raptor/contig_extraction/task_000 \
  --profile branching \
  --k 11 \
  --component-count 6 \
  --primary-reads-per-component 6 \
  --alternate-reads-per-component 2 \
  --seed 7
```

Generate a fixed panel:

```bash
cargo run -- component-bench prepare-contig-extraction-panel \
  --output artifacts/autoresearch_raptor/contig_extraction/panel_small \
  --tasks 8 \
  --profile mixed \
  --seed 7 \
  --seed-step 1
```

Evaluate a task or panel root:

```bash
cargo run -- component-bench contig-extraction \
  --task artifacts/autoresearch_raptor/contig_extraction/panel_small \
  --json
```

Useful ablations:

```bash
cargo run -- component-bench contig-extraction \
  --task artifacts/autoresearch_raptor/contig_extraction/panel_small \
  --disable-prefer-high-count-seeds \
  --json

cargo run -- component-bench contig-extraction \
  --task artifacts/autoresearch_raptor/contig_extraction/panel_small \
  --disable-repeat-seed-completion \
  --json

cargo run -- component-bench contig-extraction \
  --task artifacts/autoresearch_raptor/contig_extraction/panel_small \
  --suppress-redundant-contigs \
  --json
```

## Task Layout

Each task directory contains:

```text
task_dir/
  counts.json
  graph.json
  truth.json
  task.json
```

- `counts.json` stores canonical k-mer counts, including extra decoys that can
  influence repeat classification without surviving in the cleaned graph.
- `graph.json` stores the surviving canonical k-mers and any branch-support
  evidence used during extraction.
- `truth.json` stores the expected contig sequences.

## Metrics

The harness reports:

- exact contig recovery rate
- truth k-mer precision / recall / F1
- contig-count agreement
- throughput in graph-nodes/second
- `score`

Current score:

```text
score =
  0.50 * truth_kmer_f1 +
  0.35 * exact_contig_rate +
  0.15 * contig_count_agreement
```

## Current Mutable Surface

The first optimizer-facing surface searches:

- `prefer_high_count_seeds`
- `prefer_non_repeat_seeds`
- `enable_repeat_seed_completion`
- `suppress_redundant_contigs`

`suppress_redundant_contigs` marks enclosed leftover branch paths as redundant
when both ends are already bounded by an assembled backbone. On the current
fixed panel, that is the first contig-local knob that removes the synthetic
over-extraction failure mode instead of only reordering seeds.
