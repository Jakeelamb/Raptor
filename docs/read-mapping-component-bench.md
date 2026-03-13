# Read-Mapping Component Bench

This is the first Raptor component harness intended for `autoresearch`-style
optimization loops.

It isolates the shared mapper used by polishing and shared scaffold+polish
postprocessing:

- mapper index build
- read-to-contig positional hit
- read-to-contig scaffold support hit

The goal is to optimize this component locally without rerunning the entire
assembler.

## Commands

Generate a synthetic task:

```bash
cargo run -- component-bench prepare-read-mapping \
  --output artifacts/autoresearch_raptor/read_mapping/small_pair_panel/task_000 \
  --paired \
  --num-contigs 4 \
  --contig-len 3000 \
  --read-len 150 \
  --reads 256 \
  --insert-size 450 \
  --repeat-len 80 \
  --error-rate 0.01 \
  --decoy-rate 0.25 \
  --ambiguous-repeat-decoy-rate 0.15 \
  --seed 7 \
  --k 11 \
  --w 5
```

Generate a repeatable panel for search:

```bash
cargo run -- component-bench prepare-read-mapping-panel \
  --output artifacts/autoresearch_raptor/read_mapping/small_pair_panel \
  --tasks 16 \
  --paired \
  --num-contigs 4 \
  --contig-len 3000 \
  --read-len 150 \
  --reads 256 \
  --insert-size 450 \
  --repeat-len 80 \
  --error-rate 0.01 \
  --decoy-rate 0.25 \
  --ambiguous-repeat-decoy-rate 0.15 \
  --seed 7 \
  --seed-step 1 \
  --k 11 \
  --w 5
```

Evaluate one task or a root containing many task directories:

```bash
cargo run -- component-bench read-mapping \
  --task artifacts/autoresearch_raptor/read_mapping/small_pair_panel \
  --k 11 \
  --w 5 \
  --min-primary-matches 3 \
  --min-scaffold-matches 2 \
  --position-tolerance 8 \
  --json \
  --output artifacts/autoresearch_raptor/read_mapping/latest_report.json
```

## Task Layout

Each task directory uses a small fixed contract:

```text
task_dir/
  contigs.fa
  reads_1.fastq
  reads_2.fastq            # optional
  truth.tsv
  task.json                # optional metadata
```

`truth.tsv` columns:

```text
read_id  contig  start  is_reverse  scaffold_contig  scaffold_is_reverse
```

Rules:

- `read_id` must match the FASTQ header without the leading `@`.
- `contig` and `scaffold_contig` must match FASTA record names.
- `start` is the expected 0-based contig start.
- `is_reverse` is the expected mapper strand for the positional hit.
- blank `contig` means no primary mapping is expected.
- blank `scaffold_contig` means no scaffold-support hit is expected.
- decoy reads use blank truth fields and exist to punish false-positive mapping.
- repeat-only ambiguous decoys should also use blank truth fields; they are
  intentionally shared across contigs so confident scaffold support is a false
  positive for this harness.

## Metrics

The harness emits:

- primary mapped rate
- primary exact rate
- primary near rate within `--position-tolerance`
- primary contig rate
- primary orientation rate
- unexpected primary mapping rate on decoy/unmapped reads
- primary specificity rate on decoy/unmapped reads
- scaffold mapped rate
- scaffold exact rate
- unexpected scaffold-hit rate on decoy/unmapped reads
- scaffold specificity rate on decoy/unmapped reads
- primary position MAE
- throughput in reads/second
- `score`

`score` is accuracy-weighted and intentionally does not include raw wall-clock
time directly:

- `0.40 * primary_near_rate`
- `0.20 * primary_exact_rate`
- `0.15 * scaffold_exact_rate`
- `0.05 * primary_mapped_rate`
- `0.10 * primary_specificity_rate`
- `0.10 * scaffold_specificity_rate`

Use that score under a fixed task panel and fixed runtime budget. That keeps
the optimizer focused on correctness first while still benefiting from faster
implementations through shorter evaluation cycles.

## Autoresearch Use

Use this component before touching end-to-end assembly search.

Recommended loop:

1. Generate a small panel of 8-32 synthetic tasks with different seeds.
2. Keep `k`, `w`, the support thresholds, and the score definition fixed for
   one search run.
3. Let `autoresearch` mutate only the mapper implementation or a tiny mapper
   parameter surface.
4. Accept candidates only if they also preserve `quick_test` scaffold/polish
   behavior on the end-to-end regression command.

Do not use this harness as the only quality signal for Raptor. It is a local
optimization surface, not the final acceptance criterion.

## Current Status

Implemented in the latest pass:

- unexpected scaffold hits are now counted and reported
- repeat-only ambiguous decoys can now be generated with
  `--ambiguous-repeat-decoy-rate`
- the default bridge panels now use those ambiguous decoys so local runs expose
  false scaffold support numerically

What remains:

- the new specificity signal surfaces the failure mode, but the previously
  rejected mapper candidate still wins the aggregate local score
- the next best task is to add a harder held-out split and require non-regressing
  held-out score or specificity before promotion

How to validate:

```bash
cargo test cli::component_bench::tests --lib
target/debug/raptor component-bench read-mapping \
  --task artifacts/autoresearch_raptor/read_mapping/tasks/val \
  --k 9 --w 4 \
  --min-primary-matches 5 \
  --min-scaffold-matches 3 \
  --position-tolerance 8 \
  --json
```
