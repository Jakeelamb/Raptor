# Error-Correction Component Bench

This component isolates phase-3 singleton rescue in
[large_genome_assembler.rs](/home/jake/Projects/Raptor/src/pipeline/large_genome_assembler.rs).

It evaluates `error_correct_kmers` on fixed synthetic count tables instead of
paying full assembly cost.

## Commands

Generate one synthetic task:

```bash
cargo run -- component-bench prepare-error-correction \
  --output artifacts/autoresearch_raptor/error_correction/task_000 \
  --k 21 \
  --trusted-roots 64 \
  --weak-roots 32 \
  --correctable-singletons-per-trusted 2 \
  --protected-singletons-per-weak 1 \
  --seed 7
```

Generate a fixed panel:

```bash
cargo run -- component-bench prepare-error-correction-panel \
  --output artifacts/autoresearch_raptor/error_correction/panel_small \
  --tasks 16 \
  --seed 7 \
  --seed-step 1
```

Evaluate a task or panel root:

```bash
cargo run -- component-bench error-correction \
  --task artifacts/autoresearch_raptor/error_correction/panel_small \
  --min-count 1 \
  --min-trusted-count 4 \
  --json
```

## Task Layout

Each task directory contains:

```text
task_dir/
  counts.json
  truth.json
  task.json
```

`counts.json` stores the synthetic input count table. `truth.json` stores:

- singleton k-mers that must be removed
- singleton k-mers that must be preserved
- trusted k-mers with their expected post-correction counts

The synthetic generator mixes:

- trusted roots with counts cycling through `4`, `6`, and `9`
- correctable singleton errors one edit away from trusted roots
- weak roots with counts `2` and `3`
- protected singleton neighbors around weak roots that should survive correction

## Metrics

The harness reports:

- correction precision / recall / F1
- protected-singleton retention
- trusted-target exact-count rate
- throughput in k-mers/second
- `score`

Current evaluator knobs:

- `min_count`
- `min_trusted_count`

Current score:

```text
score =
  0.50 * correction_f1 +
  0.25 * preserved_retention_rate +
  0.25 * trusted_target_exact_rate
```

## Current Status

The optimizer-facing surface is now wide enough to punish both overly lenient
and overly strict singleton rescue.

The first Raptor-side autoresearch bridge now exists:

- [scripts/autoresearch_bridge/error_correction_prepare.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/error_correction_prepare.py)
- [scripts/autoresearch_bridge/error_correction_train.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/error_correction_train.py)
- [docs/autoresearch-error-correction-bridge.md](/home/jake/Projects/Raptor/docs/autoresearch-error-correction-bridge.md)

Latest heldout checks on March 13, 2026:

- baseline `min_count=1, min_trusted_count=4`: score `1.0`, correction F1
  `1.0`, preserved retention `1.0`, trusted exact `1.0`
- overly lenient `min_count=1, min_trusted_count=3`: score `0.7800`,
  correction F1 `0.8101`, preserved retention `0.5`, trusted exact `1.0`
- overly strict `min_count=1, min_trusted_count=6`: score `0.8114`,
  correction F1 `0.7942`, preserved retention `1.0`, trusted exact `0.6573`

Interpretation:

- the panel now gives a real two-sided local gradient around the live baseline
- overly lenient trusted floors are penalized on preserved rare-variant
  retention before full assembly runs
- overly strict trusted floors are penalized on correction recall and exact
  trusted-target recovery
- the heldout-aware train/val/heldout bridge still converges to the current
  baseline `min_count=1, min_trusted_count=4`
- the next step is to add an ambiguity guard or harder near-trusted heldout
  cases before spending time on a promotion script
