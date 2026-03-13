# Scaffold+Polish Component Bench

This component isolates the shared paired-end postprocess path in
[scaffolder.rs](/home/jake/Projects/Raptor/src/pipeline/scaffolder.rs).

It exercises `scaffold_and_polish_contigs` directly on small fixed synthetic
tasks instead of paying full assembly cost.

## Commands

Generate one synthetic task:

```bash
cargo run -- component-bench prepare-scaffold-polish \
  --output artifacts/autoresearch_raptor/scaffold_polish/task_000 \
  --scaffold-groups 2 \
  --contigs-per-scaffold 3 \
  --contig-len 400 \
  --read-len 80 \
  --insert-size 200 \
  --internal-pairs-per-contig 12 \
  --true-link-pairs 5 \
  --decoy-link-pairs 2 \
  --mutations-per-contig 1 \
  --seed 7
```

Generate a fixed panel:

```bash
cargo run -- component-bench prepare-scaffold-polish-panel \
  --output artifacts/autoresearch_raptor/scaffold_polish/panel_small \
  --tasks 16 \
  --seed 7 \
  --seed-step 1
```

Evaluate a task or panel root:

```bash
cargo run -- component-bench scaffold-polish \
  --task artifacts/autoresearch_raptor/scaffold_polish/panel_small \
  --min-scaffold-links 3 \
  --min-primary-matches 2 \
  --min-scaffold-matches 4 \
  --json
```

## Task Layout

Each task directory contains:

```text
task_dir/
  contigs.fa
  reads_1.fastq
  reads_2.fastq
  truth.json
  task.json
```

`contigs.fa` contains mutated input contigs. `truth.json` stores:

- expected scaffold contig order plus orientation
- expected polished contig sequences

The reads include:

- same-contig pairs for insert-size estimation and polishing coverage
- true cross-contig pairs for real scaffold joins
- false cross-contig decoy pairs below the conservative support threshold

## Metrics

The harness reports:

- scaffold link precision / recall / F1
- exact scaffold recovery rate
- polished exact-contig recovery rate
- polished base accuracy
- throughput in pairs/second
- `score`

Current evaluator knobs:

- `min_scaffold_links`
- `min_primary_matches`
- `min_scaffold_matches`

Current score:

```text
score = 0.6 * scaffold_link_f1 + 0.4 * polished_exact_rate
```

That keeps the local objective anchored on postprocess correctness, with a
slight bias toward scaffold-link quality.

## Current Status

The optimizer bridge is now live end to end:

- `scripts/autoresearch_bridge/scaffold_polish_prepare.py`
- `scripts/autoresearch_bridge/scaffold_polish_train.py`
- `scripts/autoresearch_bridge/promote_scaffold_polish_candidate.py`

The fixed panel now mixes:

- balanced support cases
- sparse true-link cases that punish over-conservative thresholds on recall
- precision-edge cases that punish overly permissive thresholds on false joins

Current local result on the fixed mixed-support panels:

- baseline/default `min_scaffold_links=3`
- train score `0.7427`
- val score `0.7479`
- heldout score `0.7569`
- heldout-eligible candidates: `1`

Initial release promotion showed why the stronger local panel was necessary:

- candidate `min_scaffold_links=5` matched baseline assembly quality but was
  slightly slower on `quick_test`
- the promotion gate now requires a meaningful quality or runtime win before a
  scaffold+polish candidate can auto-promote

Current interpretation:

- the scaffold+polish optimizer surface can now reject bad default changes
  locally
- the existing production default `min_scaffold_links=3` remains the correct
  setting on the current fixed panels

## Guardrails

Do not promote a local winner on this task alone.

A shared scaffold+polish candidate is only acceptable if:

1. it improves scaffold-link and/or polished-contig validation score
2. it does not regress `quick_test`
3. it does not materially worsen polished-assembly reference metrics
