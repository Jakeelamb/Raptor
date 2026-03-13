# Autoresearch Branch-Resolution Bridge

This bridge keeps branch-choice optimization inside Raptor.

Unlike read mapping, the initial branch-resolution search space is small enough
to evaluate exhaustively instead of mutating stochastically.

## Files

- [scripts/autoresearch_bridge/branch_resolution_prepare.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/branch_resolution_prepare.py)
- [scripts/autoresearch_bridge/branch_resolution_train.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/branch_resolution_train.py)

## Prepare Panels

```bash
python3 scripts/autoresearch_bridge/branch_resolution_prepare.py \
  --binary target/debug/raptor \
  --output-root artifacts/autoresearch_raptor/branch_resolution/tasks \
  --train-tasks 16 \
  --val-tasks 8 \
  --heldout-tasks 8 \
  --cases-per-task 64 \
  --heldout-cases-per-task 96
```

This writes:

```text
artifacts/autoresearch_raptor/branch_resolution/tasks/
  train/
  val/
  heldout/
  metadata.json
```

## Run The Search Loop

```bash
python3 scripts/autoresearch_bridge/branch_resolution_train.py \
  --binary target/release/raptor \
  --tasks-root artifacts/autoresearch_raptor/branch_resolution/tasks \
  --require-heldout
```

## Promote A Winner

```bash
python3 scripts/autoresearch_bridge/promote_branch_resolution_candidate.py \
  --binary target/release/raptor \
  --tasks-root artifacts/autoresearch_raptor/branch_resolution/tasks \
  --train-result artifacts/autoresearch_raptor/branch_resolution/latest_result.json
```

That reruns the candidate on the fixed branch-resolution train/val/heldout
panels and, unless disabled, runs the current `quick_test` end-to-end
assembler command.

The promotion summary now includes:

- the candidate metrics
- the default baseline metrics
- heldout scenario deltas for the stress cases
- minimap2-based evaluation of the final polished assembly
- a final verdict of `promote_default`, `review`, or `reject`
- explicit reasons when a candidate is held for review

## Current Mutable Surface

This bridge currently searches:

- `branch_support_min_win`
- `branch_support_min_margin`
- `prefer_non_repeat`

That is enough to vary the production heuristic meaningfully while keeping the
promotion surface small and auditable.

## Current Status

This bridge evaluates the real branch-resolution path inside
[large_genome_assembler.rs](/home/jake/Projects/Raptor/src/pipeline/large_genome_assembler.rs),
not a Python reimplementation.

Latest autonomous run:

- promoted default: `branch_support_min_win=1`,
  `branch_support_min_margin=2`, `prefer_non_repeat=true`
- synthetic train/val/heldout score: `0.9658` / `0.9590` / `1.0000`
- heldout stress panel improved `support_margin_gate` exact rate from `0.0` to
  `1.0` without regressing the other guarded scenarios
- promotion summary:
  [summary.json](/home/jake/Projects/Raptor/artifacts/autoresearch_raptor/branch_resolution/promotions/branch_resolution_20260311_184053/summary.json)
- promoted candidate `quick_test`: `97.53s`, `22` contigs, contig `N50 120054 bp`,
  `19` scaffolds, scaffold `N50 360080 bp`, `369` polish corrections
- polished-assembly reference coverage: `0.8953175` (matches baseline)
- polished matched-base fraction: `0.9926612` vs baseline `0.9926444`
- live default `quick_test`: `98.40s`, `22` contigs, contig `N50 120054 bp`,
  `19` scaffolds, scaffold `N50 360080 bp`, `369` polish corrections
- latest verdict: `promote_default`

Next improvements for this bridge:

1. keep the promoted branch defaults fixed unless a later end-to-end
   regression points back to this heuristic
2. decide whether `support_dominates` needs its own guarded heldout split if a
   future regression appears there
3. add peak RSS and per-phase timing capture to promotion summaries
4. move the next optimizer loop to the shared scaffold + polish evidence layer
