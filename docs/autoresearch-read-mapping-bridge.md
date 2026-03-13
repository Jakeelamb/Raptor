# Autoresearch Read-Mapping Bridge

This bridge lives inside Raptor because the current external
`/home/jake/Projects/autoresearch` repo does not yet have a plugin interface for
arbitrary component scorers.

It mirrors the shape of `autoresearch`:

- a fixed train/val task root
- a small mutable search surface
- a scalar score
- a time-budgeted search loop
- an optional multi-fidelity successive-halving policy for the large mapper grid

## Files

- [scripts/autoresearch_bridge/read_mapping_prepare.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/read_mapping_prepare.py)
- [scripts/autoresearch_bridge/read_mapping_train.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/read_mapping_train.py)

## Prepare Panels

```bash
python3 scripts/autoresearch_bridge/read_mapping_prepare.py \
  --binary target/debug/raptor \
  --output-root artifacts/autoresearch_raptor/read_mapping/tasks \
  --train-tasks 16 \
  --val-tasks 8 \
  --ambiguous-repeat-decoy-rate 0.15 \
  --paired
```

Default panel profile now deliberately includes:

- shorter reads
- longer planted repeats
- substitution noise
- random decoy reads that should stay unmapped
- repeat-only ambiguous decoys that should not produce confident scaffold hits

This writes:

```text
artifacts/autoresearch_raptor/read_mapping/tasks/
  train/
  val/
  heldout/
  metadata.json
```

## Run The Search Loop

```bash
python3 scripts/autoresearch_bridge/read_mapping_train.py \
  --binary target/debug/raptor \
  --tasks-root artifacts/autoresearch_raptor/read_mapping/tasks \
  --time-budget 120
```

Candidate-level racing on the large grid:

```bash
python3 scripts/autoresearch_bridge/read_mapping_train.py \
  --binary target/debug/raptor \
  --tasks-root artifacts/autoresearch_raptor/read_mapping/tasks \
  --require-heldout \
  --search-policy successive-halving \
  --halving-factor 4 \
  --final-heldout-candidates 32
```

Developer smoke:

```bash
python3 scripts/autoresearch_bridge/read_mapping_train.py \
  --binary target/debug/raptor \
  --tasks-root /tmp/raptor_autoresearch_mapper \
  --smoke
```

## Promote A Winner

```bash
python3 scripts/autoresearch_bridge/promote_read_mapping_candidate.py \
  --binary target/release/raptor \
  --tasks-root artifacts/autoresearch_raptor/read_mapping/tasks \
  --train-result artifacts/autoresearch_raptor/read_mapping/latest_result.json
```

That reruns the winner on the fixed read-mapping train/val panels and, unless
disabled, runs the paired `quick_test` assembler command with the same mapper
knobs threaded into live post-processing.

The promotion summary now includes:

- the candidate metrics
- the default baseline metrics
- heldout metrics for candidate and baseline when available
- candidate-minus-baseline deltas
- minimap2-based reference evaluation against `quick_test/reference.fa`
- a final verdict of `promote_default`, `review`, or `reject`
- explicit reasons when a candidate is held for review

## Current Mutable Surface

Right now the bridge searches:

- minimizer `k`
- minimizer `w`
- minimum primary-hit support
- minimum scaffold-hit support

That is still intentionally small. The goal is to optimize a real tradeoff
surface without letting the search sprawl across the whole mapper.

## Current Status

The initial easy smoke panel saturated too often. The bridge now defaults to a
harder panel so configs separate on both sensitivity and specificity.

The local objective is still synthetic, but the promotion path now uses the
same `k`, `w`, primary-hit threshold, and scaffold-hit threshold in:

- standalone polishing
- standalone scaffolding
- the shared scaffold+polish fast path
- `assemble-large` `quick_test` promotion runs

Latest autonomous comparison:

- promoted default: `k=15`, `w=4`, `min_primary_matches=2`,
  `min_scaffold_matches=4`
- synthetic train/val/heldout score: `0.8407` / `0.8293` / `0.8025`
- heldout scaffold specificity held steady at `0.5488`
- heldout pair-link contig correctness improved from `0.7695` to `0.8242`
- live default `quick_test`: `96.38s`, `22` contigs, contig `N50 120054 bp`,
  `19` scaffolds, scaffold `N50 360080 bp`, `370` polish corrections
- polished-assembly reference eval: coverage `0.8953175`,
  matched-base fraction `0.9926444`
- baseline polished-assembly reference eval: coverage `0.8953175`,
  matched-base fraction `0.9925343`
- latest verdict: `promote_default`

Latest implementation update:

- `read_mapping_train.py` now runs a deterministic heldout-aware grid search
  instead of a stochastic local search
- `read_mapping_train.py` now also supports
  `--search-policy successive-halving` so the largest mapper grid does not have
  to pay full `train + val + heldout` cost for every config
- the read-mapping component bench now records paired-link correctness and
  paired-link specificity, not just per-read metrics
- `promote_read_mapping_candidate.py` now evaluates the final polished assembly
  and only treats lower correction counts as suspicious when polished quality
  actually regresses
- production read-mapping defaults in the CLI and `ReadMappingConfig::default()`
  now match the promoted config
- the default campaign manifest now opts the mapper into successive halving
  with `--halving-factor 4` and `--final-heldout-candidates 32`

Latest successive-halving smoke:

- manifest:
  [read_mapping_halving_smoke.json](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/manifests/read_mapping_halving_smoke.json)
- smoke grid size: `81` configs
- split evaluations under successive halving: `81 train`, `27 val`, `12 heldout`
- total split evaluations: `120` instead of the `243` evaluations a full
  smoke grid would pay
- heldout-gated selected config on that smoke run remained the safe baseline
  `k=15`, `w=10`, `min_primary_matches=3`, `min_scaffold_matches=2`

What remains:

- `quick_test` is still the final promotion gate; local mapper score is good
  enough to guide search, not to replace end-to-end validation
- `branch_resolution` is now promoted too, so the next component loop should
  move to shared scaffold + polish evidence
- if downstream scaffold/polish quality becomes the next blocker, we should add
  a dedicated shared-evidence component harness instead of forcing mapper panels
  to model the full phenotype

Next best implementation task for this bridge:

1. keep the mapper bridge stable and only return to it if the next end-to-end
   regression traces back to post-assembly read mapping

How to validate this bridge now:

```bash
python3 scripts/autoresearch_bridge/read_mapping_prepare.py \
  --binary target/debug/raptor \
  --output-root artifacts/autoresearch_raptor/read_mapping/tasks \
  --paired
target/debug/raptor component-bench read-mapping \
  --task artifacts/autoresearch_raptor/read_mapping/tasks/val \
  --k 15 --w 4 \
  --min-primary-matches 2 \
  --min-scaffold-matches 4 \
  --position-tolerance 8 \
  --json
python3 scripts/autoresearch_bridge/read_mapping_train.py \
  --binary target/release/raptor \
  --tasks-root artifacts/autoresearch_raptor/read_mapping/tasks \
  --require-heldout \
  --search-policy successive-halving \
  --halving-factor 4 \
  --final-heldout-candidates 32
python3 scripts/autoresearch_bridge/promote_read_mapping_candidate.py \
  --binary target/release/raptor \
  --tasks-root artifacts/autoresearch_raptor/read_mapping/tasks \
  --train-result artifacts/autoresearch_raptor/read_mapping/latest_result.json
```

Next improvements for this bridge:

1. keep the current promoted mapper defaults fixed until a later component gate
   says otherwise
2. use the same train/heldout/promote discipline for the shared scaffold +
   polish evidence layer
3. only return to mapper work if a later end-to-end regression traces back to
   post-assembly mapping again
