# Autoresearch Contig-Extraction Bridge

This bridge lives inside Raptor because the current external
`/home/jake/Projects/autoresearch` repo does not yet have a plugin interface for
arbitrary component scorers.

It mirrors the existing Raptor-side component bridges:

- a fixed train/val task root
- a held-out split for non-regression checks
- a tiny mutable search surface
- a scalar local score

## Files

- [scripts/autoresearch_bridge/contig_extraction_prepare.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/contig_extraction_prepare.py)
- [scripts/autoresearch_bridge/contig_extraction_train.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/contig_extraction_train.py)
- [scripts/autoresearch_bridge/promote_contig_extraction_candidate.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/promote_contig_extraction_candidate.py)
- [docs/contig-extraction-component-bench.md](/home/jake/Projects/Raptor/docs/contig-extraction-component-bench.md)

## Prepare Panels

```bash
python3 scripts/autoresearch_bridge/contig_extraction_prepare.py \
  --binary target/debug/raptor \
  --output-root artifacts/autoresearch_raptor/contig_extraction/tasks
```

This writes:

```text
artifacts/autoresearch_raptor/contig_extraction/tasks/
  train/
  val/
  heldout/
  metadata.json
```

Default panel profiles deliberately mix:

- branching tasks that reward high-count seed ranking
- close-margin branching tasks that punish low-count-first extraction
- repeat-completion tasks that require a second repeat-seed pass
- repeat-fallback tasks where only repeat-classified seeds survive cleanup
- heldout repeat-priority tasks where a high-count repeat branch can steal the
  shared suffix unless non-repeat seeds are processed first under redundant
  contig suppression

## Run The Search Loop

```bash
python3 scripts/autoresearch_bridge/contig_extraction_train.py \
  --binary target/debug/raptor \
  --tasks-root artifacts/autoresearch_raptor/contig_extraction/tasks \
  --require-heldout
```

## Run The Promotion Gate

```bash
python3 scripts/autoresearch_bridge/promote_contig_extraction_candidate.py \
  --binary target/release/raptor \
  --tasks-root artifacts/autoresearch_raptor/contig_extraction/tasks \
  --train-result artifacts/autoresearch_raptor/contig_extraction/latest_result.json \
  --threads 8
```

This promotion pass compares the heldout-selected candidate against the
baseline on:

- fixed train/val/heldout panels
- `quick_test`
- a cached repeat-heavy `drosophila` subset under
  `artifacts/autoresearch_raptor/contig_extraction/repeat_heavy_subset/`

The full contig-extraction surface is now wired into live `assemble-large`, so
promotion runs can test seed-policy changes end to end without code edits.

## Current Mutable Surface

Right now the bridge searches:

- `prefer_high_count_seeds`
- `prefer_non_repeat_seeds`
- `enable_repeat_seed_completion`
- `suppress_redundant_contigs`

The selection tie-break now prefers the smallest safe change over gratuitous
knob churn when multiple candidates tie on heldout quality.

The heldout split is stricter than train/val now. It specifically includes
repeat-priority tasks that expose a blind spot in the earlier panel:
`prefer_non_repeat_seeds=false` could still look harmless after
`suppress_redundant_contigs=true`, even though repeat-first seeding can claim a
shared suffix and prevent the true non-repeat contig from being extracted.

Current checked result:

- heldout-selected winner: keep `prefer_non_repeat_seeds=true` and add
  `suppress_redundant_contigs=true`
- train-only winner is now stricter and wrong in exactly the way we wanted the
  panel to expose:
  `prefer_non_repeat_seeds=false`, `suppress_redundant_contigs=true`
- heldout-selected score is `1.0 / 1.0 / 1.0`
- the repeat-first ablation on heldout now drops to score `0.9695`, with
  aggregate exact contig rate falling to `0.9394`
- the two new `repeat_priority` heldout tasks each collapse from exact-contig
  `1.0` to `0.0` and truth-kmer `F1` `1.0` to `0.3125` when
  `prefer_non_repeat_seeds` is disabled under suppression
- baseline heldout score before suppression is still `0.8691`
- baseline heldout over-extraction before suppression is still `62` observed
  contigs vs `34` expected
- the latest promotion artifact is
  [summary.json](/home/jake/Projects/Raptor/artifacts/autoresearch_raptor/contig_extraction/promotions/contig_extraction_20260313_120134/summary.json)
- promotion status is `review`, not `promote_default`
- heldout deltas versus baseline are real: score `+0.1369`, truth-kmer `F1`
  `+0.1310`, count-agreement `+0.4762`
- `quick_test` remains output-identical to baseline on this dataset (`22`
  contigs, `19` scaffolds, `369` polish corrections), with candidate
  wall-clock `78.45s` vs baseline `75.75s`
- the cached repeat-heavy subset reduces final contigs by `87` (`20092` vs
  `20179`) but regresses runtime by about `297s` (`626.91s` vs `329.76s`)
- the candidate repeat-heavy run reports `Redundant-contig suppression removed
  0 contigs`, so the current stress suite still is not exercising the local
  failure mode directly enough

Interpretation:

- the local failure mode is real and now specifically benchable on heldout
- the current `quick_test` dataset still does not exercise it strongly enough
  to justify a production default flip by itself
- the new promotion gate is doing its job: it preserves the local win artifact
  but keeps the live default unchanged when repeat-heavy end-to-end evidence is
  ambiguous
- the next high-value work is to add phase-level timing or a better
  repeat-heavy suite that triggers explicit redundant-contig suppression
