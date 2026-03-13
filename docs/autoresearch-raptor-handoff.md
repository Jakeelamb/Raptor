# Autoresearch x Raptor Handoff

Date: 2026-03-13

This document is the interruption-proof execution plan for decomposing Raptor
into optimizer-friendly components and wiring an `autoresearch`-style loop
around them.

## Goal

Use small, fixed, component-local task panels so an autonomous search loop can
improve Raptor one subsystem at a time without paying full end-to-end assembly
cost on every candidate.

The near-term goal is not unrestricted code mutation. It is:

1. define stable component contracts
2. generate deterministic task panels
3. expose scalar local scores
4. guard all local wins with fixed end-to-end regression runs

## Current State

Already done:

- stage benches for large-genome branch threading, contig extraction, and error
  correction
- shared paired-end scaffold+polish fast path
- `component-bench read-mapping`
- `component-bench branch-resolution`
- `component-bench scaffold-polish`
- `component-bench error-correction`
- `component-bench contig-extraction`
- synthetic read-mapping task and panel generation
- branch-resolution task and panel generation
- error-correction train/val/heldout panel generation
- contig-extraction train/val/heldout panel generation
- Raptor-side read-mapping bridge scripts
- Raptor-side branch-resolution bridge scripts
- Raptor-side error-correction bridge scripts
- Raptor-side contig-extraction bridge scripts
- a SQLite-backed campaign runner with isolated per-campaign artifacts,
  manifests, status output, and resume support
- read-mapping promotion script for fixed-panel and `quick_test` rechecks
- branch-resolution promotion script for fixed-panel and `quick_test` rechecks
- mapper knobs threaded into standalone polishing, standalone scaffolding, and
  the shared scaffold+polish fast path
- JSON output suitable for external orchestration

Still missing:

- end-to-end quality gate automation for component winners beyond mapper,
  branch resolution, and shared scaffold+polish
- promotion parity for error correction, which still stops at heldout-aware
  local search
- candidate-level multi-fidelity scheduling beyond the mapper grid
- cross-campaign lineage or leaderboard views
- biological quality evidence strong enough to flip later component defaults
  with confidence

## Latest Autonomous Result

Latest checked autonomous component result:

- component: `contig_extraction`
- tested knobs:
  `prefer_high_count_seeds`, `prefer_non_repeat_seeds`,
  `enable_repeat_seed_completion`, `suppress_redundant_contigs`
- current train/val/heldout local winner:
  `prefer_high_count_seeds=true`, `prefer_non_repeat_seeds=true`,
  `enable_repeat_seed_completion=true`, `suppress_redundant_contigs=true`
- eligible candidates after heldout non-regression: `6`
- current local scores:
  train `1.0`, val `1.0`, heldout `1.0`
- train-only throughput winner:
  `prefer_high_count_seeds=true`, `prefer_non_repeat_seeds=false`,
  `enable_repeat_seed_completion=false`, `suppress_redundant_contigs=true`
- baseline heldout before suppression:
  score `0.8691`, `62` observed contigs vs `34` expected
- heldout comparison at `prefer_high_count_seeds=false`:
  score `0.5809`, exact-contig rate `0.1765`, truth-kmer `F1 0.8737`
- current `quick_test` baseline vs candidate:
  output-identical on this dataset (`22` contigs, `19` scaffolds,
  `369` polish corrections), candidate `78.45s` vs baseline `75.75s`
- latest contig-extraction train artifact:
  [latest_result.json](/home/jake/Projects/Raptor/artifacts/autoresearch_raptor/contig_extraction/latest_result.json)
- latest contig-extraction promotion artifact:
  [summary.json](/home/jake/Projects/Raptor/artifacts/autoresearch_raptor/contig_extraction/promotions/contig_extraction_20260313_120134/summary.json)
- promotion verdict: `review`
- repeat-heavy `drosophila` subset result:
  `20092` contigs vs baseline `20179`, but candidate runtime `626.91s` vs
  baseline `329.76s`

Interpretation:

- the contig-extraction bridge is now working end to end
- mapper and branch-resolution defaults are both promoted and live
- shared scaffold+polish already has a guarded bridge and still converges to
  the production default
- error correction now has a real two-sided local gradient, but it still does
  not need a promotion script because the widened local surface still selects
  the baseline `min_count=1, min_trusted_count=4`
- contig extraction now has a real local winner: add
  `suppress_redundant_contigs=true` while keeping the existing seed defaults
- contig extraction now also has a real promotion gate, and it correctly held
  the candidate at `review` instead of flipping the live default
- the current repeat-heavy suite still is not aligned enough with the local
  failure mode, because the candidate run reports `Redundant-contig
  suppression removed 0 contigs`
- the next high-value step is to refine the repeat-heavy suite or add
  phase-level timing so the runtime regression can be attributed before
  spending more promotion cycles on this surface

## Latest Implementation Update

Date: 2026-03-13

Implemented in this pass:

- [run_campaign.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/run_campaign.py)
  now runs resumable campaigns across existing component bridges
- [campaign_db.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/campaign_db.py)
  now persists `campaigns`, `component_runs`, `stage_runs`, and `artifacts`
- [promote_contig_extraction_candidate.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/promote_contig_extraction_candidate.py)
  now compares contig candidates against the baseline on train/val/heldout,
  `quick_test`, and a cached repeat-heavy `drosophila` subset
- `assemble-large` now exposes the full contig-extraction surface:
  `prefer_high_count_seeds`, `prefer_non_repeat_seeds`,
  `enable_repeat_seed_completion`, and `suppress_redundant_contigs`
- campaign-local manifests now exist in
  [default_campaign.json](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/manifests/default_campaign.json)
  and
  [smoke_campaign.json](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/manifests/smoke_campaign.json)
- each campaign now writes isolated artifacts under
  `artifacts/autoresearch_raptor/campaigns/<name>_<timestamp>/` instead of
  overwriting the component-global `latest_result.json`
- the runner now marks stale `running` stages as failed on resume so an
  interrupted unattended pass can be restarted safely
- [read_mapping_train.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/read_mapping_train.py)
  now supports `--search-policy successive-halving` for the heaviest bridge
  grid
- the default campaign manifest now opts `read_mapping` into successive
  halving, and a dedicated smoke manifest now exists in
  [read_mapping_halving_smoke.json](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/manifests/read_mapping_halving_smoke.json)
- the new runbook now exists in
  [docs/autoresearch-campaign-runner.md](/home/jake/Projects/Raptor/docs/autoresearch-campaign-runner.md)

Validation run in this pass:

- `python3 -m py_compile scripts/autoresearch_bridge/bridge_common.py scripts/autoresearch_bridge/campaign_db.py scripts/autoresearch_bridge/run_campaign.py`
- `python3 scripts/autoresearch_bridge/run_campaign.py run --manifest scripts/autoresearch_bridge/manifests/smoke_campaign.json`
- `python3 scripts/autoresearch_bridge/run_campaign.py status --campaign-id 1`
- `python3 scripts/autoresearch_bridge/run_campaign.py resume --campaign-id 1`
- `python3 -m py_compile scripts/autoresearch_bridge/read_mapping_train.py`
- `python3 scripts/autoresearch_bridge/read_mapping_prepare.py --binary target/debug/raptor --output-root /tmp/raptor_mapper_halving_smoke --train-tasks 4 --val-tasks 2 --heldout-tasks 1 --paired`
- `python3 scripts/autoresearch_bridge/read_mapping_train.py --binary target/debug/raptor --tasks-root /tmp/raptor_mapper_halving_smoke --require-heldout --smoke --search-policy successive-halving --halving-factor 3 --final-heldout-candidates 9 --output /tmp/raptor_mapper_halving_result.json`
- `python3 scripts/autoresearch_bridge/run_campaign.py run --manifest scripts/autoresearch_bridge/manifests/read_mapping_halving_smoke.json`

Interpretation:

- the repo can now run unattended multi-component bridge campaigns without
  clobbering global artifacts
- the orchestration layer is stable enough for repeated local use and resume
  after interruption
- the heaviest bridge grid now has real candidate-level multi-fidelity racing:
  the validated mapper smoke run narrowed `81` configs to `27` `val`
  candidates and `12` `heldout` candidates, cutting split evaluations from
  `243` to `120`
- the current limitation is no longer campaign plumbing; it is expanding that
  pattern only where it is justified and building stronger end-to-end gold
  suites
- the best next actions are to harden contig heldout evidence and then teach
  more components to use multi-fidelity search only if their train grids become
  large enough to justify the extra machinery

## Architecture Decisions

Keep these fixed:

- `component-bench` owns task contracts and scoring
- component tasks live under `artifacts/autoresearch_raptor/`
- local search loops only mutate a tiny parameter surface at first
- end-to-end `quick_test` remains the promotion gate

Do not do these yet:

- mutate the whole assembler
- let the optimizer rewrite arbitrary Rust files
- optimize tiny helpers with no biological score
- treat local component score as sufficient for acceptance

## Component Roadmap

### 1. Read Mapping

Status: implemented.

Files:

- [src/cli/component_bench.rs](/home/jake/Projects/Raptor/src/cli/component_bench.rs)
- [docs/read-mapping-component-bench.md](/home/jake/Projects/Raptor/docs/read-mapping-component-bench.md)

Remaining work:

- keep mapper defaults fixed unless a later end-to-end regression points back
  here
- only add more mapper knobs if a future failure mode cannot be explained with
  the current surface
- only return to mapper work if a future regression points back to post-assembly
  mapping

### 2. Branch Resolution

Status: implemented.

Current contract:

- input: branch-choice cases with `current_count` and candidate extensions
- output: chosen branch
- score: exact decision rate under fixed cases
- guardrail: no increase in held-out misassemblies and no `quick_test`
  regression

Files:

- [src/pipeline/large_genome_assembler.rs](/home/jake/Projects/Raptor/src/pipeline/large_genome_assembler.rs)
- [src/cli/component_bench.rs](/home/jake/Projects/Raptor/src/cli/component_bench.rs)
- [docs/branch-resolution-component-bench.md](/home/jake/Projects/Raptor/docs/branch-resolution-component-bench.md)

Remaining work:

- keep the promoted defaults fixed unless a future regression points back to
  branch choice
- decide later whether `support_dominates` needs its own guarded heldout split
  if a real assembly regression appears there

### 3. Error Correction

Status: prepare/train bridge implemented; promotion loop not built.

Needed contract:

- input: cached k-mer count tables and truth-derived trusted targets
- output: corrected count table
- score: precision/recall plus runtime

Files:

- [src/cli/component_bench.rs](/home/jake/Projects/Raptor/src/cli/component_bench.rs)
- [docs/error-correction-component-bench.md](/home/jake/Projects/Raptor/docs/error-correction-component-bench.md)
- [scripts/autoresearch_bridge/error_correction_prepare.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/error_correction_prepare.py)
- [scripts/autoresearch_bridge/error_correction_train.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/error_correction_train.py)
- [docs/autoresearch-error-correction-bridge.md](/home/jake/Projects/Raptor/docs/autoresearch-error-correction-bridge.md)

Remaining work:

- decide whether the local panel needs harder near-trusted heldout cases or an
  ambiguity-aware singleton guard before promotion exists
- only add an error-correction promotion script if a non-trivial local winner
  appears

### 4. Shared Scaffold + Polish Evidence

Status: implemented and guarded; current local search converges to the existing
default.

Needed contract:

- input: contigs, paired reads, truth targets
- output: scaffold links and polish corrections
- score: scaffold correctness, correction precision, runtime

Files:

- [src/cli/component_bench.rs](/home/jake/Projects/Raptor/src/cli/component_bench.rs)
- [docs/scaffold-polish-component-bench.md](/home/jake/Projects/Raptor/docs/scaffold-polish-component-bench.md)
- [scripts/autoresearch_bridge/scaffold_polish_prepare.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/scaffold_polish_prepare.py)
- [scripts/autoresearch_bridge/scaffold_polish_train.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/scaffold_polish_train.py)
- [scripts/autoresearch_bridge/promote_scaffold_polish_candidate.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/promote_scaffold_polish_candidate.py)

### 5. Contig Extraction / Path Selection

Status: prepare/train/promotion implemented; current promotion verdict is
`review`, so the live default stays unchanged.

Needed contract:

- input: cached graph state, branch-support evidence, and truth contigs
- output: extracted contigs
- score: exact contig recovery, truth-kmer F1, contig-count agreement, runtime

Files:

- [src/pipeline/large_genome_assembler.rs](/home/jake/Projects/Raptor/src/pipeline/large_genome_assembler.rs)
- [src/cli/component_bench.rs](/home/jake/Projects/Raptor/src/cli/component_bench.rs)
- [docs/contig-extraction-component-bench.md](/home/jake/Projects/Raptor/docs/contig-extraction-component-bench.md)
- [scripts/autoresearch_bridge/contig_extraction_prepare.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/contig_extraction_prepare.py)
- [scripts/autoresearch_bridge/contig_extraction_train.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/contig_extraction_train.py)
- [scripts/autoresearch_bridge/promote_contig_extraction_candidate.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/promote_contig_extraction_candidate.py)
- [docs/autoresearch-contig-extraction-bridge.md](/home/jake/Projects/Raptor/docs/autoresearch-contig-extraction-bridge.md)

Remaining work:

- keep `suppress_redundant_contigs` as the current local winner, but do not
  flip the live default while the promotion verdict stays `review`
- refine the repeat-heavy end-to-end suite so it triggers explicit
  redundant-contig suppression instead of only a runtime regression
- add phase-level timing to the promotion path so the repeat-heavy slowdown can
  be attributed before another promotion decision

## Raptor-Side Autoresearch Bridge

Because the current external `autoresearch` repo has no plugin interface, the
practical bridge should live inside Raptor first.

Bridge shape:

- `scripts/autoresearch_bridge/read_mapping_prepare.py`
- `scripts/autoresearch_bridge/read_mapping_train.py`
- `scripts/autoresearch_bridge/promote_read_mapping_candidate.py`
- `scripts/autoresearch_bridge/branch_resolution_prepare.py`
- `scripts/autoresearch_bridge/branch_resolution_train.py`
- `scripts/autoresearch_bridge/promote_branch_resolution_candidate.py`
- `scripts/autoresearch_bridge/error_correction_prepare.py`
- `scripts/autoresearch_bridge/error_correction_train.py`
- `scripts/autoresearch_bridge/contig_extraction_prepare.py`
- `scripts/autoresearch_bridge/contig_extraction_train.py`
- `scripts/autoresearch_bridge/promote_contig_extraction_candidate.py`

Rules:

- use the existing `target/debug/raptor` or `target/release/raptor` binary
- call `component-bench ... --json`
- keep train/val task roots fixed
- keep the score definition fixed for a run
- reserve a held-out validation panel

## Mutable Parameter Surfaces

Start narrow.

### Read Mapping

Allowed first:

- minimizer `k`
- minimizer `w`
- minimum primary-hit support threshold
- minimum scaffold-hit support threshold

Allowed next:

- optional cap on repetitive minimizer postings

### Branch Resolution

Allowed first:

- minimum winning read-support margin
- minimum absolute read-support threshold
- repeat fallback priority

### Error Correction

Allowed first:

- explicit `min_count`

Allowed next:

- trusted-threshold multiplier
- ambiguity guard for competing trusted neighbors

### Contig Extraction

Allowed first:

- seed ranking
- repeat-seed fallback priority
- repeat-seed completion pass

Allowed next:

- path termination policy
- duplicate / contained-contig suppression
- branch-support-aware path selection

## Promotion Policy

A component-local winner is not a production winner until it passes:

1. component train score improvement
2. component validation score improvement
3. no `quick_test` wall-clock regression beyond agreed tolerance
4. no component-specific repeat-heavy regression beyond agreed tolerance
5. no obvious assembly-quality regression on `quick_test`

The promotion scripts now also record:

- the default baseline config and metrics
- candidate-minus-baseline deltas
- minimap2-based reference coverage and aligned-query fractions
- a final status of `promote_default`, `review`, or `reject`
- explicit reasons when a candidate is held for review

## Concrete Next Steps

1. Keep `suppress_redundant_contigs` as the current contig-local winner, but
   keep the live default unchanged while the contig promotion verdict remains
   `review`.
2. Refine the contig repeat-heavy suite or add phase-level timing so the
   runtime regression is attributable.
3. Revisit `error_correction` by adding an ambiguity guard or harder
   near-trusted heldout cases if we want a non-trivial winner.
4. Decide later whether patching the external
   `/home/jake/Projects/autoresearch` repo is worth the complexity.

## If The Conversation Is Interrupted

Resume in this order:

1. Read [docs/autonomous-optimization-plan.md](/home/jake/Projects/Raptor/docs/autonomous-optimization-plan.md).
2. Read this file.
3. Read [docs/autoresearch-component-plan.md](/home/jake/Projects/Raptor/docs/autoresearch-component-plan.md).
4. Read [docs/autoresearch-read-mapping-bridge.md](/home/jake/Projects/Raptor/docs/autoresearch-read-mapping-bridge.md).
5. Read [docs/autoresearch-branch-resolution-bridge.md](/home/jake/Projects/Raptor/docs/autoresearch-branch-resolution-bridge.md).
6. Read [docs/autoresearch-error-correction-bridge.md](/home/jake/Projects/Raptor/docs/autoresearch-error-correction-bridge.md).
7. Read [docs/autoresearch-contig-extraction-bridge.md](/home/jake/Projects/Raptor/docs/autoresearch-contig-extraction-bridge.md).
8. Read [docs/scaffold-polish-component-bench.md](/home/jake/Projects/Raptor/docs/scaffold-polish-component-bench.md).
9. Read [docs/contig-extraction-component-bench.md](/home/jake/Projects/Raptor/docs/contig-extraction-component-bench.md).
10. Inspect the latest contig promotion summary under
   `artifacts/autoresearch_raptor/contig_extraction/promotions/`.
11. Inspect the latest mapper promotion summaries under
   `artifacts/autoresearch_raptor/read_mapping/promotions/`.
12. Inspect current CLI surface in [src/cli_main.rs](/home/jake/Projects/Raptor/src/cli_main.rs).
13. Continue with contig stress-suite refinement or error-correction hardening.

Validation commands for this implementation:

```bash
cargo build
python3 -m py_compile \
  scripts/autoresearch_bridge/bridge_common.py \
  scripts/autoresearch_bridge/contig_extraction_prepare.py \
  scripts/autoresearch_bridge/contig_extraction_train.py \
  scripts/autoresearch_bridge/promote_contig_extraction_candidate.py
cargo test cli::component_bench::tests --lib
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all
python3 scripts/autoresearch_bridge/contig_extraction_prepare.py \
  --binary target/debug/raptor \
  --output-root artifacts/autoresearch_raptor/contig_extraction/tasks
python3 scripts/autoresearch_bridge/contig_extraction_train.py \
  --binary target/debug/raptor \
  --tasks-root artifacts/autoresearch_raptor/contig_extraction/tasks \
  --require-heldout
target/debug/raptor component-bench contig-extraction \
  --task artifacts/autoresearch_raptor/contig_extraction/tasks/heldout \
  --disable-prefer-high-count-seeds \
  --json
python3 scripts/autoresearch_bridge/promote_contig_extraction_candidate.py \
  --binary target/release/raptor \
  --tasks-root artifacts/autoresearch_raptor/contig_extraction/tasks \
  --train-result artifacts/autoresearch_raptor/contig_extraction/latest_result.json \
  --threads 8
```

First command to re-run for context:

```bash
python3 scripts/autoresearch_bridge/promote_contig_extraction_candidate.py \
  --binary target/release/raptor \
  --tasks-root artifacts/autoresearch_raptor/contig_extraction/tasks \
  --train-result artifacts/autoresearch_raptor/contig_extraction/latest_result.json \
  --threads 8
```
