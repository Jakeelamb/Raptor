# Autonomous Optimization Plan

Date: 2026-03-13

This is the execution plan for continuing Raptor optimization without waiting
for interactive guidance.

## Objective

Improve Raptor one component at a time with a loop that is:

1. cheap enough to run often
2. narrow enough to avoid overfitting
3. strict enough that local gains do not silently damage full assembly output

## Non-Negotiable Rules

- Never treat a synthetic component score as sufficient for production changes.
- Keep mutable parameter surfaces small and auditable.
- Reuse cached task panels and benchmark artifacts whenever possible.
- Only promote a component-local winner if it survives fixed validation and
  `quick_test` plus any component-specific end-to-end stress suite.
- Prefer removing duplicated work over tuning around it.

## Campaign Orchestration Status

- the repo now has a stage-level campaign runner:
  [run_campaign.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/run_campaign.py)
- campaign state now persists in SQLite through
  [campaign_db.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/campaign_db.py)
- each campaign now writes isolated artifacts under
  `artifacts/autoresearch_raptor/campaigns/<name>_<timestamp>/`
- default and smoke manifests now exist in
  [scripts/autoresearch_bridge/manifests/default_campaign.json](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/manifests/default_campaign.json)
  and
  [scripts/autoresearch_bridge/manifests/smoke_campaign.json](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/manifests/smoke_campaign.json)
- the smoke campaign succeeded locally as campaign `1`, covering
  `error_correction` and `contig_extraction` end to end through prepare/train
  plus `status` and no-op `resume`
- current scope is still stage-level orchestration around the existing
  prepare/train/promote scripts, not candidate-level successive halving inside
  the raw parameter grids

## Current Component Order

Work in this order unless new evidence changes the bottleneck:

1. `read_mapping`
2. `branch_resolution`
3. shared scaffold + polish evidence
4. error correction optimizer harness
5. contig extraction / path selection

Rationale:

- `read_mapping` is still on the hot path for post-assembly work.
- `branch_resolution` is already componentized and promotable.
- scaffold/polish should be optimized after mapper promotion because they now
  share the same evidence stream.

## Autonomous Loop

For each component:

1. Prepare fixed train and validation panels.
2. Search only the approved parameter surface.
3. Persist the latest local winner under
   `artifacts/autoresearch_raptor/<component>/latest_result.json`.
4. Re-run the winner on held-out validation.
5. Promote the winner through `quick_test` and any component-specific
   end-to-end stress suite.
6. Only consider changing defaults after the promotion result is stable.

## Acceptance Gates

Promote a candidate only if all of these hold:

1. train score improves or remains tied with better throughput
2. validation score improves or remains tied with better throughput
3. `quick_test` wall-clock does not regress beyond noise tolerance
4. any component-specific stress suite does not regress beyond its configured
   guardrails
5. `quick_test` core output does not obviously regress

For now, "obvious regression" means any of:

- much worse contig or scaffold count
- materially worse N50
- failure to complete
- unexpected drop in mapped-read counts or correction yield

Read-mapping promotions now emit one of three states:

- `promote_default`: candidate beats baseline and clears all configured gates
- `review`: candidate improved important signals but still changed quality or
  correction behavior enough to require inspection
- `reject`: candidate failed validation or clearly regressed protected metrics

Current read-mapping status:

- the mapper bridge now runs a deterministic heldout-aware grid search across
  the discrete knob surface
- the mapper bridge now also supports candidate-level successive halving for
  the large grid, and the default campaign manifest uses that policy with
  `--halving-factor 4` and `--final-heldout-candidates 32`
- the local mapper objective now includes pair-link correctness/specificity in
  addition to per-read specificity
- promotion now evaluates the final polished assembly, not the pre-polish
  scaffold file
- the promoted mapper default is now `k=15`, `w=4`,
  `min_primary_matches=2`, `min_scaffold_matches=4`
- the live default `assemble-large --scaffold --polish` `quick_test` run now
  matches the promotion artifact: `19` scaffolds, scaffold `N50 360080 bp`,
  `370` polish corrections, about `97s` wall-clock
- the validated halving smoke run reduced mapper split evaluations from the
  full-grid `243` down to `120` (`81 train`, `27 val`, `12 heldout`) while the
  heldout gate still selected the safe baseline on that panel

Current branch-resolution status:

- the branch bridge now uses fixed train/val plus a heldout stress split that
  emphasizes support-floor, support-margin, coverage-closeness, and
  non-repeat-preference cases
- promotion now compares candidate and baseline on train/val/heldout, runs the
  final polished assembly through minimap2, and emits an explicit verdict
- the promoted branch default is now `branch_support_min_win=1`,
  `branch_support_min_margin=2`, `prefer_non_repeat=true`
- the promotion artifact is
  [summary.json](/home/jake/Projects/Raptor/artifacts/autoresearch_raptor/branch_resolution/promotions/branch_resolution_20260311_184053/summary.json)
- the live default `assemble-large --scaffold --polish` `quick_test` run now
  matches the promoted branch artifact: `22` contigs, contig `N50 120054 bp`,
  `19` scaffolds, scaffold `N50 360080 bp`, `369` polish corrections, about
  `98s` wall-clock

Current error-correction status:

- the error-correction bridge now has fixed `train/val/heldout` panel
  generation plus a deterministic local search loop
- the latest local artifact is
  [latest_result.json](/home/jake/Projects/Raptor/artifacts/autoresearch_raptor/error_correction/latest_result.json)
- the mutable surface now includes explicit `min_count` plus
  `min_trusted_count`
- the current local winner remains `min_count=1, min_trusted_count=4`
- train/val/heldout score is still `1.0 / 1.0 / 1.0`
- a lenient heldout comparison at `min_count=1, min_trusted_count=3` drops
  score to `0.7800` by cutting preserved retention to `0.5`
- a strict heldout comparison at `min_count=1, min_trusted_count=6` drops
  score to `0.8114` by cutting correction recall to `0.6586` and trusted exact
  to `0.6573`
- do not add a promotion script for this component until the local surface
  produces a non-trivial winner beyond the current baseline

Current contig-extraction status:

- the contig-extraction bridge now has fixed `train/val/heldout` panel
  generation across branching, repeat-completion, and repeat-fallback cases
- the latest local artifact is
  [latest_result.json](/home/jake/Projects/Raptor/artifacts/autoresearch_raptor/contig_extraction/latest_result.json)
- the heldout-aware local winner now keeps the existing seed defaults and adds
  `suppress_redundant_contigs=true`
- a contig-extraction promotion script now exists at
  [promote_contig_extraction_candidate.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/promote_contig_extraction_candidate.py)
- local train/val/heldout score is now `1.0 / 1.0 / 1.0`
- baseline heldout before suppression was
  `62` observed contigs vs `34` expected, with score `0.8691`
- disabling `prefer_high_count_seeds` on heldout drops exact-contig recovery
  from `1.0` to `0.1765` and score from `0.8691` to `0.5809`
- the latest promotion artifact is
  [summary.json](/home/jake/Projects/Raptor/artifacts/autoresearch_raptor/contig_extraction/promotions/contig_extraction_20260313_120134/summary.json)
- promotion status is `review`: no hard quality regressions, but the
  repeat-heavy runtime gate failed
- `quick_test` is still output-identical to baseline on this dataset: `22`
  contigs, contig `N50 120054 bp`, `19` scaffolds, scaffold `N50 360080 bp`,
  `369` polish corrections
- `quick_test` runtime delta is only `+2.70s` (`78.45s` vs `75.75s`)
- the cached repeat-heavy `drosophila` subset reduces final contigs by `87`
  (`20092` vs `20179`) but regresses runtime by about `297s`
- the repeat-heavy candidate run reported `Redundant-contig suppression removed
  0 contigs`, so the current stress suite still needs refinement before this
  knob should flip by default
- keep the new surface and promotion artifact, but keep the live default
  unchanged until the stress suite exercises a real suppression win

Validation for this pass:

- `cargo fmt --all --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all`
- `python3 -m py_compile scripts/autoresearch_bridge/branch_resolution_prepare.py scripts/autoresearch_bridge/branch_resolution_train.py scripts/autoresearch_bridge/promote_branch_resolution_candidate.py`
- `python3 scripts/autoresearch_bridge/read_mapping_train.py --binary target/release/raptor --tasks-root artifacts/autoresearch_raptor/read_mapping/tasks --require-heldout`
- `python3 scripts/autoresearch_bridge/promote_read_mapping_candidate.py --binary target/release/raptor --tasks-root artifacts/autoresearch_raptor/read_mapping/tasks --train-result artifacts/autoresearch_raptor/read_mapping/latest_result.json`
- `target/release/raptor assemble-large -i bench/genome_assembly/data/quick_test/reads_1.fastq.gz --input2 bench/genome_assembly/data/quick_test/reads_2.fastq.gz -o /tmp/raptor_default_mapper_quick.fa -t 8 --min-count 0 --scaffold --polish --polish-iterations 1 --compress-buckets`
- `python3 scripts/autoresearch_bridge/branch_resolution_prepare.py --binary target/debug/raptor --output-root artifacts/autoresearch_raptor/branch_resolution/tasks`
- `python3 scripts/autoresearch_bridge/branch_resolution_train.py --binary target/release/raptor --tasks-root artifacts/autoresearch_raptor/branch_resolution/tasks --require-heldout`
- `python3 scripts/autoresearch_bridge/promote_branch_resolution_candidate.py --binary target/release/raptor --tasks-root artifacts/autoresearch_raptor/branch_resolution/tasks --train-result artifacts/autoresearch_raptor/branch_resolution/latest_result.json --threads 8`
- `target/release/raptor assemble-large -i bench/genome_assembly/data/quick_test/reads_1.fastq.gz --input2 bench/genome_assembly/data/quick_test/reads_2.fastq.gz -o /tmp/raptor_branch_default_quick.fa -t 8 --min-count 0 --scaffold --polish --polish-iterations 1 --compress-buckets`
- `python3 -m py_compile scripts/autoresearch_bridge/bridge_common.py scripts/autoresearch_bridge/error_correction_prepare.py scripts/autoresearch_bridge/error_correction_train.py`
- `python3 scripts/autoresearch_bridge/error_correction_prepare.py --binary target/debug/raptor --output-root artifacts/autoresearch_raptor/error_correction/tasks`
- `python3 scripts/autoresearch_bridge/error_correction_train.py --binary target/debug/raptor --tasks-root artifacts/autoresearch_raptor/error_correction/tasks --require-heldout`
- `target/debug/raptor component-bench error-correction --task artifacts/autoresearch_raptor/error_correction/tasks/heldout --min-count 1 --min-trusted-count 3 --json`
- `target/debug/raptor component-bench error-correction --task artifacts/autoresearch_raptor/error_correction/tasks/heldout --min-count 1 --min-trusted-count 6 --json`
- `python3 -m py_compile scripts/autoresearch_bridge/bridge_common.py scripts/autoresearch_bridge/contig_extraction_prepare.py scripts/autoresearch_bridge/contig_extraction_train.py scripts/autoresearch_bridge/promote_contig_extraction_candidate.py`
- `python3 scripts/autoresearch_bridge/contig_extraction_prepare.py --binary target/debug/raptor --output-root artifacts/autoresearch_raptor/contig_extraction/tasks`
- `python3 scripts/autoresearch_bridge/contig_extraction_train.py --binary target/debug/raptor --tasks-root artifacts/autoresearch_raptor/contig_extraction/tasks --require-heldout`
- `target/debug/raptor component-bench contig-extraction --task artifacts/autoresearch_raptor/contig_extraction/tasks/heldout --disable-prefer-high-count-seeds --json`
- `target/debug/raptor component-bench contig-extraction --task artifacts/autoresearch_raptor/contig_extraction/tasks/heldout --suppress-redundant-contigs --json`
- `python3 scripts/autoresearch_bridge/promote_contig_extraction_candidate.py --binary target/release/raptor --tasks-root artifacts/autoresearch_raptor/contig_extraction/tasks --train-result artifacts/autoresearch_raptor/contig_extraction/latest_result.json --threads 8`
- `python3 -m py_compile scripts/autoresearch_bridge/bridge_common.py scripts/autoresearch_bridge/campaign_db.py scripts/autoresearch_bridge/run_campaign.py`
- `python3 scripts/autoresearch_bridge/run_campaign.py run --manifest scripts/autoresearch_bridge/manifests/smoke_campaign.json`
- `python3 scripts/autoresearch_bridge/run_campaign.py status --campaign-id 1`
- `python3 scripts/autoresearch_bridge/run_campaign.py resume --campaign-id 1`
- `python3 -m py_compile scripts/autoresearch_bridge/read_mapping_train.py`
- `python3 scripts/autoresearch_bridge/read_mapping_prepare.py --binary target/debug/raptor --output-root /tmp/raptor_mapper_halving_smoke --train-tasks 4 --val-tasks 2 --heldout-tasks 1 --paired`
- `python3 scripts/autoresearch_bridge/read_mapping_train.py --binary target/debug/raptor --tasks-root /tmp/raptor_mapper_halving_smoke --require-heldout --smoke --search-policy successive-halving --halving-factor 3 --final-heldout-candidates 9 --output /tmp/raptor_mapper_halving_result.json`
- `python3 scripts/autoresearch_bridge/run_campaign.py run --manifest scripts/autoresearch_bridge/manifests/read_mapping_halving_smoke.json`

## Current Execution Order

### Phase 1: Productionize Mapper Promotion

Needed:

- thread mapper knobs into live `assemble-large` post-processing
- add a read-mapping promotion script
- record train, validation, and `quick_test` artifacts for mapper winners

Deliverable:

- `read_mapping` can be optimized locally and then evaluated on the real
  assembly command without code edits

### Phase 2: Harden Mapper Panels

Needed:

- add a tougher held-out validation panel with higher ambiguous-repeat pressure
- make finalist selection or promotion require non-regressing held-out score or
  specificity versus baseline
- decide whether scaffold specificity needs a larger explicit weight after the
  held-out split exists

Deliverable:

- mapper search rejects repeat-driven over-scaffolding candidates before
  `quick_test`

### Phase 3: Add Better Promotion Metrics

Needed:

- capture per-phase timings from `quick_test`
- capture peak RSS
- add assembly-quality comparison beyond simple stdout summaries
- tighten the new mapper decision gate with better biological quality evidence

Deliverable:

- promotion decisions are based on better evidence than wall-clock alone

### Phase 4: Extend The Same Pattern

Apply the same structure to:

- shared scaffold + polish evidence
- error correction
- contig extraction

Current next target:

- keep `suppress_redundant_contigs` as the current contig-local winner, but do
  not flip the live default while the promotion verdict stays `review`
- refine the repeat-heavy suite so it triggers explicit redundant-contig
  suppression instead of only showing an end-to-end runtime regression
- add phase-level timing to the contig promotion path so the repeat-heavy
  slowdown can be attributed before spending more promotion cycles on this
  surface
- keep `error_correction` parked until an ambiguity guard or harder
  near-trusted heldout split yields a non-trivial winner
- extend the same candidate-level multi-fidelity pattern beyond
  `read_mapping` only if another component’s train grid becomes large enough to
  justify the added complexity

## Artifact Layout

Keep artifacts here:

```text
artifacts/autoresearch_raptor/
  read_mapping/
    tasks/
    latest_result.json
    promotions/
  branch_resolution/
    tasks/
    latest_result.json
    promotions/
  error_correction/
    tasks/
    latest_result.json
  contig_extraction/
    tasks/
    latest_result.json
    promotions/
```

## Command Skeleton

### Read Mapping

Prepare:

```bash
python3 scripts/autoresearch_bridge/read_mapping_prepare.py \
  --binary target/debug/raptor \
  --output-root artifacts/autoresearch_raptor/read_mapping/tasks \
  --paired
```

Train:

```bash
python3 scripts/autoresearch_bridge/read_mapping_train.py \
  --binary target/debug/raptor \
  --tasks-root artifacts/autoresearch_raptor/read_mapping/tasks
```

Promote:

```bash
python3 scripts/autoresearch_bridge/promote_read_mapping_candidate.py \
  --binary target/release/raptor \
  --tasks-root artifacts/autoresearch_raptor/read_mapping/tasks \
  --train-result artifacts/autoresearch_raptor/read_mapping/latest_result.json
```

### Branch Resolution

Prepare:

```bash
python3 scripts/autoresearch_bridge/branch_resolution_prepare.py \
  --binary target/debug/raptor \
  --output-root artifacts/autoresearch_raptor/branch_resolution/tasks
```

Train:

```bash
python3 scripts/autoresearch_bridge/branch_resolution_train.py \
  --binary target/debug/raptor \
  --tasks-root artifacts/autoresearch_raptor/branch_resolution/tasks
```

Promote:

```bash
python3 scripts/autoresearch_bridge/promote_branch_resolution_candidate.py \
  --binary target/release/raptor \
  --tasks-root artifacts/autoresearch_raptor/branch_resolution/tasks \
  --train-result artifacts/autoresearch_raptor/branch_resolution/latest_result.json
```

### Error Correction

Prepare:

```bash
python3 scripts/autoresearch_bridge/error_correction_prepare.py \
  --binary target/debug/raptor \
  --output-root artifacts/autoresearch_raptor/error_correction/tasks
```

Train:

```bash
python3 scripts/autoresearch_bridge/error_correction_train.py \
  --binary target/debug/raptor \
  --tasks-root artifacts/autoresearch_raptor/error_correction/tasks \
  --require-heldout
```

### Contig Extraction

Prepare:

```bash
python3 scripts/autoresearch_bridge/contig_extraction_prepare.py \
  --binary target/debug/raptor \
  --output-root artifacts/autoresearch_raptor/contig_extraction/tasks
```

Train:

```bash
python3 scripts/autoresearch_bridge/contig_extraction_train.py \
  --binary target/debug/raptor \
  --tasks-root artifacts/autoresearch_raptor/contig_extraction/tasks \
  --require-heldout
```

## Decision Rules While Unattended

- If a local winner fails `quick_test`, keep the artifact but do not change
  production defaults.
- If multiple candidates tie on score, prefer the faster one.
- If a candidate helps only train but not validation, reject it.
- If the local objective saturates, make the panel harder before adding more
  knobs.
- If a bottleneck is caused by duplicated passes, fix architecture before
  micro-optimizing code.
- If a promotion result lands in `review`, keep the artifact and continue
  improving the gate before changing defaults.
- If a promotion result lands in `reject`, improve the objective or held-out
  panels before spending time on more promotion runs for the same surface.
- If a new local metric exposes the failure mode but does not yet change the
  winner, add a harder held-out split before retuning metric weights again.
- If a component default is promoted and the live default `quick_test` matches
  it, move to the next component instead of continuing to retune the promoted
  surface.

## Resume Order If Interrupted

1. Read this file.
2. Read [docs/autoresearch-raptor-handoff.md](/home/jake/Projects/Raptor/docs/autoresearch-raptor-handoff.md).
3. Check the latest promotion summaries under
   `artifacts/autoresearch_raptor/*/promotions/`.
4. Continue with the first incomplete phase above.
