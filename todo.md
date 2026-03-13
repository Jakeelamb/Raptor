# Nightly Autonomous Todo

Updated: 2026-03-13

## Current Focus

- [x] Strengthen read-mapping promotion with an explicit baseline comparison.
- [x] Add a machine-readable promotion verdict instead of manual inspection.
- [x] Validate the new promotion flow on the current best mapper candidate.
- [x] Update handoff docs with the new decision policy and latest result.
- [x] Validate the new reference-aware mapper promotion gate.
- [x] Make the local read-mapping objective sensitive to repeat-driven false
  scaffold support by scoring scaffold specificity.
- [x] Add repeat-only ambiguous decoys to the synthetic mapper panels and
  thread the new knob through the bridge prep script.
- [x] Rebuild the default mapper panels and confirm the rejected mapper now
  exposes worse scaffold specificity locally, even though it still wins the
  scalar score.
- [x] Add an explicit read-mapping `heldout/` stress split with higher ambiguous
  repeat-decoy pressure and require candidates to beat the baseline on held-out
  score or specificity before `quick_test`.
- [x] Wire `heldout_metrics` into `read_mapping_train` result output when
  `heldout/` is present.
- [x] Validate heldout-gated promotion against a fresh candidate/baseline pair
  in smoke mode.

## Next

- [x] Tune the mapper search surface so finalist selection respects heldout
  scaffold-specificity non-regression instead of converging back to the same
  rejected `k=9, w=4, min_scaffold_matches=4` family.
- [x] Validate the new deterministic heldout-aware grid search end to end and
  rerun release promotion on the stronger selector.
- [x] Add pair-link correctness/specificity metrics to the read-mapping
  component bench so the local score can penalize false scaffold-join pressure,
  not just per-read placement errors.
- [x] Validate the new pair-aware mapper metric in Rust and rerun the mapper
  search.
- [x] Promote the new pair-aware mapper winner and confirm it removes the
  previous reference-coverage hard fail.
- [x] Switch promotion reference evaluation to the final polished assembly so
  the gate judges the real output rather than the pre-polish scaffold file.
- [x] Rerun promotion on the current mapper winner under the polished-output
  gate and confirm it clears to `promote_default`.
- [x] Flip the production read-mapping defaults to the promoted
  `k=15, w=4, min_primary_matches=2, min_scaffold_matches=4` config.
- [x] Rebuild and run the live `assemble-large` `quick_test` command with
  defaults to confirm the code path now matches the promotion artifact.
- [x] Update the runbooks with the new promoted mapper baseline.
- [x] Start the next component loop on `branch_resolution`.
- [x] Extend the branch-resolution component bench with scenario-level metrics
  and a heldout stress-profile generator so the next bridge can gate on actual
  failure modes, not only overall exact rate.
- [x] Fix branch-resolution panel summary aggregation so heldout scenario
  metrics survive task folding.
- [x] Upgrade the branch-resolution promotion script to baseline-aware
  train/val/heldout gating plus polished-assembly minimap2 evaluation.
- [x] Validate the new branch bench/profile in Rust and the branch Python
  bridge.
- [x] Run branch-resolution prepare/train/promote end to end in release mode
  and promote the heldout-cleared winner.
- [x] Flip the production branch-resolution defaults to the promoted
  `min_win=1`, `min_margin=2`, `prefer_non_repeat=true` config.
- [x] Rebuild and rerun the live default `assemble-large --scaffold --polish`
  `quick_test` path on the promoted branch defaults.
- [x] Update the runbooks and handoff docs with the promoted branch artifact
  and the new live baseline.
- [x] Start the shared scaffold + polish evidence optimizer harness by adding
  a first `component-bench` task/eval surface around
  `scaffold_and_polish_contigs`.
- [x] Compile and validate the new scaffold+polish harness, and add its docs.
- [x] Build the scaffold+polish prepare/train/promote bridge on top of the new
  component bench, keeping the mutable search surface narrow and baseline-aware.
- Validate the new scaffold+polish bridge with panel generation and heldout
  search smoke, then decide whether to run a full promotion gate.
- Keep using this file after each implementation update: mark the completed
  item, then replace this `Next` section with the immediate follow-on task.

## Next

- [x] Harden the autoresearch bridge subprocess parser so component benches can
  recover JSON summaries even when the Rust CLI emits tracing logs on stdout.
- [x] Run the scaffold+polish prepare/train loop end to end, inspect whether
  the heldout winner is materially different from the baseline, and only then
  spend release `quick_test` time on promotion.
- [x] Tighten the scaffold+polish promotion gate so a default flip requires a
  meaningful quality or runtime win, not just non-regression.
- [x] Run the scaffold+polish promotion gate in release mode, compare the
  promoted `min_scaffold_links` candidate against the live baseline on
  `quick_test`, and decide whether to flip the production default.
- [x] Make the scaffold+polish panel mix threshold-edge support profiles so the
  local objective can distinguish `min_scaffold_links` precision/recall
  tradeoffs instead of only measuring throughput.
- [x] Add sparse-true-link scaffold+polish cases so the local panel can punish
  over-conservative `min_scaffold_links` settings on recall instead of keeping
  `3` and `4` quality-identical.
- [x] Rebuild the scaffold+polish panels under the new mixed-support profiles,
  rerun the heldout-aware search, and confirm the local winner converges back
  to the production default `min_scaffold_links=3`.
- [x] Start the `error_correction` optimizer harness by defining the component
  contract, fixed synthetic task shape, and the first `component-bench`
  prepare/evaluate surface.
- [x] Validate the new `error_correction` component bench with Rust tests and a
  CLI smoke run, then document the first local score gradient.
- [x] Build the Python prepare/train bridge for `error_correction`, keeping the
  mutable search surface narrow around `min_count` before adding a promotion
  script.
- [x] Run the new `error_correction` prepare/train loop end to end, confirm the
  generated panel currently saturates at `min_count=1`, and capture the latest
  local artifact.
- [x] Move directly to `contig_extraction` by building the component bench,
  fixed train/val/heldout panels, and the heldout-aware prepare/train bridge.
- [x] Run the new `contig_extraction` prepare/train loop end to end, confirm
  the heldout-selected winner stays at the live default surface, and verify
  the panel rejects `prefer_high_count_seeds=false`.
- [x] Widen `contig_extraction` into path selection / duplicate suppression by
  adding the `suppress_redundant_contigs` surface and rerunning the heldout
  bridge.
- [x] Thread `suppress_redundant_contigs` into `assemble-large` and compare the
  baseline vs candidate on `quick_test`.
- [x] Add a resumable campaign runner with SQLite state and isolated
  per-campaign artifacts around the existing prepare/train/promote bridges.
- [x] Add a contig-extraction promotion gate with repeat-heavy end-to-end
  evidence and use it to decide whether `suppress_redundant_contigs` should
  flip by default. Current verdict: `review`, so keep the live default
  unchanged.
- [x] Make the contig heldout panel harder so it still distinguishes
  `prefer_non_repeat_seeds` and `enable_repeat_seed_completion` after
  suppression is enabled.
- [ ] Refine the contig repeat-heavy promotion suite so it triggers explicit
  redundant-contig suppression instead of only a runtime regression.
- [ ] Add phase-level timing to the contig promotion path so the repeat-heavy
  slowdown is attributable before another default decision.
- [x] Revisit `error_correction` by adding a second meaningful knob
  (`min_trusted_count`), rerun the heldout-aware bridge, and confirm the
  widened surface still selects the baseline while exposing two-sided local
  regressions.
- [ ] Only spend more time on `error_correction` if an ambiguity guard or a
  harder near-trusted heldout split looks likely to produce a non-trivial
  winner.
- [x] Extend the campaign runner from stage-level orchestration into
  candidate-level multi-fidelity racing for the heaviest bridge grid
  (`read_mapping`).
- [ ] Only extend multi-fidelity racing to another component if its train grid
  grows large enough to pay back the extra complexity.
