# PLAN

## Goal

Build Raptor into a full-fledged Rust replacement for Trinity with stage-by-stage biological benchmark evidence.

## Current Strategy

Start by building the parity map and benchmark harness before changing core algorithms. No parity claim is allowed until Trinity and Raptor are compared on a frozen benchmark panel.

## Phases

- [x] Create Trinity parity map.
- [ ] Freeze benchmark panel.
- [ ] Build Trinity-vs-Raptor benchmark harness. Initial tiny fixture scaffold exists; full Trinity comparison and frozen panel still missing.
- [ ] Close normalization parity.
- [ ] Close Inchworm-equivalent contig construction parity.
- [ ] Close Chrysalis-equivalent clustering/graph partitioning parity.
- [ ] Close Butterfly-equivalent isoform reconstruction parity.
- [ ] Close paired-end evidence handling parity.
- [ ] Close output/evaluation/reporting parity.
- [ ] Run final full benchmark panel.
- [ ] Update docs and parity claim only after evidence passes.

## Open Decisions

- Exact public RNA-seq benchmark datasets need to be selected and frozen.
- Metric thresholds and tolerance for "comparable biological outputs" need to be encoded after the first benchmark harness draft.
- CUDA backend work is deferred until CPU/OpenCL correctness baselines are stable and comparable.

## Current Findings

- `docs/trinity-parity-map.md` is the current stage inventory.
- Overall Trinity replacement readiness is currently 0 because multiple Trinity-stage equivalents and the frozen comparison artifacts are missing or unproven.
- `bench/trinity_parity/run_tiny_fixture.py` now generates a tiny truth-known alternative-isoform fixture and runs current `raptor assemble`.
- Current tiny fixture result after read-overlap rescue and non-repetitive fixture correction: Raptor exits 0 and emits 2 contigs/transcripts with N50 252 against truth transcript lengths 252 and 240. Mean and minimum best truth coverage are 1.0. Previous baseline was 51 contigs with N50 27. This is Inchworm-direction progress, not parity.
- Normalization has a single-codepath problem: CLI knobs and standalone binaries do not currently converge on one parameterized implementation/default set.
