# PLAN

## Goal

Build Raptor into a full-fledged Rust replacement for Trinity with stage-by-stage biological benchmark evidence.

## Current Strategy

Start by building the parity map and benchmark harness before changing core algorithms. No parity claim is allowed until Trinity and Raptor are compared on a frozen benchmark panel.

## Phases

- [ ] Create Trinity parity map.
- [ ] Freeze benchmark panel.
- [ ] Build Trinity-vs-Raptor benchmark harness.
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

