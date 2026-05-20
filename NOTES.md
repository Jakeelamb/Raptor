# NOTES

## Chronological Notes

- 2026-05-20 Raptor was restored to `rescue/gpu-trinity`; current head includes OpenCL rescue smoke and CPU-vs-OpenCL k-mer baseline.
- 2026-05-20 Current evidence proves GPU k-mer acceleration and general test health, not Trinity replacement parity.
- 2026-05-20 Trinity parity requires explicit normalization, Inchworm, Chrysalis, Butterfly, paired-end evidence, output/evaluation, and end-to-end workflow gates.
- 2026-05-20 Official Trinity docs describe normalization, Inchworm, Chrysalis, Butterfly, paired-end/strand-specific handling, optional Jaccard clipping, staged execution, and final transcript outputs. Raptor currently has pieces for several of these, but no frozen Trinity-vs-Raptor biological benchmark.
- 2026-05-20 The biggest immediate engineering risk is single-codepath drift in normalization: `raptor normalize`, `normalize_reads`, and `normalize_paired_reads` do not present one parameterized implementation.
- 2026-05-20 First tiny alternative-isoform harness is `bench/trinity_parity/run_tiny_fixture.py`. Current Raptor result is deliberately recorded as a failing proxy: 51 short contigs, N50 27, while truth transcripts are 252 and 240 bp.
- 2026-05-20 Bounded read-overlap rescue plus a non-repetitive tiny fixture moved the fast benchmark to 2 contigs, N50 252, mean best truth coverage 1.0, and min best truth coverage 1.0. This is progress, but the next proof needs Trinity/oracle comparison, paired-end evidence, harder fixtures, and public datasets.
- 2026-05-20 Tiny frozen oracle is checked in at `bench/trinity_parity/oracles/tiny_alt_isoform.fa`. `run_tiny_fixture.py` now gates generated-truth and oracle recovery separately at default minimum best coverage 0.95.
- 2026-05-20 `raptor assemble` now accepts `--input2` for paired-end reads. Current implementation uses R1 plus reverse-complemented R2 as sequence evidence; it is not yet full Trinity/Butterfly-style paired path support.
