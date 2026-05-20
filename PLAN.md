# PLAN

## Goal

Build Raptor into a full-fledged Rust replacement for Trinity with stage-by-stage biological benchmark evidence.

## Current Strategy

Start by building the parity map and benchmark harness before changing core algorithms. No parity claim is allowed until Trinity and Raptor are compared on a frozen benchmark panel.

## Phases

- [x] Create Trinity parity map.
- [ ] Freeze benchmark panel. Current deterministic candidate panel exists in `bench/trinity_parity/panel.json`; final public dataset list and thresholds still need approval.
- [ ] Build Trinity-vs-Raptor benchmark harness. Initial fixture scaffold and candidate-panel runner exist and can capture paired-end Trinity output when Trinity is installed; frozen public panel still missing.
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
- `bench/trinity_parity/run_tiny_fixture.py` now generates a tiny truth-known alternative-isoform fixture, runs current `raptor assemble`, checks both generated truth and the frozen oracle at `bench/trinity_parity/oracles/tiny_alt_isoform.fa`, can run paired-end Trinity with optional oracle freeze when Trinity is installed, can exercise Raptor normalize -> assemble with kept-pair metrics, and can run the single-command `raptor trinity` workflow. `bench/trinity_parity/run_panel.py` now runs the current candidate panel from `bench/trinity_parity/panel.json` and writes one machine-readable panel report with elapsed time, RSS, and output-footprint metrics.
- Current tiny paired-end fixture result after read-overlap rescue and non-repetitive fixture correction: Raptor exits 0 with `assemble --input R1 --input2 R2` across insert sweep 110,140,160,180 and emits 2 contigs/transcripts with N50 252 against truth/oracle transcript lengths 252 and 240 for every insert. Mean and minimum best truth/oracle coverage are 1.0. Previous baseline was 51 contigs with N50 27. This is Inchworm/paired-ingestion progress, not Trinity parity.
- Normalization now has a shared `NormalizeConfig` pipeline path used by `raptor normalize`, `normalize_reads`, and `normalize_paired_reads`; CLI `coverage_target`, `max_reads`, and GPU request state are no longer silently dropped. The tiny harness can now run Raptor normalization before assembly and records kept-pair metrics. This is not Trinity normalization parity until kept-read behavior is compared to Trinity.
- `raptor trinity` now provides a minimal Trinity-like end-to-end CLI over the production normalize and assemble paths, writing final FASTA plus a JSON workflow report. Tiny fixture workflow output is [252,240] with truth/oracle min coverage 1.0. This is not final end-to-end Trinity parity until Trinity is run and the frozen panel exists.
- `raptor trinity` now also emits `raptor_components.json`, `raptor_component_graphs.json`, `raptor_component_paths.fasta`, `raptor_component_selected_isoforms.fasta`, `raptor_component_selected_isoforms.json`, and `raptor_component_isoform_candidates.json`. On the tiny paired-end workflow fixture it reports `sequence_or_read_kmer` component clustering, 1 component containing contigs [0,1] with 492 total bases, 30 assigned reads, 15 assigned pairs, a 2-node/1-edge component graph with a 91 bp edge supported by 10 reads, 10 pairs, 133 observed shared k-mers, plus a read-derived k=25 graph with 622 serialized k-mer node records, 622 serialized k-mer edge records, 8 non-branching graph paths covering 622 edges, and a capped sample of 32 read k-mer edges. The component candidate FASTA contains [252,240,91,90,121,109,121,109,91,90] and recovers truth/oracle min coverage 1.0; the production selected component isoform FASTA contains [252,240] and also recovers truth/oracle min coverage 1.0 with precision/recall/F1 1.0 across insert sweep 110,140,160,180. Selected isoform evidence JSON has 2 records with `component_contig_evidence_score_v2` dense ranks, descending scores, insert-sweep total evidence scores [258,288,221,233], max scores [138,150,115,127], direct read support [50,55,40,42], direct pair support [30,35,25,27], read-kmer path support [12,12,12,12], and max path support [6,6,5,5]. Candidate JSON now proves 10 scored candidates, 2 selected contigs, and 8 rejected read-kmer path fragments for every insert in [110,140,160,180]. This is first read/k-mer-backed Chrysalis/Butterfly-shaped graph evidence, not Trinity parity: candidate scoring still needs fusion stress and Trinity comparison.
- The harness now has `--fixture ambiguous_paralog`. At insert 160 it generates truth transcript lengths [240,234,240], and `raptor trinity` emits selected lengths [240,234,240] with truth min coverage 1.0, selected precision 1.0, recall 1.0, F1 1.0, 0 false positives, 0 false negatives, 1 component, 5 component graph edges, 20 scored isoform candidates, 3 selected contig candidates, and 17 rejected candidates including 16 rejected read-kmer path candidates. This adds real stress beyond the tiny two-isoform case and now suppresses the prior partial false-positive contig.
- The harness now has `--fixture compact_fusion`. At insert 160 it generates two overlapping 216 bp truth transcripts; `raptor trinity` now emits selected lengths [216,216] with selected precision 1.0, recall 1.0, F1 1.0, 0 false positives, 0 false negatives, 1 component, 1 component graph edge, 4 isoform candidates, 2 selected contig candidates, and 2 rejected candidates. The previous 336 bp fused contig is split by paired-start gap evidence.
