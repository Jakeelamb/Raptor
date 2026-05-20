# SPEC: Raptor As A Rust Trinity Replacement

## Objective

Build Raptor into a full-fledged Rust replacement for Trinity RNA-seq de novo transcriptome assembly. The goal does not stop at feature resemblance. Raptor must demonstrate comparable biological outputs and benchmark evidence for every major Trinity pipeline stage.

## Product Target

Raptor should become a production-grade transcriptome assembler for Illumina RNA-seq reads with:

- Trinity-like de novo paired-end transcriptome assembly.
- Stage-level equivalents for normalization, Inchworm, Chrysalis, and Butterfly behavior.
- Rust-native performance, safety, observability, and reproducibility.
- Optional GPU acceleration where it improves measured throughput without changing biological results.
- Clear CLI workflows, benchmark artifacts, and documentation that make claims auditable.

## Non-Goals

- Do not preserve genome-assembly features at the expense of transcriptome correctness.
- Do not claim Trinity parity from synthetic smoke tests alone.
- CUDA, OpenCL, and other acceleration paths are allowed when tied to a measured bottleneck, but they must preserve outputs and beat a saved baseline before they are used in parity claims.
- Do not rewrite the whole repository blindly. Work stage by stage, with tests and benchmark artifacts.
- Do not tune metrics by filtering away hard cases.

## Architecture Mapping

Raptor must explicitly map and implement/test equivalents for:

- Input preparation and in silico normalization.
- Inchworm: k-mer driven transcript contig construction, dominant isoform recovery, and unique-region reporting for alternative isoforms.
- Chrysalis: clustering contigs into gene/locus-like components, constructing component de Bruijn graphs, and assigning reads to graph components.
- Butterfly: graph path tracing using read and read-pair support to recover full-length isoforms and separate paralog-derived transcripts.
- Post-assembly outputs: transcript FASTA, graph exports, GTF/GFF3 where appropriate, abundance/count matrix support, and evaluation reports.
- Optional genome-guided mode only after de novo parity is credible.

## Current Repo Context

Important starting files:

- `README.md`
- `docs/rescue-baseline.md`
- `bench/gpu_kmer_baseline.sh`
- `scripts/rescue_smoke.sh`
- `src/cli_main.rs`
- `src/main.rs`
- `src/pipeline/assemble.rs`
- `src/pipeline/normalize.rs`
- `src/pipeline/large_genome_assembler.rs`
- `src/pipeline/isoform_processor.rs`
- `src/graph/*`
- `src/kmer/*`
- `src/gpu/*`
- `tests/*`
- `benches/*`

## Scorecard

Primary score: Trinity replacement readiness, evaluated as the minimum score across pipeline stages. A single weak stage blocks completion.

Stage readiness levels:

- 0: absent or only stubbed.
- 1: compiles and has unit tests.
- 2: passes representative synthetic fixtures.
- 3: matches Trinity on controlled simulated RNA-seq fixtures.
- 4: matches Trinity on at least three real public RNA-seq datasets by biological metrics.
- 5: matches or beats Trinity on biological metrics and runtime/resource usage across the approved benchmark panel.

Completion requires every core stage at level 5:

- normalization
- Inchworm-equivalent contig construction
- Chrysalis-equivalent clustering/graph partitioning
- Butterfly-equivalent isoform reconstruction
- paired-end evidence handling
- output/evaluation/reporting
- end-to-end CLI workflow

## Biological Benchmark Requirements

Raptor must be compared to Trinity on a frozen benchmark panel. The panel must include:

- at least one small synthetic truth-known RNA-seq fixture for fast iteration
- at least one medium simulated transcriptome with known isoforms and expression
- at least three real public Illumina RNA-seq datasets with accepted references or benchmark workflows
- paired-end data
- at least one strand-specific dataset if supported
- a compact-genome or overlapping-UTR stress case for fusion behavior

Required metrics:

- transcript count and length distribution
- N50 and ExN-style or expression-aware continuity where available
- read representation/mapping rate back to assembled transcripts
- BUSCO or equivalent completeness where appropriate
- full-length transcript recovery on truth-known datasets
- isoform precision/recall/F1 where truth is known
- paralog/fusion error indicators
- runtime, peak RSS, disk usage, and GPU usage when enabled

## Feedback Loops

Fast loop:

- `cargo check`
- focused unit/property tests for the stage under active development
- a tiny synthetic RNA-seq fixture that runs Raptor and Trinity or a frozen Trinity output oracle

Medium loop:

- stage-specific benchmark panel for the active stage
- `cargo test --features gpu`
- `./scripts/rescue_smoke.sh`
- `./bench/gpu_kmer_baseline.sh`

Final loop:

- full frozen benchmark panel
- Trinity vs Raptor report with raw commands, versions, hardware, metrics, and output artifacts
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-targets --all-features`

## Done When

The goal is complete only when all of these are true:

- `docs/trinity-parity-map.md` maps Trinity stages to Raptor modules and marks every core stage complete with evidence links.
- `bench/trinity_parity/` contains reproducible scripts for downloading/preparing fixtures, running Trinity, running Raptor, and comparing outputs.
- `docs/trinity-parity-report.md` records the frozen benchmark panel, commands, versions, hardware, metrics, and pass/fail conclusions.
- `bench/trinity_parity/run_public_panel.py` passes the frozen public panel with final Raptor-vs-Trinity reciprocal FASTA F1 at or above 0.9 for each gated dataset.
- Every core stage reaches scorecard level 5.
- Raptor end-to-end outputs are biologically comparable to Trinity across the approved benchmark panel.
- Raptor runtime and peak memory are no worse than Trinity by more than the approved tolerance on the benchmark panel, and at least one major stage is measurably faster or lower-memory.
- GPU acceleration, if used in final claims, produces equivalent biological outputs to CPU and has measured speed/resource benefit.
- `cargo clippy --all-targets --all-features -- -D warnings` passes.
- `cargo test --all-targets --all-features` passes.
- No README or docs claim Trinity replacement status without linking to the parity report.

## Human Approval Gates

Require Jake's explicit approval before:

- changing the benchmark panel after it is frozen
- weakening biological metric thresholds
- deleting major existing pipeline functionality
- declaring parity complete
