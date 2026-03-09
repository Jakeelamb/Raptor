# RAPTOR Architecture Roadmap

Date: 2026-03-09

This document describes the staged redesign of RAPTOR from a single large-genome assembler into a multi-engine assembly platform.

## Goal

RAPTOR should not try to beat every assembler with one graph and one heuristic stack.

RAPTOR should become:

- a HiFi-first engine for accurate long reads
- an ONT-first engine for ultra-long noisy reads
- a hybrid engine that resolves repeats and phases haplotypes with multiple evidence types
- a graph-native system that treats assembly graphs and phased paths as first-class outputs
- a heterogeneous compute platform that uses CPU, GPU, and NVMe for the stages they are best at

## Current Direction

The current codebase has:

- Rust as a strong systems-language control plane
- optional GPU kernels for a subset of operations
- a large-genome de Bruijn pipeline with disk-based k-mer counting

The immediate structural weakness is that the large-genome path still loses too much evidence between counting and final contig extraction.

The current greedy contig extractor is also seed-order-sensitive. That is a useful signal: RAPTOR needs a stronger graph model and branch-resolution stack, not just different seed ordering heuristics.

## Staged Plan

### Phase 1: Weighted Unitig Core

Objective: replace ad hoc seed walking with an explicit weighted unitig graph.

Deliverables:

- weighted bidirected graph representation over cleaned k-mers
- maximal non-branching path extraction
- per-edge support tracking
- graph-native diagnostics for indegree, outdegree, unitig count, and dead ends

Status:

- weighted graph infrastructure and diagnostics landed in the large-genome assembler
- the production contig builder remains the adaptive coverage-guided traversal until branch resolution is strong enough to preserve current benchmark quality

### Phase 2: Read-Threaded Branch Resolution

Objective: stop making branch decisions from local coverage alone.

Deliverables:

- paired-end read threading through graph edges
- long-read threading through unitigs
- branch scoring from multiple evidence types
- path confidence scores and ambiguity annotations

Status:

- initial short-read branch threading landed for `assemble-large`
- ambiguous edges now collect directed read support before contig extraction
- on `quick_test`, the first pass improved raw contig N50 substantially while keeping misassemblies flat, but runtime increased and NGA50 still needs work

Success criteria:

- improved NGA50 without increasing QUAST misassemblies
- fewer broken unitigs in repetitive regions

### Phase 3: Modality-Specific Engines

Objective: use the right graph model for each sequencing technology.

Deliverables:

- `raptor hifi`: minimizer-space or sparse/multiplex de Bruijn graph engine
- `raptor ont`: run-length and marker/overlap graph engine
- `raptor hybrid`: HiFi graph plus ultra-long and Hi-C/Pore-C evidence

Success criteria:

- no single generic graph representation forced across all read types
- explicit benchmarks per modality

### Phase 4: Heterogeneous Compute Backends

Objective: move from partial acceleration to deliberate CPU/GPU partitioning.

GPU targets:

- homopolymer compression
- syncmer/minimizer extraction
- k-mer counting and filtering
- seed generation
- chaining
- banded consensus kernels
- graph frontier operations

CPU targets:

- graph surgery
- repeat disentangling
- haplotype decisions
- final consensus orchestration

Backend policy:

- CUDA/HIP for peak throughput
- `wgpu` for portable compute
- common Rust orchestration and buffer interfaces

### Phase 5: Graph-Native Products

Objective: make FASTA a derived artifact, not the only artifact.

Deliverables:

- GFA output with weighted edges and evidence tags
- phased path output
- ambiguity and confidence annotations
- pangenome export path from assemblies

## Benchmark Discipline

Every structural change should be evaluated on:

- `quick_test` as the fast regression gate
- at least one larger public real-world dataset
- SPAdes and one modern long-read/hybrid comparator relevant to the data type

Primary metrics:

- runtime
- peak RSS
- N50 and NGA50
- genome fraction
- misassemblies
- duplication ratio

## Near-Term Execution Order

1. Finish replacing the old contig walk with a richer unitig graph and graph diagnostics.
2. Add read-threaded branch scoring.
3. Add an ONT run-length path and a HiFi graph path as separate engines.
4. Move the regular kernels to a real heterogeneous backend.
5. Promote graph-native outputs and pangenome workflows.
