# Raptor x Autoresearch Component Plan

Date: 2026-03-12

## Summary

`autoresearch` is not currently a drop-in optimizer for arbitrary Raptor internals.

Its current loop is built around:

- a fixed cached benchmark
- a small mutable search surface
- a pure function shaped like `assemble_fn(reads) -> contigs`
- a single scalar score for fast repeated evaluation

That matches end-to-end local assembly experiments well.
It does not directly match internal Raptor stages like branch threading, scaffolding, or polishing.

The right adaptation is:

1. keep `autoresearch` for small, repeatable task panels
2. expose Raptor components as stable mini-programs with tiny input/output contracts
3. give each component its own local score
4. keep a separate end-to-end regression gate so locally optimized changes do not damage global quality

This is a "Frankenstein" strategy, but it only works if each component is optimized against the right local objective and guarded by end-to-end checks.

## What Autoresearch Expects Today

From the current local repo:

- benchmark tasks are cached local windows with bounded read subsets
- the mutable surface is small and explicit
- evaluation is deterministic and cheap enough to run repeatedly

Important implication:

- do not start by asking `autoresearch` to mutate the whole Raptor codebase
- first define a small mutable surface for each component
- then give `autoresearch` a task cache and evaluator specific to that component

## Recommended Decomposition

Do not split Raptor into too many tiny pieces.
Some stages are too coupled to optimize independently without creating bad local incentives.

Use these component families instead.

### 1. Error Correction

Candidate scope:

- singleton rescue
- trusted-neighbor policy
- threshold policy
- batching / sparsification strategy

Current anchor:

- `error_correct_kmers` in [src/pipeline/large_genome_assembler.rs](/home/jake/Projects/Raptor/src/pipeline/large_genome_assembler.rs#L1261)

Task input:

- cached k-mer count tables
- small read-derived truth panels

Task output:

- corrected k-mer count table

Primary local score:

- downstream filtered k-mer precision/recall against truth-derived trusted set
- runtime
- memory

Guardrail:

- no drop in downstream contig quality on fixed local tasks

### 2. Branch Evidence + Local Traversal

Treat these together, not separately.

Reason:

- read threading and branch choice are tightly coupled
- optimizing thread support without branch choice can produce useless local wins

Candidate scope:

- ambiguous-edge selection
- support thresholds
- branch tie-break logic
- evidence aggregation policy

Current anchors:

- branch support collection in [src/pipeline/large_genome_assembler.rs](/home/jake/Projects/Raptor/src/pipeline/large_genome_assembler.rs#L2130)
- branch choice in [src/pipeline/large_genome_assembler.rs](/home/jake/Projects/Raptor/src/pipeline/large_genome_assembler.rs#L3127)

Task input:

- cached local graph fixture
- ambiguous-edge set
- reads overlapping the local branch neighborhood

Task output:

- chosen branch path or local contig extension

Primary local score:

- exact path recovery on local repeat windows
- edge-support precision/recall
- runtime

Guardrail:

- no increase in misassemblies on held-out end-to-end tasks

### 3. Contig Extraction

Candidate scope:

- seed ranking
- traversal termination
- repeat fallback policy
- unitig/path extraction logic

Current first mutable surface:

- `prefer_high_count_seeds`
- `prefer_non_repeat_seeds`
- `enable_repeat_seed_completion`
- `suppress_redundant_contigs`

Current anchor:

- `build_contigs_from_graph` in [src/pipeline/large_genome_assembler.rs](/home/jake/Projects/Raptor/src/pipeline/large_genome_assembler.rs#L1916)

Task input:

- cached cleaned local graph
- branch support map
- local truth sequence

Task output:

- contigs

Primary local score:

- local NGA50-like score
- truth k-mer F1
- contig count penalty

Guardrail:

- deterministic output under insertion-order perturbation

### 4. Shared Read Mapping

This is now a real independent target.

Candidate scope:

- minimizer parameters
- hit aggregation
- orientation handling
- chaining / best-hit choice
- batching / parallelization policy

Current anchor:

- `map_read_with_scaffold_support` in [src/pipeline/polisher.rs](/home/jake/Projects/Raptor/src/pipeline/polisher.rs#L185)

Task input:

- contigs
- reads
- truth alignments or synthetic planted mappings

Task output:

- best positional hit
- best scaffolding hit

Primary local score:

- mapping accuracy
- mapped-read rate
- runtime

Guardrail:

- scaffold statistics and polishing corrections remain stable on fixed tasks

### 5. Scaffolding + Polishing Postprocess

Treat these as one component family when using paired reads.

Reason:

- they now share the same evidence stream
- optimizing them separately duplicates work and misaligns the objective

Current anchor:

- `scaffold_and_polish_contigs` in [src/pipeline/scaffolder.rs](/home/jake/Projects/Raptor/src/pipeline/scaffolder.rs#L549)

Task input:

- contigs
- paired reads
- local truth or truth-derived scaffold targets

Task output:

- scaffolds
- polished contigs

Primary local score:

- scaffold correctness
- polished contig identity
- runtime

Guardrail:

- no regression in scaffold N50 / correction quality on held-out tasks

## What Not To Split Yet

Do not make these standalone autoresearch targets first:

- disk bucketing / low-level FASTQ I/O
- tiny graph-cleaning primitives in isolation
- individual helper functions with no direct biological score

Those are better handled with Criterion benches, profiling, and manual engineering.

## Benchmark Shape For Each Component

Each component should expose the same five things:

1. fixed cached fixtures
2. a tiny mutable parameter/program surface
3. a pure input/output function
4. a scalar local score
5. an end-to-end guardrail set

If any component cannot be described this way, it is not ready for autoresearch.

## Proposed Task Layout

Create a parallel benchmark root for Raptor component tasks:

```text
artifacts/autoresearch_raptor/
  error_correction/
  branch_resolution/
  contig_extraction/
  read_mapping/
  scaffold_polish/
```

Each task directory should contain only what the component needs:

- `input.*`
- `truth.*`
- `metadata.json`

Avoid full pipeline context unless the component genuinely needs it.

## Mutable Surface Policy

Do not expose all of Raptor to the agent.

For each component, choose one of these mutation modes:

- parameter search only
- constrained Rust function body mutation
- strategy selection among a few implementations

Preferred order:

1. parameters
2. strategy switches
3. narrow function mutation

Avoid broad code mutation until the component benchmark is stable and reproducible.

## Scoring Policy

Every component score should be a weighted combination of:

- correctness
- runtime
- memory or allocation pressure when relevant

Use hard rejects for:

- nondeterminism
- invalid output
- output violating interface invariants

Do not let the agent trade correctness for speed unless the top-level benchmark explicitly allows it.

## Infrastructure Needed In Raptor

Before serious autoresearch integration, Raptor should add:

- cached fixture loaders for each target component
- per-component benchmark CLI entrypoints
- machine-readable output for scores and timings
- deterministic seed control
- explicit phase checkpointing

Minimum useful interface:

```text
raptor component-bench <component> --task <dir> --json
```

The command should:

- load one component task
- run exactly one component implementation
- emit metrics as JSON
- avoid unrelated pipeline work

## Recommended Execution Order

1. Start with `read_mapping`.
   Reason: strong local objective, current real bottleneck, shared by scaffolding and polishing.
2. Then `branch_resolution`.
   Reason: most likely quality limiter in repetitive regions.
3. Then `error_correction`.
   Reason: easy to cache and cheap to evaluate.
4. Then `scaffold_polish`.
   Reason: now has a unified evidence path.
5. Last, `contig_extraction`.
   Reason: highest leverage, but easiest place to overfit a local score.

## Immediate Next Steps

1. Add a `component-bench` CLI for `read_mapping`.
2. Define a tiny cached task format for mapping accuracy + speed.
3. Add one held-out end-to-end regression suite that every component mutation must pass.
4. Only then adapt `autoresearch` to drive that component.

## Bottom Line

Yes, Raptor should be decomposed for autoresearch.

But not into arbitrary source files and not all at once.

The right unit is:

- biologically meaningful
- small enough to iterate fast
- large enough that the local score still matters

That means:

- not "optimize the whole assembler"
- not "optimize every helper"
- optimize a few evidence-bearing component families with hard end-to-end guardrails
