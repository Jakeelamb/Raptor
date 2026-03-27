# Assembler Campaign Matrix

This document tracks the multi-agent optimization campaign started on
2026-03-26 US/Pacific.

## Baseline

Current promotion baseline comes from the best measured `quick_test` run on
2026-03-25 UTC:

- core `assemble-large`: `18.794s`
- end-to-end `assemble-large --scaffold --polish`: `30.59s`
- core phase timings:
  - `distribute=2.348s`
  - `count=2.836s`
  - `error_correct=7.309s`
  - `graph_clean=0.753s`
  - `graph_analyze=0.567s`
  - `branch_thread=3.475s`
  - `contig_extract=1.301s`
  - `write=0.010s`
- shared scaffold+polish:
  - `index=0.139s`
  - `map=11.485s`
  - `total=11.653s`

Expected stable output on the same dataset:

- 22 raw contigs
- 1,817,070 assembled bases
- N50 `120,054 bp`
- 19 scaffolds
- scaffold N50 `360,080 bp`
- 369 polishing corrections

Primary source of record: `docs/assembler-fast-path-log.md`

## Candidate Lanes

Each lane should be developed and benchmarked independently first.

### Lane A: Shared Mapper

Scope:

- `src/pipeline/scaffolder.rs`
- `src/pipeline/polisher.rs`

Target:

- reduce shared `map` time

Primary KPI:

- end-to-end wall-clock on `quick_test`

### Lane B: Error Correction

Scope:

- `src/pipeline/large_genome_assembler.rs`

Target:

- reduce `error_correct`

Primary KPI:

- core `assemble-large` wall-clock

### Lane C: Branch Evidence / Threading

Scope:

- `src/pipeline/large_genome_assembler.rs`
- possibly `src/io/fastq.rs`

Target:

- reduce or eliminate second read pass
- reduce `branch_thread`

Primary KPI:

- core `assemble-large` wall-clock

### Lane D: Counting Backend

Scope:

- `src/kmer/disk_counting_v2.rs`
- `src/kmer/disk_counting_optimized.rs`

Target:

- reduce `count`

Primary KPI:

- core `assemble-large` wall-clock

### Lane E: Graph Core Rewrite

Scope:

- `src/pipeline/large_genome_assembler.rs`

Target:

- dense unitig graph or sparse/minimizer-space migration

Primary KPI:

- core `assemble-large` wall-clock and peak complexity reduction

## Promotion Gates

Every candidate must pass the following before combination testing:

1. targeted unit/property tests for the changed subsystem
2. relevant stage bench:
   - branch threading: `cargo bench --bench large_genome_stages branch_threading -- --noplot`
   - error correction: `cargo bench --bench large_genome_stages error_correction -- --noplot`
   - parsing/distribution/graph analysis: `cargo bench --bench assembly_hot_path -- --noplot`
3. release `quick_test` run with phase timings captured
4. no regression in the stable output summary

Hard rejection rules:

- worse end-to-end `quick_test` without a compensating measurable win in a
  strategically dominant lane
- lower scaffold N50, fewer expected corrections, or obvious output drift
- statistically significant synthetic regression without a larger real-data win

## Combination Matrix

Only combine candidates with disjoint or low-conflict scopes:

- `A + B`
- `A + D`
- `B + D`
- `A + B + D`

Do not combine early prototypes that both rewrite `large_genome_assembler.rs`
control flow until they have individually cleared `quick_test`.

Combination evaluation order:

1. best single candidate by end-to-end wall-clock
2. best orthogonal second candidate
3. full pairwise matrix for surviving candidates
4. one triple only if both pairs are non-regressive

## Current Campaign State

- lane A is now the dominant accepted win:
  - live default mapper config promoted to `k=15, w=10, min_primary=3, min_scaffold=2`
  - current default `quick_test` frontier is `22.13s` end-to-end
  - shared `map` timing is down to `5.843s`
- lane B error-correction masked-signature prefilter was prototyped and
  rejected after `quick_test` regressed to `39.08s` with
  `error_correct=16.428s`
- lane C branch-read spool was prototyped and rejected from production after
  the stage bench showed the spool path slower than reread across all measured
  fixture sizes
- lane D counting backend did not produce a faster kernel, but the campaign
  kept a dedicated `disk_bucket_count` benchmark and targeted count-path tests
- the synthetic-task mapper winner `k=9, w=4, min_primary=5, min_scaffold=4`
  was rejected on real `quick_test` after regressing to `47.03s`
- campaign tooling is now in place under `scripts/autoresearch_bridge`
- next high-value work is still:
  - a different error-correction acceleration strategy that does not add a large
    per-run index build
  - branch-threading evidence capture that wins on real data, not just on paper
  - a more radical graph-core rewrite rather than more local singleton-rescue
    heuristics
