# Assembler Fast-Path Optimization Log

Date: 2026-03-24

This log records the concrete performance changes landed in the large-genome
assembler hot path, the benchmark commands used to evaluate them, and the
measured deltas on the checked-in benchmark harness.

## Changes Landed

### 1. Persistent bucket writers

Files:

- `src/kmer/disk_counting_v2.rs`

Change:

- keep bucket files open across repeated `distribute()` batches
- flush/finalize once before `count_all()`
- report disk usage from tracked bucket counts instead of filesystem metadata

Why:

- the old code reopened every bucket file on every batch
- the assembler feeds the disk counter in many batches, so repeated open/close
  work added avoidable I/O overhead before counting even started

### 2. Sequence-only FASTQ hot-path reader

Files:

- `src/io/fastq.rs`
- `src/pipeline/large_genome_assembler.rs`
- `src/pipeline/streaming_assembly.rs`

Change:

- added checked sequence-only FASTQ callbacks with reusable buffers
- switched large-genome distribution and branch-threading paths to sequence-only
  parsing
- switched streaming assembly ingestion to the same path

Why:

- the previous hot path allocated full `FastqRecord` values even when it only
  needed the sequence
- this keeps strict FASTQ validation while reducing per-record allocation and
  string churn

### 3. Direct graph diagnostics / ambiguous-edge scan

Files:

- `src/pipeline/large_genome_assembler.rs`

Change:

- replaced full weighted-graph materialization in the `assemble-large` hot path
  with a direct scan over adjacency for:
  - node / edge / source / sink diagnostics
  - ambiguous branch-edge collection
- kept the legacy weighted-graph implementation for tests and reference
  benchmarking

Why:

- the production path built a throwaway weighted graph only to compute
  diagnostics and ambiguous edges before branch threading
- the new scan returns the same result without allocating the full intermediate
  graph

### 4. Phase-level assembly timing

Files:

- `src/pipeline/large_genome_assembler.rs`

Change:

- added per-phase timings to `AssemblyStats` and its text rendering:
  - distribute
  - count
  - error correction
  - graph clean
  - graph analyze
  - branch threading
  - contig extraction
  - write output
  - total

Why:

- wall-clock alone was not enough to tell which stage regressed
- the assembler now exposes phase attribution directly in normal runs

### 5. Reproducible hot-path benchmark suite

Files:

- `benches/assembly_hot_path.rs`
- `Cargo.toml`

Change:

- added a Criterion bench that compares:
  - full checked FASTQ parsing vs sequence-only checked parsing
  - legacy reopen-per-batch disk distribution vs persistent writers
  - legacy weighted-graph analysis vs direct adjacency scan

### 6. O(1) canonical-neighbor search for singleton error correction

Files:

- `src/pipeline/large_genome_assembler.rs`

Change:

- replaced per-variant `KmerU64::canonical()` construction in
  `find_trusted_neighbor()` with direct bit-twiddled canonicalization
- compute the original reverse complement once, then derive each mutated
  reverse-complement candidate in O(1)

Why:

- singleton rescue was the dominant core-assembly phase on the refreshed
  `quick_test` run
- the old path recomputed a full reverse complement for every one-base variant
  of every singleton k-mer
- for `k=31`, that meant `93` candidate probes per singleton and each probe
  performed an additional `31`-step reverse-complement loop

### 7. Rolling minimizer mapper and sequence-only scaffold/polish ingestion

Files:

- `src/pipeline/polisher.rs`
- `src/pipeline/scaffolder.rs`

Change:

- replaced the naive minimizer scan with a rolling 2-bit k-mer + monotonic
  queue implementation for `k <= 32`
- kept the legacy nested minimizer scan as a fallback reference for unusual
  parameter ranges
- switched polisher and scaffolder hot paths from full-record FASTQ parsing to
  sequence-only checked callbacks
- added explicit shared scaffold+polish phase timings for:
  - index build
  - mapping
  - scaffold write
  - consensus
  - polished FASTA write

Why:

- the refreshed `quick_test` run showed the shared scaffold/polish tail was
  dominated by read mapping, not by consensus or output
- the old mapper recomputed every candidate k-mer hash in every minimizer
  window, twice per read orientation pair, and the scaffolder/polisher still
  allocated full FASTQ records they never used
- on `quick_test`, that old shared mapping phase consumed `79.511s` by itself

## Benchmark Commands

Primary microbench:

```bash
cargo bench --bench assembly_hot_path -- --noplot
```

Error-correction stage bench:

```bash
cargo bench --bench large_genome_stages error_correction -- --noplot
```

Phase-timing smoke run:

```bash
cargo run --release -- assemble-large \
  -i sample.fastq \
  -o /tmp/raptor_phase_timing.fa \
  -k 11 \
  -c 1 \
  --min-contig 10 \
  --threads 1
```

Validation:

```bash
cargo fmt --all
cargo clippy --all-targets --all-features -- -D warnings
cargo test
```

`quick_test` runtime refresh:

```bash
time -p target/release/raptor assemble-large \
  -i bench/genome_assembly/data/quick_test/reads_1.fastq.gz \
  --input2 bench/genome_assembly/data/quick_test/reads_2.fastq.gz \
  -o /tmp/raptor_quick_test_post_error_opt.fa \
  -t 8 \
  -k 31 \
  --min-count 0 \
  --scaffold \
  --polish \
  --compress-buckets
```

Shared scaffold+polish baseline and optimized runs used the same command above;
the stage split comes from the new timing line emitted by
`scaffold_and_polish_contigs()`.

## Benchmark Results

All numbers below are from the checked-in `assembly_hot_path` bench after the
changes landed. Values shown are the benchmark center estimates from the final
Criterion run.

### FASTQ parsing

| Fixture | Full checked record parse | Sequence-only checked parse | Improvement |
|---------|---------------------------|-----------------------------|-------------|
| 10,000 reads | 1.289 ms | 0.714 ms | 1.81x faster |
| 50,000 reads | 6.324 ms | 3.530 ms | 1.79x faster |

### Disk distribution

| Fixture | Legacy reopen per batch | Persistent writers | Improvement |
|---------|--------------------------|--------------------|-------------|
| 10,000 reads | 20.831 ms | 16.665 ms | 1.25x faster |
| 40,000 reads | 79.844 ms | 60.959 ms | 1.31x faster |

### Graph analysis

| Fixture | Legacy weighted graph | Direct scan | Improvement |
|---------|------------------------|-------------|-------------|
| 128 components | 2.070 ms | 0.565 ms | 3.66x faster |
| 512 components | 11.478 ms | 2.579 ms | 4.45x faster |
| 2048 components | 104.900 ms | 16.467 ms | 6.37x faster |

### Error correction

| Fixture | Before | After | Improvement |
|---------|--------|-------|-------------|
| 4,096 trusted kmers | 5.823 ms | 2.948 ms | 1.98x faster |
| 16,384 trusted kmers | 10.542 ms | 6.153 ms | 1.71x faster |
| 65,536 trusted kmers | 26.603 ms | 15.875 ms | 1.68x faster |

## Phase Timing Example

The phase-timing smoke run now prints a line like:

```text
Phase timings (s): distribute=0.062, count=0.004, error_correct=0.000, graph_clean=0.000, graph_analyze=0.000, branch_thread=0.000, contig_extract=0.000, write=0.000, total=0.066
```

That run used:

- command date/time: 2026-03-25 02:59 UTC
- input: `sample.fastq`
- threads: `1`

This was only a phase-timing smoke check on a tiny input, not a publication
benchmark.

## `quick_test` Refresh

The earlier checked-in benchmark report was stale relative to the current code.
I reran the same `quick_test` dataset before and after the singleton-neighbor
optimization on 2026-03-25 to capture the live phase breakdown.

### Core `assemble-large` phase timings

| Phase | Before | After | Improvement |
|------|--------|-------|-------------|
| distribute | 2.174 s | 1.965 s | 1.11x faster |
| count | 2.757 s | 2.738 s | flat |
| error_correct | 12.008 s | 7.690 s | 1.56x faster |
| graph_clean | 0.717 s | 0.737 s | flat |
| graph_analyze | 0.593 s | 0.632 s | flat |
| branch_thread | 3.544 s | 3.759 s | flat |
| contig_extract | 1.258 s | 1.296 s | flat |
| write | 0.008 s | 0.008 s | flat |
| total | 23.250 s | 19.013 s | 1.22x faster |

### End-to-end pipeline (`--scaffold --polish`)

| Run | Wall-clock |
|-----|------------|
| Before hotspot fix | 94.45 s |
| After hotspot fix | 90.10 s |

### Assembly output stability

The refreshed `quick_test` runs produced the same user-visible assembly summary
before and after the optimization:

- 22 raw contigs
- 1,817,070 assembled bases
- N50 120,054 bp
- 19 scaffolds with scaffold N50 360,080 bp
- 369 polishing corrections

## `quick_test` Shared Scaffold+Polish Refresh

After the singleton-rescue optimization landed, I reran `quick_test` with
shared scaffold+polish timing enabled and then optimized the mapper. The
before/after numbers below are from the same host, dataset, and command line.

### Shared scaffold+polish phase timings

| Phase | Before | After | Improvement |
|------|--------|-------|-------------|
| index | 0.831 s | 0.147 s | 5.65x faster |
| map | 79.511 s | 11.962 s | 6.65x faster |
| scaffold_write | 0.001 s | 0.001 s | flat |
| consensus | 0.020 s | 0.019 s | flat |
| polish_write | 0.001 s | 0.001 s | flat |
| total | 80.371 s | 12.137 s | 6.62x faster |

### Full pipeline (`assemble-large --scaffold --polish`)

| Run | Core `assemble-large` | End-to-end wall-clock |
|-----|-----------------------|-----------------------|
| Before mapping optimization | 18.880 s | 99.41 s |
| After mapping optimization | 19.472 s | 31.76 s |

### Output stability

The optimized run preserved the same observed scaffold/polish output on
`quick_test`:

- 19 scaffolds
- scaffold N50 360,080 bp
- 369 polishing corrections
- 1,817,073 scaffold span

## Correctness Guardrails

The direct graph scan is verified against the legacy weighted-graph path by:

- unit test: `test_direct_graph_scan_matches_weighted_graph_reference`

The new sequence-only FASTQ path is verified by:

- unit tests in `src/io/fastq.rs`
- full workspace test suite

The singleton-neighbor optimization is verified by:

- `cargo test trusted_neighbor -- --nocapture`
- `cargo test error_correct_kmers_matches_serial_reference -- --nocapture`
- `cargo bench --bench large_genome_stages error_correction -- --noplot`
- refreshed `quick_test` before/after runs on the same host and dataset

The rolling minimizer and sequence-only scaffold/polish path is verified by:

- `cargo test fast_minimizer_path_matches_reference -- --nocapture`
- `cargo test scaffold_and_polish_contigs_matches_separate_passes_on_simple_dataset -- --nocapture`
- refreshed `quick_test` shared scaffold+polish before/after runs on the same
  host and dataset

## Next Real Benchmark

The next high-value runtime target is now back inside core assembly. After this
patch, the biggest remaining measured phases on `quick_test` are:

- error correction, at `7.849s`
- branch threading, at `3.897s`
- bucket counting, at `2.493s`

The old shared mapper is no longer the whole-pipeline bottleneck. With shared
scaffold+polish down to `12.137s`, the remaining performance work should focus
on:

- further shrinking error correction
- replacing the second FASTQ read pass for branch threading
- reducing the hash-heavy branch-threading inner loop

## Batch-Packing Experiment (Reverted)

On 2026-03-25 I tried packing short-read batches into one contiguous byte
buffer with range metadata in the large-genome distribution and branch-thread
paths. The goal was to eliminate per-read `Vec<u8>` allocation churn.

### Branch-threading stage bench

`cargo bench --bench large_genome_stages branch_threading -- --noplot`

| Fixture | Result |
|---------|--------|
| 32 | no statistically significant change |
| 128 | regressed to `480.67-565.89 us` (`p < 0.05`) |
| 512 | no statistically significant change |

### `quick_test` regression

Using the same `assemble-large --scaffold --polish` command as earlier runs on
2026-03-25:

| Run | Core `assemble-large` | End-to-end wall-clock |
|-----|-----------------------|-----------------------|
| Pre-experiment baseline | 19.472 s | 31.76 s |
| Contiguous batch experiment | 20.395 s | 32.89 s |

The regression showed up in the core timing line as well:

| Phase | Baseline | Experiment |
|------|----------|------------|
| distribute | 2.256 s | 2.367 s |
| count | 2.493 s | 3.033 s |
| error_correct | 7.849 s | 8.027 s |
| branch_thread | 3.897 s | 4.034 s |

Outcome: reverted. The idea was directionally reasonable, but the measured
result was worse on both the synthetic branch-thread benchmark and the real
dataset.

## Zero-Mark Singleton Cleanup + mmap Bucket Reads

After reverting the batch-packing experiment, I made two small hot-path changes
on 2026-03-25:

- corrected singleton kmers are now zeroed and removed in one `retain()` pass
  instead of millions of individual `AHashMap::remove()` calls
- `disk_counting_v2` now memory-maps bucket files during counting instead of
  doing per-kmer `read_exact()` calls

### Targeted validation

- `cargo test test_error_correct_kmers_matches_serial_reference -- --nocapture`
- `cargo test trusted_neighbor -- --nocapture`
- `cargo test disk_counter -- --nocapture`
- `cargo test count_all_rejects_truncated_bucket_files -- --nocapture`
- `cargo bench --bench large_genome_stages error_correction -- --noplot`

The error-correction stage bench stayed statistically flat on the synthetic
fixture:

| Fixture | Result |
|---------|--------|
| 4,096 trusted kmers | no statistically significant change |
| 16,384 trusted kmers | no statistically significant change |
| 65,536 trusted kmers | no statistically significant change |

### Refreshed `quick_test`

I then reran the same release command before and after the new changes on
2026-03-25:

| Run | Core `assemble-large` | End-to-end wall-clock |
|-----|-----------------------|-----------------------|
| Restored baseline after revert | 19.141 s | 31.01 s |
| Zero-mark + mmap path | 18.794 s | 30.59 s |

### Core phase timings

| Phase | Baseline | After | Improvement |
|------|----------|-------|-------------|
| distribute | 2.335 s | 2.348 s | flat |
| count | 2.828 s | 2.836 s | flat |
| error_correct | 7.859 s | 7.309 s | 1.08x faster |
| graph_clean | 0.741 s | 0.753 s | flat |
| graph_analyze | 0.520 s | 0.567 s | flat |
| branch_thread | 3.414 s | 3.475 s | flat |
| contig_extract | 1.246 s | 1.301 s | flat |
| write | 0.009 s | 0.010 s | flat |
| total | 19.141 s | 18.794 s | 1.02x faster |

### Shared scaffold+polish timings

| Phase | Baseline | After |
|------|----------|-------|
| index | 0.136 s | 0.139 s |
| map | 11.559 s | 11.485 s |
| total | 11.728 s | 11.653 s |

### Output stability

The faster run preserved the same observed assembly result:

- 22 raw contigs
- 1,817,070 assembled bases
- N50 120,054 bp
- 19 scaffolds with scaffold N50 360,080 bp
- 369 polishing corrections

## Multi-Agent Campaign Sweep (2026-03-26)

On 2026-03-26 US/Pacific I ran a parallel optimization campaign with separate
lanes for error correction, branch threading / read spooling, counting, and
evaluation tooling. The result was mixed: the new benchmarking surface is worth
keeping, but the two new production-path candidates were explicitly rejected.

### New evaluation tooling kept

Added:

- `scripts/autoresearch_bridge/profile_assembler_quick_test.py`
- `scripts/autoresearch_bridge/run_patch_matrix.py`
- `docs/assembler-campaign-matrix.md`
- dedicated `disk_bucket_count` Criterion coverage in
  `benches/assembly_hot_path.rs`

These do not materially change the production fast path, but they do make it
far easier to test future candidates against the real `quick_test` gate.

### Rejected: masked-signature singleton rescue

I tried a per-position masked-signature prefilter for singleton correction in
`large_genome_assembler.rs`, with exact neighbor probing retained as the final
verifier.

The synthetic rejection-heavy fixture was not convincing, and the real dataset
rejected it outright. On `quick_test` with that path enabled:

- end-to-end wall-clock regressed from the `2026-03-25` baseline `30.59s` to
  `39.08s`
- core `assemble-large` regressed from `18.794s` to `28.430s`
- `error_correct` exploded from `7.309s` to `16.428s`

Output shape stayed stable, but the runtime regression was too large to accept.
Outcome: reverted from the production path.

Reference run:

- `artifacts/autoresearch_raptor/campaign_runs/20260326/combined_20260326_current/summary.json`

### Rejected: production branch-read spool

I also tried writing a compact branch-read spool during phase-1 distribution and
replaying it for branch threading instead of rereading FASTQ.

The dedicated stage bench in `large_genome_stages` showed the spool path losing
to reread on all measured fixture sizes:

| Fixture | reread | spool | Outcome |
|---------|--------|-------|---------|
| 32 | `92.904-94.259 us` | `106.00-128.30 us` | slower |
| 128 | `330.14-392.69 us` | `590.91-594.65 us` | slower |
| 512 | `1.6043-3.0684 ms` | `2.4091-2.4342 ms` | slower |

On the real `quick_test` run with both experimental lanes enabled,
`branch_thread` was only `3.518s`, which was not enough to justify the large
`error_correct` regression and the extra phase-1 work. Outcome: reverted from
the production path, but the spool module and parity bench were kept for future
experiments.

### Counting lane outcome

The counting lane produced useful direct measurement but no faster kernel yet.
The new `disk_bucket_count` Criterion group measured the attempted radix path as
substantially worse than the current `sort_unstable` path:

| Fixture | `sort_unstable` | radix candidate | Outcome |
|---------|------------------|----------------|---------|
| 10,000 seqs | `2.4327-2.6391 ms` | `13.818-15.176 ms` | rejected |
| 40,000 seqs | `6.9599-7.5132 ms` | `45.422-49.430 ms` | rejected |

Outcome: keep the benchmark and helper coverage, keep the production count
kernel unchanged.

### Restored production path after rejects

After removing the masked-signature and production spool paths, I rebuilt and
reran `quick_test` on 2026-03-26:

| Run | Core `assemble-large` | End-to-end wall-clock |
|-----|-----------------------|-----------------------|
| 2026-03-25 accepted baseline | 18.794 s | 30.59 s |
| 2026-03-26 restored production path | 19.445 s | 30.24 s |

Core phase timings for the restored path:

| Phase | 2026-03-25 baseline | 2026-03-26 restored |
|------|----------------------|---------------------|
| distribute | 2.348 s | 2.255 s |
| count | 2.836 s | 2.828 s |
| error_correct | 7.309 s | 7.724 s |
| graph_clean | 0.753 s | 0.727 s |
| graph_analyze | 0.567 s | 0.609 s |
| branch_thread | 3.475 s | 3.719 s |
| contig_extract | 1.301 s | 1.353 s |
| write | 0.010 s | 0.009 s |
| total | 18.794 s | 19.445 s |

Shared scaffold+polish timings for the restored path:

| Phase | 2026-03-25 baseline | 2026-03-26 restored |
|------|----------------------|---------------------|
| index | 0.139 s | 0.138 s |
| map | 11.485 s | 10.477 s |
| total | 11.653 s | 10.644 s |

The observed assembly output remained stable:

- 22 raw contigs
- 1,817,070 assembled bases
- N50 120,054 bp
- 19 scaffolds with scaffold N50 360,080 bp
- 369 polishing corrections

Important note: the campaign tooling smoke test earlier in the day showed
visible run-to-run noise even with the same binary, so the `30.24s` restored
run should be treated as a promising current measurement, not a definitive
algorithmic gain attributable to the new campaign code alone.

Reference run:

- `artifacts/autoresearch_raptor/campaign_runs/20260326/production_after_rejects_20260326/summary.json`

## Mapper Default Promotion (2026-03-26 Late Session)

After the multi-agent sweep, I stayed on the shared mapper lane and found the
first truly dominant improvement of the night: the old speed-oriented read
mapping defaults were substantially faster on real `quick_test` than the
current production defaults, with only a very small reference-alignment tradeoff.

### Kept code changes

- `polisher.rs` now uses a dense reusable contig/strand scratch buffer and a
  single pass over positional hits to choose both the primary mapping and the
  scaffold-support hit.
- live read-mapping defaults were changed to:
  - `k=15`
  - `w=10`
  - `min_primary_matches=3`
  - `min_scaffold_matches=2`

### Default production `quick_test`

With those defaults promoted into the normal no-flag production path on
2026-03-26:

| Run | Core `assemble-large` | End-to-end wall-clock |
|-----|-----------------------|-----------------------|
| 2026-03-26 pre-promotion default | 19.445 s | 30.24 s |
| 2026-03-26 promoted mapper default | 16.138 s | 22.13 s |

Core phase timings:

| Phase | Pre-promotion default | Promoted default |
|------|------------------------|------------------|
| distribute | 2.255 s | 1.651 s |
| count | 2.828 s | 2.504 s |
| error_correct | 7.724 s | 6.402 s |
| graph_clean | 0.727 s | 0.604 s |
| graph_analyze | 0.609 s | 0.509 s |
| branch_thread | 3.719 s | 3.136 s |
| contig_extract | 1.353 s | 1.174 s |
| write | 0.009 s | 0.007 s |
| total | 19.445 s | 16.138 s |

Shared scaffold+polish timings:

| Phase | Pre-promotion default | Promoted default |
|------|------------------------|------------------|
| index | 0.138 s | 0.069 s |
| map | 10.477 s | 5.843 s |
| total | 10.644 s | 5.934 s |

Observed assembly output for the promoted default:

- 22 raw contigs
- 1,817,070 assembled bases
- N50 120,054 bp
- 19 scaffolds with scaffold N50 360,080 bp
- 520 polishing corrections

Reference run:

- `artifacts/autoresearch_raptor/campaign_runs/20260326/default_after_mapper_default_promotion_20260326/summary.json`

### Reference-eval tradeoff

The read-mapping promotion script also compared the faster `w=10/3/2` default
against the slower `w=4/2/4` mapping configuration on `quick_test` reference
alignment:

- reference coverage fraction stayed the same
- aligned query fraction dropped by only `0.00055`
- matched-base fraction dropped by only `0.00018`

That is a tiny quality hit for a large runtime win, so the faster mapper
configuration was promoted into the live defaults.

Reference promotion summary:

- `artifacts/autoresearch_raptor/read_mapping/promotions/read_mapping_20260326_220903/summary.json`

### Rejected mapper candidate from synthetic search

The synthetic-task winner `k=9, w=4, min_primary=5, min_scaffold=4` failed
hard on real `quick_test`:

- end-to-end wall-clock regressed to `47.03s`
- shared mapping exploded to `26.150s`
- output shape changed to 18 scaffolds with scaffold N50 `1,045,669 bp`

This is a good example of why no synthetic win gets promoted without the real
dataset gate.

Reference run:

- `artifacts/autoresearch_raptor/campaign_runs/20260326/mapper_config_k9_w4_mp5_ms4_20260326/summary.json`
