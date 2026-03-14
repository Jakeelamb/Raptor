# RepeatMasker Rust Port

This repo now carries a first RepeatMasker Rust port prototype aimed at the
postprocessing path rather than the external search engines.

## Scope

Current prototype:

- vendors the upstream RepeatMasker source under `third_party/repeatmasker/`
- reads RepeatMasker `.out` annotations
- reads RepeatMasker `.cat` raw annotation files plus footer metadata
- rewrites normalized `.out` output
- emits GFF3 and machine-readable TSV
- emits repeat family/class summary TSV and enriched JSON summary with occupied-base accounting
- emits a first RepeatMasker-style `.tbl` report
- applies a lightweight `ProcessRepeats`-style adjudication pass for report outputs by suppressing obvious contained lower-score duplicates and tiny fragments
- parses RepeatMasker `.align` files and can use their alignment-level metrics to break ties for report-level adjudication
- reconstructs conservative exact-family fragment chains from adjudicated annotations and can refine those chains from `.align` block terminals for `.tbl`, JSON, and TSV export
- exports `.align`-derived CIGAR / CAF / block summaries as TSV
- uses source FASTA, when supplied, to compute exact `.tbl`/summary header stats without requiring masked FASTA output
- merges overlapping mask intervals per query sequence
- masks FASTA in `N`, `X`, or lowercase mode
- writes masked FASTA and optional JSON stats
- ships a vendored-fixture parity/runtime harness for the production `.out` path

Entrypoint:

- [repeatmasker_rs.rs](/home/jake/Projects/Raptor/crates/repeatmasker-rs/src/bin/repeatmasker_rs.rs)
- [process_repeats_rs.rs](/home/jake/Projects/Raptor/crates/repeatmasker-rs/src/bin/process_repeats_rs.rs)
- [lib.rs](/home/jake/Projects/Raptor/crates/repeatmasker-rs/src/lib.rs)

## Why This Slice First

RepeatMasker itself is mostly orchestration around external engines like
RMBlast, cross_match, and nhmmer. Rewriting those wrappers alone would not
change the dominant search cost much.

The tractable Rust slice is the Perl postprocessing path:

- `.out` parsing
- `.cat` parsing and footer metadata ingestion
- `.out` normalization / serialization
- GFF3 / tabular export
- repeat class/family summary export
- `.tbl` report generation from raw `.cat` annotations and footer metadata
- lightweight duplicate/fragment adjudication for `.tbl`, summary TSV, JSON stats, and masked FASTA
- optional `.align`-guided adjudication for same-family tie breaking
- conservative fragment-chain reconstruction for exact-family multi-hit elements, with `.align`-guided terminal refinement
- `.align`-derived TSV export with block/CIGAR/CAF encodings
- overlap consolidation
- sequence masking
- summary statistics

That gives us a compatible tool surface we can benchmark immediately and widen
later into more of `ProcessRepeats`.

## Production Status

The current publication-ready scope is the `.out` postprocessing path, not the
full legacy `.cat` adjudication surface.

Current status:

- `process_repeats_rs --annotations *.out --tbl --fasta --masked-output` is the primary compatibility target
- the default vendored-fixture gate matches gold outputs on `4/4` `.out` cases
- those cases cover `small-1` default masking, `small-1 --xsmall`, `small-1 --x`, and `hum-1` default masking
- the latest primary report is written under `artifacts/repeatmasker_fixture_benchmark/`

Explicitly non-final areas:

- old-format `is1` `.out` clipping semantics still differ from the legacy gold output
- old-format `.cat` adjudication remains experimental
- the explicit legacy sweep currently lands at `7/10` cases under `artifacts/repeatmasker_fixture_benchmark_legacy/`

That is enough to position the Rust port as a fast, production-usable
RepeatMasker `.out` postprocessor and masker, while keeping old-format `.cat`
replacement clearly marked as incomplete.

## Crate Packaging

The implementation now lives in the standalone package
[repeatmasker-rs](/home/jake/Projects/Raptor/crates/repeatmasker-rs/Cargo.toml).
Raptor re-exports it for in-repo callers, but the publishable crate source of
truth is the package under `crates/repeatmasker-rs/`.

## Current Limitations

This is not yet a full RepeatMasker replacement.

It does not yet:

- run RMBlast, nhmmer, cross_match, or ABBlast
- parse `.cat` refinement alignments or reconstruct the full upstream fragment/linkage graph
- reproduce the full species-aware `ProcessRepeats` `.tbl` variants
- reproduce taxonomy-aware adjudication from `ProcessRepeats`
- expose both raw and adjudicated annotation streams in the summary TSV or `.tbl` output files
- use `.align` records to rebuild full fragment chains with upstream family-specific heuristics

## Benchmark Gate

Primary publication gate:

```bash
python3 bench/repeatmasker/run_fixture_benchmark.py --repeat-count 3 --build-release
```

Explicit legacy sweep:

```bash
python3 bench/repeatmasker/run_fixture_benchmark.py \
  --repeat-count 1 \
  --include-legacy \
  --report-dir artifacts/repeatmasker_fixture_benchmark_legacy
```

The harness writes JSON and Markdown reports and fails if any selected case
misses the gold `.tbl` or masked FASTA output.

## Example

```bash
cargo run --bin repeatmasker_rs -- \
  --annotations sample.out \
  --fasta sample.fa \
  --output sample.fa.masked \
  --mask lowercase \
  --stats-json sample.mask.stats.json
```

```bash
cargo run --bin process_repeats_rs -- \
  --annotations sample.cat \
  --align sample.align \
  --align-tsv sample.align.tsv \
  --chain-tsv sample.chains.tsv \
  --out sample.normalized.out \
  --gff sample.out.gff \
  --tsv sample.annotations.tsv \
  --summary-tsv sample.summary.tsv \
  --tbl sample.tbl \
  --fasta sample.fa \
  --stats-json sample.process.stats.json
```

To also emit masked FASTA:

```bash
cargo run --bin process_repeats_rs -- \
  --annotations sample.cat \
  --align sample.align \
  --tbl sample.tbl \
  --fasta sample.fa \
  --masked-output sample.fa.masked \
  --xsmall \
  --stats-json sample.process.stats.json
```
