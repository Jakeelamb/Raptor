# RepeatMasker Fixture Benchmark

This directory contains the reproducible publication gate for the Rust
RepeatMasker port.

The packaged crate source lives under `crates/repeatmasker-rs/`.

## Primary Gate

Run the production `.out` compatibility and runtime sweep:

```bash
python3 bench/repeatmasker/run_fixture_benchmark.py --repeat-count 3 --build-release
```

This defaults to the production-ready `.out` postprocessing path and currently
targets:

- `small-1` default masking
- `small-1 --xsmall`
- `small-1 --x`
- `hum-1` default masking

Reports are written to `artifacts/repeatmasker_fixture_benchmark/`:

- `latest_report.json`
- `latest_report.md`

The command exits non-zero if any selected case misses the gold `.tbl` or
masked FASTA output.

## Legacy Sweep

Run the explicit legacy compatibility report:

```bash
python3 bench/repeatmasker/run_fixture_benchmark.py \
  --repeat-count 1 \
  --include-legacy \
  --report-dir artifacts/repeatmasker_fixture_benchmark_legacy
```

This includes the remaining old-format `.out` / `.cat` cases so they stay
visible without weakening the main publication claim.
