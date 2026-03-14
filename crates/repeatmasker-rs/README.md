# repeatmasker-rs

`repeatmasker-rs` is a fast Rust postprocessor for RepeatMasker annotations.

Current production scope:

- parse RepeatMasker `.out`, `.cat`, and `.align`
- emit normalized `.out`, GFF3, TSV, summary TSV, JSON, and RepeatMasker-style `.tbl`
- mask FASTA in `N`, `X`, or lowercase mode
- compute exact FASTA-driven header stats
- reconstruct conservative fragment chains and refine them with `.align`

Current publication-ready compatibility target:

- the `.out` postprocessing path
- vendored fixture gate currently matches gold outputs on `4/4` production `.out` cases

Still experimental:

- old-format `.cat` adjudication
- old-format `is1` clipping semantics

## CLI

`process_repeats_rs`:

```bash
process_repeats_rs \
  --annotations sample.out \
  --tbl sample.tbl \
  --fasta sample.fa \
  --masked-output sample.masked.fa \
  --stats-json sample.stats.json
```

`repeatmasker_rs`:

```bash
repeatmasker_rs \
  --annotations sample.out \
  --fasta sample.fa \
  --output sample.masked.fa \
  --mask lowercase
```

## Library

```rust,no_run
use repeatmasker_rs::{load_repeatmasker_catalog_path, write_repeatmasker_tbl_path};

let catalog = load_repeatmasker_catalog_path("sample.out".as_ref())?;
write_repeatmasker_tbl_path("sample.tbl".as_ref(), "sample.fa", &catalog, None)?;
# Ok::<(), repeatmasker_rs::RepeatMaskerError>(())
```

## Benchmark Gate

From the repo root:

```bash
python3 bench/repeatmasker/run_fixture_benchmark.py --repeat-count 3 --build-release
```
