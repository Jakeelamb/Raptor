# Trinity Parity Harness

This directory is the reproducible comparison harness for Raptor versus Trinity.

Current status: scaffold only. It intentionally starts with a tiny truth-known fixture so failures are cheap and visible. This is not the frozen public benchmark panel yet.

## Candidate Panel

`panel.json` defines the current deterministic candidate panel. It is not the
final frozen public Trinity parity panel, but it is the active fast regression
gate across the current truth-known stress cases:

- `tiny_alt_isoform` across insert sizes `110,140,160,180`
- `ambiguous_paralog` at insert size `160`
- `compact_fusion` at insert size `160`
- `antisense_overlap` at insert size `160`
- `partial_antisense_overlap` at insert size `160`
- `high_depth_normalization` at insert size `160`

Run the full candidate panel:

```bash
python3 bench/trinity_parity/run_panel.py
```

The runner writes `target/trinity_parity/candidate_panel/panel_report.json`
with per-fixture commands, pass/fail state, selected-isoform precision/recall/F1,
component counts, graph counts, lengths, elapsed time, resource usage, output
file counts, output byte counts, GPU telemetry, Trinity metrics when requested,
and reciprocal Raptor-vs-Trinity selected-isoform precision/recall/F1 whenever
both FASTAs exist. Resource usage is captured with GNU
`/usr/bin/time -v` when available; otherwise the harness records peak observed
process-group RSS by polling `/proc`. GPU telemetry is captured with
`nvidia-smi` when available.

When Trinity is enabled, `panel.json` can also set
`min_trinity_selected_f1` per fixture. The current candidate panel enforces
Raptor-vs-Trinity selected-output F1 `0.95` for `tiny_alt_isoform`,
`antisense_overlap`, and `partial_antisense_overlap`. `ambiguous_paralog` and
`compact_fusion` remain recorded divergence probes until their policy is
explicitly decided.

The panel also runs Raptor paired-read normalization and, when Trinity is
enabled, records Trinity's `insilico_read_normalization` kept-pair counts. The
small fixtures require exact normalized kept-pair retention. The
`high_depth_normalization` fixture forces read reduction with a Raptor
normalization target of `200`; current measured retention is Raptor `1703/5100`
pairs and Trinity `1697/5100` pairs, so the candidate gate allows kept-pair
fraction delta up to `0.02`. This is a retention-count gate, not finished normalization
parity. The harness also records retained original-pair overlap; current
high-depth overlap is `570` pairs with Jaccard `0.201413`, so retained-read
identity remains an explicit gap. The same fixture now assembles from
Raptor-normalized reads and requires assembly FASTA F1 `0.95` against Trinity;
current Raptor and Trinity outputs both recover one `900` bp transcript with
assembly F1 `1.0`.

If Trinity is installed on `PATH`, the same panel can also run Trinity:

```bash
python3 bench/trinity_parity/run_panel.py --run-trinity
```

If Trinity is installed somewhere else, pass an explicit binary or set
`TRINITY_BIN`:

```bash
python3 bench/trinity_parity/run_panel.py --run-trinity --trinity-bin /path/to/Trinity
TRINITY_BIN=/path/to/Trinity python3 bench/trinity_parity/run_panel.py --run-trinity
```

This repo also includes a Docker wrapper for the upstream Trinity image:

```bash
docker pull trinityrnaseq/trinityrnaseq:2.15.2
python3 bench/trinity_parity/run_panel.py --run-trinity --trinity-bin scripts/trinity_docker.sh
```

Override the image with `TRINITY_DOCKER_IMAGE` if a different Trinity tag is
needed. The wrapper mounts the current working directory at the same path inside
the container, which matches the harness because it writes all inputs and
outputs under this repo.

Use `--require-trinity` when Trinity output is mandatory for the gate.

Current Trinity-backed candidate evidence with
`trinityrnaseq/trinityrnaseq:2.15.2`:

- `tiny_alt_isoform`: Raptor and Trinity both emit `[252,240]`, reciprocal selected F1 `1.0` across inserts `110,140,160,180`.
- `ambiguous_paralog`: Raptor emits `[240,234,240]`; Trinity emits `[234,240]` with truth min coverage `0.3`; reciprocal selected F1 is `0.8`.
- `compact_fusion`: Raptor emits `[216,216]`; Trinity emits `[336]`; reciprocal selected F1 is `0.0`.
- `antisense_overlap` and `partial_antisense_overlap`: Trinity runs with `--SS_lib_type RF` on Trinity-safe renamed RF FASTQs and emits the expected `[264,264]` and `[288,288]` truth-covering transcript pairs.

The panel also runs a GPU-requested Raptor workflow with
`cargo run --features gpu -- trinity --gpu` for each case and compares its
selected isoforms back to the CPU-requested workflow. This exercises the OpenCL
k-mer counter when available, but graph build still runs on CPU and these tiny
fixtures are an output-equivalence gate, not a speedup claim.

The panel also gates single-end `raptor trinity` output for fixtures where
single-end evidence is biologically sufficient. `compact_fusion` is paired-only
because that case is designed to require paired-start evidence to split a compact
overlap.

The panel also gates Trinity-style samples-file input. The current fixture uses
the tab-delimited `condition replicate left right` form and verifies the selected
isoforms match the same biological thresholds as direct paired input. It gates
both one-row files and multi-row files that are merged into workflow-owned
left/right FASTQs before normalization.

The panel also gates Trinity-style comma-separated direct read inputs. The
current fixture passes duplicated left/right FASTQ paths through `--input1` and
`--input2` comma lists, verifies they are materialized into workflow-owned merged
FASTQs, and applies the same selected-isoform thresholds as direct paired input.

The panel also gates the current strand-specific workflow path with
`--SS_lib_type RF`. The fixture uses RF-oriented read pairs, verifies the
workflow materializes strand-oriented FASTQs before assembly, and requires
forward-strand selected-isoform precision/F1. The current panel includes an
`antisense_overlap` fixture with a reverse-complement transcript pair. This is
paired with `partial_antisense_overlap`, where only the middle segment is in
opposite-strand orientation. This is still not final strand-specific parity
because it needs Trinity comparison and broader real-data coverage.

The panel also runs a malformed FASTQ negative check. `raptor trinity` must fail
with a nonzero exit code, avoid writing an assembly FASTA, and report FASTQ
record context for the validation error.

## Tiny Fixture

`run_tiny_fixture.py` generates:

- two related transcript sequences with a shared exon and one alternative exon
- deterministic single-end and paired-end FASTQ reads
- truth FASTA
- truth metadata JSON

It then runs the current Raptor production CLI on paired-end reads with a default non-overlapping insert size of `160` bp and compares the output against both generated truth and the checked-in frozen oracle at `oracles/tiny_alt_isoform.fa`:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py
```

For the current paired-insert stress gate:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py --insert-sweep 110,140,160,180
```

To exercise the current Raptor normalize -> assemble path:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py --normalize-raptor --assemble-normalized
```

To exercise the single-command Raptor Trinity-like workflow:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py --run-raptor-workflow --skip-raptor
```

To exercise single-end input through the same workflow:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py --run-raptor-workflow-single --skip-raptor
```

To exercise Trinity-style samples-file input:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py --run-raptor-workflow-samples-file --skip-raptor
```

To exercise multi-row samples-file merging:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py --run-raptor-workflow-samples-file-multi --skip-raptor
```

To exercise comma-separated direct read-list merging:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py --run-raptor-workflow-comma-lists --skip-raptor
```

To exercise the current strand-specific RF workflow contract:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py --run-raptor-workflow-stranded-rf --skip-raptor
```

To exercise malformed FASTQ rejection:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py --run-malformed-fastq-checks --skip-raptor
```

To run the first harder ambiguous isoform/paralog stress fixture:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py --fixture ambiguous_paralog --run-raptor-workflow --skip-raptor
```

To run the current compact-overlap fusion stress fixture:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py --fixture compact_fusion --run-raptor-workflow --skip-raptor
```

This fixture is a compact-overlap fusion-control gate. It previously exposed one
fused 336 bp contig from two 216 bp truth transcripts; the current workflow
splits that compact overlap with paired-start evidence and recovers both truth
transcripts.

To run the current partial antisense-overlap stress fixture:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py --fixture partial_antisense_overlap --run-raptor-workflow-stranded-rf --skip-raptor
```

If Trinity is installed on `PATH`, the same fixture can run Trinity on the same
paired-end reads:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py --run-trinity --require-trinity
```

Use `--trinity-bin /path/to/Trinity` or `TRINITY_BIN=/path/to/Trinity` when the
executable is not on `PATH`.
`scripts/trinity_docker.sh` can be used here after pulling the Docker image.

To replace the tiny frozen oracle with a freshly captured Trinity output:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py --freeze-trinity-oracle
```

By default the harness fails if minimum best truth coverage is below `0.95`.
It also fails if minimum best oracle coverage is below `0.95`. Selected
component isoforms also need reciprocal truth matches at `0.95` coverage, with
precision at least `0.95` and F1 at least `0.95`. Those precision defaults are
current synthetic-fixture gates, not final public-panel parity thresholds.

Outputs are written under `target/trinity_parity/tiny_alt_isoform/`.
Sweep mode writes one report per insert plus `insert_sweep_report.json`.

The script records:

- exact command lines
- tool versions where available
- transcript length/count metrics
- truth and oracle recovery metrics
- Raptor normalization command and kept-pair metrics when requested
- Raptor `trinity` workflow command, input mode, report path, component clustering mode, component JSON path/count, component graph JSON path/count, component transcript candidate and selected-isoform FASTA metrics, selected-isoform precision/recall/F1 metrics, selected-isoform evidence JSON support metrics, scored isoform candidate JSON selection/rejection metrics, graph node/edge counts, edge read/pair/k-mer support, read k-mer node/edge counts, serialized read k-mer node/edge record counts, reconstructed read k-mer path counts/support, capped read k-mer edge sample count, component assigned read/pair counts, and recovery metrics when requested
- Trinity-style samples-file path and reported sample count when requested
- comma-separated direct input mode and reported input-group count when requested
- Trinity-style `SS_lib_type` value and actual assembly input paths when requested
- malformed FASTQ rejection command, exit status, output absence, and stderr context
- Raptor and Trinity normalized kept-pair counts/fractions and retained-pair overlap when requested
- direct Raptor assembly versus Trinity FASTA precision/recall/F1 when requested
- Trinity paired-end command, exit status, output metrics, and truth recovery when requested
- whether Trinity was available, which executable was resolved, and `Trinity --version` output
- the current Raptor limitation that paired-end evidence is only used as reverse-complemented mate sequence evidence, not yet full Butterfly-style pair path constraints

## Rules

- Do not call this the frozen benchmark panel until Jake approves the dataset list and thresholds.
- Do not claim Trinity parity from this fixture.
- Keep raw generated data and outputs under `target/`, not git.
- Keep scripts deterministic.
