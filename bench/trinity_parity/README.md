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

Run the full candidate panel:

```bash
python3 bench/trinity_parity/run_panel.py
```

The runner writes `target/trinity_parity/candidate_panel/panel_report.json`
with per-fixture commands, pass/fail state, selected-isoform precision/recall/F1,
component counts, graph counts, lengths, elapsed time, resource usage, output
file counts, output byte counts, GPU telemetry, and Trinity metrics when
requested. Resource usage is captured with GNU
`/usr/bin/time -v` when available; otherwise the harness records peak observed
process-group RSS by polling `/proc`. GPU telemetry is captured with
`nvidia-smi` when available.

If Trinity is installed on `PATH`, the same panel can also run Trinity:

```bash
python3 bench/trinity_parity/run_panel.py --run-trinity
```

Use `--require-trinity` when Trinity output is mandatory for the gate.

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

If Trinity is installed on `PATH`, the same fixture can run Trinity on the same
paired-end reads:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py --run-trinity --require-trinity
```

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
- Raptor `trinity` workflow command, report path, component clustering mode, component JSON path/count, component graph JSON path/count, component transcript candidate and selected-isoform FASTA metrics, selected-isoform precision/recall/F1 metrics, selected-isoform evidence JSON support metrics, scored isoform candidate JSON selection/rejection metrics, graph node/edge counts, edge read/pair/k-mer support, read k-mer node/edge counts, serialized read k-mer node/edge record counts, reconstructed read k-mer path counts/support, capped read k-mer edge sample count, component assigned read/pair counts, and recovery metrics when requested
- Trinity paired-end command, exit status, output metrics, and truth recovery when requested
- whether Trinity was available on `PATH`
- the current Raptor limitation that paired-end evidence is only used as reverse-complemented mate sequence evidence, not yet full Butterfly-style pair path constraints

## Rules

- Do not call this the frozen benchmark panel until Jake approves the dataset list and thresholds.
- Do not claim Trinity parity from this fixture.
- Keep raw generated data and outputs under `target/`, not git.
- Keep scripts deterministic.
