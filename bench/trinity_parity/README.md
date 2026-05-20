# Trinity Parity Harness

This directory is the reproducible comparison harness for Raptor versus Trinity.

Current status: scaffold only. It intentionally starts with a tiny truth-known fixture so failures are cheap and visible. This is not the frozen public benchmark panel yet.

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
It also fails if minimum best oracle coverage is below `0.95`.

Outputs are written under `target/trinity_parity/tiny_alt_isoform/`.
Sweep mode writes one report per insert plus `insert_sweep_report.json`.

The script records:

- exact command lines
- tool versions where available
- transcript length/count metrics
- truth and oracle recovery metrics
- Raptor normalization command and kept-pair metrics when requested
- Raptor `trinity` workflow command, report path, component JSON path/count, component graph JSON path/count, graph node/edge counts, edge read/pair/k-mer support, read k-mer node/edge counts, component assigned read/pair counts, and recovery metrics when requested
- Trinity paired-end command, exit status, output metrics, and truth recovery when requested
- whether Trinity was available on `PATH`
- the current Raptor limitation that paired-end evidence is only used as reverse-complemented mate sequence evidence, not yet full Butterfly-style pair path constraints

## Rules

- Do not call this the frozen benchmark panel until Jake approves the dataset list and thresholds.
- Do not claim Trinity parity from this fixture.
- Keep raw generated data and outputs under `target/`, not git.
- Keep scripts deterministic.
