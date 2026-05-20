# Trinity Parity Harness

This directory is the reproducible comparison harness for Raptor versus Trinity.

Current status: scaffold only. It intentionally starts with a tiny truth-known fixture so failures are cheap and visible. This is not the frozen public benchmark panel yet.

## Tiny Fixture

`run_tiny_fixture.py` generates:

- two related transcript sequences with a shared exon and one alternative exon
- deterministic single-end and paired-end FASTQ reads
- truth FASTA
- truth metadata JSON

It then runs the current Raptor production CLI on the single-end reads:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py
```

By default the harness fails if minimum best truth coverage is below `0.95`.

Outputs are written under `target/trinity_parity/tiny_alt_isoform/`.

The script records:

- exact command lines
- tool versions where available
- transcript length/count metrics
- whether Trinity was available on `PATH`
- the current Raptor limitation that the normal `assemble` CLI accepts one input FASTQ and does not yet consume paired-end evidence

## Rules

- Do not call this the frozen benchmark panel until Jake approves the dataset list and thresholds.
- Do not claim Trinity parity from this fixture.
- Keep raw generated data and outputs under `target/`, not git.
- Keep scripts deterministic.
