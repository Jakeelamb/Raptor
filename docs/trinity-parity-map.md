# Trinity Parity Map

Status date: 2026-05-20
Branch: `rescue/gpu-trinity`

This document is the control map for turning Raptor into a credible Rust replacement for Trinity. It is intentionally conservative. A stage is not complete because the repo has similar-sounding code; it is complete only when the behavior is mapped, tested, benchmarked, and biologically comparable to Trinity.

## External Trinity Reference

Primary reference:

- Trinity wiki home: https://github.com/trinityrnaseq/trinityrnaseq/wiki
- Trinity running guide: https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity

Trinity's documented de novo pipeline has these core pieces:

- in silico read normalization, enabled by default in modern Trinity runs
- Inchworm: k-mer-driven construction of transcript sequences, often recovering the dominant isoform and unique portions of alternative isoforms
- Chrysalis: clustering Inchworm contigs into locus-like clusters, constructing complete de Bruijn graph components, and partitioning reads across those disjoint graph components
- Butterfly: tracing graph paths using reads and read-pair support to report full-length alternatively spliced isoforms and separate paralog-derived transcripts
- paired-end handling, strand-specific mode, sample-file support, and optional Jaccard clipping for compact/gene-dense genomes
- final transcript FASTA and downstream evaluation/quantification workflows

## Readiness Scale

- 0: absent or only stubbed
- 1: compiles and has unit tests
- 2: passes representative synthetic fixtures
- 3: matches Trinity on controlled simulated RNA-seq fixtures
- 4: matches Trinity on at least three real public RNA-seq datasets by biological metrics
- 5: matches or beats Trinity on biological metrics and runtime/resource usage across the approved benchmark panel

The current overall readiness is the minimum stage score: **0**. That is not a value judgment; it means several Trinity-required stages and the frozen comparison artifacts do not exist yet.

## Current Stage Inventory

| Trinity responsibility | Raptor evidence | Current status | Readiness | Next required proof |
| --- | --- | --- | --- | --- |
| Input FASTQ handling | `src/io/fastq.rs`, checked streaming readers, gzip-capable callers | Real infrastructure exists; needs transcriptome benchmark fixtures and CLI compatibility review | 1 | Add parity fixtures covering single-end, paired-end, gzip, malformed FASTQ, and Trinity-style input combinations |
| In silico normalization | `src/pipeline/normalize.rs`, `src/kmer/normalize.rs`, `src/bin/normalize_reads.rs`, `src/bin/normalize_paired_reads.rs`, `normalize` CLI, `bench/trinity_parity/run_tiny_fixture.py` | Two-pass CMS/ntHash normalization exists for single and paired reads. `NormalizeConfig` now centralizes k, target coverage, min abundance, max read/pair limit, and GPU request state for the CLI and standalone binaries. The tiny harness can run Raptor normalize -> assemble and record kept-pair metrics. GPU is still a logged CPU fallback, and no Trinity kept-read comparison exists yet. | 1 | Compare kept-read behavior to Trinity normalization on fixtures, then add CPU/GPU output-equivalence if GPU normalization becomes real |
| GPU/OpenCL k-mer acceleration | `src/gpu/kmer_gpu.rs`, `src/bin/count_gpu.rs`, `bench/gpu_kmer_baseline.sh`, `docs/rescue-baseline.md` | OpenCL k-mer smoke and CPU/OpenCL baseline exist. This is useful acceleration infrastructure, not Trinity parity. | 2 | Add CPU-vs-GPU output-equivalence checks on normalization/assembly inputs before using acceleration in any parity claim |
| Inchworm-equivalent k-mer contig construction | `src/graph/assembler.rs`, `src/kmer/*`, `src/pipeline/assemble.rs`, `src/pipeline/large_genome_assembler.rs`, `tests/greedy.rs`, `benches/greedy_assembly.rs` | Greedy/adaptive k-mer assembly and a larger de Bruijn-style genome path exist. They are not yet mapped to Inchworm's transcript behavior: dominant isoform recovery plus unique alternative segments. | 1 | Build controlled transcript fixtures with shared exons/alternative exons and compare contig outputs to Trinity Inchworm-stage outputs or a frozen Trinity oracle |
| Chrysalis-equivalent clustering and graph partitioning | `src/graph/partition.rs`, `src/dist/partition.rs`, `src/graph/builder.rs`, `src/pipeline/large_genome_assembler.rs` weighted graph functions, `raptor_components.json` and `raptor_component_graphs.json` from `raptor trinity` | First component graph artifact exists. Workflow components now use `sequence_or_read_kmer` clustering, so observed read k-mers can define component boundaries when sequence overlap alone is below threshold. Reads/pairs are assigned back to components, graph edges report supporting reads/pairs/observed shared k-mers, and each component graph reports read-derived k-mer node/edge counts plus a capped read k-mer edge sample from assigned reads. This is still not Trinity Chrysalis-style complete de Bruijn graph construction from reads. | 1 | Promote read/k-mer graph construction from a boundary signal, summary, and capped sample into full component de Bruijn graph artifacts, then prove component membership/read assignment against Trinity |
| Butterfly-equivalent isoform reconstruction | `src/pipeline/isoform_processor.rs`, `src/graph/isoform_graph.rs`, `src/graph/isoform_traverse.rs`, `src/graph/transcript.rs`, `src/io/transcript_io.rs`, `tests/isoform_expression_properties.rs`, `benches/isoform_paths.rs` | Directed path traversal and transcript writers exist. Current path support is coverage/link based and does not yet prove read/read-pair-supported full-length isoform recovery or paralog separation. README language calling this "Butterfly-like" should remain provisional until benchmarked. | 1 | Add truth-known alternative-splicing fixtures, read-pair path evidence, paralog/fusion stress cases, and compare isoform precision/recall/F1 to Trinity |
| Paired-end evidence handling | `src/io/fastq.rs`, `normalize_paired`, `assemble --input2`, `assemble-large --input2`, branch/scaffold code in `src/pipeline/large_genome_assembler.rs`, component bench prepare/evaluate paths | Paired reads are parsed and `raptor assemble --input R1 --input2 R2` now feeds R1 plus reverse-complemented R2 sequence evidence through the production path. This proves paired ingestion and mismatch validation, but not full Butterfly-style read-pair path constraints. | 1 | Feed paired evidence into component graph/path selection with tests for insert-size, gap-bridging, and paralog/isoform disambiguation |
| Strand-specific RNA-seq | No clear `SS_lib_type` equivalent found | No mapped strand-specific semantics. | 0 | Add explicit strand-orientation model and fixtures before claiming strand-specific support |
| Jaccard clipping / compact genome fusion control | No mapped transcriptome implementation found | Missing. | 0 | Add compact-genome overlapping-UTR stress fixture and either implement support or document an explicit non-goal before parity threshold is frozen |
| Output transcript FASTA/GFA/GTF/GFF3/counts | `src/io/transcript_io.rs`, `src/io/gfa.rs`, `src/io/gff3.rs`, `src/io/gtf.rs`, `assemble --isoforms --gtf --gff3 --counts-matrix`, `isoform` CLI | Writers exist. Current GTF writer emits one exon per transcript and uses transcript IDs as seqnames, so it is not a faithful transcript annotation model for reference-based evaluation. | 1 | Define output contract for de novo transcript FASTA and evaluation-side reference mappings; fix GTF/GFF3 semantics or restrict them to synthetic truth fixtures |
| Evaluation and reporting | `src/eval/metrics.rs`, `src/eval/gtf_compare.rs`, `gtf-compare`, `eval`, `stats`, `bench/README.md`, `bench/trinity_parity/run_tiny_fixture.py`, `src/pipeline/trinity_workflow.rs` | Length metrics, GTF precision/recall, a tiny paired-end Raptor/Trinity harness, a JSON report from `raptor trinity`, and workflow component-count reporting exist. The harness can freeze a live Trinity output as oracle when Trinity is installed. There is still no frozen public panel, BUSCO/read-representation/full-length recovery integration, or resource capture. | 1 | Expand `bench/trinity_parity/` into exact versions, commands, metrics, runtime/RSS/disk/GPU capture, and reproducible public-panel reports |
| End-to-end Trinity-like CLI | `raptor trinity`, `src/pipeline/trinity_workflow.rs`, `raptor normalize`, `raptor assemble`, `bench/trinity_parity/run_tiny_fixture.py --run-raptor-workflow` | A minimal single-command de novo workflow now accepts paired FASTQ, runs production normalization and assembly, writes final transcript FASTA, assembly sidecars, `raptor_components.json`, `raptor_component_graphs.json`, and `raptor_trinity_report.json`. It passes the tiny fixture, but has not been compared directly to Trinity or the frozen panel. | 1 | Compare the single-command workflow directly to Trinity on the tiny fixture, then expand it across the approved benchmark panel |
| Genome-guided Trinity mode | `assemble-large` is genome assembly oriented, not genome-guided transcriptome assembly | Missing for Trinity parity. It remains optional until de novo parity is credible. | 0 | Defer until de novo stages reach at least readiness 4 |

## Current Evidence Worth Keeping

- Rescue branch and GPU/OpenCL smoke are documented in `docs/rescue-baseline.md`.
- `bench/gpu_kmer_baseline.sh` provides a useful starting pattern for captured command/version/hardware baselines.
- `bench/trinity_parity/run_tiny_fixture.py` now generates a deterministic tiny alternative-isoform fixture and runs current `raptor assemble` through the production CLI.
- Property tests already exist for normalization, graph navigation, graph stats, expression/TPM, metrics, and polishing. These are good fast-loop guards, but they are not biological parity evidence.
- `bench/README.md` describes transcript/isoform evaluation intent, but it references scripts that are not currently present in `bench/` and should be replaced by the new parity harness.

## First Tiny Fixture Result

Command:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py
```

Result file:

- `target/trinity_parity/tiny_alt_isoform/report.json`

Initial rescue-baseline Raptor output on the generated single-end reads:

- exit code: 0
- transcript/contig count: 51
- total output bases: 1364
- N50: 27
- truth transcripts: 252 bp and 240 bp
- Trinity executable on PATH: false

After adding bounded read-overlap rescue for small fragmented transcriptome inputs:

- exit code: 0
- transcript/contig count: 2
- total output bases: 492
- N50: 252
- mean best truth coverage: 1.0
- minimum best truth coverage: 1.0
- mean best frozen-oracle coverage: 1.0
- minimum best frozen-oracle coverage: 1.0
- truth transcripts: 252 bp and 240 bp
- Trinity executable on PATH: false

The harness now has a live Trinity capture path:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py --run-trinity --require-trinity
python3 bench/trinity_parity/run_tiny_fixture.py --freeze-trinity-oracle
```

Raptor normalization can now be included before assembly:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py --normalize-raptor --assemble-normalized
```

On the default insert 160 fixture, this kept 15/15 pairs with matched R1/R2 counts and still recovered [252,240] with truth/oracle minimum coverage 1.0.

The single-command workflow can also be exercised:

```bash
python3 bench/trinity_parity/run_tiny_fixture.py --run-raptor-workflow --skip-raptor
```

On the default insert 160 fixture, `raptor trinity` emits [252,240], N50 252, and truth/oracle minimum coverage 1.0 while writing `raptor_trinity_report.json`.

The workflow also emits a first component artifact:

- `raptor_components.json`
- `raptor_component_graphs.json`
- component count: 1
- component clustering: `sequence_or_read_kmer`
- component graph count: 1
- component 0 contigs: [0, 1]
- component 0 total bases: 492
- component 0 assigned reads: 30
- component 0 assigned pairs: 15
- component graph 0 nodes: 2
- component graph 0 edges: 1
- component graph 0 shared edge length: 91 bp
- component graph 0 shared edge read support: 10
- component graph 0 shared edge pair support: 10
- component graph 0 observed shared k-mers: 133
- component graph 0 read k: 25
- component graph 0 read k-mer nodes: 622
- component graph 0 read k-mer edges: 622
- component graph 0 sampled read k-mer edges: 32

Interpretation: the two related assembled transcripts are grouped into one deterministic component graph on the tiny fixture using `sequence_or_read_kmer` clustering, every workflow read pair is assigned back to that component, the component edge has direct read/pair/k-mer support, and the component exposes a read-derived k-mer graph summary plus a capped edge sample. This is a useful first Chrysalis-shaped graph artifact, but it is not yet Trinity Chrysalis parity because full component de Bruijn graph artifacts and Trinity component membership comparison are still missing.

Interpretation: current `raptor assemble --input R1 --input2 R2` now recovers the two truth transcripts exactly on this tiny non-repetitive two-isoform fixture across insert sweep 110,140,160,180 and matches the checked-in frozen oracle at `bench/trinity_parity/oracles/tiny_alt_isoform.fa`. This is useful normalization/Inchworm/paired-ingestion progress, but it is still not Trinity parity: Trinity is not installed on this machine's `PATH` yet, paired-end evidence is not yet used as full path/scaffold constraints, and broader isoform correctness is not measured.

## Misleading Or Risky Areas

- The README currently presents Raptor as a graph-based RNA-seq assembler with "Butterfly-like traversal." That is acceptable as aspiration only if docs do not imply Trinity replacement status.
- `normalize` CLI and standalone binaries now route through `NormalizeConfig`, but GPU normalization is still an explicit CPU fallback rather than accelerated normalization.
- Normalization parity remains unproven until Raptor kept-read behavior is compared to Trinity on transcriptome fixtures.
- Much of the strongest graph work is in `large_genome_assembler`; it may be reusable, but it is not automatically transcriptome/Trinity-equivalent.
- Current GTF output is too simplified for strong biological claims.

## Immediate Next Work

1. Capture a real Trinity tiny oracle on a machine with Trinity installed, or install Trinity locally.
2. Promote read/k-mer graph construction from a boundary signal, summary, and capped sample into full component de Bruijn graph artifacts, then compare membership/read assignment against Trinity.
3. Compare Raptor normalization kept-read behavior to Trinity normalization on a small transcriptome fixture.
4. Replace stale `bench/README.md` claims about missing scripts with the new harness contract.
5. Only after the fast fixture is stable, decide the frozen public benchmark panel and require approval before changing it.
