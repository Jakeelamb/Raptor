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
| In silico normalization | `src/pipeline/normalize.rs`, `src/kmer/normalize.rs`, `src/bin/normalize_reads.rs`, `src/bin/normalize_paired_reads.rs`, `normalize` CLI | Two-pass CMS/ntHash normalization exists for single and paired reads. CLI exposes knobs, but `pipeline::normalize` currently hardcodes k=25, target=50, min_abundance=2 and ignores `_use_gpu`; standalone binaries use different defaults. | 1 | Collapse to one production normalization path, wire CLI parameters into the implementation, then compare kept-read behavior to Trinity normalization on fixtures |
| GPU/OpenCL k-mer acceleration | `src/gpu/kmer_gpu.rs`, `src/bin/count_gpu.rs`, `bench/gpu_kmer_baseline.sh`, `docs/rescue-baseline.md` | OpenCL k-mer smoke and CPU/OpenCL baseline exist. This is useful acceleration infrastructure, not Trinity parity. | 2 | Add CPU-vs-GPU output-equivalence checks on normalization/assembly inputs before using acceleration in any parity claim |
| Inchworm-equivalent k-mer contig construction | `src/graph/assembler.rs`, `src/kmer/*`, `src/pipeline/assemble.rs`, `src/pipeline/large_genome_assembler.rs`, `tests/greedy.rs`, `benches/greedy_assembly.rs` | Greedy/adaptive k-mer assembly and a larger de Bruijn-style genome path exist. They are not yet mapped to Inchworm's transcript behavior: dominant isoform recovery plus unique alternative segments. | 1 | Build controlled transcript fixtures with shared exons/alternative exons and compare contig outputs to Trinity Inchworm-stage outputs or a frozen Trinity oracle |
| Chrysalis-equivalent clustering and graph partitioning | `src/graph/partition.rs`, `src/dist/partition.rs`, `src/graph/builder.rs`, `src/pipeline/large_genome_assembler.rs` weighted graph functions | Generic graph and partitioning pieces exist, plus large-genome unitig graph logic. No evidence yet that Inchworm-like contigs are clustered into locus-level transcript graph components or that reads are partitioned among those disjoint components. | 0 | Define `RaptorComponent` semantics, emit component graph artifacts, and prove component membership/read assignment against Trinity on small fixtures |
| Butterfly-equivalent isoform reconstruction | `src/pipeline/isoform_processor.rs`, `src/graph/isoform_graph.rs`, `src/graph/isoform_traverse.rs`, `src/graph/transcript.rs`, `src/io/transcript_io.rs`, `tests/isoform_expression_properties.rs`, `benches/isoform_paths.rs` | Directed path traversal and transcript writers exist. Current path support is coverage/link based and does not yet prove read/read-pair-supported full-length isoform recovery or paralog separation. README language calling this "Butterfly-like" should remain provisional until benchmarked. | 1 | Add truth-known alternative-splicing fixtures, read-pair path evidence, paralog/fusion stress cases, and compare isoform precision/recall/F1 to Trinity |
| Paired-end evidence handling | `src/io/fastq.rs`, `normalize_paired`, `assemble-large --input2`, branch/scaffold code in `src/pipeline/large_genome_assembler.rs`, component bench prepare/evaluate paths | Paired reads are parsed and used in some normalization/genome branch/scaffold paths. The normal `assemble` transcriptome CLI accepts only one input file and Butterfly-equivalent isoform traversal is not yet driven by read-pair paths. | 1 | Add paired-end transcriptome CLI path and feed paired evidence into component graph/path selection with tests for insert-size and mismatched-pair behavior |
| Strand-specific RNA-seq | No clear `SS_lib_type` equivalent found | No mapped strand-specific semantics. | 0 | Add explicit strand-orientation model and fixtures before claiming strand-specific support |
| Jaccard clipping / compact genome fusion control | No mapped transcriptome implementation found | Missing. | 0 | Add compact-genome overlapping-UTR stress fixture and either implement support or document an explicit non-goal before parity threshold is frozen |
| Output transcript FASTA/GFA/GTF/GFF3/counts | `src/io/transcript_io.rs`, `src/io/gfa.rs`, `src/io/gff3.rs`, `src/io/gtf.rs`, `assemble --isoforms --gtf --gff3 --counts-matrix`, `isoform` CLI | Writers exist. Current GTF writer emits one exon per transcript and uses transcript IDs as seqnames, so it is not a faithful transcript annotation model for reference-based evaluation. | 1 | Define output contract for de novo transcript FASTA and evaluation-side reference mappings; fix GTF/GFF3 semantics or restrict them to synthetic truth fixtures |
| Evaluation and reporting | `src/eval/metrics.rs`, `src/eval/gtf_compare.rs`, `gtf-compare`, `eval`, `stats`, `bench/README.md` | Length metrics and GTF precision/recall exist, but there is no Trinity-vs-Raptor harness, no frozen dataset panel, no BUSCO/read-representation/full-length recovery integration, and no resource capture. | 1 | Create `bench/trinity_parity/` harness with exact versions, commands, metrics, runtime/RSS/disk/GPU capture, and reproducible reports |
| End-to-end Trinity-like CLI | `raptor normalize`, `raptor assemble`, `raptor isoform`, `raptor eval` | Pieces are exposed as separate commands. There is no single Trinity-equivalent workflow that accepts paired FASTQ and produces final de novo transcript FASTA plus report. | 0 | Add or script an end-to-end workflow over the same production code paths and compare directly to a Trinity command |
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
- truth transcripts: 252 bp and 240 bp
- Trinity executable on PATH: false

Interpretation: current `raptor assemble` now recovers the two truth transcripts exactly on this tiny non-repetitive two-isoform fixture. This is useful Inchworm-direction progress, but it is still not Trinity parity: Trinity is not yet run or frozen as an oracle, paired-end evidence is not consumed by this path, and broader isoform correctness is not measured.

## Misleading Or Risky Areas

- The README currently presents Raptor as a graph-based RNA-seq assembler with "Butterfly-like traversal." That is acceptable as aspiration only if docs do not imply Trinity replacement status.
- `normalize` CLI exposes `coverage_target`, `max_reads`, and `gpu`, but the current `pipeline::normalize` implementation does not honor those knobs. That is a correctness and operator-trust issue.
- There are separate normalization binaries and pipeline functions with different defaults. This violates the single-codepath requirement for a Trinity replacement.
- Much of the strongest graph work is in `large_genome_assembler`; it may be reusable, but it is not automatically transcriptome/Trinity-equivalent.
- Current GTF output is too simplified for strong biological claims.

## Immediate Next Work

1. Extend `bench/trinity_parity/run_tiny_fixture.py` to compare against Trinity when Trinity is installed or a frozen Trinity oracle is supplied.
2. Replace stale `bench/README.md` claims about missing scripts with the new harness contract.
3. Fix normalization single-codepath drift: one implementation should own k, target coverage, min abundance, paired/single behavior, and GPU/no-GPU semantics.
4. Use the tiny fixture to drive Inchworm-equivalent improvements until the output matches truth/transcript-oracle identity and not just length scale.
5. Only after the fast fixture is stable, decide the frozen public benchmark panel and require approval before changing it.
