# Raptor

**A blazing-fast, parallel, graph-based RNA-Seq assembler.**  
_K-mer powered. Isoform aware. Built for scale._

[![Build Status](https://img.shields.io/github/actions/workflow/status/Jakeelamb/Raptor/build.yml?branch=main)](https://github.com/Jakeelamb/Raptor/actions)  
[![Crates.io](https://img.shields.io/crates/v/raptor)](https://crates.io/crates/raptor)  
[![GPU Accelerated](https://img.shields.io/badge/GPU-accelerated-green)](#gpu-acceleration)  
[![License](https://img.shields.io/github/license/Jakeelamb/Raptor)](./LICENSE)

---

Raptor is a modern RNA-Seq assembler built for performance and biological accuracy. Inspired by Trinity, bbnorm, and SeqKit, Raptor supports:

- **Greedy k-mer extension** (adaptive k, canonical hashing)
- **Parallel assembly** (Rayon, SIMD)
- **GPU-accelerated k-mer normalization and overlap detection through OpenCL**
- **Graph-based isoform stitching** (Butterfly-like traversal)
- **Streaming input and low-RAM support**
- **Isoform filtering, polishing, quantification**
- **PCA, heatmaps, TPM matrices, and GTF export**
- **HPC-oriented shell workflows**

---

## Installation

> Requires Rust 1.72+. GPU support currently uses OpenCL through the optional `gpu` feature.

```bash
git clone https://github.com/Jakeelamb/Raptor.git
cd raptor
cargo build --release

# Optional: compile with GPU support:
cargo build --release --features "gpu"

# HPC environments with module system:
./compile_hpc.sh         # Default build path for HPC environments
./compile_hpc.sh --gpu   # With GPU support
```

## Known-Good Rescue Checks

This branch was rescued on May 20, 2026 as `rescue/gpu-trinity`. The current restart checks are:

```bash
cargo check
cargo check --features gpu
cargo test
cargo test --features gpu
./scripts/rescue_smoke.sh
```

The GPU backend is OpenCL today. CUDA is a valid future backend for NVIDIA hardware, but it should be added only after the OpenCL path has a saved baseline to beat.

## Quick Start

```bash
# Normalize reads using GPU-accelerated counting when an OpenCL device is available
raptor normalize \
  -i sample_R1.fastq.gz \
  -o norm.fastq.gz \
  --gpu \
  --streaming

# Assemble transcriptome from normalized reads
raptor assemble \
  -i norm.fastq.gz \
  -o my_assembly \
  --threads 16 \
  --gfa --isoforms \
  --json-metadata metadata.json \
  --min-confidence 0.75

# Visualize transcript diversity
raptor stats --input my_assembly_isoform.counts.matrix --pca pca.png --heatmap heatmap.png
```

## Example Outputs

| File | Description |
|------|-------------|
| my_assembly.fasta | Assembled contigs |
| my_assembly.gfa | GFA1 graph of overlaps |
| my_assembly_isoforms.fasta | Inferred transcripts |
| my_assembly_isoforms.gfa | Graph with isoform P lines |
| my_assembly_isoforms.gtf | GTF format annotation |
| my_assembly_isoform.counts.matrix | TPM + confidence scores |
| heatmap.png, pca.png | Visual TPM analysis |

## Key Features

- Adaptive k-mer selection
- Paired-end support
- Long-read polishing
- Splicing-aware path inference
- GFA2 + BandageNG annotations
- JSON/TSV metadata export
- Differential isoform comparison via GTF
- Optional OpenCL acceleration for k-mer counting and overlap detection
- HPC-ready shell scripts for build, submit, and monitoring workflows

## Benchmarking

Raptor now ships with a reproducible genome-assembly benchmark workflow under [bench/genome_assembly](./bench/genome_assembly).

Use it to:

- download or generate benchmark datasets
- run Raptor and comparator assemblers with fixed commands
- capture runtime, peak memory, and assembly statistics
- generate machine-readable summaries and markdown reports

Quick entry points:

```bash
# Install benchmark dependencies
./bench/genome_assembly/setup_environment.sh

# Generate simulated data or download public datasets
./bench/genome_assembly/download_data.sh simulated

# Run a benchmark
./bench/genome_assembly/run_benchmark.sh simulated 8

# Aggregate all benchmark runs into CSV/Markdown summaries
python3 ./bench/genome_assembly/summarize_results.py
```

Methodology and publishing guidance live in [bench/genome_assembly/README.md](./bench/genome_assembly/README.md).

The current checked-in baseline report is in [docs/genome-benchmark-report.md](./docs/genome-benchmark-report.md).

The generated cross-tool comparison view is in [docs/genome-benchmark-comparison.md](./docs/genome-benchmark-comparison.md).

The next-generation assembler redesign plan is in [docs/raptor-architecture-roadmap.md](./docs/raptor-architecture-roadmap.md).

## HPC Support

Raptor includes scripts specifically designed for high-performance computing environments:

- `compile_hpc.sh` - Easy compilation with module detection
- `raptor_hpc.sh` - Job submission script for SLURM/PBS/SGE
- `monitor_hpc.sh` - Job monitoring tool for resource usage
- `test_features.sh` - Validate builds with different features

For detailed HPC setup instructions, see [HPC_INSTRUCTIONS.md](./HPC_INSTRUCTIONS.md).

## Citations & References

If you use Raptor in your research, please cite the tool (citation coming soon) and the underlying software inspirations:

- Grabherr et al. "Full-length transcriptome assembly from RNA-Seq data without a reference genome." Nat Biotech (2011)

- BBTools (BBNorm): https://jgi.doe.gov/data-and-tools/bbtools

- SeqKit: https://bioinf.shenwei.me/seqkit/

## Contributing

PRs welcome! Run `cargo fmt && cargo clippy` before submitting.
See CONTRIBUTING.md for details.

## License

MIT © 2024 Jacob Lamb / Mueller Lab 
