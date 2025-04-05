# 🦖 Raptor

**A blazing-fast, parallel, graph-based RNA-Seq assembler.**  
_K-mer powered. Isoform aware. Built for scale._

[![Build Status](https://img.shields.io/github/actions/workflow/status/Jakeelamb/Raptor/build.yml?branch=main)](https://github.com/Jakeelamb/Raptor/actions)  
[![Crates.io](https://img.shields.io/crates/v/raptor)](https://crates.io/crates/raptor)  
[![GPU Accelerated](https://img.shields.io/badge/GPU-accelerated-green)](#gpu-acceleration)  
[![License](https://img.shields.io/github/license/Jakeelamb/Raptor)](./LICENSE)

---

Raptor is a modern RNA-Seq assembler built for performance and biological accuracy. Inspired by Trinity, bbnorm, and SeqKit, Raptor supports:

- 🧠 **Greedy k-mer extension** (adaptive k, canonical hashing)
- ⚙️ **Parallel assembly** (Rayon, SIMD)
- ⚡ **GPU-accelerated k-mer normalization**
- 🔗 **Graph-based isoform stitching** (Butterfly-like traversal)
- 💾 **Streaming input and low-RAM support**
- 🧬 **Isoform filtering, polishing, quantification**
- 📈 **PCA, heatmaps, TPM matrices, and GTF export**

---

## 📦 Installation

> Requires Rust 1.72+ and optionally CUDA for GPU support.

```bash
git clone https://github.com/Jakeelamb/Raptor.git
cd raptor
cargo build --release

# Optional: compile with CUDA:
cargo build --release --features "gpu"
```

## 🚀 Quick Start

```bash
# Normalize reads using GPU-accelerated CMS
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

## 🧪 Example Outputs

| File | Description |
|------|-------------|
| my_assembly.fasta | Assembled contigs |
| my_assembly.gfa | GFA1 graph of overlaps |
| my_assembly_isoforms.fasta | Inferred transcripts |
| my_assembly_isoforms.gfa | Graph with isoform P lines |
| my_assembly_isoforms.gtf | GTF format annotation |
| my_assembly_isoform.counts.matrix | TPM + confidence scores |
| heatmap.png, pca.png | Visual TPM analysis |

## 🛠️ Key Features

✅ Adaptive k-mer selection  
✅ Paired-end support  
✅ Long-read polishing  
✅ Splicing-aware path inference  
✅ GFA2 + BandageNG annotations  
✅ JSON/TSV metadata export  
✅ Differential isoform comparison via GTF  

## 📚 Citations & References

If you use Raptor in your research, please cite the tool (citation coming soon) and the underlying software inspirations:

- Grabherr et al. "Full-length transcriptome assembly from RNA-Seq data without a reference genome." Nat Biotech (2011)

- BBTools (BBNorm): https://jgi.doe.gov/data-and-tools/bbtools

- SeqKit: https://bioinf.shenwei.me/seqkit/

## 🤝 Contributing

PRs welcome! Run `cargo fmt && cargo clippy` before submitting.
See CONTRIBUTING.md for details.

## 🧠 License

MIT © 2024 Jacob Lamb / Mueller Lab 