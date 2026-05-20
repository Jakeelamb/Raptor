# Raptor Rescue Baseline

Date: 2026-05-20
Branch: `rescue/gpu-trinity`
Baseline commits:

- `299c294 Revive GPU count smoke CLI`
- `954a255 Polish rescued GPU baseline`

## What Was Rescued

The restored line is the pre-rewrite Rust assembler with:

- OpenCL GPU k-mer counting and overlap detection behind `--features gpu`
- Rayon/SIMD CPU paths
- isoform graph and transcript workflows
- disk-backed large-genome assembly pipeline
- benchmark and HPC shell workflows

The abandoned tiny rewrite was preserved before switching back:

- `stash@{0}: preserve foundation rewrite before GPU rescue`
- `/tmp/raptor-foundation-rewrite-2026-05-20.patch`

## Current Truth

GPU support currently means OpenCL through the `ocl` crate. The README previously said CUDA; that was stale. CUDA can be added later as a separate NVIDIA-focused backend, but the OpenCL path is the current working implementation and should be the baseline to beat.

The repository does not currently define an `mpi-support` Cargo feature. HPC scripts exist, but MPI feature claims should not be treated as active Cargo functionality until the feature is reintroduced and tested.

## Baseline Commands

```bash
cargo check
cargo check --features gpu
cargo test
cargo test --features gpu
./scripts/rescue_smoke.sh
```

## Last Verified Results

- `cargo check`: passed
- `cargo check --features gpu`: passed
- `cargo test`: passed
- `cargo test --features gpu`: passed
- `./scripts/rescue_smoke.sh`: passed with `gpu_smoke=ok`
- Current host GPU runtime status: `nvidia-smi` reports `NVIDIA GeForce RTX 5070 Laptop GPU`, driver `595.71.05`, CUDA runtime `13.2`, and 8151 MiB VRAM.
- OpenCL k-mer smoke completed and reported GPU counting.

## CPU vs OpenCL Baseline

Run:

```bash
./bench/gpu_kmer_baseline.sh
```

Default input shape:

- reads: 20000
- read length: 150
- synthetic genome length: 50000
- k: 25

The script writes:

- `target/rescue-baseline/kmer_cpu_opencl.tsv`
- `target/rescue-baseline/cpu.log`
- `target/rescue-baseline/opencl.log`

First measured run on the RTX 5070 Laptop GPU:

| Backend | Status | End-to-end seconds | Unique k-mers |
|---|---:|---:|---:|
| CPU exact | ok | 3.424227 | 49974 |
| OpenCL | ok | 0.648744 | 49974 |

The OpenCL binary also reported kernel-side counting time of `231.493299ms`.

## CUDA Note

Jake has NVIDIA hardware, so CUDA is a reasonable future backend. Do not add it as a blind rewrite. The right path is:

1. Keep `nvidia-smi` and OpenCL smoke green.
2. Save CPU and OpenCL baselines with `./bench/gpu_kmer_baseline.sh`.
3. Add a CUDA backend behind a separate feature, such as `cuda`, with the same public counting contract.
4. Compare OpenCL, CUDA, and CPU on the same generated FASTQ inputs.

## Next Polish Targets

1. Keep GPU/OpenCL docs and command examples honest.
2. Add measured CPU-vs-GPU k-mer counting baselines on generated FASTQ data.
3. Fix remaining production-path `unwrap()` and `expect()` calls in binaries.
4. Audit direct dependencies and feature flags before upgrading versions.
5. Evaluate CUDA only after the OpenCL baseline is saved with comparable data.
