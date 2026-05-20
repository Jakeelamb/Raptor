# ATTEMPTS

| Time | Attempt | Evidence | Result | Next Adjustment |
| --- | --- | --- | --- | --- |
| 2026-05-20 | Rescued old GPU-enabled branch and established OpenCL baseline. | `cargo test --features gpu`; `./scripts/rescue_smoke.sh`; `./bench/gpu_kmer_baseline.sh` | Passed; OpenCL baseline was faster than CPU exact on small synthetic k-mer run. | Build Trinity parity map and benchmark harness. |

