# Performance Frontier

Status date: 2026-05-20
Branch: `rescue/gpu-trinity`

This note records the current Trinity-replacement performance frontier so future
optimization work does not repeat local wins that failed the public biological
gate. The public benchmark lane is the Trinity source `test_Trinity_Assembly`
dataset materialized under `data/trinity_public_panel/`.

## Current Public Frontier

Command shape:

```bash
python3 bench/trinity_parity/run_public_panel.py \
  --dataset trinity_source_test_assembly \
  --timeout-seconds 180 \
  --out-root target/trinity_parity/public_panel_oriented_kmers_capped

python3 bench/trinity_parity/run_public_panel.py \
  --dataset trinity_source_test_assembly \
  --skip-raptor \
  --run-trinity \
  --timeout-seconds 120 \
  --out-root target/trinity_parity/public_panel_oriented_kmers_capped
```

Current Raptor output:

| Metric | Value |
| --- | ---: |
| Contigs | 31 |
| Total bases | 50,507 |
| N50 | 3,532 bp |
| Max length | 8,756 bp |
| Contigs >=1kb | 12 |
| Bases >=1kb | 43,274 |
| App-level workflow elapsed | 1.20 s |
| Harness wall time with release compile | 14.99 s |

Trinity in the same output root emits 74 transcripts, 132,171 bases, N50 3,697
bp, and max length 8,973 bp.

Strict parity still fails. Reciprocal FASTA F1 is 0.0 at the frozen 0.95
coverage threshold because no reciprocal match reaches that threshold. The best
observed reciprocal coverage is 0.917342, so the current problem is path
correctness/isoform boundary accuracy, not gross contiguity or raw speed.

Current diagnostic threshold sweep:

| Min coverage | Raptor matches | Trinity matches | F1 |
| --- | ---: | ---: | ---: |
| 0.900 | 1 / 31 | 1 / 74 | 0.019048 |
| 0.925 | 0 / 31 | 0 / 74 | 0.0 |
| 0.950 | 0 / 31 | 0 / 74 | 0.0 |
| 0.975 | 0 / 31 | 0 / 74 | 0.0 |

## Accepted Changes

| Change | Public effect | Status |
| --- | --- | --- |
| Default `raptor trinity --min-len` raised from 50 to 200 | Reduced obvious short-fragment output while preserving synthetic gates | Kept |
| Component artifact cap above 1,000 contigs | Public fragmented runs finish instead of hanging in exhaustive component graph emission | Kept |
| Prefix index for exact read-overlap rescue | Public runtime improved from 209.67 s to 115.09 s with unchanged 4,392-contig output | Kept as a fallback speed fix |
| Paired-start fusion split bounded to compact contigs | Public output improved to 3,974 contigs, N50 338, max 1,052 without breaking compact-fusion fixture | Kept |
| Oriented k-mer path fallback before read-overlap rescue | Public output jumped to 31 contigs, N50 3,532, max 8,756; best reciprocal coverage improved from 0.634298 to 0.917342 | Current frontier |
| Component artifact cap above 20,000 contig bases | Lets public oriented-contig output finish and be scored while full graph artifacts remain expensive | Kept as a bounded reporting guard |

## Limited Or Negative Results

- Prefix-index read-overlap rescue was a runtime win only. It did not improve
  public biological F1 or max contig length.
- Bounding paired-start fusion splitting helped contiguity, but only modestly.
  It did not move reciprocal FASTA F1 off 0.0.
- Oriented k-mer paths fixed the gross canonical-orientation bug and made output
  Trinity-scale, but still miss the strict 0.95 reciprocal coverage gate.
- Full component graph artifact generation is not yet public-scale. Current
  caps are honest bounded-reporting guards, not a Chrysalis parity solution.

## Next Serious Moves

1. Compare Raptor oriented k-mer contigs directly against Trinity
   `inchworm.DS.fa`, not only final `Trinity.fasta`.
2. Use the public report threshold sweep to track near-match progress at 0.90
   and 0.925 coverage while the hard 0.95 gate remains intact.
3. Tighten oriented path selection around the current 0.917 best match:
   investigate whether the miss is truncation, over-extension, reverse-complement
   orientation, or isoform boundary selection.
4. Replace capped public component artifacts with incremental/read-assignment
   summaries before claiming Chrysalis-scale performance.
5. Keep CUDA/OpenCL work isolated behind CPU equivalence gates. GPU speed is not
   useful until biological FASTA output matches the CPU path and public Trinity
   comparison remains stable.

## Stop Rule

Do not make a performance claim from a local microbenchmark alone. A change only
counts if it preserves the synthetic candidate panel and improves or explains
the public Trinity source case under `run_public_panel.py`.
