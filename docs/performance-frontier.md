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
  --out-root target/trinity_parity/public_panel_oriented_k25_probe

python3 bench/trinity_parity/run_public_panel.py \
  --dataset trinity_source_test_assembly \
  --skip-raptor \
  --run-trinity \
  --timeout-seconds 120 \
  --out-root target/trinity_parity/public_panel_oriented_k25_probe
```

Current Raptor output:

| Metric | Value |
| --- | ---: |
| Contigs | 30 |
| Total bases | 50,612 |
| N50 | 3,580 bp |
| Max length | 8,756 bp |
| Contigs >=1kb | 12 |
| Bases >=1kb | 43,391 |
| App-level workflow elapsed | 1.35 s |
| Harness wall time with release compile | 18.40 s |

Trinity in the same output root emits 76 transcripts, 155,664 bases, N50 5,399
bp, and max length 8,973 bp.

Strict parity still fails. Reciprocal FASTA F1 is 0.018868 at the frozen 0.95
coverage threshold, far below the required 0.9. The best observed reciprocal
coverage remains 0.999886 for one nearly full-length match, while the next
public boundary miss remains 0.917342. The current problem is still transcript
selection and boundary accuracy, not gross contiguity or raw speed.

Current diagnostic threshold sweep:

| Min coverage | Raptor matches | Trinity matches | F1 |
| --- | ---: | ---: | ---: |
| 0.900 | 2 / 30 | 2 / 76 | 0.037736 |
| 0.925 | 1 / 30 | 1 / 76 | 0.018868 |
| 0.950 | 1 / 30 | 1 / 76 | 0.018868 |
| 0.975 | 1 / 30 | 1 / 76 | 0.018868 |

Against Trinity's `inchworm.DS.fa`, Raptor's current contigs have stronger
stage-level overlap:

| Min coverage | Raptor matches | Trinity Inchworm matches | F1 |
| --- | ---: | ---: | ---: |
| 0.900 | 3 / 30 | 3 / 1,278 | 0.004586 |
| 0.925 | 2 / 30 | 2 / 1,278 | 0.003058 |
| 0.950 | 1 / 30 | 1 / 1,278 | 0.001528 |
| 0.975 | 1 / 30 | 1 / 1,278 | 0.001528 |

Interpretation: the oriented k-mer fallback now tries the adaptive k and
Trinity-style k=25, then keeps the candidate with the best N50/longest/total
bases key. This slightly improves final Trinity FASTA matching, but it trades
off some Inchworm-stage reciprocal hits on the fresh Trinity run. The next
blocker is likely selection/boundary reconstruction across
Chrysalis/Butterfly-style stages, not only k selection or raw contig length.

## Accepted Changes

| Change | Public effect | Status |
| --- | --- | --- |
| Default `raptor trinity --min-len` raised from 50 to 200 | Reduced obvious short-fragment output while preserving synthetic gates | Kept |
| Component artifact cap above 1,000 contigs | Public fragmented runs finish instead of hanging in exhaustive component graph emission | Kept |
| Prefix index for exact read-overlap rescue | Public runtime improved from 209.67 s to 115.09 s with unchanged 4,392-contig output | Kept as a fallback speed fix |
| Paired-start fusion split bounded to compact contigs | Public output improved to 3,974 contigs, N50 338, max 1,052 without breaking compact-fusion fixture | Kept |
| Oriented k-mer path fallback before read-overlap rescue | Public output jumped to 31 contigs, N50 3,532, max 8,756; best reciprocal coverage improved from 0.634298 to 0.917342 | Current frontier |
| Adaptive oriented k plus Trinity-style k=25 candidate selection | Public final FASTA F1 moved from 0.0 to 0.018868 at 0.95 coverage; output is 30 contigs, N50 3,580, max 8,756 | Current frontier |
| Component artifact cap above 20,000 contig bases | Lets public oriented-contig output finish and be scored while full graph artifacts remain expensive | Kept as a bounded reporting guard |

## Limited Or Negative Results

- Prefix-index read-overlap rescue was a runtime win only. It did not improve
  public biological F1 or max contig length.
- Bounding paired-start fusion splitting helped contiguity, but only modestly.
  It did not move reciprocal FASTA F1 off 0.0.
- Oriented k-mer paths fixed the gross canonical-orientation bug and made output
  Trinity-scale, but still miss the strict 0.95 reciprocal coverage gate.
- Trying Trinity-style k=25 improves final FASTA reciprocal matching, but it
  does not solve the boundary problem and reduces fresh-run Inchworm-stage
  matches from the prior 31-contig frontier.
- Full component graph artifact generation is not yet public-scale. Current
  caps are honest bounded-reporting guards, not a Chrysalis parity solution.

## Next Serious Moves

1. Use the Raptor-vs-Trinity-Inchworm diagnostics to separate raw contig
   construction progress from final transcript selection failures.
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
