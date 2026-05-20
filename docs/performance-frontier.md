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
  --out-root target/trinity_parity/public_panel_oriented_boundary_probe

python3 bench/trinity_parity/run_public_panel.py \
  --dataset trinity_source_test_assembly \
  --skip-raptor \
  --run-trinity \
  --timeout-seconds 120 \
  --out-root target/trinity_parity/public_panel_oriented_boundary_probe
```

Current Raptor output:

| Metric | Value |
| --- | ---: |
| Contigs | 30 |
| Total bases | 96,750 |
| N50 | 5,379 bp |
| Max length | 8,756 bp |
| Contigs >=1kb | 22 |
| Bases >=1kb | 93,247 |
| App-level workflow elapsed | 1.14 s |
| Harness wall time with release compile | 18.51 s |

Trinity in the same output root emits 74 transcripts, 149,325 bases, N50 5,399
bp, and max length 8,973 bp.

Strict parity still fails. Reciprocal FASTA F1 is 0.019231 at the frozen 0.95
coverage threshold, far below the required 0.9. The best observed reciprocal
coverage remains 0.999886 for one nearly full-length match, while the next
public boundary miss remains 0.917342. The current problem is still transcript
selection and boundary accuracy, not gross contiguity or raw speed.

Current diagnostic threshold sweep:

| Min coverage | Raptor matches | Trinity matches | F1 |
| --- | ---: | ---: | ---: |
| 0.900 | 2 / 30 | 2 / 74 | 0.038462 |
| 0.925 | 1 / 30 | 1 / 74 | 0.019231 |
| 0.950 | 1 / 30 | 1 / 74 | 0.019231 |
| 0.975 | 1 / 30 | 1 / 74 | 0.019231 |

Against Trinity's `inchworm.DS.fa`, Raptor's current contigs have stronger
stage-level overlap:

| Min coverage | Raptor matches | Trinity Inchworm matches | F1 |
| --- | ---: | ---: | ---: |
| 0.900 | 2 / 30 | 2 / 1,275 | 0.003066 |
| 0.925 | 1 / 30 | 1 / 1,275 | 0.001532 |
| 0.950 | 1 / 30 | 1 / 1,275 | 0.001532 |
| 0.975 | 1 / 30 | 1 / 1,275 | 0.001532 |

Interpretation: the oriented k-mer fallback now extends contig boundaries inside
the same oriented k-mer graph before ranking adaptive-k and Trinity-style k=25
candidates. This makes the public output much more Trinity-scale by contiguity
and slightly improves final FASTA F1, but it still leaves most public
transcripts unmatched. The next blocker is path selection/isoform reconstruction
across Chrysalis/Butterfly-style stages, not raw contig length alone.

## Accepted Changes

| Change | Public effect | Status |
| --- | --- | --- |
| Default `raptor trinity --min-len` raised from 50 to 200 | Reduced obvious short-fragment output while preserving synthetic gates | Kept |
| Component artifact cap above 1,000 contigs | Public fragmented runs finish instead of hanging in exhaustive component graph emission | Kept |
| Prefix index for exact read-overlap rescue | Public runtime improved from 209.67 s to 115.09 s with unchanged 4,392-contig output | Kept as a fallback speed fix |
| Paired-start fusion split bounded to compact contigs | Public output improved to 3,974 contigs, N50 338, max 1,052 without breaking compact-fusion fixture | Kept |
| Oriented k-mer path fallback before read-overlap rescue | Public output jumped to 31 contigs, N50 3,532, max 8,756; best reciprocal coverage improved from 0.634298 to 0.917342 | Current frontier |
| Adaptive oriented k plus Trinity-style k=25 candidate selection | Public final FASTA F1 moved from 0.0 to 0.018868 at 0.95 coverage; output is 30 contigs, N50 3,580, max 8,756 | Kept |
| Oriented graph boundary extension | Public output improved to 30 contigs, 96,750 bases, N50 5,379; final FASTA F1 moved to 0.019231 | Current frontier |
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
- Naive exact read-backed end extension is a trap. It lengthened 25 public
  contigs and raised total bases to 51,804, but final FASTA F1 regressed to 0.0,
  best final coverage dropped to 0.634298, and two one-base `N`s appeared in the
  output. Do not revive that path without a stronger graph/evidence gate.
- Unique-only oriented boundary extension is also a trap. It preserved the same
  strict final FASTA F1 (`0.019231`) and the same 0.90 sweep F1 (`0.038462`),
  but cut output from 96,750 bases and N50 5,379 down to 51,951 bases and N50
  3,580. The current best-count unused-neighbor boundary extension is the better
  frontier until path selection has a stronger evidence model.
- Full component graph artifact generation is not yet public-scale. Current
  caps are honest bounded-reporting guards, not a Chrysalis parity solution.

## Boundary Diagnostics

`run_public_panel.py` now records `best_matches_top20` and a
`containment_summary` for each FASTA comparison direction. Each match record
includes query/reference IDs, orientation, lengths, whether the query contains
the reference, containment offsets, and missing-end sizes. The summary counts
how many queries have any exact containment, how many are contained in the
reference, how many contain the reference, how many are near-full-length at
`>=0.95`, and the mean missing-end sizes.

Current oriented-boundary frontier diagnostics:

- `contig_1` matches `TRINITY_DN10_c0_g1_i3` at 8,756 / 8,757 bp
  (`0.999886` coverage), missing 1 bp on the left.
- `contig_14` matches `TRINITY_DN7_c0_g1_i1` at 3,507 / 3,823 bp
  (`0.917342` coverage), reverse-complemented in final Trinity, with 271 bp
  missing on one side and 45 bp on the other.
- `contig_9` is a new near-match to `TRINITY_DN3_c0_g1_i1` at 4,853 / 5,494 bp
  (`0.883327` coverage), missing 461 bp and 180 bp at the ends.
- Against Trinity Inchworm, the 3,507 bp near-match maps to `a301;19` with the
  same coverage and missing-end sizes, so this is still an Inchworm/path-boundary
  problem before final transcript selection.

Containment summary from
`target/trinity_parity/public_panel_oriented_boundary_probe/containment_summary_report.json`:

| Comparison | Direction | Any containment | Query contained in reference | Query contains reference | Near full length | Mean left missing | Mean right missing |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Trinity final | Raptor -> Trinity | 6 | 6 | 0 | 1 | 270.00 | 139.00 |
| Trinity final | Trinity -> Raptor | 7 | 0 | 7 | 1 | 475.57 | 216.43 |
| Trinity Inchworm | Raptor -> Inchworm | 16 | 5 | 11 | 1 | 578.81 | 1501.88 |
| Trinity Inchworm | Inchworm -> Raptor | 18 | 13 | 5 | 1 | 181.28 | 2163.78 |

Interpretation: the final Trinity comparison says most Raptor long contigs have
no exact containment in final Trinity transcripts, while the few hits are
mostly under-extended relative to Trinity. The Inchworm comparison shows both
under-contained and over-contained relationships, which points to graph path
selection and isoform boundary decisions rather than a simple "extend all ends"
fix.

## Next Serious Moves

1. Use the Raptor-vs-Trinity-Inchworm diagnostics to separate raw contig
   construction progress from final transcript selection failures.
2. Use the public report threshold sweep to track near-match progress at 0.90
   and 0.925 coverage while the hard 0.95 gate remains intact.
3. Tighten oriented path selection around the current 0.917 and 0.883 near-matches:
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
