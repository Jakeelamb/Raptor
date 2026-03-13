# Autoresearch Error-Correction Bridge

This bridge lives inside Raptor because the current external
`/home/jake/Projects/autoresearch` repo does not yet have a plugin interface for
arbitrary component scorers.

It mirrors the existing Raptor-side component bridges:

- a fixed train/val task root
- a held-out split for non-regression checks
- a tiny mutable search surface
- a scalar local score

## Files

- [scripts/autoresearch_bridge/error_correction_prepare.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/error_correction_prepare.py)
- [scripts/autoresearch_bridge/error_correction_train.py](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/error_correction_train.py)

## Prepare Panels

```bash
python3 scripts/autoresearch_bridge/error_correction_prepare.py \
  --binary target/debug/raptor \
  --output-root artifacts/autoresearch_raptor/error_correction/tasks
```

This writes:

```text
artifacts/autoresearch_raptor/error_correction/tasks/
  train/
  val/
  heldout/
  metadata.json
```

Default panel profiles deliberately mix:

- balanced correction/retention tasks
- correction-heavy tasks with more correctable singleton errors
- retention-heavy tasks with more protected rare variants
- held-out tasks with higher weak-root pressure so strict thresholds regress
  on trusted-target recovery instead of silently looking fine on train/val

## Run The Search Loop

```bash
python3 scripts/autoresearch_bridge/error_correction_train.py \
  --binary target/debug/raptor \
  --tasks-root artifacts/autoresearch_raptor/error_correction/tasks \
  --require-heldout
```

Developer smoke:

```bash
python3 scripts/autoresearch_bridge/error_correction_prepare.py \
  --binary target/debug/raptor \
  --output-root /tmp/raptor_autoresearch_error_correction
python3 scripts/autoresearch_bridge/error_correction_train.py \
  --binary target/debug/raptor \
  --tasks-root /tmp/raptor_autoresearch_error_correction \
  --smoke \
  --require-heldout
```

## Current Mutable Surface

Right now the bridge searches:

- explicit `min_count` values `1` through `5`
- explicit `min_trusted_count` floors `3` through `6`

That is intentionally narrow. This first pass is meant to prove the component
loop shape and expose whether the current synthetic panel yields a meaningful
local tradeoff before adding more knobs.

## Current Status

The first error-correction bridge is local-search only.

It does not yet ship with a promotion script because:

- the full assembler still uses adaptive `min_count=0` downstream
- the isolated error-correction stage does not yet expose enough knobs to make
  a default flip defensible on its own
- we should only add a release promotion gate after the local panel shows a
  real, non-trivial winner

The trainer still writes the latest local result to:

- [latest_result.json](/home/jake/Projects/Raptor/artifacts/autoresearch_raptor/error_correction/latest_result.json)

Latest checked local result on March 13, 2026:

- local winner: `min_count=1, min_trusted_count=4`
- train/val/heldout score: `1.0` / `1.0` / `1.0`
- eligible candidates after heldout non-regression: `1`
- lenient heldout comparison `min_count=1, min_trusted_count=3`: score
  `0.7800`, correction `F1 0.8101`, trusted exact `1.0`, preserved retention
  `0.5`
- strict heldout comparison `min_count=1, min_trusted_count=6`: score
  `0.8114`, correction `F1 0.7942`, trusted exact `0.6573`,
  preserved retention `1.0`

Interpretation:

- the bridge shape is working correctly
- the current heldout panel is now strong enough to reject both overly lenient
  and overly strict trusted-neighbor thresholds
- the widened surface still selects the live baseline cleanly, so promotion
  would still be premature without either a richer panel or an ambiguity-aware
  singleton-rescue guard

## Next Best Step

1. run the fixed error-correction train/val/heldout loop and see whether the
   search saturates at the most lenient threshold
2. if it still saturates, add an ambiguity guard or harder near-trusted heldout
   cases before spending time on a release promotion path
3. if a non-trivial winner appears, add a guarded promotion script and then
   compare the candidate against the live `assemble-large` path
