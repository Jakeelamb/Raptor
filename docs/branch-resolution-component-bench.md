# Branch-Resolution Component Bench

This component isolates the branch-choice heuristic in the large-genome
assembler.

It targets the decision logic behind:

- strong read-support wins
- weak support noise should not override closer coverage
- non-repeat candidates are preferred when support does not decide
- repeat fallback remains deterministic

## Commands

Generate one synthetic task:

```bash
cargo run -- component-bench prepare-branch-resolution \
  --output artifacts/autoresearch_raptor/branch_resolution/task_000 \
  --cases 64 \
  --seed 7
```

Generate a fixed panel:

```bash
cargo run -- component-bench prepare-branch-resolution-panel \
  --output artifacts/autoresearch_raptor/branch_resolution/panel_small \
  --tasks 16 \
  --cases-per-task 64 \
  --seed 7 \
  --seed-step 1
```

Evaluate a task or panel root:

```bash
cargo run -- component-bench branch-resolution \
  --task artifacts/autoresearch_raptor/branch_resolution/panel_small \
  --branch-support-min-win 2 \
  --branch-support-min-margin 1 \
  --json
```

## Task Layout

Each task directory contains:

```text
task_dir/
  cases.json
  task.json
```

`cases.json` stores a list of branch-choice cases. Each case contains:

- `scenario`
- `current_count`
- `expected_base_idx`
- `expected_next_kmer`
- `candidates[]`

Each candidate contains:

- `base_idx`
- `next_kmer`
- `count`
- `is_repeat`
- `read_support`

## Metrics

The harness reports:

- exact choice rate
- any-choice rate
- throughput in cases/second
- `score`

Mutable evaluator knobs:

- `branch_support_min_win`
- `branch_support_min_margin`
- `prefer_non_repeat_branches` via `--disable-prefer-non-repeat`

For this component, `score = exact_choice_rate`.

That is intentionally strict. This component is small enough that approximate
matches are not useful.

## Guardrails

Do not promote a local winner on this task alone.

A branch-resolution candidate is only acceptable if:

1. it improves branch-resolution validation score
2. it does not regress `quick_test`
3. it does not obviously worsen misassembly-oriented whole-pipeline metrics
