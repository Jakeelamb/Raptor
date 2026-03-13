# Autoresearch Campaign Runner

Date: 2026-03-13

Raptor now has a resumable campaign runner around the existing component
bridges.

The goal is not to replace the current per-component prepare/train/promote
scripts.
It is to make them run unattended for long local campaigns with:

- one manifest
- one SQLite state store
- isolated per-campaign artifacts
- resumable stage execution
- optional multi-component parallelism

## Entry Point

Use:

```bash
python3 scripts/autoresearch_bridge/run_campaign.py run \
  --manifest scripts/autoresearch_bridge/manifests/default_campaign.json
```

For the higher-budget unattended path, use:

```bash
python3 scripts/autoresearch_bridge/run_campaign.py run \
  --manifest scripts/autoresearch_bridge/manifests/overnight_campaign.json
```

Useful commands:

```bash
python3 scripts/autoresearch_bridge/run_campaign.py status
python3 scripts/autoresearch_bridge/run_campaign.py status --campaign-id 3
python3 scripts/autoresearch_bridge/run_campaign.py resume --campaign-id 3
```

The default and overnight manifests now also run contig-extraction promotion,
so campaign roots may include both `quick_test/` and repeat-heavy stress-suite
artifacts under `contig_extraction/promote/`.

## Manifest Shape

The manifest is JSON.

Top-level keys:

- `name`
- `max_workers`
- `fail_fast`
- `default_debug_binary` (optional)
- `default_release_binary` (optional)
- `campaigns_root` (optional)
- `database_path` (optional)
- `components`

Each component entry supports:

- `component`: one of `read_mapping`, `branch_resolution`,
  `scaffold_polish`, `error_correction`, `contig_extraction`
- `tasks_root`: optional external task root to reuse instead of campaign-local
  tasks
- `prepare`
- `train`
- `promote`

Each stage object supports:

- `enabled`
- `binary`
- `args`

The runner owns these flags and rejects them if they appear in `args`:

- prepare: `--binary`, `--output-root`
- train: `--binary`, `--tasks-root`, `--output`
- promote: `--binary`, `--tasks-root`, `--output-root`, `--train-result`

Reference manifests:

- [default_campaign.json](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/manifests/default_campaign.json)
- [overnight_campaign.json](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/manifests/overnight_campaign.json)
- [smoke_campaign.json](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/manifests/smoke_campaign.json)
- [read_mapping_halving_smoke.json](/home/jake/Projects/Raptor/scripts/autoresearch_bridge/manifests/read_mapping_halving_smoke.json)

## Artifact Layout

Each run gets an isolated root under:

```text
artifacts/autoresearch_raptor/campaigns/<name>_<timestamp>/
  manifest.json
  campaign.json
  <component>/
    tasks/
    prepare/
      attempt_001.stdout.log
      attempt_001.stderr.log
    train/
      latest_result.json
      attempt_001.stdout.log
      attempt_001.stderr.log
    promote/
      attempt_001.stdout.log
      attempt_001.stderr.log
      <component>_<timestamp>/
        summary.json
        quick_test/
        repeat_heavy/
```

That means a campaign never overwrites the global component
`latest_result.json` files.

## SQLite State

The runner stores campaign state in:

- default: `artifacts/autoresearch_raptor/campaigns/campaigns.sqlite3`

Tracked tables:

- `campaigns`
- `component_runs`
- `stage_runs`
- `artifacts`

This is enough to answer:

- what campaign is running
- which stage failed
- which result file was produced
- whether a resume needs to rerun a stale stage

If a campaign is interrupted, `resume` marks any stale `running` stages as
failed and schedules them again.

## Current Scope

This is intentionally mostly stage-level orchestration.

Current exception:

- `read_mapping` now supports candidate-level successive halving inside its
  train stage, and the default campaign manifest enables that policy
- the validated smoke run narrowed an `81`-config mapper grid to `27` `val`
  candidates and `12` `heldout` candidates, reducing split evaluations from
  `243` to `120`

Today it automates:

1. component task preparation
2. component-local search
3. optional promotion through heldout plus `quick_test`
4. component-specific end-to-end stress suites when a promotion script defines
   them

It does not yet do:

- candidate-level successive halving across the raw parameter grid
- multi-agent git lineage
- automatic default flips
- cross-campaign leaderboard views

That is the right next layer once the stage-level system has proven stable.
