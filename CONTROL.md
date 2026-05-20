# CONTROL

## Status Contract

status_file: PLAN.md
attempt_log: ATTEMPTS.md
durable_notes: NOTES.md
update_memory_after: every_experiment
check_control_before: phase_change, strategic_pivot, expensive_step, dependency_change, accelerator_backend_change, parity_claim

## Human Priorities

primary_priority: biological_correctness
secondary_priority: evidence_quality

## Scope Knobs

allowed_files:
- src/
- tests/
- benches/
- bench/
- scripts/
- docs/
- README.md
- Cargo.toml
- Cargo.lock

protected_files:
- none

max_blast_radius: stage_by_stage_changes_only

## Resource Knobs

max_runtime_per_step: none
max_parallel_jobs: host_reasonable
network_allowed: true
external_api_allowed: false

## Decision Gates

require_approval_for:
- strategic_pivot
- destructive_change
- dependency_change
- public_api_change
- benchmark_panel_change_after_freeze
- metric_threshold_weakening
- parity_completion_claim

## Sidecar Inputs

sidecar_apply_cadence: before_phase_change
nudge_file: none
human_overlay_file: none
review_queue_file: none

## Latest Human Nudge

Pursue Raptor to the end as a full Trinity replacement using the better infrastructure now in place. CUDA is allowed on this NVIDIA host when it is justified by a measured bottleneck and preserves biological output equivalence.
