#!/usr/bin/env python3
"""Run resumable autonomous autoresearch campaigns across Raptor components."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import time
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from bridge_common import extract_last_json_document, format_command
from campaign_db import CampaignDb, decode_json


REPO_ROOT = Path(__file__).resolve().parents[2]
BRIDGE_ROOT = REPO_ROOT / "scripts" / "autoresearch_bridge"
DEFAULT_CAMPAIGNS_ROOT = REPO_ROOT / "artifacts" / "autoresearch_raptor" / "campaigns"
DEFAULT_DB_PATH = DEFAULT_CAMPAIGNS_ROOT / "campaigns.sqlite3"
DEFAULT_DEBUG_BINARY = REPO_ROOT / "target" / "debug" / "raptor"
DEFAULT_RELEASE_BINARY = REPO_ROOT / "target" / "release" / "raptor"


@dataclass(frozen=True)
class ComponentSpec:
    name: str
    prepare_script: str
    train_script: str
    promote_script: str | None


@dataclass(frozen=True)
class StageConfig:
    enabled: bool
    args: tuple[str, ...]
    binary: Path | None


@dataclass(frozen=True)
class ComponentConfig:
    component: str
    tasks_root: Path | None
    prepare: StageConfig
    train: StageConfig
    promote: StageConfig


@dataclass(frozen=True)
class CampaignConfig:
    name: str
    manifest_path: Path
    campaigns_root: Path
    db_path: Path
    max_workers: int
    fail_fast: bool
    default_debug_binary: Path
    default_release_binary: Path
    components: tuple[ComponentConfig, ...]


@dataclass(frozen=True)
class ResolvedComponent:
    config: ComponentConfig
    spec: ComponentSpec
    artifact_root: Path
    tasks_root: Path
    prepare_log_root: Path
    train_log_root: Path
    train_result_path: Path
    promote_log_root: Path


@dataclass(frozen=True)
class StageExecution:
    stage_run_id: int
    component_run_id: int
    component_name: str
    stage_name: str
    attempt: int
    command: list[str]
    workdir: Path
    stdout_path: Path
    stderr_path: Path
    tasks_root: Path
    train_result_path: Path
    promote_log_root: Path


@dataclass(frozen=True)
class StageArtifact:
    artifact_type: str
    path: Path
    metadata: dict[str, Any] | None = None


@dataclass(frozen=True)
class StageOutcome:
    stage_run_id: int
    component_run_id: int
    component_name: str
    stage_name: str
    success: bool
    exit_code: int
    duration_seconds: float
    result_path: Path | None
    result_json: dict[str, Any] | None
    error_message: str | None
    artifacts: tuple[StageArtifact, ...]


COMPONENT_SPECS: dict[str, ComponentSpec] = {
    "read_mapping": ComponentSpec(
        name="read_mapping",
        prepare_script="read_mapping_prepare.py",
        train_script="read_mapping_train.py",
        promote_script="promote_read_mapping_candidate.py",
    ),
    "branch_resolution": ComponentSpec(
        name="branch_resolution",
        prepare_script="branch_resolution_prepare.py",
        train_script="branch_resolution_train.py",
        promote_script="promote_branch_resolution_candidate.py",
    ),
    "scaffold_polish": ComponentSpec(
        name="scaffold_polish",
        prepare_script="scaffold_polish_prepare.py",
        train_script="scaffold_polish_train.py",
        promote_script="promote_scaffold_polish_candidate.py",
    ),
    "error_correction": ComponentSpec(
        name="error_correction",
        prepare_script="error_correction_prepare.py",
        train_script="error_correction_train.py",
        promote_script=None,
    ),
    "contig_extraction": ComponentSpec(
        name="contig_extraction",
        prepare_script="contig_extraction_prepare.py",
        train_script="contig_extraction_train.py",
        promote_script="promote_contig_extraction_candidate.py",
    ),
}

RESERVED_STAGE_FLAGS = {
    "prepare": {"--binary", "--output-root"},
    "train": {"--binary", "--tasks-root", "--output"},
    "promote": {"--binary", "--tasks-root", "--output-root", "--train-result"},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run or resume resumable autoresearch campaigns with SQLite state tracking"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Start a new campaign from a manifest")
    run_parser.add_argument("--manifest", type=Path, required=True)
    run_parser.add_argument("--db", type=Path)
    run_parser.add_argument("--campaigns-root", type=Path)
    run_parser.add_argument("--max-workers", type=int)
    run_parser.add_argument("--fail-fast", action="store_true")
    run_parser.add_argument("--no-fail-fast", action="store_true")

    resume_parser = subparsers.add_parser("resume", help="Resume an existing campaign by id")
    resume_parser.add_argument("--campaign-id", type=int, required=True)
    resume_parser.add_argument("--db", type=Path, default=DEFAULT_DB_PATH)
    resume_parser.add_argument("--max-workers", type=int)
    resume_parser.add_argument("--fail-fast", action="store_true")
    resume_parser.add_argument("--no-fail-fast", action="store_true")

    status_parser = subparsers.add_parser("status", help="Show campaign status")
    status_parser.add_argument("--campaign-id", type=int)
    status_parser.add_argument("--db", type=Path, default=DEFAULT_DB_PATH)
    status_parser.add_argument("--limit", type=int, default=10)

    return parser.parse_args()


def slugify(text: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", text.strip()).strip("_").lower()
    return slug or "campaign"


def resolve_path(raw: str | Path, *, base_dir: Path) -> Path:
    path = Path(raw)
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    return path


def normalize_args(value: Any, *, stage_name: str) -> tuple[str, ...]:
    if value is None:
        return ()
    if not isinstance(value, list):
        raise ValueError(f"{stage_name}.args must be a JSON array")
    args = tuple(str(item) for item in value)
    reserved = RESERVED_STAGE_FLAGS[stage_name]
    for item in args:
        if item in reserved or any(item.startswith(f"{flag}=") for flag in reserved):
            raise ValueError(
                f"{stage_name}.args may not include reserved flag {item}; the campaign runner owns it"
            )
    return args


def parse_stage_config(
    payload: dict[str, Any] | None,
    *,
    default_enabled: bool,
    stage_name: str,
    base_dir: Path,
) -> StageConfig:
    if payload is None:
        return StageConfig(enabled=default_enabled, args=(), binary=None)
    if not isinstance(payload, dict):
        raise ValueError(f"{stage_name} must be an object")
    binary = payload.get("binary")
    return StageConfig(
        enabled=bool(payload.get("enabled", default_enabled)),
        args=normalize_args(payload.get("args"), stage_name=stage_name),
        binary=resolve_path(binary, base_dir=base_dir) if binary is not None else None,
    )


def load_manifest(
    manifest_path: Path,
    *,
    db_path_override: Path | None,
    campaigns_root_override: Path | None,
    max_workers_override: int | None,
    fail_fast_override: bool | None,
) -> CampaignConfig:
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    base_dir = manifest_path.parent.resolve()
    components_payload = payload.get("components")
    if not isinstance(components_payload, list) or not components_payload:
        raise ValueError("manifest must define a non-empty components array")

    default_debug_binary = resolve_path(
        payload.get("default_debug_binary", DEFAULT_DEBUG_BINARY),
        base_dir=base_dir,
    )
    default_release_binary = resolve_path(
        payload.get("default_release_binary", DEFAULT_RELEASE_BINARY),
        base_dir=base_dir,
    )

    components: list[ComponentConfig] = []
    seen: set[str] = set()
    for entry in components_payload:
        if not isinstance(entry, dict):
            raise ValueError("each components entry must be an object")
        component_name = entry.get("component")
        if not isinstance(component_name, str) or component_name not in COMPONENT_SPECS:
            raise ValueError(
                f"unknown component {component_name!r}; expected one of {sorted(COMPONENT_SPECS)}"
            )
        if component_name in seen:
            raise ValueError(f"duplicate component entry for {component_name}")
        seen.add(component_name)

        if not bool(entry.get("enabled", True)):
            continue

        tasks_root = entry.get("tasks_root")
        prepare = parse_stage_config(
            entry.get("prepare"),
            default_enabled=True,
            stage_name="prepare",
            base_dir=base_dir,
        )
        train = parse_stage_config(
            entry.get("train"),
            default_enabled=True,
            stage_name="train",
            base_dir=base_dir,
        )
        promote_default_enabled = COMPONENT_SPECS[component_name].promote_script is not None
        promote = parse_stage_config(
            entry.get("promote"),
            default_enabled=promote_default_enabled,
            stage_name="promote",
            base_dir=base_dir,
        )
        if promote.enabled and COMPONENT_SPECS[component_name].promote_script is None:
            raise ValueError(f"{component_name} does not have a promotion script yet")
        if promote.enabled and not train.enabled:
            raise ValueError(f"{component_name} cannot enable promote while train is disabled")

        components.append(
            ComponentConfig(
                component=component_name,
                tasks_root=resolve_path(tasks_root, base_dir=base_dir) if tasks_root else None,
                prepare=prepare,
                train=train,
                promote=promote,
            )
        )

    if not components:
        raise ValueError("manifest does not enable any components")

    fail_fast = bool(payload.get("fail_fast", True))
    if fail_fast_override is not None:
        fail_fast = fail_fast_override

    max_workers = int(payload.get("max_workers", 1))
    if max_workers_override is not None:
        max_workers = max_workers_override
    if max_workers < 1:
        raise ValueError("max_workers must be at least 1")

    campaigns_root = (
        campaigns_root_override.resolve()
        if campaigns_root_override is not None
        else resolve_path(payload.get("campaigns_root", DEFAULT_CAMPAIGNS_ROOT), base_dir=base_dir)
    )
    db_path = (
        db_path_override.resolve()
        if db_path_override is not None
        else resolve_path(payload.get("database_path", DEFAULT_DB_PATH), base_dir=base_dir)
    )

    return CampaignConfig(
        name=str(payload.get("name", manifest_path.stem)),
        manifest_path=manifest_path.resolve(),
        campaigns_root=campaigns_root,
        db_path=db_path,
        max_workers=max_workers,
        fail_fast=fail_fast,
        default_debug_binary=default_debug_binary,
        default_release_binary=default_release_binary,
        components=tuple(components),
    )


def manifest_to_jsonable(config: CampaignConfig) -> dict[str, Any]:
    return {
        "name": config.name,
        "manifest_path": str(config.manifest_path),
        "campaigns_root": str(config.campaigns_root),
        "database_path": str(config.db_path),
        "max_workers": config.max_workers,
        "fail_fast": config.fail_fast,
        "default_debug_binary": str(config.default_debug_binary),
        "default_release_binary": str(config.default_release_binary),
        "components": [
            {
                "component": component.component,
                "tasks_root": str(component.tasks_root) if component.tasks_root else None,
                "prepare": {
                    "enabled": component.prepare.enabled,
                    "binary": str(component.prepare.binary) if component.prepare.binary else None,
                    "args": list(component.prepare.args),
                },
                "train": {
                    "enabled": component.train.enabled,
                    "binary": str(component.train.binary) if component.train.binary else None,
                    "args": list(component.train.args),
                },
                "promote": {
                    "enabled": component.promote.enabled,
                    "binary": str(component.promote.binary) if component.promote.binary else None,
                    "args": list(component.promote.args),
                },
            }
            for component in config.components
        ],
    }


def resolve_component_layout(config: CampaignConfig, campaign_root: Path) -> list[ResolvedComponent]:
    resolved: list[ResolvedComponent] = []
    for component in config.components:
        spec = COMPONENT_SPECS[component.component]
        artifact_root = campaign_root / component.component
        tasks_root = component.tasks_root or artifact_root / "tasks"
        resolved.append(
            ResolvedComponent(
                config=component,
                spec=spec,
                artifact_root=artifact_root,
                tasks_root=tasks_root,
                prepare_log_root=artifact_root / "prepare",
                train_log_root=artifact_root / "train",
                train_result_path=artifact_root / "train" / "latest_result.json",
                promote_log_root=artifact_root / "promote",
            )
        )
    return resolved


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def find_latest_summary_json(output_root: Path) -> Path | None:
    candidates = [path for path in output_root.glob("*/summary.json") if path.is_file()]
    if not candidates:
        return None
    return max(candidates, key=lambda path: path.stat().st_mtime)


def parse_json_fallback(stdout_path: Path, stderr_path: Path) -> dict[str, Any] | None:
    streams: list[str] = []
    for path in (stdout_path, stderr_path):
        if path.exists():
            streams.append(path.read_text(encoding="utf-8", errors="replace"))
    combined = "\n".join(streams)
    if not combined.strip():
        return None
    try:
        value = extract_last_json_document(combined)
    except json.JSONDecodeError:
        return None
    return value if isinstance(value, dict) else None


def execute_stage(execution: StageExecution) -> StageOutcome:
    execution.stdout_path.parent.mkdir(parents=True, exist_ok=True)
    execution.stderr_path.parent.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    completed: subprocess.CompletedProcess[str] | None = None
    artifacts = [
        StageArtifact("stdout_log", execution.stdout_path),
        StageArtifact("stderr_log", execution.stderr_path),
    ]
    try:
        with execution.stdout_path.open("w", encoding="utf-8") as stdout_handle, execution.stderr_path.open(
            "w", encoding="utf-8"
        ) as stderr_handle:
            completed = subprocess.run(
                execution.command,
                cwd=execution.workdir,
                text=True,
                stdout=stdout_handle,
                stderr=stderr_handle,
                check=False,
            )
        duration_seconds = time.monotonic() - start

        result_path: Path | None = None
        result_json: dict[str, Any] | None = None
        error_message: str | None = None

        if completed.returncode == 0:
            if execution.stage_name == "prepare":
                metadata_path = execution.tasks_root / "metadata.json"
                if metadata_path.exists():
                    result_path = metadata_path
                    result_json = json.loads(metadata_path.read_text(encoding="utf-8"))
                    artifacts.append(
                        StageArtifact(
                            "prepare_metadata",
                            metadata_path,
                            {"tasks_root": str(execution.tasks_root)},
                        )
                    )
                    artifacts.append(StageArtifact("tasks_root", execution.tasks_root))
            elif execution.stage_name == "train":
                if execution.train_result_path.exists():
                    result_path = execution.train_result_path
                    result_json = json.loads(execution.train_result_path.read_text(encoding="utf-8"))
                    artifacts.append(StageArtifact("train_result", execution.train_result_path))
            elif execution.stage_name == "promote":
                summary_path = find_latest_summary_json(execution.promote_log_root)
                if summary_path is not None:
                    result_path = summary_path
                    result_json = json.loads(summary_path.read_text(encoding="utf-8"))
                    artifacts.append(StageArtifact("promotion_summary", summary_path))

            if result_json is None:
                result_json = parse_json_fallback(execution.stdout_path, execution.stderr_path)
            success = True
        else:
            stderr_text = execution.stderr_path.read_text(encoding="utf-8", errors="replace")
            stdout_text = execution.stdout_path.read_text(encoding="utf-8", errors="replace")
            trimmed_stderr = stderr_text.strip().splitlines()[-20:]
            trimmed_stdout = stdout_text.strip().splitlines()[-10:]
            error_lines = []
            if trimmed_stderr:
                error_lines.append("\n".join(trimmed_stderr))
            if trimmed_stdout:
                error_lines.append("\n".join(trimmed_stdout))
            error_message = "\n".join(part for part in error_lines if part).strip() or (
                f"{execution.component_name}:{execution.stage_name} failed"
            )
            success = False

        return StageOutcome(
            stage_run_id=execution.stage_run_id,
            component_run_id=execution.component_run_id,
            component_name=execution.component_name,
            stage_name=execution.stage_name,
            success=success,
            exit_code=completed.returncode,
            duration_seconds=duration_seconds,
            result_path=result_path,
            result_json=result_json,
            error_message=error_message,
            artifacts=tuple(artifacts),
        )
    except Exception as exc:
        duration_seconds = time.monotonic() - start
        return StageOutcome(
            stage_run_id=execution.stage_run_id,
            component_run_id=execution.component_run_id,
            component_name=execution.component_name,
            stage_name=execution.stage_name,
            success=False,
            exit_code=-1 if completed is None else completed.returncode,
            duration_seconds=duration_seconds,
            result_path=None,
            result_json=None,
            error_message=str(exc),
            artifacts=tuple(artifacts),
        )


def build_prepare_command(component: ResolvedComponent, campaign: CampaignConfig) -> list[str]:
    binary = component.config.prepare.binary or campaign.default_debug_binary
    return [
        "python3",
        str(BRIDGE_ROOT / component.spec.prepare_script),
        "--binary",
        str(binary),
        "--output-root",
        str(component.tasks_root),
        *component.config.prepare.args,
    ]


def build_train_command(component: ResolvedComponent, campaign: CampaignConfig) -> list[str]:
    binary = component.config.train.binary or campaign.default_debug_binary
    return [
        "python3",
        str(BRIDGE_ROOT / component.spec.train_script),
        "--binary",
        str(binary),
        "--tasks-root",
        str(component.tasks_root),
        "--output",
        str(component.train_result_path),
        *component.config.train.args,
    ]


def build_promote_command(component: ResolvedComponent, campaign: CampaignConfig) -> list[str]:
    if component.spec.promote_script is None:
        raise ValueError(f"{component.config.component} has no promotion script")
    binary = component.config.promote.binary or campaign.default_release_binary
    return [
        "python3",
        str(BRIDGE_ROOT / component.spec.promote_script),
        "--binary",
        str(binary),
        "--tasks-root",
        str(component.tasks_root),
        "--output-root",
        str(component.promote_log_root),
        "--train-result",
        str(component.train_result_path),
        *component.config.promote.args,
    ]


def create_stage_execution(
    *,
    db: CampaignDb,
    component_run_id: int,
    component: ResolvedComponent,
    stage_name: str,
    campaign: CampaignConfig,
) -> StageExecution:
    attempt = db.next_stage_attempt(component_run_id, stage_name)
    log_root = {
        "prepare": component.prepare_log_root,
        "train": component.train_log_root,
        "promote": component.promote_log_root,
    }[stage_name]
    stdout_path = log_root / f"attempt_{attempt:03d}.stdout.log"
    stderr_path = log_root / f"attempt_{attempt:03d}.stderr.log"
    if stage_name == "prepare":
        command = build_prepare_command(component, campaign)
    elif stage_name == "train":
        command = build_train_command(component, campaign)
    else:
        command = build_promote_command(component, campaign)

    stage_run_id = db.start_stage_run(
        component_run_id=component_run_id,
        stage_name=stage_name,
        attempt=attempt,
        command=command,
        workdir=REPO_ROOT,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
    )
    return StageExecution(
        stage_run_id=stage_run_id,
        component_run_id=component_run_id,
        component_name=component.config.component,
        stage_name=stage_name,
        attempt=attempt,
        command=command,
        workdir=REPO_ROOT,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        tasks_root=component.tasks_root,
        train_result_path=component.train_result_path,
        promote_log_root=component.promote_log_root,
    )


def latest_stage_status(db: CampaignDb, component_run_id: int, stage_name: str) -> str | None:
    row = db.latest_stage_run(component_run_id, stage_name)
    if row is None:
        return None
    return str(row["status"])


def next_stage_for_component(
    *,
    db: CampaignDb,
    component_run_id: int,
    component: ResolvedComponent,
) -> str | None:
    prepare_status = latest_stage_status(db, component_run_id, "prepare")
    train_status = latest_stage_status(db, component_run_id, "train")
    promote_status = latest_stage_status(db, component_run_id, "promote")

    if component.config.prepare.enabled:
        if prepare_status in (None, "failed"):
            return "prepare"
        if prepare_status != "succeeded":
            return None
    elif not component.tasks_root.exists():
        raise FileNotFoundError(
            f"{component.config.component} prepare is disabled but tasks root does not exist: {component.tasks_root}"
        )

    if component.config.train.enabled:
        if train_status in (None, "failed"):
            return "train"
        if train_status != "succeeded":
            return None
    else:
        return None

    if component.config.promote.enabled:
        if promote_status in (None, "failed"):
            return "promote"
        if promote_status != "succeeded":
            return None

    return None


def component_complete(db: CampaignDb, component_run_id: int, component: ResolvedComponent) -> bool:
    if component.config.promote.enabled:
        return latest_stage_status(db, component_run_id, "promote") == "succeeded"
    if component.config.train.enabled:
        return latest_stage_status(db, component_run_id, "train") == "succeeded"
    if component.config.prepare.enabled:
        return latest_stage_status(db, component_run_id, "prepare") == "succeeded"
    return component.tasks_root.exists()


def stage_summary_text(stage_name: str, result_json: dict[str, Any] | None) -> str | None:
    if result_json is None:
        return None
    if stage_name == "train":
        best_val = result_json.get("best_val_config")
        val_metrics = result_json.get("val_metrics") or {}
        heldout_metrics = result_json.get("heldout_metrics") or {}
        score = val_metrics.get("score")
        heldout_score = heldout_metrics.get("score")
        return (
            f"best_val={best_val} val_score={score} heldout_score={heldout_score}"
            if best_val is not None
            else None
        )
    if stage_name == "promote":
        decision = result_json.get("decision") or result_json.get("promotion_decision")
        quick_test = result_json.get("quick_test") or {}
        elapsed_seconds = quick_test.get("elapsed_seconds")
        return (
            f"decision={decision} quick_test_seconds={elapsed_seconds}"
            if decision is not None
            else None
        )
    if stage_name == "prepare":
        train_tasks = result_json.get("train_tasks")
        val_tasks = result_json.get("val_tasks")
        heldout_tasks = result_json.get("heldout_tasks")
        return (
            f"train={train_tasks} val={val_tasks} heldout={heldout_tasks}"
            if train_tasks is not None
            else None
        )
    return None


def finalize_stage(db: CampaignDb, outcome: StageOutcome) -> None:
    db.finish_stage_run(
        outcome.stage_run_id,
        status="succeeded" if outcome.success else "failed",
        exit_code=outcome.exit_code,
        duration_seconds=outcome.duration_seconds,
        result_path=outcome.result_path,
        result_json=outcome.result_json,
        error_message=outcome.error_message,
    )
    for artifact in outcome.artifacts:
        db.add_artifact(
            stage_run_id=outcome.stage_run_id,
            artifact_type=artifact.artifact_type,
            path=artifact.path,
            metadata=artifact.metadata,
        )

    if outcome.success:
        update_kwargs: dict[str, Any] = {"status": "running", "last_error": None}
        if outcome.stage_name == "prepare" and outcome.result_path is not None:
            update_kwargs["latest_prepare_metadata_path"] = outcome.result_path
        elif outcome.stage_name == "train" and outcome.result_path is not None:
            update_kwargs["latest_train_result_path"] = outcome.result_path
        elif outcome.stage_name == "promote" and outcome.result_path is not None:
            update_kwargs["latest_promote_summary_path"] = outcome.result_path
        db.update_component_run(outcome.component_run_id, **update_kwargs)
    else:
        db.update_component_run(
            outcome.component_run_id,
            status="failed",
            last_error=outcome.error_message,
        )


def initialize_campaign_rows(
    db: CampaignDb,
    *,
    campaign_id: int,
    components: list[ResolvedComponent],
) -> dict[str, int]:
    component_run_ids: dict[str, int] = {}
    for component in components:
        component_run_id = db.create_component_run(
            campaign_id=campaign_id,
            component_name=component.config.component,
                artifact_root=component.artifact_root,
                tasks_root=component.tasks_root,
                manifest={
                    "component": component.config.component,
                    "tasks_root": str(component.tasks_root),
                    "prepare": {
                        "enabled": component.config.prepare.enabled,
                        "args": list(component.config.prepare.args),
                        "binary": (
                            str(component.config.prepare.binary)
                            if component.config.prepare.binary is not None
                            else None
                        ),
                    },
                    "train": {
                        "enabled": component.config.train.enabled,
                        "args": list(component.config.train.args),
                        "binary": (
                            str(component.config.train.binary)
                            if component.config.train.binary is not None
                            else None
                        ),
                    },
                    "promote": {
                        "enabled": component.config.promote.enabled,
                        "args": list(component.config.promote.args),
                        "binary": (
                            str(component.config.promote.binary)
                            if component.config.promote.binary is not None
                            else None
                        ),
                    },
                },
            )
        component_run_ids[component.config.component] = component_run_id
        metadata_path = component.tasks_root / "metadata.json"
        if metadata_path.exists():
            db.update_component_run(
                component_run_id,
                latest_prepare_metadata_path=metadata_path,
            )
    return component_run_ids


def run_campaign(
    *,
    db: CampaignDb,
    campaign_id: int,
    campaign: CampaignConfig,
    campaign_root: Path,
) -> int:
    resolved_components = resolve_component_layout(campaign, campaign_root)
    component_run_ids = initialize_campaign_rows(
        db,
        campaign_id=campaign_id,
        components=resolved_components,
    )
    stale = db.mark_running_stages_failed(
        campaign_id,
        error_message="campaign resumed after interruption; rerun stage",
    )
    if stale:
        print(f"[campaign {campaign_id}] marked {stale} stale running stage(s) as failed", flush=True)

    running: dict[Future[StageOutcome], StageExecution] = {}
    fail_fast_triggered = False
    any_failure = False

    with ThreadPoolExecutor(max_workers=campaign.max_workers) as executor:
        while True:
            scheduled_work = False
            if not fail_fast_triggered:
                for component in resolved_components:
                    if any(
                        execution.component_name == component.config.component
                        for execution in running.values()
                    ):
                        continue

                    component_run_id = component_run_ids[component.config.component]
                    if component_complete(db, component_run_id, component):
                        db.update_component_run(component_run_id, status="succeeded", last_error=None)
                        continue

                    try:
                        stage_name = next_stage_for_component(
                            db=db,
                            component_run_id=component_run_id,
                            component=component,
                        )
                    except Exception as exc:
                        any_failure = True
                        fail_fast_triggered = campaign.fail_fast
                        db.update_component_run(
                            component_run_id,
                            status="failed",
                            last_error=str(exc),
                        )
                        continue
                    if stage_name is None:
                        continue

                    execution = create_stage_execution(
                        db=db,
                        component_run_id=component_run_id,
                        component=component,
                        stage_name=stage_name,
                        campaign=campaign,
                    )
                    print(
                        f"[campaign {campaign_id}] start {component.config.component}:{stage_name} "
                        f"attempt={execution.attempt} cmd={format_command(execution.command)}",
                        flush=True,
                    )
                    future = executor.submit(execute_stage, execution)
                    running[future] = execution
                    scheduled_work = True
                    if len(running) >= campaign.max_workers:
                        break

            if not running and not scheduled_work:
                break

            if running:
                completed_futures, _ = wait(running.keys(), return_when=FIRST_COMPLETED)
                for future in completed_futures:
                    execution = running.pop(future)
                    outcome = future.result()
                    finalize_stage(db, outcome)
                    summary = stage_summary_text(outcome.stage_name, outcome.result_json)
                    if outcome.success:
                        print(
                            f"[campaign {campaign_id}] done {outcome.component_name}:{outcome.stage_name} "
                            f"seconds={outcome.duration_seconds:.2f}"
                            + (f" {summary}" if summary else ""),
                            flush=True,
                        )
                    else:
                        any_failure = True
                        if campaign.fail_fast:
                            fail_fast_triggered = True
                        print(
                            f"[campaign {campaign_id}] fail {outcome.component_name}:{outcome.stage_name} "
                            f"seconds={outcome.duration_seconds:.2f} error={outcome.error_message}",
                            flush=True,
                        )

    campaign_status = "failed" if any_failure else "succeeded"
    db.update_campaign_status(campaign_id, campaign_status)
    for component in resolved_components:
        component_run_id = component_run_ids[component.config.component]
        if component_complete(db, component_run_id, component):
            db.update_component_run(component_run_id, status="succeeded", last_error=None)
        elif latest_stage_status(db, component_run_id, "promote") == "failed" or latest_stage_status(
            db, component_run_id, "train"
        ) == "failed" or latest_stage_status(db, component_run_id, "prepare") == "failed":
            db.update_component_run(
                component_run_id,
                status="failed",
                last_error=db.get_component_run(campaign_id, component.config.component)["last_error"],
            )
        else:
            db.update_component_run(component_run_id, status="pending")
    return 0 if campaign_status == "succeeded" else 1


def print_campaign_status(db: CampaignDb, campaign_id: int) -> None:
    campaign = db.get_campaign(campaign_id)
    print(
        f"campaign_id={campaign.id} name={campaign.name} status={campaign.status} "
        f"workers={campaign.max_workers} fail_fast={campaign.fail_fast}"
    )
    print(f"campaign_root={campaign.campaign_root}")
    print(f"manifest={campaign.manifest_path}")
    if campaign.last_error:
        print(f"last_error={campaign.last_error}")
    for component_row in db.list_component_runs(campaign_id):
        component_manifest = decode_json(component_row["manifest_json"]) or {}
        stage_enabled = {
            "prepare": bool((component_manifest.get("prepare") or {}).get("enabled", False)),
            "train": bool((component_manifest.get("train") or {}).get("enabled", False)),
            "promote": bool((component_manifest.get("promote") or {}).get("enabled", False)),
        }
        print(
            f"- {component_row['component_name']}: status={component_row['status']} "
            f"tasks_root={component_row['tasks_root']}"
        )
        stage_rows = {row["stage_name"]: row for row in db.list_stage_runs(int(component_row["id"]))}
        for stage_name in ("prepare", "train", "promote"):
            row = stage_rows.get(stage_name)
            if row is None:
                if not stage_enabled[stage_name]:
                    print(f"  {stage_name}: disabled")
                else:
                    print(f"  {stage_name}: pending")
                continue
            summary = stage_summary_text(stage_name, decode_json(row["result_json"]))
            summary_suffix = f" {summary}" if summary else ""
            print(
                f"  {stage_name}: status={row['status']} attempt={row['attempt']} "
                f"seconds={row['duration_seconds']}{summary_suffix}"
            )


def start_new_campaign(args: argparse.Namespace) -> int:
    if args.fail_fast and args.no_fail_fast:
        raise ValueError("cannot pass both --fail-fast and --no-fail-fast")
    fail_fast_override = None
    if args.fail_fast:
        fail_fast_override = True
    elif args.no_fail_fast:
        fail_fast_override = False

    manifest = load_manifest(
        args.manifest.resolve(),
        db_path_override=args.db,
        campaigns_root_override=args.campaigns_root,
        max_workers_override=args.max_workers,
        fail_fast_override=fail_fast_override,
    )
    manifest_json = manifest_to_jsonable(manifest)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    campaign_root = manifest.campaigns_root / f"{slugify(manifest.name)}_{timestamp}"
    campaign_root.mkdir(parents=True, exist_ok=True)
    manifest_copy_path = campaign_root / "manifest.json"
    write_json(manifest_copy_path, manifest_json)

    db = CampaignDb(manifest.db_path)
    try:
        campaign_id = db.create_campaign(
            name=manifest.name,
            manifest_path=manifest_copy_path,
            campaign_root=campaign_root,
            max_workers=manifest.max_workers,
            fail_fast=manifest.fail_fast,
            metadata=manifest_json,
        )
        write_json(
            campaign_root / "campaign.json",
            {
                "campaign_id": campaign_id,
                "database_path": str(manifest.db_path),
                "campaign_root": str(campaign_root),
            },
        )
        exit_code = run_campaign(
            db=db,
            campaign_id=campaign_id,
            campaign=manifest,
            campaign_root=campaign_root,
        )
        print_campaign_status(db, campaign_id)
        return exit_code
    finally:
        db.close()


def resume_campaign(args: argparse.Namespace) -> int:
    if args.fail_fast and args.no_fail_fast:
        raise ValueError("cannot pass both --fail-fast and --no-fail-fast")
    db = CampaignDb(args.db)
    try:
        campaign_record = db.get_campaign(args.campaign_id)
        manifest_path = Path(campaign_record.manifest_path)
        fail_fast_override = None
        if args.fail_fast:
            fail_fast_override = True
        elif args.no_fail_fast:
            fail_fast_override = False
        manifest = load_manifest(
            manifest_path,
            db_path_override=args.db,
            campaigns_root_override=Path(campaign_record.campaign_root).parent,
            max_workers_override=args.max_workers,
            fail_fast_override=fail_fast_override,
        )
        db.update_campaign_status(args.campaign_id, "running")
        exit_code = run_campaign(
            db=db,
            campaign_id=args.campaign_id,
            campaign=manifest,
            campaign_root=Path(campaign_record.campaign_root),
        )
        print_campaign_status(db, args.campaign_id)
        return exit_code
    finally:
        db.close()


def show_status(args: argparse.Namespace) -> int:
    db = CampaignDb(args.db)
    try:
        if args.campaign_id is not None:
            print_campaign_status(db, args.campaign_id)
            return 0
        latest_id = db.latest_campaign_id()
        if latest_id is None:
            print("no campaigns found")
            return 0
        print_campaign_status(db, latest_id)
        print()
        print("recent campaigns:")
        for row in db.list_campaigns(limit=args.limit):
            print(
                f"- id={row['id']} name={row['name']} status={row['status']} "
                f"created_at={row['created_at']} root={row['campaign_root']}"
            )
        return 0
    finally:
        db.close()


def main() -> int:
    args = parse_args()
    if args.command == "run":
        return start_new_campaign(args)
    if args.command == "resume":
        return resume_campaign(args)
    if args.command == "status":
        return show_status(args)
    raise ValueError(f"unknown command {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
