#!/usr/bin/env python3
"""Materialize declared public Trinity-parity datasets into the local data root."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from plan_public_panel import DEFAULT_DATA_ROOT, build_plan
from validate_public_panel import DEFAULT_MANIFEST, load_json, validate_manifest


def input_relative_paths(dataset: dict[str, object]) -> list[str]:
    inputs = dataset.get("inputs", {})
    if not isinstance(inputs, dict):
        return []
    paths: list[str] = []
    for raw_paths in inputs.values():
        values = raw_paths if isinstance(raw_paths, list) else [raw_paths]
        paths.extend(str(value) for value in values)
    return paths


def copy_required_inputs(
    source_root: Path,
    target_root: Path,
    relative_paths: list[str],
    force: bool,
) -> list[str]:
    copied: list[str] = []
    for relative_path in relative_paths:
        source = source_root / relative_path
        target = target_root / relative_path
        if not source.exists():
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and not force:
            continue
        shutil.copy2(source, target)
        copied.append(str(target))
    return copied


def run_git(command: list[str], cwd: Path | None = None) -> None:
    completed = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "command failed: "
            + " ".join(command)
            + "\nstdout:\n"
            + completed.stdout
            + "\nstderr:\n"
            + completed.stderr
        )


def materialize_git_sparse(
    dataset: dict[str, object],
    materialization: dict[str, object],
    target_root: Path,
    force: bool,
    dry_run: bool,
) -> dict[str, object]:
    repo = str(materialization["repo"])
    ref = str(materialization["ref"])
    sparse_paths = [str(path) for path in materialization.get("paths", [])]
    required_paths = input_relative_paths(dataset)
    if dry_run:
        return {
            "type": "git_sparse",
            "repo": repo,
            "ref": ref,
            "paths": sparse_paths,
            "copied": [],
            "dry_run": True,
        }

    with tempfile.TemporaryDirectory(prefix="raptor-public-panel-") as temp_name:
        temp_root = Path(temp_name) / "source"
        run_git(
            [
                "git",
                "clone",
                "--quiet",
                "--depth",
                "1",
                "--filter=blob:none",
                "--sparse",
                repo,
                str(temp_root),
            ]
        )
        run_git(["git", "fetch", "--quiet", "--depth", "1", "origin", ref], cwd=temp_root)
        run_git(["git", "checkout", "--quiet", ref], cwd=temp_root)
        if sparse_paths:
            run_git(["git", "sparse-checkout", "set", *sparse_paths], cwd=temp_root)
        copied = copy_required_inputs(temp_root, target_root, required_paths, force)
    return {
        "type": "git_sparse",
        "repo": repo,
        "ref": ref,
        "paths": sparse_paths,
        "copied": copied,
        "dry_run": False,
    }


def materialize_dataset(
    dataset: dict[str, object],
    data_root: Path,
    force: bool,
    dry_run: bool,
) -> dict[str, object]:
    dataset_id = str(dataset["id"])
    target_root = data_root / dataset_id
    materialization = dataset.get("materialization", {})
    if not isinstance(materialization, dict):
        materialization = {"type": "manual", "note": "No materialization recipe declared."}
    materialization_type = materialization.get("type")
    if materialization_type == "git_sparse":
        result = materialize_git_sparse(
            dataset,
            materialization,
            target_root,
            force,
            dry_run,
        )
    else:
        result = {
            "type": materialization_type or "manual",
            "note": materialization.get("note", "Manual staging required."),
            "copied": [],
            "dry_run": dry_run,
        }
    return {
        "id": dataset_id,
        "target_root": str(target_root),
        "materialization": result,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument(
        "--dataset",
        action="append",
        default=[],
        help="Materialize only this dataset id. Can be repeated.",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing files.")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    manifest = load_json(args.manifest)
    failures = validate_manifest(manifest)
    if failures:
        print(json.dumps({"passed": False, "failures": failures}, indent=2))
        return 1

    datasets = [dataset for dataset in manifest["datasets"] if isinstance(dataset, dict)]
    selected = set(args.dataset)
    if selected:
        datasets = [dataset for dataset in datasets if dataset.get("id") in selected]
        missing_ids = selected.difference(str(dataset.get("id")) for dataset in datasets)
        if missing_ids:
            print(
                json.dumps(
                    {
                        "passed": False,
                        "failures": [f"unknown dataset id: {dataset_id}" for dataset_id in sorted(missing_ids)],
                    },
                    indent=2,
                )
            )
            return 1

    data_root = args.data_root.resolve()
    if not args.dry_run:
        data_root.mkdir(parents=True, exist_ok=True)

    results = []
    errors = []
    for dataset in datasets:
        try:
            results.append(
                materialize_dataset(dataset, data_root, args.force, args.dry_run)
            )
        except Exception as exc:  # noqa: BLE001 - report all materialization failures.
            errors.append(f"{dataset.get('id')}: {exc}")

    plan = build_plan(
        manifest,
        data_root,
        Path("target/trinity_parity/public_panel").resolve(),
        "scripts/trinity_docker.sh",
    )
    report = {
        "data_root": str(data_root),
        "dry_run": args.dry_run,
        "results": results,
        "ready_dataset_count": plan["ready_dataset_count"],
        "missing_inputs": plan["missing_inputs"],
        "passed": not errors,
        "errors": errors,
    }
    print(json.dumps(report, indent=2))
    if errors:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
