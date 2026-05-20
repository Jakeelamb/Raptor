#!/usr/bin/env python3
"""Build a reproducible run plan for the public Trinity-parity panel."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from validate_public_panel import DEFAULT_MANIFEST, load_json, validate_manifest


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA_ROOT = ROOT / "data" / "trinity_public_panel"
DEFAULT_OUT = ROOT / "target" / "trinity_parity" / "public_panel_run_plan.json"


def as_string_list(value: object) -> list[str]:
    if isinstance(value, list):
        return [str(item) for item in value]
    return []


def dataset_root(data_root: Path, dataset_id: str) -> Path:
    return data_root / dataset_id


def resolve_inputs(dataset: dict[str, object], data_root: Path) -> dict[str, object]:
    dataset_id = str(dataset["id"])
    root = dataset_root(data_root, dataset_id)
    inputs = dataset.get("inputs", {})
    if not isinstance(inputs, dict):
        return {"root": str(root), "files": [], "missing": []}

    files: list[dict[str, object]] = []
    for role, raw_paths in inputs.items():
        paths = raw_paths if isinstance(raw_paths, list) else [raw_paths]
        for raw_path in paths:
            relative_path = str(raw_path)
            absolute_path = root / relative_path
            files.append(
                {
                    "role": role,
                    "relative_path": relative_path,
                    "path": str(absolute_path),
                    "exists": absolute_path.exists(),
                    "bytes": absolute_path.stat().st_size
                    if absolute_path.exists()
                    else None,
                }
            )
    return {
        "root": str(root),
        "files": files,
        "missing": [entry["path"] for entry in files if not entry["exists"]],
    }


def comma_paths(resolved_inputs: dict[str, object], role: str) -> str:
    files = resolved_inputs.get("files", [])
    if not isinstance(files, list):
        return ""
    paths = [
        str(entry["path"])
        for entry in files
        if isinstance(entry, dict) and entry.get("role") == role
    ]
    return ",".join(paths)


def raptor_command(
    dataset: dict[str, object],
    resolved_inputs: dict[str, object],
    out_root: Path,
) -> list[str]:
    dataset_id = str(dataset["id"])
    output_dir = out_root / dataset_id / "raptor"
    command = [
        "cargo",
        "run",
        "--release",
        "--",
        "trinity",
        "--input1",
        comma_paths(resolved_inputs, "left"),
        "--output-dir",
        str(output_dir),
    ]
    right = comma_paths(resolved_inputs, "right")
    if right:
        command.extend(["--input2", right])
    library = dataset.get("library", {})
    ss_lib_type = library.get("ss_lib_type") if isinstance(library, dict) else None
    if ss_lib_type:
        command.extend(["--SS_lib_type", str(ss_lib_type)])
    return command


def trinity_command(
    dataset: dict[str, object],
    resolved_inputs: dict[str, object],
    out_root: Path,
    trinity_bin: str,
) -> list[str]:
    dataset_id = str(dataset["id"])
    output_dir = out_root / dataset_id / "trinity"
    command = [
        trinity_bin,
        "--left",
        comma_paths(resolved_inputs, "left"),
        "--output",
        str(output_dir),
    ]
    right = comma_paths(resolved_inputs, "right")
    if right:
        command.extend(["--right", right])
    command.extend(as_string_list(dataset.get("trinity_args")))
    return command


def dataset_plan(
    dataset: dict[str, object],
    data_root: Path,
    out_root: Path,
    trinity_bin: str,
) -> dict[str, object]:
    resolved_inputs = resolve_inputs(dataset, data_root)
    missing = resolved_inputs["missing"]
    return {
        "id": dataset["id"],
        "organism": dataset.get("organism"),
        "source": dataset.get("source"),
        "description": dataset.get("description"),
        "data_root": resolved_inputs["root"],
        "ready": not missing,
        "missing_inputs": missing,
        "inputs": resolved_inputs["files"],
        "required_stages": dataset.get("required_stages", []),
        "raptor_command": raptor_command(dataset, resolved_inputs, out_root),
        "trinity_command": trinity_command(
            dataset,
            resolved_inputs,
            out_root,
            trinity_bin,
        ),
    }


def build_plan(
    manifest: dict[str, object],
    data_root: Path,
    out_root: Path,
    trinity_bin: str,
) -> dict[str, object]:
    datasets = manifest.get("datasets", [])
    dataset_plans = [
        dataset_plan(dataset, data_root, out_root, trinity_bin)
        for dataset in datasets
        if isinstance(dataset, dict)
    ]
    missing = [
        f"{plan['id']}: {path}"
        for plan in dataset_plans
        for path in plan["missing_inputs"]
    ]
    return {
        "manifest": manifest.get("name"),
        "status": manifest.get("status"),
        "data_root": str(data_root),
        "out_root": str(out_root),
        "dataset_count": len(dataset_plans),
        "ready_dataset_count": sum(1 for plan in dataset_plans if plan["ready"]),
        "ready": not missing,
        "missing_inputs": missing,
        "defaults": manifest.get("defaults", {}),
        "metrics": manifest.get("metrics", []),
        "datasets": dataset_plans,
    }


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument(
        "--out-root",
        type=Path,
        default=ROOT / "target" / "trinity_parity" / "public_panel",
    )
    parser.add_argument("--report", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--trinity-bin", default="scripts/trinity_docker.sh")
    parser.add_argument(
        "--require-data",
        action="store_true",
        help="Exit nonzero when any declared input file is missing.",
    )
    args = parser.parse_args()

    manifest = load_json(args.manifest)
    failures = validate_manifest(manifest)
    if failures:
        print(json.dumps({"passed": False, "failures": failures}, indent=2))
        return 1

    plan = build_plan(
        manifest,
        args.data_root.resolve(),
        args.out_root.resolve(),
        args.trinity_bin,
    )
    write_json(args.report, plan)
    print(json.dumps(plan, indent=2))
    if args.require_data and plan["missing_inputs"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
