#!/usr/bin/env python3
"""Validate the public Trinity-parity panel manifest.

This intentionally does not download data. It verifies that the manifest has
enough machine-readable provenance, inputs, thresholds, and required-stage
coverage before any expensive public benchmark run is allowed to claim parity.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = ROOT / "bench" / "trinity_parity" / "public_panel.json"
REQUIRED_TOP_LEVEL = {"name", "status", "defaults", "metrics", "datasets"}
REQUIRED_DEFAULTS = {
    "seq_type",
    "threads",
    "max_memory",
    "min_raptor_vs_trinity_selected_f1",
    "min_reference_full_length_recovery",
    "min_read_representation",
    "max_wall_time_ratio_vs_trinity",
    "max_peak_rss_ratio_vs_trinity",
}
REQUIRED_METRICS = {
    "selected_transcript_reciprocal_f1",
    "reference_full_length_recovery",
    "read_representation",
    "wall_time_seconds",
    "peak_rss_kb",
}
REQUIRED_STAGES = {
    "input_fastq",
    "normalization",
    "inchworm",
    "chrysalis",
    "butterfly",
    "reporting",
}


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def require(condition: bool, failures: list[str], message: str) -> None:
    if not condition:
        failures.append(message)


def nonempty_string(value: object) -> bool:
    return isinstance(value, str) and bool(value.strip())


def string_list(value: object) -> bool:
    return (
        isinstance(value, list)
        and bool(value)
        and all(nonempty_string(item) for item in value)
    )


def validate_thresholds(defaults: dict[str, object], failures: list[str]) -> None:
    for key in REQUIRED_DEFAULTS:
        require(key in defaults, failures, f"defaults missing {key}")
    require(defaults.get("seq_type") == "fq", failures, "defaults.seq_type must be fq")
    require(
        isinstance(defaults.get("threads"), int) and defaults["threads"] > 0,
        failures,
        "defaults.threads must be a positive integer",
    )
    for key in [
        "min_raptor_vs_trinity_selected_f1",
        "min_reference_full_length_recovery",
        "min_read_representation",
    ]:
        value = defaults.get(key)
        require(
            isinstance(value, (int, float)) and 0.0 <= float(value) <= 1.0,
            failures,
            f"defaults.{key} must be in [0,1]",
        )
    for key in ["max_wall_time_ratio_vs_trinity", "max_peak_rss_ratio_vs_trinity"]:
        value = defaults.get(key)
        require(
            isinstance(value, (int, float)) and float(value) >= 1.0,
            failures,
            f"defaults.{key} must be >= 1",
        )


def validate_dataset(dataset: object, seen_ids: set[str], failures: list[str]) -> None:
    require(isinstance(dataset, dict), failures, "dataset entry must be an object")
    if not isinstance(dataset, dict):
        return

    dataset_id = dataset.get("id")
    require(nonempty_string(dataset_id), failures, "dataset missing id")
    if isinstance(dataset_id, str):
        require(dataset_id not in seen_ids, failures, f"duplicate dataset id {dataset_id}")
        seen_ids.add(dataset_id)

    require(
        nonempty_string(dataset.get("organism")),
        failures,
        f"{dataset_id}: missing organism",
    )
    require(
        nonempty_string(dataset.get("description")),
        failures,
        f"{dataset_id}: missing description",
    )
    source = dataset.get("source")
    require(isinstance(source, dict), failures, f"{dataset_id}: source must be an object")
    if isinstance(source, dict):
        for key in ["name", "url", "documentation_url"]:
            require(nonempty_string(source.get(key)), failures, f"{dataset_id}: source missing {key}")

    library = dataset.get("library")
    require(isinstance(library, dict), failures, f"{dataset_id}: library must be an object")
    if isinstance(library, dict):
        require(
            library.get("layout") in {"paired", "single"},
            failures,
            f"{dataset_id}: library.layout must be paired or single",
        )
        ss_lib_type = library.get("ss_lib_type")
        require(
            ss_lib_type in {None, "F", "R", "FR", "RF"},
            failures,
            f"{dataset_id}: invalid ss_lib_type",
        )

    inputs = dataset.get("inputs")
    require(isinstance(inputs, dict), failures, f"{dataset_id}: inputs must be an object")
    if isinstance(inputs, dict):
        require(
            string_list(inputs.get("left")),
            failures,
            f"{dataset_id}: inputs.left must be a nonempty string list",
        )
        if isinstance(library, dict) and library.get("layout") == "paired":
            require(
                string_list(inputs.get("right")),
                failures,
                f"{dataset_id}: paired dataset requires inputs.right",
            )
            if string_list(inputs.get("left")) and string_list(inputs.get("right")):
                require(
                    len(inputs["left"]) == len(inputs["right"]),
                    failures,
                    f"{dataset_id}: left/right input counts differ",
                )

    require(
        string_list(dataset.get("trinity_args")),
        failures,
        f"{dataset_id}: trinity_args must be a nonempty string list",
    )
    stages = dataset.get("required_stages")
    require(string_list(stages), failures, f"{dataset_id}: required_stages must be a nonempty string list")
    if isinstance(stages, list):
        missing = sorted(REQUIRED_STAGES.difference(stages))
        require(not missing, failures, f"{dataset_id}: missing required stages {','.join(missing)}")


def validate_manifest(manifest: dict[str, object]) -> list[str]:
    failures: list[str] = []
    missing = sorted(REQUIRED_TOP_LEVEL.difference(manifest))
    require(not missing, failures, f"manifest missing top-level keys {','.join(missing)}")

    defaults = manifest.get("defaults")
    require(isinstance(defaults, dict), failures, "defaults must be an object")
    if isinstance(defaults, dict):
        validate_thresholds(defaults, failures)

    metrics = manifest.get("metrics")
    require(string_list(metrics), failures, "metrics must be a nonempty string list")
    if isinstance(metrics, list):
        missing_metrics = sorted(REQUIRED_METRICS.difference(metrics))
        require(not missing_metrics, failures, f"metrics missing {','.join(missing_metrics)}")

    datasets = manifest.get("datasets")
    require(
        isinstance(datasets, list) and len(datasets) >= 3,
        failures,
        "public panel needs at least 3 datasets",
    )
    seen_ids: set[str] = set()
    if isinstance(datasets, list):
        for dataset in datasets:
            validate_dataset(dataset, seen_ids, failures)
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    args = parser.parse_args()

    manifest = load_json(args.manifest)
    failures = validate_manifest(manifest)
    report = {
        "manifest": str(args.manifest),
        "dataset_count": len(manifest.get("datasets", []))
        if isinstance(manifest.get("datasets"), list)
        else 0,
        "passed": not failures,
        "failures": failures,
    }
    print(json.dumps(report, indent=2))
    if failures:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
