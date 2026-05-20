#!/usr/bin/env python3
"""Run ready datasets from the public Trinity-parity panel."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import signal
import subprocess
import time
from pathlib import Path

from plan_public_panel import DEFAULT_DATA_ROOT, build_plan
from validate_public_panel import DEFAULT_MANIFEST, load_json, validate_manifest


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT_ROOT = ROOT / "target" / "trinity_parity" / "public_panel"
DEFAULT_REPORT = ROOT / "target" / "trinity_parity" / "public_panel_report.json"
DEFAULT_COVERAGE_SWEEP = [0.90, 0.925, 0.95, 0.975]


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def read_fasta_records(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    if not path.exists():
        return records
    current_name: str | None = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    records[current_name] = "".join(chunks).upper()
                current_name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current_name is not None:
        records[current_name] = "".join(chunks).upper()
    return records


def n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    half = sum(lengths) / 2
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= half:
            return length
    return 0


def fasta_stats(path: Path) -> dict[str, object]:
    records = read_fasta_records(path)
    lengths = sorted((len(seq) for seq in records.values()), reverse=True)
    return {
        "path": str(path),
        "exists": path.exists(),
        "transcript_count": len(lengths),
        "total_bases": sum(lengths),
        "n50": n50(lengths),
        "max_length": lengths[0] if lengths else 0,
        "lengths_top20": lengths[:20],
    }


def best_reciprocal_coverage(query: str, references: list[str]) -> float:
    best = 0.0
    rc_query = reverse_complement(query)
    for reference in references:
        for candidate in (query, rc_query):
            if candidate in reference or reference in candidate:
                coverage = min(len(candidate), len(reference)) / max(
                    len(candidate), len(reference)
                )
                best = max(best, coverage)
    return best


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


def fasta_match_metrics(
    query_fasta: Path,
    reference_fasta: Path,
    min_coverage: float,
) -> dict[str, object]:
    query = read_fasta_records(query_fasta)
    reference = read_fasta_records(reference_fasta)
    reference_sequences = list(reference.values())
    matched = 0
    coverages: list[float] = []
    for sequence in query.values():
        coverage = best_reciprocal_coverage(sequence, reference_sequences)
        coverages.append(round(coverage, 6))
        if coverage >= min_coverage:
            matched += 1
    precision = matched / len(query) if query else 0.0
    return {
        "query_count": len(query),
        "reference_count": len(reference),
        "matched": matched,
        "precision": round(precision, 6),
        "min_match_coverage": min_coverage,
        "best_coverages_top20": sorted(coverages, reverse=True)[:20],
    }


def reciprocal_fasta_metrics(
    raptor_fasta: Path,
    trinity_fasta: Path,
    min_coverage: float,
) -> dict[str, object]:
    raptor_to_trinity = fasta_match_metrics(raptor_fasta, trinity_fasta, min_coverage)
    trinity_to_raptor = fasta_match_metrics(trinity_fasta, raptor_fasta, min_coverage)
    precision = float(raptor_to_trinity["precision"])
    recall = float(trinity_to_raptor["precision"])
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "raptor_to_trinity": raptor_to_trinity,
        "trinity_to_raptor": trinity_to_raptor,
        "precision": round(precision, 6),
        "recall": round(recall, 6),
        "f1": round(f1, 6),
    }


def reciprocal_fasta_threshold_sweep(
    raptor_fasta: Path,
    trinity_fasta: Path,
    thresholds: list[float],
) -> list[dict[str, object]]:
    sweep = []
    for threshold in thresholds:
        metrics = reciprocal_fasta_metrics(raptor_fasta, trinity_fasta, threshold)
        sweep.append(
            {
                "min_match_coverage": threshold,
                "precision": metrics["precision"],
                "recall": metrics["recall"],
                "f1": metrics["f1"],
                "raptor_to_trinity_matched": metrics["raptor_to_trinity"]["matched"],
                "trinity_to_raptor_matched": metrics["trinity_to_raptor"]["matched"],
            }
        )
    return sweep


def command_resource_usage(path: Path | None) -> dict[str, object]:
    if path is None or not path.exists():
        return {
            "available": False,
            "source": None,
            "max_rss_kb": None,
            "user_seconds": None,
            "system_seconds": None,
        }
    values: dict[str, object] = {
        "available": True,
        "source": "gnu_time",
        "max_rss_kb": None,
        "user_seconds": None,
        "system_seconds": None,
    }
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if ":" not in line:
            continue
        key, raw_value = line.split(":", 1)
        value = raw_value.strip()
        if key == "Maximum resident set size (kbytes)":
            values["max_rss_kb"] = int(value)
        elif key == "User time (seconds)":
            values["user_seconds"] = float(value)
        elif key == "System time (seconds)":
            values["system_seconds"] = float(value)
    return values


def run_command(
    command: list[str],
    cwd: Path,
    work_dir: Path,
    timeout_seconds: int | None,
) -> dict[str, object]:
    work_dir.mkdir(parents=True, exist_ok=True)
    time_bin = Path("/usr/bin/time")
    time_output = work_dir / "time.txt"
    measured = command
    if time_bin.exists():
        measured = [str(time_bin), "-v", "-o", str(time_output), *command]
    start = time.monotonic()
    process = subprocess.Popen(
        measured,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
    )
    timed_out = False
    try:
        stdout, stderr = process.communicate(timeout=timeout_seconds)
    except subprocess.TimeoutExpired:
        timed_out = True
        os.killpg(process.pid, signal.SIGTERM)
        try:
            stdout, stderr = process.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            stdout, stderr = process.communicate()
    elapsed = time.monotonic() - start
    return {
        "command": command,
        "exit_code": process.returncode,
        "timed_out": timed_out,
        "timeout_seconds": timeout_seconds,
        "elapsed_seconds": round(elapsed, 6),
        "stdout": stdout,
        "stderr": stderr,
        "resource_usage": command_resource_usage(time_output if time_bin.exists() else None),
    }


def raptor_fasta_path(dataset_id: str, out_root: Path) -> Path:
    return out_root / dataset_id / "raptor" / "raptor_trinity.fasta.gz"


def trinity_fasta_path(dataset_id: str, out_root: Path) -> Path:
    dataset_root = out_root / dataset_id
    candidates = [
        dataset_root / "trinity" / "Trinity.fasta",
        dataset_root / "trinity.Trinity.fasta",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def trinity_inchworm_fasta_path(dataset_id: str, out_root: Path) -> Path:
    return out_root / dataset_id / "trinity" / "inchworm.DS.fa"


def selected_dataset_plans(
    plan: dict[str, object],
    dataset_ids: list[str],
    ready_only: bool,
) -> tuple[list[dict[str, object]], list[str]]:
    datasets = [dataset for dataset in plan.get("datasets", []) if isinstance(dataset, dict)]
    selected = set(dataset_ids)
    if selected:
        datasets = [dataset for dataset in datasets if dataset.get("id") in selected]
        missing_ids = selected.difference(str(dataset.get("id")) for dataset in datasets)
    else:
        missing_ids = set()
    if ready_only:
        datasets = [dataset for dataset in datasets if dataset.get("ready")]
    return datasets, [f"unknown dataset id: {dataset_id}" for dataset_id in sorted(missing_ids)]


def run_dataset(
    dataset: dict[str, object],
    out_root: Path,
    run_raptor: bool,
    run_trinity: bool,
    min_match_coverage: float,
    coverage_sweep: list[float],
    dry_run: bool,
    timeout_seconds: int | None,
    min_fasta_f1: float | None,
) -> dict[str, object]:
    dataset_id = str(dataset["id"])
    dataset_out = out_root / dataset_id
    result: dict[str, object] = {
        "id": dataset_id,
        "ready": dataset.get("ready"),
        "missing_inputs": dataset.get("missing_inputs", []),
        "raptor": None,
        "trinity": None,
        "metrics": {},
    }
    if dry_run:
        result["raptor"] = {"command": dataset["raptor_command"], "dry_run": True}
        result["trinity"] = {"command": dataset["trinity_command"], "dry_run": True}
        return result
    if run_raptor:
        result["raptor"] = run_command(
            [str(part) for part in dataset["raptor_command"]],
            ROOT,
            dataset_out / "raptor_command",
            timeout_seconds,
        )
    if run_trinity:
        result["trinity"] = run_command(
            [str(part) for part in dataset["trinity_command"]],
            ROOT,
            dataset_out / "trinity_command",
            timeout_seconds,
        )

    raptor_fasta = raptor_fasta_path(dataset_id, out_root)
    trinity_fasta = trinity_fasta_path(dataset_id, out_root)
    trinity_inchworm_fasta = trinity_inchworm_fasta_path(dataset_id, out_root)
    metrics = {
        "raptor_fasta": fasta_stats(raptor_fasta),
        "trinity_fasta": fasta_stats(trinity_fasta),
        "trinity_inchworm_fasta": fasta_stats(trinity_inchworm_fasta),
    }
    if raptor_fasta.exists() and trinity_fasta.exists():
        metrics["raptor_trinity_fasta_match"] = reciprocal_fasta_metrics(
            raptor_fasta,
            trinity_fasta,
            min_match_coverage,
        )
        metrics["raptor_trinity_fasta_threshold_sweep"] = reciprocal_fasta_threshold_sweep(
            raptor_fasta,
            trinity_fasta,
            coverage_sweep,
        )
        metrics["min_raptor_vs_trinity_selected_f1"] = min_fasta_f1
    if raptor_fasta.exists() and trinity_inchworm_fasta.exists():
        metrics["raptor_trinity_inchworm_fasta_match"] = reciprocal_fasta_metrics(
            raptor_fasta,
            trinity_inchworm_fasta,
            min_match_coverage,
        )
        metrics["raptor_trinity_inchworm_fasta_threshold_sweep"] = (
            reciprocal_fasta_threshold_sweep(
                raptor_fasta,
                trinity_inchworm_fasta,
                coverage_sweep,
            )
        )
    result["metrics"] = metrics
    return result


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--trinity-bin", default="scripts/trinity_docker.sh")
    parser.add_argument("--dataset", action="append", default=[])
    parser.add_argument("--include-unready", action="store_true")
    parser.add_argument("--skip-raptor", action="store_true")
    parser.add_argument("--run-trinity", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--min-match-coverage", type=float, default=0.95)
    parser.add_argument(
        "--coverage-sweep",
        default=",".join(str(value) for value in DEFAULT_COVERAGE_SWEEP),
        help="Comma-separated reciprocal FASTA coverage thresholds to report as diagnostics.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=None,
        help="Terminate each Raptor/Trinity command after this many seconds.",
    )
    args = parser.parse_args()

    manifest = load_json(args.manifest)
    failures = validate_manifest(manifest)
    if failures:
        print(json.dumps({"passed": False, "failures": failures}, indent=2))
        return 1

    out_root = args.out_root.resolve()
    try:
        coverage_sweep = [
            float(value)
            for value in args.coverage_sweep.split(",")
            if value.strip()
        ]
    except ValueError as err:
        print(json.dumps({"passed": False, "failures": [str(err)]}, indent=2))
        return 1
    plan = build_plan(
        manifest,
        args.data_root.resolve(),
        out_root,
        args.trinity_bin,
    )
    datasets, selection_failures = selected_dataset_plans(
        plan,
        args.dataset,
        ready_only=not args.include_unready,
    )
    results = [
        run_dataset(
            dataset,
            out_root,
            run_raptor=not args.skip_raptor,
            run_trinity=args.run_trinity,
            min_match_coverage=args.min_match_coverage,
            coverage_sweep=coverage_sweep,
            dry_run=args.dry_run,
            timeout_seconds=args.timeout_seconds,
            min_fasta_f1=plan.get("defaults", {}).get(
                "min_raptor_vs_trinity_selected_f1"
            ),
        )
        for dataset in datasets
    ]
    failures.extend(selection_failures)
    for result in results:
        for command_key in ["raptor", "trinity"]:
            command_result = result.get(command_key)
            if isinstance(command_result, dict) and command_result.get("exit_code") not in {
                None,
                0,
            }:
                if command_result.get("timed_out"):
                    failures.append(f"{result['id']}: {command_key} command timed out")
                else:
                    failures.append(f"{result['id']}: {command_key} command failed")
        metrics = result.get("metrics", {})
        if isinstance(metrics, dict):
            fasta_match = metrics.get("raptor_trinity_fasta_match")
            min_f1 = metrics.get("min_raptor_vs_trinity_selected_f1")
            if isinstance(fasta_match, dict) and isinstance(min_f1, (int, float)):
                f1 = float(fasta_match.get("f1", 0.0))
                if f1 < float(min_f1):
                    failures.append(
                        f"{result['id']}: Raptor-vs-Trinity FASTA F1 {f1} < {float(min_f1)}"
                    )
    report = {
        "manifest": plan["manifest"],
        "status": plan["status"],
        "data_root": plan["data_root"],
        "out_root": str(out_root),
        "dry_run": args.dry_run,
        "run_trinity": args.run_trinity,
        "dataset_count": len(results),
        "ready_dataset_count": plan["ready_dataset_count"],
        "passed": not failures,
        "failures": failures,
        "datasets": results,
    }
    write_json(args.report, report)
    print(json.dumps(report, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
