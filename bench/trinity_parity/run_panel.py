#!/usr/bin/env python3
"""Run the deterministic Raptor Trinity-parity candidate panel."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PANEL = ROOT / "bench" / "trinity_parity" / "panel.json"
DEFAULT_OUT = ROOT / "target" / "trinity_parity" / "candidate_panel"
FIXTURE_RUNNER = ROOT / "bench" / "trinity_parity" / "run_tiny_fixture.py"


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def run_command(command: list[str], cwd: Path) -> dict[str, object]:
    time_bin = Path("/usr/bin/time")
    if not time_bin.exists():
        return run_command_with_procfs_resource_poll(command, cwd)

    with tempfile.NamedTemporaryFile(prefix="raptor-panel-time-", delete=False) as handle:
        time_output_path = Path(handle.name)
    measured_command = [
        str(time_bin),
        "-v",
        "-o",
        str(time_output_path),
        *command,
    ]

    start = time.monotonic()
    gpu_samples = [sample_gpu_usage()]
    completed = subprocess.run(
        measured_command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    gpu_samples.append(sample_gpu_usage())
    result = {
        "command": command,
        "exit_code": completed.returncode,
        "elapsed_seconds": round(time.monotonic() - start, 6),
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "resource_usage": command_resource_usage(time_output_path),
        "gpu_usage": summarize_gpu_samples(gpu_samples),
    }
    time_output_path.unlink(missing_ok=True)
    return result


def run_command_with_procfs_resource_poll(
    command: list[str], cwd: Path
) -> dict[str, object]:
    with tempfile.NamedTemporaryFile(prefix="raptor-panel-stdout-", delete=False) as stdout_handle:
        stdout_path = Path(stdout_handle.name)
    with tempfile.NamedTemporaryFile(prefix="raptor-panel-stderr-", delete=False) as stderr_handle:
        stderr_path = Path(stderr_handle.name)

    start = time.monotonic()
    max_rss_kb = 0
    gpu_samples = [sample_gpu_usage()]
    with stdout_path.open("w", encoding="utf-8") as stdout_file, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr_file:
        process = subprocess.Popen(
            command,
            cwd=cwd,
            text=True,
            stdout=stdout_file,
            stderr=stderr_file,
            start_new_session=True,
        )
        next_gpu_sample = time.monotonic() + 0.2
        while process.poll() is None:
            max_rss_kb = max(max_rss_kb, process_group_rss_kb(process.pid))
            now = time.monotonic()
            if now >= next_gpu_sample:
                gpu_samples.append(sample_gpu_usage())
                next_gpu_sample = now + 0.2
            time.sleep(0.01)
        max_rss_kb = max(max_rss_kb, process_group_rss_kb(process.pid))
        gpu_samples.append(sample_gpu_usage())

    elapsed = time.monotonic() - start
    stdout = stdout_path.read_text(encoding="utf-8", errors="replace")
    stderr = stderr_path.read_text(encoding="utf-8", errors="replace")
    stdout_path.unlink(missing_ok=True)
    stderr_path.unlink(missing_ok=True)
    return {
        "command": command,
        "exit_code": process.returncode,
        "elapsed_seconds": round(elapsed, 6),
        "stdout": stdout,
        "stderr": stderr,
        "resource_usage": {
            "available": max_rss_kb > 0,
            "source": "procfs_process_group_poll",
            "max_rss_kb": max_rss_kb if max_rss_kb > 0 else None,
            "user_seconds": None,
            "system_seconds": None,
        },
        "gpu_usage": summarize_gpu_samples(gpu_samples),
    }


def process_group_rss_kb(process_group_id: int) -> int:
    total = 0
    proc_root = Path("/proc")
    for status_path in proc_root.glob("[0-9]*/status"):
        try:
            status = status_path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        current_group = None
        current_rss = 0
        for line in status.splitlines():
            if line.startswith("NSpgid:"):
                fields = line.split()
                if fields:
                    current_group = int(fields[-1])
            elif line.startswith("VmRSS:"):
                fields = line.split()
                if len(fields) >= 2:
                    current_rss = int(fields[1])
        if current_group == process_group_id:
            total += current_rss
    return total


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


def sample_gpu_usage() -> dict[str, object]:
    nvidia_smi = shutil.which("nvidia-smi")
    if not nvidia_smi:
        return {"available": False, "source": None, "gpus": []}
    command = [
        nvidia_smi,
        "--query-gpu=index,name,memory.used,utilization.gpu,power.draw",
        "--format=csv,noheader,nounits",
    ]
    completed = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        return {
            "available": False,
            "source": "nvidia-smi",
            "error": completed.stderr.strip(),
            "gpus": [],
        }

    gpus: list[dict[str, object]] = []
    for line in completed.stdout.splitlines():
        parts = [part.strip() for part in line.split(",")]
        if len(parts) < 5:
            continue
        gpus.append(
            {
                "index": int(parts[0]),
                "name": parts[1],
                "memory_used_mib": parse_optional_int(parts[2]),
                "utilization_gpu_percent": parse_optional_int(parts[3]),
                "power_draw_watts": parse_optional_float(parts[4]),
            }
        )
    return {"available": bool(gpus), "source": "nvidia-smi", "gpus": gpus}


def summarize_gpu_samples(samples: list[dict[str, object]]) -> dict[str, object]:
    available_samples = [sample for sample in samples if sample.get("available")]
    if not available_samples:
        return {"available": False, "source": None, "sample_count": len(samples), "gpus": []}

    by_index: dict[int, dict[str, object]] = {}
    for sample in available_samples:
        for gpu in sample.get("gpus", []):
            if not isinstance(gpu, dict):
                continue
            index = int(gpu["index"])
            summary = by_index.setdefault(
                index,
                {
                    "index": index,
                    "name": gpu.get("name"),
                    "max_memory_used_mib": 0,
                    "max_utilization_gpu_percent": 0,
                    "max_power_draw_watts": 0.0,
                },
            )
            memory = gpu.get("memory_used_mib")
            utilization = gpu.get("utilization_gpu_percent")
            power = gpu.get("power_draw_watts")
            if isinstance(memory, int):
                summary["max_memory_used_mib"] = max(
                    int(summary["max_memory_used_mib"]), memory
                )
            if isinstance(utilization, int):
                summary["max_utilization_gpu_percent"] = max(
                    int(summary["max_utilization_gpu_percent"]), utilization
                )
            if isinstance(power, float):
                summary["max_power_draw_watts"] = max(
                    float(summary["max_power_draw_watts"]), power
                )

    return {
        "available": True,
        "source": "nvidia-smi",
        "sample_count": len(samples),
        "gpus": [by_index[index] for index in sorted(by_index)],
    }


def parse_optional_int(value: str) -> int | None:
    try:
        return int(value)
    except ValueError:
        return None


def parse_optional_float(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None


def first_gpu_metric(gpu_usage: dict[str, object], key: str) -> object:
    gpus = gpu_usage.get("gpus", [])
    if not isinstance(gpus, list) or not gpus:
        return None
    first = gpus[0]
    if not isinstance(first, dict):
        return None
    return first.get(key)


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def directory_footprint(path: Path) -> dict[str, object]:
    file_count = 0
    total_bytes = 0
    if not path.exists():
        return {"exists": False, "file_count": 0, "total_bytes": 0}
    for child in path.rglob("*"):
        if child.is_file():
            file_count += 1
            total_bytes += child.stat().st_size
    return {"exists": True, "file_count": file_count, "total_bytes": total_bytes}


def fixture_report_path(out_dir: Path, inserts: list[int]) -> Path:
    if len(inserts) == 1:
        return out_dir / "report.json"
    return out_dir / "insert_sweep_report.json"


def fixture_summary(report_path: Path) -> dict[str, object]:
    payload = load_json(report_path)
    if "summaries" in payload:
        summaries = payload["summaries"]
        failures = payload.get("failures", [])
        return {
            "report_path": str(report_path),
            "insert_sweep": payload.get("insert_sweep", []),
            "failures": failures,
            "passed": not failures,
            "runs": summaries,
        }
    return {
        "report_path": str(report_path),
        "insert_sweep": [payload["fixture"]["insert"]],
        "failures": [],
        "passed": True,
        "runs": [summarize_single_fixture_report(payload)],
    }


def summarize_single_fixture_report(payload: dict[str, object]) -> dict[str, object]:
    raptor_workflow = payload.get("raptor_workflow") or {}
    raptor_workflow_gpu = payload.get("raptor_workflow_gpu") or {}
    raptor_workflow_single = payload.get("raptor_workflow_single") or {}
    raptor_workflow_samples_file = payload.get("raptor_workflow_samples_file") or {}
    malformed_fastq_check = payload.get("malformed_fastq_check") or {}
    workflow_metrics = raptor_workflow.get("metrics", {})
    gpu_workflow_metrics = raptor_workflow_gpu.get("metrics", {})
    single_workflow_metrics = raptor_workflow_single.get("metrics", {})
    samples_workflow_metrics = raptor_workflow_samples_file.get("metrics", {})
    malformed_metrics = malformed_fastq_check.get("metrics", {})
    trinity = payload.get("trinity") or {}
    trinity_result = trinity.get("result", {})
    trinity_metrics = trinity_result.get("metrics", {})
    workflow_resources = raptor_workflow.get("resource_usage", {})
    gpu_workflow_resources = raptor_workflow_gpu.get("resource_usage", {})
    single_workflow_resources = raptor_workflow_single.get("resource_usage", {})
    samples_workflow_resources = raptor_workflow_samples_file.get("resource_usage", {})
    trinity_resources = trinity_result.get("resource_usage", {})
    workflow_gpu = raptor_workflow.get("gpu_usage", {})
    gpu_workflow_gpu = raptor_workflow_gpu.get("gpu_usage", {})
    trinity_gpu = trinity_result.get("gpu_usage", {})
    return {
        "insert": payload["fixture"]["insert"],
        "paired_end_pairs": payload["fixture"]["paired_end_pairs"],
        "report_path": payload.get("report_path"),
        "raptor_workflow_exit_code": raptor_workflow.get("exit_code"),
        "raptor_workflow_elapsed_seconds": raptor_workflow.get("elapsed_seconds"),
        "raptor_workflow_max_rss_kb": workflow_resources.get("max_rss_kb"),
        "raptor_workflow_user_seconds": workflow_resources.get("user_seconds"),
        "raptor_workflow_system_seconds": workflow_resources.get("system_seconds"),
        "raptor_workflow_gpu_available": workflow_gpu.get("available"),
        "raptor_workflow_gpu_max_memory_used_mib": first_gpu_metric(
            workflow_gpu, "max_memory_used_mib"
        ),
        "raptor_workflow_gpu_max_utilization_percent": first_gpu_metric(
            workflow_gpu, "max_utilization_gpu_percent"
        ),
        "raptor_workflow_gpu_max_power_draw_watts": first_gpu_metric(
            workflow_gpu, "max_power_draw_watts"
        ),
        "raptor_gpu_workflow_exit_code": raptor_workflow_gpu.get("exit_code"),
        "raptor_gpu_workflow_elapsed_seconds": raptor_workflow_gpu.get("elapsed_seconds"),
        "raptor_gpu_workflow_max_rss_kb": gpu_workflow_resources.get("max_rss_kb"),
        "raptor_gpu_workflow_gpu_available": gpu_workflow_gpu.get("available"),
        "raptor_gpu_workflow_gpu_max_memory_used_mib": first_gpu_metric(
            gpu_workflow_gpu, "max_memory_used_mib"
        ),
        "raptor_gpu_workflow_gpu_max_utilization_percent": first_gpu_metric(
            gpu_workflow_gpu, "max_utilization_gpu_percent"
        ),
        "raptor_gpu_workflow_gpu_max_power_draw_watts": first_gpu_metric(
            gpu_workflow_gpu, "max_power_draw_watts"
        ),
        "raptor_gpu_workflow_selected_cpu_match_precision": gpu_workflow_metrics.get(
            "cpu_selected_isoform_match", {}
        ).get("precision"),
        "raptor_gpu_workflow_selected_cpu_match_f1": gpu_workflow_metrics.get(
            "cpu_selected_isoform_match", {}
        ).get("f1"),
        "raptor_gpu_workflow_selected_isoform_lengths": gpu_workflow_metrics.get(
            "component_selected_isoform_lengths"
        ),
        "raptor_single_workflow_exit_code": raptor_workflow_single.get("exit_code"),
        "raptor_single_workflow_elapsed_seconds": raptor_workflow_single.get(
            "elapsed_seconds"
        ),
        "raptor_single_workflow_max_rss_kb": single_workflow_resources.get("max_rss_kb"),
        "raptor_single_workflow_input_mode": single_workflow_metrics.get("input_mode"),
        "raptor_single_workflow_selected_precision": single_workflow_metrics.get(
            "component_selected_truth_precision", {}
        ).get("precision"),
        "raptor_single_workflow_selected_f1": single_workflow_metrics.get(
            "component_selected_truth_precision", {}
        ).get("f1"),
        "raptor_single_workflow_selected_isoform_lengths": single_workflow_metrics.get(
            "component_selected_isoform_lengths"
        ),
        "raptor_samples_file_workflow_exit_code": raptor_workflow_samples_file.get("exit_code"),
        "raptor_samples_file_workflow_elapsed_seconds": raptor_workflow_samples_file.get(
            "elapsed_seconds"
        ),
        "raptor_samples_file_workflow_max_rss_kb": samples_workflow_resources.get("max_rss_kb"),
        "raptor_samples_file_workflow_input_mode": samples_workflow_metrics.get("input_mode"),
        "raptor_samples_file_workflow_sample_count": samples_workflow_metrics.get("sample_count"),
        "raptor_samples_file_workflow_selected_precision": samples_workflow_metrics.get(
            "component_selected_truth_precision", {}
        ).get("precision"),
        "raptor_samples_file_workflow_selected_f1": samples_workflow_metrics.get(
            "component_selected_truth_precision", {}
        ).get("f1"),
        "raptor_samples_file_workflow_selected_isoform_lengths": samples_workflow_metrics.get(
            "component_selected_isoform_lengths"
        ),
        "malformed_fastq_exit_code": malformed_fastq_check.get("exit_code"),
        "malformed_fastq_failed": malformed_metrics.get("failed"),
        "malformed_fastq_stderr_mentions_fastq": malformed_metrics.get("stderr_mentions_fastq"),
        "malformed_fastq_stderr_mentions_quality": malformed_metrics.get(
            "stderr_mentions_quality"
        ),
        "raptor_workflow_output_fasta_bytes": workflow_metrics.get("output_fasta_bytes"),
        "raptor_workflow_output_dir_bytes": workflow_metrics.get(
            "output_dir_footprint", {}
        ).get("total_bytes"),
        "raptor_workflow_output_file_count": workflow_metrics.get(
            "output_dir_footprint", {}
        ).get("file_count"),
        "workflow_lengths": workflow_metrics.get("lengths"),
        "workflow_n50": workflow_metrics.get("n50"),
        "workflow_component_count": workflow_metrics.get("component_count"),
        "workflow_component_graph_count": workflow_metrics.get("component_graph_count"),
        "workflow_component_selected_isoform_lengths": workflow_metrics.get(
            "component_selected_isoform_lengths"
        ),
        "workflow_component_selected_precision": workflow_metrics.get(
            "component_selected_truth_precision", {}
        ).get("precision"),
        "workflow_component_selected_recall": workflow_metrics.get(
            "component_selected_truth_precision", {}
        ).get("recall"),
        "workflow_component_selected_f1": workflow_metrics.get(
            "component_selected_truth_precision", {}
        ).get("f1"),
        "workflow_component_selected_false_positive": workflow_metrics.get(
            "component_selected_truth_precision", {}
        ).get("false_positive"),
        "workflow_component_selected_false_negative": workflow_metrics.get(
            "component_selected_truth_precision", {}
        ).get("false_negative"),
        "workflow_truth_min_coverage": workflow_metrics.get("truth_recovery", {}).get(
            "min_best_coverage"
        ),
        "trinity_available": trinity.get("available"),
        "trinity_ran": trinity.get("ran"),
        "trinity_exit_code": trinity_result.get("exit_code"),
        "trinity_elapsed_seconds": trinity_result.get("elapsed_seconds"),
        "trinity_max_rss_kb": trinity_resources.get("max_rss_kb"),
        "trinity_user_seconds": trinity_resources.get("user_seconds"),
        "trinity_system_seconds": trinity_resources.get("system_seconds"),
        "trinity_gpu_available": trinity_gpu.get("available"),
        "trinity_gpu_max_memory_used_mib": first_gpu_metric(
            trinity_gpu, "max_memory_used_mib"
        ),
        "trinity_gpu_max_utilization_percent": first_gpu_metric(
            trinity_gpu, "max_utilization_gpu_percent"
        ),
        "trinity_gpu_max_power_draw_watts": first_gpu_metric(
            trinity_gpu, "max_power_draw_watts"
        ),
        "trinity_output_fasta_bytes": trinity_metrics.get("output_fasta_bytes"),
        "trinity_output_dir_bytes": trinity_metrics.get("output_dir_footprint", {}).get(
            "total_bytes"
        ),
        "trinity_output_file_count": trinity_metrics.get("output_dir_footprint", {}).get(
            "file_count"
        ),
        "trinity_lengths": trinity_metrics.get("lengths"),
        "trinity_truth_min_coverage": trinity_metrics.get("truth_recovery", {}).get(
            "min_best_coverage"
        ),
    }


def run_fixture(
    fixture: dict[str, object],
    defaults: dict[str, object],
    out_dir: Path,
    run_trinity: bool,
    require_trinity: bool,
) -> dict[str, object]:
    name = str(fixture["name"])
    inserts = [int(insert) for insert in fixture["inserts"]]
    fixture_out = out_dir / name
    run_raptor_workflow = bool(
        fixture.get("run_raptor_workflow", defaults.get("run_raptor_workflow", True))
    )
    run_raptor_workflow_gpu = bool(
        fixture.get("run_raptor_workflow_gpu", defaults.get("run_raptor_workflow_gpu", False))
    )
    run_raptor_workflow_single = bool(
        fixture.get(
            "run_raptor_workflow_single",
            defaults.get("run_raptor_workflow_single", False),
        )
    )
    run_malformed_fastq_checks = bool(
        fixture.get(
            "run_malformed_fastq_checks",
            defaults.get("run_malformed_fastq_checks", False),
        )
    )
    run_raptor_workflow_samples_file = bool(
        fixture.get(
            "run_raptor_workflow_samples_file",
            defaults.get("run_raptor_workflow_samples_file", False),
        )
    )
    skip_raptor = bool(fixture.get("skip_raptor", defaults.get("skip_raptor", True)))
    command = [
        sys.executable,
        str(FIXTURE_RUNNER),
        "--fixture",
        name,
        "--out-dir",
        str(fixture_out),
        "--insert-sweep",
        ",".join(str(insert) for insert in inserts),
        "--min-truth-coverage",
        str(defaults["min_truth_coverage"]),
        "--min-selected-match-coverage",
        str(defaults["min_selected_match_coverage"]),
        "--min-selected-precision",
        str(defaults["min_selected_precision"]),
        "--min-selected-f1",
        str(defaults["min_selected_f1"]),
    ]
    if run_raptor_workflow:
        command.append("--run-raptor-workflow")
    if run_raptor_workflow_gpu:
        command.append("--run-raptor-workflow-gpu")
    if run_raptor_workflow_single:
        command.append("--run-raptor-workflow-single")
    if run_raptor_workflow_samples_file:
        command.append("--run-raptor-workflow-samples-file")
    if run_malformed_fastq_checks:
        command.append("--run-malformed-fastq-checks")
    if skip_raptor:
        command.append("--skip-raptor")
    if run_trinity:
        command.append("--run-trinity")
    if require_trinity:
        command.append("--require-trinity")

    result = run_command(command, ROOT)
    report_path = fixture_report_path(fixture_out, inserts)
    summary = fixture_summary(report_path) if report_path.exists() else None
    output_footprint = directory_footprint(fixture_out)
    failures = [] if summary is None else list(summary.get("failures", []))
    if result["exit_code"] != 0 and not failures:
        failures.append("fixture runner failed before writing threshold failures")
    return {
        "fixture": name,
        "description": fixture.get("description"),
        "inserts": inserts,
        "command": command,
        "exit_code": result["exit_code"],
        "elapsed_seconds": result["elapsed_seconds"],
        "resource_usage": result["resource_usage"],
        "gpu_usage": result["gpu_usage"],
        "output_footprint": output_footprint,
        "report_path": str(report_path),
        "stdout": result["stdout"],
        "stderr": result["stderr"],
        "summary": summary,
        "failures": failures,
        "passed": result["exit_code"] == 0 and not failures,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, default=DEFAULT_PANEL)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--run-trinity", action="store_true")
    parser.add_argument("--require-trinity", action="store_true")
    args = parser.parse_args()

    panel = load_json(args.panel)
    defaults = panel.get("defaults", {})
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    fixture_results = [
        run_fixture(
            fixture,
            defaults,
            out_dir,
            args.run_trinity or args.require_trinity,
            args.require_trinity,
        )
        for fixture in panel.get("fixtures", [])
    ]
    failures = [
        f"{result['fixture']}: {failure}"
        for result in fixture_results
        for failure in result["failures"]
    ]
    report = {
        "panel": panel,
        "panel_path": str(args.panel),
        "out_dir": str(out_dir),
        "run_trinity": args.run_trinity or args.require_trinity,
        "require_trinity": args.require_trinity,
        "passed": not failures,
        "failures": failures,
        "fixtures": fixture_results,
    }
    report_path = out_dir / "panel_report.json"
    write_text(report_path, json.dumps(report, indent=2) + "\n")
    print(f"wrote {report_path}")
    if failures:
        for failure in failures:
            print(failure, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
