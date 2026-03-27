#!/usr/bin/env python3
"""Run a baseline-plus-candidates patch matrix using existing promotion tooling."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from bridge_common import format_command
from profile_assembler_quick_test import run_profile


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "artifacts" / "autoresearch_raptor" / "patch_matrix"
DEFAULT_THREADS = 8
PROMOTE_SCRIPTS = {
    "read-mapping": REPO_ROOT / "scripts" / "autoresearch_bridge" / "promote_read_mapping_candidate.py",
    "branch-resolution": REPO_ROOT
    / "scripts"
    / "autoresearch_bridge"
    / "promote_branch_resolution_candidate.py",
    "scaffold-polish": REPO_ROOT
    / "scripts"
    / "autoresearch_bridge"
    / "promote_scaffold_polish_candidate.py",
    "contig-extraction": REPO_ROOT
    / "scripts"
    / "autoresearch_bridge"
    / "promote_contig_extraction_candidate.py",
}


@dataclass(frozen=True)
class LabeledBinary:
    label: str
    path: Path


def parse_labeled_binary(raw: str) -> LabeledBinary:
    if "=" not in raw:
        raise argparse.ArgumentTypeError("expected LABEL=/path/to/binary")
    label, raw_path = raw.split("=", 1)
    label = label.strip()
    path = Path(raw_path).expanduser().resolve()
    if not label:
        raise argparse.ArgumentTypeError("candidate label must be non-empty")
    return LabeledBinary(label=label, path=path)


def delta(candidate: float | int | None, baseline: float | int | None) -> float | None:
    if candidate is None or baseline is None:
        return None
    return float(candidate) - float(baseline)


def ratio_change(candidate: float | int | None, baseline: float | int | None) -> float | None:
    if candidate is None or baseline in (None, 0):
        return None
    return (float(candidate) - float(baseline)) / float(baseline)


def compare_quick_profiles(candidate: dict[str, Any], baseline: dict[str, Any]) -> dict[str, Any]:
    comparison = {
        "elapsed_seconds_delta": delta(candidate.get("elapsed_seconds"), baseline.get("elapsed_seconds")),
        "elapsed_ratio_change": ratio_change(
            candidate.get("elapsed_seconds"), baseline.get("elapsed_seconds")
        ),
        "contigs_delta": delta(candidate.get("contigs"), baseline.get("contigs")),
        "contig_n50_bp_delta": delta(candidate.get("contig_n50_bp"), baseline.get("contig_n50_bp")),
        "scaffolds_delta": delta(candidate.get("scaffolds"), baseline.get("scaffolds")),
        "scaffold_n50_bp_delta": delta(
            candidate.get("scaffold_n50_bp"), baseline.get("scaffold_n50_bp")
        ),
        "polish_corrections_delta": delta(
            candidate.get("polish_corrections"), baseline.get("polish_corrections")
        ),
        "core_phase_timing_deltas": {},
        "shared_phase_timing_deltas": {},
    }
    for key, baseline_value in (baseline.get("core_phase_timings") or {}).items():
        comparison["core_phase_timing_deltas"][key] = delta(
            (candidate.get("core_phase_timings") or {}).get(key), baseline_value
        )
    for key, baseline_value in (baseline.get("shared_phase_timings") or {}).items():
        comparison["shared_phase_timing_deltas"][key] = delta(
            (candidate.get("shared_phase_timings") or {}).get(key), baseline_value
        )
    return comparison


def run_promotion(
    component: str,
    binary: Path,
    output_root: Path,
    threads: int,
    promote_args: list[str],
    baseline_summary: Path | None,
) -> dict[str, Any]:
    script = PROMOTE_SCRIPTS[component]
    output_root.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(script),
        "--binary",
        str(binary),
        "--output-root",
        str(output_root),
        "--threads",
        str(threads),
    ]
    if baseline_summary is not None:
        cmd.extend(["--baseline-summary", str(baseline_summary)])
    cmd.extend(promote_args)
    completed = subprocess.run(cmd, check=True, capture_output=True, text=True, cwd=REPO_ROOT)
    stdout_path = output_root / "promotion.stdout.log"
    stderr_path = output_root / "promotion.stderr.log"
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")
    summary_path = output_root / "summary.json"
    result = json.loads(summary_path.read_text(encoding="utf-8"))
    result["command"] = cmd
    result["command_pretty"] = format_command(cmd)
    result["stdout_log"] = str(stdout_path)
    result["stderr_log"] = str(stderr_path)
    result["summary_path"] = str(summary_path)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=parse_labeled_binary)
    parser.add_argument("--candidate", action="append", required=True, type=parse_labeled_binary)
    parser.add_argument(
        "--component",
        choices=sorted(PROMOTE_SCRIPTS),
        help="Optional component promotion script to run alongside quick-test profiling.",
    )
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--quick-extra-arg",
        action="append",
        default=[],
        help="Extra assemble-large arg for the quick-test profiler. Repeat for multiple args.",
    )
    parser.add_argument(
        "--promote-arg",
        action="append",
        default=[],
        help="Extra arg forwarded verbatim to the component promotion script. Repeat as needed.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)

    baseline_dir = args.output_root / args.baseline.label
    baseline_quick = run_profile(
        args.baseline.path,
        args.threads,
        args.baseline.label,
        args.output_root / "quick_profiles",
        args.quick_extra_arg,
    )

    baseline_promotion = None
    baseline_summary_path = None
    if args.component is not None:
        baseline_promotion = run_promotion(
            args.component,
            args.baseline.path,
            baseline_dir / "promotion",
            args.threads,
            args.promote_arg,
            baseline_summary=None,
        )
        baseline_summary_path = Path(baseline_promotion["summary_path"])

    candidates: list[dict[str, Any]] = []
    for candidate in args.candidate:
        candidate_dir = args.output_root / candidate.label
        quick = run_profile(
            candidate.path,
            args.threads,
            candidate.label,
            args.output_root / "quick_profiles",
            args.quick_extra_arg,
        )
        promotion = None
        if args.component is not None and baseline_summary_path is not None:
            promotion = run_promotion(
                args.component,
                candidate.path,
                candidate_dir / "promotion",
                args.threads,
                args.promote_arg,
                baseline_summary=baseline_summary_path,
            )
        candidates.append(
            {
                "label": candidate.label,
                "binary": str(candidate.path),
                "quick_profile": quick,
                "quick_vs_baseline": compare_quick_profiles(quick, baseline_quick),
                "promotion": promotion,
            }
        )

    result = {
        "baseline": {
            "label": args.baseline.label,
            "binary": str(args.baseline.path),
            "quick_profile": baseline_quick,
            "promotion": baseline_promotion,
        },
        "component": args.component,
        "threads": args.threads,
        "quick_extra_args": args.quick_extra_arg,
        "promote_args": args.promote_arg,
        "candidates": candidates,
    }

    summary_path = args.output_root / "matrix_summary.json"
    summary_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    print(f"Wrote patch matrix summary to {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
