#!/usr/bin/env python3
"""Promote a branch-resolution candidate through heldout eval and quick_test."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from bridge_common import run_json_command


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BINARY = REPO_ROOT / "target" / "release" / "raptor"
DEFAULT_TASKS_ROOT = (
    REPO_ROOT / "artifacts" / "autoresearch_raptor" / "branch_resolution" / "tasks"
)
DEFAULT_OUTPUT_ROOT = (
    REPO_ROOT / "artifacts" / "autoresearch_raptor" / "branch_resolution" / "promotions"
)
DEFAULT_TRAIN_RESULT = (
    REPO_ROOT
    / "artifacts"
    / "autoresearch_raptor"
    / "branch_resolution"
    / "latest_result.json"
)
QUICK_TEST_READS_1 = (
    REPO_ROOT / "bench" / "genome_assembly" / "data" / "quick_test" / "reads_1.fastq.gz"
)
QUICK_TEST_READS_2 = (
    REPO_ROOT / "bench" / "genome_assembly" / "data" / "quick_test" / "reads_2.fastq.gz"
)
QUICK_TEST_REFERENCE = (
    REPO_ROOT / "bench" / "genome_assembly" / "data" / "quick_test" / "reference.fa"
)
DEFAULT_BASELINE_CONFIG = {
    "branch_support_min_win": 2,
    "branch_support_min_margin": 1,
    "prefer_non_repeat": True,
}
HELDOUT_STRESS_SCENARIOS = (
    "support_floor_gate",
    "support_margin_gate",
    "coverage_closeness",
    "non_repeat_preferred",
)


@dataclass(frozen=True)
class BranchConfig:
    branch_support_min_win: int = 2
    branch_support_min_margin: int = 1
    prefer_non_repeat: bool = True


@dataclass(frozen=True)
class PromotionGate:
    min_val_score_delta: float = 0.0
    min_heldout_score_delta: float = 0.0
    max_heldout_scenario_drop: float = 0.0
    max_runtime_regression_ratio: float = 0.05
    max_runtime_regression_seconds: float = 5.0
    max_contig_n50_drop_ratio: float = 0.0
    max_scaffold_n50_drop_ratio: float = 0.0
    max_polish_correction_drop_ratio: float = 0.15
    max_reference_coverage_drop: float = 0.002
    max_aligned_query_fraction_drop: float = 0.01


def load_config(args: argparse.Namespace) -> BranchConfig:
    if args.train_result is not None and args.train_result.exists():
        result = json.loads(args.train_result.read_text(encoding="utf-8"))
        winner = result.get("best_val_config") or result.get("best_train_config")
        if winner:
            return BranchConfig(
                branch_support_min_win=winner["branch_support_min_win"],
                branch_support_min_margin=winner["branch_support_min_margin"],
                prefer_non_repeat=winner["prefer_non_repeat"],
            )

    return BranchConfig(
        branch_support_min_win=args.branch_support_min_win,
        branch_support_min_margin=args.branch_support_min_margin,
        prefer_non_repeat=not args.disable_prefer_non_repeat,
    )


def load_baseline_summary(summary_path: Path) -> dict[str, Any]:
    result = json.loads(summary_path.read_text(encoding="utf-8"))
    return {
        "config": result["config"],
        "train_metrics": result["train_metrics"],
        "val_metrics": result["val_metrics"],
        "heldout_metrics": result.get("heldout_metrics"),
        "quick_test": result.get("quick_test"),
        "summary_path": str(summary_path),
    }


def evaluate_branch_panel(binary: Path, task_root: Path, config: BranchConfig) -> dict[str, Any]:
    cmd = [
        str(binary),
        "component-bench",
        "branch-resolution",
        "--task",
        str(task_root),
        "--branch-support-min-win",
        str(config.branch_support_min_win),
        "--branch-support-min-margin",
        str(config.branch_support_min_margin),
        "--json",
    ]
    if not config.prefer_non_repeat:
        cmd.append("--disable-prefer-non-repeat")
    return run_json_command(cmd)["summary"]


def parse_first_int(pattern: str, text: str) -> int | None:
    match = re.search(pattern, text, re.MULTILINE)
    if match is None:
        return None
    return int(match.group(1))


def delta(candidate: float | int | None, baseline: float | int | None) -> float | None:
    if candidate is None or baseline is None:
        return None
    return float(candidate) - float(baseline)


def ratio_change(candidate: float | int | None, baseline: float | int | None) -> float | None:
    if candidate is None or baseline in (None, 0):
        return None
    return (float(candidate) - float(baseline)) / float(baseline)


def safe_scenario_rate(report: dict[str, Any] | None, scenario: str) -> float | None:
    if report is None:
        return None
    scenario_metrics = report.get("scenario_metrics") or {}
    scenario_report = scenario_metrics.get(scenario)
    if scenario_report is None:
        return None
    exact_rate = scenario_report.get("exact_rate")
    if exact_rate is None:
        return None
    return float(exact_rate)


def read_fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current_name: str | None = None
    current_len = 0
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    lengths[current_name] = current_len
                current_name = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line)
    if current_name is not None:
        lengths[current_name] = current_len
    return lengths


def merged_interval_bp(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    current_start, current_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
            continue
        total += current_end - current_start
        current_start, current_end = start, end
    total += current_end - current_start
    return total


def run_reference_eval(
    reference_path: Path,
    assembly_path: Path,
    output_dir: Path,
    minimap2_path: str | None,
    threads: int,
    min_mapq: int,
    min_aligned_bp: int,
) -> dict[str, Any]:
    if minimap2_path is None:
        return {"skipped": "minimap2 not found"}
    if not reference_path.exists():
        return {"skipped": f"reference not found at {reference_path}"}
    if not assembly_path.exists():
        return {"skipped": f"assembly not found at {assembly_path}"}

    output_dir.mkdir(parents=True, exist_ok=True)
    paf_path = output_dir / "alignments.paf"
    stderr_path = output_dir / "minimap2.stderr.log"
    cmd = [
        minimap2_path,
        "-x",
        "asm5",
        "-t",
        str(max(1, threads)),
        str(reference_path),
        str(assembly_path),
    ]
    completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
    paf_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")

    reference_lengths = read_fasta_lengths(reference_path)
    assembly_lengths = read_fasta_lengths(assembly_path)
    reference_intervals: dict[str, list[tuple[int, int]]] = {}
    query_intervals: dict[str, list[tuple[int, int]]] = {}
    query_records: dict[str, int] = {}
    matched_bases = 0
    aligned_bases = 0
    mapq_sum = 0
    kept_records = 0

    for line in completed.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        qname = fields[0]
        qstart = int(fields[2])
        qend = int(fields[3])
        tname = fields[5]
        tstart = int(fields[7])
        tend = int(fields[8])
        nmatch = int(fields[9])
        block_len = int(fields[10])
        mapq = int(fields[11])
        if mapq < min_mapq or (qend - qstart) < min_aligned_bp:
            continue
        kept_records += 1
        matched_bases += nmatch
        aligned_bases += block_len
        mapq_sum += mapq
        reference_intervals.setdefault(tname, []).append((tstart, tend))
        query_intervals.setdefault(qname, []).append((qstart, qend))
        query_records[qname] = query_records.get(qname, 0) + 1

    reference_total_bp = sum(reference_lengths.values())
    assembly_total_bp = sum(assembly_lengths.values())
    reference_covered_bp = sum(
        merged_interval_bp(intervals) for intervals in reference_intervals.values()
    )
    aligned_query_bp = sum(merged_interval_bp(intervals) for intervals in query_intervals.values())
    queries_with_alignment = len(query_intervals)
    queries_with_multiple_alignments = sum(
        1 for records in query_records.values() if records > 1
    )

    return {
        "tool": "minimap2",
        "reference_path": str(reference_path),
        "assembly_path": str(assembly_path),
        "paf_path": str(paf_path),
        "stderr_log": str(stderr_path),
        "reference_total_bp": reference_total_bp,
        "assembly_total_bp": assembly_total_bp,
        "reference_covered_bp": reference_covered_bp,
        "reference_coverage_fraction": (
            reference_covered_bp / reference_total_bp if reference_total_bp else None
        ),
        "aligned_query_bp": aligned_query_bp,
        "aligned_query_fraction": (
            aligned_query_bp / assembly_total_bp if assembly_total_bp else None
        ),
        "aligned_records": kept_records,
        "queries_with_alignment": queries_with_alignment,
        "queries_with_multiple_alignments": queries_with_multiple_alignments,
        "multi_alignment_query_fraction": (
            queries_with_multiple_alignments / queries_with_alignment
            if queries_with_alignment
            else None
        ),
        "mean_mapq": (mapq_sum / kept_records if kept_records else None),
        "matched_base_fraction": (matched_bases / aligned_bases if aligned_bases else None),
    }


def run_quick_test(
    binary: Path,
    config: BranchConfig,
    threads: int,
    output_dir: Path,
    reference_path: Path,
    minimap2_path: str | None,
    reference_eval_min_mapq: int,
    reference_eval_min_aligned_bp: int,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    contigs_path = output_dir / "contigs.fa"
    scaffold_path = output_dir / "contigs.scaffolds.fa"
    polished_path = output_dir / "contigs.polished.fa"
    cmd = [
        str(binary),
        "assemble-large",
        "-i",
        str(QUICK_TEST_READS_1),
        "--input2",
        str(QUICK_TEST_READS_2),
        "-o",
        str(contigs_path),
        "-t",
        str(threads),
        "--min-count",
        "0",
        "--scaffold",
        "--polish",
        "--polish-iterations",
        "1",
        "--compress-buckets",
        "--branch-support-min-win",
        str(config.branch_support_min_win),
        "--branch-support-min-margin",
        str(config.branch_support_min_margin),
    ]
    if not config.prefer_non_repeat:
        cmd.append("--disable-prefer-non-repeat")

    start = time.time()
    completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
    elapsed = time.time() - start

    stdout_path = output_dir / "stdout.log"
    stderr_path = output_dir / "stderr.log"
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")

    combined = completed.stdout + "\n" + completed.stderr
    result = {
        "elapsed_seconds": elapsed,
        "contigs": parse_first_int(r"^Contigs:\s+(\d+)$", combined),
        "total_length_bp": parse_first_int(r"^Total length:\s+(\d+)\s+bp$", combined),
        "contig_n50_bp": parse_first_int(r"^N50:\s+(\d+)\s+bp$", combined),
        "scaffolds": parse_first_int(r"^\s+Scaffolds:\s+(\d+)$", combined),
        "scaffold_n50_bp": parse_first_int(r"^\s+N50:\s+(\d+)\s+bp$", combined),
        "polish_corrections": parse_first_int(r"^\s+Corrections:\s+(\d+)$", combined),
        "stdout_log": str(stdout_path),
        "stderr_log": str(stderr_path),
        "contigs_path": str(contigs_path),
        "scaffold_path": str(scaffold_path),
        "polished_path": str(polished_path),
    }
    if polished_path.exists():
        reference_eval_input = polished_path
        reference_eval_target = "polished"
    elif scaffold_path.exists():
        reference_eval_input = scaffold_path
        reference_eval_target = "scaffold"
    else:
        reference_eval_input = contigs_path
        reference_eval_target = "contigs"
    result["reference_eval"] = run_reference_eval(
        reference_path,
        reference_eval_input,
        output_dir / "reference_eval",
        minimap2_path,
        threads,
        reference_eval_min_mapq,
        reference_eval_min_aligned_bp,
    )
    result["reference_eval_target"] = reference_eval_target
    return result


def compare_candidate_against_baseline(
    candidate: dict[str, Any], baseline: dict[str, Any]
) -> dict[str, Any]:
    candidate_quick_test = candidate.get("quick_test") or {}
    baseline_quick_test = baseline.get("quick_test") or {}
    candidate_reference_eval = candidate_quick_test.get("reference_eval") or {}
    baseline_reference_eval = baseline_quick_test.get("reference_eval") or {}
    candidate_heldout = candidate.get("heldout_metrics") or {}
    baseline_heldout = baseline.get("heldout_metrics") or {}

    comparison = {
        "train_score_delta": delta(
            candidate["train_metrics"]["score"], baseline["train_metrics"]["score"]
        ),
        "val_score_delta": delta(candidate["val_metrics"]["score"], baseline["val_metrics"]["score"]),
        "heldout_score_delta": delta(
            candidate_heldout.get("score"), baseline_heldout.get("score")
        ),
        "train_cases_per_second_delta": delta(
            candidate["train_metrics"]["cases_per_second"],
            baseline["train_metrics"]["cases_per_second"],
        ),
        "val_cases_per_second_delta": delta(
            candidate["val_metrics"]["cases_per_second"],
            baseline["val_metrics"]["cases_per_second"],
        ),
        "quick_test_elapsed_seconds_delta": delta(
            candidate_quick_test.get("elapsed_seconds"),
            baseline_quick_test.get("elapsed_seconds"),
        ),
        "quick_test_elapsed_ratio_change": ratio_change(
            candidate_quick_test.get("elapsed_seconds"),
            baseline_quick_test.get("elapsed_seconds"),
        ),
        "quick_test_contig_n50_delta": delta(
            candidate_quick_test.get("contig_n50_bp"),
            baseline_quick_test.get("contig_n50_bp"),
        ),
        "quick_test_contig_n50_ratio_change": ratio_change(
            candidate_quick_test.get("contig_n50_bp"),
            baseline_quick_test.get("contig_n50_bp"),
        ),
        "quick_test_scaffold_n50_delta": delta(
            candidate_quick_test.get("scaffold_n50_bp"),
            baseline_quick_test.get("scaffold_n50_bp"),
        ),
        "quick_test_scaffold_n50_ratio_change": ratio_change(
            candidate_quick_test.get("scaffold_n50_bp"),
            baseline_quick_test.get("scaffold_n50_bp"),
        ),
        "quick_test_polish_corrections_delta": delta(
            candidate_quick_test.get("polish_corrections"),
            baseline_quick_test.get("polish_corrections"),
        ),
        "quick_test_polish_corrections_ratio_change": ratio_change(
            candidate_quick_test.get("polish_corrections"),
            baseline_quick_test.get("polish_corrections"),
        ),
        "quick_test_total_length_delta": delta(
            candidate_quick_test.get("total_length_bp"),
            baseline_quick_test.get("total_length_bp"),
        ),
        "reference_coverage_fraction_delta": delta(
            candidate_reference_eval.get("reference_coverage_fraction"),
            baseline_reference_eval.get("reference_coverage_fraction"),
        ),
        "aligned_query_fraction_delta": delta(
            candidate_reference_eval.get("aligned_query_fraction"),
            baseline_reference_eval.get("aligned_query_fraction"),
        ),
        "multi_alignment_query_fraction_delta": delta(
            candidate_reference_eval.get("multi_alignment_query_fraction"),
            baseline_reference_eval.get("multi_alignment_query_fraction"),
        ),
        "matched_base_fraction_delta": delta(
            candidate_reference_eval.get("matched_base_fraction"),
            baseline_reference_eval.get("matched_base_fraction"),
        ),
    }
    for scenario in HELDOUT_STRESS_SCENARIOS:
        comparison[f"heldout_{scenario}_delta"] = delta(
            safe_scenario_rate(candidate_heldout, scenario),
            safe_scenario_rate(baseline_heldout, scenario),
        )
    return comparison


def decide_promotion(
    candidate: dict[str, Any],
    baseline: dict[str, Any],
    comparison: dict[str, Any],
    gate: PromotionGate,
) -> dict[str, Any]:
    hard_fail_reasons: list[str] = []
    review_reasons: list[str] = []

    if (
        comparison["val_score_delta"] is not None
        and comparison["val_score_delta"] < gate.min_val_score_delta
    ):
        hard_fail_reasons.append("validation score did not beat the baseline")

    heldout_score_delta = comparison["heldout_score_delta"]
    if heldout_score_delta is None:
        hard_fail_reasons.append("heldout panel evaluation missing for branch gate")
    elif heldout_score_delta < gate.min_heldout_score_delta:
        hard_fail_reasons.append("heldout score regressed versus baseline")

    for scenario in HELDOUT_STRESS_SCENARIOS:
        scenario_delta = comparison[f"heldout_{scenario}_delta"]
        if scenario_delta is None:
            hard_fail_reasons.append(f"heldout scenario metric missing for {scenario}")
            continue
        if scenario_delta < -gate.max_heldout_scenario_drop:
            hard_fail_reasons.append(
                f"heldout scenario {scenario} regressed beyond the allowed threshold"
            )

    contig_n50_ratio = comparison["quick_test_contig_n50_ratio_change"]
    if contig_n50_ratio is not None and contig_n50_ratio < -gate.max_contig_n50_drop_ratio:
        hard_fail_reasons.append("contig N50 regressed beyond the allowed threshold")

    scaffold_n50_ratio = comparison["quick_test_scaffold_n50_ratio_change"]
    if scaffold_n50_ratio is not None and scaffold_n50_ratio < -gate.max_scaffold_n50_drop_ratio:
        hard_fail_reasons.append("scaffold N50 regressed beyond the allowed threshold")

    reference_coverage_delta = comparison["reference_coverage_fraction_delta"]
    if (
        reference_coverage_delta is not None
        and reference_coverage_delta < -gate.max_reference_coverage_drop
    ):
        hard_fail_reasons.append("reference coverage regressed beyond the allowed threshold")

    aligned_query_fraction_delta = comparison["aligned_query_fraction_delta"]
    if (
        aligned_query_fraction_delta is not None
        and aligned_query_fraction_delta < -gate.max_aligned_query_fraction_drop
    ):
        hard_fail_reasons.append(
            "aligned assembly fraction regressed beyond the allowed threshold"
        )

    runtime_ratio = comparison["quick_test_elapsed_ratio_change"]
    runtime_delta = comparison["quick_test_elapsed_seconds_delta"]
    if runtime_ratio is not None and runtime_delta is not None:
        if (
            runtime_ratio > gate.max_runtime_regression_ratio
            and runtime_delta > gate.max_runtime_regression_seconds
        ):
            review_reasons.append("quick_test runtime regressed beyond the allowed threshold")

    correction_ratio = comparison["quick_test_polish_corrections_ratio_change"]
    matched_base_fraction_delta = comparison["matched_base_fraction_delta"]
    if (
        correction_ratio is not None
        and correction_ratio < -gate.max_polish_correction_drop_ratio
        and (matched_base_fraction_delta is None or matched_base_fraction_delta < 0.0)
    ):
        review_reasons.append("polish corrections dropped materially versus baseline")

    multi_alignment_delta = comparison["multi_alignment_query_fraction_delta"]
    if multi_alignment_delta is not None and multi_alignment_delta > 0.0:
        review_reasons.append("reference alignment fragmentation increased versus baseline")

    if hard_fail_reasons:
        status = "reject"
    elif review_reasons:
        status = "review"
    else:
        status = "promote_default"

    return {
        "status": status,
        "promote_default": status == "promote_default",
        "hard_fail_reasons": hard_fail_reasons,
        "review_reasons": review_reasons,
        "gate": asdict(gate),
        "candidate_summary": {
            "config": candidate["config"],
            "heldout_score": candidate.get("heldout_metrics", {}).get("score"),
            "heldout_scenario_rates": {
                scenario: safe_scenario_rate(candidate.get("heldout_metrics"), scenario)
                for scenario in HELDOUT_STRESS_SCENARIOS
            },
            "val_score": candidate["val_metrics"]["score"],
            "quick_test_elapsed_seconds": (candidate.get("quick_test") or {}).get(
                "elapsed_seconds"
            ),
            "quick_test_contig_n50_bp": (candidate.get("quick_test") or {}).get("contig_n50_bp"),
            "quick_test_scaffold_n50_bp": (candidate.get("quick_test") or {}).get(
                "scaffold_n50_bp"
            ),
        },
        "baseline_summary": {
            "config": baseline["config"],
            "heldout_score": baseline.get("heldout_metrics", {}).get("score"),
            "heldout_scenario_rates": {
                scenario: safe_scenario_rate(baseline.get("heldout_metrics"), scenario)
                for scenario in HELDOUT_STRESS_SCENARIOS
            },
            "val_score": baseline["val_metrics"]["score"],
            "quick_test_elapsed_seconds": (baseline.get("quick_test") or {}).get(
                "elapsed_seconds"
            ),
            "quick_test_contig_n50_bp": (baseline.get("quick_test") or {}).get("contig_n50_bp"),
            "quick_test_scaffold_n50_bp": (baseline.get("quick_test") or {}).get(
                "scaffold_n50_bp"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Promote a branch-resolution candidate through heldout eval and quick_test"
    )
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--tasks-root", type=Path, default=DEFAULT_TASKS_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--train-result", type=Path, default=DEFAULT_TRAIN_RESULT)
    parser.add_argument("--branch-support-min-win", type=int, default=2)
    parser.add_argument("--branch-support-min-margin", type=int, default=1)
    parser.add_argument("--disable-prefer-non-repeat", action="store_true")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--skip-quick-test", action="store_true")
    parser.add_argument("--baseline-summary", type=Path)
    parser.add_argument("--skip-baseline-quick-test", action="store_true")
    parser.add_argument("--min-heldout-score-delta", type=float, default=0.0)
    parser.add_argument("--max-heldout-scenario-drop", type=float, default=0.0)
    parser.add_argument("--max-runtime-regression-ratio", type=float, default=0.05)
    parser.add_argument("--max-runtime-regression-seconds", type=float, default=5.0)
    parser.add_argument("--max-contig-n50-drop-ratio", type=float, default=0.0)
    parser.add_argument("--max-scaffold-n50-drop-ratio", type=float, default=0.0)
    parser.add_argument("--max-polish-correction-drop-ratio", type=float, default=0.15)
    parser.add_argument("--max-reference-coverage-drop", type=float, default=0.002)
    parser.add_argument("--max-aligned-query-fraction-drop", type=float, default=0.01)
    parser.add_argument("--reference", type=Path, default=QUICK_TEST_REFERENCE)
    parser.add_argument("--minimap2", type=str, default=shutil.which("minimap2"))
    parser.add_argument("--reference-eval-min-mapq", type=int, default=20)
    parser.add_argument("--reference-eval-min-aligned-bp", type=int, default=500)
    args = parser.parse_args()

    if not args.binary.exists():
        raise FileNotFoundError(
            f"missing raptor binary at {args.binary}; build it first with cargo build --release"
        )

    train_root = args.tasks_root / "train"
    val_root = args.tasks_root / "val"
    heldout_root = args.tasks_root / "heldout"
    if not train_root.exists() or not val_root.exists() or not heldout_root.exists():
        raise FileNotFoundError(
            f"expected train/val/heldout panels under {args.tasks_root}; run branch_resolution_prepare.py first"
        )

    config = load_config(args)
    gate = PromotionGate(
        min_heldout_score_delta=args.min_heldout_score_delta,
        max_heldout_scenario_drop=args.max_heldout_scenario_drop,
        max_runtime_regression_ratio=args.max_runtime_regression_ratio,
        max_runtime_regression_seconds=args.max_runtime_regression_seconds,
        max_contig_n50_drop_ratio=args.max_contig_n50_drop_ratio,
        max_scaffold_n50_drop_ratio=args.max_scaffold_n50_drop_ratio,
        max_polish_correction_drop_ratio=args.max_polish_correction_drop_ratio,
        max_reference_coverage_drop=args.max_reference_coverage_drop,
        max_aligned_query_fraction_drop=args.max_aligned_query_fraction_drop,
    )
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = args.output_root / f"branch_resolution_{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Binary: {args.binary}")
    print(f"Tasks root: {args.tasks_root}")
    print(f"Promotion output: {output_dir}")
    print(f"Config: {asdict(config)}")

    candidate_result = {
        "config": asdict(config),
        "train_metrics": evaluate_branch_panel(args.binary, train_root, config),
        "val_metrics": evaluate_branch_panel(args.binary, val_root, config),
        "heldout_metrics": evaluate_branch_panel(args.binary, heldout_root, config),
    }

    if args.baseline_summary is not None and args.baseline_summary.exists():
        baseline_result = load_baseline_summary(args.baseline_summary)
    else:
        baseline_config = BranchConfig(**DEFAULT_BASELINE_CONFIG)
        baseline_result = {
            "config": asdict(baseline_config),
            "train_metrics": evaluate_branch_panel(args.binary, train_root, baseline_config),
            "val_metrics": evaluate_branch_panel(args.binary, val_root, baseline_config),
            "heldout_metrics": evaluate_branch_panel(args.binary, heldout_root, baseline_config),
            "quick_test": None,
        }

    if baseline_result.get("heldout_metrics") is None:
        baseline_config = BranchConfig(**baseline_result["config"])
        baseline_result["heldout_metrics"] = evaluate_branch_panel(
            args.binary, heldout_root, baseline_config
        )
    else:
        baseline_config = BranchConfig(**baseline_result["config"])

    quick_data_available = QUICK_TEST_READS_1.exists() and QUICK_TEST_READS_2.exists()

    if args.skip_quick_test:
        candidate_result["quick_test"] = None
        baseline_result["quick_test"] = None
    elif quick_data_available:
        if baseline_result.get("quick_test") is None:
            baseline_result["quick_test"] = {"skipped": "pending heldout stress gate"}
        candidate_result["quick_test"] = {"skipped": "pending heldout stress gate"}
    else:
        candidate_result["quick_test"] = {"skipped": "quick_test data not found"}
        if baseline_result.get("quick_test") is None:
            baseline_result["quick_test"] = {"skipped": "quick_test data not found"}

    comparison = compare_candidate_against_baseline(candidate_result, baseline_result)
    pre_decision = decide_promotion(candidate_result, baseline_result, comparison, gate)
    heldout_gate_failed = any(
        reason.startswith("heldout") for reason in pre_decision["hard_fail_reasons"]
    )

    if args.skip_quick_test:
        pass
    elif pre_decision["status"] in {"promote_default", "review"}:
        candidate_result["quick_test"] = run_quick_test(
            args.binary,
            config,
            args.threads,
            output_dir / "quick_test",
            args.reference,
            args.minimap2,
            args.reference_eval_min_mapq,
            args.reference_eval_min_aligned_bp,
        )
        if not args.skip_baseline_quick_test and baseline_result.get("quick_test", {}).get(
            "skipped"
        ) == "pending heldout stress gate":
            baseline_result["quick_test"] = run_quick_test(
                args.binary,
                baseline_config,
                args.threads,
                output_dir / "baseline_quick_test",
                args.reference,
                args.minimap2,
                args.reference_eval_min_mapq,
                args.reference_eval_min_aligned_bp,
            )
        elif args.skip_baseline_quick_test:
            baseline_result["quick_test"] = None
    elif pre_decision["status"] == "reject":
        if heldout_gate_failed:
            candidate_result["quick_test"] = {"skipped": "heldout stress gate failed"}
            if baseline_result.get("quick_test", {}).get("skipped") == "pending heldout stress gate":
                baseline_result["quick_test"] = {"skipped": "heldout stress gate failed"}
        elif candidate_result.get("quick_test", {}).get("skipped") == "pending heldout stress gate":
            candidate_result["quick_test"] = {"skipped": "preliminary promotion rejection"}
            if baseline_result.get("quick_test", {}).get("skipped") == "pending heldout stress gate":
                baseline_result["quick_test"] = {"skipped": "preliminary promotion rejection"}

    comparison = compare_candidate_against_baseline(candidate_result, baseline_result)
    promotion_decision = decide_promotion(candidate_result, baseline_result, comparison, gate)

    summary = {
        "config": candidate_result["config"],
        "train_metrics": candidate_result["train_metrics"],
        "val_metrics": candidate_result["val_metrics"],
        "heldout_metrics": candidate_result["heldout_metrics"],
        "baseline": baseline_result,
        "quick_test": candidate_result.get("quick_test"),
        "comparison": comparison,
        "promotion_decision": promotion_decision,
    }

    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Wrote promotion summary to {summary_path}")


if __name__ == "__main__":
    main()
