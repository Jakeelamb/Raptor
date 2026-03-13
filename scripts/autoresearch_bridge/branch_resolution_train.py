#!/usr/bin/env python3
"""Heldout-aware exhaustive search for branch-resolution component tasks."""

from __future__ import annotations

import argparse
import itertools
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path

from bridge_common import run_json_command


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BINARY = REPO_ROOT / "target" / "debug" / "raptor"
DEFAULT_TASKS_ROOT = (
    REPO_ROOT / "artifacts" / "autoresearch_raptor" / "branch_resolution" / "tasks"
)
DEFAULT_OUTPUT = (
    REPO_ROOT / "artifacts" / "autoresearch_raptor" / "branch_resolution" / "latest_result.json"
)

MIN_WIN_VALUES = (1, 2, 3, 4, 5, 6)
MIN_MARGIN_VALUES = (1, 2, 3, 4)
PREFER_NON_REPEAT_VALUES = (True, False)
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


def evaluate_config(binary: Path, task_root: Path, config: BranchConfig) -> dict:
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
    report = run_json_command(cmd)
    return report["summary"]


def metrics_key(metrics: dict) -> tuple[float, float]:
    return (metrics["score"], metrics["cases_per_second"])


def safe_scenario_rate(metrics: dict | None, scenario: str) -> float | None:
    if metrics is None:
        return None
    scenario_metrics = metrics.get("scenario_metrics") or {}
    scenario_report = scenario_metrics.get(scenario)
    if scenario_report is None:
        return None
    exact_rate = scenario_report.get("exact_rate")
    if exact_rate is None:
        return None
    return float(exact_rate)


def enumerate_configs() -> list[BranchConfig]:
    return [
        BranchConfig(
            branch_support_min_win=min_win,
            branch_support_min_margin=min_margin,
            prefer_non_repeat=prefer_non_repeat,
        )
        for min_win, min_margin, prefer_non_repeat in itertools.product(
            MIN_WIN_VALUES,
            MIN_MARGIN_VALUES,
            PREFER_NON_REPEAT_VALUES,
        )
    ]


def selection_key(row: dict) -> tuple[float, float, float, float, float, float]:
    val_metrics = row["val_metrics"]
    heldout_metrics = row.get("heldout_metrics")
    val_score = val_metrics["score"]
    heldout_score = heldout_metrics["score"] if heldout_metrics is not None else val_score
    robust_score = min(val_score, heldout_score)
    support_margin_rate = safe_scenario_rate(heldout_metrics, "support_margin_gate") or 0.0
    support_floor_rate = safe_scenario_rate(heldout_metrics, "support_floor_gate") or 0.0
    coverage_rate = safe_scenario_rate(heldout_metrics, "coverage_closeness") or 0.0
    non_repeat_rate = safe_scenario_rate(heldout_metrics, "non_repeat_preferred") or 0.0
    return (
        robust_score,
        support_margin_rate,
        support_floor_rate,
        coverage_rate,
        non_repeat_rate,
        val_metrics["cases_per_second"],
    )


def evaluate_grid(
    binary: Path,
    train_root: Path,
    val_root: Path,
    heldout_root: Path | None,
) -> tuple[list[dict], float]:
    start = time.time()
    evaluated: list[dict] = []

    for index, config in enumerate(enumerate_configs(), start=1):
        train_metrics = evaluate_config(binary, train_root, config)
        val_metrics = evaluate_config(binary, val_root, config)
        heldout_metrics = (
            evaluate_config(binary, heldout_root, config) if heldout_root is not None else None
        )
        evaluated.append(
            {
                "config": config,
                "train_metrics": train_metrics,
                "val_metrics": val_metrics,
                "heldout_metrics": heldout_metrics,
            }
        )
        if index % 8 == 0:
            leader = max(evaluated, key=selection_key)
            heldout_margin = safe_scenario_rate(
                leader.get("heldout_metrics"), "support_margin_gate"
            )
            print(
                f"candidates={index:03d} | "
                f"best_val_score={leader['val_metrics']['score']:.4f} | "
                f"best_heldout_margin={heldout_margin if heldout_margin is not None else float('nan'):.4f} | "
                f"best={asdict(leader['config'])}",
                flush=True,
            )

    return evaluated, time.time() - start


def heldout_non_regressing(candidate: dict, baseline: dict) -> bool:
    candidate_heldout = candidate.get("heldout_metrics")
    baseline_heldout = baseline.get("heldout_metrics")
    if candidate_heldout is None or baseline_heldout is None:
        return True
    if candidate_heldout["score"] < baseline_heldout["score"]:
        return False
    for scenario in HELDOUT_STRESS_SCENARIOS:
        candidate_rate = safe_scenario_rate(candidate_heldout, scenario)
        baseline_rate = safe_scenario_rate(baseline_heldout, scenario)
        if candidate_rate is None or baseline_rate is None:
            continue
        if candidate_rate < baseline_rate:
            return False
    return True


def choose_finalists(evaluated: list[dict]) -> tuple[dict, dict, dict, int]:
    baseline = next(row for row in evaluated if row["config"] == BranchConfig())
    best_train = max(
        evaluated,
        key=lambda row: metrics_key(row["train_metrics"]),
    )
    best_unconstrained = max(evaluated, key=selection_key)
    eligible = [row for row in evaluated if heldout_non_regressing(row, baseline)]
    if not eligible:
        eligible = [baseline]
    best_selected = max(eligible, key=selection_key)
    return best_train, best_selected, best_unconstrained, len(eligible)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the branch-resolution autoresearch bridge")
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--tasks-root", type=Path, default=DEFAULT_TASKS_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--skip-heldout", action="store_true")
    parser.add_argument("--require-heldout", action="store_true")
    args = parser.parse_args()

    if not args.binary.exists():
        raise FileNotFoundError(
            f"missing raptor binary at {args.binary}; build it first with cargo build"
        )

    train_root = args.tasks_root / "train"
    val_root = args.tasks_root / "val"
    heldout_root = args.tasks_root / "heldout"
    if not train_root.exists() or not val_root.exists():
        raise FileNotFoundError(
            f"expected train/val panels under {args.tasks_root}; run branch_resolution_prepare.py first"
        )
    has_heldout = heldout_root.exists()
    if args.require_heldout and not has_heldout:
        raise FileNotFoundError(
            f"expected heldout panel under {heldout_root}; pass --skip-heldout if you only want train/val"
        )

    print(f"Binary: {args.binary}")
    print(f"Train root: {train_root}")
    print(f"Val root: {val_root}")
    if has_heldout:
        print(f"Heldout root: {heldout_root}")
    else:
        print("Heldout root: missing (using train/val only)")
    print(f"Candidate count: {len(enumerate_configs())}")

    heldout_eval_root = None if args.skip_heldout or not has_heldout else heldout_root
    evaluated, search_seconds = evaluate_grid(args.binary, train_root, val_root, heldout_eval_root)
    best_train, best_selected, best_unconstrained, eligible_candidates = choose_finalists(evaluated)
    baseline = next(row for row in evaluated if row["config"] == BranchConfig())

    result = {
        "best_train_config": asdict(best_train["config"]),
        "best_train_metrics": best_train["train_metrics"],
        "best_val_config": asdict(best_selected["config"]),
        "val_metrics": best_selected["val_metrics"],
        "heldout_metrics": best_selected.get("heldout_metrics"),
        "baseline_config": asdict(baseline["config"]),
        "baseline_train_metrics": baseline["train_metrics"],
        "baseline_val_metrics": baseline["val_metrics"],
        "baseline_heldout_metrics": baseline.get("heldout_metrics"),
        "best_unconstrained_config": asdict(best_unconstrained["config"]),
        "best_unconstrained_val_metrics": best_unconstrained["val_metrics"],
        "best_unconstrained_heldout_metrics": best_unconstrained.get("heldout_metrics"),
        "eligible_candidates": eligible_candidates,
        "selection_policy": {
            "mode": "deterministic_grid_with_heldout_non_regression",
            "heldout_stress_scenarios": list(HELDOUT_STRESS_SCENARIOS),
        },
        "search_seconds": search_seconds,
        "evaluated_candidates": len(evaluated),
    }

    print(json.dumps(result, indent=2, sort_keys=True))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote branch-resolution result to {args.output}")


if __name__ == "__main__":
    main()
