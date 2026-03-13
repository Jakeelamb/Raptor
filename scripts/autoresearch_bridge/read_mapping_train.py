#!/usr/bin/env python3
"""Autoresearch-style search loop for Raptor read-mapping component tasks."""

from __future__ import annotations

import argparse
import itertools
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable

from bridge_common import run_json_command


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BINARY = REPO_ROOT / "target" / "debug" / "raptor"
DEFAULT_TASKS_ROOT = REPO_ROOT / "artifacts" / "autoresearch_raptor" / "read_mapping" / "tasks"
DEFAULT_OUTPUT = (
    REPO_ROOT / "artifacts" / "autoresearch_raptor" / "read_mapping" / "latest_result.json"
)

K_VALUES = (9, 11, 13, 15, 17, 19, 21)
W_VALUES = (4, 5, 6, 8, 10, 12)
MIN_PRIMARY_MATCH_VALUES = (2, 3, 4, 5)
MIN_SCAFFOLD_MATCH_VALUES = (1, 2, 3, 4)
SMOKE_K_VALUES = (9, 11, 15)
SMOKE_W_VALUES = (4, 5, 10)
SMOKE_MIN_PRIMARY_MATCH_VALUES = (3, 4, 5)
SMOKE_MIN_SCAFFOLD_MATCH_VALUES = (2, 3, 4)
SEARCH_POLICY_CHOICES = ("grid", "successive-halving")
DEFAULT_HALVING_FACTOR = 4
DEFAULT_FINAL_HELDOUT_CANDIDATES = 32
DEFAULT_SMOKE_FINAL_HELDOUT_CANDIDATES = 12


@dataclass(frozen=True)
class MapperConfig:
    k: int = 15
    w: int = 10
    min_primary_matches: int = 3
    min_scaffold_matches: int = 2


def evaluate_config(
    binary: Path,
    task_root: Path,
    config: MapperConfig,
    position_tolerance: int,
) -> dict[str, Any]:
    cmd = [
        str(binary),
        "component-bench",
        "read-mapping",
        "--task",
        str(task_root),
        "--k",
        str(config.k),
        "--w",
        str(config.w),
        "--min-primary-matches",
        str(config.min_primary_matches),
        "--min-scaffold-matches",
        str(config.min_scaffold_matches),
        "--position-tolerance",
        str(position_tolerance),
        "--json",
    ]
    report = run_json_command(cmd)
    return report["summary"]


def ordered_grid(smoke: bool) -> list[MapperConfig]:
    if smoke:
        k_values = SMOKE_K_VALUES
        w_values = SMOKE_W_VALUES
        min_primary_values = SMOKE_MIN_PRIMARY_MATCH_VALUES
        min_scaffold_values = SMOKE_MIN_SCAFFOLD_MATCH_VALUES
    else:
        k_values = K_VALUES
        w_values = W_VALUES
        min_primary_values = MIN_PRIMARY_MATCH_VALUES
        min_scaffold_values = MIN_SCAFFOLD_MATCH_VALUES

    configs = [
        MapperConfig(
            k=k,
            w=w,
            min_primary_matches=min_primary_matches,
            min_scaffold_matches=min_scaffold_matches,
        )
        for k, w, min_primary_matches, min_scaffold_matches in itertools.product(
            k_values, w_values, min_primary_values, min_scaffold_values
        )
    ]

    default = MapperConfig()
    if default in configs:
        configs.remove(default)
    return [default, *configs]


def safe_metric(report: dict[str, Any] | None, key: str) -> float | None:
    if report is None:
        return None
    value = report.get(key)
    if value is not None:
        return float(value)
    if key == "scaffold_specificity_rate":
        value = report.get("primary_specificity_rate")
        if value is not None:
            return float(value)
    return None


def empty_row(config: MapperConfig) -> dict[str, Any]:
    return {
        "config": config,
        "train_metrics": None,
        "val_metrics": None,
        "heldout_metrics": None,
    }


def train_stage_key(row: dict[str, Any]) -> tuple[float, float]:
    metrics = row["train_metrics"]
    if metrics is None:
        return (float("-inf"), float("-inf"))
    return (metrics["score"], metrics["reads_per_second"])


def val_stage_key(row: dict[str, Any]) -> tuple[float, float, float]:
    metrics = row["val_metrics"]
    if metrics is None:
        return (float("-inf"), float("-inf"), float("-inf"))
    return (
        metrics["score"],
        safe_metric(metrics, "scaffold_specificity_rate") or 0.0,
        metrics["reads_per_second"],
    )


def selection_key(row: dict[str, Any]) -> tuple[float, float, float, float, float]:
    val_metrics = row.get("val_metrics")
    if val_metrics is None:
        return (float("-inf"), float("-inf"), float("-inf"), float("-inf"), float("-inf"))

    val_score = val_metrics["score"]
    val_reads_per_second = val_metrics["reads_per_second"]
    heldout_metrics = row.get("heldout_metrics")
    heldout_score = safe_metric(heldout_metrics, "score")
    heldout_specificity = safe_metric(heldout_metrics, "scaffold_specificity_rate")

    if heldout_score is None:
        return (
            val_score,
            safe_metric(val_metrics, "scaffold_specificity_rate") or 0.0,
            val_score,
            0.0,
            val_reads_per_second,
        )

    robust_score = min(val_score, heldout_score)
    return (
        robust_score,
        heldout_specificity or 0.0,
        heldout_score,
        val_score,
        val_reads_per_second,
    )


def stage_key(split_name: str) -> Callable[[dict[str, Any]], tuple]:
    if split_name == "train":
        return train_stage_key
    if split_name == "val":
        return val_stage_key
    return selection_key


def keep_top_rows(
    rows: list[dict[str, Any]],
    *,
    keep_count: int,
    key_fn: Callable[[dict[str, Any]], tuple],
    baseline_config: MapperConfig,
) -> list[dict[str, Any]]:
    if not rows:
        return []
    keep_count = max(1, min(len(rows), keep_count))
    ordered = sorted(rows, key=key_fn, reverse=True)
    survivors = ordered[:keep_count]
    baseline_row = next((row for row in rows if row["config"] == baseline_config), None)
    if baseline_row is not None and all(row["config"] != baseline_config for row in survivors):
        survivors[-1] = baseline_row
    deduped: list[dict[str, Any]] = []
    seen: set[MapperConfig] = set()
    for row in survivors:
        if row["config"] in seen:
            continue
        deduped.append(row)
        seen.add(row["config"])
    return deduped


def print_stage_progress(
    *,
    split_name: str,
    completed: int,
    total: int,
    candidates: list[dict[str, Any]],
    key_fn: Callable[[dict[str, Any]], tuple],
) -> None:
    ready = [row for row in candidates if row.get(f"{split_name}_metrics") is not None]
    if not ready:
        return
    leader = max(ready, key=key_fn)
    val_metrics = leader.get("val_metrics")
    heldout_specificity = safe_metric(leader.get("heldout_metrics"), "scaffold_specificity_rate")
    if split_name == "train":
        leader_score = leader["train_metrics"]["score"]
        detail = f"best_train_score={leader_score:.4f}"
    elif split_name == "val":
        leader_score = val_metrics["score"] if val_metrics is not None else float("nan")
        detail = f"best_val_score={leader_score:.4f}"
    else:
        leader_score = val_metrics["score"] if val_metrics is not None else float("nan")
        detail = (
            f"best_val_score={leader_score:.4f} | "
            f"best_heldout_spec={heldout_specificity if heldout_specificity is not None else float('nan'):.4f}"
        )
    print(
        f"stage={split_name} | candidates={completed:03d}/{total:03d} | "
        f"{detail} | best={asdict(leader['config'])}",
        flush=True,
    )


def evaluate_split_rows(
    *,
    binary: Path,
    task_root: Path,
    split_name: str,
    rows: list[dict[str, Any]],
    start_time: float,
    budget_seconds: float,
    position_tolerance: int,
) -> tuple[bool, int]:
    key_fn = stage_key(split_name)
    metrics_key = f"{split_name}_metrics"
    evaluated_now = 0
    total = len(rows)

    for index, row in enumerate(rows, start=1):
        if budget_seconds > 0.0 and (time.time() - start_time) >= budget_seconds:
            return False, evaluated_now
        if row[metrics_key] is not None:
            continue
        row[metrics_key] = evaluate_config(binary, task_root, row["config"], position_tolerance)
        evaluated_now += 1
        if index % 12 == 0 or index == total:
            print_stage_progress(
                split_name=split_name,
                completed=index,
                total=total,
                candidates=rows,
                key_fn=key_fn,
            )

    return True, evaluated_now


def describe_stage(
    *,
    split_name: str,
    evaluated_rows: list[dict[str, Any]],
    survivors: list[dict[str, Any]],
) -> dict[str, Any]:
    key_fn = stage_key(split_name)
    leader = max(evaluated_rows, key=key_fn) if evaluated_rows else None
    return {
        "stage": split_name,
        "evaluated_candidates": len(evaluated_rows),
        "survivor_candidates": len(survivors),
        "best_config": asdict(leader["config"]) if leader is not None else None,
        "best_metrics": leader.get(f"{split_name}_metrics") if leader is not None else None,
    }


def evaluate_grid(
    *,
    binary: Path,
    train_root: Path,
    val_root: Path,
    heldout_root: Path | None,
    configs: list[MapperConfig],
    budget_seconds: float,
    position_tolerance: int,
) -> tuple[list[dict[str, Any]], float, bool, dict[str, int], list[dict[str, Any]]]:
    start = time.time()
    rows = [empty_row(config) for config in configs]
    completed, train_evals = evaluate_split_rows(
        binary=binary,
        task_root=train_root,
        split_name="train",
        rows=rows,
        start_time=start,
        budget_seconds=budget_seconds,
        position_tolerance=position_tolerance,
    )
    if not completed:
        return rows, time.time() - start, False, {"train": train_evals, "val": 0, "heldout": 0}, []

    completed, val_evals = evaluate_split_rows(
        binary=binary,
        task_root=val_root,
        split_name="val",
        rows=rows,
        start_time=start,
        budget_seconds=budget_seconds,
        position_tolerance=position_tolerance,
    )
    stage_records = [
        describe_stage(split_name="train", evaluated_rows=rows, survivors=rows),
        describe_stage(split_name="val", evaluated_rows=rows, survivors=rows),
    ]
    if not completed:
        return rows, time.time() - start, False, {"train": train_evals, "val": val_evals, "heldout": 0}, stage_records

    heldout_evals = 0
    if heldout_root is not None:
        completed, heldout_evals = evaluate_split_rows(
            binary=binary,
            task_root=heldout_root,
            split_name="heldout",
            rows=rows,
            start_time=start,
            budget_seconds=budget_seconds,
            position_tolerance=position_tolerance,
        )
        stage_records.append(
            describe_stage(split_name="heldout", evaluated_rows=rows, survivors=rows)
        )
    return (
        rows,
        time.time() - start,
        completed,
        {"train": train_evals, "val": val_evals, "heldout": heldout_evals},
        stage_records,
    )


def evaluate_successive_halving(
    *,
    binary: Path,
    train_root: Path,
    val_root: Path,
    heldout_root: Path | None,
    configs: list[MapperConfig],
    budget_seconds: float,
    position_tolerance: int,
    halving_factor: int,
    final_heldout_candidates: int,
) -> tuple[list[dict[str, Any]], float, bool, dict[str, int], list[dict[str, Any]]]:
    start = time.time()
    rows = [empty_row(config) for config in configs]
    baseline_config = MapperConfig()

    completed, train_evals = evaluate_split_rows(
        binary=binary,
        task_root=train_root,
        split_name="train",
        rows=rows,
        start_time=start,
        budget_seconds=budget_seconds,
        position_tolerance=position_tolerance,
    )
    if not completed:
        return rows, time.time() - start, False, {"train": train_evals, "val": 0, "heldout": 0}, []

    train_survivor_count = max(
        final_heldout_candidates,
        (len(rows) + max(2, halving_factor) - 1) // max(2, halving_factor),
    )
    train_survivors = keep_top_rows(
        rows,
        keep_count=train_survivor_count,
        key_fn=train_stage_key,
        baseline_config=baseline_config,
    )

    completed, val_evals = evaluate_split_rows(
        binary=binary,
        task_root=val_root,
        split_name="val",
        rows=train_survivors,
        start_time=start,
        budget_seconds=budget_seconds,
        position_tolerance=position_tolerance,
    )
    stage_records = [
        describe_stage(split_name="train", evaluated_rows=rows, survivors=train_survivors),
    ]
    if not completed:
        return (
            rows,
            time.time() - start,
            False,
            {"train": train_evals, "val": val_evals, "heldout": 0},
            stage_records,
        )

    if heldout_root is None:
        stage_records.append(
            describe_stage(split_name="val", evaluated_rows=train_survivors, survivors=train_survivors)
        )
        return (
            rows,
            time.time() - start,
            True,
            {"train": train_evals, "val": val_evals, "heldout": 0},
            stage_records,
        )

    heldout_survivors = keep_top_rows(
        train_survivors,
        keep_count=final_heldout_candidates,
        key_fn=val_stage_key,
        baseline_config=baseline_config,
    )
    completed, heldout_evals = evaluate_split_rows(
        binary=binary,
        task_root=heldout_root,
        split_name="heldout",
        rows=heldout_survivors,
        start_time=start,
        budget_seconds=budget_seconds,
        position_tolerance=position_tolerance,
    )
    stage_records.append(
        describe_stage(split_name="val", evaluated_rows=train_survivors, survivors=heldout_survivors)
    )
    stage_records.append(
        describe_stage(
            split_name="heldout",
            evaluated_rows=heldout_survivors,
            survivors=heldout_survivors,
        )
    )
    return (
        rows,
        time.time() - start,
        completed,
        {"train": train_evals, "val": val_evals, "heldout": heldout_evals},
        stage_records,
    )


def choose_finalist(
    evaluated: list[dict[str, Any]],
    baseline: dict[str, Any],
    max_heldout_scaffold_specificity_drop: float,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], int]:
    best_train = max(
        (row for row in evaluated if row["train_metrics"] is not None),
        key=train_stage_key,
    )
    val_ready = [row for row in evaluated if row["val_metrics"] is not None]
    if not val_ready:
        raise RuntimeError("read-mapping search did not evaluate any validation candidates")

    constrained = val_ready
    baseline_heldout_specificity = safe_metric(
        baseline.get("heldout_metrics"), "scaffold_specificity_rate"
    )
    if baseline_heldout_specificity is not None:
        constrained = []
        for row in val_ready:
            heldout_specificity = safe_metric(row.get("heldout_metrics"), "scaffold_specificity_rate")
            if heldout_specificity is None:
                continue
            if heldout_specificity + max_heldout_scaffold_specificity_drop < baseline_heldout_specificity:
                continue
            constrained.append(row)
        if not constrained:
            constrained = [baseline]

    best_selected = max(constrained, key=selection_key)
    best_unconstrained = max(val_ready, key=selection_key)
    return best_train, best_selected, best_unconstrained, len(constrained)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the read-mapping autoresearch bridge")
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--tasks-root", type=Path, default=DEFAULT_TASKS_ROOT)
    parser.add_argument("--time-budget", type=float, default=300.0)
    parser.add_argument("--position-tolerance", type=int, default=8)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--smoke", action="store_true", help="Short developer run")
    parser.add_argument("--skip-heldout", action="store_true", help="Do not evaluate heldout split")
    parser.add_argument(
        "--search-policy",
        choices=SEARCH_POLICY_CHOICES,
        default="grid",
        help="Use full deterministic grid or multi-fidelity successive halving",
    )
    parser.add_argument(
        "--halving-factor",
        type=int,
        default=DEFAULT_HALVING_FACTOR,
        help="Successive-halving factor for narrowing the train winner set",
    )
    parser.add_argument(
        "--final-heldout-candidates",
        type=int,
        default=DEFAULT_FINAL_HELDOUT_CANDIDATES,
        help="Maximum number of configs that reach heldout under successive halving",
    )
    parser.add_argument(
        "--max-heldout-scaffold-specificity-drop",
        type=float,
        default=0.0,
        help="Reject finalist configs that regress heldout scaffold specificity more than this amount",
    )
    parser.add_argument(
        "--require-heldout",
        action="store_true",
        help="Fail if heldout panel is missing",
    )
    args = parser.parse_args()

    if not args.binary.exists():
        raise FileNotFoundError(
            f"missing raptor binary at {args.binary}; build it first with cargo build"
        )
    if args.halving_factor < 2:
        raise ValueError("--halving-factor must be at least 2")
    if args.final_heldout_candidates < 1:
        raise ValueError("--final-heldout-candidates must be at least 1")

    train_root = args.tasks_root / "train"
    val_root = args.tasks_root / "val"
    heldout_root = args.tasks_root / "heldout"
    if not train_root.exists() or not val_root.exists():
        raise FileNotFoundError(
            f"expected train/val panels under {args.tasks_root}; run read_mapping_prepare.py first"
        )
    has_heldout = heldout_root.exists()
    if args.require_heldout and not has_heldout:
        raise FileNotFoundError(
            f"expected heldout panel under {heldout_root}; pass --skip-heldout if you only want train/val"
        )

    budget_seconds = 45.0 if args.smoke else args.time_budget
    configs = ordered_grid(args.smoke)
    heldout_eval_root = None if args.skip_heldout or not has_heldout else heldout_root
    final_heldout_candidates = (
        min(DEFAULT_SMOKE_FINAL_HELDOUT_CANDIDATES, len(configs))
        if args.smoke
        else min(args.final_heldout_candidates, len(configs))
    )

    print(f"Binary: {args.binary}")
    print(f"Train root: {train_root}")
    print(f"Val root: {val_root}")
    if has_heldout:
        print(f"Heldout root: {heldout_root}")
    else:
        print("Heldout root: missing (using train/val only)")
    print(f"Budget seconds: {budget_seconds:.1f}")
    print(f"Position tolerance: {args.position_tolerance}")
    print(f"Search policy: {args.search_policy}")
    print(f"Grid candidates: {len(configs)}")
    if args.search_policy == "successive-halving":
        print(f"Halving factor: {args.halving_factor}")
        print(f"Final heldout candidates: {final_heldout_candidates}")
    print(
        "Max heldout scaffold specificity drop: "
        f"{args.max_heldout_scaffold_specificity_drop:.4f}"
    )

    if args.search_policy == "grid":
        evaluated, search_seconds, completed_search, split_evaluations, stage_records = evaluate_grid(
            binary=args.binary,
            train_root=train_root,
            val_root=val_root,
            heldout_root=heldout_eval_root,
            configs=configs,
            budget_seconds=budget_seconds,
            position_tolerance=args.position_tolerance,
        )
        completed_grid = completed_search
    else:
        evaluated, search_seconds, completed_search, split_evaluations, stage_records = (
            evaluate_successive_halving(
                binary=args.binary,
                train_root=train_root,
                val_root=val_root,
                heldout_root=heldout_eval_root,
                configs=configs,
                budget_seconds=budget_seconds,
                position_tolerance=args.position_tolerance,
                halving_factor=args.halving_factor,
                final_heldout_candidates=final_heldout_candidates,
            )
        )
        completed_grid = False

    if not any(row["train_metrics"] is not None for row in evaluated):
        raise RuntimeError("read-mapping search did not evaluate any candidates")
    if not completed_search and not args.smoke:
        raise RuntimeError(
            "read-mapping search did not finish within the time budget; rerun with a higher "
            "--time-budget or a more selective search policy"
        )

    baseline = next(row for row in evaluated if row["config"] == MapperConfig())
    if baseline["val_metrics"] is None:
        raise RuntimeError("baseline mapper config did not survive to validation")
    if heldout_eval_root is not None and baseline["heldout_metrics"] is None:
        raise RuntimeError("baseline mapper config did not survive to heldout")

    best_train, best_selected, best_unconstrained, eligible_candidates = choose_finalist(
        evaluated=evaluated,
        baseline=baseline,
        max_heldout_scaffold_specificity_drop=args.max_heldout_scaffold_specificity_drop,
    )

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
        "grid_candidates": len(configs),
        "eligible_candidates": eligible_candidates,
        "completed_grid": completed_grid,
        "completed_search": completed_search,
        "selection_policy": {
            "mode": (
                "deterministic_grid_with_heldout_specificity_gate"
                if args.search_policy == "grid"
                else "successive_halving_with_heldout_specificity_gate"
            ),
            "search_policy": args.search_policy,
            "halving_factor": args.halving_factor if args.search_policy == "successive-halving" else None,
            "final_heldout_candidates": (
                final_heldout_candidates if args.search_policy == "successive-halving" else None
            ),
            "max_heldout_scaffold_specificity_drop": args.max_heldout_scaffold_specificity_drop,
        },
        "stage_records": stage_records,
        "search_seconds": search_seconds,
        "evaluated_candidates": sum(1 for row in evaluated if row["train_metrics"] is not None),
        "split_evaluations": split_evaluations,
        "total_split_evaluations": sum(split_evaluations.values()),
    }

    print(json.dumps(result, indent=2, sort_keys=True))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote read-mapping result to {args.output}")


if __name__ == "__main__":
    main()
