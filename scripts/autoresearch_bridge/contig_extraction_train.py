#!/usr/bin/env python3
"""Heldout-aware exhaustive search for contig-extraction component tasks."""

from __future__ import annotations

import argparse
import json
import itertools
import time
from dataclasses import asdict, dataclass
from pathlib import Path

from bridge_common import run_json_command


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BINARY = REPO_ROOT / "target" / "debug" / "raptor"
DEFAULT_TASKS_ROOT = (
    REPO_ROOT / "artifacts" / "autoresearch_raptor" / "contig_extraction" / "tasks"
)
DEFAULT_OUTPUT = (
    REPO_ROOT / "artifacts" / "autoresearch_raptor" / "contig_extraction" / "latest_result.json"
)


@dataclass(frozen=True)
class ContigExtractionConfig:
    prefer_high_count_seeds: bool = True
    prefer_non_repeat_seeds: bool = True
    enable_repeat_seed_completion: bool = True
    suppress_redundant_contigs: bool = False


PROTECTED_HELDOUT_METRICS = (
    "exact_contig_rate",
    "truth_kmer_f1",
)


def safe_metric(metrics: dict | None, key: str) -> float | None:
    if metrics is None:
        return None
    value = metrics.get(key)
    if value is None:
        return None
    return float(value)


def baseline_distance(config: ContigExtractionConfig) -> int:
    baseline = ContigExtractionConfig()
    return sum(
        (
            config.prefer_high_count_seeds != baseline.prefer_high_count_seeds,
            config.prefer_non_repeat_seeds != baseline.prefer_non_repeat_seeds,
            config.enable_repeat_seed_completion != baseline.enable_repeat_seed_completion,
            config.suppress_redundant_contigs != baseline.suppress_redundant_contigs,
        )
    )


def evaluate_config(binary: Path, task_root: Path, config: ContigExtractionConfig) -> dict:
    cmd = [
        str(binary),
        "component-bench",
        "contig-extraction",
        "--task",
        str(task_root),
        "--json",
    ]
    if not config.prefer_high_count_seeds:
        cmd.append("--disable-prefer-high-count-seeds")
    if not config.prefer_non_repeat_seeds:
        cmd.append("--disable-prefer-non-repeat-seeds")
    if not config.enable_repeat_seed_completion:
        cmd.append("--disable-repeat-seed-completion")
    if config.suppress_redundant_contigs:
        cmd.append("--suppress-redundant-contigs")
    return run_json_command(cmd)["summary"]


def enumerate_configs() -> list[ContigExtractionConfig]:
    configs = [
        ContigExtractionConfig(
            prefer_high_count_seeds=prefer_high_count_seeds,
            prefer_non_repeat_seeds=prefer_non_repeat_seeds,
            enable_repeat_seed_completion=enable_repeat_seed_completion,
            suppress_redundant_contigs=suppress_redundant_contigs,
        )
        for prefer_high_count_seeds, prefer_non_repeat_seeds, enable_repeat_seed_completion, suppress_redundant_contigs in itertools.product(
            (True, False),
            (True, False),
            (True, False),
            (False, True),
        )
    ]
    default = ContigExtractionConfig()
    configs.remove(default)
    return [default, *configs]


def selection_key(row: dict) -> tuple[float, float, float, int, float]:
    val_metrics = row["val_metrics"]
    heldout_metrics = row.get("heldout_metrics")
    val_score = val_metrics["score"]
    heldout_score = heldout_metrics["score"] if heldout_metrics is not None else val_score
    robust_score = min(val_score, heldout_score)
    heldout_exact = safe_metric(heldout_metrics, "exact_contig_rate") or 0.0
    heldout_kmer_f1 = safe_metric(heldout_metrics, "truth_kmer_f1") or 0.0
    minimal_change_bias = -baseline_distance(row["config"])
    return (
        robust_score,
        heldout_exact,
        heldout_kmer_f1,
        minimal_change_bias,
        val_metrics["graph_nodes_per_second"],
    )


def evaluate_grid(
    binary: Path,
    train_root: Path,
    val_root: Path,
    heldout_root: Path | None,
) -> tuple[list[dict], float]:
    start = time.time()
    evaluated: list[dict] = []

    configs = enumerate_configs()
    for index, config in enumerate(configs, start=1):
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
        leader = max(evaluated, key=selection_key)
        heldout_exact = safe_metric(leader.get("heldout_metrics"), "exact_contig_rate")
        print(
            f"candidates={index:03d}/{len(configs):03d} | "
            f"best_val_score={leader['val_metrics']['score']:.4f} | "
            f"best_heldout_exact={heldout_exact if heldout_exact is not None else float('nan'):.4f} | "
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
    for metric in PROTECTED_HELDOUT_METRICS:
        candidate_value = safe_metric(candidate_heldout, metric)
        baseline_value = safe_metric(baseline_heldout, metric)
        if candidate_value is None or baseline_value is None:
            continue
        if candidate_value < baseline_value:
            return False
    return True


def choose_finalists(evaluated: list[dict]) -> tuple[dict, dict, dict, int]:
    baseline = next(row for row in evaluated if row["config"] == ContigExtractionConfig())
    best_train = max(
        evaluated,
        key=lambda row: (
            row["train_metrics"]["score"],
            row["train_metrics"]["graph_nodes_per_second"],
        ),
    )
    best_unconstrained = max(evaluated, key=selection_key)
    eligible = [row for row in evaluated if heldout_non_regressing(row, baseline)]
    if not eligible:
        eligible = [baseline]
    best_selected = max(eligible, key=selection_key)
    return best_train, best_selected, best_unconstrained, len(eligible)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the contig-extraction autoresearch bridge")
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
            f"expected train/val panels under {args.tasks_root}; run contig_extraction_prepare.py first"
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
    baseline = next(row for row in evaluated if row["config"] == ContigExtractionConfig())

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
            "protected_metrics": list(PROTECTED_HELDOUT_METRICS),
        },
        "search_seconds": search_seconds,
        "evaluated_candidates": len(evaluated),
    }

    print(json.dumps(result, indent=2, sort_keys=True))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote contig-extraction result to {args.output}")


if __name__ == "__main__":
    main()
