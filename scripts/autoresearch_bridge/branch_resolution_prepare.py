#!/usr/bin/env python3
"""Prepare fixed branch-resolution train/val/heldout panels for the Raptor bridge."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BINARY = REPO_ROOT / "target" / "debug" / "raptor"
DEFAULT_OUTPUT_ROOT = (
    REPO_ROOT / "artifacts" / "autoresearch_raptor" / "branch_resolution" / "tasks"
)


def run_prepare_panel(
    binary: Path,
    output_dir: Path,
    tasks: int,
    cases_per_task: int,
    seed: int,
    seed_step: int,
    scenario_profile: str,
) -> None:
    cmd = [
        str(binary),
        "component-bench",
        "prepare-branch-resolution-panel",
        "--output",
        str(output_dir),
        "--tasks",
        str(tasks),
        "--cases-per-task",
        str(cases_per_task),
        "--seed",
        str(seed),
        "--seed-step",
        str(seed_step),
        "--scenario-profile",
        scenario_profile,
    ]
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare branch-resolution train/val/heldout panels")
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--train-tasks", type=int, default=16)
    parser.add_argument("--val-tasks", type=int, default=8)
    parser.add_argument("--heldout-tasks", type=int, default=8)
    parser.add_argument("--cases-per-task", type=int, default=64)
    parser.add_argument("--heldout-cases-per-task", type=int, default=96)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--seed-step", type=int, default=1)
    args = parser.parse_args()

    if not args.binary.exists():
        raise FileNotFoundError(
            f"missing raptor binary at {args.binary}; build it first with cargo build"
        )

    train_root = args.output_root / "train"
    val_root = args.output_root / "val"
    heldout_root = args.output_root / "heldout"
    train_root.mkdir(parents=True, exist_ok=True)
    val_root.mkdir(parents=True, exist_ok=True)
    heldout_root.mkdir(parents=True, exist_ok=True)

    run_prepare_panel(
        binary=args.binary,
        output_dir=train_root,
        tasks=args.train_tasks,
        cases_per_task=args.cases_per_task,
        seed=args.seed,
        seed_step=args.seed_step,
        scenario_profile="balanced",
    )
    run_prepare_panel(
        binary=args.binary,
        output_dir=val_root,
        tasks=args.val_tasks,
        cases_per_task=args.cases_per_task,
        seed=args.seed + args.train_tasks * args.seed_step + 1000,
        seed_step=args.seed_step,
        scenario_profile="balanced",
    )
    run_prepare_panel(
        binary=args.binary,
        output_dir=heldout_root,
        tasks=args.heldout_tasks,
        cases_per_task=args.heldout_cases_per_task,
        seed=args.seed + (args.train_tasks + args.val_tasks) * args.seed_step + 2000,
        seed_step=args.seed_step,
        scenario_profile="support_stress",
    )

    metadata = {
        "binary": str(args.binary),
        "train_root": str(train_root),
        "val_root": str(val_root),
        "heldout_root": str(heldout_root),
        "train_tasks": args.train_tasks,
        "val_tasks": args.val_tasks,
        "heldout_tasks": args.heldout_tasks,
        "cases_per_task": args.cases_per_task,
        "heldout_cases_per_task": args.heldout_cases_per_task,
        "seed": args.seed,
        "seed_step": args.seed_step,
    }
    metadata_path = args.output_root / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote branch-resolution train/val/heldout panels under {args.output_root}")


if __name__ == "__main__":
    main()
