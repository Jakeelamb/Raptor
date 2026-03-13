#!/usr/bin/env python3
"""Prepare fixed error-correction train/val/heldout panels for the Raptor bridge."""

from __future__ import annotations

import argparse
import json
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BINARY = REPO_ROOT / "target" / "debug" / "raptor"
DEFAULT_OUTPUT_ROOT = (
    REPO_ROOT / "artifacts" / "autoresearch_raptor" / "error_correction" / "tasks"
)


@dataclass(frozen=True)
class TaskProfile:
    name: str
    trusted_roots: int
    weak_roots: int
    correctable_singletons_per_trusted: int
    protected_singletons_per_weak: int


TRAIN_VAL_PROFILES = (
    TaskProfile(
        name="balanced",
        trusted_roots=64,
        weak_roots=32,
        correctable_singletons_per_trusted=2,
        protected_singletons_per_weak=1,
    ),
    TaskProfile(
        name="correction_heavy",
        trusted_roots=72,
        weak_roots=24,
        correctable_singletons_per_trusted=3,
        protected_singletons_per_weak=1,
    ),
    TaskProfile(
        name="retention_heavy",
        trusted_roots=48,
        weak_roots=56,
        correctable_singletons_per_trusted=1,
        protected_singletons_per_weak=2,
    ),
    TaskProfile(
        name="mixed_pressure",
        trusted_roots=80,
        weak_roots=40,
        correctable_singletons_per_trusted=2,
        protected_singletons_per_weak=2,
    ),
)

HELDOUT_PROFILES = (
    TaskProfile(
        name="heldout_retention_stress",
        trusted_roots=40,
        weak_roots=72,
        correctable_singletons_per_trusted=1,
        protected_singletons_per_weak=3,
    ),
    TaskProfile(
        name="heldout_tradeoff",
        trusted_roots=56,
        weak_roots=64,
        correctable_singletons_per_trusted=2,
        protected_singletons_per_weak=2,
    ),
    TaskProfile(
        name="heldout_correction_stress",
        trusted_roots=88,
        weak_roots=28,
        correctable_singletons_per_trusted=4,
        protected_singletons_per_weak=1,
    ),
)


def run_prepare_task(
    binary: Path,
    output_dir: Path,
    k: int,
    profile: TaskProfile,
    seed: int,
) -> None:
    cmd = [
        str(binary),
        "component-bench",
        "prepare-error-correction",
        "--output",
        str(output_dir),
        "--k",
        str(k),
        "--trusted-roots",
        str(profile.trusted_roots),
        "--weak-roots",
        str(profile.weak_roots),
        "--correctable-singletons-per-trusted",
        str(profile.correctable_singletons_per_trusted),
        "--protected-singletons-per-weak",
        str(profile.protected_singletons_per_weak),
        "--seed",
        str(seed),
    ]
    subprocess.run(cmd, check=True)


def prepare_split(
    binary: Path,
    output_dir: Path,
    tasks: int,
    k: int,
    seed: int,
    seed_step: int,
    profiles: tuple[TaskProfile, ...],
) -> list[dict[str, int | str]]:
    tasks_written: list[dict[str, int | str]] = []
    output_dir.mkdir(parents=True, exist_ok=True)

    for task_index in range(tasks):
        profile = profiles[task_index % len(profiles)]
        task_seed = seed + task_index * seed_step
        task_dir = output_dir / f"task_{task_index:03d}"
        run_prepare_task(
            binary=binary,
            output_dir=task_dir,
            k=k,
            profile=profile,
            seed=task_seed,
        )
        tasks_written.append(
            {
                "task_name": task_dir.name,
                "profile": profile.name,
                "trusted_roots": profile.trusted_roots,
                "weak_roots": profile.weak_roots,
                "correctable_singletons_per_trusted": profile.correctable_singletons_per_trusted,
                "protected_singletons_per_weak": profile.protected_singletons_per_weak,
                "seed": task_seed,
            }
        )

    return tasks_written


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare error-correction train/val/heldout panels")
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--train-tasks", type=int, default=16)
    parser.add_argument("--val-tasks", type=int, default=8)
    parser.add_argument("--heldout-tasks", type=int, default=8)
    parser.add_argument("--k", type=int, default=21)
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

    train_tasks = prepare_split(
        binary=args.binary,
        output_dir=train_root,
        tasks=args.train_tasks,
        k=args.k,
        seed=args.seed,
        seed_step=args.seed_step,
        profiles=TRAIN_VAL_PROFILES,
    )
    val_tasks = prepare_split(
        binary=args.binary,
        output_dir=val_root,
        tasks=args.val_tasks,
        k=args.k,
        seed=args.seed + args.train_tasks * args.seed_step + 1000,
        seed_step=args.seed_step,
        profiles=TRAIN_VAL_PROFILES,
    )
    heldout_tasks = prepare_split(
        binary=args.binary,
        output_dir=heldout_root,
        tasks=args.heldout_tasks,
        k=args.k,
        seed=args.seed + (args.train_tasks + args.val_tasks) * args.seed_step + 2000,
        seed_step=args.seed_step,
        profiles=HELDOUT_PROFILES,
    )

    metadata = {
        "binary": str(args.binary),
        "train_root": str(train_root),
        "val_root": str(val_root),
        "heldout_root": str(heldout_root),
        "train_tasks": args.train_tasks,
        "val_tasks": args.val_tasks,
        "heldout_tasks": args.heldout_tasks,
        "k": args.k,
        "seed": args.seed,
        "seed_step": args.seed_step,
        "train_val_profiles": [asdict(profile) for profile in TRAIN_VAL_PROFILES],
        "heldout_profiles": [asdict(profile) for profile in HELDOUT_PROFILES],
        "train_tasks_written": train_tasks,
        "val_tasks_written": val_tasks,
        "heldout_tasks_written": heldout_tasks,
    }
    metadata_path = args.output_root / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote error-correction train/val/heldout panels under {args.output_root}")


if __name__ == "__main__":
    main()
