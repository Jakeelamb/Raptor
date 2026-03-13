#!/usr/bin/env python3
"""Prepare fixed contig-extraction train/val/heldout panels for the Raptor bridge."""

from __future__ import annotations

import argparse
import json
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BINARY = REPO_ROOT / "target" / "debug" / "raptor"
DEFAULT_OUTPUT_ROOT = (
    REPO_ROOT / "artifacts" / "autoresearch_raptor" / "contig_extraction" / "tasks"
)


@dataclass(frozen=True)
class TaskProfile:
    name: str
    profile: str
    component_count: int
    primary_reads_per_component: int
    alternate_reads_per_component: int


TRAIN_VAL_PROFILES = (
    TaskProfile(
        name="branching_dominant",
        profile="branching",
        component_count=6,
        primary_reads_per_component=6,
        alternate_reads_per_component=2,
    ),
    TaskProfile(
        name="branching_close",
        profile="branching",
        component_count=6,
        primary_reads_per_component=4,
        alternate_reads_per_component=3,
    ),
    TaskProfile(
        name="repeat_completion",
        profile="repeat_completion",
        component_count=4,
        primary_reads_per_component=4,
        alternate_reads_per_component=3,
    ),
    TaskProfile(
        name="repeat_fallback",
        profile="repeat_fallback",
        component_count=4,
        primary_reads_per_component=4,
        alternate_reads_per_component=3,
    ),
)

HELDOUT_PROFILES = (
    TaskProfile(
        name="heldout_branching_close",
        profile="branching",
        component_count=8,
        primary_reads_per_component=4,
        alternate_reads_per_component=3,
    ),
    TaskProfile(
        name="heldout_branching_edge",
        profile="branching",
        component_count=6,
        primary_reads_per_component=5,
        alternate_reads_per_component=4,
    ),
    TaskProfile(
        name="heldout_repeat_completion",
        profile="repeat_completion",
        component_count=4,
        primary_reads_per_component=4,
        alternate_reads_per_component=3,
    ),
    TaskProfile(
        name="heldout_repeat_fallback",
        profile="repeat_fallback",
        component_count=4,
        primary_reads_per_component=4,
        alternate_reads_per_component=3,
    ),
    TaskProfile(
        name="heldout_repeat_priority",
        profile="repeat_priority",
        component_count=4,
        primary_reads_per_component=4,
        alternate_reads_per_component=3,
    ),
    TaskProfile(
        name="heldout_repeat_priority_stress",
        profile="repeat_priority",
        component_count=4,
        primary_reads_per_component=5,
        alternate_reads_per_component=4,
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
        "prepare-contig-extraction",
        "--output",
        str(output_dir),
        "--profile",
        profile.profile,
        "--k",
        str(k),
        "--component-count",
        str(profile.component_count),
        "--primary-reads-per-component",
        str(profile.primary_reads_per_component),
        "--alternate-reads-per-component",
        str(profile.alternate_reads_per_component),
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
                "profile_name": profile.name,
                "profile": profile.profile,
                "component_count": profile.component_count,
                "primary_reads_per_component": profile.primary_reads_per_component,
                "alternate_reads_per_component": profile.alternate_reads_per_component,
                "seed": task_seed,
            }
        )

    return tasks_written


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare contig-extraction train/val/heldout panels")
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--train-tasks", type=int, default=16)
    parser.add_argument("--val-tasks", type=int, default=8)
    parser.add_argument("--heldout-tasks", type=int, default=8)
    parser.add_argument("--k", type=int, default=11)
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
    print(f"Wrote contig-extraction train/val/heldout panels under {args.output_root}")


if __name__ == "__main__":
    main()
