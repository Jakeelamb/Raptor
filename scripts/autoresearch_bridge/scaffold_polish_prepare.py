#!/usr/bin/env python3
"""Prepare fixed scaffold+polish train/val/heldout panels for the Raptor bridge."""

from __future__ import annotations

import argparse
import json
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BINARY = REPO_ROOT / "target" / "debug" / "raptor"
DEFAULT_OUTPUT_ROOT = (
    REPO_ROOT / "artifacts" / "autoresearch_raptor" / "scaffold_polish" / "tasks"
)


@dataclass(frozen=True)
class TaskProfile:
    name: str
    true_link_pairs: int
    decoy_link_pairs: int
    scaffold_groups: int | None = None


TRAIN_VAL_PROFILES = (
    TaskProfile(name="balanced", true_link_pairs=5, decoy_link_pairs=2),
    TaskProfile(name="sparse_truth", true_link_pairs=3, decoy_link_pairs=2),
    TaskProfile(name="recall_edge", true_link_pairs=4, decoy_link_pairs=2),
    TaskProfile(name="precision_edge", true_link_pairs=5, decoy_link_pairs=4),
    TaskProfile(name="tradeoff", true_link_pairs=4, decoy_link_pairs=4),
)

HELDOUT_PROFILES = (
    TaskProfile(name="heldout_balanced", true_link_pairs=5, decoy_link_pairs=3, scaffold_groups=3),
    TaskProfile(name="heldout_sparse_truth", true_link_pairs=3, decoy_link_pairs=2, scaffold_groups=3),
    TaskProfile(name="heldout_recall_edge", true_link_pairs=4, decoy_link_pairs=3, scaffold_groups=3),
    TaskProfile(name="heldout_precision_edge", true_link_pairs=5, decoy_link_pairs=4, scaffold_groups=3),
    TaskProfile(name="heldout_tradeoff", true_link_pairs=4, decoy_link_pairs=4, scaffold_groups=3),
)


def run_prepare_task(
    binary: Path,
    output_dir: Path,
    scaffold_groups: int,
    contigs_per_scaffold: int,
    contig_len: int,
    read_len: int,
    insert_size: int,
    internal_pairs_per_contig: int,
    true_link_pairs: int,
    decoy_link_pairs: int,
    mutations_per_contig: int,
    seed: int,
) -> None:
    cmd = [
        str(binary),
        "component-bench",
        "prepare-scaffold-polish",
        "--output",
        str(output_dir),
        "--scaffold-groups",
        str(scaffold_groups),
        "--contigs-per-scaffold",
        str(contigs_per_scaffold),
        "--contig-len",
        str(contig_len),
        "--read-len",
        str(read_len),
        "--insert-size",
        str(insert_size),
        "--internal-pairs-per-contig",
        str(internal_pairs_per_contig),
        "--true-link-pairs",
        str(true_link_pairs),
        "--decoy-link-pairs",
        str(decoy_link_pairs),
        "--mutations-per-contig",
        str(mutations_per_contig),
        "--seed",
        str(seed),
    ]
    subprocess.run(cmd, check=True)


def prepare_split(
    binary: Path,
    output_dir: Path,
    tasks: int,
    scaffold_groups: int,
    contigs_per_scaffold: int,
    contig_len: int,
    read_len: int,
    insert_size: int,
    internal_pairs_per_contig: int,
    mutations_per_contig: int,
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
            scaffold_groups=profile.scaffold_groups or scaffold_groups,
            contigs_per_scaffold=contigs_per_scaffold,
            contig_len=contig_len,
            read_len=read_len,
            insert_size=insert_size,
            internal_pairs_per_contig=internal_pairs_per_contig,
            true_link_pairs=profile.true_link_pairs,
            decoy_link_pairs=profile.decoy_link_pairs,
            mutations_per_contig=mutations_per_contig,
            seed=task_seed,
        )
        tasks_written.append(
            {
                "task_name": task_dir.name,
                "profile": profile.name,
                "scaffold_groups": profile.scaffold_groups or scaffold_groups,
                "true_link_pairs": profile.true_link_pairs,
                "decoy_link_pairs": profile.decoy_link_pairs,
                "seed": task_seed,
            }
        )
    return tasks_written


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare scaffold+polish train/val/heldout panels")
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--train-tasks", type=int, default=16)
    parser.add_argument("--val-tasks", type=int, default=8)
    parser.add_argument("--heldout-tasks", type=int, default=8)
    parser.add_argument("--scaffold-groups", type=int, default=2)
    parser.add_argument("--contigs-per-scaffold", type=int, default=3)
    parser.add_argument("--contig-len", type=int, default=400)
    parser.add_argument("--read-len", type=int, default=80)
    parser.add_argument("--insert-size", type=int, default=200)
    parser.add_argument("--internal-pairs-per-contig", type=int, default=12)
    parser.add_argument("--true-link-pairs", type=int, default=5)
    parser.add_argument("--decoy-link-pairs", type=int, default=2)
    parser.add_argument("--heldout-scaffold-groups", type=int, default=3)
    parser.add_argument("--heldout-true-link-pairs", type=int, default=5)
    parser.add_argument("--heldout-decoy-link-pairs", type=int, default=2)
    parser.add_argument("--mutations-per-contig", type=int, default=1)
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
        scaffold_groups=args.scaffold_groups,
        contigs_per_scaffold=args.contigs_per_scaffold,
        contig_len=args.contig_len,
        read_len=args.read_len,
        insert_size=args.insert_size,
        internal_pairs_per_contig=args.internal_pairs_per_contig,
        mutations_per_contig=args.mutations_per_contig,
        seed=args.seed,
        seed_step=args.seed_step,
        profiles=TRAIN_VAL_PROFILES,
    )
    val_tasks = prepare_split(
        binary=args.binary,
        output_dir=val_root,
        tasks=args.val_tasks,
        scaffold_groups=args.scaffold_groups,
        contigs_per_scaffold=args.contigs_per_scaffold,
        contig_len=args.contig_len,
        read_len=args.read_len,
        insert_size=args.insert_size,
        internal_pairs_per_contig=args.internal_pairs_per_contig,
        mutations_per_contig=args.mutations_per_contig,
        seed=args.seed + args.train_tasks * args.seed_step + 1000,
        seed_step=args.seed_step,
        profiles=TRAIN_VAL_PROFILES,
    )
    heldout_tasks = prepare_split(
        binary=args.binary,
        output_dir=heldout_root,
        tasks=args.heldout_tasks,
        scaffold_groups=args.heldout_scaffold_groups,
        contigs_per_scaffold=args.contigs_per_scaffold,
        contig_len=args.contig_len,
        read_len=args.read_len,
        insert_size=args.insert_size,
        internal_pairs_per_contig=args.internal_pairs_per_contig,
        mutations_per_contig=args.mutations_per_contig,
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
        "scaffold_groups": args.scaffold_groups,
        "contigs_per_scaffold": args.contigs_per_scaffold,
        "heldout_scaffold_groups": args.heldout_scaffold_groups,
        "mutations_per_contig": args.mutations_per_contig,
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
    print(f"Wrote scaffold+polish train/val/heldout panels under {args.output_root}")


if __name__ == "__main__":
    main()
