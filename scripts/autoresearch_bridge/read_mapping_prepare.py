#!/usr/bin/env python3
"""Prepare fixed read-mapping train/val/heldout panels for the Raptor autoresearch bridge."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BINARY = REPO_ROOT / "target" / "debug" / "raptor"
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "artifacts" / "autoresearch_raptor" / "read_mapping" / "tasks"


def run_prepare_panel(
    binary: Path,
    output_dir: Path,
    tasks: int,
    paired: bool,
    num_contigs: int,
    contig_len: int,
    read_len: int,
    reads: int,
    insert_size: int,
    repeat_len: int,
    error_rate: float,
    decoy_rate: float,
    ambiguous_repeat_decoy_rate: float,
    seed: int,
    seed_step: int,
    minimizer_k: int,
    minimizer_w: int,
) -> None:
    cmd = [
        str(binary),
        "component-bench",
        "prepare-read-mapping-panel",
        "--output",
        str(output_dir),
        "--tasks",
        str(tasks),
        "--num-contigs",
        str(num_contigs),
        "--contig-len",
        str(contig_len),
        "--read-len",
        str(read_len),
        "--reads",
        str(reads),
        "--insert-size",
        str(insert_size),
        "--repeat-len",
        str(repeat_len),
        "--error-rate",
        str(error_rate),
        "--decoy-rate",
        str(decoy_rate),
        "--ambiguous-repeat-decoy-rate",
        str(ambiguous_repeat_decoy_rate),
        "--seed",
        str(seed),
        "--seed-step",
        str(seed_step),
        "--k",
        str(minimizer_k),
        "--w",
        str(minimizer_w),
    ]
    if paired:
        cmd.append("--paired")
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare read-mapping train/val/heldout panels")
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--train-tasks", type=int, default=16)
    parser.add_argument("--val-tasks", type=int, default=8)
    parser.add_argument("--heldout-tasks", type=int, default=4)
    parser.add_argument("--paired", action="store_true")
    parser.add_argument("--num-contigs", type=int, default=5)
    parser.add_argument("--contig-len", type=int, default=2200)
    parser.add_argument("--read-len", type=int, default=90)
    parser.add_argument("--reads", type=int, default=64)
    parser.add_argument("--insert-size", type=int, default=250)
    parser.add_argument("--repeat-len", type=int, default=120)
    parser.add_argument("--error-rate", type=float, default=0.04)
    parser.add_argument("--decoy-rate", type=float, default=0.25)
    parser.add_argument("--ambiguous-repeat-decoy-rate", type=float, default=0.15)
    parser.add_argument("--heldout-decoy-rate", type=float, default=0.35)
    parser.add_argument("--heldout-ambiguous-repeat-decoy-rate", type=float, default=0.30)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--seed-step", type=int, default=1)
    parser.add_argument("--k", type=int, default=11)
    parser.add_argument("--w", type=int, default=5)
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

    val_seed = args.seed + args.train_tasks * args.seed_step + 1000
    heldout_seed = args.seed + (args.train_tasks + args.val_tasks) * args.seed_step + 2000

    run_prepare_panel(
        binary=args.binary,
        output_dir=train_root,
        tasks=args.train_tasks,
        paired=args.paired,
        num_contigs=args.num_contigs,
        contig_len=args.contig_len,
        read_len=args.read_len,
        reads=args.reads,
        insert_size=args.insert_size,
        repeat_len=args.repeat_len,
        error_rate=args.error_rate,
        decoy_rate=args.decoy_rate,
        ambiguous_repeat_decoy_rate=args.ambiguous_repeat_decoy_rate,
        seed=args.seed,
        seed_step=args.seed_step,
        minimizer_k=args.k,
        minimizer_w=args.w,
    )
    run_prepare_panel(
        binary=args.binary,
        output_dir=val_root,
        tasks=args.val_tasks,
        paired=args.paired,
        num_contigs=args.num_contigs,
        contig_len=args.contig_len,
        read_len=args.read_len,
        reads=args.reads,
        insert_size=args.insert_size,
        repeat_len=args.repeat_len,
        error_rate=args.error_rate,
        decoy_rate=args.decoy_rate,
        ambiguous_repeat_decoy_rate=args.ambiguous_repeat_decoy_rate,
        seed=val_seed,
        seed_step=args.seed_step,
        minimizer_k=args.k,
        minimizer_w=args.w,
    )
    run_prepare_panel(
        binary=args.binary,
        output_dir=heldout_root,
        tasks=args.heldout_tasks,
        paired=args.paired,
        num_contigs=args.num_contigs,
        contig_len=args.contig_len,
        read_len=args.read_len,
        reads=args.reads,
        insert_size=args.insert_size,
        repeat_len=args.repeat_len,
        error_rate=args.error_rate,
        decoy_rate=args.heldout_decoy_rate,
        ambiguous_repeat_decoy_rate=args.heldout_ambiguous_repeat_decoy_rate,
        seed=heldout_seed,
        seed_step=args.seed_step,
        minimizer_k=args.k,
        minimizer_w=args.w,
    )

    metadata = {
        "binary": str(args.binary),
        "train_root": str(train_root),
        "val_root": str(val_root),
        "heldout_root": str(heldout_root),
        "train_tasks": args.train_tasks,
        "val_tasks": args.val_tasks,
        "heldout_tasks": args.heldout_tasks,
        "paired": args.paired,
        "num_contigs": args.num_contigs,
        "contig_len": args.contig_len,
        "read_len": args.read_len,
        "reads": args.reads,
        "insert_size": args.insert_size,
        "repeat_len": args.repeat_len,
        "error_rate": args.error_rate,
        "decoy_rate": args.decoy_rate,
        "ambiguous_repeat_decoy_rate": args.ambiguous_repeat_decoy_rate,
        "heldout_decoy_rate": args.heldout_decoy_rate,
        "heldout_ambiguous_repeat_decoy_rate": args.heldout_ambiguous_repeat_decoy_rate,
        "seed": args.seed,
        "seed_step": args.seed_step,
        "task_k": args.k,
        "task_w": args.w,
    }
    metadata_path = args.output_root / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote read-mapping train/val/heldout panels under {args.output_root}")


if __name__ == "__main__":
    main()
