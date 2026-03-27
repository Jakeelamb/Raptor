#!/usr/bin/env python3
"""Run assemble-large on quick_test and emit a structured timing summary."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import time
from pathlib import Path
from typing import Any

from bridge_common import format_command


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BINARY = REPO_ROOT / "target" / "release" / "raptor"
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "artifacts" / "autoresearch_raptor" / "campaign_runs"
QUICK_TEST_READS_1 = REPO_ROOT / "bench" / "genome_assembly" / "data" / "quick_test" / "reads_1.fastq.gz"
QUICK_TEST_READS_2 = REPO_ROOT / "bench" / "genome_assembly" / "data" / "quick_test" / "reads_2.fastq.gz"

CORE_TIMING_RE = re.compile(
    r"Phase timings \(s\): distribute=(?P<distribute>[0-9.]+), "
    r"count=(?P<count>[0-9.]+), error_correct=(?P<error_correct>[0-9.]+), "
    r"graph_clean=(?P<graph_clean>[0-9.]+), graph_analyze=(?P<graph_analyze>[0-9.]+), "
    r"branch_thread=(?P<branch_thread>[0-9.]+), contig_extract=(?P<contig_extract>[0-9.]+), "
    r"write=(?P<write>[0-9.]+), total=(?P<total>[0-9.]+)"
)
SHARED_TIMING_RE = re.compile(
    r"Shared scaffold\+polish timings \(s\): index=(?P<index>[0-9.]+), "
    r"map=(?P<map>[0-9.]+), scaffold_write=(?P<scaffold_write>[0-9.]+), "
    r"consensus=(?P<consensus>[0-9.]+), polish_write=(?P<polish_write>[0-9.]+), "
    r"total=(?P<total>[0-9.]+)"
)


def parse_first_int(pattern: str, text: str) -> int | None:
    match = re.search(pattern, text, re.MULTILINE)
    if match is None:
        return None
    return int(match.group(1))


def parse_named_float_block(regex: re.Pattern[str], text: str) -> dict[str, float] | None:
    matches = list(regex.finditer(text))
    if not matches:
        return None
    match = matches[-1]
    return {name: float(value) for name, value in match.groupdict().items()}


def build_command(binary: Path, threads: int, output_dir: Path, extra_args: list[str]) -> tuple[list[str], Path]:
    contigs_path = output_dir / "contigs.fa"
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
        "-k",
        "31",
        "--min-count",
        "0",
        "--scaffold",
        "--polish",
        "--compress-buckets",
    ]
    cmd.extend(extra_args)
    return cmd, contigs_path


def run_profile(binary: Path, threads: int, label: str, output_root: Path, extra_args: list[str]) -> dict[str, Any]:
    run_dir = output_root / label
    run_dir.mkdir(parents=True, exist_ok=True)
    cmd, contigs_path = build_command(binary, threads, run_dir, extra_args)

    start = time.time()
    completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
    elapsed = time.time() - start

    stdout_path = run_dir / "stdout.log"
    stderr_path = run_dir / "stderr.log"
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")

    combined = completed.stdout + "\n" + completed.stderr
    result = {
        "label": label,
        "command": cmd,
        "command_pretty": format_command(cmd),
        "elapsed_seconds": elapsed,
        "core_phase_timings": parse_named_float_block(CORE_TIMING_RE, combined),
        "shared_phase_timings": parse_named_float_block(SHARED_TIMING_RE, combined),
        "contigs": parse_first_int(r"^Contigs:\s+(\d+)$", combined),
        "total_length_bp": parse_first_int(r"^Total length:\s+(\d+)\s+bp$", combined),
        "contig_n50_bp": parse_first_int(r"^N50:\s+(\d+)\s+bp$", combined),
        "scaffolds": parse_first_int(r"^\s+Scaffolds:\s+(\d+)$", combined),
        "scaffold_n50_bp": parse_first_int(r"^\s+N50:\s+(\d+)\s+bp$", combined),
        "polish_corrections": parse_first_int(r"^\s+Corrections:\s+(\d+)$", combined),
        "stdout_log": str(stdout_path),
        "stderr_log": str(stderr_path),
        "contigs_path": str(contigs_path),
    }

    summary_path = run_dir / "summary.json"
    summary_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    result["summary_path"] = str(summary_path)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--label", required=True)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--extra-arg",
        action="append",
        default=[],
        help="Extra argument to append to the assemble-large command. Repeat for multiple args.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = run_profile(args.binary, args.threads, args.label, args.output_root, args.extra_arg)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
