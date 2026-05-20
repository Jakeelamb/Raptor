#!/usr/bin/env python3
"""Run a tiny truth-known transcriptome fixture through current Raptor.

This is a fast-loop harness, not a parity claim. It creates a small fixture with
two related isoforms, runs the production Raptor CLI, and writes machine-readable
metrics that expose current gaps.
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "target" / "trinity_parity" / "tiny_alt_isoform"
DEFAULT_ORACLE = ROOT / "bench" / "trinity_parity" / "oracles" / "tiny_alt_isoform.fa"


def deterministic_dna(label: str, length: int) -> str:
    state = 0x9E3779B97F4A7C15
    for byte in label.encode("utf-8"):
        state ^= byte
        state = (state * 0xBF58476D1CE4E5B9) & ((1 << 64) - 1)

    bases = "ACGT"
    out: list[str] = []
    while len(out) < length:
        state ^= (state >> 12) & ((1 << 64) - 1)
        state ^= (state << 25) & ((1 << 64) - 1)
        state ^= (state >> 27) & ((1 << 64) - 1)
        value = (state * 0x2545F4914F6CDD1D) & ((1 << 64) - 1)
        base = bases[value & 3]
        if len(out) >= 3 and out[-1] == out[-2] == out[-3] == base:
            base = bases[(bases.index(base) + 1) & 3]
        out.append(base)
    return "".join(out)


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1].upper()


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_fastq_gz(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for name, seq in records:
            handle.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def tiled_starts(sequence_len: int, window_len: int, step: int) -> list[int]:
    if sequence_len < window_len:
        return []
    last = sequence_len - window_len
    starts = list(range(0, last + 1, step))
    if starts[-1] != last:
        starts.append(last)
    return starts


def read_fasta_lengths(path: Path) -> list[int]:
    opener = gzip.open if path.suffix == ".gz" else open
    lengths: list[int] = []
    current = 0
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line)
    if current:
        lengths.append(current)
    return lengths


def read_fasta_records(path: Path) -> dict[str, str]:
    opener = gzip.open if path.suffix == ".gz" else open
    records: dict[str, str] = {}
    name: str | None = None
    parts: list[str] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())
    if name is not None:
        records[name] = "".join(parts)
    return records


def longest_common_substring_len(a: str, b: str) -> int:
    if not a or not b:
        return 0
    previous = [0] * (len(b) + 1)
    best = 0
    for base_a in a:
        current = [0] * (len(b) + 1)
        for idx, base_b in enumerate(b, start=1):
            if base_a == base_b:
                value = previous[idx - 1] + 1
                current[idx] = value
                if value > best:
                    best = value
        previous = current
    return best


def fasta_recovery_metrics(reference_fasta: Path, assembly_fasta: Path) -> dict[str, object]:
    reference = read_fasta_records(reference_fasta)
    assembled = read_fasta_records(assembly_fasta)
    per_reference: dict[str, object] = {}
    for reference_name, reference_seq in reference.items():
        best = 0
        best_contig = None
        for contig_name, contig_seq in assembled.items():
            forward = longest_common_substring_len(reference_seq, contig_seq)
            reverse = longest_common_substring_len(reference_seq, revcomp(contig_seq))
            observed = max(forward, reverse)
            if observed > best:
                best = observed
                best_contig = contig_name
        coverage = best / len(reference_seq) if reference_seq else 0.0
        per_reference[reference_name] = {
            "reference_length": len(reference_seq),
            "best_matching_bases": best,
            "best_contig": best_contig,
            "best_coverage": round(coverage, 6),
        }

    coverages = [
        entry["best_coverage"]
        for entry in per_reference.values()
        if isinstance(entry, dict)
    ]
    return {
        "reference_transcript_count": len(reference),
        "assembled_record_count": len(assembled),
        "mean_best_coverage": round(sum(coverages) / len(coverages), 6)
        if coverages
        else 0.0,
        "min_best_coverage": min(coverages) if coverages else 0.0,
        "per_reference": per_reference,
    }


def truth_recovery_metrics(truth_fasta: Path, assembly_fasta: Path) -> dict[str, object]:
    metrics = fasta_recovery_metrics(truth_fasta, assembly_fasta)
    metrics["truth_transcript_count"] = metrics.pop("reference_transcript_count")
    metrics["per_transcript"] = metrics.pop("per_reference")
    for entry in metrics["per_transcript"].values():
        if isinstance(entry, dict):
            entry["truth_length"] = entry.pop("reference_length")
    return metrics


def n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    total = sum(lengths)
    acc = 0
    for length in sorted(lengths, reverse=True):
        acc += length
        if acc * 2 >= total:
            return length
    return 0


def run_command(command: list[str], cwd: Path) -> dict[str, object]:
    start = time.monotonic()
    completed = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    elapsed = time.monotonic() - start
    return {
        "command": command,
        "exit_code": completed.returncode,
        "elapsed_seconds": round(elapsed, 6),
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def generate_fixture(out_dir: Path) -> dict[str, object]:
    exon_a = deterministic_dna("shared_exon_a", 90)
    exon_b = deterministic_dna("dominant_exon_b", 72)
    exon_alt = deterministic_dna("alternative_exon", 60)
    exon_c = deterministic_dna("shared_exon_c", 90)

    tx1 = exon_a + exon_b + exon_c
    tx2 = exon_a + exon_alt + exon_c
    transcripts = {"tx_dominant": tx1, "tx_alt": tx2}

    truth_fasta = out_dir / "truth" / "transcripts.fa"
    truth_lines = []
    for name, seq in transcripts.items():
        truth_lines.append(f">{name} length={len(seq)}")
        truth_lines.extend(seq[i : i + 80] for i in range(0, len(seq), 80))
    write_text(truth_fasta, "\n".join(truth_lines) + "\n")

    single_records: list[tuple[str, str]] = []
    r1_records: list[tuple[str, str]] = []
    r2_records: list[tuple[str, str]] = []
    read_len = 75
    insert = 160
    step = 24

    for tx_name, seq in transcripts.items():
        coverage_rounds = 3 if tx_name == "tx_dominant" else 2
        for round_idx in range(coverage_rounds):
            for start in tiled_starts(len(seq), read_len, step):
                read = seq[start : start + read_len]
                single_records.append((f"{tx_name}_se_{round_idx}_{start}", read))
            for start in tiled_starts(len(seq), insert, step * 2):
                frag = seq[start : start + insert]
                r1_records.append((f"{tx_name}_pe_{round_idx}_{start}/1", frag[:read_len]))
                r2_records.append((f"{tx_name}_pe_{round_idx}_{start}/2", revcomp(frag[-read_len:])))

    single_fastq = out_dir / "reads" / "single.fastq.gz"
    r1_fastq = out_dir / "reads" / "R1.fastq.gz"
    r2_fastq = out_dir / "reads" / "R2.fastq.gz"
    write_fastq_gz(single_fastq, single_records)
    write_fastq_gz(r1_fastq, r1_records)
    write_fastq_gz(r2_fastq, r2_records)

    metadata = {
        "fixture": "tiny_alt_isoform",
        "read_len": read_len,
        "insert": insert,
        "truth_transcripts": {name: len(seq) for name, seq in transcripts.items()},
        "single_end_reads": len(single_records),
        "paired_end_pairs": len(r1_records),
        "paths": {
            "truth_fasta": str(truth_fasta),
            "single_fastq": str(single_fastq),
            "r1_fastq": str(r1_fastq),
            "r2_fastq": str(r2_fastq),
        },
    }
    write_text(out_dir / "truth" / "metadata.json", json.dumps(metadata, indent=2) + "\n")
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--skip-raptor", action="store_true")
    parser.add_argument("--run-trinity", action="store_true")
    parser.add_argument("--oracle-fasta", type=Path, default=DEFAULT_ORACLE)
    parser.add_argument("--min-truth-coverage", type=float, default=0.95)
    parser.add_argument("--min-oracle-coverage", type=float, default=0.95)
    args = parser.parse_args()

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    fixture = generate_fixture(out_dir)

    report: dict[str, object] = {
        "fixture": fixture,
        "repo": str(ROOT),
        "raptor": None,
        "trinity": {
            "available": shutil.which("Trinity") is not None,
            "ran": False,
            "note": "Trinity is optional for this scaffold until the benchmark panel is frozen.",
        },
        "oracle": {
            "path": str(args.oracle_fasta),
            "available": args.oracle_fasta.exists(),
            "note": "Checked-in tiny oracle is a frozen stand-in until Trinity is installed or a Trinity oracle is captured.",
        },
        "known_limitations": [
            "Current raptor assemble CLI accepts one input FASTQ, so this scaffold runs single-end reads only.",
            "No Trinity parity claim is made from this tiny fixture.",
        ],
    }

    if not args.skip_raptor:
        output_fasta = out_dir / "raptor" / "assembly.fa.gz"
        output_fasta.parent.mkdir(parents=True, exist_ok=True)
        command = [
            "cargo",
            "run",
            "--quiet",
            "--",
            "assemble",
            "--input",
            fixture["paths"]["single_fastq"],
            "--output",
            str(output_fasta),
            "--min-len",
            "25",
            "--threads",
            "1",
        ]
        result = run_command(command, ROOT)
        metrics = {"output_exists": output_fasta.exists()}
        if output_fasta.exists():
            lengths = read_fasta_lengths(output_fasta)
            metrics.update(
                {
                    "transcript_count": len(lengths),
                    "total_bases": sum(lengths),
                    "n50": n50(lengths),
                    "lengths": lengths,
                }
            )
            metrics["truth_recovery"] = truth_recovery_metrics(
                Path(fixture["paths"]["truth_fasta"]), output_fasta
            )
            if args.oracle_fasta.exists():
                metrics["oracle_recovery"] = fasta_recovery_metrics(
                    args.oracle_fasta, output_fasta
                )
        result["metrics"] = metrics
        report["raptor"] = result

    if args.run_trinity:
        trinity_bin = shutil.which("Trinity")
        if not trinity_bin:
            report["trinity"]["ran"] = False
            report["trinity"]["error"] = "Trinity executable not found on PATH"
        else:
            trinity_out = out_dir / "trinity"
            command = [
                trinity_bin,
                "--seqType",
                "fq",
                "--single",
                fixture["paths"]["single_fastq"],
                "--CPU",
                "1",
                "--max_memory",
                "2G",
                "--output",
                str(trinity_out),
            ]
            report["trinity"]["ran"] = True
            report["trinity"]["result"] = run_command(command, ROOT)

    report_path = out_dir / "report.json"
    write_text(report_path, json.dumps(report, indent=2) + "\n")

    print(f"wrote {report_path}")
    if report["raptor"] and report["raptor"]["exit_code"] != 0:
        print("raptor assemble failed; see report.json", file=sys.stderr)
        return 1
    if report["raptor"]:
        metrics = report["raptor"].get("metrics", {})
        truth_recovery = metrics.get("truth_recovery", {})
        min_coverage = truth_recovery.get("min_best_coverage", 0.0)
        if min_coverage < args.min_truth_coverage:
            print(
                "truth recovery below threshold: "
                f"{min_coverage} < {args.min_truth_coverage}",
                file=sys.stderr,
            )
            return 1
        oracle_recovery = metrics.get("oracle_recovery")
        if oracle_recovery:
            min_oracle_coverage = oracle_recovery.get("min_best_coverage", 0.0)
            if min_oracle_coverage < args.min_oracle_coverage:
                print(
                    "oracle recovery below threshold: "
                    f"{min_oracle_coverage} < {args.min_oracle_coverage}",
                    file=sys.stderr,
                )
                return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
