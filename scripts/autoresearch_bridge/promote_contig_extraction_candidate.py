#!/usr/bin/env python3
"""Promote a contig-extraction candidate through heldout eval and repeat-heavy end-to-end suites."""

from __future__ import annotations

import argparse
import gzip
import json
import re
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from bridge_common import run_json_command


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BINARY = REPO_ROOT / "target" / "release" / "raptor"
DEFAULT_TASKS_ROOT = (
    REPO_ROOT / "artifacts" / "autoresearch_raptor" / "contig_extraction" / "tasks"
)
DEFAULT_OUTPUT_ROOT = (
    REPO_ROOT / "artifacts" / "autoresearch_raptor" / "contig_extraction" / "promotions"
)
DEFAULT_TRAIN_RESULT = (
    REPO_ROOT / "artifacts" / "autoresearch_raptor" / "contig_extraction" / "latest_result.json"
)
QUICK_TEST_READS_1 = (
    REPO_ROOT / "bench" / "genome_assembly" / "data" / "quick_test" / "reads_1.fastq.gz"
)
QUICK_TEST_READS_2 = (
    REPO_ROOT / "bench" / "genome_assembly" / "data" / "quick_test" / "reads_2.fastq.gz"
)
QUICK_TEST_REFERENCE = (
    REPO_ROOT / "bench" / "genome_assembly" / "data" / "quick_test" / "reference.fa"
)
REPEAT_HEAVY_SOURCE_READS_1 = (
    REPO_ROOT / "bench" / "genome_assembly" / "data" / "drosophila" / "reads_1.fastq.gz"
)
REPEAT_HEAVY_SOURCE_READS_2 = (
    REPO_ROOT / "bench" / "genome_assembly" / "data" / "drosophila" / "reads_2.fastq.gz"
)
REPEAT_HEAVY_REFERENCE = (
    REPO_ROOT / "bench" / "genome_assembly" / "data" / "drosophila" / "reference.fa"
)
REPEAT_HEAVY_SUBSET_ROOT = (
    REPO_ROOT
    / "artifacts"
    / "autoresearch_raptor"
    / "contig_extraction"
    / "repeat_heavy_subset"
)
DEFAULT_REPEAT_HEAVY_MAX_PAIRS = 2_000_000
DEFAULT_BASELINE_CONFIG = {
    "prefer_high_count_seeds": True,
    "prefer_non_repeat_seeds": True,
    "enable_repeat_seed_completion": True,
    "suppress_redundant_contigs": False,
}


@dataclass(frozen=True)
class ContigExtractionConfig:
    prefer_high_count_seeds: bool = True
    prefer_non_repeat_seeds: bool = True
    enable_repeat_seed_completion: bool = True
    suppress_redundant_contigs: bool = False


@dataclass(frozen=True)
class DatasetSpec:
    name: str
    reads_1: Path
    reads_2: Path
    reference: Path


@dataclass(frozen=True)
class PromotionGate:
    min_val_score_delta: float = 0.0
    min_heldout_score_delta: float = 0.0
    max_heldout_exact_contig_rate_drop: float = 0.0
    max_heldout_truth_kmer_f1_drop: float = 0.0
    max_heldout_count_agreement_drop: float = 0.0
    max_quick_test_runtime_regression_ratio: float = 0.05
    max_quick_test_runtime_regression_seconds: float = 5.0
    max_quick_test_contig_n50_drop_ratio: float = 0.0
    max_quick_test_scaffold_n50_drop_ratio: float = 0.0
    max_quick_test_reference_coverage_drop: float = 0.002
    max_quick_test_aligned_query_fraction_drop: float = 0.01
    max_repeat_heavy_runtime_regression_ratio: float = 0.10
    max_repeat_heavy_runtime_regression_seconds: float = 30.0
    max_repeat_heavy_contig_n50_drop_ratio: float = 0.05
    max_repeat_heavy_scaffold_n50_drop_ratio: float = 0.05
    max_repeat_heavy_reference_coverage_drop: float = 0.002
    max_repeat_heavy_aligned_query_fraction_drop: float = 0.01
    min_repeat_heavy_contig_reduction: int = 1
    min_repeat_heavy_multi_alignment_improvement: float = 0.01


def stream_fastq_record(handle: Any) -> list[str] | None:
    record = [handle.readline() for _ in range(4)]
    if not record[0]:
        return None
    if any(line == "" for line in record):
        raise ValueError("truncated FASTQ record while materializing repeat-heavy subset")
    return record


def ensure_repeat_heavy_subset(max_pairs: int) -> tuple[Path, Path]:
    subset_dir = REPEAT_HEAVY_SUBSET_ROOT / f"drosophila_{max_pairs:07d}_pairs"
    reads_1 = subset_dir / "reads_1.fastq.gz"
    reads_2 = subset_dir / "reads_2.fastq.gz"
    metadata_path = subset_dir / "metadata.json"
    if reads_1.exists() and reads_2.exists() and metadata_path.exists():
        return reads_1, reads_2

    subset_dir.mkdir(parents=True, exist_ok=True)
    tmp_reads_1 = subset_dir / "reads_1.fastq.gz.tmp"
    tmp_reads_2 = subset_dir / "reads_2.fastq.gz.tmp"
    written_pairs = 0

    with (
        gzip.open(REPEAT_HEAVY_SOURCE_READS_1, "rt", encoding="utf-8") as source_1,
        gzip.open(REPEAT_HEAVY_SOURCE_READS_2, "rt", encoding="utf-8") as source_2,
        gzip.open(tmp_reads_1, "wt", encoding="utf-8") as out_1,
        gzip.open(tmp_reads_2, "wt", encoding="utf-8") as out_2,
    ):
        while written_pairs < max_pairs:
            record_1 = stream_fastq_record(source_1)
            record_2 = stream_fastq_record(source_2)
            if record_1 is None or record_2 is None:
                break
            out_1.writelines(record_1)
            out_2.writelines(record_2)
            written_pairs += 1

    tmp_reads_1.replace(reads_1)
    tmp_reads_2.replace(reads_2)
    metadata = {
        "source_reads_1": str(REPEAT_HEAVY_SOURCE_READS_1),
        "source_reads_2": str(REPEAT_HEAVY_SOURCE_READS_2),
        "subset_pairs": written_pairs,
        "requested_pairs": max_pairs,
        "reference": str(REPEAT_HEAVY_REFERENCE),
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Prepared repeat-heavy subset with {written_pairs} read pairs under {subset_dir}")
    return reads_1, reads_2


def load_config(args: argparse.Namespace) -> ContigExtractionConfig:
    if args.train_result is not None and args.train_result.exists():
        result = json.loads(args.train_result.read_text(encoding="utf-8"))
        winner = result.get("best_val_config") or result.get("best_train_config")
        if winner:
            return ContigExtractionConfig(
                prefer_high_count_seeds=winner["prefer_high_count_seeds"],
                prefer_non_repeat_seeds=winner["prefer_non_repeat_seeds"],
                enable_repeat_seed_completion=winner["enable_repeat_seed_completion"],
                suppress_redundant_contigs=winner["suppress_redundant_contigs"],
            )

    return ContigExtractionConfig(
        prefer_high_count_seeds=not args.disable_prefer_high_count_seeds,
        prefer_non_repeat_seeds=not args.disable_prefer_non_repeat_seeds,
        enable_repeat_seed_completion=not args.disable_repeat_seed_completion,
        suppress_redundant_contigs=args.suppress_redundant_contigs,
    )


def load_baseline_summary(summary_path: Path) -> dict[str, Any]:
    result = json.loads(summary_path.read_text(encoding="utf-8"))
    return {
        "config": result["config"],
        "train_metrics": result["train_metrics"],
        "val_metrics": result["val_metrics"],
        "heldout_metrics": result.get("heldout_metrics"),
        "quick_test": result.get("quick_test"),
        "repeat_heavy": result.get("repeat_heavy"),
        "summary_path": str(summary_path),
    }


def evaluate_panel(binary: Path, task_root: Path, config: ContigExtractionConfig) -> dict[str, Any]:
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


def parse_first_int(pattern: str, text: str) -> int | None:
    match = re.search(pattern, text, re.MULTILINE)
    if match is None:
        return None
    return int(match.group(1))


def delta(candidate: float | int | None, baseline: float | int | None) -> float | None:
    if candidate is None or baseline is None:
        return None
    return float(candidate) - float(baseline)


def ratio_change(candidate: float | int | None, baseline: float | int | None) -> float | None:
    if candidate is None or baseline in (None, 0):
        return None
    return (float(candidate) - float(baseline)) / float(baseline)


def read_fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current_name: str | None = None
    current_len = 0
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    lengths[current_name] = current_len
                current_name = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line)
    if current_name is not None:
        lengths[current_name] = current_len
    return lengths


def merged_interval_bp(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    current_start, current_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
            continue
        total += current_end - current_start
        current_start, current_end = start, end
    total += current_end - current_start
    return total


def run_reference_eval(
    reference_path: Path,
    assembly_path: Path,
    output_dir: Path,
    minimap2_path: str | None,
    threads: int,
    min_mapq: int,
    min_aligned_bp: int,
) -> dict[str, Any]:
    if minimap2_path is None:
        return {"skipped": "minimap2 not found"}
    if not reference_path.exists():
        return {"skipped": f"reference not found at {reference_path}"}
    if not assembly_path.exists():
        return {"skipped": f"assembly not found at {assembly_path}"}

    output_dir.mkdir(parents=True, exist_ok=True)
    paf_path = output_dir / "alignments.paf"
    stderr_path = output_dir / "minimap2.stderr.log"
    cmd = [
        minimap2_path,
        "-x",
        "asm5",
        "-t",
        str(max(1, threads)),
        str(reference_path),
        str(assembly_path),
    ]
    completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
    paf_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")

    reference_lengths = read_fasta_lengths(reference_path)
    assembly_lengths = read_fasta_lengths(assembly_path)
    reference_intervals: dict[str, list[tuple[int, int]]] = {}
    query_intervals: dict[str, list[tuple[int, int]]] = {}
    query_records: dict[str, int] = {}
    matched_bases = 0
    aligned_bases = 0
    mapq_sum = 0
    kept_records = 0

    for line in completed.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        qname = fields[0]
        qstart = int(fields[2])
        qend = int(fields[3])
        tname = fields[5]
        tstart = int(fields[7])
        tend = int(fields[8])
        nmatch = int(fields[9])
        block_len = int(fields[10])
        mapq = int(fields[11])
        if mapq < min_mapq or (qend - qstart) < min_aligned_bp:
            continue
        kept_records += 1
        matched_bases += nmatch
        aligned_bases += block_len
        mapq_sum += mapq
        reference_intervals.setdefault(tname, []).append((tstart, tend))
        query_intervals.setdefault(qname, []).append((qstart, qend))
        query_records[qname] = query_records.get(qname, 0) + 1

    reference_total_bp = sum(reference_lengths.values())
    assembly_total_bp = sum(assembly_lengths.values())
    reference_covered_bp = sum(
        merged_interval_bp(intervals) for intervals in reference_intervals.values()
    )
    aligned_query_bp = sum(merged_interval_bp(intervals) for intervals in query_intervals.values())
    queries_with_alignment = len(query_intervals)
    queries_with_multiple_alignments = sum(
        1 for records in query_records.values() if records > 1
    )

    return {
        "tool": "minimap2",
        "reference_path": str(reference_path),
        "assembly_path": str(assembly_path),
        "paf_path": str(paf_path),
        "stderr_log": str(stderr_path),
        "reference_total_bp": reference_total_bp,
        "assembly_total_bp": assembly_total_bp,
        "reference_covered_bp": reference_covered_bp,
        "reference_coverage_fraction": (
            reference_covered_bp / reference_total_bp if reference_total_bp else None
        ),
        "aligned_query_bp": aligned_query_bp,
        "aligned_query_fraction": (
            aligned_query_bp / assembly_total_bp if assembly_total_bp else None
        ),
        "aligned_records": kept_records,
        "queries_with_alignment": queries_with_alignment,
        "queries_with_multiple_alignments": queries_with_multiple_alignments,
        "multi_alignment_query_fraction": (
            queries_with_multiple_alignments / queries_with_alignment
            if queries_with_alignment
            else None
        ),
        "mean_mapq": (mapq_sum / kept_records if kept_records else None),
        "matched_base_fraction": (matched_bases / aligned_bases if aligned_bases else None),
    }


def build_assembly_command(
    binary: Path,
    dataset: DatasetSpec,
    config: ContigExtractionConfig,
    threads: int,
    output_dir: Path,
) -> tuple[list[str], Path, Path, Path]:
    contigs_path = output_dir / "contigs.fa"
    scaffold_path = output_dir / "contigs.scaffolds.fa"
    polished_path = output_dir / "contigs.polished.fa"
    cmd = [
        str(binary),
        "assemble-large",
        "-i",
        str(dataset.reads_1),
        "--input2",
        str(dataset.reads_2),
        "-o",
        str(contigs_path),
        "-t",
        str(threads),
        "--min-count",
        "0",
        "--scaffold",
        "--polish",
        "--polish-iterations",
        "1",
        "--compress-buckets",
    ]
    if not config.prefer_high_count_seeds:
        cmd.append("--disable-prefer-high-count-seeds")
    if not config.prefer_non_repeat_seeds:
        cmd.append("--disable-prefer-non-repeat-seeds")
    if not config.enable_repeat_seed_completion:
        cmd.append("--disable-repeat-seed-completion")
    if config.suppress_redundant_contigs:
        cmd.append("--suppress-redundant-contigs")
    return cmd, contigs_path, scaffold_path, polished_path


def run_end_to_end_dataset(
    binary: Path,
    dataset: DatasetSpec,
    config: ContigExtractionConfig,
    threads: int,
    output_dir: Path,
    minimap2_path: str | None,
    reference_eval_min_mapq: int,
    reference_eval_min_aligned_bp: int,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd, contigs_path, scaffold_path, polished_path = build_assembly_command(
        binary, dataset, config, threads, output_dir
    )

    start = time.time()
    completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
    elapsed = time.time() - start

    stdout_path = output_dir / "stdout.log"
    stderr_path = output_dir / "stderr.log"
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")

    combined = completed.stdout + "\n" + completed.stderr
    if polished_path.exists():
        reference_eval_input = polished_path
        reference_eval_target = "polished"
    elif scaffold_path.exists():
        reference_eval_input = scaffold_path
        reference_eval_target = "scaffold"
    else:
        reference_eval_input = contigs_path
        reference_eval_target = "contigs"

    result = {
        "dataset": dataset.name,
        "command": cmd,
        "elapsed_seconds": elapsed,
        "contigs": parse_first_int(r"^Contigs:\s+(\d+)$", combined),
        "total_length_bp": parse_first_int(r"^Total length:\s+(\d+)\s+bp$", combined),
        "contig_n50_bp": parse_first_int(r"^N50:\s+(\d+)\s+bp$", combined),
        "scaffolds": parse_first_int(r"^\s+Scaffolds:\s+(\d+)$", combined),
        "scaffold_n50_bp": parse_first_int(r"^\s+N50:\s+(\d+)\s+bp$", combined),
        "polish_corrections": parse_first_int(r"^\s+Corrections:\s+(\d+)$", combined),
        "stdout_log": str(stdout_path),
        "stderr_log": str(stderr_path),
        "contigs_path": str(contigs_path),
        "scaffold_path": str(scaffold_path),
        "polished_path": str(polished_path),
        "reference_eval_target": reference_eval_target,
    }
    result["reference_eval"] = run_reference_eval(
        dataset.reference,
        reference_eval_input,
        output_dir / "reference_eval",
        minimap2_path,
        threads,
        reference_eval_min_mapq,
        reference_eval_min_aligned_bp,
    )
    return result


def compare_dataset_results(
    prefix: str,
    candidate_result: dict[str, Any] | None,
    baseline_result: dict[str, Any] | None,
) -> dict[str, Any]:
    candidate_result = candidate_result or {}
    baseline_result = baseline_result or {}
    candidate_reference = candidate_result.get("reference_eval") or {}
    baseline_reference = baseline_result.get("reference_eval") or {}
    return {
        f"{prefix}_elapsed_seconds_delta": delta(
            candidate_result.get("elapsed_seconds"), baseline_result.get("elapsed_seconds")
        ),
        f"{prefix}_elapsed_ratio_change": ratio_change(
            candidate_result.get("elapsed_seconds"), baseline_result.get("elapsed_seconds")
        ),
        f"{prefix}_contigs_delta": delta(
            candidate_result.get("contigs"), baseline_result.get("contigs")
        ),
        f"{prefix}_total_length_delta": delta(
            candidate_result.get("total_length_bp"), baseline_result.get("total_length_bp")
        ),
        f"{prefix}_contig_n50_delta": delta(
            candidate_result.get("contig_n50_bp"), baseline_result.get("contig_n50_bp")
        ),
        f"{prefix}_contig_n50_ratio_change": ratio_change(
            candidate_result.get("contig_n50_bp"), baseline_result.get("contig_n50_bp")
        ),
        f"{prefix}_scaffold_n50_delta": delta(
            candidate_result.get("scaffold_n50_bp"), baseline_result.get("scaffold_n50_bp")
        ),
        f"{prefix}_scaffold_n50_ratio_change": ratio_change(
            candidate_result.get("scaffold_n50_bp"), baseline_result.get("scaffold_n50_bp")
        ),
        f"{prefix}_reference_coverage_fraction_delta": delta(
            candidate_reference.get("reference_coverage_fraction"),
            baseline_reference.get("reference_coverage_fraction"),
        ),
        f"{prefix}_aligned_query_fraction_delta": delta(
            candidate_reference.get("aligned_query_fraction"),
            baseline_reference.get("aligned_query_fraction"),
        ),
        f"{prefix}_multi_alignment_query_fraction_delta": delta(
            candidate_reference.get("multi_alignment_query_fraction"),
            baseline_reference.get("multi_alignment_query_fraction"),
        ),
        f"{prefix}_matched_base_fraction_delta": delta(
            candidate_reference.get("matched_base_fraction"),
            baseline_reference.get("matched_base_fraction"),
        ),
    }


def compare_candidate_against_baseline(
    candidate: dict[str, Any], baseline: dict[str, Any]
) -> dict[str, Any]:
    candidate_heldout = candidate.get("heldout_metrics") or {}
    baseline_heldout = baseline.get("heldout_metrics") or {}
    comparison = {
        "train_score_delta": delta(
            candidate["train_metrics"]["score"], baseline["train_metrics"]["score"]
        ),
        "val_score_delta": delta(
            candidate["val_metrics"]["score"], baseline["val_metrics"]["score"]
        ),
        "heldout_score_delta": delta(
            candidate_heldout.get("score"), baseline_heldout.get("score")
        ),
        "heldout_exact_contig_rate_delta": delta(
            candidate_heldout.get("exact_contig_rate"),
            baseline_heldout.get("exact_contig_rate"),
        ),
        "heldout_truth_kmer_f1_delta": delta(
            candidate_heldout.get("truth_kmer_f1"),
            baseline_heldout.get("truth_kmer_f1"),
        ),
        "heldout_count_agreement_delta": delta(
            candidate_heldout.get("contig_count_agreement"),
            baseline_heldout.get("contig_count_agreement"),
        ),
    }
    comparison.update(
        compare_dataset_results(
            "quick_test", candidate.get("quick_test"), baseline.get("quick_test")
        )
    )
    comparison.update(
        compare_dataset_results(
            "repeat_heavy", candidate.get("repeat_heavy"), baseline.get("repeat_heavy")
        )
    )
    return comparison


def decide_promotion(
    candidate: dict[str, Any],
    baseline: dict[str, Any],
    comparison: dict[str, Any],
    gate: PromotionGate,
) -> dict[str, Any]:
    hard_fail_reasons: list[str] = []
    review_reasons: list[str] = []

    if (
        comparison["val_score_delta"] is not None
        and comparison["val_score_delta"] < gate.min_val_score_delta
    ):
        hard_fail_reasons.append("validation score did not beat the baseline")

    heldout_score_delta = comparison["heldout_score_delta"]
    if heldout_score_delta is None:
        hard_fail_reasons.append("heldout panel evaluation missing for contig gate")
    elif heldout_score_delta < gate.min_heldout_score_delta:
        hard_fail_reasons.append("heldout score regressed versus baseline")

    if (
        comparison["heldout_exact_contig_rate_delta"] is not None
        and comparison["heldout_exact_contig_rate_delta"]
        < -gate.max_heldout_exact_contig_rate_drop
    ):
        hard_fail_reasons.append("heldout exact contig rate regressed beyond the allowed threshold")
    if (
        comparison["heldout_truth_kmer_f1_delta"] is not None
        and comparison["heldout_truth_kmer_f1_delta"] < -gate.max_heldout_truth_kmer_f1_drop
    ):
        hard_fail_reasons.append("heldout truth-kmer F1 regressed beyond the allowed threshold")
    if (
        comparison["heldout_count_agreement_delta"] is not None
        and comparison["heldout_count_agreement_delta"] < -gate.max_heldout_count_agreement_drop
    ):
        hard_fail_reasons.append("heldout contig-count agreement regressed beyond the allowed threshold")

    def check_suite_quality(prefix: str, contig_n50_drop_ratio: float, scaffold_n50_drop_ratio: float,
                            reference_coverage_drop: float, aligned_query_fraction_drop: float) -> None:
        contig_n50_ratio = comparison[f"{prefix}_contig_n50_ratio_change"]
        if contig_n50_ratio is not None and contig_n50_ratio < -contig_n50_drop_ratio:
            hard_fail_reasons.append(
                f"{prefix} contig N50 regressed beyond the allowed threshold"
            )

        scaffold_n50_ratio = comparison[f"{prefix}_scaffold_n50_ratio_change"]
        if scaffold_n50_ratio is not None and scaffold_n50_ratio < -scaffold_n50_drop_ratio:
            hard_fail_reasons.append(
                f"{prefix} scaffold N50 regressed beyond the allowed threshold"
            )

        reference_coverage_delta = comparison[f"{prefix}_reference_coverage_fraction_delta"]
        if (
            reference_coverage_delta is not None
            and reference_coverage_delta < -reference_coverage_drop
        ):
            hard_fail_reasons.append(
                f"{prefix} reference coverage regressed beyond the allowed threshold"
            )

        aligned_query_delta = comparison[f"{prefix}_aligned_query_fraction_delta"]
        if (
            aligned_query_delta is not None
            and aligned_query_delta < -aligned_query_fraction_drop
        ):
            hard_fail_reasons.append(
                f"{prefix} aligned assembly fraction regressed beyond the allowed threshold"
            )

    check_suite_quality(
        "quick_test",
        gate.max_quick_test_contig_n50_drop_ratio,
        gate.max_quick_test_scaffold_n50_drop_ratio,
        gate.max_quick_test_reference_coverage_drop,
        gate.max_quick_test_aligned_query_fraction_drop,
    )
    check_suite_quality(
        "repeat_heavy",
        gate.max_repeat_heavy_contig_n50_drop_ratio,
        gate.max_repeat_heavy_scaffold_n50_drop_ratio,
        gate.max_repeat_heavy_reference_coverage_drop,
        gate.max_repeat_heavy_aligned_query_fraction_drop,
    )

    def check_runtime(prefix: str, max_ratio: float, max_seconds: float) -> None:
        runtime_ratio = comparison[f"{prefix}_elapsed_ratio_change"]
        runtime_delta = comparison[f"{prefix}_elapsed_seconds_delta"]
        if runtime_ratio is not None and runtime_delta is not None:
            if runtime_ratio > max_ratio and runtime_delta > max_seconds:
                review_reasons.append(
                    f"{prefix} runtime regressed beyond the allowed threshold"
                )

    check_runtime(
        "quick_test",
        gate.max_quick_test_runtime_regression_ratio,
        gate.max_quick_test_runtime_regression_seconds,
    )
    check_runtime(
        "repeat_heavy",
        gate.max_repeat_heavy_runtime_regression_ratio,
        gate.max_repeat_heavy_runtime_regression_seconds,
    )

    if candidate.get("repeat_heavy") is None or baseline.get("repeat_heavy") is None:
        review_reasons.append("repeat-heavy end-to-end suite missing")
    else:
        contig_reduction = comparison["repeat_heavy_contigs_delta"]
        multi_alignment_delta = comparison["repeat_heavy_multi_alignment_query_fraction_delta"]
        has_contig_reduction = (
            contig_reduction is not None
            and contig_reduction <= -gate.min_repeat_heavy_contig_reduction
        )
        has_multi_alignment_improvement = (
            multi_alignment_delta is not None
            and multi_alignment_delta <= -gate.min_repeat_heavy_multi_alignment_improvement
        )
        if not has_contig_reduction and not has_multi_alignment_improvement:
            review_reasons.append(
                "repeat-heavy suite did not show a meaningful redundant-contig reduction"
            )
        elif multi_alignment_delta is not None and multi_alignment_delta > 0.0:
            review_reasons.append(
                "repeat-heavy reference alignment fragmentation increased versus baseline"
            )

    if hard_fail_reasons:
        status = "reject"
    elif review_reasons:
        status = "review"
    else:
        status = "promote_default"

    candidate_repeat = candidate.get("repeat_heavy") or {}
    baseline_repeat = baseline.get("repeat_heavy") or {}
    return {
        "status": status,
        "promote_default": status == "promote_default",
        "hard_fail_reasons": hard_fail_reasons,
        "review_reasons": review_reasons,
        "gate": asdict(gate),
        "candidate_summary": {
            "config": candidate["config"],
            "heldout_score": (candidate.get("heldout_metrics") or {}).get("score"),
            "heldout_exact_contig_rate": (candidate.get("heldout_metrics") or {}).get("exact_contig_rate"),
            "heldout_truth_kmer_f1": (candidate.get("heldout_metrics") or {}).get("truth_kmer_f1"),
            "quick_test_elapsed_seconds": (candidate.get("quick_test") or {}).get("elapsed_seconds"),
            "repeat_heavy_contigs": candidate_repeat.get("contigs"),
            "repeat_heavy_multi_alignment_query_fraction": (
                (candidate_repeat.get("reference_eval") or {}).get("multi_alignment_query_fraction")
            ),
        },
        "baseline_summary": {
            "config": baseline["config"],
            "heldout_score": (baseline.get("heldout_metrics") or {}).get("score"),
            "heldout_exact_contig_rate": (baseline.get("heldout_metrics") or {}).get("exact_contig_rate"),
            "heldout_truth_kmer_f1": (baseline.get("heldout_metrics") or {}).get("truth_kmer_f1"),
            "quick_test_elapsed_seconds": (baseline.get("quick_test") or {}).get("elapsed_seconds"),
            "repeat_heavy_contigs": baseline_repeat.get("contigs"),
            "repeat_heavy_multi_alignment_query_fraction": (
                (baseline_repeat.get("reference_eval") or {}).get("multi_alignment_query_fraction")
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Promote a contig-extraction candidate through heldout eval and repeat-heavy suites"
    )
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--tasks-root", type=Path, default=DEFAULT_TASKS_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--train-result", type=Path, default=DEFAULT_TRAIN_RESULT)
    parser.add_argument("--disable-prefer-high-count-seeds", action="store_true")
    parser.add_argument("--disable-prefer-non-repeat-seeds", action="store_true")
    parser.add_argument("--disable-repeat-seed-completion", action="store_true")
    parser.add_argument("--suppress-redundant-contigs", action="store_true")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--skip-quick-test", action="store_true")
    parser.add_argument("--skip-repeat-heavy", action="store_true")
    parser.add_argument("--baseline-summary", type=Path)
    parser.add_argument("--skip-baseline-end-to-end", action="store_true")
    parser.add_argument("--repeat-heavy-max-pairs", type=int, default=DEFAULT_REPEAT_HEAVY_MAX_PAIRS)
    parser.add_argument("--min-heldout-score-delta", type=float, default=0.0)
    parser.add_argument("--max-heldout-exact-contig-rate-drop", type=float, default=0.0)
    parser.add_argument("--max-heldout-truth-kmer-f1-drop", type=float, default=0.0)
    parser.add_argument("--max-heldout-count-agreement-drop", type=float, default=0.0)
    parser.add_argument("--max-quick-test-runtime-regression-ratio", type=float, default=0.05)
    parser.add_argument("--max-quick-test-runtime-regression-seconds", type=float, default=5.0)
    parser.add_argument("--max-quick-test-contig-n50-drop-ratio", type=float, default=0.0)
    parser.add_argument("--max-quick-test-scaffold-n50-drop-ratio", type=float, default=0.0)
    parser.add_argument("--max-quick-test-reference-coverage-drop", type=float, default=0.002)
    parser.add_argument("--max-quick-test-aligned-query-fraction-drop", type=float, default=0.01)
    parser.add_argument("--max-repeat-heavy-runtime-regression-ratio", type=float, default=0.10)
    parser.add_argument("--max-repeat-heavy-runtime-regression-seconds", type=float, default=30.0)
    parser.add_argument("--max-repeat-heavy-contig-n50-drop-ratio", type=float, default=0.05)
    parser.add_argument("--max-repeat-heavy-scaffold-n50-drop-ratio", type=float, default=0.05)
    parser.add_argument("--max-repeat-heavy-reference-coverage-drop", type=float, default=0.002)
    parser.add_argument("--max-repeat-heavy-aligned-query-fraction-drop", type=float, default=0.01)
    parser.add_argument("--min-repeat-heavy-contig-reduction", type=int, default=1)
    parser.add_argument("--min-repeat-heavy-multi-alignment-improvement", type=float, default=0.01)
    parser.add_argument("--minimap2", type=str, default=shutil.which("minimap2"))
    parser.add_argument("--reference-eval-min-mapq", type=int, default=20)
    parser.add_argument("--reference-eval-min-aligned-bp", type=int, default=500)
    args = parser.parse_args()

    if not args.binary.exists():
        raise FileNotFoundError(
            f"missing raptor binary at {args.binary}; build it first with cargo build --release"
        )

    train_root = args.tasks_root / "train"
    val_root = args.tasks_root / "val"
    heldout_root = args.tasks_root / "heldout"
    if not train_root.exists() or not val_root.exists() or not heldout_root.exists():
        raise FileNotFoundError(
            f"expected train/val/heldout panels under {args.tasks_root}; run contig_extraction_prepare.py first"
        )

    config = load_config(args)
    gate = PromotionGate(
        min_heldout_score_delta=args.min_heldout_score_delta,
        max_heldout_exact_contig_rate_drop=args.max_heldout_exact_contig_rate_drop,
        max_heldout_truth_kmer_f1_drop=args.max_heldout_truth_kmer_f1_drop,
        max_heldout_count_agreement_drop=args.max_heldout_count_agreement_drop,
        max_quick_test_runtime_regression_ratio=args.max_quick_test_runtime_regression_ratio,
        max_quick_test_runtime_regression_seconds=args.max_quick_test_runtime_regression_seconds,
        max_quick_test_contig_n50_drop_ratio=args.max_quick_test_contig_n50_drop_ratio,
        max_quick_test_scaffold_n50_drop_ratio=args.max_quick_test_scaffold_n50_drop_ratio,
        max_quick_test_reference_coverage_drop=args.max_quick_test_reference_coverage_drop,
        max_quick_test_aligned_query_fraction_drop=args.max_quick_test_aligned_query_fraction_drop,
        max_repeat_heavy_runtime_regression_ratio=args.max_repeat_heavy_runtime_regression_ratio,
        max_repeat_heavy_runtime_regression_seconds=args.max_repeat_heavy_runtime_regression_seconds,
        max_repeat_heavy_contig_n50_drop_ratio=args.max_repeat_heavy_contig_n50_drop_ratio,
        max_repeat_heavy_scaffold_n50_drop_ratio=args.max_repeat_heavy_scaffold_n50_drop_ratio,
        max_repeat_heavy_reference_coverage_drop=args.max_repeat_heavy_reference_coverage_drop,
        max_repeat_heavy_aligned_query_fraction_drop=args.max_repeat_heavy_aligned_query_fraction_drop,
        min_repeat_heavy_contig_reduction=args.min_repeat_heavy_contig_reduction,
        min_repeat_heavy_multi_alignment_improvement=args.min_repeat_heavy_multi_alignment_improvement,
    )
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = args.output_root / f"contig_extraction_{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=True)

    quick_test_dataset = DatasetSpec(
        name="quick_test",
        reads_1=QUICK_TEST_READS_1,
        reads_2=QUICK_TEST_READS_2,
        reference=QUICK_TEST_REFERENCE,
    )
    if args.repeat_heavy_max_pairs > 0:
        repeat_heavy_reads_1, repeat_heavy_reads_2 = ensure_repeat_heavy_subset(
            args.repeat_heavy_max_pairs
        )
        repeat_heavy_name = f"drosophila_subset_{args.repeat_heavy_max_pairs:07d}_pairs"
    else:
        repeat_heavy_reads_1, repeat_heavy_reads_2 = (
            REPEAT_HEAVY_SOURCE_READS_1,
            REPEAT_HEAVY_SOURCE_READS_2,
        )
        repeat_heavy_name = "drosophila_full"
    repeat_heavy_dataset = DatasetSpec(
        name=repeat_heavy_name,
        reads_1=repeat_heavy_reads_1,
        reads_2=repeat_heavy_reads_2,
        reference=REPEAT_HEAVY_REFERENCE,
    )

    print(f"Binary: {args.binary}")
    print(f"Tasks root: {args.tasks_root}")
    print(f"Promotion output: {output_dir}")
    print(f"Config: {asdict(config)}")

    candidate_result = {
        "config": asdict(config),
        "train_metrics": evaluate_panel(args.binary, train_root, config),
        "val_metrics": evaluate_panel(args.binary, val_root, config),
        "heldout_metrics": evaluate_panel(args.binary, heldout_root, config),
    }

    if args.baseline_summary is not None and args.baseline_summary.exists():
        baseline_result = load_baseline_summary(args.baseline_summary)
        baseline_config = ContigExtractionConfig(**baseline_result["config"])
    else:
        baseline_config = ContigExtractionConfig(**DEFAULT_BASELINE_CONFIG)
        baseline_result = {
            "config": asdict(baseline_config),
            "train_metrics": evaluate_panel(args.binary, train_root, baseline_config),
            "val_metrics": evaluate_panel(args.binary, val_root, baseline_config),
            "heldout_metrics": evaluate_panel(args.binary, heldout_root, baseline_config),
            "quick_test": None,
            "repeat_heavy": None,
        }

    quick_test_available = (
        quick_test_dataset.reads_1.exists()
        and quick_test_dataset.reads_2.exists()
        and quick_test_dataset.reference.exists()
    )
    repeat_heavy_available = (
        repeat_heavy_dataset.reads_1.exists()
        and repeat_heavy_dataset.reads_2.exists()
        and repeat_heavy_dataset.reference.exists()
    )

    if args.skip_quick_test:
        candidate_result["quick_test"] = None
        if args.skip_baseline_end_to_end:
            baseline_result["quick_test"] = None
    elif quick_test_available:
        candidate_result["quick_test"] = run_end_to_end_dataset(
            args.binary,
            quick_test_dataset,
            config,
            args.threads,
            output_dir / "quick_test",
            args.minimap2,
            args.reference_eval_min_mapq,
            args.reference_eval_min_aligned_bp,
        )
        if not args.skip_baseline_end_to_end and baseline_result.get("quick_test") is None:
            baseline_result["quick_test"] = run_end_to_end_dataset(
                args.binary,
                quick_test_dataset,
                baseline_config,
                args.threads,
                output_dir / "baseline_quick_test",
                args.minimap2,
                args.reference_eval_min_mapq,
                args.reference_eval_min_aligned_bp,
            )
    else:
        candidate_result["quick_test"] = {"skipped": "quick_test data not found"}
        if baseline_result.get("quick_test") is None:
            baseline_result["quick_test"] = {"skipped": "quick_test data not found"}

    if args.skip_repeat_heavy:
        candidate_result["repeat_heavy"] = None
        if args.skip_baseline_end_to_end:
            baseline_result["repeat_heavy"] = None
    elif repeat_heavy_available:
        candidate_result["repeat_heavy"] = run_end_to_end_dataset(
            args.binary,
            repeat_heavy_dataset,
            config,
            args.threads,
            output_dir / "repeat_heavy",
            args.minimap2,
            args.reference_eval_min_mapq,
            args.reference_eval_min_aligned_bp,
        )
        if not args.skip_baseline_end_to_end and baseline_result.get("repeat_heavy") is None:
            baseline_result["repeat_heavy"] = run_end_to_end_dataset(
                args.binary,
                repeat_heavy_dataset,
                baseline_config,
                args.threads,
                output_dir / "baseline_repeat_heavy",
                args.minimap2,
                args.reference_eval_min_mapq,
                args.reference_eval_min_aligned_bp,
            )
    else:
        candidate_result["repeat_heavy"] = {"skipped": "repeat-heavy data not found"}
        if baseline_result.get("repeat_heavy") is None:
            baseline_result["repeat_heavy"] = {"skipped": "repeat-heavy data not found"}

    comparison = compare_candidate_against_baseline(candidate_result, baseline_result)
    decision = decide_promotion(candidate_result, baseline_result, comparison, gate)

    result = {
        "config": candidate_result["config"],
        "train_metrics": candidate_result["train_metrics"],
        "val_metrics": candidate_result["val_metrics"],
        "heldout_metrics": candidate_result["heldout_metrics"],
        "quick_test": candidate_result.get("quick_test"),
        "repeat_heavy": candidate_result.get("repeat_heavy"),
        "baseline": baseline_result,
        "comparison": comparison,
        "decision": decision,
        "gate": asdict(gate),
    }

    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    print(f"Wrote promotion summary to {summary_path}")


if __name__ == "__main__":
    main()
