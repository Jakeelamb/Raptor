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
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "target" / "trinity_parity" / "tiny_alt_isoform"
DEFAULT_ORACLE = ROOT / "bench" / "trinity_parity" / "oracles" / "tiny_alt_isoform.fa"
DEFAULT_FIXTURE = "tiny_alt_isoform"


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


def file_size_bytes(path: Path) -> int | None:
    if not path.exists() or not path.is_file():
        return None
    return path.stat().st_size


def directory_footprint(path: Path) -> dict[str, object]:
    file_count = 0
    total_bytes = 0
    if not path.exists():
        return {"exists": False, "file_count": 0, "total_bytes": 0}
    for child in path.rglob("*"):
        if child.is_file():
            file_count += 1
            total_bytes += child.stat().st_size
    return {"exists": True, "file_count": file_count, "total_bytes": total_bytes}


def component_evidence_metrics(path: Path) -> dict[str, object]:
    if not path.exists():
        return {"component_evidence_exists": False}
    components = json.loads(path.read_text(encoding="utf-8"))
    return {
        "component_evidence_exists": True,
        "component_assigned_read_count": sum(
            int(component.get("assigned_read_count", 0)) for component in components
        ),
        "component_assigned_pair_count": sum(
            int(component.get("assigned_pair_count", 0)) for component in components
        ),
    }


def component_graph_metrics(path: Path) -> dict[str, object]:
    if not path.exists():
        return {"component_graphs_exist": False}
    graphs = json.loads(path.read_text(encoding="utf-8"))
    edges = [
        edge
        for graph in graphs
        for edge in graph.get("edges", [])
        if isinstance(edge, dict)
    ]
    read_kmer_edge_samples = [
        edge
        for graph in graphs
        for edge in graph.get("read_kmer_edges_sample", [])
        if isinstance(edge, dict)
    ]
    read_kmer_nodes = [
        node
        for graph in graphs
        for node in graph.get("read_kmer_nodes", [])
        if isinstance(node, str)
    ]
    read_kmer_edges = [
        edge
        for graph in graphs
        for edge in graph.get("read_kmer_edges", [])
        if isinstance(edge, dict)
    ]
    read_kmer_paths = [
        path
        for graph in graphs
        for path in graph.get("read_kmer_paths", [])
        if isinstance(path, dict)
    ]
    read_kmer_node_count = sum(int(graph.get("read_kmer_node_count", 0)) for graph in graphs)
    read_kmer_edge_count = sum(int(graph.get("read_kmer_edge_count", 0)) for graph in graphs)
    read_kmer_path_count = sum(int(graph.get("read_kmer_path_count", 0)) for graph in graphs)
    return {
        "component_graphs_exist": True,
        "component_graph_count": len(graphs),
        "component_graph_node_count": sum(int(graph.get("node_count", 0)) for graph in graphs),
        "component_graph_edge_count": sum(int(graph.get("edge_count", 0)) for graph in graphs),
        "component_graph_edge_read_support": sum(
            int(edge.get("shared_read_count", 0)) for edge in edges
        ),
        "component_graph_edge_pair_support": sum(
            int(edge.get("shared_pair_count", 0)) for edge in edges
        ),
        "component_graph_edge_observed_kmers": sum(
            int(edge.get("shared_observed_kmer_count", 0)) for edge in edges
        ),
        "component_graph_read_kmer_node_count": read_kmer_node_count,
        "component_graph_read_kmer_edge_count": read_kmer_edge_count,
        "component_graph_read_kmer_node_record_count": len(read_kmer_nodes),
        "component_graph_read_kmer_edge_record_count": len(read_kmer_edges),
        "component_graph_read_kmer_edge_sample_count": len(read_kmer_edge_samples),
        "component_graph_read_kmer_path_count": read_kmer_path_count,
        "component_graph_read_kmer_path_record_count": len(read_kmer_paths),
        "component_graph_read_kmer_path_edge_count": sum(
            int(path.get("edge_count", 0)) for path in read_kmer_paths
        ),
        "component_graph_read_kmer_path_min_support": min(
            (int(path.get("min_support", 0)) for path in read_kmer_paths),
            default=0,
        ),
    }


def selected_isoform_evidence_metrics(path: Path) -> dict[str, object]:
    if not path.exists():
        return {"component_selected_isoforms_json_exists": False}
    records = json.loads(path.read_text(encoding="utf-8"))
    methods = sorted(
        {
            str(record.get("selection_method", ""))
            for record in records
            if record.get("selection_method")
        }
    )
    ranks_by_component: dict[int, list[int]] = {}
    scores_by_component: dict[int, list[int]] = {}
    for record in records:
        component_id = int(record.get("component_id", -1))
        ranks_by_component.setdefault(component_id, []).append(
            int(record.get("component_rank", 0))
        )
        scores_by_component.setdefault(component_id, []).append(
            int(record.get("evidence_score", 0))
        )
    ranks_are_dense = all(
        sorted(ranks) == list(range(1, len(ranks) + 1))
        for ranks in ranks_by_component.values()
    )
    scores_are_descending = all(
        scores == sorted(scores, reverse=True) for scores in scores_by_component.values()
    )
    return {
        "component_selected_isoforms_json_exists": True,
        "component_selected_isoform_record_count": len(records),
        "component_selected_selection_methods": methods,
        "component_selected_component_ranks_are_dense": ranks_are_dense,
        "component_selected_scores_are_descending": scores_are_descending,
        "component_selected_evidence_score": sum(
            int(record.get("evidence_score", 0)) for record in records
        ),
        "component_selected_max_evidence_score": max(
            (int(record.get("evidence_score", 0)) for record in records),
            default=0,
        ),
        "component_selected_direct_read_support": sum(
            int(record.get("direct_read_support", 0)) for record in records
        ),
        "component_selected_direct_pair_support": sum(
            int(record.get("direct_pair_support", 0)) for record in records
        ),
        "component_selected_read_kmer_path_support": sum(
            int(record.get("overlapping_read_kmer_path_count", 0)) for record in records
        ),
        "component_selected_max_read_kmer_path_support": max(
            (
                int(record.get("max_overlapping_read_kmer_path_support", 0))
                for record in records
            ),
            default=0,
        ),
    }


def isoform_candidate_metrics(path: Path) -> dict[str, object]:
    if not path.exists():
        return {"component_isoform_candidates_json_exists": False}
    records = json.loads(path.read_text(encoding="utf-8"))
    selected = [record for record in records if record.get("selected") is True]
    rejected = [record for record in records if record.get("selected") is not True]
    path_candidates = [
        record for record in records if record.get("source_kind") == "read_kmer_path"
    ]
    selected_contigs = [
        record
        for record in selected
        if record.get("source_kind") == "contig"
    ]
    ranks_by_component: dict[int, list[int]] = {}
    scores_by_component: dict[int, list[int]] = {}
    for record in records:
        component_id = int(record.get("component_id", -1))
        ranks_by_component.setdefault(component_id, []).append(
            int(record.get("candidate_rank", 0))
        )
        scores_by_component.setdefault(component_id, []).append(
            int(record.get("evidence_score", 0))
        )
    return {
        "component_isoform_candidates_json_exists": True,
        "component_isoform_candidate_count": len(records),
        "component_isoform_selected_candidate_count": len(selected),
        "component_isoform_rejected_candidate_count": len(rejected),
        "component_isoform_path_candidate_count": len(path_candidates),
        "component_isoform_rejected_path_candidate_count": sum(
            1 for record in path_candidates if record.get("selected") is not True
        ),
        "component_isoform_selected_contig_candidate_count": len(selected_contigs),
        "component_isoform_candidate_ranks_are_dense": all(
            sorted(ranks) == list(range(1, len(ranks) + 1))
            for ranks in ranks_by_component.values()
        ),
        "component_isoform_candidate_scores_are_descending": all(
            scores == sorted(scores, reverse=True)
            for scores in scores_by_component.values()
        ),
        "component_isoform_candidate_max_score": max(
            (int(record.get("evidence_score", 0)) for record in records),
            default=0,
        ),
        "component_isoform_rejected_candidate_max_score": max(
            (int(record.get("evidence_score", 0)) for record in rejected),
            default=0,
        ),
    }


def count_fastq_records(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    lines = 0
    with opener(path, "rt", encoding="utf-8") as handle:
        for lines, _ in enumerate(handle, start=1):
            pass
    if lines % 4 != 0:
        raise ValueError(f"FASTQ line count is not divisible by 4: {path}")
    return lines // 4


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


def component_candidate_metrics(
    path: Path,
    selected_path: Path | None,
    truth_fasta: Path,
    oracle_fasta: Path,
    min_match_coverage: float,
) -> dict[str, object]:
    if not path.exists():
        return {"component_paths_fasta_exists": False}
    records = read_fasta_records(path)
    metrics: dict[str, object] = {
        "component_paths_fasta_exists": True,
        "component_path_transcript_count": len(records),
        "component_path_total_bases": sum(len(seq) for seq in records.values()),
        "component_path_n50": n50([len(seq) for seq in records.values()]),
        "component_path_lengths": [len(seq) for seq in records.values()],
    }
    metrics["component_candidate_truth_recovery"] = truth_recovery_metrics(truth_fasta, path)
    if oracle_fasta.exists():
        metrics["component_candidate_oracle_recovery"] = fasta_recovery_metrics(
            oracle_fasta, path
        )
    if selected_path is not None:
        metrics["component_selected_isoforms_fasta"] = str(selected_path)
        metrics["component_selected_isoforms_fasta_exists"] = selected_path.exists()
    if selected_path is not None and selected_path.exists():
        selected = read_fasta_records(selected_path)
        metrics["component_selected_isoform_count"] = len(selected)
        metrics["component_selected_isoform_lengths"] = [len(seq) for seq in selected.values()]
        metrics["component_selected_truth_recovery"] = truth_recovery_metrics(
            truth_fasta, selected_path
        )
        metrics["component_selected_truth_precision"] = fasta_precision_recall_metrics(
            truth_fasta, selected_path, min_match_coverage
        )
        if oracle_fasta.exists():
            metrics["component_selected_oracle_recovery"] = fasta_recovery_metrics(
                oracle_fasta, selected_path
            )
    return metrics


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


def fasta_precision_recall_metrics(
    truth_fasta: Path,
    assembly_fasta: Path,
    min_match_coverage: float,
) -> dict[str, object]:
    truth = read_fasta_records(truth_fasta)
    assembled = read_fasta_records(assembly_fasta)
    candidates: list[tuple[float, str, str, int]] = []
    for truth_name, truth_seq in truth.items():
        for assembled_name, assembled_seq in assembled.items():
            forward = longest_common_substring_len(truth_seq, assembled_seq)
            reverse = longest_common_substring_len(truth_seq, revcomp(assembled_seq))
            matching_bases = max(forward, reverse)
            truth_coverage = matching_bases / len(truth_seq) if truth_seq else 0.0
            assembled_coverage = matching_bases / len(assembled_seq) if assembled_seq else 0.0
            score = min(truth_coverage, assembled_coverage)
            if score >= min_match_coverage:
                candidates.append((score, truth_name, assembled_name, matching_bases))

    candidates.sort(key=lambda item: (-item[0], item[1], item[2]))
    matched_truth: set[str] = set()
    matched_assembled: set[str] = set()
    matches: list[dict[str, object]] = []
    for score, truth_name, assembled_name, matching_bases in candidates:
        if truth_name in matched_truth or assembled_name in matched_assembled:
            continue
        matched_truth.add(truth_name)
        matched_assembled.add(assembled_name)
        matches.append(
            {
                "truth": truth_name,
                "assembled": assembled_name,
                "matching_bases": matching_bases,
                "match_coverage": round(score, 6),
            }
        )

    true_positive = len(matches)
    false_positive = len(assembled) - true_positive
    false_negative = len(truth) - true_positive
    precision = true_positive / len(assembled) if assembled else 0.0
    recall = true_positive / len(truth) if truth else 0.0
    f1 = (
        2 * precision * recall / (precision + recall)
        if precision + recall > 0.0
        else 0.0
    )
    return {
        "match_min_coverage": min_match_coverage,
        "truth_transcript_count": len(truth),
        "assembled_record_count": len(assembled),
        "true_positive": true_positive,
        "false_positive": false_positive,
        "false_negative": false_negative,
        "precision": round(precision, 6),
        "recall": round(recall, 6),
        "f1": round(f1, 6),
        "matches": matches,
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
    time_bin = Path("/usr/bin/time")
    if not time_bin.exists():
        return run_command_with_procfs_resource_poll(command, cwd)

    with tempfile.NamedTemporaryFile(prefix="raptor-time-", delete=False) as handle:
        time_output_path = Path(handle.name)
    measured_command = [
        str(time_bin),
        "-v",
        "-o",
        str(time_output_path),
        *command,
    ]

    start = time.monotonic()
    gpu_samples = [sample_gpu_usage()]
    completed = subprocess.run(
        measured_command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    gpu_samples.append(sample_gpu_usage())
    elapsed = time.monotonic() - start
    result = {
        "command": command,
        "exit_code": completed.returncode,
        "elapsed_seconds": round(elapsed, 6),
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "resource_usage": command_resource_usage(time_output_path),
        "gpu_usage": summarize_gpu_samples(gpu_samples),
    }
    time_output_path.unlink(missing_ok=True)
    return result


def run_command_with_procfs_resource_poll(
    command: list[str], cwd: Path
) -> dict[str, object]:
    with tempfile.NamedTemporaryFile(prefix="raptor-stdout-", delete=False) as stdout_handle:
        stdout_path = Path(stdout_handle.name)
    with tempfile.NamedTemporaryFile(prefix="raptor-stderr-", delete=False) as stderr_handle:
        stderr_path = Path(stderr_handle.name)

    start = time.monotonic()
    max_rss_kb = 0
    gpu_samples = [sample_gpu_usage()]
    with stdout_path.open("w", encoding="utf-8") as stdout_file, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr_file:
        process = subprocess.Popen(
            command,
            cwd=cwd,
            text=True,
            stdout=stdout_file,
            stderr=stderr_file,
            start_new_session=True,
        )
        next_gpu_sample = time.monotonic() + 0.2
        while process.poll() is None:
            max_rss_kb = max(max_rss_kb, process_group_rss_kb(process.pid))
            now = time.monotonic()
            if now >= next_gpu_sample:
                gpu_samples.append(sample_gpu_usage())
                next_gpu_sample = now + 0.2
            time.sleep(0.01)
        max_rss_kb = max(max_rss_kb, process_group_rss_kb(process.pid))
        gpu_samples.append(sample_gpu_usage())

    elapsed = time.monotonic() - start
    stdout = stdout_path.read_text(encoding="utf-8", errors="replace")
    stderr = stderr_path.read_text(encoding="utf-8", errors="replace")
    stdout_path.unlink(missing_ok=True)
    stderr_path.unlink(missing_ok=True)
    return {
        "command": command,
        "exit_code": process.returncode,
        "elapsed_seconds": round(elapsed, 6),
        "stdout": stdout,
        "stderr": stderr,
        "resource_usage": {
            "available": max_rss_kb > 0,
            "source": "procfs_process_group_poll",
            "max_rss_kb": max_rss_kb if max_rss_kb > 0 else None,
            "user_seconds": None,
            "system_seconds": None,
        },
        "gpu_usage": summarize_gpu_samples(gpu_samples),
    }


def process_group_rss_kb(process_group_id: int) -> int:
    total = 0
    proc_root = Path("/proc")
    for status_path in proc_root.glob("[0-9]*/status"):
        try:
            status = status_path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        current_group = None
        current_rss = 0
        for line in status.splitlines():
            if line.startswith("NSpgid:"):
                fields = line.split()
                if fields:
                    current_group = int(fields[-1])
            elif line.startswith("VmRSS:"):
                fields = line.split()
                if len(fields) >= 2:
                    current_rss = int(fields[1])
        if current_group == process_group_id:
            total += current_rss
    return total


def command_resource_usage(path: Path | None) -> dict[str, object]:
    if path is None or not path.exists():
        return {
            "available": False,
            "source": None,
            "max_rss_kb": None,
            "user_seconds": None,
            "system_seconds": None,
        }

    values: dict[str, object] = {
        "available": True,
        "source": "gnu_time",
        "max_rss_kb": None,
        "user_seconds": None,
        "system_seconds": None,
    }
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if ":" not in line:
            continue
        key, raw_value = line.split(":", 1)
        value = raw_value.strip()
        if key == "Maximum resident set size (kbytes)":
            values["max_rss_kb"] = int(value)
        elif key == "User time (seconds)":
            values["user_seconds"] = float(value)
        elif key == "System time (seconds)":
            values["system_seconds"] = float(value)
    return values


def sample_gpu_usage() -> dict[str, object]:
    nvidia_smi = shutil.which("nvidia-smi")
    if not nvidia_smi:
        return {"available": False, "source": None, "gpus": []}
    command = [
        nvidia_smi,
        "--query-gpu=index,name,memory.used,utilization.gpu,power.draw",
        "--format=csv,noheader,nounits",
    ]
    completed = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        return {
            "available": False,
            "source": "nvidia-smi",
            "error": completed.stderr.strip(),
            "gpus": [],
        }

    gpus: list[dict[str, object]] = []
    for line in completed.stdout.splitlines():
        parts = [part.strip() for part in line.split(",")]
        if len(parts) < 5:
            continue
        gpus.append(
            {
                "index": int(parts[0]),
                "name": parts[1],
                "memory_used_mib": parse_optional_int(parts[2]),
                "utilization_gpu_percent": parse_optional_int(parts[3]),
                "power_draw_watts": parse_optional_float(parts[4]),
            }
        )
    return {"available": bool(gpus), "source": "nvidia-smi", "gpus": gpus}


def summarize_gpu_samples(samples: list[dict[str, object]]) -> dict[str, object]:
    available_samples = [sample for sample in samples if sample.get("available")]
    if not available_samples:
        return {"available": False, "source": None, "sample_count": len(samples), "gpus": []}

    by_index: dict[int, dict[str, object]] = {}
    for sample in available_samples:
        for gpu in sample.get("gpus", []):
            if not isinstance(gpu, dict):
                continue
            index = int(gpu["index"])
            summary = by_index.setdefault(
                index,
                {
                    "index": index,
                    "name": gpu.get("name"),
                    "max_memory_used_mib": 0,
                    "max_utilization_gpu_percent": 0,
                    "max_power_draw_watts": 0.0,
                },
            )
            memory = gpu.get("memory_used_mib")
            utilization = gpu.get("utilization_gpu_percent")
            power = gpu.get("power_draw_watts")
            if isinstance(memory, int):
                summary["max_memory_used_mib"] = max(
                    int(summary["max_memory_used_mib"]), memory
                )
            if isinstance(utilization, int):
                summary["max_utilization_gpu_percent"] = max(
                    int(summary["max_utilization_gpu_percent"]), utilization
                )
            if isinstance(power, float):
                summary["max_power_draw_watts"] = max(
                    float(summary["max_power_draw_watts"]), power
                )

    return {
        "available": True,
        "source": "nvidia-smi",
        "sample_count": len(samples),
        "gpus": [by_index[index] for index in sorted(by_index)],
    }


def parse_optional_int(value: str) -> int | None:
    try:
        return int(value)
    except ValueError:
        return None


def parse_optional_float(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None


def first_gpu_metric(gpu_usage: dict[str, object], key: str) -> object:
    gpus = gpu_usage.get("gpus", [])
    if not isinstance(gpus, list) or not gpus:
        return None
    first = gpus[0]
    if not isinstance(first, dict):
        return None
    return first.get(key)


def parse_insert_sweep(text: str) -> list[int]:
    inserts: list[int] = []
    for part in text.split(","):
        value = part.strip()
        if not value:
            continue
        insert = int(value)
        if insert <= 0:
            raise ValueError(f"insert sizes must be positive, got {insert}")
        inserts.append(insert)
    if not inserts:
        raise ValueError("insert sweep must contain at least one insert size")
    return inserts


def write_fixture_files(
    out_dir: Path,
    fixture_name: str,
    transcripts: dict[str, str],
    insert: int,
    coverage_rounds_by_transcript: dict[str, int],
) -> dict[str, object]:
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
    step = 24

    for tx_name, seq in transcripts.items():
        coverage_rounds = coverage_rounds_by_transcript[tx_name]
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
        "fixture": fixture_name,
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


def tiny_alt_isoform_transcripts() -> tuple[dict[str, str], dict[str, int]]:
    exon_a = deterministic_dna("shared_exon_a", 90)
    exon_b = deterministic_dna("dominant_exon_b", 72)
    exon_alt = deterministic_dna("alternative_exon", 60)
    exon_c = deterministic_dna("shared_exon_c", 90)

    tx1 = exon_a + exon_b + exon_c
    tx2 = exon_a + exon_alt + exon_c
    return {"tx_dominant": tx1, "tx_alt": tx2}, {"tx_dominant": 3, "tx_alt": 2}


def ambiguous_paralog_transcripts() -> tuple[dict[str, str], dict[str, int]]:
    shared_a = deterministic_dna("ambiguous_shared_exon_a", 84)
    shared_c = deterministic_dna("ambiguous_shared_exon_c", 84)
    dominant_mid = deterministic_dna("ambiguous_dominant_mid", 72)
    alt_mid = deterministic_dna("ambiguous_alt_mid", 66)
    paralog_a = deterministic_dna("ambiguous_paralog_a", 84)
    paralog_c = deterministic_dna("ambiguous_paralog_c", 84)

    transcripts = {
        "tx_major": shared_a + dominant_mid + shared_c,
        "tx_alt": shared_a + alt_mid + shared_c,
        "tx_paralog": paralog_a + dominant_mid + paralog_c,
    }
    coverage = {"tx_major": 4, "tx_alt": 2, "tx_paralog": 2}
    return transcripts, coverage


def compact_fusion_transcripts() -> tuple[dict[str, str], dict[str, int]]:
    left_unique = deterministic_dna("fusion_left_unique", 120)
    shared_overlap = deterministic_dna("fusion_shared_overlap", 96)
    right_unique = deterministic_dna("fusion_right_unique", 120)

    transcripts = {
        "tx_left": left_unique + shared_overlap,
        "tx_right": shared_overlap + right_unique,
    }
    coverage = {"tx_left": 3, "tx_right": 3}
    return transcripts, coverage


def generate_fixture(out_dir: Path, fixture_name: str, insert: int) -> dict[str, object]:
    if fixture_name == "tiny_alt_isoform":
        transcripts, coverage = tiny_alt_isoform_transcripts()
    elif fixture_name == "ambiguous_paralog":
        transcripts, coverage = ambiguous_paralog_transcripts()
    elif fixture_name == "compact_fusion":
        transcripts, coverage = compact_fusion_transcripts()
    else:
        raise ValueError(f"unknown fixture: {fixture_name}")
    return write_fixture_files(out_dir, fixture_name, transcripts, insert, coverage)


def run_raptor_workflow_case(
    out_dir: Path,
    fixture: dict[str, object],
    oracle_fasta: Path,
    min_match_coverage: float,
    use_gpu: bool,
) -> dict[str, object]:
    workflow_dir = out_dir / ("raptor_workflow_gpu" if use_gpu else "raptor_workflow")
    workflow_fasta = workflow_dir / "raptor_trinity.fasta.gz"
    command = [
        "cargo",
        "run",
        "--quiet",
    ]
    if use_gpu:
        command.extend(["--features", "gpu"])
    command.extend(
        [
        "--",
        "trinity",
        "--input1",
        fixture["paths"]["r1_fastq"],
        "--input2",
        fixture["paths"]["r2_fastq"],
        "--output-dir",
        str(workflow_dir),
        "--output-fasta",
        str(workflow_fasta),
        "--min-len",
        "25",
        ]
    )
    if use_gpu:
        command.append("--gpu")
    result = run_command(command, ROOT)
    metrics = {
        "gpu_requested": use_gpu,
        "output_exists": workflow_fasta.exists(),
        "output_fasta_bytes": file_size_bytes(workflow_fasta),
        "output_dir_footprint": directory_footprint(workflow_dir),
        "workflow_report": str(workflow_dir / "raptor_trinity_report.json"),
    }
    workflow_report = workflow_dir / "raptor_trinity_report.json"
    if workflow_report.exists():
        workflow_payload = json.loads(workflow_report.read_text(encoding="utf-8"))
        metrics["component_count"] = workflow_payload.get("component_count")
        metrics["components_json"] = workflow_payload.get("components_json")
        metrics["component_clustering"] = workflow_payload.get("component_clustering")
        metrics["component_graph_count"] = workflow_payload.get("component_graph_count")
        metrics["component_graphs_json"] = workflow_payload.get("component_graphs_json")
        metrics["component_paths_fasta"] = workflow_payload.get("component_paths_fasta")
        metrics["component_selected_isoforms_fasta"] = workflow_payload.get(
            "component_selected_isoforms_fasta"
        )
        metrics["component_selected_isoforms_json"] = workflow_payload.get(
            "component_selected_isoforms_json"
        )
        metrics["component_isoform_candidates_json"] = workflow_payload.get(
            "component_isoform_candidates_json"
        )
        components_json = workflow_payload.get("components_json")
        if components_json:
            metrics.update(component_evidence_metrics(Path(components_json)))
        component_graphs_json = workflow_payload.get("component_graphs_json")
        if component_graphs_json:
            metrics.update(component_graph_metrics(Path(component_graphs_json)))
        selected_isoforms_json = workflow_payload.get("component_selected_isoforms_json")
        if selected_isoforms_json:
            metrics.update(selected_isoform_evidence_metrics(Path(selected_isoforms_json)))
        isoform_candidates_json = workflow_payload.get("component_isoform_candidates_json")
        if isoform_candidates_json:
            metrics.update(isoform_candidate_metrics(Path(isoform_candidates_json)))
        component_paths_fasta = workflow_payload.get("component_paths_fasta")
        if component_paths_fasta:
            metrics.update(
                component_candidate_metrics(
                    Path(component_paths_fasta),
                    Path(workflow_payload["component_selected_isoforms_fasta"])
                    if workflow_payload.get("component_selected_isoforms_fasta")
                    else None,
                    Path(fixture["paths"]["truth_fasta"]),
                    oracle_fasta,
                    min_match_coverage,
                )
            )
    if workflow_fasta.exists():
        lengths = read_fasta_lengths(workflow_fasta)
        metrics.update(
            {
                "transcript_count": len(lengths),
                "total_bases": sum(lengths),
                "n50": n50(lengths),
                "lengths": lengths,
            }
        )
        metrics["truth_recovery"] = truth_recovery_metrics(
            Path(fixture["paths"]["truth_fasta"]), workflow_fasta
        )
        if oracle_fasta.exists():
            metrics["oracle_recovery"] = fasta_recovery_metrics(oracle_fasta, workflow_fasta)
    result["metrics"] = metrics
    return result


def add_cpu_gpu_workflow_comparison(report: dict[str, object]) -> None:
    cpu = report.get("raptor_workflow") or {}
    gpu = report.get("raptor_workflow_gpu") or {}
    cpu_metrics = cpu.get("metrics", {})
    gpu_metrics = gpu.get("metrics", {})
    cpu_selected = cpu_metrics.get("component_selected_isoforms_fasta")
    gpu_selected = gpu_metrics.get("component_selected_isoforms_fasta")
    if not cpu_selected or not gpu_selected:
        return
    cpu_path = Path(cpu_selected)
    gpu_path = Path(gpu_selected)
    if not cpu_path.exists() or not gpu_path.exists():
        return
    gpu_metrics["cpu_selected_isoform_match"] = fasta_precision_recall_metrics(
        cpu_path, gpu_path, 0.95
    )
    cpu_metrics["gpu_selected_isoform_match"] = fasta_precision_recall_metrics(
        gpu_path, cpu_path, 0.95
    )


def run_one_fixture(
    out_dir: Path,
    fixture_name: str,
    insert: int,
    oracle_fasta: Path,
    normalize_raptor: bool,
    assemble_normalized: bool,
    run_raptor_workflow: bool,
    run_raptor_workflow_gpu: bool,
    run_trinity: bool,
    require_trinity: bool,
    freeze_trinity_oracle: bool,
    skip_raptor: bool,
    min_match_coverage: float,
) -> dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    fixture = generate_fixture(out_dir, fixture_name, insert)

    report: dict[str, object] = {
        "fixture": fixture,
        "repo": str(ROOT),
        "raptor_normalize": None,
        "raptor": None,
        "raptor_workflow": None,
        "raptor_workflow_gpu": None,
        "trinity": {
            "available": shutil.which("Trinity") is not None,
            "ran": False,
            "required": require_trinity,
            "freeze_oracle_requested": freeze_trinity_oracle,
            "note": "Trinity is optional unless --require-trinity or --freeze-trinity-oracle is set.",
        },
        "oracle": {
            "path": str(oracle_fasta),
            "available": oracle_fasta.exists(),
            "note": "Checked-in tiny oracle is a frozen stand-in until Trinity is installed or a Trinity oracle is captured.",
        },
        "known_limitations": [
            "Tiny fixture now exercises paired-end input for both Raptor and Trinity, but downstream Raptor path scoring still needs stronger read-pair constraints.",
            "No Trinity parity claim is made from this tiny fixture.",
        ],
    }

    raptor_r1 = Path(fixture["paths"]["r1_fastq"])
    raptor_r2 = Path(fixture["paths"]["r2_fastq"])
    if normalize_raptor:
        normalized_prefix = out_dir / "raptor_normalized" / "reads"
        normalized_prefix.parent.mkdir(parents=True, exist_ok=True)
        command = [
            "cargo",
            "run",
            "--quiet",
            "--",
            "normalize",
            "--input1",
            fixture["paths"]["r1_fastq"],
            "--input2",
            fixture["paths"]["r2_fastq"],
            "--output",
            str(normalized_prefix),
            "--coverage-target",
            "500",
            "--max-reads",
            "5000000",
        ]
        result = run_command(command, ROOT)
        normalized_r1 = Path(f"{normalized_prefix}_R1.fastq.gz")
        normalized_r2 = Path(f"{normalized_prefix}_R2.fastq.gz")
        metrics: dict[str, object] = {
            "input_pairs": fixture["paired_end_pairs"],
            "output_r1": str(normalized_r1),
            "output_r2": str(normalized_r2),
            "output_exists": normalized_r1.exists() and normalized_r2.exists(),
            "output_r1_bytes": file_size_bytes(normalized_r1),
            "output_r2_bytes": file_size_bytes(normalized_r2),
            "output_dir_footprint": directory_footprint(normalized_prefix.parent),
        }
        if normalized_r1.exists() and normalized_r2.exists():
            r1_count = count_fastq_records(normalized_r1)
            r2_count = count_fastq_records(normalized_r2)
            metrics.update(
                {
                    "kept_r1_reads": r1_count,
                    "kept_r2_reads": r2_count,
                    "kept_pairs": min(r1_count, r2_count),
                    "kept_pair_fraction": round(
                        min(r1_count, r2_count) / fixture["paired_end_pairs"], 6
                    )
                    if fixture["paired_end_pairs"]
                    else 0.0,
                    "paired_counts_match": r1_count == r2_count,
                }
            )
            if assemble_normalized:
                raptor_r1 = normalized_r1
                raptor_r2 = normalized_r2
        result["metrics"] = metrics
        report["raptor_normalize"] = result

    if not skip_raptor:
        output_fasta = out_dir / "raptor" / "assembly.fa.gz"
        output_fasta.parent.mkdir(parents=True, exist_ok=True)
        command = [
            "cargo",
            "run",
            "--quiet",
            "--",
            "assemble",
            "--input",
            str(raptor_r1),
            "--input2",
            str(raptor_r2),
            "--output",
            str(output_fasta),
            "--min-len",
            "25",
            "--threads",
            "1",
        ]
        result = run_command(command, ROOT)
        metrics = {
            "output_exists": output_fasta.exists(),
            "output_fasta_bytes": file_size_bytes(output_fasta),
            "output_dir_footprint": directory_footprint(output_fasta.parent),
            "assembled_from_normalized_reads": assemble_normalized and normalize_raptor,
            "input_r1": str(raptor_r1),
            "input_r2": str(raptor_r2),
        }
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
            if oracle_fasta.exists():
                metrics["oracle_recovery"] = fasta_recovery_metrics(oracle_fasta, output_fasta)
        result["metrics"] = metrics
        report["raptor"] = result

    if run_raptor_workflow:
        report["raptor_workflow"] = run_raptor_workflow_case(
            out_dir, fixture, oracle_fasta, min_match_coverage, False
        )

    if run_raptor_workflow_gpu:
        report["raptor_workflow_gpu"] = run_raptor_workflow_case(
            out_dir, fixture, oracle_fasta, min_match_coverage, True
        )

    add_cpu_gpu_workflow_comparison(report)

    if run_trinity:
        trinity_bin = shutil.which("Trinity")
        if not trinity_bin:
            report["trinity"]["ran"] = False
            report["trinity"]["error"] = "Trinity executable not found on PATH"
        else:
            trinity_out = out_dir / "trinity"
            if trinity_out.exists():
                shutil.rmtree(trinity_out)
            command = [
                trinity_bin,
                "--seqType",
                "fq",
                "--left",
                fixture["paths"]["r1_fastq"],
                "--right",
                fixture["paths"]["r2_fastq"],
                "--CPU",
                "1",
                "--max_memory",
                "2G",
                "--output",
                str(trinity_out),
            ]
            report["trinity"]["ran"] = True
            result = run_command(command, ROOT)
            output_fasta = trinity_out / "Trinity.fasta"
            metrics = {
                "output_exists": output_fasta.exists(),
                "output_fasta_bytes": file_size_bytes(output_fasta),
                "output_dir_footprint": directory_footprint(trinity_out),
            }
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
                if freeze_trinity_oracle and result["exit_code"] == 0:
                    oracle_fasta.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(output_fasta, oracle_fasta)
                    report["oracle"]["path"] = str(oracle_fasta)
                    report["oracle"]["available"] = True
                    report["oracle"]["frozen_from_trinity"] = str(output_fasta)
            result["metrics"] = metrics
            report["trinity"]["result"] = result

    report_path = out_dir / "report.json"
    write_text(report_path, json.dumps(report, indent=2) + "\n")
    report["report_path"] = str(report_path)
    return report


def check_report_thresholds(
    report: dict[str, object],
    min_truth_coverage: float,
    min_oracle_coverage: float,
    min_selected_precision: float,
    min_selected_f1: float,
) -> list[str]:
    failures: list[str] = []
    trinity = report.get("trinity", {})
    if trinity.get("required") and not trinity.get("available"):
        failures.append("Trinity required but executable was not found on PATH")
    if trinity.get("required") and trinity.get("available") and not trinity.get("ran"):
        failures.append("Trinity required but run was not requested")
    if trinity.get("required") and trinity.get("ran"):
        result = trinity.get("result", {})
        if result.get("exit_code") != 0:
            failures.append("required Trinity run failed")
        elif not result.get("metrics", {}).get("output_exists"):
            failures.append("required Trinity run did not produce Trinity.fasta")

    raptor = report.get("raptor")
    raptor_normalize = report.get("raptor_normalize")
    raptor_workflow = report.get("raptor_workflow")
    raptor_workflow_gpu = report.get("raptor_workflow_gpu")
    if raptor_normalize:
        if raptor_normalize.get("exit_code") != 0:
            failures.append("raptor normalize failed")
        normalize_metrics = raptor_normalize.get("metrics", {})
        if not normalize_metrics.get("output_exists"):
            failures.append("raptor normalize did not produce paired outputs")
        if normalize_metrics.get("paired_counts_match") is False:
            failures.append("raptor normalize produced mismatched R1/R2 record counts")

    if not raptor:
        if not raptor_workflow and not raptor_workflow_gpu:
            return failures
    if raptor_workflow:
        if raptor_workflow.get("exit_code") != 0:
            failures.append("raptor trinity workflow failed")
        workflow_metrics = raptor_workflow.get("metrics", {})
        if not workflow_metrics.get("output_exists"):
            failures.append("raptor trinity workflow did not produce assembly output")
        if workflow_metrics.get("component_count", 0) < 1:
            failures.append("raptor trinity workflow did not emit transcript components")
        if workflow_metrics.get("component_clustering") != "sequence_or_read_kmer":
            failures.append("raptor trinity workflow did not use read-kmer-aware clustering")
        if workflow_metrics.get("component_graph_count", 0) < 1:
            failures.append("raptor trinity workflow did not emit component graphs")
        if workflow_metrics.get("component_graph_node_count", 0) < 1:
            failures.append("raptor trinity workflow emitted empty component graphs")
        if workflow_metrics.get("component_graph_edge_count", 0) < 1:
            failures.append("raptor trinity workflow did not emit component graph edges")
        if workflow_metrics.get("component_graph_edge_read_support", 0) < 1:
            failures.append("raptor trinity workflow component graph edges lack read support")
        if workflow_metrics.get("component_graph_edge_pair_support", 0) < 1:
            failures.append("raptor trinity workflow component graph edges lack pair support")
        if workflow_metrics.get("component_graph_edge_observed_kmers", 0) < 1:
            failures.append("raptor trinity workflow component graph edges lack observed k-mer support")
        if workflow_metrics.get("component_graph_read_kmer_node_count", 0) < 1:
            failures.append("raptor trinity workflow component graphs lack read k-mer nodes")
        if workflow_metrics.get("component_graph_read_kmer_edge_count", 0) < 1:
            failures.append("raptor trinity workflow component graphs lack read k-mer edges")
        if workflow_metrics.get(
            "component_graph_read_kmer_node_record_count", 0
        ) != workflow_metrics.get("component_graph_read_kmer_node_count", 0):
            failures.append(
                "raptor trinity workflow component graph read k-mer node records are incomplete"
            )
        if workflow_metrics.get(
            "component_graph_read_kmer_edge_record_count", 0
        ) != workflow_metrics.get("component_graph_read_kmer_edge_count", 0):
            failures.append(
                "raptor trinity workflow component graph read k-mer edge records are incomplete"
            )
        if workflow_metrics.get("component_graph_read_kmer_edge_sample_count", 0) < 1:
            failures.append("raptor trinity workflow component graphs lack read k-mer edge samples")
        if workflow_metrics.get("component_graph_read_kmer_path_count", 0) < 1:
            failures.append("raptor trinity workflow component graphs lack read k-mer paths")
        if workflow_metrics.get(
            "component_graph_read_kmer_path_record_count", 0
        ) != workflow_metrics.get("component_graph_read_kmer_path_count", 0):
            failures.append(
                "raptor trinity workflow component graph read k-mer path records are incomplete"
            )
        if workflow_metrics.get("component_graph_read_kmer_path_edge_count", 0) < 1:
            failures.append("raptor trinity workflow component graph paths lack edges")
        if workflow_metrics.get("component_graph_read_kmer_path_min_support", 0) < 1:
            failures.append("raptor trinity workflow component graph paths lack read support")
        if not workflow_metrics.get("component_paths_fasta_exists"):
            failures.append("raptor trinity workflow did not emit component path transcript FASTA")
        if workflow_metrics.get("component_path_transcript_count", 0) < 1:
            failures.append("raptor trinity workflow component path FASTA is empty")
        if not workflow_metrics.get("component_selected_isoforms_fasta_exists"):
            failures.append("raptor trinity workflow did not emit selected isoform FASTA")
        if workflow_metrics.get("component_selected_isoform_count", 0) < 1:
            failures.append("raptor trinity workflow did not select component isoforms")
        if not workflow_metrics.get("component_selected_isoforms_json_exists"):
            failures.append("raptor trinity workflow did not emit selected isoform evidence JSON")
        if workflow_metrics.get("component_selected_isoform_record_count", 0) != workflow_metrics.get(
            "component_selected_isoform_count", 0
        ):
            failures.append("selected isoform evidence records do not match selected FASTA records")
        if "component_contig_evidence_score_v2" not in workflow_metrics.get(
            "component_selected_selection_methods", []
        ):
            failures.append("selected isoforms do not report the expected scoring method")
        if not workflow_metrics.get("component_selected_component_ranks_are_dense"):
            failures.append("selected isoform component ranks are not dense")
        if not workflow_metrics.get("component_selected_scores_are_descending"):
            failures.append("selected isoforms are not sorted by descending evidence score")
        if workflow_metrics.get("component_selected_evidence_score", 0) < 1:
            failures.append("selected isoforms lack positive evidence score")
        if workflow_metrics.get("component_selected_direct_read_support", 0) < 1:
            failures.append("selected isoforms lack direct read support")
        if workflow_metrics.get("component_selected_direct_pair_support", 0) < 1:
            failures.append("selected isoforms lack direct pair support")
        if workflow_metrics.get("component_selected_read_kmer_path_support", 0) < 1:
            failures.append("selected isoforms lack read k-mer path support")
        if not workflow_metrics.get("component_isoform_candidates_json_exists"):
            failures.append("raptor trinity workflow did not emit isoform candidate JSON")
        if workflow_metrics.get("component_isoform_candidate_count", 0) <= workflow_metrics.get(
            "component_selected_isoform_count", 0
        ):
            failures.append("isoform candidate set does not exceed selected isoform count")
        if workflow_metrics.get("component_isoform_selected_candidate_count", 0) != workflow_metrics.get(
            "component_selected_isoform_count", 0
        ):
            failures.append("selected isoform candidates do not match selected isoform count")
        if workflow_metrics.get("component_isoform_path_candidate_count", 0) < 1:
            failures.append("isoform candidate set lacks read k-mer path candidates")
        if workflow_metrics.get("component_isoform_rejected_path_candidate_count", 0) < 1:
            failures.append("isoform scoring did not reject read k-mer path candidates")
        if not workflow_metrics.get("component_isoform_candidate_ranks_are_dense"):
            failures.append("isoform candidate ranks are not dense")
        if not workflow_metrics.get("component_isoform_candidate_scores_are_descending"):
            failures.append("isoform candidates are not sorted by descending evidence score")
        if workflow_metrics.get("component_assigned_read_count", 0) < 1:
            failures.append("raptor trinity workflow did not assign reads to components")
        if workflow_metrics.get("component_assigned_pair_count", 0) < 1:
            failures.append("raptor trinity workflow did not assign pairs to components")
        truth_recovery = workflow_metrics.get("truth_recovery", {})
        min_coverage = truth_recovery.get("min_best_coverage", 0.0)
        if min_coverage < min_truth_coverage:
            failures.append(
                f"workflow truth recovery below threshold: {min_coverage} < {min_truth_coverage}"
            )
        path_truth_recovery = workflow_metrics.get("component_candidate_truth_recovery", {})
        path_min_coverage = path_truth_recovery.get("min_best_coverage", 0.0)
        if path_min_coverage < min_truth_coverage:
            failures.append(
                f"component path truth recovery below threshold: {path_min_coverage} < {min_truth_coverage}"
            )
        selected_truth_recovery = workflow_metrics.get("component_selected_truth_recovery", {})
        selected_min_coverage = selected_truth_recovery.get("min_best_coverage", 0.0)
        if selected_min_coverage < min_truth_coverage:
            failures.append(
                f"selected component isoform truth recovery below threshold: {selected_min_coverage} < {min_truth_coverage}"
            )
        selected_truth_precision = workflow_metrics.get("component_selected_truth_precision", {})
        selected_precision = selected_truth_precision.get("precision", 0.0)
        selected_f1 = selected_truth_precision.get("f1", 0.0)
        if selected_precision < min_selected_precision:
            failures.append(
                f"selected component isoform precision below threshold: {selected_precision} < {min_selected_precision}"
            )
        if selected_f1 < min_selected_f1:
            failures.append(
                f"selected component isoform F1 below threshold: {selected_f1} < {min_selected_f1}"
            )
        oracle_recovery = workflow_metrics.get("oracle_recovery")
        if oracle_recovery:
            min_oracle = oracle_recovery.get("min_best_coverage", 0.0)
            if min_oracle < min_oracle_coverage:
                failures.append(
                    f"workflow oracle recovery below threshold: {min_oracle} < {min_oracle_coverage}"
                )
        path_oracle_recovery = workflow_metrics.get("component_candidate_oracle_recovery")
        if path_oracle_recovery:
            path_min_oracle = path_oracle_recovery.get("min_best_coverage", 0.0)
            if path_min_oracle < min_oracle_coverage:
                failures.append(
                    f"component path oracle recovery below threshold: {path_min_oracle} < {min_oracle_coverage}"
                )
        selected_oracle_recovery = workflow_metrics.get("component_selected_oracle_recovery")
        if selected_oracle_recovery:
            selected_min_oracle = selected_oracle_recovery.get("min_best_coverage", 0.0)
            if selected_min_oracle < min_oracle_coverage:
                failures.append(
                    f"selected component isoform oracle recovery below threshold: {selected_min_oracle} < {min_oracle_coverage}"
                )

    if raptor_workflow_gpu:
        if raptor_workflow_gpu.get("exit_code") != 0:
            failures.append("raptor trinity GPU-requested workflow failed")
        gpu_workflow_metrics = raptor_workflow_gpu.get("metrics", {})
        if not gpu_workflow_metrics.get("output_exists"):
            failures.append("raptor trinity GPU-requested workflow did not produce assembly output")
        selected_match = gpu_workflow_metrics.get("cpu_selected_isoform_match", {})
        selected_match_precision = selected_match.get("precision", 0.0)
        selected_match_f1 = selected_match.get("f1", 0.0)
        if selected_match_precision < min_selected_precision:
            failures.append(
                f"GPU-requested selected isoform CPU-match precision below threshold: {selected_match_precision} < {min_selected_precision}"
            )
        if selected_match_f1 < min_selected_f1:
            failures.append(
                f"GPU-requested selected isoform CPU-match F1 below threshold: {selected_match_f1} < {min_selected_f1}"
            )

    if not raptor:
        return failures

    if raptor.get("exit_code") != 0:
        failures.append("raptor assemble failed")
        return failures

    metrics = raptor.get("metrics", {})
    truth_recovery = metrics.get("truth_recovery", {})
    min_coverage = truth_recovery.get("min_best_coverage", 0.0)
    if min_coverage < min_truth_coverage:
        failures.append(
            f"truth recovery below threshold: {min_coverage} < {min_truth_coverage}"
        )

    oracle_recovery = metrics.get("oracle_recovery")
    if oracle_recovery:
        min_oracle = oracle_recovery.get("min_best_coverage", 0.0)
        if min_oracle < min_oracle_coverage:
            failures.append(
                f"oracle recovery below threshold: {min_oracle} < {min_oracle_coverage}"
            )

    return failures


def summarize_report(report: dict[str, object]) -> dict[str, object]:
    raptor = report.get("raptor") or {}
    raptor_normalize = report.get("raptor_normalize") or {}
    raptor_workflow = report.get("raptor_workflow") or {}
    raptor_workflow_gpu = report.get("raptor_workflow_gpu") or {}
    trinity = report.get("trinity") or {}
    metrics = raptor.get("metrics", {})
    normalize_metrics = raptor_normalize.get("metrics", {})
    workflow_metrics = raptor_workflow.get("metrics", {})
    gpu_workflow_metrics = raptor_workflow_gpu.get("metrics", {})
    trinity_result = trinity.get("result", {})
    trinity_metrics = trinity_result.get("metrics", {})
    raptor_resources = raptor.get("resource_usage", {})
    normalize_resources = raptor_normalize.get("resource_usage", {})
    workflow_resources = raptor_workflow.get("resource_usage", {})
    gpu_workflow_resources = raptor_workflow_gpu.get("resource_usage", {})
    trinity_resources = trinity_result.get("resource_usage", {})
    raptor_gpu = raptor.get("gpu_usage", {})
    normalize_gpu = raptor_normalize.get("gpu_usage", {})
    workflow_gpu = raptor_workflow.get("gpu_usage", {})
    gpu_workflow_gpu = raptor_workflow_gpu.get("gpu_usage", {})
    trinity_gpu = trinity_result.get("gpu_usage", {})
    return {
        "insert": report["fixture"]["insert"],
        "paired_end_pairs": report["fixture"]["paired_end_pairs"],
        "report_path": report.get("report_path"),
        "raptor_exit_code": raptor.get("exit_code"),
        "raptor_elapsed_seconds": raptor.get("elapsed_seconds"),
        "raptor_max_rss_kb": raptor_resources.get("max_rss_kb"),
        "raptor_user_seconds": raptor_resources.get("user_seconds"),
        "raptor_system_seconds": raptor_resources.get("system_seconds"),
        "raptor_gpu_available": raptor_gpu.get("available"),
        "raptor_gpu_max_memory_used_mib": first_gpu_metric(
            raptor_gpu, "max_memory_used_mib"
        ),
        "raptor_gpu_max_utilization_percent": first_gpu_metric(
            raptor_gpu, "max_utilization_gpu_percent"
        ),
        "raptor_normalize_exit_code": raptor_normalize.get("exit_code"),
        "raptor_normalize_elapsed_seconds": raptor_normalize.get("elapsed_seconds"),
        "raptor_normalize_max_rss_kb": normalize_resources.get("max_rss_kb"),
        "raptor_normalize_user_seconds": normalize_resources.get("user_seconds"),
        "raptor_normalize_system_seconds": normalize_resources.get("system_seconds"),
        "raptor_normalize_gpu_available": normalize_gpu.get("available"),
        "raptor_normalize_gpu_max_memory_used_mib": first_gpu_metric(
            normalize_gpu, "max_memory_used_mib"
        ),
        "raptor_normalize_gpu_max_utilization_percent": first_gpu_metric(
            normalize_gpu, "max_utilization_gpu_percent"
        ),
        "raptor_normalize_output_bytes": (
            (normalize_metrics.get("output_r1_bytes") or 0)
            + (normalize_metrics.get("output_r2_bytes") or 0)
        )
        if raptor_normalize
        else None,
        "normalized_kept_pairs": normalize_metrics.get("kept_pairs"),
        "normalized_kept_pair_fraction": normalize_metrics.get("kept_pair_fraction"),
        "assembled_from_normalized_reads": metrics.get("assembled_from_normalized_reads"),
        "raptor_workflow_exit_code": raptor_workflow.get("exit_code"),
        "raptor_workflow_elapsed_seconds": raptor_workflow.get("elapsed_seconds"),
        "raptor_workflow_max_rss_kb": workflow_resources.get("max_rss_kb"),
        "raptor_workflow_user_seconds": workflow_resources.get("user_seconds"),
        "raptor_workflow_system_seconds": workflow_resources.get("system_seconds"),
        "raptor_workflow_gpu_available": workflow_gpu.get("available"),
        "raptor_workflow_gpu_max_memory_used_mib": first_gpu_metric(
            workflow_gpu, "max_memory_used_mib"
        ),
        "raptor_workflow_gpu_max_utilization_percent": first_gpu_metric(
            workflow_gpu, "max_utilization_gpu_percent"
        ),
        "raptor_workflow_gpu_max_power_draw_watts": first_gpu_metric(
            workflow_gpu, "max_power_draw_watts"
        ),
        "raptor_gpu_workflow_exit_code": raptor_workflow_gpu.get("exit_code"),
        "raptor_gpu_workflow_elapsed_seconds": raptor_workflow_gpu.get("elapsed_seconds"),
        "raptor_gpu_workflow_max_rss_kb": gpu_workflow_resources.get("max_rss_kb"),
        "raptor_gpu_workflow_gpu_available": gpu_workflow_gpu.get("available"),
        "raptor_gpu_workflow_gpu_max_memory_used_mib": first_gpu_metric(
            gpu_workflow_gpu, "max_memory_used_mib"
        ),
        "raptor_gpu_workflow_gpu_max_utilization_percent": first_gpu_metric(
            gpu_workflow_gpu, "max_utilization_gpu_percent"
        ),
        "raptor_gpu_workflow_gpu_max_power_draw_watts": first_gpu_metric(
            gpu_workflow_gpu, "max_power_draw_watts"
        ),
        "raptor_gpu_workflow_selected_cpu_match_precision": gpu_workflow_metrics.get(
            "cpu_selected_isoform_match", {}
        ).get("precision"),
        "raptor_gpu_workflow_selected_cpu_match_recall": gpu_workflow_metrics.get(
            "cpu_selected_isoform_match", {}
        ).get("recall"),
        "raptor_gpu_workflow_selected_cpu_match_f1": gpu_workflow_metrics.get(
            "cpu_selected_isoform_match", {}
        ).get("f1"),
        "raptor_gpu_workflow_selected_isoform_lengths": gpu_workflow_metrics.get(
            "component_selected_isoform_lengths"
        ),
        "raptor_workflow_output_fasta_bytes": workflow_metrics.get("output_fasta_bytes"),
        "raptor_workflow_output_dir_bytes": workflow_metrics.get(
            "output_dir_footprint", {}
        ).get("total_bytes"),
        "raptor_workflow_output_file_count": workflow_metrics.get(
            "output_dir_footprint", {}
        ).get("file_count"),
        "workflow_lengths": workflow_metrics.get("lengths"),
        "workflow_n50": workflow_metrics.get("n50"),
        "workflow_component_count": workflow_metrics.get("component_count"),
        "workflow_component_clustering": workflow_metrics.get("component_clustering"),
        "workflow_component_graph_count": workflow_metrics.get("component_graph_count"),
        "workflow_component_graph_nodes": workflow_metrics.get("component_graph_node_count"),
        "workflow_component_graph_edges": workflow_metrics.get("component_graph_edge_count"),
        "workflow_component_graph_edge_read_support": workflow_metrics.get(
            "component_graph_edge_read_support"
        ),
        "workflow_component_graph_edge_pair_support": workflow_metrics.get(
            "component_graph_edge_pair_support"
        ),
        "workflow_component_graph_edge_observed_kmers": workflow_metrics.get(
            "component_graph_edge_observed_kmers"
        ),
        "workflow_component_graph_read_kmer_nodes": workflow_metrics.get(
            "component_graph_read_kmer_node_count"
        ),
        "workflow_component_graph_read_kmer_edges": workflow_metrics.get(
            "component_graph_read_kmer_edge_count"
        ),
        "workflow_component_graph_read_kmer_node_records": workflow_metrics.get(
            "component_graph_read_kmer_node_record_count"
        ),
        "workflow_component_graph_read_kmer_edge_records": workflow_metrics.get(
            "component_graph_read_kmer_edge_record_count"
        ),
        "workflow_component_graph_read_kmer_edge_samples": workflow_metrics.get(
            "component_graph_read_kmer_edge_sample_count"
        ),
        "workflow_component_graph_read_kmer_paths": workflow_metrics.get(
            "component_graph_read_kmer_path_count"
        ),
        "workflow_component_graph_read_kmer_path_records": workflow_metrics.get(
            "component_graph_read_kmer_path_record_count"
        ),
        "workflow_component_graph_read_kmer_path_edges": workflow_metrics.get(
            "component_graph_read_kmer_path_edge_count"
        ),
        "workflow_component_graph_read_kmer_path_min_support": workflow_metrics.get(
            "component_graph_read_kmer_path_min_support"
        ),
        "workflow_component_paths_fasta": workflow_metrics.get("component_paths_fasta"),
        "workflow_component_path_transcript_count": workflow_metrics.get(
            "component_path_transcript_count"
        ),
        "workflow_component_path_lengths": workflow_metrics.get("component_path_lengths"),
        "workflow_component_selected_isoforms_fasta": workflow_metrics.get(
            "component_selected_isoforms_fasta"
        ),
        "workflow_component_selected_isoforms_json": workflow_metrics.get(
            "component_selected_isoforms_json"
        ),
        "workflow_component_isoform_candidates_json": workflow_metrics.get(
            "component_isoform_candidates_json"
        ),
        "workflow_component_isoform_candidate_count": workflow_metrics.get(
            "component_isoform_candidate_count"
        ),
        "workflow_component_isoform_selected_candidate_count": workflow_metrics.get(
            "component_isoform_selected_candidate_count"
        ),
        "workflow_component_isoform_rejected_candidate_count": workflow_metrics.get(
            "component_isoform_rejected_candidate_count"
        ),
        "workflow_component_isoform_path_candidate_count": workflow_metrics.get(
            "component_isoform_path_candidate_count"
        ),
        "workflow_component_isoform_rejected_path_candidate_count": workflow_metrics.get(
            "component_isoform_rejected_path_candidate_count"
        ),
        "workflow_component_isoform_candidate_ranks_dense": workflow_metrics.get(
            "component_isoform_candidate_ranks_are_dense"
        ),
        "workflow_component_isoform_candidate_scores_descending": workflow_metrics.get(
            "component_isoform_candidate_scores_are_descending"
        ),
        "workflow_component_isoform_candidate_max_score": workflow_metrics.get(
            "component_isoform_candidate_max_score"
        ),
        "workflow_component_isoform_rejected_candidate_max_score": workflow_metrics.get(
            "component_isoform_rejected_candidate_max_score"
        ),
        "workflow_component_selected_isoform_count": workflow_metrics.get(
            "component_selected_isoform_count"
        ),
        "workflow_component_selected_isoform_records": workflow_metrics.get(
            "component_selected_isoform_record_count"
        ),
        "workflow_component_selected_selection_methods": workflow_metrics.get(
            "component_selected_selection_methods"
        ),
        "workflow_component_selected_ranks_dense": workflow_metrics.get(
            "component_selected_component_ranks_are_dense"
        ),
        "workflow_component_selected_scores_descending": workflow_metrics.get(
            "component_selected_scores_are_descending"
        ),
        "workflow_component_selected_evidence_score": workflow_metrics.get(
            "component_selected_evidence_score"
        ),
        "workflow_component_selected_max_evidence_score": workflow_metrics.get(
            "component_selected_max_evidence_score"
        ),
        "workflow_component_selected_isoform_lengths": workflow_metrics.get(
            "component_selected_isoform_lengths"
        ),
        "workflow_component_selected_direct_read_support": workflow_metrics.get(
            "component_selected_direct_read_support"
        ),
        "workflow_component_selected_direct_pair_support": workflow_metrics.get(
            "component_selected_direct_pair_support"
        ),
        "workflow_component_selected_read_kmer_path_support": workflow_metrics.get(
            "component_selected_read_kmer_path_support"
        ),
        "workflow_component_selected_max_read_kmer_path_support": workflow_metrics.get(
            "component_selected_max_read_kmer_path_support"
        ),
        "workflow_component_path_truth_min_coverage": workflow_metrics.get(
            "component_candidate_truth_recovery", {}
        ).get("min_best_coverage"),
        "workflow_component_path_oracle_min_coverage": workflow_metrics.get(
            "component_candidate_oracle_recovery", {}
        ).get("min_best_coverage"),
        "workflow_component_selected_truth_min_coverage": workflow_metrics.get(
            "component_selected_truth_recovery", {}
        ).get("min_best_coverage"),
        "workflow_component_selected_precision": workflow_metrics.get(
            "component_selected_truth_precision", {}
        ).get("precision"),
        "workflow_component_selected_recall": workflow_metrics.get(
            "component_selected_truth_precision", {}
        ).get("recall"),
        "workflow_component_selected_f1": workflow_metrics.get(
            "component_selected_truth_precision", {}
        ).get("f1"),
        "workflow_component_selected_true_positive": workflow_metrics.get(
            "component_selected_truth_precision", {}
        ).get("true_positive"),
        "workflow_component_selected_false_positive": workflow_metrics.get(
            "component_selected_truth_precision", {}
        ).get("false_positive"),
        "workflow_component_selected_false_negative": workflow_metrics.get(
            "component_selected_truth_precision", {}
        ).get("false_negative"),
        "workflow_component_selected_oracle_min_coverage": workflow_metrics.get(
            "component_selected_oracle_recovery", {}
        ).get("min_best_coverage"),
        "workflow_component_assigned_reads": workflow_metrics.get(
            "component_assigned_read_count"
        ),
        "workflow_component_assigned_pairs": workflow_metrics.get(
            "component_assigned_pair_count"
        ),
        "workflow_truth_min_coverage": workflow_metrics.get("truth_recovery", {}).get(
            "min_best_coverage"
        ),
        "workflow_oracle_min_coverage": workflow_metrics.get("oracle_recovery", {}).get(
            "min_best_coverage"
        ),
        "lengths": metrics.get("lengths"),
        "n50": metrics.get("n50"),
        "truth_min_coverage": metrics.get("truth_recovery", {}).get("min_best_coverage"),
        "oracle_min_coverage": metrics.get("oracle_recovery", {}).get("min_best_coverage"),
        "trinity_available": trinity.get("available"),
        "trinity_ran": trinity.get("ran"),
        "trinity_exit_code": trinity_result.get("exit_code"),
        "trinity_elapsed_seconds": trinity_result.get("elapsed_seconds"),
        "trinity_max_rss_kb": trinity_resources.get("max_rss_kb"),
        "trinity_user_seconds": trinity_resources.get("user_seconds"),
        "trinity_system_seconds": trinity_resources.get("system_seconds"),
        "trinity_gpu_available": trinity_gpu.get("available"),
        "trinity_gpu_max_memory_used_mib": first_gpu_metric(
            trinity_gpu, "max_memory_used_mib"
        ),
        "trinity_gpu_max_utilization_percent": first_gpu_metric(
            trinity_gpu, "max_utilization_gpu_percent"
        ),
        "trinity_gpu_max_power_draw_watts": first_gpu_metric(
            trinity_gpu, "max_power_draw_watts"
        ),
        "trinity_output_fasta_bytes": trinity_metrics.get("output_fasta_bytes"),
        "trinity_output_dir_bytes": trinity_metrics.get("output_dir_footprint", {}).get(
            "total_bytes"
        ),
        "trinity_output_file_count": trinity_metrics.get("output_dir_footprint", {}).get(
            "file_count"
        ),
        "trinity_lengths": trinity_metrics.get("lengths"),
        "trinity_n50": trinity_metrics.get("n50"),
        "trinity_truth_min_coverage": trinity_metrics.get("truth_recovery", {}).get(
            "min_best_coverage"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--fixture",
        choices=["tiny_alt_isoform", "ambiguous_paralog", "compact_fusion"],
        default=DEFAULT_FIXTURE,
    )
    parser.add_argument("--skip-raptor", action="store_true")
    parser.add_argument(
        "--normalize-raptor",
        action="store_true",
        help="Run raptor normalize on the paired FASTQs before assembly and record kept-pair metrics",
    )
    parser.add_argument(
        "--assemble-normalized",
        action="store_true",
        help="Assemble Raptor output from --normalize-raptor instead of raw fixture reads",
    )
    parser.add_argument(
        "--run-raptor-workflow",
        action="store_true",
        help="Run the raptor trinity end-to-end CLI and record its output metrics",
    )
    parser.add_argument(
        "--run-raptor-workflow-gpu",
        action="store_true",
        help="Run the raptor trinity end-to-end CLI with --gpu and compare to the CPU-requested workflow when both are present",
    )
    parser.add_argument("--run-trinity", action="store_true")
    parser.add_argument("--oracle-fasta", type=Path, default=DEFAULT_ORACLE)
    parser.add_argument("--insert", type=int, default=160)
    parser.add_argument(
        "--require-trinity",
        action="store_true",
        help="Fail when --run-trinity cannot execute Trinity and produce Trinity.fasta",
    )
    parser.add_argument(
        "--freeze-trinity-oracle",
        action="store_true",
        help="After a successful --run-trinity run, copy Trinity.fasta to --oracle-fasta",
    )
    parser.add_argument(
        "--insert-sweep",
        default=None,
        help="Comma-separated paired insert sizes to run as a sweep, e.g. 110,160,200",
    )
    parser.add_argument("--min-truth-coverage", type=float, default=0.95)
    parser.add_argument("--min-oracle-coverage", type=float, default=0.95)
    parser.add_argument(
        "--min-selected-match-coverage",
        type=float,
        default=0.95,
        help="Minimum reciprocal coverage for a selected isoform to count as a truth match",
    )
    parser.add_argument(
        "--min-selected-precision",
        type=float,
        default=0.95,
        help="Minimum selected isoform precision against truth transcripts",
    )
    parser.add_argument(
        "--min-selected-f1",
        type=float,
        default=0.95,
        help="Minimum selected isoform F1 against truth transcripts",
    )
    args = parser.parse_args()

    out_dir = args.out_dir.resolve()
    if args.fixture != DEFAULT_FIXTURE and args.out_dir == DEFAULT_OUT:
        out_dir = ROOT / "target" / "trinity_parity" / args.fixture
    oracle_fasta = args.oracle_fasta
    if args.fixture != DEFAULT_FIXTURE and args.oracle_fasta == DEFAULT_ORACLE:
        oracle_fasta = (
            ROOT / "bench" / "trinity_parity" / "oracles" / f"{args.fixture}.fa"
        )
    if args.assemble_normalized and not args.normalize_raptor:
        print("--assemble-normalized requires --normalize-raptor", file=sys.stderr)
        return 2
    run_trinity = args.run_trinity or args.require_trinity or args.freeze_trinity_oracle
    require_trinity = args.require_trinity or args.freeze_trinity_oracle

    inserts = parse_insert_sweep(args.insert_sweep) if args.insert_sweep else [args.insert]
    reports = []
    failures = []
    for insert in inserts:
        run_dir = out_dir if len(inserts) == 1 else out_dir / f"insert_{insert}"
        report = run_one_fixture(
            run_dir,
            args.fixture,
            insert,
            oracle_fasta,
            args.normalize_raptor,
            args.assemble_normalized,
            args.run_raptor_workflow,
            args.run_raptor_workflow_gpu,
            run_trinity,
            require_trinity,
            args.freeze_trinity_oracle,
            args.skip_raptor,
            args.min_selected_match_coverage,
        )
        reports.append(report)
        failures.extend(
            f"insert {insert}: {failure}"
            for failure in check_report_thresholds(
                report,
                args.min_truth_coverage,
                args.min_oracle_coverage,
                args.min_selected_precision,
                args.min_selected_f1,
            )
        )
        print(f"wrote {report['report_path']}")

    if len(reports) > 1:
        summary = {
            "insert_sweep": inserts,
            "summaries": [summarize_report(report) for report in reports],
            "failures": failures,
        }
        summary_path = out_dir / "insert_sweep_report.json"
        write_text(summary_path, json.dumps(summary, indent=2) + "\n")
        print(f"wrote {summary_path}")

    if failures:
        for failure in failures:
            print(failure, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
