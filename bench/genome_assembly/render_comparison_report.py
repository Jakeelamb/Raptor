#!/usr/bin/env python3

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SUMMARY_CSV = ROOT / "benchmark_summary.csv"
OUTPUT_MD = ROOT.parent.parent / "docs" / "genome-benchmark-comparison.md"


def load_rows() -> list[dict[str, str]]:
    if not SUMMARY_CSV.exists():
        return []
    with SUMMARY_CSV.open() as handle:
        return list(csv.DictReader(handle))


def val(row: dict[str, str], key: str) -> str:
    value = row.get(key, "N/A").strip()
    return value or "N/A"


def build_table(rows: list[dict[str, str]]) -> list[str]:
    lines = [
        "| Tool | Source Run | Time (s) | Peak RSS (KB) | Contigs | N50 | Genome Fraction % | Misassemblies | NGA50 | Duplication Ratio |",
        "|------|------------|----------|---------------|---------|-----|-------------------|---------------|-------|-------------------|",
    ]
    for row in rows:
        lines.append(
            "| {tool} | {run} | {time_seconds} | {peak_memory_kb} | {num_contigs} | {n50} | "
            "{genome_fraction_percent} | {misassemblies} | {nga50} | {duplication_ratio} |".format(
                **row
            )
        )
    return lines


def metric_ready(row: dict[str, str]) -> bool:
    return any(
        val(row, key) != "N/A"
        for key in ("time_seconds", "num_contigs", "n50", "genome_fraction_percent", "nga50")
    )


def select_latest_rows(dataset_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    latest: dict[str, dict[str, str]] = {}
    for row in sorted(dataset_rows, key=lambda r: (val(r, "run"), metric_ready(r))):
        tool = val(row, "tool")
        current = latest.get(tool)
        if current is None:
            latest[tool] = row
            continue
        if val(row, "run") > val(current, "run"):
            latest[tool] = row
            continue
        if val(row, "run") == val(current, "run") and metric_ready(row) and not metric_ready(current):
            latest[tool] = row

    return sorted(latest.values(), key=lambda r: (val(r, "tool") != "raptor", val(r, "tool")))


def main() -> None:
    rows = load_rows()
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[val(row, "dataset")].append(row)

    lines = [
        "# Genome Benchmark Comparison",
        "",
        "This document is generated from `bench/genome_assembly/benchmark_summary.csv`.",
        "",
        "It is the publication-facing comparison view for Raptor's genome assembly benchmarks.",
        "",
        "## Current Status",
        "",
        "- Raptor baseline runs are available in the current repo state.",
        "- Comparator conclusions should only be made for datasets where the same summary contains both Raptor and comparator rows with non-`N/A` metrics.",
        "",
    ]

    if not grouped:
        lines.extend(
            [
                "## No Benchmark Rows",
                "",
                "Run `python3 ./bench/genome_assembly/summarize_results.py` after benchmark execution to populate this report.",
                "",
            ]
        )
    else:
        for dataset in sorted(grouped):
            dataset_rows = select_latest_rows(grouped[dataset])
            tools = ", ".join(val(r, "tool") for r in dataset_rows)
            comparator_ready = any(
                val(r, "tool") != "raptor" and metric_ready(r) for r in dataset_rows
            ) and any(val(r, "tool") == "raptor" and metric_ready(r) for r in dataset_rows)

            lines.extend(
                [
                    f"## Dataset: `{dataset}`",
                    "",
                    f"Latest complete rows by tool: {tools}",
                    "",
                ]
            )

            lines.extend(build_table(dataset_rows))
            lines.append("")

            if comparator_ready:
                lines.append("Comparator metrics are present for this dataset.")
            else:
                lines.append("Comparator metrics are not yet complete for this dataset; treat this as an execution record, not a winner table.")
            lines.append("")

    lines.extend(
        [
            "## Update Workflow",
            "",
            "```bash",
            "python3 ./bench/genome_assembly/summarize_results.py",
            "python3 ./bench/genome_assembly/render_comparison_report.py",
            "```",
            "",
        ]
    )

    OUTPUT_MD.write_text("\n".join(lines) + "\n")
    print(f"Wrote comparison report to {OUTPUT_MD}")


if __name__ == "__main__":
    main()
