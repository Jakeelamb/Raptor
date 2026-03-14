#!/usr/bin/env python3
"""Run RepeatMasker Rust fixture parity and runtime checks.

The default suite targets the production `.out` postprocessing path. Legacy
old-format `.out` / `.cat` fixtures can be included explicitly to track
remaining compatibility work without blocking the primary publication story.
"""

from __future__ import annotations

import argparse
import dataclasses
import difflib
import json
import pathlib
import shutil
import statistics
import subprocess
import time


ROOT = pathlib.Path(__file__).resolve().parents[2]
REPEATMASKER_ROOT = ROOT / "third_party" / "repeatmasker" / "RepeatMasker-master"
DEFAULT_REPORT_DIR = ROOT / "artifacts" / "repeatmasker_fixture_benchmark"
RUST_BINARY = ROOT / "target" / "release" / "process_repeats_rs"


@dataclasses.dataclass(frozen=True)
class FixtureCase:
    name: str
    suite: str
    annotations: pathlib.Path
    fasta: pathlib.Path
    expected_tbl_cmp: pathlib.Path
    expected_masked: pathlib.Path
    mask_mode: str
    ignore_blank_mask_lines: bool = False


@dataclasses.dataclass
class CaseResult:
    name: str
    suite: str
    mask_mode: str
    total_bases: int
    median_runtime_ms: float
    runtimes_ms: list[float]
    tbl_match: bool
    masked_match: bool
    tbl_diff_preview: list[str]
    masked_diff_preview: list[str]
    return_code: int
    stdout: str
    stderr: str


def fixture_cases(include_legacy: bool) -> list[FixtureCase]:
    primary = [
        FixtureCase(
            name="small-1-default-out",
            suite="primary-out",
            annotations=REPEATMASKER_ROOT / "t/seqs/general/small-1-rcmp-1/small-1.fa.out",
            fasta=REPEATMASKER_ROOT / "t/seqs/general/small-1.fa",
            expected_tbl_cmp=REPEATMASKER_ROOT
            / "t/seqs/general/small-1-rcmp-1/small-1.fa.tbl.cmp",
            expected_masked=REPEATMASKER_ROOT
            / "t/seqs/general/small-1-rcmp-1/small-1.fa.masked",
            mask_mode="n",
        ),
        FixtureCase(
            name="small-1-xsmall-out",
            suite="primary-out",
            annotations=REPEATMASKER_ROOT / "t/seqs/general/small-1-rcmp-2/small-1.fa.out",
            fasta=REPEATMASKER_ROOT / "t/seqs/general/small-1.fa",
            expected_tbl_cmp=REPEATMASKER_ROOT
            / "t/seqs/general/small-1-rcmp-2/small-1.fa.tbl.cmp",
            expected_masked=REPEATMASKER_ROOT
            / "t/seqs/general/small-1-rcmp-2/small-1.fa.masked",
            mask_mode="xsmall",
        ),
        FixtureCase(
            name="small-1-x-out",
            suite="primary-out",
            annotations=REPEATMASKER_ROOT / "t/seqs/general/small-1-rcmp-3/small-1.fa.out",
            fasta=REPEATMASKER_ROOT / "t/seqs/general/small-1.fa",
            expected_tbl_cmp=REPEATMASKER_ROOT
            / "t/seqs/general/small-1-rcmp-3/small-1.fa.tbl.cmp",
            expected_masked=REPEATMASKER_ROOT
            / "t/seqs/general/small-1-rcmp-3/small-1.fa.masked",
            mask_mode="x",
        ),
        FixtureCase(
            name="hum-1-default-out",
            suite="primary-out",
            annotations=REPEATMASKER_ROOT / "t/seqs/general/hum-1-rcmp-1/hum-1.fa.out",
            fasta=REPEATMASKER_ROOT / "t/seqs/general/hum-1.fa",
            expected_tbl_cmp=REPEATMASKER_ROOT
            / "t/seqs/general/hum-1-rcmp-1/hum-1.fa.tbl.cmp",
            expected_masked=REPEATMASKER_ROOT
            / "t/seqs/general/hum-1-rcmp-1/hum-1.fa.masked",
            mask_mode="n",
            ignore_blank_mask_lines=True,
        ),
    ]
    if not include_legacy:
        return primary

    legacy = [
        FixtureCase(
            name="is1-default-out",
            suite="legacy-out",
            annotations=REPEATMASKER_ROOT / "t/seqs/ISElements/is1-rcmp-1/is1.fa.out",
            fasta=REPEATMASKER_ROOT / "t/seqs/ISElements/is1.fa",
            expected_tbl_cmp=REPEATMASKER_ROOT
            / "t/seqs/ISElements/is1-rcmp-1/is1.fa.tbl.cmp",
            expected_masked=REPEATMASKER_ROOT
            / "t/seqs/ISElements/is1-rcmp-1/is1.fa.masked",
            mask_mode="n",
        ),
        FixtureCase(
            name="small-1-default-cat",
            suite="legacy-cat",
            annotations=REPEATMASKER_ROOT / "t/seqs/general/small-1-rcmp-1/small-1.fa.cat",
            fasta=REPEATMASKER_ROOT / "t/seqs/general/small-1.fa",
            expected_tbl_cmp=REPEATMASKER_ROOT
            / "t/seqs/general/small-1-rcmp-1/small-1.fa.tbl.cmp",
            expected_masked=REPEATMASKER_ROOT
            / "t/seqs/general/small-1-rcmp-1/small-1.fa.masked",
            mask_mode="n",
        ),
        FixtureCase(
            name="small-1-xsmall-cat",
            suite="legacy-cat",
            annotations=REPEATMASKER_ROOT / "t/seqs/general/small-1-rcmp-2/small-1.fa.cat",
            fasta=REPEATMASKER_ROOT / "t/seqs/general/small-1.fa",
            expected_tbl_cmp=REPEATMASKER_ROOT
            / "t/seqs/general/small-1-rcmp-2/small-1.fa.tbl.cmp",
            expected_masked=REPEATMASKER_ROOT
            / "t/seqs/general/small-1-rcmp-2/small-1.fa.masked",
            mask_mode="xsmall",
        ),
        FixtureCase(
            name="small-1-x-cat",
            suite="legacy-cat",
            annotations=REPEATMASKER_ROOT / "t/seqs/general/small-1-rcmp-3/small-1.fa.cat",
            fasta=REPEATMASKER_ROOT / "t/seqs/general/small-1.fa",
            expected_tbl_cmp=REPEATMASKER_ROOT
            / "t/seqs/general/small-1-rcmp-3/small-1.fa.tbl.cmp",
            expected_masked=REPEATMASKER_ROOT
            / "t/seqs/general/small-1-rcmp-3/small-1.fa.masked",
            mask_mode="x",
        ),
        FixtureCase(
            name="hum-1-default-cat",
            suite="legacy-cat",
            annotations=REPEATMASKER_ROOT / "t/seqs/general/hum-1-rcmp-1/hum-1.fa.cat",
            fasta=REPEATMASKER_ROOT / "t/seqs/general/hum-1.fa",
            expected_tbl_cmp=REPEATMASKER_ROOT
            / "t/seqs/general/hum-1-rcmp-1/hum-1.fa.tbl.cmp",
            expected_masked=REPEATMASKER_ROOT
            / "t/seqs/general/hum-1-rcmp-1/hum-1.fa.masked",
            mask_mode="n",
            ignore_blank_mask_lines=True,
        ),
        FixtureCase(
            name="is1-default-cat",
            suite="legacy-cat",
            annotations=REPEATMASKER_ROOT / "t/seqs/ISElements/is1-rcmp-1/is1.fa.cat",
            fasta=REPEATMASKER_ROOT / "t/seqs/ISElements/is1.fa",
            expected_tbl_cmp=REPEATMASKER_ROOT
            / "t/seqs/ISElements/is1-rcmp-1/is1.fa.tbl.cmp",
            expected_masked=REPEATMASKER_ROOT
            / "t/seqs/ISElements/is1-rcmp-1/is1.fa.masked",
            mask_mode="n",
        ),
    ]
    return primary + legacy


def ensure_release_binary(force_build: bool) -> pathlib.Path:
    if force_build or not RUST_BINARY.exists():
        subprocess.run(
            ["cargo", "build", "--release", "-p", "repeatmasker-rs", "--bin", "process_repeats_rs"],
            cwd=ROOT,
            check=True,
        )
    return RUST_BINARY


def fasta_total_bases(path: pathlib.Path) -> int:
    total = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def normalize_tbl_cmp(text: str) -> str:
    lines = []
    for line in text.splitlines():
        if any(
            token in line
            for token in (
                "RepeatMasker",
                "cross_match",
                "RepBase",
                "The query species was assumed",
            )
        ):
            continue
        lines.append(" ".join(line.rstrip().split()))
    return "\n".join(lines).strip()


def normalize_masked(text: str, ignore_blank_lines: bool) -> str:
    lines = text.splitlines()
    if ignore_blank_lines:
        lines = [line for line in lines if line.strip()]
    return "\n".join(lines).strip()


def preview_diff(expected: str, observed: str, label: str) -> list[str]:
    diff = list(
        difflib.unified_diff(
            expected.splitlines(),
            observed.splitlines(),
            fromfile=f"expected/{label}",
            tofile=f"observed/{label}",
            n=3,
        )
    )
    return diff[:20]


def rust_command(binary: pathlib.Path, case: FixtureCase, out_dir: pathlib.Path) -> list[str]:
    cmd = [
        str(binary),
        "--annotations",
        str(case.annotations),
        "--tbl",
        str(out_dir / f"{case.name}.tbl"),
        "--fasta",
        str(case.fasta),
        "--masked-output",
        str(out_dir / f"{case.name}.masked"),
        "--stats-json",
        str(out_dir / f"{case.name}.json"),
    ]
    if case.mask_mode == "xsmall":
        cmd.append("--xsmall")
    elif case.mask_mode == "x":
        cmd.append("--x")
    return cmd


def run_case(binary: pathlib.Path, case: FixtureCase, repeat_count: int, out_dir: pathlib.Path) -> CaseResult:
    case_dir = out_dir / case.name
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True)

    runtimes_ms: list[float] = []
    stdout = ""
    stderr = ""
    return_code = 0
    for iteration in range(repeat_count):
        run_dir = case_dir / f"run_{iteration + 1}"
        run_dir.mkdir()
        start = time.perf_counter()
        result = subprocess.run(
            rust_command(binary, case, run_dir),
            cwd=ROOT,
            text=True,
            capture_output=True,
        )
        runtimes_ms.append((time.perf_counter() - start) * 1000.0)
        stdout = result.stdout
        stderr = result.stderr
        return_code = result.returncode
        if result.returncode != 0:
            break

    latest_run = case_dir / f"run_{min(repeat_count, len(runtimes_ms))}"
    tbl_path = latest_run / f"{case.name}.tbl"
    masked_path = latest_run / f"{case.name}.masked"

    expected_tbl = normalize_tbl_cmp(case.expected_tbl_cmp.read_text(encoding="utf-8"))
    observed_tbl = normalize_tbl_cmp(tbl_path.read_text(encoding="utf-8")) if tbl_path.exists() else ""
    expected_masked = normalize_masked(
        case.expected_masked.read_text(encoding="utf-8"),
        case.ignore_blank_mask_lines,
    )
    observed_masked = (
        normalize_masked(masked_path.read_text(encoding="utf-8"), case.ignore_blank_mask_lines)
        if masked_path.exists()
        else ""
    )

    return CaseResult(
        name=case.name,
        suite=case.suite,
        mask_mode=case.mask_mode,
        total_bases=fasta_total_bases(case.fasta),
        median_runtime_ms=statistics.median(runtimes_ms) if runtimes_ms else float("nan"),
        runtimes_ms=runtimes_ms,
        tbl_match=expected_tbl == observed_tbl,
        masked_match=expected_masked == observed_masked,
        tbl_diff_preview=[] if expected_tbl == observed_tbl else preview_diff(expected_tbl, observed_tbl, "tbl"),
        masked_diff_preview=[]
        if expected_masked == observed_masked
        else preview_diff(expected_masked, observed_masked, "masked"),
        return_code=return_code,
        stdout=stdout,
        stderr=stderr,
    )


def suite_summaries(results: list[CaseResult]) -> dict[str, dict[str, int]]:
    suites = sorted({result.suite for result in results})
    return {
        suite: {
            "total_cases": len([result for result in results if result.suite == suite]),
            "passing_cases": sum(
                1
                for result in results
                if result.suite == suite
                and result.return_code == 0
                and result.tbl_match
                and result.masked_match
            ),
        }
        for suite in suites
    }


def write_report(report_dir: pathlib.Path, results: list[CaseResult]) -> tuple[pathlib.Path, pathlib.Path]:
    report_dir.mkdir(parents=True, exist_ok=True)
    json_path = report_dir / "latest_report.json"
    md_path = report_dir / "latest_report.md"

    aggregate = {
        "total_cases": len(results),
        "passing_cases": sum(
            1
            for result in results
            if result.return_code == 0 and result.tbl_match and result.masked_match
        ),
        "failing_cases": [
            result.name
            for result in results
            if result.return_code != 0 or not result.tbl_match or not result.masked_match
        ],
        "median_runtime_ms_across_cases": statistics.median(
            [result.median_runtime_ms for result in results]
        )
        if results
        else float("nan"),
        "suite_summaries": suite_summaries(results),
    }
    payload = {
        "aggregate": aggregate,
        "results": [dataclasses.asdict(result) for result in results],
    }
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# RepeatMasker Fixture Benchmark",
        "",
        f"- total cases: `{aggregate['total_cases']}`",
        f"- passing cases: `{aggregate['passing_cases']}`",
        f"- median runtime across cases: `{aggregate['median_runtime_ms_across_cases']:.3f} ms`",
        "",
    ]
    for suite, summary in aggregate["suite_summaries"].items():
        lines.append(f"- `{suite}`: `{summary['passing_cases']}/{summary['total_cases']}` cases matched")
    lines.extend(
        [
            "",
            "| case | suite | mode | bases | median ms | tbl | masked |",
            "| --- | --- | --- | ---: | ---: | --- | --- |",
        ]
    )
    for result in results:
        lines.append(
            f"| {result.name} | `{result.suite}` | `{result.mask_mode}` | {result.total_bases} | "
            f"{result.median_runtime_ms:.3f} | "
            f"{'ok' if result.tbl_match else 'diff'} | "
            f"{'ok' if result.masked_match else 'diff'} |"
        )
    failing = [result for result in results if not (result.tbl_match and result.masked_match)]
    if failing:
        lines.extend(["", "## Mismatches", ""])
        for result in failing:
            lines.append(f"### {result.name}")
            if result.tbl_diff_preview:
                lines.extend(["```diff", *result.tbl_diff_preview, "```", ""])
            if result.masked_diff_preview:
                lines.extend(["```diff", *result.masked_diff_preview, "```", ""])
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return json_path, md_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeat-count", type=int, default=3, help="Runs per fixture")
    parser.add_argument(
        "--report-dir",
        type=pathlib.Path,
        default=DEFAULT_REPORT_DIR,
        help="Output directory for JSON and Markdown reports",
    )
    parser.add_argument(
        "--build-release",
        action="store_true",
        help="Force a fresh release build of process_repeats_rs before benchmarking",
    )
    parser.add_argument(
        "--include-legacy",
        action="store_true",
        help="Also benchmark legacy `.cat` and old-format `.out` fixtures",
    )
    args = parser.parse_args()

    binary = ensure_release_binary(force_build=args.build_release)
    report_dir = args.report_dir
    report_dir.mkdir(parents=True, exist_ok=True)
    results = [
        run_case(binary, case, args.repeat_count, report_dir)
        for case in fixture_cases(include_legacy=args.include_legacy)
    ]
    json_path, md_path = write_report(report_dir, results)

    passing = sum(
        1
        for result in results
        if result.return_code == 0 and result.tbl_match and result.masked_match
    )
    print(
        f"Wrote RepeatMasker fixture report to {json_path} and {md_path} "
        f"({passing}/{len(results)} cases matched gold outputs)"
    )
    return 0 if passing == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
