#!/usr/bin/env python3
"""Shared helpers for robust autoresearch bridge subprocess handling."""

from __future__ import annotations

import json
import shlex
import subprocess
from pathlib import Path
from typing import Any, Sequence


def extract_last_json_document(text: str) -> Any:
    decoder = json.JSONDecoder()
    best_value: Any | None = None
    best_end = -1
    search_start = 0

    while True:
        candidate_start = text.find("{", search_start)
        if candidate_start == -1:
            break
        try:
            value, end = decoder.raw_decode(text, candidate_start)
        except json.JSONDecodeError:
            search_start = candidate_start + 1
            continue
        if end > best_end:
            best_value = value
            best_end = end
        search_start = candidate_start + 1

    if best_value is None:
        raise json.JSONDecodeError("no JSON document found in subprocess output", text, 0)
    return best_value


def format_command(cmd: Sequence[str]) -> str:
    return " ".join(shlex.quote(part) for part in cmd)


def run_json_command(cmd: Sequence[str], cwd: Path | None = None) -> Any:
    completed = subprocess.run(
        list(cmd),
        check=True,
        capture_output=True,
        text=True,
        cwd=str(cwd) if cwd is not None else None,
    )
    for stream in (completed.stdout, completed.stderr, f"{completed.stdout}\n{completed.stderr}"):
        if not stream.strip():
            continue
        try:
            return extract_last_json_document(stream)
        except json.JSONDecodeError:
            continue
    raise json.JSONDecodeError(
        "failed to recover JSON document from subprocess output",
        f"stdout:\n{completed.stdout}\n\nstderr:\n{completed.stderr}",
        0,
    )
