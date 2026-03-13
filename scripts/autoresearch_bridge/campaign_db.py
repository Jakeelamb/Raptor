#!/usr/bin/env python3
"""SQLite-backed state tracking for autonomous autoresearch campaigns."""

from __future__ import annotations

import json
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


CAMPAIGN_TERMINAL_STATUSES = {"succeeded", "failed", "cancelled"}
STAGE_TERMINAL_STATUSES = {"succeeded", "failed", "skipped"}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def encode_json(value: Any | None) -> str | None:
    if value is None:
        return None
    return json.dumps(value, indent=2, sort_keys=True)


def decode_json(raw: str | None) -> Any | None:
    if raw is None or raw == "":
        return None
    return json.loads(raw)


@dataclass(frozen=True)
class CampaignRecord:
    id: int
    name: str
    status: str
    manifest_path: str
    campaign_root: str
    max_workers: int
    fail_fast: bool
    metadata: dict[str, Any] | None
    last_error: str | None
    created_at: str
    updated_at: str
    finished_at: str | None


class CampaignDb:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(self.path)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute("PRAGMA foreign_keys = ON")
        self._init_schema()

    def close(self) -> None:
        self.conn.close()

    def _init_schema(self) -> None:
        self.conn.executescript(
            """
            CREATE TABLE IF NOT EXISTS campaigns (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT NOT NULL,
                status TEXT NOT NULL,
                manifest_path TEXT NOT NULL,
                campaign_root TEXT NOT NULL,
                max_workers INTEGER NOT NULL,
                fail_fast INTEGER NOT NULL,
                metadata_json TEXT,
                last_error TEXT,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL,
                finished_at TEXT
            );

            CREATE TABLE IF NOT EXISTS component_runs (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                campaign_id INTEGER NOT NULL REFERENCES campaigns(id) ON DELETE CASCADE,
                component_name TEXT NOT NULL,
                status TEXT NOT NULL,
                artifact_root TEXT NOT NULL,
                tasks_root TEXT NOT NULL,
                manifest_json TEXT NOT NULL,
                latest_prepare_metadata_path TEXT,
                latest_train_result_path TEXT,
                latest_promote_summary_path TEXT,
                last_error TEXT,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL,
                UNIQUE(campaign_id, component_name)
            );

            CREATE TABLE IF NOT EXISTS stage_runs (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                component_run_id INTEGER NOT NULL REFERENCES component_runs(id) ON DELETE CASCADE,
                stage_name TEXT NOT NULL,
                attempt INTEGER NOT NULL,
                status TEXT NOT NULL,
                command_json TEXT NOT NULL,
                workdir TEXT NOT NULL,
                stdout_path TEXT NOT NULL,
                stderr_path TEXT NOT NULL,
                result_path TEXT,
                result_json TEXT,
                exit_code INTEGER,
                duration_seconds REAL,
                error_message TEXT,
                started_at TEXT NOT NULL,
                finished_at TEXT,
                UNIQUE(component_run_id, stage_name, attempt)
            );

            CREATE TABLE IF NOT EXISTS artifacts (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                stage_run_id INTEGER NOT NULL REFERENCES stage_runs(id) ON DELETE CASCADE,
                artifact_type TEXT NOT NULL,
                path TEXT NOT NULL,
                metadata_json TEXT,
                created_at TEXT NOT NULL,
                UNIQUE(stage_run_id, artifact_type, path)
            );
            """
        )
        self.conn.commit()

    def create_campaign(
        self,
        *,
        name: str,
        manifest_path: Path,
        campaign_root: Path,
        max_workers: int,
        fail_fast: bool,
        metadata: dict[str, Any] | None,
    ) -> int:
        now = utc_now()
        cursor = self.conn.execute(
            """
            INSERT INTO campaigns (
                name,
                status,
                manifest_path,
                campaign_root,
                max_workers,
                fail_fast,
                metadata_json,
                created_at,
                updated_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                name,
                "running",
                str(manifest_path),
                str(campaign_root),
                max_workers,
                int(fail_fast),
                encode_json(metadata),
                now,
                now,
            ),
        )
        self.conn.commit()
        return int(cursor.lastrowid)

    def get_campaign(self, campaign_id: int) -> CampaignRecord:
        row = self.conn.execute("SELECT * FROM campaigns WHERE id = ?", (campaign_id,)).fetchone()
        if row is None:
            raise KeyError(f"unknown campaign id {campaign_id}")
        return CampaignRecord(
            id=int(row["id"]),
            name=str(row["name"]),
            status=str(row["status"]),
            manifest_path=str(row["manifest_path"]),
            campaign_root=str(row["campaign_root"]),
            max_workers=int(row["max_workers"]),
            fail_fast=bool(row["fail_fast"]),
            metadata=decode_json(row["metadata_json"]),
            last_error=row["last_error"],
            created_at=str(row["created_at"]),
            updated_at=str(row["updated_at"]),
            finished_at=row["finished_at"],
        )

    def latest_campaign_id(self) -> int | None:
        row = self.conn.execute("SELECT id FROM campaigns ORDER BY id DESC LIMIT 1").fetchone()
        if row is None:
            return None
        return int(row["id"])

    def list_campaigns(self, limit: int = 20) -> list[sqlite3.Row]:
        cursor = self.conn.execute(
            "SELECT * FROM campaigns ORDER BY id DESC LIMIT ?",
            (limit,),
        )
        return list(cursor.fetchall())

    def update_campaign_status(
        self,
        campaign_id: int,
        status: str,
        *,
        last_error: str | None = None,
    ) -> None:
        now = utc_now()
        finished_at = now if status in CAMPAIGN_TERMINAL_STATUSES else None
        self.conn.execute(
            """
            UPDATE campaigns
            SET status = ?,
                last_error = ?,
                updated_at = ?,
                finished_at = ?
            WHERE id = ?
            """,
            (status, last_error, now, finished_at, campaign_id),
        )
        self.conn.commit()

    def create_component_run(
        self,
        *,
        campaign_id: int,
        component_name: str,
        artifact_root: Path,
        tasks_root: Path,
        manifest: dict[str, Any],
    ) -> int:
        existing = self.get_component_run(campaign_id, component_name)
        if existing is not None:
            return int(existing["id"])

        now = utc_now()
        cursor = self.conn.execute(
            """
            INSERT INTO component_runs (
                campaign_id,
                component_name,
                status,
                artifact_root,
                tasks_root,
                manifest_json,
                created_at,
                updated_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                campaign_id,
                component_name,
                "pending",
                str(artifact_root),
                str(tasks_root),
                encode_json(manifest),
                now,
                now,
            ),
        )
        self.conn.commit()
        return int(cursor.lastrowid)

    def get_component_run(self, campaign_id: int, component_name: str) -> sqlite3.Row | None:
        return self.conn.execute(
            """
            SELECT * FROM component_runs
            WHERE campaign_id = ? AND component_name = ?
            """,
            (campaign_id, component_name),
        ).fetchone()

    def list_component_runs(self, campaign_id: int) -> list[sqlite3.Row]:
        cursor = self.conn.execute(
            """
            SELECT * FROM component_runs
            WHERE campaign_id = ?
            ORDER BY component_name
            """,
            (campaign_id,),
        )
        return list(cursor.fetchall())

    def update_component_run(
        self,
        component_run_id: int,
        *,
        status: str | None = None,
        latest_prepare_metadata_path: Path | None = None,
        latest_train_result_path: Path | None = None,
        latest_promote_summary_path: Path | None = None,
        last_error: str | None = None,
    ) -> None:
        row = self.conn.execute(
            "SELECT * FROM component_runs WHERE id = ?",
            (component_run_id,),
        ).fetchone()
        if row is None:
            raise KeyError(f"unknown component run id {component_run_id}")
        now = utc_now()
        self.conn.execute(
            """
            UPDATE component_runs
            SET status = COALESCE(?, status),
                latest_prepare_metadata_path = COALESCE(?, latest_prepare_metadata_path),
                latest_train_result_path = COALESCE(?, latest_train_result_path),
                latest_promote_summary_path = COALESCE(?, latest_promote_summary_path),
                last_error = ?,
                updated_at = ?
            WHERE id = ?
            """,
            (
                status,
                str(latest_prepare_metadata_path) if latest_prepare_metadata_path else None,
                str(latest_train_result_path) if latest_train_result_path else None,
                str(latest_promote_summary_path) if latest_promote_summary_path else None,
                last_error,
                now,
                component_run_id,
            ),
        )
        self.conn.commit()

    def next_stage_attempt(self, component_run_id: int, stage_name: str) -> int:
        row = self.conn.execute(
            """
            SELECT MAX(attempt) AS max_attempt
            FROM stage_runs
            WHERE component_run_id = ? AND stage_name = ?
            """,
            (component_run_id, stage_name),
        ).fetchone()
        max_attempt = row["max_attempt"] if row is not None else None
        return 1 if max_attempt is None else int(max_attempt) + 1

    def start_stage_run(
        self,
        *,
        component_run_id: int,
        stage_name: str,
        attempt: int,
        command: list[str],
        workdir: Path,
        stdout_path: Path,
        stderr_path: Path,
    ) -> int:
        now = utc_now()
        cursor = self.conn.execute(
            """
            INSERT INTO stage_runs (
                component_run_id,
                stage_name,
                attempt,
                status,
                command_json,
                workdir,
                stdout_path,
                stderr_path,
                started_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                component_run_id,
                stage_name,
                attempt,
                "running",
                encode_json(command),
                str(workdir),
                str(stdout_path),
                str(stderr_path),
                now,
            ),
        )
        self.conn.execute(
            """
            UPDATE component_runs
            SET status = ?, updated_at = ?, last_error = NULL
            WHERE id = ?
            """,
            ("running", now, component_run_id),
        )
        self.conn.commit()
        return int(cursor.lastrowid)

    def finish_stage_run(
        self,
        stage_run_id: int,
        *,
        status: str,
        exit_code: int,
        duration_seconds: float,
        result_path: Path | None,
        result_json: dict[str, Any] | None,
        error_message: str | None,
    ) -> None:
        now = utc_now()
        self.conn.execute(
            """
            UPDATE stage_runs
            SET status = ?,
                exit_code = ?,
                duration_seconds = ?,
                result_path = ?,
                result_json = ?,
                error_message = ?,
                finished_at = ?
            WHERE id = ?
            """,
            (
                status,
                exit_code,
                duration_seconds,
                str(result_path) if result_path else None,
                encode_json(result_json),
                error_message,
                now,
                stage_run_id,
            ),
        )
        self.conn.commit()

    def latest_stage_run(self, component_run_id: int, stage_name: str) -> sqlite3.Row | None:
        return self.conn.execute(
            """
            SELECT *
            FROM stage_runs
            WHERE component_run_id = ? AND stage_name = ?
            ORDER BY attempt DESC
            LIMIT 1
            """,
            (component_run_id, stage_name),
        ).fetchone()

    def list_stage_runs(self, component_run_id: int) -> list[sqlite3.Row]:
        cursor = self.conn.execute(
            """
            SELECT *
            FROM stage_runs
            WHERE component_run_id = ?
            ORDER BY stage_name, attempt
            """,
            (component_run_id,),
        )
        return list(cursor.fetchall())

    def add_artifact(
        self,
        *,
        stage_run_id: int,
        artifact_type: str,
        path: Path,
        metadata: dict[str, Any] | None = None,
    ) -> None:
        self.conn.execute(
            """
            INSERT OR REPLACE INTO artifacts (
                stage_run_id,
                artifact_type,
                path,
                metadata_json,
                created_at
            ) VALUES (?, ?, ?, ?, ?)
            """,
            (
                stage_run_id,
                artifact_type,
                str(path),
                encode_json(metadata),
                utc_now(),
            ),
        )
        self.conn.commit()

    def list_artifacts(self, stage_run_id: int) -> list[sqlite3.Row]:
        cursor = self.conn.execute(
            """
            SELECT *
            FROM artifacts
            WHERE stage_run_id = ?
            ORDER BY artifact_type, path
            """,
            (stage_run_id,),
        )
        return list(cursor.fetchall())

    def mark_running_stages_failed(self, campaign_id: int, *, error_message: str) -> int:
        rows = self.conn.execute(
            """
            SELECT sr.id AS stage_run_id, sr.component_run_id
            FROM stage_runs sr
            INNER JOIN component_runs cr ON cr.id = sr.component_run_id
            WHERE cr.campaign_id = ? AND sr.status = 'running'
            """,
            (campaign_id,),
        ).fetchall()
        if not rows:
            return 0

        now = utc_now()
        self.conn.execute(
            """
            UPDATE stage_runs
            SET status = 'failed',
                exit_code = COALESCE(exit_code, -1),
                duration_seconds = COALESCE(duration_seconds, 0.0),
                error_message = ?,
                finished_at = ?
            WHERE id IN (
                SELECT sr.id
                FROM stage_runs sr
                INNER JOIN component_runs cr ON cr.id = sr.component_run_id
                WHERE cr.campaign_id = ? AND sr.status = 'running'
            )
            """,
            (error_message, now, campaign_id),
        )
        self.conn.execute(
            """
            UPDATE component_runs
            SET status = 'failed',
                last_error = ?,
                updated_at = ?
            WHERE campaign_id = ? AND status = 'running'
            """,
            (error_message, now, campaign_id),
        )
        self.conn.commit()
        return len(rows)
