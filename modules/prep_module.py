"""Colab back end for prep (TASK-004): run bookkeeping around :func:`simprep.prep.run.run_prep`.

Layout (CLAUDE.md, Colab conventions): ``<drive_root>/<yymm>/run_<yymmdd_hhmmss>/`` holds
``config.json`` (the config, the manifest, the simprep commit and the ``run_id``) and
the prep outputs. Work happens under ``work_root`` (local disk on Colab) and finished
outputs are copied to the run directory in one go. Re-running with ``run_dir`` set to an
existing run resumes it: if the recorded outputs are all present with matching SHA-256
the run is skipped, otherwise prep runs again and the outputs are replaced.

Locally (tests) the same function runs with temporary directories in place of Drive.
"""

from __future__ import annotations

import json
import random
import shutil
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from simprep.manifest import write_json
from simprep.prep.run import RECORD_FILE, PrepRequest, run_prep
from simprep.provenance import sha256_file, simprep_provenance

CONFIG_FILE = "config.json"
REQUIRED_KEYS = ("drive_root", "work_root", "structure", "manifest", "seed")


class RunConfigError(ValueError):
    """The config dict cannot drive a run, or disagrees with the run it resumes."""


@dataclass(frozen=True)
class RunResult:
    run_dir: Path
    run_id: str
    skipped: bool


def run_prep_job(config: dict, now: datetime | None = None) -> RunResult:
    """Run (or resume) one prep run described by the notebook's config dict."""
    missing = [key for key in REQUIRED_KEYS if key not in config]
    if missing:
        raise RunConfigError(f"config is missing {', '.join(missing)}")
    random.seed(config["seed"])
    run_dir = Path(config["run_dir"]) if config.get("run_dir") else new_run_dir(config, now)
    run_id = run_dir.name
    _write_or_check_config(run_dir, config)
    if outputs_complete(run_dir):
        return RunResult(run_dir, run_id, skipped=True)
    work_dir = Path(config["work_root"]) / run_id
    shutil.rmtree(work_dir, ignore_errors=True)
    run_prep(PrepRequest(Path(config["structure"]), Path(config["manifest"]), work_dir))
    shutil.copytree(work_dir, run_dir, dirs_exist_ok=True)
    return RunResult(run_dir, run_id, skipped=False)


def new_run_dir(config: dict, now: datetime | None = None) -> Path:
    """``<drive_root>/<yymm>/run_<yymmdd_hhmmss>`` for ``now`` (UTC)."""
    moment = now or datetime.now(UTC)
    return (
        Path(config["drive_root"]) / moment.strftime("%y%m") / moment.strftime("run_%y%m%d_%H%M%S")
    )


def _write_or_check_config(run_dir: Path, config: dict) -> None:
    """Write config.json for a new run; on resume, refuse a config that changed."""
    document = {
        "run_id": run_dir.name,
        "config": {key: config[key] for key in sorted(config) if key != "run_dir"},
        "manifest": json.loads(Path(config["manifest"]).read_text()),
        "input_sha256": sha256_file(Path(config["structure"])),
    }
    path = run_dir / CONFIG_FILE
    if path.exists():
        recorded = json.loads(path.read_text())
        if {k: recorded[k] for k in document} != document:
            raise RunConfigError(
                f"{path} records a different config, manifest or input; start a new run "
                "(leave run_dir unset) instead of resuming this one"
            )
        return
    run_dir.mkdir(parents=True, exist_ok=True)
    write_json(document | {"simprep": simprep_provenance()}, path)


def outputs_complete(run_dir: Path) -> bool:
    """True when prep_record.json exists and every file it lists is present and matches."""
    record_path = run_dir / RECORD_FILE
    if not record_path.exists():
        return False
    record = json.loads(record_path.read_text())
    files = [f for system in record["systems"] for f in system["files"]]
    return all(
        (run_dir / f["path"]).exists() and sha256_file(run_dir / f["path"]) == f["sha256"]
        for f in files
    )
