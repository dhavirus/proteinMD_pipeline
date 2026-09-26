"""Provenance at the I/O edge: file hashes, simprep version and git commit, timestamps."""

from __future__ import annotations

import hashlib
import subprocess
from datetime import UTC, datetime
from pathlib import Path

from simprep import __version__
from simprep.paths import REPO_ROOT

HASH_CHUNK_BYTES = 1 << 20


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(HASH_CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_commit(repo: Path = REPO_ROOT) -> str | None:
    """HEAD commit of the simprep checkout, suffixed ``-dirty`` if tracked files changed.

    Returns None when simprep does not run from a git checkout (recorded as such).
    """
    try:
        head = _git(repo, "rev-parse", "HEAD")
        dirty = _git(repo, "status", "--porcelain", "--untracked-files=no")
    except (OSError, subprocess.CalledProcessError):
        return None
    return f"{head}-dirty" if dirty else head


def _git(repo: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", "-C", str(repo), *args], capture_output=True, text=True, check=True
    )
    return result.stdout.strip()


def simprep_provenance() -> dict:
    return {"version": __version__, "git_commit": git_commit()}


def utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")
