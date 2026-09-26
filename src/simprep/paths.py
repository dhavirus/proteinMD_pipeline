"""Locations of the repository-level data directories (schemas and knowledge base).

simprep is used from a checkout (editable install or a pinned clone on Colab); the
schemas and rules live beside ``src/`` and are versioned with the code.
"""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_DIR = REPO_ROOT / "schema"
KNOWLEDGE_DIR = REPO_ROOT / "knowledge"


def require_dir(path: Path) -> Path:
    """Return ``path`` if it is a directory, else raise with a remedy."""
    if not path.is_dir():
        raise FileNotFoundError(
            f"{path} not found. simprep must run from a repository checkout "
            "(pip install -e .) so that schema/ and knowledge/ are available."
        )
    return path
