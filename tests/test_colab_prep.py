"""Colab prep module and notebook (TASK-004 §5), run locally with temp dirs for Drive."""

import json
from datetime import UTC, datetime
from pathlib import Path

import nbformat
import pytest
import yaml

from modules.prep_module import RunConfigError, run_prep_job
from tests.prep_helpers import decided_manifest

REPO_ROOT = Path(__file__).parent.parent
PANEL_DIR = REPO_ROOT / "tests" / "panel"
NOTEBOOK = REPO_ROOT / "notebooks" / "01_prep.ipynb"
NOW = datetime(2026, 9, 26, 14, 5, 9, tzinfo=UTC)


@pytest.fixture
def config(tmp_path):
    fixture = yaml.safe_load((PANEL_DIR / "prep_fixtures" / "6OIM.yaml").read_text())
    structure = PANEL_DIR / "6OIM.cif.gz"
    manifest = decided_manifest(structure, fixture["decisions"], tmp_path / "manifest.json")
    return {
        "repo_commit": "0" * 40,
        "drive_root": str(tmp_path / "drive"),
        "work_root": str(tmp_path / "work"),
        "structure": str(structure),
        "manifest": str(manifest),
        "run_dir": None,
        "seed": 7,
    }


def test_new_run_follows_the_drive_layout_and_records_its_config(config):
    result = run_prep_job(config, NOW)
    assert result.run_dir == Path(config["drive_root"]) / "2609" / "run_260926_140509"
    assert not result.skipped
    recorded = json.loads((result.run_dir / "config.json").read_text())
    assert recorded["run_id"] == "run_260926_140509"
    assert recorded["config"]["seed"] == 7
    assert recorded["manifest"] == json.loads(Path(config["manifest"]).read_text())
    assert "git_commit" in recorded["simprep"]
    assert (result.run_dir / "system.cif").is_file()
    assert (result.run_dir / "prep_record.json").is_file()


def test_resume_skips_a_complete_run_and_redoes_a_damaged_one(config):
    first = run_prep_job(config, NOW)
    resume = config | {"run_dir": str(first.run_dir)}
    assert run_prep_job(resume).skipped
    cif = first.run_dir / "system.cif"
    original = cif.read_bytes()
    cif.write_text("truncated")
    assert not run_prep_job(resume).skipped
    assert cif.read_bytes() == original


def test_resume_refuses_a_changed_config(config):
    first = run_prep_job(config, NOW)
    with pytest.raises(RunConfigError, match="different config"):
        run_prep_job(config | {"run_dir": str(first.run_dir), "seed": 8})


def test_notebook_is_valid_with_one_config_cell_and_no_autoreload():
    notebook = nbformat.read(NOTEBOOK, as_version=4)
    nbformat.validate(notebook)
    code = [cell.source for cell in notebook.cells if cell.cell_type == "code"]
    assert code[0].count("CONFIG = {") == 1
    assert all("CONFIG = {" not in source for source in code[1:])
    assert "importlib.reload(prep_module)" in "\n".join(code)
    assert "autoreload" not in "\n".join(code)
    assert all("\n" in source for source in code), "cells must keep their newlines"
