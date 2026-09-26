"""Colab variants module (TASK-005 decision 6), run locally with temp dirs for Drive."""

import json
from datetime import UTC, datetime
from pathlib import Path

import pytest
import yaml

from modules.variants_module import run_variants_job
from simprep.knowledge import load_ruleset
from simprep.prep.run import PrepRequest
from simprep.variants.run import run_variants
from tests.prep_helpers import decided_manifest

PANEL_DIR = Path(__file__).parent / "panel"
FIXTURE = yaml.safe_load((PANEL_DIR / "variant_fixtures" / "5FQL.yaml").read_text())
PREP = yaml.safe_load((PANEL_DIR / "prep_fixtures" / FIXTURE["prep_fixture"]).read_text())
STRUCTURE = PANEL_DIR / FIXTURE["structure"]
NOW = datetime(2026, 9, 26, 15, 0, 0, tzinfo=UTC)


def config(tmp_path, manifest):
    return {
        "repo_commit": "0" * 40,
        "drive_root": str(tmp_path / "drive"),
        "work_root": str(tmp_path / "work"),
        "structure": str(STRUCTURE),
        "manifest": str(manifest),
        "run_dir": None,
        "seed": 3,
    }


@pytest.fixture
def undecided(tmp_path):
    decisions = PREP["decisions"] | FIXTURE["prep_overrides"]
    path = decided_manifest(STRUCTURE, decisions, tmp_path / "m.json")
    document = json.loads(path.read_text())
    document["variants"] = [v for v in FIXTURE["variants"] if v["name"] == "R468Q"]
    document["relaxation"] = {**load_ruleset().relaxation, "enabled": False}  # speed
    path.write_text(json.dumps(document))
    return path


def decided(undecided, tmp_path):
    """Pass 1 locally, then the recommended rotamer as a recorded decision."""
    run_variants(PrepRequest(STRUCTURE, undecided, tmp_path / "pass1"))
    document = json.loads((tmp_path / "pass1" / "manifest.json").read_text())
    (finding,) = document["variant_snapshot"]["findings"]
    rotamer = next(e["value"] for e in finding["evidence"] if e["key"] == "recommended_rotamer")
    document["decisions"].append(
        {
            "finding_id": finding["id"],
            "option_id": "choose_rotamer",
            "parameters": {"rotamer": rotamer},
            "rationale": "test",
            "decided_by": "test",
            "timestamp": "2026-09-26T00:00:00Z",
        }
    )
    path = tmp_path / "decided.json"
    path.write_text(json.dumps(document))
    return path


def test_first_run_stops_with_the_candidates_in_the_run_directory(tmp_path, undecided):
    result = run_variants_job(config(tmp_path, undecided), NOW)
    assert result.status == "needs_decisions"
    assert "variant_snapshot" in json.loads((result.run_dir / "manifest.json").read_text())


def test_decided_run_builds_and_resumes(tmp_path, undecided):
    settings = config(tmp_path, decided(undecided, tmp_path))
    result = run_variants_job(settings, NOW)
    assert result.status == "complete" and not result.skipped
    for path in ("wt/system.cif", "R468Q/system.cif", "variant_record.json", "config.json"):
        assert (result.run_dir / path).is_file(), path
    resume = settings | {"run_dir": str(result.run_dir)}
    assert run_variants_job(resume).skipped
    (result.run_dir / "R468Q" / "system.cif").write_text("damaged")
    again = run_variants_job(resume)
    assert not again.skipped and again.status == "complete"
