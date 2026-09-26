"""`simprep prep` end to end on a small panel structure (6OIM): outputs, gate, input check."""

import json
from pathlib import Path

import yaml

from simprep.cli import main
from tests.prep_helpers import decided_manifest

PANEL_DIR = Path(__file__).parent / "panel"
STRUCTURE = PANEL_DIR / "6OIM.cif.gz"
BLOCKING_GROUP = "unrecognized/group/A:303"


def decisions() -> dict:
    return yaml.safe_load((PANEL_DIR / "prep_fixtures" / "6OIM.yaml").read_text())["decisions"]


def run_cli(*args) -> int:
    return main([str(a) for a in args])


def test_prep_writes_system_files_record_and_report(tmp_path):
    manifest = decided_manifest(STRUCTURE, decisions(), tmp_path / "m.json")
    out = tmp_path / "out"
    assert run_cli("prep", STRUCTURE, "--manifest", manifest, "--out", out) == 0
    assert sorted(p.name for p in out.iterdir()) == [
        "prep_record.json",
        "prep_report.md",
        "system.cif",
        "system.pdb",
    ]
    record = json.loads((out / "prep_record.json").read_text())
    assert record["input"]["sha256"] == json.loads(manifest.read_text())["input"]["sha256"]
    assert "## Work order for later stages" in (out / "prep_report.md").read_text()


def test_expert_review_blocks_prep_and_names_the_finding(tmp_path, capsys):
    chosen = decisions() | {BLOCKING_GROUP: {"option": "expert_review", "rationale": "open"}}
    manifest = decided_manifest(STRUCTURE, chosen, tmp_path / "m.json")
    out = tmp_path / "out"
    assert run_cli("prep", STRUCTURE, "--manifest", manifest, "--out", out) == 2
    assert f"{BLOCKING_GROUP}: expert_review is not a final decision" in capsys.readouterr().err
    assert not out.exists()
    assert run_cli("manifest", "status", manifest) == 3


def test_undecided_blocking_finding_blocks_prep(tmp_path, capsys):
    chosen = {k: v for k, v in decisions().items() if k != BLOCKING_GROUP}
    manifest = decided_manifest(STRUCTURE, chosen, tmp_path / "m.json")
    assert run_cli("prep", STRUCTURE, "--manifest", manifest, "--out", tmp_path / "out") == 2
    assert f"[blocking] {BLOCKING_GROUP}" in capsys.readouterr().err


def test_prep_refuses_a_different_structure(tmp_path, capsys):
    manifest = decided_manifest(STRUCTURE, decisions(), tmp_path / "m.json")
    other = PANEL_DIR / "3KS3.cif.gz"
    assert run_cli("prep", other, "--manifest", manifest, "--out", tmp_path / "out") == 2
    assert "manifest was built for input SHA-256" in capsys.readouterr().err
