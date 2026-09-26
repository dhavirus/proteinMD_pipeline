"""Cross-language contract: manifests built by the front end's JavaScript modules must pass
simprep's schema, snapshot-hash and decision checks, and drive `simprep` as intended.

Runs the modules under Node (no browser, no network). Node is part of CI (TASK-002).
"""

import json
import shutil
import subprocess
import sys

import pytest

from simprep.cli import main
from simprep.manifest import load_manifest
from simprep.paths import REPO_ROOT
from simprep.schemas import validation_errors

FRONTEND = REPO_ROOT / "frontend"
DEMO_FINDINGS = FRONTEND / "demo" / "5FQL.findings.json"
STRUCTURE = REPO_ROOT / "tests" / "panel" / "5FQL.cif.gz"
ACTIVE_SITE = "A:45 A:46 A:84 A:334-335 A:1551"
NODE = shutil.which("node")

pytestmark = pytest.mark.skipif(
    NODE is None, reason="Node.js is required for front-end contract tests"
)


@pytest.fixture(scope="module")
def exported(tmp_path_factory):
    out = tmp_path_factory.mktemp("manifests")
    script = FRONTEND / "tests" / "export_manifests.mjs"
    subprocess.run([NODE, str(script), str(DEMO_FINDINGS), str(out), ACTIVE_SITE], check=True)
    return out


def cli(*args) -> int:
    return main([str(a) for a in args])


def test_every_exported_manifest_passes_python_checks(exported):
    paths = sorted(exported.glob("*.json"))
    assert [p.stem for p in paths] == [
        "blocking_decided",
        "explicit_choice",
        "region_all_decided",
        "undecided",
    ]
    for path in paths:
        manifest = json.loads(path.read_text())
        assert validation_errors(manifest, "manifest") == [], path.name
        load_manifest(path)  # snapshot hash (RFC 8785) and decisions


def test_status_gate_follows_decisions(exported, capsys):
    assert cli("manifest", "status", exported / "undecided.json", "--structure", STRUCTURE) == 3
    assert (
        cli("manifest", "status", exported / "blocking_decided.json", "--structure", STRUCTURE) == 0
    )
    assert cli("manifest", "status", exported / "explicit_choice.json") == 0
    assert "explicit-choice option" in capsys.readouterr().out


def test_region_from_front_end_escalates_on_reaudit(exported, tmp_path):
    out = tmp_path / "reaudit"
    assert (
        cli("audit", STRUCTURE, "--manifest", exported / "region_all_decided.json", "--out", out)
        == 0
    )
    findings = {f["id"]: f for f in json.loads((out / "findings.json").read_text())["findings"]}
    metal = findings["metals/A:1551"]
    assert (metal["base_severity"], metal["effective_severity"]) == ("warn", "blocking")
    assert metal["recommended_option"] == "bonded_mcpb"
    updated = load_manifest(out / "manifest.json")
    assert len(updated["decisions"]) == len(findings)


def test_demo_findings_are_current():
    result = subprocess.run(
        [sys.executable, str(FRONTEND / "demo" / "build_demo.py"), "--check"],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    assert result.returncode == 0, result.stderr
