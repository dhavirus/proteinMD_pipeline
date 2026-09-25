"""Regression tests on the committed panel structures.

Expected findings come from tests/panel/expected_findings/*.yaml, derived from each
file's own annotations (derive_expected.py) plus hand-written publication checks, never
from detector output. Any mismatch fails CI.
"""

from functools import cache
from pathlib import Path

import pytest
import yaml

from simprep.audit import ReportContext, build_findings_report, run_audit
from simprep.config import default_config
from simprep.knowledge import load_ruleset
from simprep.schemas import validation_errors
from simprep.structure.parse import read_structure

PANEL_DIR = Path(__file__).parent / "panel"
EXPECTED_DIR = PANEL_DIR / "expected_findings"
PANEL_SIZE = 5
FAMILIES = ("metals", "nonstandard_residues", "altlocs", "missing_residues", "unrecognized")
EXPECTED = sorted(EXPECTED_DIR.glob("*.yaml"))


def expected_param(path):
    return pytest.param(yaml.safe_load(path.read_text()), id=path.stem)


@cache
def audit(file_name: str) -> dict:
    """findings.json for a panel file (default config: no regions)."""
    ruleset = load_ruleset()
    path = PANEL_DIR / file_name
    structure = read_structure(path)
    findings = run_audit(structure, ruleset, default_config(ruleset.audit_defaults))
    context = ReportContext(
        {"path": f"tests/panel/{file_name}", "sha256": "0" * 64, "format": "mmcif"},
        ruleset,
        default_config(ruleset.audit_defaults),
        {"version": "test", "git_commit": None},
        "2026-01-01T00:00:00Z",
    )
    return build_findings_report(structure, findings, context)


def by_id(report: dict) -> dict[str, dict]:
    return {finding["id"]: finding for finding in report["findings"]}


def test_panel_is_complete():
    assert len(EXPECTED) == PANEL_SIZE, f"expected {PANEL_SIZE} panel entries, got {EXPECTED}"
    for path in EXPECTED:
        entry = yaml.safe_load(path.read_text())
        assert (PANEL_DIR / entry["file"]).is_file(), entry["file"]
        assert entry["reviewed"] is False or entry["reviewed"] is True


@pytest.mark.parametrize("expected", [expected_param(p) for p in EXPECTED])
def test_findings_json_is_schema_valid(expected):
    assert validation_errors(audit(expected["file"]), "findings_report") == []


@pytest.mark.parametrize("family", FAMILIES)
@pytest.mark.parametrize("expected", [expected_param(p) for p in EXPECTED])
def test_family_findings_match_annotations(expected, family):
    actual = {f["id"]: f for f in audit(expected["file"])["findings"] if f["rule_family"] == family}
    wanted = {entry["id"]: entry for entry in expected["families"][family]}
    assert sorted(actual) == sorted(wanted)
    for finding_id, entry in wanted.items():
        if "rule_id" in entry:
            assert actual[finding_id]["rule_id"] == entry["rule_id"], finding_id


def shell_atoms(finding: dict) -> set[tuple[str, str]]:
    return {
        (f"{a['chain']}:{a['seq_num']}{a['ins_code']}", a["atom_name"])
        for item in finding["evidence"]
        if item["key"] == "ligand_distance"
        for a in item["atoms"][1:]
    }


@pytest.mark.parametrize("expected", [expected_param(p) for p in EXPECTED])
def test_annotated_metal_ligands_within_cutoff_are_in_shell(expected):
    findings = by_id(audit(expected["file"]))
    for entry in expected["families"]["metals"]:
        finding = findings[entry["id"]]
        cutoff = next(i["value"] for i in finding["evidence"] if i["key"] == "coordination_cutoff")
        within = {
            (lig["residue"], lig["atom"])
            for lig in entry["annotated_ligands"]
            if lig["distance_angstrom"] <= cutoff
        }
        assert within <= shell_atoms(finding), entry["id"]


@pytest.mark.parametrize("expected", [expected_param(p) for p in EXPECTED])
def test_record_accounting_balances(expected):
    for row in audit(expected["file"])["counts"]:
        assert row["in"] == row["passed_through"] + row["flagged"] + row["excluded"], row


PUBLICATION = [
    pytest.param(e, id=p.stem)
    for p in EXPECTED
    if "publication_checks" in (e := yaml.safe_load(p.read_text()))
]


@pytest.mark.parametrize("expected", PUBLICATION)
def test_publication_checks(expected):
    findings = by_id(audit(expected["file"]))
    for check in expected["publication_checks"]:
        if "count" in check:
            prefix = check["count"]["id_prefix"]
            matching = [f for f in findings if f.startswith(prefix)]
            assert len(matching) >= check["count"]["min"], check["source"]
            continue
        finding = findings[check["finding_id"]]
        if "shell_residues" in check:
            residues = {residue for residue, _ in shell_atoms(finding)}
            assert set(check["shell_residues"]) <= residues, check["source"]
        for key, value in check.get("evidence", {}).items():
            actual = next(i.get("value") for i in finding["evidence"] if i["key"] == key)
            assert actual == value, (check["finding_id"], key)
