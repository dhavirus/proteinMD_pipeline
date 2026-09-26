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
CURATED_FAMILIES = ("covalent_contacts",)
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


@pytest.mark.parametrize("family", CURATED_FAMILIES)
@pytest.mark.parametrize("expected", [expected_param(p) for p in EXPECTED])
def test_curated_family_findings(expected, family):
    """Families the annotations cannot give are listed by hand in curated_families."""
    actual = {f["id"]: f for f in audit(expected["file"])["findings"] if f["rule_family"] == family}
    wanted = {entry["id"]: entry for entry in expected["curated_families"][family]}
    assert sorted(actual) == sorted(wanted)
    for finding_id, entry in wanted.items():
        assert actual[finding_id]["rule_id"] == entry["rule_id"]


@pytest.mark.parametrize("expected", [expected_param(p) for p in EXPECTED])
def test_glycosylation_candidates_sit_on_derived_sequons(expected):
    """Independent check: n_glycosylation candidates start an N-X-S/T sequon according to
    the file's own entity sequence (derive_expected.py), and are recommended add_link."""
    sequon_residues = {r for chain in expected["sequon_asparagines"].values() for r in chain}
    findings = by_id(audit(expected["file"]))
    for entry in expected["curated_families"]["covalent_contacts"]:
        if entry["pattern"] != "n_glycosylation":
            continue
        assert entry["polymer_residue"] in sequon_residues, entry["id"]
        finding = findings[entry["id"]]
        evidence = {item["key"]: item.get("value") for item in finding["evidence"]}
        assert evidence["in_sequon"] is True and finding["recommended_option"] == "add_link"
        assert evidence["ligand_group_finding"].startswith("unrecognized/group/")


def shell_atoms(finding: dict, keys: tuple[str, ...] = ("ligand_distance",)) -> set[tuple]:
    return {
        (f"{a['chain']}:{a['seq_num']}{a['ins_code']}", a["atom_name"])
        for item in finding["evidence"]
        if item["key"] in keys
        for a in item["atoms"][1:]
    }


@pytest.mark.parametrize("expected", [expected_param(p) for p in EXPECTED])
def test_annotated_metal_ligands_within_cutoff_are_in_shell(expected):
    findings = by_id(audit(expected["file"]))
    for entry in expected["families"]["metals"]:
        finding = findings[entry["id"]]
        rule = next(r for r in load_ruleset().rules if r.rule_id == entry["rule_id"])
        cutoff = rule.matcher["coordination_cutoff_angstrom"]
        donors = set(rule.matcher["donor_elements"])
        within = {
            (lig["residue"], lig["atom"])
            for lig in entry["annotated_ligands"]
            if lig["distance_angstrom"] <= cutoff and lig["element"] in donors
        }
        # Annotated ligands appear in the shell, counted or excluded for partial occupancy.
        reported = shell_atoms(finding, ("ligand_distance", "excluded_partial_donor"))
        assert within <= reported, entry["id"]


@pytest.mark.parametrize("expected", [expected_param(p) for p in EXPECTED])
def test_record_accounting_balances(expected):
    for row in audit(expected["file"])["counts"]:
        assert row["in"] == row["passed_through"] + row["flagged"] + row["excluded"], row


CURATED = [
    pytest.param(e, id=p.stem)
    for p in EXPECTED
    if "curated_checks" in (e := yaml.safe_load(p.read_text()))
]


@pytest.mark.parametrize("expected", CURATED)
def test_curated_checks(expected):
    findings = by_id(audit(expected["file"]))
    for check in expected["curated_checks"]:
        if "absent" in check:
            assert not set(check["absent"]) & set(findings), check["source"]
        elif "glycans_attached_to" in check:
            assert set(check["glycans_attached_to"]) == attachment_residues(findings)
        else:
            check_finding(findings[check["finding_id"]], check)


def attachment_residues(findings: dict[str, dict]) -> set[str]:
    """Polymer residues that unrecognized groups are covalently attached to."""
    return {
        f"{atom['chain']}:{atom['seq_num']}{atom['ins_code']}"
        for finding in findings.values()
        if finding["rule_id"] == "unrecognized.chemistry_group"
        for item in finding["evidence"]
        if item["key"] == "attachment"
        for atom in item["atoms"]
        if atom["res_name"] == "ASN"
    }


def check_finding(finding: dict, check: dict) -> None:
    if "shell_residues" in check:
        residues = {residue for residue, _ in shell_atoms(finding)}
        assert set(check["shell_residues"]) == residues, check["source"]
    for key, value in check.get("evidence", {}).items():
        actual = next(i.get("value") for i in finding["evidence"] if i["key"] == key)
        assert actual == value, (finding["id"], key)
