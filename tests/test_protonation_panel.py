"""Protonation of 5FQL at pH 7.2 (TASK-009): `simprep protonate` on the prepared wild type
and `simprep variants` protonating the wild type and R468Q. Expected values come from
tests/panel/protonation_fixtures/5FQL.yaml (written by hand from the file's records, the
CCD and the spec's measurements). Needs PROPKA and OpenMM (simprep[relax])."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from simprep.structure.model import ResidueId
from simprep.structure.parse import read_structure
from tests.test_variants_panel import FIXTURE as VARIANTS
from tests.test_variants_panel import (
    PANEL_DIR,
    cli,
    decide,
    evidence,
    manifest_with,
    variant_findings,
)

pytest.importorskip("propka")

FIXTURE = yaml.safe_load((PANEL_DIR / "protonation_fixtures" / "5FQL.yaml").read_text())
STRUCTURE = PANEL_DIR / FIXTURE["structure"]
EXPECTED = FIXTURE["expected"]


def label(text: str) -> ResidueId:
    chain, num = text.split(":")
    return ResidueId(chain, int(num))


def with_protonation(path: Path, variants: list[dict] | None = None) -> Path:
    manifest_with(variants or [], path)
    document = json.loads(path.read_text())
    document["protonation"] = FIXTURE["protonation"]
    if not variants:
        del document["variants"]
    path.write_text(json.dumps(document))
    return path


@pytest.fixture(scope="module")
def protonated(tmp_path_factory):
    root = tmp_path_factory.mktemp("protonate")
    manifest = with_protonation(root / "m.json")
    assert cli("protonate", STRUCTURE, "--manifest", manifest, "--out", root / "out") == 0
    record = json.loads((root / "out" / "protonation_record.json").read_text())
    return record, root


def states(record: dict, system: str = "wt") -> dict:
    entry = next(s for s in record["systems"] if s["name"] == system)
    return {ResidueId(s["chain"], s["seq_num"], s["ins_code"]): s for s in entry["states"]}


def atom_names(structure, text: str) -> set[str]:
    return {a.name for a in structure.residue(label(text)).atoms}


def test_metal_ligands_and_disulfides_follow_the_rules(protonated):
    record, root = protonated
    by_residue = states(record)
    structure = read_structure(root / "out" / "wt_protonated" / "system.cif")
    for text, variant in EXPECTED["metal_rule"].items():
        assert by_residue[label(text)]["variant"] == variant
        assert by_residue[label(text)]["basis"].startswith("metal rule")
    assert "HD1" in atom_names(structure, "A:335") and "HE2" not in atom_names(structure, "A:335")
    assert not {"HD1", "HD2"} & atom_names(structure, "A:45")
    for text in EXPECTED["disulfide"]:
        assert by_residue[label(text)]["variant"] == "CYX"
        assert "HG" not in atom_names(structure, text)


def test_estimates_near_ignored_chemistry_are_flagged_with_the_standard_default(protonated):
    record, _ = protonated
    (system,) = record["systems"]
    unreliable = [
        f for f in system["findings"] if f["rule_id"] == "protonation.unreliable_estimate"
    ]
    assert sorted(f["id"].split("/")[-1] for f in unreliable) == EXPECTED["unreliable"]
    assert all(f["recommended_option"] == EXPECTED["unreliable_default"] for f in unreliable)
    by_residue = states(record)
    for text in EXPECTED["unreliable"]:
        assert by_residue[label(text)]["basis"].startswith(EXPECTED["unreliable_default"])


def test_ddz_gets_its_hydrogens_and_heavy_atoms_are_unchanged(protonated):
    _, root = protonated
    structure = read_structure(root / "out" / "wt_protonated" / "system.cif")
    hydrogens = {n for n in atom_names(structure, "A:84") if n.startswith("H")}
    assert hydrogens == set(EXPECTED["ddz_hydrogens"])
    prepared = read_structure(root / "out" / "wt" / "system.cif")
    heavy = {(r.id, a.name, a.position) for r in prepared.residues for a in r.heavy_atoms}
    assert heavy == {(r.id, a.name, a.position) for r in structure.residues for a in r.heavy_atoms}
    polymer = structure.polymer_residues("A")
    assert all(any(a.element == "H" for a in r.atoms) for r in polymer)


def test_record_has_the_ph_and_no_protonation_work_left(protonated):
    record, _ = protonated
    assert record["ph"] == FIXTURE["protonation"]["ph"]
    assert all(w["stage"] != "protonation" for w in record["work_order"])


def test_protonation_is_deterministic(protonated, tmp_path):
    _, root = protonated
    assert cli("protonate", STRUCTURE, "--manifest", root / "m.json", "--out", tmp_path) == 0
    for name in ("system.cif", "system.pdb"):
        first = (root / "out" / "wt_protonated" / name).read_bytes()
        assert first == (tmp_path / "wt_protonated" / name).read_bytes()


def test_a_missing_ph_stops_before_any_work(tmp_path, capsys):
    manifest = manifest_with([], tmp_path / "m.json")
    assert cli("protonate", STRUCTURE, "--manifest", manifest, "--out", tmp_path / "out") == 2
    assert "states no pH" in capsys.readouterr().err
    assert not (tmp_path / "out").exists()


def test_variants_protonate_every_system_and_compare_with_the_wild_type(tmp_path):
    variant = next(v for v in VARIANTS["variants"] if v["name"] == "R468Q")
    start = with_protonation(tmp_path / "m.json", [variant])
    assert cli("variants", STRUCTURE, "--manifest", start, "--out", tmp_path / "p1") == 3
    finding = variant_findings(tmp_path / "p1" / "manifest.json")["R468Q"]
    rotamer = evidence(finding, "recommended_rotamer")
    decided = decide(
        tmp_path / "p1" / "manifest.json",
        {"R468Q": ("choose_rotamer", rotamer)},
        tmp_path / "d.json",
    )
    assert cli("variants", STRUCTURE, "--manifest", decided, "--out", tmp_path / "out") == 0
    record = json.loads((tmp_path / "out" / "protonation_record.json").read_text())
    assert [s["name"] for s in record["systems"]] == ["wt", "R468Q"]
    mutant = read_structure(tmp_path / "out" / "R468Q_protonated" / "system.cif")
    gln = mutant.residue(ResidueId("A", 468))
    assert gln.name == "GLN" and {"HE21", "HE22"} <= {a.name for a in gln.atoms}
