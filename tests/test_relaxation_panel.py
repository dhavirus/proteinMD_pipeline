"""Relaxation of 5FQL R468Q / R468W and the matched wild type (TASK-006). Expected
outcomes are the hand-written limits in tests/panel/variant_fixtures/5FQL.yaml (from the
TASK-006 spike), not values read back from this implementation. Needs OpenMM
(simprep[relax])."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from simprep.structure.model import ResidueId
from simprep.structure.parse import read_structure
from tests.test_variants_panel import (
    FIXTURE,
    STRUCTURE,
    cli,
    decide,
    manifest_with,
    variant_findings,
)

pytest.importorskip("openmm")

EXPECTED = FIXTURE["relaxation"]
SITE = ResidueId("A", 468)


@pytest.fixture(scope="module")
def relaxed(tmp_path_factory):
    """Both variants with their recommended rotamers, relaxed (default protocol)."""
    root = tmp_path_factory.mktemp("relaxed")
    manifest = manifest_with(FIXTURE["variants"], root / "m.json", relaxation=True)
    assert cli("variants", STRUCTURE, "--manifest", manifest, "--out", root / "p1") == 3
    findings = variant_findings(root / "p1" / "manifest.json")
    choices = {
        name: (
            "choose_rotamer",
            next(e["value"] for e in f["evidence"] if e["key"] == "recommended_rotamer"),
        )
        for name, f in findings.items()
    }
    decided = decide(root / "p1" / "manifest.json", choices, root / "decided.json")
    assert cli("variants", STRUCTURE, "--manifest", decided, "--out", root / "out") == 0
    return json.loads((root / "out" / "variant_record.json").read_text()), root / "out"


def entry(record: dict, name: str) -> dict:
    return next(v for v in record["variants"] if v["name"] == name)


def test_outputs_and_protocol_are_recorded(relaxed):
    record, out = relaxed
    assert record["relaxation_protocol"]["platform"] == "Reference"
    (wild_type,) = record["relaxed_wild_types"]
    assert wild_type["directory"] == "wt_relaxed_A468"
    for path in (
        "wt/system.cif",
        "wt_relaxed_A468/system.cif",
        "R468Q/system.cif",
        "R468W/system.cif",
        "R468W/unrelaxed/system.cif",
    ):
        assert (out / path).is_file(), path


def test_r468w_clashes_drop_and_what_remains_is_a_warning(relaxed):
    record, _ = relaxed
    relaxation = entry(record, "R468W")["relaxation"]
    assert relaxation["clashes_before"] == EXPECTED["R468W"]["clashes_before"]
    assert relaxation["clashes_after"] <= EXPECTED["R468W"]["clashes_after_at_most"]
    findings = entry(record, "R468W")["findings"]
    assert bool(findings) is bool(relaxation["clashes_after"])
    assert all(f["effective_severity"] == "warn" for f in findings)


def test_r468q_stays_clash_free(relaxed):
    record, _ = relaxed
    assert (
        entry(record, "R468Q")["relaxation"]["clashes_after"]
        <= EXPECTED["R468Q"]["clashes_after_at_most"]
    )


def heavy_positions(path: Path) -> dict:
    return {
        (r.id, a.name): a.position
        for r in read_structure(path).residues
        for a in r.atoms
        if a.is_heavy
    }


def test_relaxed_wild_type_moves_little(relaxed):
    _, out = relaxed
    before, after = (
        heavy_positions(out / "wt/system.cif"),
        heavy_positions(out / "wt_relaxed_A468/system.cif"),
    )
    largest = max(math.dist(before[k], after[k]) for k in after)
    assert 0 < largest <= EXPECTED["wild_type_max_displacement_angstrom"]


def test_relaxed_wild_type_and_variant_differ_only_inside_the_shell(relaxed):
    record, out = relaxed
    wild_type, mutant = (
        read_structure(out / p) for p in ("wt_relaxed_A468/system.cif", "R468W/system.cif")
    )
    differing = {a.id for a, b in zip(wild_type.residues, mutant.residues, strict=True) if a != b}
    moved = {
        ResidueId(r["chain"], r["seq_num"], r["ins_code"])
        for system in (record["relaxed_wild_types"][0], entry(record, "R468W"))
        for r in system["relaxation"]["moved_residues"]
    }
    assert SITE in differing
    assert differing <= moved | {SITE}


def test_moved_residues_lose_hydrogens_and_get_protonation_items(relaxed):
    record, out = relaxed
    mutant = read_structure(out / "R468W/system.cif")
    moved = entry(record, "R468W")["relaxation"]["moved_residues"]
    for ref in moved:
        residue = mutant.residue(ResidueId(ref["chain"], ref["seq_num"], ref["ins_code"]))
        assert all(atom.is_heavy for atom in residue.atoms), residue.id
    protonated = {
        tuple(i["residues"][0].values())
        for i in entry(record, "R468W")["work_order"]
        if i["stage"] == "protonation"
    }
    assert {tuple(r.values()) for r in moved} <= protonated


def test_frozen_and_left_out_chemistry_is_recorded(relaxed):
    record, _ = relaxed
    relaxation = entry(record, "R468W")["relaxation"]
    # A:84 is the gem-diol DDZ after prep (TASK-008); neither it nor Ca2+ has a template.
    assert {r["res_name"] for r in relaxation["left_out"]} >= {"DDZ", "CA"}
    assert all(f["near"].startswith("DDZ") for f in relaxation["frozen"])
    assert relaxation["frozen"] == record["relaxed_wild_types"][0]["relaxation"]["frozen"]
