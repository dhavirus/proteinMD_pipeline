"""Relaxation building blocks (TASK-006): shell selection and freezing (pure), and the
OpenMM path on the prepared 3KS3 panel structure (small, with a template-less Zn2+)."""

from __future__ import annotations

import json
import math
import tempfile
from functools import cache
from pathlib import Path

import pytest
import yaml

from simprep.knowledge import load_ruleset
from simprep.prep.apply import apply_plan
from simprep.prep.plan import build_plan
from simprep.relax.shell import RelaxError, build_shell, heavy, segments, shell_residues
from simprep.structure.model import ResidueId
from simprep.structure.parse import read_structure
from tests.prep_helpers import decided_manifest
from tests.unit import builders as b

PANEL = Path(__file__).parents[1] / "panel"
ZINC = ResidueId("A", 262)
ZINC_LIGAND = ResidueId("A", 94)  # His94 NE2 coordinates the zinc (2.0 A in the file)


def peptide(*gaps):
    """Three residues along x; a gap after residue number n moves the next one away."""
    residues, x = [], 0.0
    for number in (1, 2, 3):
        atoms = [
            b.atom("N", "N", (x, 0, 0)),
            b.atom("CA", "C", (x + 1.46, 0, 0)),
            b.atom("C", "C", (x + 2.0, 1.4, 0)),
            b.atom("CB", "C", (x + 1.46, -1.5, 0)),
        ]
        residues.append(b.residue("A", number, "ALA", b.P, atoms))
        x += 3.33 + (5.0 if number in gaps else 0.0)
    return b.structure(residues)


def test_segments_split_where_the_peptide_bond_is_missing():
    assert [len(s) for s in segments(list(peptide().residues))] == [3]
    assert [len(s) for s in segments(list(peptide(1).residues))] == [1, 2]


def test_shell_moves_side_chains_and_frees_only_the_site():
    structure = peptide()
    sites = frozenset({ResidueId("A", 2)})
    shell = build_shell(structure, sites, shell_residues(structure, sites, 4.0), 1.0)
    assert shell.residues == {ResidueId("A", n) for n in (1, 2, 3)}
    assert shell.mobile == {(ResidueId("A", n), "CB") for n in (1, 2, 3)}
    assert shell.unrestrained == {(ResidueId("A", 2), "CB")}


@cache
def prepared_3ks3():
    """3KS3 after its prep fixture (altlocs collapsed, glycerol excluded, Zn kept)."""
    fixture = yaml.safe_load((PANEL / "prep_fixtures" / "3KS3.yaml").read_text())
    path = PANEL / "3KS3.cif.gz"
    with tempfile.TemporaryDirectory() as tmp:
        manifest = json.loads(
            decided_manifest(path, fixture["decisions"], Path(tmp) / "m.json").read_text()
        )
    structure = read_structure(path)
    (system,) = apply_plan(structure, build_plan(structure, manifest, load_ruleset()))
    return system.structure


def far_site(structure):
    """A lysine far from the zinc: its shell stays clear of template-less chemistry."""
    zinc = structure.residue(ZINC).atoms[0].position
    lysines = [r for r in structure.residues if r.name == "LYS"]
    return max(lysines, key=lambda r: min(math.dist(a.position, zinc) for a in heavy(r))).id


def test_shell_residues_near_the_zinc_are_frozen():
    structure = prepared_3ks3()
    sites = frozenset({ResidueId("A", 199)})  # Thr199, next to the zinc site
    shell = build_shell(structure, sites, shell_residues(structure, sites, 6.0), 1.0)
    frozen = {rid for rid, near, _ in shell.frozen}
    assert frozen and all(near.startswith("ZN") for _, near, _ in shell.frozen)
    assert not any(rid in frozen for rid, _ in shell.mobile)


openmm = pytest.importorskip("openmm")
from simprep.relax.openmm_run import relax  # noqa: E402


def test_a_site_next_to_template_less_chemistry_stops_relaxation():
    structure = prepared_3ks3()
    sites = frozenset({ZINC_LIGAND})
    shell = build_shell(structure, sites, shell_residues(structure, sites, 6.0), 1.0)
    with pytest.raises(RelaxError, match=r"ZN A:262 .* Parameterize them first"):
        relax(structure, shell, load_ruleset().relaxation)


def test_altlocs_are_refused():
    structure = read_structure(PANEL / "3KS3.cif.gz")
    sites = frozenset({far_site(structure)})
    shell = build_shell(structure, sites, shell_residues(structure, sites, 6.0), 1.0)
    with pytest.raises(RelaxError, match="one conformer per residue"):
        relax(structure, shell, load_ruleset().relaxation)


@cache
def relaxed_far_site():
    structure = prepared_3ks3()
    sites = frozenset({far_site(structure)})
    shell = build_shell(structure, sites, shell_residues(structure, sites, 6.0), 1.0)
    return structure, shell, relax(structure, shell, load_ruleset().relaxation)


def test_relaxation_is_deterministic():
    structure, shell, first = relaxed_far_site()
    second = relax(structure, shell, load_ruleset().relaxation)
    assert second.structure == first.structure


def test_only_moved_residues_change_and_they_lose_hydrogens():
    structure, shell, outcome = relaxed_far_site()
    changed = {
        a.id
        for a, b_ in zip(structure.residues, outcome.structure.residues, strict=True)
        if a != b_
    }
    assert changed == set(outcome.moved_residues) and changed <= shell.residues
    for rid in changed:
        assert all(atom.is_heavy for atom in outcome.structure.residue(rid).atoms)

    def largest_shift(residue_ids):
        return max(
            math.dist(a.position, b_.position)
            for rid in residue_ids
            for a, b_ in zip(
                heavy(structure.residue(rid)), outcome.structure.residue(rid).atoms, strict=True
            )
        )

    # The site side chain is unrestrained; everything else is held by strong restraints.
    assert largest_shift(changed - shell.sites) < largest_shift(shell.sites)
    assert {r["res_name"] for r in outcome.left_out} == {"ZN"}


def test_disulfides_get_cyx_templates():
    from simprep.relax.openmm_run import _disulfide_templates, _model

    structure = read_structure(PANEL / "1HZH.cif.gz")
    heavy_only = b.structure([r for r in structure.residues[:230] if r.residue_class is b.P])
    model = _model(heavy_only, openmm_modules(), load_ruleset().relaxation)
    templates = _disulfide_templates(model.topology)
    assert templates and set(templates.values()) == {"CYX"}
    assert all(residue.name == "CYS" for residue in templates)


def openmm_modules():
    from openmm import app, unit

    return openmm, app, unit
