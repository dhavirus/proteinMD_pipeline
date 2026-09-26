"""Loop modelling, pure parts (TASK-007): which gaps are built, how the loop is inserted
and marked, the torsion restraints, and the geometry and contact checks. Structures are
hand-built; the backbone uses textbook trans-peptide internal coordinates."""

import json
from dataclasses import replace

import pytest

from simprep.model.checks import criterion as clash_criterion
from simprep.model.checks import gap_findings, geometry_failures, loop_clashes
from simprep.model.loops import (
    ModelError,
    gaps_to_model,
    loop_shell,
    torsion_restraints,
    with_loops,
)
from simprep.model.run import LoopRejected, Modelled, write_model
from simprep.prep.apply import apply_plan
from simprep.prep.plan import build_plan
from simprep.prep.run import PrepInputs
from simprep.relax.openmm_run import RelaxOutcome
from simprep.severity import apply_context
from simprep.structure.geometry import place_atom
from simprep.structure.model import ResidueId, UnobservedResidue
from simprep.variants.relaxation import Relaxed
from tests.unit.builders import P, atom, config, residue, structure

CHECKS = {
    "peptide_bond_angstrom": 1.33,
    "peptide_bond_tolerance_angstrom": 0.05,
    "ca_ca_trans_angstrom": 3.8,
    "ca_ca_trans_tolerance_angstrom": 0.1,
}
MARKING = {"occupancy": 0.0, "b_iso_angstrom2": 0.0}
# (bond A, angle degree) of N, CA, C and the omega/phi/psi used to place them.
N_CA, CA_C, C_N = 1.458, 1.525, 1.329
ANGLES = {"N": 121.7, "CA": 111.2, "C": 116.2}
SEQUENCE = ("GLY", "ALA", "ALA", "SER", "GLY")


def backbone(count: int, chirality: dict | None = None):
    """``count`` residues of an extended trans chain: {index: -1} mirrors that CB."""
    points = [(0.0, 1.4, 0.0), (0.0, 0.0, 0.0), (1.5, 0.0, 0.0)]  # N, CA, C of residue 1
    for _ in range(count - 1):
        for bond, angle, torsion in (
            (C_N, ANGLES["C"], 180.0),
            (N_CA, ANGLES["N"], 180.0),
            (CA_C, ANGLES["CA"], -120.0),
        ):
            points.append(place_atom(tuple(points[-3:]), (bond, angle, torsion)))
    chirality = chirality or {}
    residues = []
    for index in range(count):
        n, ca, c = points[3 * index : 3 * index + 3]
        cb = place_atom((c, n, ca), (1.53, 110.5, -122.5 * chirality.get(index, 1)))
        atoms = [atom("N", "N", n), atom("CA", "C", ca), atom("C", "C", c), atom("CB", "C", cb)]
        residues.append(residue("A", index + 1, SEQUENCE[index % 5], P, atoms, index + 1))
    return residues


def gapped(names=("ALA", "SER")):
    """A:1, A:2, A:5 observed; A:3-4 unobserved with ``names``; sequence SEQUENCE."""
    chain = backbone(5)
    observed = [chain[0], chain[1], chain[4]]
    unobserved = tuple(
        UnobservedResidue(ResidueId("A", n), name, n, True)
        for n, name in zip((3, 4), names, strict=True)
    )
    return structure(
        observed, unobserved_residues=unobserved, polymer_sequences=(("A", SEQUENCE),)
    ), chain


def gap_manifest(option: str, evidence=()):
    options = [
        {"id": "model_loop", "label": "Model", "prep_action": "apply", "prep_stage": "modelling"},
        {"id": "truncate", "label": "Truncate", "prep_action": "record", "prep_stage": "topology"},
    ]
    finding = {
        "id": "missing_residues/A:3-4",
        "rule_id": "missing_residues.internal",
        "title": "Internal gap: A:3-4 (2 residues)",
        "effective_severity": "warn",
        "options": options,
        "locus": {"extent": [{"chain": "A", "seq_num": n, "ins_code": ""} for n in (3, 4)]},
        "anchor_residues": [{"chain": "A", "seq_num": n, "ins_code": ""} for n in (2, 5)],
        "evidence": list(evidence),
    }
    decision = {
        "finding_id": finding["id"],
        "option_id": option,
        "rationale": "t",
        "decided_by": "t",
        "timestamp": "2026-09-26T00:00:00Z",
    }
    return {"findings_snapshot": {"findings": [finding]}, "decisions": [decision]}


def test_model_loop_is_left_to_the_modelling_stage(ruleset):
    site, _ = gapped()
    plan = build_plan(site, gap_manifest("model_loop"), ruleset)
    (action,) = plan.actions
    assert (action.prep_action, action.prep_stage) == ("apply", "modelling")
    assert [item.stage for item in plan.work_order] == ["modelling"]
    (system,) = apply_plan(site, plan)
    assert system.structure == site  # prep never invents coordinates


def test_truncate_names_the_charged_termini(ruleset):
    evidence = [
        {"key": "position", "type": "label", "value": "internal"},
        {"key": "flank_before", "type": "label", "value": "ALAA:2"},
        {"key": "flank_after", "type": "label", "value": "GLYA:5"},
    ]
    site, _ = gapped()
    (item,) = build_plan(site, gap_manifest("truncate", evidence), ruleset).work_order
    assert item.description.endswith(
        "chain break: charged C-terminus (COO-) at ALA A:2, charged N-terminus (NH3+) at GLY A:5"
    )


def test_the_gap_is_inserted_marked_and_removed_from_the_annotation(ruleset):
    site, chain = gapped()
    (gap,) = gaps_to_model(site, build_plan(site, gap_manifest("model_loop"), ruleset))
    assert gap.flanks == (ResidueId("A", 2), ResidueId("A", 5))
    placed = {r.id: r.atoms for r in chain[2:4]}
    built = with_loops(site, (gap,), (placed, MARKING))
    assert [r.id.seq_num for r in built.residues] == [1, 2, 3, 4, 5]
    assert [r.label_seq for r in built.residues] == [1, 2, 3, 4, 5]
    assert all(a.occupancy == 0.0 for r in built.residues[2:4] for a in r.atoms)
    assert built.unobserved_residues == ()


def test_a_sequence_that_disagrees_with_the_entity_is_refused(ruleset):
    site, _ = gapped(names=("ALA", "TRP"))
    with pytest.raises(ModelError, match="TRP in the unobserved-residue annotation"):
        gaps_to_model(site, build_plan(site, gap_manifest("model_loop"), ruleset))


def test_torsion_restraints_cover_every_peptide_bond_and_non_glycine_residue(ruleset):
    site, chain = gapped()
    (gap,) = gaps_to_model(site, build_plan(site, gap_manifest("model_loop"), ruleset))
    built = with_loops(site, (gap,), ({r.id: r.atoms for r in chain[2:4]}, MARKING))
    restraints = torsion_restraints(built, (gap,), 34.0)
    omegas = [r for r in restraints if r.target_degree == 180.0]
    assert len(omegas) == 3  # A:2-3, A:3-4, A:4-5
    assert [r.atoms[0][0].seq_num for r in restraints if r.target_degree == 34.0] == [3, 4]


def modelled(ruleset, chain):
    site, _ = gapped()
    (gap,) = gaps_to_model(site, build_plan(site, gap_manifest("model_loop"), ruleset))
    return with_loops(site, (gap,), ({r.id: r.atoms for r in chain[2:4]}, MARKING)), gap


def test_ideal_trans_geometry_passes_and_distortions_are_named(ruleset):
    built, gap = modelled(ruleset, backbone(5))
    assert geometry_failures(built, gap, CHECKS) == []
    built, gap = modelled(ruleset, backbone(5, chirality={3: -1}))
    (failure,) = geometry_failures(built, gap, CHECKS)
    assert failure["key"] == "ca_improper" and failure["atoms"][0]["seq_num"] == 4
    stretched = backbone(5)
    moved = [
        replace(a, position=(a.position[0] + 0.3, *a.position[1:])) for a in stretched[3].atoms
    ]
    stretched[3] = replace(stretched[3], atoms=tuple(moved))
    built, gap = modelled(ruleset, stretched)
    assert {f["key"] for f in geometry_failures(built, gap, CHECKS)} >= {"c_n_distance"}


def test_contacts_skip_sequence_neighbours_but_not_the_rest(ruleset):
    chain = backbone(5)
    built, gap = modelled(ruleset, chain)
    assert loop_clashes(built, gap, clash_criterion(ruleset)) == []
    clash_point = next(a.position for a in built.residue(ResidueId("A", 4)).atoms if a.name == "CB")
    intruder = residue("A", 900, "HOH", P, [atom("O", "O", clash_point)])
    crowded = replace(built, residues=(*built.residues, intruder))
    found = loop_clashes(crowded, gap, clash_criterion(ruleset))
    assert {c.residue.id for _, c in found} == {ResidueId("A", 900)}
    assert ResidueId("A", 4) in {owner.id for owner, _ in found}


def test_a_rejected_loop_writes_its_record_and_no_system(ruleset, tmp_path):
    """A blocking geometry finding: record with status rejected, no wt_modelled/."""
    pytest.importorskip("pdbfixer")  # the record names the installed PDBFixer version
    built, gap = modelled(ruleset, backbone(5, chirality={3: -1}))
    raw = gap_findings(built, gap, (ruleset, CHECKS))
    findings = [f.to_dict() for f in apply_context(raw, built, config())]
    outcome = RelaxOutcome(built, (), 0, 0.0, 0.0, (0.0, 0.0), 0.0, (), "test")
    relaxed = Relaxed("wt_modelled", outcome, loop_shell(built, (gap,)), (), ())
    (tmp_path / "wt").mkdir()
    (tmp_path / "wt" / "prep_record.json").write_text("{}")
    source = {"path": "test.cif", "format": "mmcif", "sha256": "0" * 64}
    inputs = PrepInputs(source, {}, "0" * 64, ruleset, built)
    kb = {"version": ruleset.version, "sha256": ruleset.sha256}
    wt_record = {"schema_version": "0.1.0", "knowledge_base": kb, "work_order": []}
    rejected = Modelled(built, (gap,), relaxed, findings, ruleset.modelling)
    with pytest.raises(LoopRejected, match="Modelled loop fails the peptide-geometry checks"):
        write_model(rejected, (inputs, wt_record), tmp_path)
    record = json.loads((tmp_path / "model_record.json").read_text())
    assert (record["status"], record["system"]) == ("rejected", None)
    assert [f["rule_id"] for f in record["findings"]] == ["modelling.junction_geometry"]
    assert not (tmp_path / "wt_modelled").exists()
