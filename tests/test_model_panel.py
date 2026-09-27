"""Loop modelling on 5FQL (TASK-007): `simprep model` for the wild type, and `simprep
variants` building R468Q / R468W on the modelled wild type. Expected values come from
tests/panel/model_fixtures/5FQL.yaml (written by hand from the spec and the file's own
records). Needs PDBFixer and OpenMM (simprep[relax]); each loop build takes about two
minutes, so the module builds it three times only: `simprep model`, variants pass 1 and
variants pass 2."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest
import yaml

from simprep.knowledge import load_ruleset
from simprep.structure.geometry import dihedral_degree
from simprep.structure.model import ResidueClass, ResidueId
from simprep.structure.parse import read_structure
from tests.prep_helpers import decided_manifest
from tests.test_variants_panel import FIXTURE as VARIANTS
from tests.test_variants_panel import PANEL_DIR, cli, decide, evidence, variant_findings

pytest.importorskip("pdbfixer")

MODEL = yaml.safe_load((PANEL_DIR / "model_fixtures" / "5FQL.yaml").read_text())
STRUCTURE = PANEL_DIR / MODEL["structure"]
PREP = yaml.safe_load((PANEL_DIR / "prep_fixtures" / MODEL["prep_fixture"]).read_text())
EXPECTED = MODEL["expected"]
LOOP = EXPECTED["loop"]
LOOP_IDS = [ResidueId(LOOP["chain"], num) for num, _ in LOOP["residues"]]
FLANKS = [ResidueId(LOOP["chain"], num) for num in LOOP["flanks"]]


def manifest(path: Path, variants: list[dict] | None = None) -> Path:
    decided_manifest(STRUCTURE, PREP["decisions"] | MODEL["decisions"], path)
    if variants:
        document = json.loads(path.read_text())
        document["variants"] = variants
        document["relaxation"] = {**load_ruleset().relaxation, "enabled": False}
        path.write_text(json.dumps(document))
    return path


@pytest.fixture(scope="module")
def modelled(tmp_path_factory):
    root = tmp_path_factory.mktemp("model")
    assert cli("model", STRUCTURE, "--manifest", manifest(root / "m.json"), "--out", root) == 0
    return json.loads((root / "model_record.json").read_text()), root


@pytest.fixture(scope="module")
def variants(tmp_path_factory):
    """R468Q and R468W with their recommended rotamers, rigid, on the modelled wild type."""
    root = tmp_path_factory.mktemp("variants")
    start = manifest(root / "m.json", VARIANTS["variants"])
    assert cli("variants", STRUCTURE, "--manifest", start, "--out", root / "p1") == 3
    findings = variant_findings(root / "p1" / "manifest.json")
    choices = {
        name: ("choose_rotamer", evidence(finding, "recommended_rotamer"))
        for name, finding in findings.items()
    }
    decided = decide(root / "p1" / "manifest.json", choices, root / "decided.json")
    assert cli("variants", STRUCTURE, "--manifest", decided, "--out", root / "out") == 0
    return root / "out"


def residues(path: Path) -> dict:
    return {r.id: r for r in read_structure(path).residues}


def position(residue, name: str):
    return next(a.position for a in residue.atoms if a.name == name)


def test_the_loop_has_the_entity_sequence_between_its_flanks(modelled):
    _, root = modelled
    structure = read_structure(root / "wt_modelled" / "system.cif")
    chain = [r for r in structure.polymer_residues(LOOP["chain"])]
    ids = [r.id for r in chain]
    start = ids.index(FLANKS[0])
    assert ids[start + 1 : start + 1 + len(LOOP_IDS)] == LOOP_IDS
    assert ids[start + 1 + len(LOOP_IDS)] == FLANKS[1]
    names = [structure.residue(rid).name for rid in LOOP_IDS]
    assert names == [name for _, name in LOOP["residues"]]
    assert sum(len(structure.residue(rid).heavy_atoms) for rid in LOOP_IDS) == LOOP["heavy_atoms"]


def test_every_peptide_bond_passes_and_every_residue_is_l(modelled):
    _, root = modelled
    index = residues(root / "wt_modelled" / "system.cif")
    segment = [index[rid] for rid in [FLANKS[0], *LOOP_IDS, FLANKS[1]]]
    low_cn, high_cn = EXPECTED["checks"]["c_n_angstrom"]
    low_ca, high_ca = EXPECTED["checks"]["ca_ca_angstrom"]
    for first, second in zip(segment, segment[1:], strict=False):
        assert low_cn <= math.dist(position(first, "C"), position(second, "N")) <= high_cn
        assert low_ca <= math.dist(position(first, "CA"), position(second, "CA")) <= high_ca
    for residue in segment[1:-1]:
        improper = dihedral_degree(*(position(residue, n) for n in ("CA", "N", "C", "CB")))
        assert improper > 0, residue.id


def test_modelled_atoms_are_marked_and_listed(modelled):
    record, root = modelled
    index = residues(root / "wt_modelled" / "system.cif")
    occupancy = EXPECTED["marking"]["occupancy"]
    assert all(a.occupancy == occupancy for rid in LOOP_IDS for a in index[rid].atoms)
    (loop,) = record["loops"]
    listed = {(r["seq_num"], tuple(r["atoms"])) for r in loop["residues"]}
    assert listed == {(rid.seq_num, tuple(a.name for a in index[rid].atoms)) for rid in LOOP_IDS}
    assert "no experimental support" in loop["note"]
    assert loop["source"].startswith("PDBFixer 1.12.0 placement, then OpenMM")


def test_atoms_outside_the_loop_and_its_flanks_are_unmoved(modelled):
    _, root = modelled
    prepared = residues(root / "wt" / "system.cif")
    modelled_ = residues(root / "wt_modelled" / "system.cif")
    changed = {rid for rid, r in prepared.items() if modelled_[rid] != r}
    assert changed <= set(FLANKS)
    assert set(modelled_) - set(prepared) == set(LOOP_IDS)


def test_accounting_and_the_unobserved_annotation(modelled):
    record, root = modelled
    rows = {row["record_type"]: row for row in record["system"]["counts"]}
    for row in rows.values():
        assert row["in"] == row["passed_through"] + row["modified"] + row["excluded"], row
        assert row["out"] == row["passed_through"] + row["modified"] + row["added"], row
    accounting = EXPECTED["accounting"]
    assert rows["polymer residues"]["added"] == accounting["polymer_residues_added"]
    unobserved = rows["unobserved polymer residues (annotation)"]
    assert (unobserved["in"], unobserved["out"]) == (
        accounting["unobserved_in"],
        accounting["unobserved_out"],
    )
    written = read_structure(root / "wt_modelled" / "system.cif")
    assert len([u for u in written.unobserved_residues if u.is_polymer]) == unobserved["out"]


def test_the_truncated_n_terminus_is_recorded_without_atoms(modelled):
    record, root = modelled
    notes = [w["description"] for w in record["work_order"] if w["stage"] == "topology"]
    assert any(EXPECTED["truncated_n_terminus"] in note for note in notes)
    index = residues(root / "wt_modelled" / "system.cif")
    assert not {ResidueId("A", num) for num in EXPECTED["absent_residues"]} & set(index)
    modelling_left = [w for w in record["work_order"] if w["stage"] == "modelling"]
    assert all(w["finding_id"] != "missing_residues/A:444-453" for w in modelling_left)


def test_modelling_is_deterministic(modelled, variants):
    _, root = modelled
    for name in ("system.cif", "system.pdb"):
        first = (root / "wt_modelled" / name).read_bytes()
        assert first == (variants / "wt_modelled" / name).read_bytes()


@pytest.mark.parametrize("name", ["R468Q", "R468W"])
def test_every_variant_carries_the_same_loop(variants, name):
    wild_type = residues(variants / "wt_modelled" / "system.cif")
    variant = residues(variants / name / "system.cif")
    for rid in [*LOOP_IDS, *FLANKS]:
        assert variant[rid] == wild_type[rid]
    assert variant[ResidueId("A", 468)].residue_class is ResidueClass.POLYMER
    record = json.loads((variants / "variant_record.json").read_text())
    assert record["wild_type"]["modelled"]["directory"] == "wt_modelled"


CHARGE = {"ASP": -1, "GLU": -1, "HIP": 1, "LYS": 1}  # formal charge per recorded variant
PROTONATION = yaml.safe_load((PANEL_DIR / "protonation_fixtures" / "5FQL.yaml").read_text())


def test_the_modelled_wild_type_protonates_and_parameterizes_completely(modelled):
    """TASK-009/010 on the modelled wild type: every atom gets parameters; the net charge
    equals the formal charges of the recorded states, Arg, the termini and the ions."""
    import openmm

    from simprep.parameterize.run import run_parameterize
    from simprep.prep.run import PrepRequest, load_inputs
    from simprep.protonate.run import protonate, write_protonation

    record, root = modelled
    document = json.loads((root / "m.json").read_text())
    document["protonation"] = PROTONATION["protonation"]
    (root / "mp.json").write_text(json.dumps(document))
    inputs = load_inputs(PrepRequest(STRUCTURE, root / "mp.json", root))
    structure = read_structure(root / "wt_modelled" / "system.cif")
    wild_type = protonate("wt", structure, ("wt_modelled", "wt_protonated"), inputs)
    protonation = write_protonation([wild_type], (inputs, record["work_order"]), root)
    parameters = run_parameterize(root, root / "mp.json")
    (system,) = parameters["systems"]
    states = protonation["systems"][0]["states"]
    arginines = sum(1 for r in structure.polymer_residues("A") if r.name == "ARG")
    ions = sum(2 if r.name == "CA" else -1 for r in structure.residues if r.name in ("CA", "CL"))
    termini = 1 - 1  # N-terminal Thr34 NH3+, C-terminal Pro550 COO-
    expected = sum(CHARGE.get(s["variant"], 0) for s in states) + arginines + ions + termini
    assert system["net_charge"] == expected
    assert system["checks"]["angles"] == system["checks"]["angles_expected"]
    (ion,) = system["ion_models"]
    assert (ion["res_name"], ion["seq_num"], ion["model"]) == ("CA", 1551, "lj1264")
    built = openmm.XmlSerializer.deserialize(
        (root / system["directory"] / "system.xml").read_text()
    )
    c4 = next(f for f in built.getForces() if isinstance(f, openmm.CustomNonbondedForce))
    water_oxygens = [
        a.index
        for r in _topology(root).residues()
        if r.name == "HOH"
        for a in r.atoms()
        if a.element.symbol == "O"
    ]
    assert c4.getParticleParameters(water_oxygens[0])[1] == pytest.approx(
        87.3 * 4.184e-4
    )  # C4(Ca, OW)


def _topology(root):
    from openmm import app

    return app.PDBFile(str(root / "wt_protonated" / "system.pdb")).topology
