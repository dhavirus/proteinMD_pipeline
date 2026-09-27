"""The parameterization stage at the I/O edge (TASK-010, ADR-0010).

``simprep parameterize DIR --manifest M`` reads DIR/protonation_record.json (written by
``simprep protonate`` or ``simprep variants``), builds an OpenMM System for each
protonated system with ff14SB + TIP3P + the generated DDZ file, gives the metal ions the
model the manifest decided (12-6-4), checks that every bond and angle has parameters, and
writes ``<name>_parameterized/system.xml`` and ``parameterization_record.json``.
"""

from __future__ import annotations

import json
import os
import xml.etree.ElementTree as ET
from pathlib import Path

from simprep.knowledge import load_ruleset
from simprep.manifest import check_input, load_manifest, write_json
from simprep.manifest.manifest import sha256_text
from simprep.parameterize.lj1264 import apply
from simprep.paths import KNOWLEDGE_DIR
from simprep.prep.plan import finding_residues
from simprep.protonate.hydrogens import load_definitions
from simprep.provenance import sha256_file, simprep_provenance, utc_now
from simprep.schemas import validate

PROTONATION_RECORD = "protonation_record.json"
RECORD_FILE = "parameterization_record.json"
SYSTEM_FILE = "system.xml"
SUFFIX = "_parameterized"
STAGE = "parameterization"
SOURCES = {
    "HOH": "TIP3P (OpenMM amber14/tip3p.xml)",
    "CL": "Cl- Joung-Cheatham 12-6 (OpenMM amber14/tip3p.xml)",
    "CA": "Ca2+ Li-Merz 12-6 (OpenMM amber14/tip3p.xml)",
    "DDZ": "DDZ: AM1-BCC + GAFF2 side chain, ff14SB backbone (knowledge/forcefield/ddz.xml)",
}
PROTEIN = "Amber ff14SB (OpenMM amber14/protein.ff14SB.xml)"
MODEL_SOURCE = {"lj1264": "12-6-4 Li-Merz TIP3P set (knowledge/forcefield/lj1264.yaml)"}


class ParameterizationError(ValueError):
    """A system cannot be parameterized completely; the message says what is missing."""


def run_parameterize(directory: Path, manifest_path: Path, knowledge: Path = KNOWLEDGE_DIR) -> dict:
    protonation = json.loads((directory / PROTONATION_RECORD).read_text())
    manifest = load_manifest(manifest_path)
    check_input(manifest, protonation["input"]["sha256"])
    ruleset = load_ruleset(knowledge)
    bundle = ruleset.parameterization
    models = _ion_models(manifest, bundle["protocol"]["metal_models"])
    systems = [
        _system(entry, (directory, bundle, models), ruleset.protonation["residue_definitions"])
        for entry in protonation["systems"]
    ]
    import openmm

    record = {
        "schema_version": "0.1.0",
        "generated_at": utc_now(),
        "simprep": simprep_provenance(),
        "knowledge_base": {"version": ruleset.version, "sha256": ruleset.sha256},
        "input": protonation["input"],
        "manifest_sha256": sha256_text(manifest_path.read_text()),
        "protonation_record_sha256": sha256_file(directory / PROTONATION_RECORD),
        "protocol": bundle["protocol"],
        "openmm_version": openmm.__version__,
        "systems": systems,
        "work_order": [w for w in protonation["work_order"] if w["stage"] != STAGE],
    }
    validate(record, "parameterization_record")
    write_json(record, directory / RECORD_FILE)
    return record


def _ion_models(manifest: dict, metal_models: dict) -> dict:
    """Residue label -> (finding id, model) for metals decided with a modelled option."""
    findings = {f["id"]: f for f in manifest["findings_snapshot"]["findings"]}
    models = {}
    for decision in manifest["decisions"]:
        model = metal_models.get(decision["option_id"])
        finding = findings.get(decision["finding_id"])
        if model and finding:
            for rid in finding_residues(finding):
                models[rid.label()] = (finding["id"], model)
    return models


def _system(entry: dict, context: tuple, definitions: dict) -> dict:
    directory, bundle, models = context
    import openmm
    from openmm import app

    load_definitions(app, definitions)
    pdb = app.PDBFile(str(directory / entry["directory"] / "system.pdb"))
    files = [*bundle["protocol"]["force_field_files"], *bundle["generated_files"]]
    forcefield = app.ForceField(*files)
    try:
        system = forcefield.createSystem(
            pdb.topology, nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False
        )
    except ValueError as error:
        raise ParameterizationError(
            f"{entry['name']}: {error}. A chain end left by `truncate` needs its terminal "
            "atoms (e.g. OXT), which the topology stage adds (TASK-011); other residues need "
            "parameters in knowledge/forcefield/"
        ) from error
    checks = _checks(pdb.topology, system, entry["name"])
    ions = _apply_ions(system, pdb.topology, (files, models, bundle["lj1264"]))
    name = f"{entry['name']}{SUFFIX}"
    (directory / name).mkdir(parents=True, exist_ok=True)
    path = directory / name / SYSTEM_FILE
    path.write_text(openmm.XmlSerializer.serialize(system))
    return {
        "name": entry["name"],
        "directory": name,
        "source": entry["directory"],
        "files": [
            {"path": SYSTEM_FILE, "format": "openmm_system_xml", "sha256": sha256_file(path)}
        ],
        "particles": system.getNumParticles(),
        "net_charge": round(_net_charge(system), 4),
        "residue_sources": _sources(pdb.topology, models),
        "ion_models": ions,
        "checks": checks,
    }


def _checks(topology, system, name: str) -> dict:
    """Every bond and every angle (from the bond graph) must have parameters: OpenMM skips
    an angle without parameters silently."""

    forces = {type(f).__name__: f for f in system.getForces()}
    neighbours = {}
    for a, b in topology.bonds():
        neighbours.setdefault(a.index, set()).add(b.index)
        neighbours.setdefault(b.index, set()).add(a.index)
    checks = {
        "bonds": forces["HarmonicBondForce"].getNumBonds(),
        "bonds_expected": topology.getNumBonds(),
        "angles": forces["HarmonicAngleForce"].getNumAngles(),
        "angles_expected": sum(len(n) * (len(n) - 1) // 2 for n in neighbours.values()),
    }
    if checks["bonds"] != checks["bonds_expected"] or checks["angles"] != checks["angles_expected"]:
        raise ParameterizationError(f"{name}: bonded terms without parameters: {checks}")
    return checks


def _classes(system, files: list[str], polarizabilities: dict) -> list[str]:
    """A class for every particle, found by its (sigma, epsilon): ParmEd's add12_6_4
    likewise assumes atom types with the same LJ terms share a polarizability, and it is
    checked here. Returns one representative class per particle."""
    import openmm
    from openmm import unit

    by_lj = {}
    for name in files:
        for tree in _trees(Path(name) if os.path.isabs(name) else _openmm_data() / name):
            type_class = {t.get("name"): t.get("class") for t in tree.iter("Type")}
            for atom in tree.iter("Atom"):
                if atom.get("sigma") is not None and atom.get("type") in type_class:
                    key = _lj_key(float(atom.get("sigma")), float(atom.get("epsilon")))
                    by_lj.setdefault(key, set()).add(type_class[atom.get("type")])
    nonbonded = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
    classes = []
    for index in range(system.getNumParticles()):
        _, sigma, epsilon = nonbonded.getParticleParameters(index)
        key = _lj_key(
            sigma.value_in_unit(unit.nanometer), epsilon.value_in_unit(unit.kilojoule_per_mole)
        )
        candidates = sorted(by_lj.get(key, ()))
        values = {
            polarizabilities[c]["value_cubic_angstrom"] for c in candidates if c in polarizabilities
        }
        if not candidates or len(values) != 1:
            raise ParameterizationError(
                f"particle {index}: LJ classes {candidates} give polarizabilities {sorted(values)}"
            )
        classes.append(next(c for c in candidates if c in polarizabilities))
    return classes


def _trees(path: Path) -> list:
    """The force-field file and every file it includes (<Include file=...>)."""
    tree = ET.parse(path)
    included = [p for i in tree.iter("Include") for p in _trees(path.parent / i.get("file"))]
    return [tree, *included]


def _lj_key(sigma: float, epsilon: float) -> tuple[float, float]:
    return round(sigma, 6), round(epsilon, 6)


def _openmm_data() -> Path:
    from openmm import app

    return Path(app.__file__).parent / "data"


def _apply_ions(system, topology, context: tuple) -> list[dict]:
    """Apply each decided ion model; one ion species per model in v0.1."""
    files, models, data = context
    applied, targets = [], {}
    for residue in topology.residues():
        label = f"{residue.chain.id}:{residue.id}"
        if label in models:
            targets.setdefault(models[label][1], []).append((residue, models[label][0]))
    for model, residues in targets.items():
        species = {r.name for r, _ in residues}
        if len(species) != 1 or next(iter(species)) not in data["ions"]:
            raise ParameterizationError(
                f"{model}: ions {sorted(species)} not in lj1264.yaml (one species per model)"
            )
        classes = _classes(system, files, data["polarizabilities"])
        indices = tuple(a.index for r, _ in residues for a in r.atoms())
        pairs = apply(system, (next(iter(species)), indices, classes), data)
        applied += [
            {
                "chain": r.chain.id,
                "seq_num": int(r.id),
                "ins_code": (r.insertionCode or "").strip(),
                "res_name": r.name,
                "finding_id": fid,
                "model": model,
                "c4_pairs": pairs,
            }
            for r, fid in residues
        ]
    return applied


def _net_charge(system) -> float:
    import openmm
    from openmm import unit

    nonbonded = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
    return sum(
        nonbonded.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
        for i in range(system.getNumParticles())
    )


def _sources(topology, models: dict) -> dict:
    counts = {}
    for residue in topology.residues():
        label = f"{residue.chain.id}:{residue.id}"
        if label in models:
            source = MODEL_SOURCE[models[label][1]]
        else:
            source = SOURCES.get(residue.name, PROTEIN if len(list(residue.atoms())) > 1 else None)
        if source is None:
            raise ParameterizationError(f"{residue.name} {label}: no parameter source recorded")
        counts[source] = counts.get(source, 0) + 1
    return dict(sorted(counts.items()))
