"""Hydrogens with OpenMM Modeller from recorded states (I/O edge; TASK-009, ADR-0009).

Every titratable residue gets an explicit variant, except neutral His (variant None),
whose tautomer Modeller picks by hydrogen bonding; Modeller does that only when its pH
argument exceeds 6.5, and the argument affects no other residue here, so it is fixed at
NEUTRAL_HIS_PH. Non-standard polymer residues need hydrogen and bond definitions in
knowledge/protonation.yaml; single-atom ions get none. Heavy atoms keep their
coordinates; each hydrogen takes its parent's occupancy and B-factor.
"""

from __future__ import annotations

import random
import tempfile
from dataclasses import replace
from pathlib import Path
from xml.sax.saxutils import quoteattr

from simprep.protonate.pka import without_hydrogens
from simprep.protonate.states import ProtonationError, State
from simprep.relax.shell import TEMPLATED, WATER, segments
from simprep.structure.model import Atom, Residue, ResidueClass, Structure

NEUTRAL_HIS_PH = 7.0
COORDINATE_DECIMALS = 3
ANGSTROM_PER_NM = 10.0


def add_hydrogens(structure: Structure, states: list[State], protonation: dict) -> Structure:
    """``structure`` with every hydrogen removed and added again for ``states``."""
    import openmm
    from openmm import app, unit

    definitions = protonation["residue_definitions"]
    _load_definitions(app, definitions)
    heavy = without_hydrogens(structure)
    polymer, waters, _ = _partition(heavy, definitions)
    topology, positions, order = _topology(app, [*segments(polymer), waters])
    quantity = unit.Quantity([openmm.Vec3(*p) for p in positions], unit.angstrom)
    topology.createDisulfideBonds(quantity)
    modeller = app.Modeller(topology, quantity)
    variants = {s.residue: s.variant for s in states}
    method = protonation["method"]
    random.seed(method["seed"])  # Modeller places new hydrogens at random positions
    modeller.addHydrogens(
        pH=NEUTRAL_HIS_PH,
        variants=[variants.get(rid) for rid in order],
        platform=openmm.Platform.getPlatformByName(method["platform"]),
    )
    added = _hydrogens(modeller, order, unit)
    residues = tuple(_with(r, added.get(r.id, ())) for r in heavy.residues)
    return replace(heavy, residues=residues)


def _partition(structure: Structure, definitions: dict) -> tuple[list, list, list]:
    """(polymer residues to protonate, waters, residues left as they are); anything else
    stops: unknown chemistry is a hard stop (CLAUDE.md rule 3)."""
    polymer, waters, left, unknown = [], [], [], []
    for residue in structure.residues:
        if residue.name == WATER:
            waters.append(residue)
        elif residue.residue_class is ResidueClass.POLYMER:
            (
                polymer if residue.name in TEMPLATED or residue.name in definitions else unknown
            ).append(residue)
        elif len(residue.heavy_atoms) == 1:
            left.append(residue)  # a monatomic ion has no hydrogens
        else:
            unknown.append(residue)
    if unknown:
        names = ", ".join(sorted({f"{r.name} {r.id.label()}" for r in unknown}))
        raise ProtonationError(
            f"no hydrogen definitions for {names}; add them to knowledge/protonation.yaml "
            "(residue_definitions) or exclude the residues in the manifest"
        )
    return polymer, waters, left


def _topology(app, groups: list[list[Residue]]) -> tuple:
    topology, positions, order = app.Topology(), [], []
    for group in groups:
        chain = topology.addChain()
        for residue in group:
            new = topology.addResidue(
                residue.name, chain, id=f"{residue.id.seq_num}{residue.id.ins_code}"
            )
            for atom in residue.heavy_atoms:
                topology.addAtom(atom.name, app.Element.getBySymbol(atom.element.title()), new)
                positions.append(atom.position)
            order.append(residue.id)
    topology.createStandardBonds()
    return topology, positions, order


def _hydrogens(modeller, order: list, unit) -> dict:
    """Residue id -> new hydrogen atoms (name, parent heavy-atom name, position in A)."""
    parents = {}
    for first, second in modeller.topology.bonds():
        for hydrogen, other in ((first, second), (second, first)):
            if hydrogen.element.symbol == "H":
                parents[hydrogen.index] = other.name
    values = modeller.positions.value_in_unit(unit.nanometer)
    added = {}
    for residue, rid in zip(modeller.topology.residues(), order, strict=True):
        added[rid] = tuple(
            (
                a.name,
                parents[a.index],
                tuple(round(c * ANGSTROM_PER_NM, COORDINATE_DECIMALS) for c in values[a.index]),
            )
            for a in residue.atoms()
            if a.element.symbol == "H"
        )
    return added


def _with(residue: Residue, hydrogens: tuple) -> Residue:
    parent = {a.name: a for a in residue.atoms}
    new = tuple(
        Atom(name, "H", position, parent[of].occupancy, parent[of].b_iso)
        for name, of, position in hydrogens
    )
    return replace(residue, atoms=residue.atoms + new)


def _load_definitions(app, definitions: dict) -> None:
    """Bond and hydrogen definitions for non-standard polymer residues, as OpenMM XML."""
    with tempfile.TemporaryDirectory() as directory:
        bonds, hydrogens = Path(directory) / "bonds.xml", Path(directory) / "hydrogens.xml"
        bonds.write_text(_xml(definitions, _bond_lines))
        hydrogens.write_text(_xml(definitions, _hydrogen_lines))
        app.Topology.loadBondDefinitions(str(bonds))
        app.Modeller.loadHydrogenDefinitions(str(hydrogens))


def _xml(definitions: dict, lines) -> str:
    body = [
        line
        for name, definition in sorted(definitions.items())
        for line in (f"<Residue name={quoteattr(name)}>", *lines(definition), "</Residue>")
    ]
    return "\n".join(["<Residues>", *body, "</Residues>"]) + "\n"


def _bond_lines(definition: dict) -> list[str]:
    return [f"<Bond from={quoteattr(a)} to={quoteattr(b)}/>" for a, b in definition["bonds"]]


def _hydrogen_lines(definition: dict) -> list[str]:
    return [
        f"<H name={quoteattr(h['name'])} parent={quoteattr(h['parent'])}"
        + (f" terminal={quoteattr(h['terminal'])}" if "terminal" in h else "")
        + "/>"
        for h in definition["hydrogens"]
    ]
