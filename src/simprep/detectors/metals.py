"""Metal-ion sites: coordination shell, geometry descriptor, classification."""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

from simprep.config import AuditConfig
from simprep.detectors import register
from simprep.detectors.common import is_standard_polymer, primary_atom, standard_polymer_residues
from simprep.findings import (
    Finding,
    atom_ref,
    ev_b_factor,
    ev_count,
    ev_distance,
    ev_label,
    ev_number,
    ev_occupancy,
    ev_source_record,
    partner_ref,
    residue_locus,
)
from simprep.rules import Rule, RuleSet
from simprep.structure.geometry import angle_degree, distance
from simprep.structure.model import Atom, Residue, ResidueClass, Structure

TETRAHEDRAL_DEGREE = math.degrees(math.acos(-1.0 / 3.0))
# Sorted L-M-L angles of ideal polyhedra (all ligand pairs).
IDEAL_ANGLES_DEGREE: dict[str, tuple[float, ...]] = {
    "linear": (180.0,),
    "trigonal_planar": (120.0,) * 3,
    "tetrahedral": (TETRAHEDRAL_DEGREE,) * 6,
    "square_planar": (90.0,) * 4 + (180.0,) * 2,
    "trigonal_bipyramidal": (90.0,) * 6 + (120.0,) * 3 + (180.0,),
    "square_pyramidal": (90.0,) * 8 + (180.0,) * 2,
    "octahedral": (90.0,) * 12 + (180.0,) * 3,
    "pentagonal_bipyramidal": (72.0,) * 5 + (90.0,) * 10 + (144.0,) * 5 + (180.0,),
}
OCCUPANCY_DECIMALS = 2  # occupancies are deposited with two decimals
LIGAND_CLASSES = ("protein", "water", "nonstandard_polymer_residue", "nonpolymer_ligand")


@dataclass(frozen=True)
class Ligand:
    """One donor atom in a metal's coordination shell."""

    residue: Residue
    atom: Atom
    distance_angstrom: float
    ligand_class: str
    shell_occupancy: float


@register("metals")
def detect_metals(structure: Structure, ruleset: RuleSet, config: AuditConfig) -> list[Finding]:
    """One finding per single-element metal-ion residue matched by a metals rule."""
    rules = ruleset.family("metals")
    standard = standard_polymer_residues(ruleset)
    findings = []
    for residue in structure.residues:
        rule = _matching_rule(residue, rules)
        if rule is not None:
            findings.append(_metal_finding(structure, residue, (rule, standard)))
    return findings


def _matching_rule(residue: Residue, rules: tuple[Rule, ...]) -> Rule | None:
    """First rule whose elements contain the residue's sole element (ions only)."""
    elements = {atom.element for atom in residue.atoms}
    if residue.residue_class is not ResidueClass.NONPOLYMER or len(elements) != 1:
        return None
    (element,) = elements
    return next((rule for rule in rules if element in rule.matcher["elements"]), None)


def coordination_shell(
    structure: Structure, metal: Residue, rule_standard: tuple[Rule, frozenset[str]]
) -> list[Ligand]:
    """Donor atoms within the rule's cutoff, sorted by distance.

    Each donor is represented by its closest altloc record; ``shell_occupancy`` sums the
    occupancies of all of that atom's records inside the cutoff.
    """
    rule, standard = rule_standard
    center = primary_atom(metal.atoms).position
    records: dict[tuple, list[tuple[float, Atom, Residue]]] = {}
    for residue in structure.residues:
        if residue.id == metal.id:
            continue
        for atom in residue.heavy_atoms:
            if atom.element not in rule.matcher["donor_elements"]:
                continue
            d = distance(center, atom.position)
            if d <= rule.matcher["coordination_cutoff_angstrom"]:
                records.setdefault((residue.id, atom.name), []).append((d, atom, residue))
    ligands = [_ligand(found, standard) for found in records.values()]
    return sorted(ligands, key=lambda lig: (lig.distance_angstrom, lig.residue.id))


def _ligand(records: list[tuple[float, Atom, Residue]], standard: frozenset[str]) -> Ligand:
    d, atom, residue = min(records, key=lambda record: (record[0], record[1].altloc))
    occupancy = sum(record[1].occupancy for record in records)
    return Ligand(residue, atom, d, _ligand_class(residue, standard), occupancy)


def split_by_occupancy(ligands: list[Ligand], rule: Rule) -> tuple[list[Ligand], list[Ligand]]:
    """(counted, excluded): donors below the rule's minimum shell occupancy do not count."""
    minimum = rule.matcher["min_donor_occupancy"]
    counted = [lig for lig in ligands if round(lig.shell_occupancy, OCCUPANCY_DECIMALS) >= minimum]
    return counted, [lig for lig in ligands if lig not in counted]


def _ligand_class(residue: Residue, standard: frozenset[str]) -> str:
    if residue.residue_class is ResidueClass.WATER:
        return "water"
    if residue.residue_class is ResidueClass.POLYMER:
        return (
            "protein" if is_standard_polymer(residue, standard) else "nonstandard_polymer_residue"
        )
    return "nonpolymer_ligand"


def geometry_descriptor(
    center: tuple, ligands: list[Ligand], rule: Rule
) -> tuple[str, float | None]:
    """Best-matching ideal polyhedron and its RMS angle deviation (degrees)."""
    candidates = rule.matcher.get("geometry_candidates", {}).get(str(len(ligands)), [])
    if not candidates:
        return "not_evaluated", None
    measured = sorted(
        angle_degree(a.atom.position, center, b.atom.position)
        for a, b in itertools.combinations(ligands, 2)
    )
    scored = [(_rms(measured, IDEAL_ANGLES_DEGREE[name]), name) for name in candidates]
    rms, name = min(scored)
    return name, rms


def _rms(measured: list[float], ideal: tuple[float, ...]) -> float:
    return math.sqrt(sum((m - i) ** 2 for m, i in zip(measured, ideal, strict=True)) / len(ideal))


def classify(ligands: list[Ligand], classification: dict) -> tuple[str, str]:
    """Apply the rule's stated criteria in order; return (class, basis)."""
    counts = {name: sum(lig.ligand_class == name for lig in ligands) for name in LIGAND_CLASSES}
    for ligand_class in classification["catalytic_if_shell_contains"]:
        if counts[ligand_class]:
            return "likely_catalytic", f"shell contains {ligand_class}"
    open_site = classification["open_site"]
    if open_site and (
        counts["protein"] <= open_site["max_protein_donors"]
        and counts["water"] >= open_site["min_water_donors"]
    ):
        return "likely_catalytic", (
            f"open site: {counts['protein']} protein donors, {counts['water']} water donors"
        )
    if counts["protein"] >= classification["structural_min_protein_donors"]:
        return (
            "likely_structural",
            f"{counts['protein']} protein donors, no catalytic criterion met",
        )
    return "ambiguous", "no criterion met"


def _metal_finding(
    structure: Structure, metal: Residue, rule_standard: tuple[Rule, frozenset[str]]
) -> Finding:
    rule = rule_standard[0]
    ligands, partial = split_by_occupancy(coordination_shell(structure, metal, rule_standard), rule)
    center_atom = primary_atom(metal.atoms)
    label, basis = classify(ligands, rule.matcher["classification"])
    return Finding(
        id=f"metals/{metal.id.label()}",
        rule=rule,
        title=f"{metal.name} site {metal.id.label()}: {label}, coordination number {len(ligands)}",
        locus=residue_locus(metal, (center_atom.name,)),
        anchor_residues=(metal.id,),
        evidence=tuple(
            _site_evidence(metal, ligands, rule)
            + [ev_label("classification", label, note=basis)]
            + _ligand_evidence(metal, center_atom, (ligands, "ligand_distance"))
            + _ligand_evidence(metal, center_atom, (partial, "excluded_partial_donor"))
            + _annotated_links(structure, metal)
            + _crystallization_record(structure)
        ),
        claims=frozenset({metal.id}),
    )


def _site_evidence(metal: Residue, ligands: list[Ligand], rule: Rule) -> list[dict]:
    center = primary_atom(metal.atoms)
    geometry, rms = geometry_descriptor(center.position, ligands, rule)
    items = [
        ev_occupancy("metal_occupancy", center.occupancy),
        ev_b_factor("metal_b_factor", center.b_iso),
        ev_number("coordination_cutoff", rule.matcher["coordination_cutoff_angstrom"], "angstrom"),
        ev_number("min_donor_occupancy", rule.matcher["min_donor_occupancy"], "fraction"),
        ev_count("coordination_number", len(ligands)),
    ]
    items += [
        ev_count(f"{name}_donors", sum(lig.ligand_class == name for lig in ligands))
        for name in LIGAND_CLASSES
    ]
    if ligands:
        mean = sum(lig.distance_angstrom for lig in ligands) / len(ligands)
        items.append(ev_number("mean_ligand_distance", mean, "angstrom"))
    items.append(ev_label("geometry", geometry))
    if rms is not None:
        items.append(ev_number("geometry_rms_angle_deviation", rms, "degree"))
    if metal.altlocs:
        items.append(
            ev_label(
                "metal_altlocs",
                ",".join(metal.altlocs),
                note="shell computed from the highest-occupancy record",
            )
        )
    return items


def _ligand_evidence(metal: Residue, center: Atom, ligands_key: tuple) -> list[dict]:
    ligands, key = ligands_key
    return [
        ev_distance(
            key,
            lig.distance_angstrom,
            [atom_ref(metal, center), atom_ref(lig.residue, lig.atom)],
            note=_ligand_note(lig),
        )
        for lig in ligands
    ]


def _annotated_links(structure: Structure, metal: Residue) -> list[dict]:
    """The file's own metalc records for this ion, verbatim, for comparison with the shell."""
    items = []
    for link in structure.links:
        if link.conn_type != "metalc" or metal.id not in (
            link.partner1.residue,
            link.partner2.residue,
        ):
            continue
        items.append(
            ev_source_record(
                "annotated_metal_link",
                "struct_conn",
                {
                    "id": link.conn_id,
                    "conn_type": link.conn_type,
                    "distance_angstrom": link.distance_angstrom,
                },
            )
            | {"atoms": [partner_ref(link.partner1), partner_ref(link.partner2)]}
        )
    return items


def _crystallization_record(structure: Structure) -> list[dict]:
    """The file's crystallization conditions, verbatim: was the ion a buffer component?"""
    if structure.crystallization_details is None:
        return []
    return [
        ev_source_record(
            "crystallization_conditions",
            "exptl_crystal_grow",
            {"pdbx_details": structure.crystallization_details},
        )
    ]


def _ligand_note(ligand: Ligand) -> str:
    """Ligand class, plus shell occupancy and altloc when the donor is not fully occupied."""
    atom = ligand.atom
    if ligand.shell_occupancy >= 1.0 and not atom.altloc:
        return ligand.ligand_class
    return f"{ligand.ligand_class}, occupancy in shell {ligand.shell_occupancy:.2f}" + (
        f", closest altloc {atom.altloc}" if atom.altloc else ""
    )
