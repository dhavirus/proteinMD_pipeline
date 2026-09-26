"""Checks on a modelled loop (pure; TASK-007): peptide geometry from flank to flank, CA
chirality, backbone dihedrals, and heavy-atom contacts with the rest of the system."""

from __future__ import annotations

from dataclasses import dataclass

from simprep.findings import Finding, Locus, atom_ref, ev_distance, ev_label, ev_number
from simprep.model.loops import GLYCINE, Gap, atom, segment
from simprep.rules import RuleSet
from simprep.structure.geometry import dihedral_degree, distance
from simprep.structure.model import Atom, Residue, Structure
from simprep.variants.check import ClashCriterion, clashes, contacts

FAMILY = "modelling"
GEOMETRY_RULE = "modelling.junction_geometry"
CONTACTS_RULE = "modelling.loop_contacts"
NO_EXPERIMENTAL_SUPPORT = (
    "Modelled: no experimental support (occupancy 0.00). Analyses that rely on the "
    "experimental structure (contact maps, native-structure models) should treat these "
    "residues as modelled."
)


@dataclass(frozen=True)
class PeptideBond:
    first: Residue
    second: Residue
    c_n_angstrom: float
    ca_ca_angstrom: float
    omega_degree: float


@dataclass(frozen=True)
class Backbone:
    """phi, psi and the CA-N-C-CB improper (None for glycine) of a modelled residue."""

    residue: Residue
    phi_degree: float
    psi_degree: float
    improper_degree: float | None


def peptide_bonds(structure: Structure, gap: Gap) -> list[PeptideBond]:
    residues = segment(structure, gap)
    return [_bond(first, second) for first, second in zip(residues, residues[1:], strict=False)]


def _bond(first: Residue, second: Residue) -> PeptideBond:
    ca1, c1, n2, ca2 = (
        atom(r, n) for r, n in ((first, "CA"), (first, "C"), (second, "N"), (second, "CA"))
    )
    return PeptideBond(
        first,
        second,
        distance(c1.position, n2.position),
        distance(ca1.position, ca2.position),
        dihedral_degree(ca1.position, c1.position, n2.position, ca2.position),
    )


def backbone(structure: Structure, gap: Gap) -> list[Backbone]:
    residues = segment(structure, gap)
    return [_backbone(*residues[i - 1 : i + 2]) for i in range(1, len(residues) - 1)]


def _backbone(previous: Residue, residue: Residue, following: Residue) -> Backbone:
    n, ca, c = (atom(residue, name).position for name in ("N", "CA", "C"))
    improper = None
    if residue.name != GLYCINE:
        improper = dihedral_degree(ca, n, c, atom(residue, "CB").position)
    return Backbone(
        residue,
        dihedral_degree(atom(previous, "C").position, n, ca, c),
        dihedral_degree(n, ca, c, atom(following, "N").position),
        improper,
    )


def geometry_failures(structure: Structure, gap: Gap, checks: dict) -> list[dict]:
    """Evidence for every peptide bond or residue outside the protocol's checks."""
    failures = []
    for bond in peptide_bonds(structure, gap):
        failures += _bond_failures(bond, checks)
    for entry in backbone(structure, gap):
        if entry.improper_degree is not None and entry.improper_degree <= 0:
            residue = entry.residue
            atoms = [atom_ref(residue, atom(residue, n)) for n in ("CA", "N", "C", "CB")]
            improper = ev_number(
                "ca_improper", entry.improper_degree, "degree", "not positive: D chirality"
            )
            failures.append(improper | {"atoms": atoms})
    return failures


def _bond_failures(bond: PeptideBond, checks: dict) -> list[dict]:
    tests = (
        ("c_n_distance", ("C", "N"), bond.c_n_angstrom, "peptide_bond"),
        ("ca_ca_distance", ("CA", "CA"), bond.ca_ca_angstrom, "ca_ca_trans"),
    )
    failures = []
    for key, names, value, check in tests:
        target, tolerance = checks[f"{check}_angstrom"], checks[f"{check}_tolerance_angstrom"]
        if abs(value - target) > tolerance:
            atoms = [
                atom_ref(bond.first, atom(bond.first, names[0])),
                atom_ref(bond.second, atom(bond.second, names[1])),
            ]
            note = f"outside {target} +/- {tolerance} A (omega {bond.omega_degree:.0f} degrees)"
            failures.append(ev_distance(key, value, atoms, note))
    return failures


def loop_clashes(structure: Structure, gap: Gap, criterion: ClashCriterion) -> list[tuple]:
    """(modelled residue, clash) for the modelled residues' heavy atoms against every other
    heavy atom within reach, except within one residue and between sequence neighbours;
    each pair once. An element without a radius raises only when it is within reach."""
    reach = 2 * max(criterion.radii_angstrom.values())
    loop = list(gap.residue_ids)
    found = []
    for position, rid in enumerate(loop):
        residue = structure.residue(rid)
        others = [
            (other, a)
            for other in structure.residues
            if _counted(residue, other, loop[: position + 1])
            for a in other.atoms
            if a.is_heavy and _within(a, residue.heavy_atoms, reach)
        ]
        found += [
            (residue, c)
            for c in clashes(contacts(residue.heavy_atoms, others, criterion), criterion)
        ]
    return sorted(found, key=lambda p: (-p[1].overlap_angstrom, p[0].id, p[1].new_atom.name))


def _within(atom: Atom, atoms: tuple[Atom, ...], reach: float) -> bool:
    return any(distance(atom.position, other.position) <= reach for other in atoms)


def _counted(residue: Residue, other: Residue, done: list) -> bool:
    if other.id in done:
        return False
    if other.id.chain == residue.id.chain and None not in (other.label_seq, residue.label_seq):
        return abs(other.label_seq - residue.label_seq) > 1
    return True


def criterion(ruleset: RuleSet) -> ClashCriterion:
    rule = ruleset.family(FAMILY)
    matcher = next(r for r in rule if r.rule_id == CONTACTS_RULE).matcher
    source = next(r for r in ruleset.rules if r.rule_id == matcher["criterion_rule"])
    return ClashCriterion.from_matcher(source.matcher)


def gap_findings(structure: Structure, gap: Gap, context: tuple) -> list[Finding]:
    """The geometry finding (if any check fails) and the contacts finding (if any clash)
    of one modelled gap; ``context`` is (ruleset, checks)."""
    ruleset, checks = context
    rules = {r.rule_id: r for r in ruleset.family(FAMILY)}
    findings = []
    failures = geometry_failures(structure, gap, checks)
    if failures:
        findings.append(_finding(structure, gap, rules[GEOMETRY_RULE], failures))
    found = loop_clashes(structure, gap, criterion(ruleset))
    if found:
        findings.append(_finding(structure, gap, rules[CONTACTS_RULE], _clash_items(found)))
    return findings


def _finding(structure: Structure, gap: Gap, rule, items: list[dict]) -> Finding:
    first = structure.residue(gap.residue_ids[0])
    return Finding(
        id=f"{rule.rule_id.replace('.', '/')}/{gap.label}",
        rule=rule,
        title=f"{rule.title}: {gap.label} ({len(items)} item(s))",
        locus=Locus(
            first.id,
            first.name,
            (),
            tuple((rid, structure.residue(rid).name) for rid in gap.residue_ids),
        ),
        anchor_residues=gap.flanks,
        evidence=(ev_label("modelled", gap.label, NO_EXPERIMENTAL_SUPPORT), *items),
    )


def _clash_items(found: list[tuple]) -> list[dict]:
    return [
        ev_distance(
            "loop_clash",
            c.distance_angstrom,
            [atom_ref(residue, c.new_atom), atom_ref(c.residue, c.atom)],
            f"{residue.name} {residue.id.label()} {c.new_atom.name} - {c.residue.name} "
            f"{c.residue.id.label()} {c.atom.name}, overlap {c.overlap_angstrom:.2f} A",
        )
        for residue, c in found
    ]
