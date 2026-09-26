"""Covalent attachments the file does not annotate, found by declared attachment chemistry.

Only residue pairs without an annotated covalent link are examined, only altloc-compatible
atoms are compared (same altloc, or either has none), and only atom pairs that match a
pattern in the knowledge base count.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from simprep.config import AuditConfig
from simprep.detectors import register
from simprep.findings import Finding, Locus, atom_ref, ev_distance, ev_flag, ev_label, ev_number
from simprep.rules import Rule, RuleSet
from simprep.structure.geometry import distance
from simprep.structure.model import Atom, Residue, ResidueClass, ResidueId, Structure

WILDCARD = "*"
LIGAND_CLASSES = (ResidueClass.NONPOLYMER, ResidueClass.BRANCHED)
SEQUON_ACCEPTORS = frozenset({"SER", "THR"})
SEQUON_BLOCKER = "PRO"


@dataclass(frozen=True)
class Contact:
    """Closest matching atom pair between one ligand residue and one polymer residue."""

    pattern: dict
    ligand: Residue
    ligand_atom: Atom
    polymer: Residue
    polymer_atom: Atom
    distance_angstrom: float


@register("covalent_contacts")
def detect_covalent_contacts(
    structure: Structure, ruleset: RuleSet, config: AuditConfig
) -> list[Finding]:
    """One finding per (ligand residue, polymer residue, pattern) with a matching contact."""
    findings = []
    for rule in ruleset.family("covalent_contacts"):
        linked = annotated_pairs(structure, rule.matcher["annotated_conn_types"])
        for pattern in rule.matcher["patterns"]:
            for contact in pattern_contacts(structure, pattern, linked):
                findings.append(_finding(structure, contact, rule))
    return findings


def annotated_pairs(structure: Structure, conn_types: list[str]) -> set[frozenset[ResidueId]]:
    return {
        frozenset((link.partner1.residue, link.partner2.residue))
        for link in structure.links
        if link.conn_type in conn_types
    }


def pattern_contacts(
    structure: Structure, pattern: dict, linked: set[frozenset[ResidueId]]
) -> list[Contact]:
    """Closest matching, altloc-compatible atom pair per unlinked residue pair."""
    ligands = [r for r in structure.residues if _ligand_matches(r, pattern["ligand"])]
    polymers = [r for r in structure.residues if _polymer_matches(r, pattern["polymer"])]
    contacts = []
    for ligand in ligands:
        for polymer in polymers:
            if frozenset((ligand.id, polymer.id)) in linked:
                continue
            contact = _closest(pattern, ligand, polymer)
            if (
                contact is not None
                and contact.distance_angstrom <= pattern["max_distance_angstrom"]
            ):
                contacts.append(contact)
    return contacts


def _ligand_matches(residue: Residue, selector: dict) -> bool:
    comp_ids = selector["comp_ids"]
    return residue.residue_class in LIGAND_CLASSES and (
        WILDCARD in comp_ids or residue.name in comp_ids
    )


def _polymer_matches(residue: Residue, selector: dict) -> bool:
    return residue.residue_class is ResidueClass.POLYMER and residue.name in selector["res_names"]


def _atom_matches(atom: Atom, selector: dict) -> bool:
    if "atom_names" in selector and atom.name not in selector["atom_names"]:
        return False
    return "elements" not in selector or atom.element in selector["elements"]


def compatible(a: Atom, b: Atom) -> bool:
    """Atoms can coexist: same altloc, or at least one has none."""
    return not a.altloc or not b.altloc or a.altloc == b.altloc


def _closest(pattern: dict, ligand: Residue, polymer: Residue) -> Contact | None:
    pairs = [
        (distance(a.position, b.position), a, b)
        for a in ligand.heavy_atoms
        if _atom_matches(a, pattern["ligand"])
        for b in polymer.heavy_atoms
        if _atom_matches(b, pattern["polymer"]) and compatible(a, b)
    ]
    if not pairs:
        return None
    d, a, b = min(pairs, key=lambda pair: (pair[0], pair[1].name, pair[1].altloc, pair[2].altloc))
    return Contact(pattern, ligand, a, polymer, b, d)


def sequon(structure: Structure, residue: Residue) -> tuple[str, ...] | None:
    """The residue and its two sequence successors, or None if they cannot be placed."""
    sequence = structure.sequence(residue.id.chain)
    if residue.label_seq is None or residue.label_seq + 2 > len(sequence):
        return None
    return sequence[residue.label_seq - 1 : residue.label_seq + 2]


def in_sequon(triplet: tuple[str, ...]) -> bool:
    """N-X-S/T with X not proline; microheterogeneous positions accept any alternative."""
    first, middle, last = (set(position.split(",")) for position in triplet)
    return "ASN" in first and SEQUON_BLOCKER not in middle and bool(last & SEQUON_ACCEPTORS)


def _finding(structure: Structure, contact: Contact, rule: Rule) -> Finding:
    ligand, polymer = contact.ligand, contact.polymer
    pair = f"{ligand.id.label()}-{polymer.id.label()}"
    return Finding(
        id=f"covalent_contacts/{contact.pattern['name']}/{pair}",
        rule=rule,
        title=f"{rule.title}: {ligand.name} {ligand.id.label()} {contact.ligand_atom.name} - "
        f"{polymer.name} {polymer.id.label()} {contact.polymer_atom.name} "
        f"({contact.distance_angstrom:.2f} A, {contact.pattern['name']})",
        locus=Locus(
            ligand.id,
            ligand.name,
            (contact.ligand_atom.name,),
            ((ligand.id, ligand.name), (polymer.id, polymer.name)),
        ),
        anchor_residues=(ligand.id, polymer.id),
        evidence=tuple(_evidence(structure, contact)),
    )


def _evidence(structure: Structure, contact: Contact) -> list[dict]:
    pattern = contact.pattern
    items = [
        ev_label("pattern", pattern["name"]),
        ev_label("bond", pattern["bond"]),
        ev_distance(
            "contact_distance",
            contact.distance_angstrom,
            [
                atom_ref(contact.ligand, contact.ligand_atom),
                atom_ref(contact.polymer, contact.polymer_atom),
            ],
        ),
        ev_number("max_distance", pattern["max_distance_angstrom"], "angstrom"),
        ev_label(
            "struct_conn_link",
            "none between these residues",
            note="checked against the annotated covalent connection types of the rule",
        ),
    ]
    if pattern["sequon_check"]:
        items += _sequon_evidence(structure, contact.polymer)
    return items


def _sequon_evidence(structure: Structure, residue: Residue) -> list[dict]:
    triplet = sequon(structure, residue)
    if triplet is None:
        return [ev_label("sequon", "not determinable: no entity sequence for this position")]
    return [
        ev_label("sequon", "-".join(triplet), note="residue and its two successors"),
        ev_flag("in_sequon", in_sequon(triplet), note="N-X-S/T with X not proline"),
    ]


def link_to_groups(findings: list[Finding]) -> list[Finding]:
    """Name, on each candidate, the unrecognized-group finding its ligand belongs to."""
    group_of = {
        residue_id: finding.id
        for finding in findings
        if finding.rule.family == "unrecognized"
        for residue_id in finding.claims
    }
    return [_with_group(finding, group_of) for finding in findings]


def _with_group(finding: Finding, group_of: dict[ResidueId, str]) -> Finding:
    if finding.rule.family != "covalent_contacts":
        return finding
    group = group_of.get(finding.anchor_residues[0])
    if group is None:
        return finding
    return replace(finding, evidence=(*finding.evidence, ev_label("ligand_group_finding", group)))
