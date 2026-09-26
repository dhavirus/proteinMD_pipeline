"""Contacts of a built side chain with the rest of the system (pure; TASK-005).

A heavy-atom pair clashes when its overlap (sum of van der Waals radii minus distance)
reaches the rule's ``overlap_angstrom``; a pair of polar atoms (possible hydrogen bond)
clashes only below ``polar_min_distance_angstrom``. Hydrogens and the mutated residue's
own atoms are ignored; altloc atoms of neighbours are each checked (conservative).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from simprep.structure.model import Atom, Residue, ResidueId, Structure

NEIGHBOUR_MARGIN_ANGSTROM = 12.0  # search radius around CB; > side-chain reach + 2 radii


class ClashDataError(ValueError):
    """A neighbouring atom's element has no van der Waals radius in the rule data."""


@dataclass(frozen=True)
class ClashCriterion:
    overlap_angstrom: float
    radii_angstrom: dict[str, float]
    polar_elements: frozenset[str]
    polar_min_distance_angstrom: float

    @classmethod
    def from_matcher(cls, matcher: dict) -> ClashCriterion:
        return cls(
            matcher["overlap_angstrom"],
            dict(matcher["radii_angstrom"]),
            frozenset(matcher["polar_elements"]),
            matcher["polar_min_distance_angstrom"],
        )

    def is_clash(self, contact: Contact) -> bool:
        polar = {contact.new_atom.element, contact.atom.element} <= self.polar_elements
        if polar:
            return contact.distance_angstrom < self.polar_min_distance_angstrom
        return contact.overlap_angstrom >= self.overlap_angstrom


@dataclass(frozen=True)
class Contact:
    """A new side-chain atom and a neighbouring atom, with distance and vdW overlap."""

    new_atom: Atom
    residue: Residue
    atom: Atom
    distance_angstrom: float
    overlap_angstrom: float


def neighbours(structure: Structure, residue_id: ResidueId) -> list[tuple[Residue, Atom]]:
    """Heavy atoms of other residues within the search radius of ``residue_id``'s CB."""
    residue = structure.residue(residue_id)
    cb = next(a.position for a in residue.atoms if a.name == "CB")
    return [
        (other, atom)
        for other in structure.residues
        if other.id != residue_id
        for atom in other.atoms
        if atom.is_heavy and math.dist(atom.position, cb) <= NEIGHBOUR_MARGIN_ANGSTROM
    ]


def contacts(new_atoms: tuple[Atom, ...], nearby: list, criterion: ClashCriterion) -> list[Contact]:
    """Every new-atom / neighbour pair with positive overlap, largest overlap first."""
    found = []
    for new in new_atoms:
        for residue, atom in nearby:
            distance = math.dist(new.position, atom.position)
            overlap = _radius(new, criterion) + _radius(atom, criterion) - distance
            if overlap > 0:
                found.append(Contact(new, residue, atom, distance, overlap))
    return sorted(found, key=lambda c: (-c.overlap_angstrom, c.residue.id, c.atom.name))


def clashes(found: list[Contact], criterion: ClashCriterion) -> list[Contact]:
    return [c for c in found if criterion.is_clash(c)]


def _radius(atom: Atom, criterion: ClashCriterion) -> float:
    try:
        return criterion.radii_angstrom[atom.element]
    except KeyError:
        raise ClashDataError(
            f"no van der Waals radius for element {atom.element} (atom {atom.name}); add it "
            "to knowledge/variant_build.yaml with its source"
        ) from None
