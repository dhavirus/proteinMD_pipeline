"""Which atoms a relaxation moves, and what the force field can describe (pure; TASK-006).

The mobile shell is every residue or water with a heavy atom within
``mobile_shell_angstrom`` of a mutated residue; side chains (not N, CA, C, O, OXT) of
shell residues and whole shell waters move, and the mutated residues' side chains are
left unrestrained. Residues without a force-field template are left out of the
calculation. So that no mobile atom feels their absence, a shell residue with a mobile
atom within the non-bonded cutoff of one is held fixed (frozen, and recorded); if a
mutated residue's own side chain is that close, relaxation stops (unknown chemistry is a
hard stop).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from simprep.structure.model import Residue, ResidueClass, ResidueId, Structure

BACKBONE = frozenset({"N", "CA", "C", "O", "OXT"})
WATER = "HOH"
# Residue names the Amber14 protein (ff14SB) and TIP3P templates describe; OpenMM picks
# the protonation variant (e.g. HID/HIE, CYX) itself.
TEMPLATED = frozenset(
    {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        WATER,
    }
)
PEPTIDE_BOND_MAX_ANGSTROM = 2.0
ANGSTROM_PER_NM = 10.0

AtomKey = tuple[ResidueId, str]


class RelaxError(ValueError):
    """Relaxation cannot run as specified; the message says why and what to do."""


@dataclass(frozen=True)
class Shell:
    sites: frozenset[ResidueId]
    residues: frozenset[ResidueId]
    mobile: frozenset[AtomKey]
    unrestrained: frozenset[AtomKey]
    frozen: tuple[tuple[ResidueId, str, float], ...] = ()  # (residue, near, distance A)


def is_templated(residue: Residue) -> bool:
    """A standard polymer residue or a water: the force field has a template for it."""
    if residue.name == WATER:
        return True
    return residue.name in TEMPLATED and residue.residue_class is ResidueClass.POLYMER


def heavy(residue: Residue):
    return [atom for atom in residue.atoms if atom.is_heavy]


def shell_residues(structure: Structure, sites: frozenset[ResidueId], radius: float) -> frozenset:
    """Sites plus every residue with a heavy atom within ``radius`` A of a site heavy atom."""
    points = [a.position for rid in sites for a in heavy(structure.residue(rid))]
    return frozenset(
        residue.id
        for residue in structure.residues
        if residue.id in sites
        or any(math.dist(a.position, p) <= radius for a in heavy(residue) for p in points)
    )


def build_shell(
    structure: Structure, sites: frozenset[ResidueId], residues: frozenset, cutoff_nm: float
) -> Shell:
    """Mobile atoms of ``residues`` (a shell, possibly a union of shells), with shell
    residues near template-less chemistry held fixed."""
    near = _near_untemplated(structure, residues - sites, cutoff_nm * ANGSTROM_PER_NM)
    mobile, unrestrained = set(), set()
    for rid in sorted(residues - {rid for rid, _, _ in near}):
        residue = structure.residue(rid)
        for atom in _movable(residue):
            mobile.add((rid, atom.name))
            if rid in sites:
                unrestrained.add((rid, atom.name))
    return Shell(sites, residues, frozenset(mobile), frozenset(unrestrained), tuple(near))


def _movable(residue: Residue):
    return [a for a in heavy(residue) if residue.name == WATER or a.name not in BACKBONE]


def _near_untemplated(structure: Structure, candidates: frozenset, limit: float) -> list:
    """(residue id, nearest template-less residue, distance) for each candidate residue
    with a movable atom within ``limit`` A of a template-less residue."""
    others = [r for r in structure.residues if not is_templated(r)]
    near = []
    for rid in sorted(candidates):
        movable = _movable(structure.residue(rid))
        pairs = [
            (math.dist(a.position, b.position), o)
            for a in movable
            for o in others
            for b in heavy(o)
        ]
        if pairs and min(pairs, key=lambda p: p[0])[0] <= limit:
            distance, other = min(pairs, key=lambda p: p[0])
            near.append((rid, f"{other.name} {other.id.label()}", round(distance, 2)))
    return near


def untemplated(structure: Structure, shell: Shell, cutoff_nm: float) -> tuple[list, list]:
    """(left out, blocking): residues without a template, each with its closest heavy-atom
    distance to a mobile atom; blocking ones lie within the non-bonded cutoff."""
    index = structure.residue_index
    mobile = [a.position for rid, name in shell.mobile for a in heavy(index[rid]) if a.name == name]
    left_out, blocking = [], []
    for residue in structure.residues:
        if is_templated(residue):
            continue
        distance = min(math.dist(a.position, p) for a in heavy(residue) for p in mobile)
        entry = (residue, round(distance, 2))
        (blocking if distance <= cutoff_nm * ANGSTROM_PER_NM else left_out).append(entry)
    return left_out, blocking


def segments(residues: list[Residue]) -> list[list[Residue]]:
    """Polymer residues split wherever the peptide bond is missing (different chain, gap,
    left-out residue, or C-N farther than PEPTIDE_BOND_MAX_ANGSTROM)."""
    result: list[list[Residue]] = []
    previous = None
    for residue in residues:
        if previous is None or not _bonded(previous, residue):
            result.append([])
        result[-1].append(residue)
        previous = residue
    return result


def _bonded(first: Residue, second: Residue) -> bool:
    if first.id.chain != second.id.chain:
        return False
    carbon = next((a for a in first.atoms if a.name == "C"), None)
    nitrogen = next((a for a in second.atoms if a.name == "N"), None)
    if carbon is None or nitrogen is None:
        return False
    return math.dist(carbon.position, nitrogen.position) <= PEPTIDE_BOND_MAX_ANGSTROM


def polymer_and_water(structure: Structure) -> tuple[list[Residue], list[Residue]]:
    """Templated residues to relax with: polymer residues in file order, then waters."""
    polymer = [
        r
        for r in structure.residues
        if r.name in TEMPLATED and r.residue_class is ResidueClass.POLYMER
    ]
    waters = [r for r in structure.residues if r.name == WATER]
    return polymer, waters
