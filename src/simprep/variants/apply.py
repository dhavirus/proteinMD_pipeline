"""Turn the prepared wild type into a variant (pure; TASK-005).

The mutated residue keeps its heavy backbone atoms and CB, loses every hydrogen
(decision 2) and every other side-chain atom, and gains the built side chain. The entity
sequence is updated at that position; connection records that referenced a removed atom
are dropped. Nothing else changes.
"""

from __future__ import annotations

from dataclasses import replace

from simprep.structure.model import Atom, Link, Residue, Structure
from simprep.variants.build import BACKBONE_KEPT

TERMINAL_OXYGEN = "OXT"


KEPT_NAMES = (*BACKBONE_KEPT, TERMINAL_OXYGEN)


def mutated_residue(residue: Residue, to: str, side_chain: tuple[Atom, ...]) -> Residue:
    kept = tuple(atom for atom in residue.atoms if atom.is_heavy and atom.name in KEPT_NAMES)
    return replace(residue, name=to, atoms=kept + side_chain)


def apply_mutations(wild_type: Structure, residues: tuple[Residue, ...]) -> Structure:
    """``wild_type`` with each of ``residues`` (already mutated) put in its place."""
    new = {residue.id: residue for residue in residues}
    kept = tuple(new.get(residue.id, residue) for residue in wild_type.residues)
    return replace(
        wild_type,
        residues=kept,
        links=tuple(link for link in wild_type.links if _link_kept(link, new)),
        polymer_sequences=_sequences(wild_type, residues),
    )


def _link_kept(link: Link, new: dict) -> bool:
    """A link to a mutated residue survives only if it names a kept backbone atom (a
    reused name on the new side chain is a different atom)."""
    return all(
        partner.residue not in new or partner.atom_name in KEPT_NAMES
        for partner in (link.partner1, link.partner2)
    )


def _sequences(wild_type: Structure, residues: tuple[Residue, ...]) -> tuple:
    sequences = {chain: list(sequence) for chain, sequence in wild_type.polymer_sequences}
    for residue in residues:
        sequence = sequences.get(residue.id.chain)
        if sequence and residue.label_seq and residue.label_seq <= len(sequence):
            sequence[residue.label_seq - 1] = residue.name
    return tuple((chain, tuple(sequence)) for chain, sequence in sequences.items())
