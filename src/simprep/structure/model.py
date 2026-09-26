"""Immutable, parser-independent structure model.

Detectors consume only this model, never gemmi objects, so they stay pure functions
and can be unit-tested on hand-built structures.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from functools import cached_property

HYDROGEN_ELEMENTS = frozenset({"H", "D"})


@dataclass(frozen=True, order=True)
class ResidueId:
    """Author-numbered residue address: chain, sequence number, insertion code."""

    chain: str
    seq_num: int
    ins_code: str = ""

    def label(self) -> str:
        """Compact label such as ``A:123`` or ``A:123B``."""
        return f"{self.chain}:{self.seq_num}{self.ins_code}"

    def to_dict(self) -> dict:
        return {"chain": self.chain, "seq_num": self.seq_num, "ins_code": self.ins_code}


class ResidueClass(StrEnum):
    """Entity class of a residue, taken from the file's entity records."""

    POLYMER = "polymer"
    NONPOLYMER = "nonpolymer"
    BRANCHED = "branched"
    WATER = "water"


@dataclass(frozen=True)
class Atom:
    """One atom record (one altloc of one atom)."""

    name: str
    element: str
    position: tuple[float, float, float]
    occupancy: float
    b_iso: float
    altloc: str = ""

    @property
    def is_heavy(self) -> bool:
        return self.element.upper() not in HYDROGEN_ELEMENTS


@dataclass(frozen=True)
class Residue:
    """A residue instance with all of its atom records (every altloc)."""

    id: ResidueId
    name: str
    residue_class: ResidueClass
    atoms: tuple[Atom, ...]
    label_seq: int | None = None

    @property
    def heavy_atoms(self) -> tuple[Atom, ...]:
        return tuple(atom for atom in self.atoms if atom.is_heavy)

    @property
    def altlocs(self) -> tuple[str, ...]:
        return tuple(sorted({atom.altloc for atom in self.atoms if atom.altloc}))


@dataclass(frozen=True)
class LinkPartner:
    """One end of an annotated connection (struct_conn / LINK / SSBOND)."""

    residue: ResidueId
    res_name: str
    atom_name: str
    altloc: str = ""


@dataclass(frozen=True)
class Link:
    """An annotated connection record, kept verbatim as evidence."""

    conn_id: str
    conn_type: str
    partner1: LinkPartner
    partner2: LinkPartner
    distance_angstrom: float | None = None


@dataclass(frozen=True)
class ModifiedResidue:
    """A modified-residue annotation (pdbx_struct_mod_residue / MODRES)."""

    residue: ResidueId
    res_name: str
    parent_res_name: str | None
    details: str | None = None


@dataclass(frozen=True)
class UnobservedResidue:
    """A residue listed as unobserved (pdbx_unobs_or_zero_occ_residues / REMARK 465)."""

    residue: ResidueId
    res_name: str
    label_seq: int | None
    is_polymer: bool


@dataclass(frozen=True)
class Structure:
    """A single model of an asymmetric unit plus the annotations detectors rely on.

    ``annotation_categories`` lists which annotation sources were present in the input,
    so that absence of a record can be distinguished from absence of the annotation.
    ``crystallization_details`` is the free-text crystallization condition, kept verbatim
    as evidence (e.g. whether a bound ion was a buffer component).
    ``polymer_sequences`` maps each chain to its full entity sequence (observed and
    unobserved residues, indexed by label_seq - 1); a position with microheterogeneity
    lists its alternatives joined by commas.
    """

    name: str
    residues: tuple[Residue, ...]
    links: tuple[Link, ...] = ()
    modified_residues: tuple[ModifiedResidue, ...] = ()
    unobserved_residues: tuple[UnobservedResidue, ...] = ()
    annotation_categories: frozenset[str] = field(default_factory=frozenset)
    crystallization_details: str | None = None
    polymer_sequences: tuple[tuple[str, tuple[str, ...]], ...] = ()

    @cached_property
    def residue_index(self) -> dict[ResidueId, Residue]:
        return {residue.id: residue for residue in self.residues}

    def residue(self, residue_id: ResidueId) -> Residue:
        """Return the residue with ``residue_id``; raise KeyError naming it if absent."""
        try:
            return self.residue_index[residue_id]
        except KeyError:
            raise KeyError(f"residue {residue_id.label()} not found in {self.name}") from None

    def sequence(self, chain: str) -> tuple[str, ...]:
        """Full entity sequence of ``chain`` (empty if the file has none)."""
        return dict(self.polymer_sequences).get(chain, ())

    def polymer_residues(self, chain: str) -> tuple[Residue, ...]:
        """Observed polymer residues of ``chain`` in file order."""
        return tuple(
            residue
            for residue in self.residues
            if residue.id.chain == chain and residue.residue_class is ResidueClass.POLYMER
        )
