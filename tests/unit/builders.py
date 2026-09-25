"""Hand-built minimal structures for detector unit tests (no parser involved)."""

from __future__ import annotations

import math

from simprep.config import AuditConfig, Region, Threshold
from simprep.structure.model import (
    Atom,
    Link,
    LinkPartner,
    Residue,
    ResidueClass,
    ResidueId,
    Structure,
)

P = ResidueClass.POLYMER
NP = ResidueClass.NONPOLYMER
BR = ResidueClass.BRANCHED
W = ResidueClass.WATER
TETRAHEDRAL = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
THRESHOLDS = (Threshold(6.0, 2), Threshold(12.0, 1))


def atom(name, element, position, occupancy=1.0, b_iso=20.0, altloc=""):
    return Atom(name, element, tuple(float(x) for x in position), occupancy, b_iso, altloc)


def residue(chain, num, name, cls, atoms, label_seq=None):
    return Residue(ResidueId(chain, num), name, cls, tuple(atoms), label_seq)


def scaled(direction, length):
    norm = math.sqrt(sum(x * x for x in direction))
    return tuple(length * x / norm for x in direction)


def offset(position, delta):
    return tuple(p + d for p, d in zip(position, delta, strict=True))


def ion(chain, num, element, position=(0, 0, 0), altloc="", occupancy=1.0):
    return residue(
        chain, num, element, NP, [atom(element, element, position, occupancy, altloc=altloc)]
    )


def water(chain, num, position, altloc="", occupancy=1.0):
    return residue(chain, num, "HOH", W, [atom("O", "O", position, occupancy, altloc=altloc)])


def his(chain, num, ne2_position, label_seq=None):
    """Histidine reduced to the atoms tests need; NE2 at ``ne2_position``."""
    away = scaled(ne2_position, 1.3)
    return residue(
        chain,
        num,
        "HIS",
        P,
        [
            atom("CA", "C", offset(ne2_position, scaled(away, 4.0))),
            atom("NE2", "N", ne2_position),
        ],
        label_seq,
    )


def cys(chain, num, sg_position, label_seq=None):
    return residue(
        chain,
        num,
        "CYS",
        P,
        [
            atom("CA", "C", offset(sg_position, scaled(sg_position, 2.5))),
            atom("SG", "S", sg_position),
        ],
        label_seq,
    )


def structure(residues, links=(), **annotations):
    return Structure("TEST", tuple(residues), tuple(links), **annotations)


def link(conn_id, conn_type, partners, distance=None):
    (r1, n1, a1), (r2, n2, a2) = partners
    return Link(conn_id, conn_type, LinkPartner(r1, n1, a1), LinkPartner(r2, n2, a2), distance)


def zinc_catalytic_site():
    """Zn with three His NE2 and one water, tetrahedral at 2.05 A."""
    ligand_positions = [scaled(d, 2.05) for d in TETRAHEDRAL]
    residues = [ion("A", 900, "ZN")]
    residues += [his("A", 10 + i, ligand_positions[i]) for i in range(3)]
    residues.append(water("A", 1000, ligand_positions[3]))
    return structure(residues)


def zinc_structural_site():
    """Zn with four Cys SG, tetrahedral at 2.33 A."""
    residues = [ion("A", 900, "ZN")]
    residues += [cys("A", 20 + i, scaled(d, 2.33)) for i, d in enumerate(TETRAHEDRAL)]
    return structure(residues)


def config(regions=()):
    return AuditConfig(tuple(regions), THRESHOLDS)


def region(name, residue_ids):
    return Region(name, f"test region {name}", frozenset(residue_ids))
