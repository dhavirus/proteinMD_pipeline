"""Build a side chain from the kept backbone and a chi set (pure; TASK-005).

Data (internal coordinates, rotamers) come from ``knowledge/side_chains.yaml``; nothing
here is residue-specific.
"""

from __future__ import annotations

from dataclasses import dataclass

from simprep.structure.geometry import place_atom
from simprep.structure.model import Atom, Residue

BACKBONE_KEPT = ("N", "CA", "C", "O", "CB")
FLIP_DEGREE = 180.0
FLIP_SUFFIX = "-flip"


class BuildError(ValueError):
    """The side chain cannot be built on this residue (missing data or atoms)."""


@dataclass(frozen=True)
class Candidate:
    """One side-chain placement to evaluate: a rotamer, possibly with its flip."""

    rotamer_id: str
    chi_degree: tuple[float, ...]
    frequency_percent: float


def candidates(data: dict) -> tuple[Candidate, ...]:
    """Every rotamer of ``data``, plus its planar-group flip when ``flip_chi`` is set."""
    result = []
    for rotamer in data["rotamers"]:
        chis = tuple(rotamer["chi_degree"])
        result.append(Candidate(rotamer["id"], chis, rotamer["frequency_percent"]))
        if "flip_chi" in data:
            flipped = list(chis)
            flipped[data["flip_chi"] - 1] = _wrap(flipped[data["flip_chi"] - 1] + FLIP_DEGREE)
            result.append(
                Candidate(rotamer["id"] + FLIP_SUFFIX, tuple(flipped), rotamer["frequency_percent"])
            )
    return tuple(result)


def _wrap(angle_degree: float) -> float:
    return (angle_degree + 180.0) % 360.0 - 180.0


def build_side_chain(
    residue: Residue, data: dict, chi_degree: tuple[float, ...]
) -> tuple[Atom, ...]:
    """Side-chain atoms beyond CB placed on ``residue``'s N, CA, CB (occupancy 1,
    B-factor of the residue's CB, no altloc)."""
    positions = _backbone(residue)
    cb = next(a for a in residue.atoms if a.name == "CB")
    atoms = []
    for spec in data["atoms"]:
        refs = tuple(positions[name] for name in spec["refs"])
        internal = (spec["bond_angstrom"], spec["angle_degree"], _dihedral(spec, chi_degree))
        positions[spec["name"]] = place_atom(refs, internal)
        atoms.append(Atom(spec["name"], spec["element"], positions[spec["name"]], 1.0, cb.b_iso))
    return tuple(atoms)


def _dihedral(spec: dict, chi_degree: tuple[float, ...]) -> float:
    dihedral = spec["dihedral"]
    if isinstance(dihedral, int | float):
        return float(dihedral)
    return chi_degree[dihedral["chi"] - 1] + dihedral.get("offset_degree", 0.0)


def _backbone(residue: Residue) -> dict[str, tuple[float, float, float]]:
    """Heavy backbone positions of a residue without altlocs; raise naming what is missing."""
    if residue.altlocs:
        raise BuildError(
            f"{residue.id.label()} still has altlocs {residue.altlocs}; decide them first"
        )
    positions = {a.name: a.position for a in residue.atoms if a.name in BACKBONE_KEPT}
    missing = [name for name in ("N", "CA", "C", "CB") if name not in positions]
    if missing:
        raise BuildError(f"{residue.id.label()} {residue.name} lacks {', '.join(missing)}")
    return positions
