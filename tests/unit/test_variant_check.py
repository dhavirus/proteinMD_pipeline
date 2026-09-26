"""Clash criterion (TASK-005 decision 3, knowledge/variant_build.yaml): unit cases and a
calibration against the deposited panel."""

from pathlib import Path

import pytest

from simprep.structure.model import ResidueId
from simprep.structure.parse import read_structure
from simprep.variants.check import (
    ClashCriterion,
    ClashDataError,
    clashes,
    contacts,
    neighbours,
)
from tests.unit import builders as b

PANEL = sorted((Path(__file__).parents[1] / "panel").glob("*.cif.gz"))
BACKBONE = {"N", "CA", "C", "O", "CB"}
POLAR_RESIDUES = ("ARG", "ASN", "ASP", "GLN", "GLU", "HIS", "LYS", "SER", "THR", "TYR")


@pytest.fixture
def criterion(ruleset):
    return ClashCriterion.from_matcher(ruleset.family("variant_build")[0].matcher)


def pair(new, other, distance_angstrom):
    """A new atom at the origin and a neighbour atom on the x axis."""
    (name1, element1), (name2, element2) = new, other
    atom = b.atom(name1, element1, (0, 0, 0))
    neighbour = b.residue("A", 9, "XXX", b.P, [b.atom(name2, element2, (distance_angstrom, 0, 0))])
    return (atom,), [(neighbour, neighbour.atoms[0])]


def found_clashes(criterion, new, other, distance_angstrom):
    atoms, nearby = pair(new, other, distance_angstrom)
    return clashes(contacts(atoms, nearby, criterion), criterion)


def test_nonpolar_pair_clashes_at_the_overlap_limit(criterion):
    # C-C: radii 1.70 + 1.70; overlap 0.4 A at 3.0 A.
    assert found_clashes(criterion, ("CG", "C"), ("CB", "C"), 2.99)
    assert not found_clashes(criterion, ("CG", "C"), ("CB", "C"), 3.01)


def test_polar_pair_clashes_only_below_the_hydrogen_bond_limit(criterion):
    assert not found_clashes(criterion, ("OE1", "O"), ("O", "O"), 2.35)
    assert found_clashes(criterion, ("OE1", "O"), ("O", "O"), 2.25)


def test_unknown_element_is_a_hard_stop(criterion):
    atoms, nearby = pair(("CG", "C"), ("ZN", "ZN"), 2.0)
    with pytest.raises(ClashDataError, match="no van der Waals radius for element ZN"):
        contacts(atoms, nearby, criterion)


def test_neighbours_skip_hydrogens_and_the_residue_itself():
    site = b.residue(
        "A", 1, "ARG", b.P, [b.atom("CB", "C", (0, 0, 0)), b.atom("CG", "C", (1, 0, 0))]
    )
    other = b.residue(
        "A", 2, "SER", b.P, [b.atom("OG", "O", (2, 0, 0)), b.atom("HG", "H", (2.5, 0, 0))]
    )
    far = b.residue("A", 3, "SER", b.P, [b.atom("OG", "O", (40, 0, 0))])
    nearby = neighbours(b.structure([site, other, far]), ResidueId("A", 1))
    assert [(r.id.seq_num, a.name) for r, a in nearby] == [(2, "OG")]


def deposited_side_chain_contacts(criterion, names):
    for path in PANEL:
        structure = read_structure(path)
        for residue in structure.residues:
            if residue.name not in names or residue.altlocs:
                continue
            side = tuple(a for a in residue.atoms if a.is_heavy and a.name not in BACKBONE)
            # Metal ions have no radius in the rule (a mutation next to one is refused);
            # the calibration leaves them out.
            nearby = [
                (r, a)
                for r, a in neighbours(structure, residue.id)
                if a.element in criterion.radii_angstrom
            ]
            yield residue, contacts(side, nearby, criterion)


def test_polar_limit_is_below_almost_every_deposited_polar_contact(criterion):
    """Recomputes the rule's empirical basis: of the deposited polar side-chain contacts
    with positive overlap, at most 0.2 % are shorter than polar_min_distance_angstrom."""
    distances = [
        c.distance_angstrom
        for _, found in deposited_side_chain_contacts(criterion, POLAR_RESIDUES)
        for c in found
        if {c.new_atom.element, c.atom.element} <= criterion.polar_elements
    ]
    assert len(distances) > 1000
    shorter = sum(d < criterion.polar_min_distance_angstrom for d in distances)
    assert shorter / len(distances) <= 0.002


def test_deposited_gln_trp_arg_side_chains_rarely_clash(criterion):
    """A criterion that flags deposited side chains would reject good placements. Measured
    in the TASK-005 session: 11 of 286 deposited GLN/TRP/ARG side chains (3.8 %), mostly
    C-O contacts of 2.7-2.8 A (heavy-atom radii, no hydrogens); the limit is 5 %."""
    flagged = [
        bool(clashes(found, criterion))
        for _, found in deposited_side_chain_contacts(criterion, ("GLN", "TRP", "ARG"))
    ]
    assert len(flagged) > 200
    assert sum(flagged) / len(flagged) <= 0.05
