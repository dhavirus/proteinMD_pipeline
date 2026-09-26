"""Side-chain builder (TASK-005): NeRF geometry, candidates, and a check against
deposited side chains, which tests the knowledge-base geometry independently of its
recalled source."""

import math
import statistics
from pathlib import Path

import pytest

from simprep.structure.geometry import angle_degree, dihedral_degree, distance, place_atom
from simprep.structure.parse import read_structure
from simprep.variants.build import BuildError, build_side_chain, candidates
from tests.unit import builders as b

PANEL = sorted((Path(__file__).parents[1] / "panel").glob("*.cif.gz"))
CHI_ATOMS = {
    "GLN": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "OE1")),
    "TRP": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")),
}
# Rebuilding a deposited side chain from its own chi angles with ideal geometry differs
# from the deposited coordinates by the refinement's deviations from ideal. Measured on
# the panel in the TASK-005 session: GLN median 0.04 A, max 0.21 A; TRP median 0.16 A,
# max 0.49 A (ring out-of-plane deviations at CG of about 1.5 degrees, amplified along
# the ring). Tolerances sit just above those values.
MEDIAN_RMSD_ANGSTROM = {"GLN": 0.10, "TRP": 0.20}
MAX_RMSD_ANGSTROM = {"GLN": 0.25, "TRP": 0.60}


@pytest.mark.parametrize("torsion", [-179.0, -65.0, 0.0, 62.0, 180.0])
def test_place_atom_reproduces_internal_coordinates(torsion):
    a, b_, c = (0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (2.0, 1.4, 0.3)
    d = place_atom((a, b_, c), (1.53, 111.0, torsion))
    assert distance(c, d) == pytest.approx(1.53)
    assert angle_degree(b_, c, d) == pytest.approx(111.0)
    assert math.cos(math.radians(dihedral_degree(a, b_, c, d) - torsion)) == pytest.approx(1.0)


def test_gln_candidates_include_amide_flips(ruleset):
    data = ruleset.side_chains["GLN"]
    ids = [c.rotamer_id for c in candidates(data)]
    assert len(ids) == 2 * len(data["rotamers"])
    flipped = {c.rotamer_id: c.chi_degree for c in candidates(data)}
    assert flipped["mt-30-flip"][2] == pytest.approx(155.0)


def test_trp_has_no_flips(ruleset):
    assert len(candidates(ruleset.side_chains["TRP"])) == len(
        ruleset.side_chains["TRP"]["rotamers"]
    )


def arg_backbone(**drop):
    atoms = [
        b.atom("N", "N", (1.46, 0, 0)),
        b.atom("CA", "C", (0, 0, 0)),
        b.atom("C", "C", (-0.55, 1.42, 0)),
        b.atom("CB", "C", (-0.53, -0.77, -1.2)),
    ]
    return b.residue("A", 5, "ARG", b.P, [a for a in atoms if a.name not in drop])


def test_builder_refuses_a_residue_without_cb(ruleset):
    with pytest.raises(BuildError, match="lacks CB"):
        build_side_chain(arg_backbone(CB=True), ruleset.side_chains["GLN"], (-65, 180, 0))


def rebuilt_rmsds(ruleset, name: str) -> list[float]:
    data, rmsds = ruleset.side_chains[name], []
    names = [spec["name"] for spec in data["atoms"]]
    for path in PANEL:
        for residue in read_structure(path).residues:
            positions = {a.name: a.position for a in residue.atoms}
            if residue.name != name or residue.altlocs or not set(names) <= positions.keys():
                continue
            chis = tuple(dihedral_degree(*(positions[x] for x in q)) for q in CHI_ATOMS[name])
            built = {a.name: a.position for a in build_side_chain(residue, data, chis)}
            squares = [distance(built[n], positions[n]) ** 2 for n in names]
            rmsds.append(math.sqrt(sum(squares) / len(squares)))
    return rmsds


@pytest.mark.parametrize("name", ["GLN", "TRP"])
def test_rebuilt_deposited_side_chains_match_the_deposited_coordinates(ruleset, name):
    rmsds = rebuilt_rmsds(ruleset, name)
    assert len(rmsds) >= 40, "the panel should provide enough deposited examples"
    assert statistics.median(rmsds) <= MEDIAN_RMSD_ANGSTROM[name]
    assert max(rmsds) <= MAX_RMSD_ANGSTROM[name]
