"""Parameter data, pure checks (TASK-010): the 12-6-4 C4 terms against ParmEd's add12_6_4
reference values (tests/panel/parameter_fixtures/lj1264_reference.yaml, computed once with
AmberTools), and the generated DDZ template against its recipe."""

import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
import yaml

from simprep.parameterize.lj1264 import A4_TO_NM4, KCAL_TO_KJ, c4_by_class, sigma_epsilon

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = yaml.safe_load(
    (ROOT / "tests/panel/parameter_fixtures/lj1264_reference.yaml").read_text()
)
DDZ = ET.parse(ROOT / "knowledge/forcefield/ddz.xml")


def test_c4_terms_match_parmed(ruleset):
    data = ruleset.parameterization["lj1264"]
    ours = c4_by_class(data, "CA")
    by_amber = {
        p["amber_type"]: ours[c] / (KCAL_TO_KJ * A4_TO_NM4)
        for c, p in data["polarizabilities"].items()
    }
    for amber_type, value in REFERENCE["c4_kcal_per_mol_a4"].items():
        if amber_type == "Ca2+":
            continue  # ion-ion pairs are not modelled (one ion per model in v0.1)
        assert by_amber[amber_type] == pytest.approx(value, abs=1e-5), amber_type


def test_ion_12_6_terms_convert_from_rmin_half(ruleset):
    sigma, epsilon = sigma_epsilon(ruleset.parameterization["lj1264"], "CA")
    assert sigma == pytest.approx(2 * 0.1642 / 2 ** (1 / 6))  # Rmin/2 1.642 A (frcmod file)
    assert epsilon == pytest.approx(0.10185975 * 4.184)


def test_ddz_is_neutral_with_the_serine_backbone():
    from openmm import app

    serine = next(
        r
        for r in ET.parse(Path(app.__file__).parent / "data/amber14/protein.ff14SB.xml").iter(
            "Residue"
        )
        if r.get("name") == "SER"
    )
    ser = {a.get("name"): (a.get("type"), float(a.get("charge"))) for a in serine.iter("Atom")}
    ddz = {
        a.get("name"): (a.get("type"), float(a.get("charge")))
        for a in DDZ.iter("Atom")
        if a.get("charge")
    }
    assert sum(q for _, q in ddz.values()) == pytest.approx(0.0, abs=1e-5)
    for name in ("N", "H", "CA", "HA", "C", "O"):
        assert ddz[name] == ser[name]
    assert {ddz[n][0] for n in ("CB", "HB", "OG1", "OG2", "HG1", "HG2")} <= {
        "ddz-CB",
        "ddz-HB",
        "ddz-OG",
        "ddz-HG",
    }


def test_ddz_has_the_gem_diol_angle_ff14sb_lacks():
    angles = {(a.get("type1"), a.get("type2"), a.get("type3")) for a in DDZ.iter("Angle")}
    assert ("ddz-OG", "ddz-CB", "ddz-OG") in angles
