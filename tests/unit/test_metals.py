import pytest

from simprep.audit import run_audit
from simprep.detectors.metals import (
    IDEAL_ANGLES_DEGREE,
    classify,
    coordination_shell,
    detect_metals,
)
from simprep.structure.model import ResidueId
from tests.unit import builders as b


def only(findings):
    assert len(findings) == 1
    return findings[0]


def test_catalytic_zinc_shell_geometry_and_class(ruleset):
    finding = only(detect_metals(b.zinc_catalytic_site(), ruleset, b.config()))
    assert finding.id == "metals/A:900"
    assert finding.rule.rule_id == "metals.zinc"
    assert finding.evidence_value("coordination_number") == 4
    assert finding.evidence_value("protein_donors") == 3
    assert finding.evidence_value("water_donors") == 1
    assert finding.evidence_value("geometry") == "tetrahedral"
    assert finding.evidence_value("geometry_rms_angle_deviation") == pytest.approx(0.0, abs=1e-3)
    assert finding.evidence_value("mean_ligand_distance") == pytest.approx(2.05)
    assert finding.evidence_value("classification") == "likely_catalytic"
    assert finding.claims == {ResidueId("A", 900)}


def test_structural_zinc_cys4(ruleset):
    finding = only(detect_metals(b.zinc_structural_site(), ruleset, b.config()))
    assert finding.evidence_value("classification") == "likely_structural"
    assert finding.evidence_value("protein_donors") == 4


def test_ligand_outside_cutoff_is_excluded(ruleset):
    far_water = b.water("A", 1001, (0.0, 0.0, 2.9))  # zinc cutoff is 2.8 A
    site = b.structure(b.zinc_structural_site().residues + (far_water,))
    finding = only(detect_metals(site, ruleset, b.config()))
    assert finding.evidence_value("coordination_number") == 4


def test_closest_altloc_of_a_ligand_counts_once(ruleset):
    site = b.zinc_catalytic_site()
    residues = list(site.residues)
    water = residues[-1]
    shifted = b.atom("O", "O", b.offset(water.atoms[0].position, (0.2, 0, 0)), 0.5, altloc="B")
    first = b.atom("O", "O", water.atoms[0].position, 0.5, altloc="A")
    residues[-1] = b.residue("A", 1000, "HOH", b.W, [first, shifted])
    zinc_rule = next(r for r in ruleset.family("metals") if r.rule_id == "metals.zinc")
    shell = coordination_shell(
        b.structure(residues), residues[0], (zinc_rule, frozenset({"HIS", "HOH"}))
    )
    assert len(shell) == 4


def test_nonstandard_ligand_makes_site_catalytic(ruleset):
    positions = [b.scaled(d, 2.4) for d in b.TETRAHEDRAL]
    residues = [b.ion("A", 900, "CA")]
    residues += [
        b.residue("A", 50 + i, "ASP", b.P, [b.atom("OD1", "O", positions[i])]) for i in range(3)
    ]
    residues.append(b.residue("A", 79, "XYZ", b.P, [b.atom("O1", "O", positions[3])]))
    finding = only(detect_metals(b.structure(residues), ruleset, b.config()))
    assert finding.rule.rule_id == "metals.calcium"
    assert finding.evidence_value("classification") == "likely_catalytic"
    assert finding.evidence_value("nonstandard_polymer_residue_donors") == 1


def test_classification_ambiguous_when_no_criterion_met(ruleset):
    calcium = next(r for r in ruleset.family("metals") if r.rule_id == "metals.calcium")
    assert classify([], calcium.matcher["classification"]) == ("ambiguous", "no criterion met")


def test_annotated_metal_links_are_reported_verbatim(ruleset):
    site = b.zinc_catalytic_site()
    annotated = b.link(
        "metalc1",
        "metalc",
        [(ResidueId("A", 900), "ZN", "ZN"), (ResidueId("A", 10), "HIS", "NE2")],
        2.05,
    )
    finding = only(detect_metals(b.structure(site.residues, [annotated]), ruleset, b.config()))
    records = [e for e in finding.evidence if e["key"] == "annotated_metal_link"]
    assert records[0]["fields"]["id"] == "metalc1"


def test_polyatomic_metal_residue_is_not_an_ion(ruleset):
    heme_like = b.residue(
        "A", 500, "HEM", b.NP, [b.atom("FE", "FE", (0, 0, 0)), b.atom("NA", "N", (2.0, 0, 0))]
    )
    assert detect_metals(b.structure([heme_like]), ruleset, b.config()) == []


def test_metal_without_rule_is_unrecognized(ruleset):
    findings = run_audit(b.structure([b.ion("A", 900, "GD")]), ruleset, b.config())
    finding = only(findings)
    assert finding.rule.rule_id == "unrecognized.chemistry_group"
    assert finding.effective_severity == "blocking"


def test_ideal_angle_tables_have_all_pairs():
    for name, angles in IDEAL_ANGLES_DEGREE.items():
        n_pairs = len(angles)
        cn = next(n for n in range(2, 9) if n * (n - 1) // 2 == n_pairs)
        assert cn >= 2, name
