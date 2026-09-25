"""Severity escalation by proximity to manifest regions (acceptance criterion)."""

from simprep.audit import run_audit
from simprep.findings import raise_severity
from simprep.structure.model import ResidueId
from tests.unit import builders as b

ZINC = ResidueId("A", 900)


def zinc_with_distant_residue():
    """Catalytic zinc site plus an alanine 30 A away."""
    distant = b.residue("B", 1, "ALA", b.P, [b.atom("CA", "C", (30.0, 0, 0))], label_seq=1)
    return b.structure(b.zinc_catalytic_site().residues + (distant,))


def metal_finding(findings):
    return next(f for f in findings if f.rule.family == "metals")


def test_same_metal_finding_escalates_only_with_containing_region(ruleset):
    site = zinc_with_distant_residue()
    distant = run_audit(site, ruleset, b.config([b.region("far", [ResidueId("B", 1)])]))
    containing = run_audit(site, ruleset, b.config([b.region("site", [ZINC])]))
    far_finding, near_finding = metal_finding(distant), metal_finding(containing)
    assert far_finding.id == near_finding.id
    assert far_finding.rule.base_severity == near_finding.rule.base_severity == "warn"
    assert far_finding.effective_severity == "warn"
    assert far_finding.context.proximity == "far"
    assert near_finding.effective_severity == "blocking"
    assert near_finding.context.proximity == "inside"
    assert near_finding.context.min_distance_angstrom == 0.0


def test_recommendation_follows_region_context(ruleset):
    site = zinc_with_distant_residue()
    far = metal_finding(run_audit(site, ruleset, b.config([b.region("r", [ResidueId("B", 1)])])))
    inside = metal_finding(run_audit(site, ruleset, b.config([b.region("r", [ZINC])])))
    unset = metal_finding(run_audit(site, ruleset, b.config()))
    assert far.recommended_option == "nonbonded_12_6"
    assert inside.recommended_option == "bonded_mcpb"
    assert inside.recommendation_basis.startswith("condition catalytic_site_in_region")
    assert unset.recommended_option == "nonbonded_12_6_4"
    assert unset.recommendation_basis == "default"


def test_near_region_uses_tightest_threshold(ruleset):
    site = b.zinc_catalytic_site()
    his_residue = ResidueId("A", 10)
    finding = metal_finding(run_audit(site, ruleset, b.config([b.region("his", [his_residue])])))
    assert finding.context.proximity == "near"
    assert finding.context.min_distance_angstrom < 6.0
    assert finding.context.raise_levels == 2


def test_no_regions_means_effective_equals_base(ruleset):
    for finding in run_audit(zinc_with_distant_residue(), ruleset, b.config()):
        assert finding.effective_severity == finding.rule.base_severity
        assert finding.context.proximity == "no_regions"


def test_raise_severity_saturates():
    assert raise_severity("info", 1) == "warn"
    assert raise_severity("info", 2) == "blocking"
    assert raise_severity("warn", 2) == "blocking"
    assert raise_severity("blocking", 0) == "blocking"
