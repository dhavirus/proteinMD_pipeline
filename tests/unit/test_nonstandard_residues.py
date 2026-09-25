from dataclasses import replace

from simprep.detectors.nonstandard_residues import detect_nonstandard_residues, modeled_form
from simprep.rules import Rule, RuleSet
from simprep.structure.model import ModifiedResidue, ResidueId
from tests.unit import builders as b


def mse(num=5):
    return b.residue(
        "A", num, "MSE", b.P, [b.atom("CA", "C", (0, 0, 0)), b.atom("SE", "SE", (2, 0, 0))]
    )


def ala(num):
    return b.residue("A", num, "ALA", b.P, [b.atom("CA", "C", (num * 3.8, 0, 0))])


def component_rule(ruleset):
    """A specific component rule for 'FGX' (synthetic id used only in this test)."""
    generic = ruleset.family("nonstandard_residues")[-1]
    return replace(
        generic,
        rule_id="nonstandard_residues.test_component",
        priority=10,
        matcher={
            "kind": "component_id",
            "comp_ids": ["FGX"],
            "residue_classes": ["polymer"],
            "modeled_forms": [
                {"form": "gem_diol", "comp_ids": ["FGX"], "required_atoms": ["O1", "O2"]},
                {
                    "form": "aldehyde",
                    "comp_ids": ["FGX"],
                    "required_atoms": ["O1"],
                    "absent_atoms": ["O2"],
                },
            ],
        },
    )


def with_rule(ruleset, rule: Rule) -> RuleSet:
    return replace(ruleset, rules=ruleset.rules + (rule,))


def test_nonstandard_by_name(ruleset):
    findings = detect_nonstandard_residues(
        b.structure([ala(4), mse(), ala(6)]), ruleset, b.config()
    )
    assert [f.id for f in findings] == ["nonstandard_residues/A:5"]
    assert findings[0].evidence_value("detected_by") == "residue_name"
    assert findings[0].recommended_option is None  # set later by the severity step


def test_nonstandard_by_annotation_and_parent_recorded(ruleset):
    annotation = ModifiedResidue(ResidueId("A", 5), "MSE", "MET", "SELENOMETHIONINE")
    site = b.structure([mse()], modified_residues=(annotation,))
    finding = detect_nonstandard_residues(site, ruleset, b.config())[0]
    assert finding.evidence_value("detected_by") == "residue_name,modification_annotation"
    record = next(e for e in finding.evidence if e["type"] == "source_record")
    assert record["fields"]["parent_comp_id"] == "MET"


def test_annotation_alone_flags_residue(ruleset):
    annotation = ModifiedResidue(ResidueId("A", 4), "ALA", "SER", "odd annotation")
    finding = detect_nonstandard_residues(
        b.structure([ala(4)], modified_residues=(annotation,)), ruleset, b.config()
    )[0]
    assert finding.evidence_value("detected_by") == "modification_annotation"


def test_attachment_site_annotation_is_not_a_nonstandard_residue(ruleset):
    annotation = ModifiedResidue(ResidueId("A", 4), "ASN", "ASN", "GLYCOSYLATION SITE")
    asn = b.residue("A", 4, "ASN", b.P, [b.atom("ND2", "N", (0, 0, 0))])
    site = b.structure([asn], modified_residues=(annotation,))
    assert detect_nonstandard_residues(site, ruleset, b.config()) == []


def test_standard_residues_and_ligands_not_flagged(ruleset):
    site = b.structure([ala(1), b.ion("A", 900, "ZN"), b.water("A", 1000, (9, 9, 9))])
    assert detect_nonstandard_residues(site, ruleset, b.config()) == []


def test_component_rule_wins_and_classifies_form(ruleset):
    rules = with_rule(ruleset, component_rule(ruleset))
    fgx = b.residue(
        "A",
        79,
        "FGX",
        b.P,
        [b.atom("CA", "C", (0, 0, 0)), b.atom("O1", "O", (1, 0, 0)), b.atom("O2", "O", (0, 1, 0))],
    )
    finding = detect_nonstandard_residues(b.structure([fgx]), rules, b.config())[0]
    assert finding.rule.rule_id == "nonstandard_residues.test_component"
    assert finding.evidence_value("modeled_form") == "gem_diol"


def test_modeled_form_aldehyde_and_unclassified(ruleset):
    matcher = component_rule(ruleset).matcher
    aldehyde = b.residue("A", 79, "FGX", b.P, [b.atom("O1", "O", (1, 0, 0))])
    bare = b.residue("A", 79, "FGX", b.P, [b.atom("CA", "C", (0, 0, 0))])
    assert modeled_form(aldehyde, matcher) == "aldehyde"
    assert modeled_form(bare, matcher) == "unclassified"


def test_revert_to_parent_is_never_recommended(ruleset):
    generic = ruleset.family("nonstandard_residues")[-1]
    explicit = [o["id"] for o in generic.options if o.get("requires_explicit_choice")]
    assert explicit == ["revert_to_parent"]
    assert generic.recommended_option["default"] != "revert_to_parent"
