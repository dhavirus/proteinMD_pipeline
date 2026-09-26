from simprep.detectors.altlocs import detect_altlocs, max_spread
from tests.unit import builders as b


def two_state_residue(num, shift, occupancies=(0.6, 0.4)):
    return b.residue(
        "A",
        num,
        "SER",
        b.P,
        [
            b.atom("CA", "C", (0, 0, 0)),
            b.atom("OG", "O", (1, 0, 0), occupancies[0], 15.0, "A"),
            b.atom("OG", "O", (1 + shift, 0, 0), occupancies[1], 25.0, "B"),
        ],
    )


def test_small_spread_is_minor_info(ruleset):
    findings = detect_altlocs(b.structure([two_state_residue(3, 0.4)]), ruleset, b.config())
    (finding,) = findings
    assert finding.rule.rule_id == "altlocs.minor_spread"
    assert finding.rule.base_severity == "info"
    assert finding.evidence_value("occupancy_a") == 0.6
    assert finding.evidence_value("b_factor_b") == 25.0
    assert finding.evidence_value("max_spread") == 0.4


def test_large_spread_is_major_warn(ruleset):
    (finding,) = detect_altlocs(b.structure([two_state_residue(3, 2.5)]), ruleset, b.config())
    assert finding.rule.rule_id == "altlocs.major_spread"
    assert finding.locus.atom_names == ("OG",)


def test_waters_are_aggregated(ruleset):
    waters = [
        b.water("A", 1000 + i, (i * 5.0, 0, 0), altloc=alt, occupancy=0.5)
        for i in range(3)
        for alt in "AB"
    ]
    merged = [
        b.residue(
            "A",
            1000 + i,
            "HOH",
            b.W,
            [a for w in waters if w.id.seq_num == 1000 + i for a in w.atoms],
        )
        for i in range(3)
    ]
    (finding,) = detect_altlocs(b.structure(merged), ruleset, b.config())
    assert finding.id == "altlocs/water"
    assert finding.evidence_value("residues_with_altlocs") == 3
    assert len(finding.anchor_residues) == 3


def test_no_altlocs_no_findings(ruleset):
    assert detect_altlocs(b.zinc_catalytic_site(), ruleset, b.config()) == []


def test_spread_without_pairs_is_zero():
    lone = b.residue("A", 1, "SER", b.P, [b.atom("OG", "O", (0, 0, 0), 0.5, altloc="A")])
    assert max_spread(lone) == (0.0, None)
