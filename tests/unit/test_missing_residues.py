import pytest

from simprep.detectors.missing_residues import detect_missing_residues
from simprep.structure.model import ResidueId, UnobservedResidue
from tests.unit import builders as b


def chain(observed_label_seqs, spacing=3.8):
    return [
        b.residue("A", n + 100, "GLY", b.P, [b.atom("CA", "C", (n * spacing, 0, 0))], n)
        for n in observed_label_seqs
    ]


def unobserved(label_seqs, chain_id="A"):
    return tuple(
        UnobservedResidue(ResidueId(chain_id, n + 100), "GLY", n, True) for n in label_seqs
    )


def findings_for(observed, missing, ruleset, spacing=3.8):
    site = b.structure(chain(observed, spacing), unobserved_residues=unobserved(missing))
    return {f.id: f for f in detect_missing_residues(site, ruleset, b.config())}


def test_terminal_and_internal_runs(ruleset):
    found = findings_for([3, 4, 5, 9, 10], [1, 2, 6, 7, 8, 11], ruleset)
    assert set(found) == {
        "missing_residues/A:101-102",
        "missing_residues/A:106-108",
        "missing_residues/A:111-111",
    }
    assert found["missing_residues/A:101-102"].rule.rule_id == "missing_residues.n_terminal"
    assert found["missing_residues/A:111-111"].rule.rule_id == "missing_residues.c_terminal"
    internal = found["missing_residues/A:106-108"]
    assert internal.rule.rule_id == "missing_residues.internal"
    assert internal.evidence_value("gap_length") == 3
    assert internal.evidence_value("flank_ca_distance") == pytest.approx(4 * 3.8)
    assert internal.evidence_value("bridgeable") is True
    assert internal.anchor_residues == (ResidueId("A", 105), ResidueId("A", 109))


def test_unbridgeable_gap_flagged(ruleset):
    found = findings_for([1, 2, 4], [3], ruleset, spacing=10.0)
    internal = found["missing_residues/A:103-103"]
    assert internal.evidence_value("bridgeable") is False


def test_whole_chain_unobserved(ruleset):
    site = b.structure([], unobserved_residues=unobserved([1, 2], "B"))
    (finding,) = detect_missing_residues(site, ruleset, b.config())
    assert finding.rule.rule_id == "missing_residues.whole_chain"


def test_nonpolymer_unobserved_ignored(ruleset):
    ligand = (UnobservedResidue(ResidueId("A", 900), "NAG", None, False),)
    site = b.structure(chain([1, 2]), unobserved_residues=ligand)
    assert detect_missing_residues(site, ruleset, b.config()) == []
