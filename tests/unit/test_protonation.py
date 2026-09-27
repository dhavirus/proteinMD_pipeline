"""Protonation states, pure parts (TASK-009): the disulfide and metal rules, pKa against
the pH, the two findings, and decisions. Structures are hand-built; pKa values are given."""

import pytest

from simprep.protonate.hydrogens import _partition
from simprep.protonate.states import Context, Estimate, ProtonationError, assign
from simprep.structure.model import ResidueId
from tests.unit.builders import NP, P, atom, config, link, residue, structure

A = {n: ResidueId("A", n) for n in (1, 2, 3, 4, 5, 6, 7, 900)}
METHOD = {"ambiguity_window_ph": 1.0, "unreliable_within_angstrom": 6.0}


def site():
    """ASP 1 and HIS 2 bound to CA 900 (NE2); CYS 3-4 disulfide; LYS 5 near the calcium;
    ASP 6 and GLU 7 far away."""

    def res(num, name, x):
        return residue("A", num, name, P, [atom("CA", "C", (x, 0, 0)), atom("CB", "C", (x, 1, 0))])

    residues = [
        res(1, "ASP", 0),
        res(2, "HIS", 2),
        res(3, "CYS", 30),
        res(4, "CYS", 32),
        res(5, "LYS", 5),
        res(6, "ASP", 60),
        res(7, "GLU", 90),
        residue("A", 900, "CA", NP, [atom("CA", "CA", (1, 2, 0))]),
    ]
    links = [
        link("metalc1", "metalc", ((A[1], "ASP", "OD1"), (A[900], "CA", "CA"))),
        link("metalc2", "metalc", ((A[2], "HIS", "NE2"), (A[900], "CA", "CA"))),
        link("disulf1", "disulf", ((A[3], "CYS", "SG"), (A[4], "CYS", "SG"))),
    ]
    return structure(residues, links=links)


def estimates(asp6=3.0, glu7=4.4):
    values = {
        1: ("ASP", 9.0),
        2: ("HIS", 9.0),
        3: ("CYS", 99.99),
        4: ("CYS", 99.99),
        5: ("LYS", 6.0),
        6: ("ASP", asp6),
        7: ("GLU", glu7),
    }
    model = {"ASP": 3.8, "HIS": 6.5, "CYS": 9.0, "LYS": 10.5, "GLU": 4.5}
    return [Estimate(A[n], name, pka, model[name]) for n, (name, pka) in values.items()]


def run(ruleset, decisions=None, **kw):
    context = Context(7.2, METHOD, ruleset, decisions or {}, "wt", config())
    states, findings = assign(site(), estimates(**kw), context)
    return {s.residue.seq_num: s for s in states}, {
        f.anchor_residues[0].seq_num: f for f in findings
    }


def test_rules_before_estimates(ruleset):
    states, findings = run(ruleset)
    assert (states[1].variant, states[2].variant) == ("ASP", "HID")  # pKa 9 ignored: metal
    assert states[3].variant == states[4].variant == "CYX"
    assert states[6].variant == "ASP" and states[7].variant == "GLU"  # pKa below 7.2
    assert set(findings) == {5}  # only the lysine near the calcium


def test_unreliable_estimate_defaults_to_the_model_pka(ruleset):
    states, findings = run(ruleset)
    assert findings[5].rule.rule_id == "protonation.unreliable_estimate"
    assert states[5].variant == "LYS"  # model pKa 10.5 > 7.2, although PROPKA said 6.0
    assert "recommended, undecided" in states[5].basis


def test_ambiguous_state_takes_the_prediction_unless_decided(ruleset):
    states, findings = run(ruleset, asp6=7.5)
    assert findings[6].rule.rule_id == "protonation.ambiguous_state"
    assert states[6].variant == "ASH"
    decision = {"option_id": "specific_state", "parameters": {"state": "ASP"}}
    states, _ = run(ruleset, {"protonation/wt/A:6": decision}, asp6=7.5)
    assert (states[6].variant, states[6].basis) == ("ASP", "specific_state (decided)")


def test_a_specific_state_the_residue_cannot_take_is_refused(ruleset):
    decision = {"option_id": "specific_state", "parameters": {"state": "HID"}}
    with pytest.raises(ProtonationError, match="one of"):
        run(ruleset, {"protonation/wt/A:6": decision}, asp6=7.5)


def test_chemistry_without_hydrogen_definitions_is_a_hard_stop():
    ligand = residue("A", 901, "XYZ", NP, [atom("C1", "C", (0, 0, 0)), atom("O1", "O", (1, 0, 0))])
    with pytest.raises(ProtonationError, match="no hydrogen definitions for XYZ A:901"):
        _partition(structure([ligand]), {})
