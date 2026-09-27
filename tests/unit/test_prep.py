"""Prep plan and apply on hand-built structures: one test per operation, the gate,
conflicts, ensembles and the mapped revert (TASK-004)."""

from dataclasses import replace
from pathlib import Path

import pytest
import yaml

from simprep.model.loops import APPLIED_OPTIONS as MODEL_OPTIONS
from simprep.paths import KNOWLEDGE_DIR, SCHEMA_DIR
from simprep.prep.apply import apply_plan
from simprep.prep.plan import OPERATIONS, PrepError, build_plan
from simprep.prep.record import system_counts
from simprep.protonate.states import APPLIED_OPTIONS as PROTONATION_OPTIONS
from simprep.structure.model import ModifiedResidue, ResidueId
from simprep.variants.findings import FAMILY as VARIANT_FAMILY
from simprep.variants.run import APPLIED_OPTIONS
from tests.unit.builders import BR, NP, P, W, atom, link, residue, structure

A1, A2, A3, L1, L2 = (ResidueId("A", n) for n in (1, 2, 3, 901, 902))
OPTION_ACTIONS = {
    "highest_occupancy": ("apply", None),
    "specific_altloc": ("apply", None),
    "keep_ensemble": ("apply", None),
    "exclude": ("apply", None),
    "exclude_group": ("apply", None),
    "add_link": ("apply", None),
    "revert_to_parent": ("apply", None),
    "model_gem_diol": ("apply", None),
    "treat_noncovalent": ("record", "topology"),
    "parameterize_manually": ("defer", "parameterization"),
    "expert_review": ("unresolved", None),
}


def split(name, element, position, occupancies):
    """One atom in two altlocs (A, B) with the given occupancies."""
    return [
        atom(name, element, position, occ, altloc=alt)
        for alt, occ in zip("AB", occupancies, strict=True)
    ]


def sample():
    """Polymer A:1-3 (A:2 and A:3 with altlocs), ligands A:901-902, water A:1001."""
    return structure(
        [
            residue("A", 1, "MSE", P, [atom("CA", "C", (0, 0, 0)), atom("SE", "SE", (1, 0, 0))]),
            residue(
                "A",
                2,
                "SER",
                P,
                [atom("CA", "C", (3, 0, 0)), *split("OG", "O", (4, 0, 0), (0.4, 0.6))],
            ),
            residue(
                "A",
                3,
                "ASN",
                P,
                [atom("CA", "C", (6, 0, 0)), *split("ND2", "N", (7, 0, 0), (0.7, 0.3))],
            ),
            residue("A", 901, "NAG", BR, [atom("C1", "C", (8, 0, 0))]),
            residue("A", 902, "NAG", BR, [atom("C1", "C", (9, 0, 0))]),
            residue("A", 1001, "HOH", W, [atom("O", "O", (20, 0, 0))]),
        ],
        links=[link("covale1", "covale", ((L1, "NAG", "O4"), (L2, "NAG", "C1")))],
    )


def ref(rid):
    return {"chain": rid.chain, "seq_num": rid.seq_num, "ins_code": rid.ins_code}


def finding(fid, residues, severity="blocking", evidence=()):
    options = [
        {"id": o, "label": o, "prep_action": a, **({"prep_stage": s} if s else {})}
        for o, (a, s) in OPTION_ACTIONS.items()
    ]
    return {
        "id": fid,
        "rule_id": f"test.{fid}",
        "title": f"title {fid}",
        "effective_severity": severity,
        "options": options,
        "locus": {"extent": [ref(r) for r in residues]},
        "evidence": list(evidence),
    }


def manifest(findings, decisions):
    return {
        "findings_snapshot": {"findings": findings},
        "decisions": [
            {
                "finding_id": f,
                "option_id": o,
                "rationale": "t",
                "decided_by": "t",
                "timestamp": "2026-09-26T00:00:00Z",
                **({"parameters": p} if p else {}),
            }
            for f, o, p in decisions
        ],
    }


def plan_for(findings, decisions, ruleset):
    return build_plan(sample(), manifest(findings, decisions), ruleset)


def only_system(plan):
    (system,) = apply_plan(sample(), plan)
    return system.structure


def test_highest_occupancy_keeps_the_major_altloc_and_clears_ids(ruleset):
    plan = plan_for([finding("alt", [A2, A3])], [("alt", "highest_occupancy", None)], ruleset)
    out = only_system(plan).residue_index
    assert [(a.name, a.altloc, a.occupancy) for a in out[A2].atoms[1:]] == [("OG", "", 0.6)]
    assert [(a.name, a.altloc, a.occupancy) for a in out[A3].atoms[1:]] == [("ND2", "", 0.7)]


def test_specific_altloc_keeps_the_named_altloc(ruleset):
    plan = plan_for([finding("alt", [A2])], [("alt", "specific_altloc", {"altloc": "A"})], ruleset)
    assert only_system(plan).residue_index[A2].atoms[1].occupancy == 0.4


def test_specific_altloc_refuses_an_absent_altloc(ruleset):
    with pytest.raises(PrepError, match="no altloc 'C'"):
        plan_for([finding("alt", [A2])], [("alt", "specific_altloc", {"altloc": "C"})], ruleset)


def test_keep_ensemble_gives_one_system_per_altloc(ruleset):
    plan = plan_for([finding("alt", [A2, A3])], [("alt", "keep_ensemble", None)], ruleset)
    systems = apply_plan(sample(), plan)
    assert [s.name for s in systems] == ["system_altloc_A", "system_altloc_B"]
    assert [s.structure.residue_index[A3].atoms[1].occupancy for s in systems] == [0.7, 0.3]


def test_keep_ensemble_needs_every_altloc_on_every_residue(ruleset):
    three = replace(
        sample(),
        residues=sample().residues[:2]
        + (residue("A", 3, "ASN", P, [atom("ND2", "N", (7, 0, 0), 1.0, altloc="C")]),)
        + sample().residues[3:],
    )
    with pytest.raises(PrepError, match="A:3 has no altloc A, B"):
        build_plan(
            three, manifest([finding("alt", [A2, A3])], [("alt", "keep_ensemble", None)]), ruleset
        )


def test_exclude_removes_residues_and_their_links(ruleset):
    plan = plan_for([finding("lig", [L1])], [("lig", "exclude", None)], ruleset)
    out = only_system(plan)
    assert L1 not in out.residue_index and out.links == ()
    rows = {r["record_type"]: r for r in system_counts(sample(), out)}
    assert rows["branched residues"]["excluded"] == 1
    assert rows["struct_conn records"]["excluded"] == 1


def test_exclude_group_removes_the_whole_ligand_group(ruleset):
    contact = finding(
        "contact", [L1, A3], evidence=[{"key": "ligand_group_finding", "value": "group"}]
    )
    group = finding("group", [L1, L2], severity="warn")
    plan = plan_for([contact, group], [("contact", "exclude_group", None)], ruleset)
    assert plan.excluded == {L1, L2}


def contact_finding():
    atoms = [
        {
            "chain": "A",
            "seq_num": 902,
            "ins_code": "",
            "res_name": "NAG",
            "atom_name": "C1",
            "altloc": "",
        },
        {
            "chain": "A",
            "seq_num": 3,
            "ins_code": "",
            "res_name": "ASN",
            "atom_name": "ND2",
            "altloc": "A",
        },
    ]
    return finding(
        "contact", [L2, A3], evidence=[{"key": "contact_distance", "value": 2.4, "atoms": atoms}]
    )


def test_add_link_writes_a_covalent_link_with_collapsed_altlocs(ruleset):
    plan = plan_for(
        [contact_finding(), finding("alt", [A3])],
        [("contact", "add_link", None), ("alt", "highest_occupancy", None)],
        ruleset,
    )
    added = only_system(plan).links[-1]
    assert (added.conn_id, added.conn_type) == ("prep_link1", "covale")
    assert (added.partner2.residue, added.partner2.atom_name, added.partner2.altloc) == (
        A3,
        "ND2",
        "",
    )


def test_add_link_to_a_dropped_altloc_is_a_conflict(ruleset):
    with pytest.raises(PrepError, match="bonds altloc A, but the altloc decision keeps B"):
        plan_for(
            [contact_finding(), finding("alt", [A3])],
            [("contact", "add_link", None), ("alt", "specific_altloc", {"altloc": "B"})],
            ruleset,
        )


def test_linked_and_excluded_residue_is_a_conflict(ruleset):
    with pytest.raises(PrepError, match="A:902 is excluded by one decision and linked"):
        plan_for(
            [contact_finding(), finding("lig", [L2])],
            [("contact", "add_link", None), ("lig", "exclude", None)],
            ruleset,
        )


def test_deferred_and_excluded_residue_is_a_conflict(ruleset):
    with pytest.raises(PrepError, match="keeps A:901, which another decision excludes"):
        plan_for(
            [finding("a", [L1]), finding("b", [L1])],
            [("a", "exclude", None), ("b", "parameterize_manually", None)],
            ruleset,
        )


def test_revert_to_parent_refuses_without_a_mapping(ruleset):
    with pytest.raises(PrepError, match="needs an atom mapping for MSE"):
        plan_for([finding("ns", [A1])], [("ns", "revert_to_parent", None)], ruleset)


def test_revert_to_parent_applies_an_explicit_mapping(ruleset):
    example = SCHEMA_DIR / "examples" / "residue_mappings" / "valid_mse_to_met.yaml"
    mapped = replace(
        ruleset, residue_mappings=tuple(yaml.safe_load(example.read_text())["mappings"])
    )
    plan = plan_for([finding("ns", [A1])], [("ns", "revert_to_parent", None)], mapped)
    reverted = only_system(plan).residue_index[A1]
    assert reverted.name == "MET"
    assert [(a.name, a.element) for a in reverted.atoms] == [("CA", "C"), ("SD", "S")]


def test_record_and_defer_become_work_items_and_change_nothing(ruleset):
    plan = plan_for(
        [finding("a", [L1]), finding("b", [A3], severity="warn")],
        [("a", "parameterize_manually", None)],
        ruleset,
    )
    assert only_system(plan) == sample()
    assert [(w.stage, w.finding_id) for w in plan.work_order] == [
        ("parameterization", "a"),
        ("decision", "b"),
    ]


def test_gate_names_undecided_blocking_findings(ruleset):
    with pytest.raises(PrepError, match=r"\[blocking\] lig"):
        plan_for([finding("lig", [L1])], [], ruleset)


def test_gate_refuses_any_non_final_decision(ruleset):
    with pytest.raises(PrepError, match="lig: expert_review is not a final decision"):
        plan_for([finding("lig", [L1], severity="info")], [("lig", "expert_review", None)], ruleset)


def test_snapshot_without_prep_action_is_refused(ruleset):
    old = finding("lig", [L1])
    old["options"] = [
        {k: v for k, v in o.items() if not k.startswith("prep_")} for o in old["options"]
    ]
    with pytest.raises(PrepError, match="predates knowledge base 0.3.0"):
        plan_for([old], [("lig", "exclude", None)], ruleset)


def apply_options(prep_families: bool, stage: str | None = None) -> set[str]:
    """Apply options of prep's rule families (or the variant family) carried out by
    ``stage`` (None: by prep itself)."""
    return {
        option["id"]
        for path in KNOWLEDGE_DIR.glob("*.yaml")
        for rule in (yaml.safe_load(Path(path).read_text()) or {}).get("rules", [])
        if (rule["family"] != VARIANT_FAMILY) == prep_families
        for option in rule["options"]
        if option.get("prep_action") == "apply" and option.get("prep_stage") == stage
    }


def test_every_apply_option_in_the_knowledge_base_has_an_operation():
    assert apply_options(prep_families=True) and apply_options(True) <= set(OPERATIONS)
    assert apply_options(prep_families=True, stage="modelling") == set(MODEL_OPTIONS)
    assert apply_options(prep_families=True, stage="protonation") == set(PROTONATION_OPTIONS)
    assert apply_options(prep_families=False) == set(APPLIED_OPTIONS)


def gem_diol_site():
    """VAL A:83 - ALS A:84 (sulfate on OS1) - CA A:1551 bound to OS1 and OS4."""
    als = residue(
        "A",
        84,
        "ALS",
        P,
        [
            atom(name, element, (float(i), 0, 0))
            for i, (name, element) in enumerate(
                [
                    ("N", "N"),
                    ("CA", "C"),
                    ("C", "C"),
                    ("O", "O"),
                    ("CB", "C"),
                    ("OG", "O"),
                    ("OS1", "O"),
                    ("S", "S"),
                    ("OS2", "O"),
                    ("OS3", "O"),
                    ("OS4", "O"),
                ]
            )
        ],
        label_seq=2,
    )
    val = residue("A", 83, "VAL", P, [atom("C", "C", (-1, 0, 0))], label_seq=1)
    ca = residue("A", 1551, "CA", NP, [atom("CA", "CA", (6, 2, 0))])
    rid, metal = ResidueId("A", 84), ResidueId("A", 1551)
    return structure(
        [val, als, ca],
        links=[
            link("covale1", "covale", ((ResidueId("A", 83), "VAL", "C"), (rid, "ALS", "N"))),
            link("metalc3", "metalc", ((rid, "ALS", "OS1"), (metal, "CA", "CA"))),
            link("metalc5", "metalc", ((rid, "ALS", "OS4"), (metal, "CA", "CA"))),
        ],
        modified_residues=(ModifiedResidue(rid, "ALS", "ALA"),),
        polymer_sequences=(("A", ("VAL", "ALS")),),
    )


def test_model_gem_diol_turns_als_into_ddz_with_links_record_and_sequence(ruleset):
    site, rid = gem_diol_site(), ResidueId("A", 84)
    manifest_ = manifest([finding("fgly", [rid])], [("fgly", "model_gem_diol", None)])
    plan = build_plan(site, manifest_, ruleset)
    assert "link metalc5 to CA A:1551 removed (metal coordination changed" in plan.actions[0].note
    (system,) = apply_plan(site, plan)
    out = system.structure
    ddz = out.residue_index[rid]
    assert ddz.name == "DDZ"
    assert [a.name for a in ddz.atoms] == ["N", "CA", "C", "O", "CB", "OG1", "OG2"]
    assert [a.position for a in ddz.atoms] == [
        a.position
        for a in site.residue_index[rid].atoms
        if a.name not in ("S", "OS2", "OS3", "OS4")
    ]
    links = {k.conn_id: k for k in out.links}
    assert set(links) == {"covale1", "metalc3"}
    assert (links["metalc3"].partner1.res_name, links["metalc3"].partner1.atom_name) == (
        "DDZ",
        "OG2",
    )
    assert links["covale1"].partner2.res_name == "DDZ"
    assert out.modified_residues == (ModifiedResidue(rid, "DDZ", "ALA"),)
    assert out.polymer_sequences == (("A", ("VAL", "DDZ")),)
